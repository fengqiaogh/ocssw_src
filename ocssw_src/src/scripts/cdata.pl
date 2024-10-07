#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';

use Getopt::Long qw(GetOptionsFromArray :config bundling no_auto_abbrev);
use Pod::Usage;

use Fcntl qw(:seek);
use Cwd qw(cwd abs_path);
use File::Basename qw(dirname basename);
use File::Copy qw(mv);
use File::Spec qw();
use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);
use File::Glob qw(bsd_glob);
use File::Find qw(find);
use Digest::SHA;
use List::MoreUtils qw(natatime);

use Data::Dumper;

BEGIN {
	$ENV{RSYNC_RSH} //= 'ssh -q';
};

my %commands = (
	'remove' => \&_remove,
	'add' => \&_add,
	'list' => \&_list,
	'push' => \&_push,
	'init' => \&_init,
	'download' => \&_download,
	'clean' => \&_clean,
);

my %aliases = (
	rm => 'remove',
	dl => 'download',
);

sub cdata {
	my @args = @_;
	if (!@args){
		pod2usage(verbose => 0); return 1;
	}
	my %options = (
		verbose => 1,
	        "rsync-path" => [
		        'analysis701:/gfs-oceanweb/www-data/oceandata.sci.gsfc.nasa.gov/htdocs/ocssw_test',
		],
		tabs_to_add => 0,
		'remote-default-from-base' => 1,
	);

	for my $file ("$ENV{HOME}/.cdata"){
		if (-f $file){
			my $do_ret = do $file;
			if ($@){
				print "Couldn't parse config $file\n";
				return 1;
			} elsif (!defined($do_ret)){
				print "Couldn't do $file\n";
				return 1;
			} elsif (!$do_ret) {
				print "Couldn't run $file\n";
				return 1;
			} elsif (!ref($do_ret) || ref($do_ret) ne 'HASH'){
				print "Hash expected in config $file\n";
				return 1;
			}
			my %c = %$do_ret;
			@options{keys(%c)} = values(%c);
		}
	}

	GetOptionsFromArray(\@args, \%options, 'man', 'help|h|?', 'debug|d',
		'verbose|v' => sub {$options{verbose}++}, 'quiet|q' => sub {$options{verbose}--},
		'id=s', 'force|f', 'create|c', 'out-file|out|o=s', 'log-file|l=s',
		'rsync-path=s@', 'remote-default-to-basename', 'remote-default-from-local',
		'remote-default-from=s', 'remote-default-from-base', 'cmake-base=s', 'base-dir=s', 'meta-file|m=s',
		'continue-on-failure',
	) || pod2usage(2);
	if ($options{'log-file'}){
		$options{'close-log-file'} = 1;
		open(my $fh, '>', $options{'log-file'}) || die("Could not open meta file: $!");
		$options{'log-file'} = $fh;
	} else {
		$options{'log-file'} = \*STDOUT;
	}
	if ($options{'help'}){
		pod2usage(output => $options{'log-file'}, noperldoc => 1, verbose => 99, sections => 'SYNOPSIS|OPTIONS');
		return clean(\%options, 1);
	} elsif ($options{'man'}){
		pod2usage(output => $options{'log-file'}, noperldoc => 1, verbose => 99, sections => 'NAME|SYNOPSIS|OPTIONS|DESCRIPTION');
		return clean(\%options, 1);
	}

	my ($search_path, @meta_files);
	if (-f 'cmake_install.cmake'){
		open(my $f, '<', 'cmake_install.cmake');
		if ($f){
			($search_path = <$f>) =~ s/^.*?: |\s+$//g;
			close($f);
		}
	} else {
		$search_path = cwd();
	}
	unless ($options{'meta-file'}){
		find(sub {/^CDataList.txt/ && !/tmp$/ && push(@meta_files, $File::Find::name)}, $search_path);
		my ($cur, $cmake_base) = (dirname($search_path));
		while (-f catfile($cur, 'CMakeLists.txt')){
			push(@meta_files, bsd_glob(catfile($cur, 'CDataList.txt*')));
			$cmake_base = $cur;
			$cur = dirname($cur);
		}
		if ($options{'remote-default-from-base'} && !$options{'cmake-base'}){
			$options{'cmake-base'} = $cmake_base;
		}
	}
	if ($options{'remote-default-from-base'} && !$options{'cmake-base'}){
		my ($cur, $cmake_base) = (dirname($search_path));
		while (-f catfile($cur, 'CMakeLists.txt')){
			$cmake_base = $cur;
			$cur = dirname($cur);
		}
		$options{'cmake-base'} = $cmake_base;
	}

	my $cmd = shift(@args);
	if (!$cmd){
		pod2usage(output => $options{'log-file'}, verbose => 1, message => "No command given");
		return clean(\%options, 1);
	}

	if ($aliases{$cmd}){
		$cmd = $aliases{$cmd};
	}
	if (!$commands{$cmd}){
		pod2usage(output => $options{'log-file'}, verbose => 1, message => "Command '$cmd' not understood");
		return clean(\%options, 1);
	}

	my $ret = 1;
	if (!$options{'id'}){
		$options{'id'} = `git rev-parse HEAD`;
		chomp($options{'id'});
	}
	if ($cmd eq 'init'){
		$ret = $commands{$cmd}(\%options, undef, $options{'meta-file'} || shift(@args), \@args);
	} elsif (!$options{'meta-file'} && !@meta_files){
		pod2usage(output => $options{'log-file'}, noperldoc => 1, verbose => 1, message => "Metadata file not given and none were found in $search_path");
		return clean(\%options, 1);
	} elsif ($options{'meta-file'}){
		my $meta;
		if (!-f $options{'meta-file'}){
			if ($options{'create'}){
				$meta = [];
			} else {
				p(\%options, 0, "Metadata file '$options{'meta-file'}' does not exist, use --create to override\n");
				return clean(\%options, 1);
			}
		}
		$meta ||= read_meta(\%options, $options{'meta-file'});

		if ($meta){
			$options{'base-dir'} ||= dirname($options{'meta-file'});
			$ret = $commands{$cmd}(\%options, $meta, $options{'meta-file'}, \@args);
		} else {
			p(\%options, 0, "Metadata file invalid\n");
		}
	} else {
		my ($remote_from, %file_map) = ($options{'remote-default-from'});
		for my $meta_file (@meta_files){
			my $meta = read_meta(\%options, $meta_file);

			if ($meta){
				if ($cmd =~ /add|remove/){
					my @files;
					for (my $i=0;$i<@args;$i++){
						if ($meta->{$args[$i]}){
							push(@files, splice(@args, $i--, 1));
							last;
						}
					}
					$file_map{$meta_file} = [$meta, \@files];
				} else {
					$file_map{$meta_file} = [$meta, undef];
				}
			} else {
				p(\%options, 0, "Metadata file invalid ($meta_file)\n");
			}
		}
		if (@args){
			if ($cmd eq 'remove'){
				p(\%options, 0, "Some files not found (bailing out):\n\t" . join("\n\t", @args) . "\n");
				return clean(\%options, 1);
			} elsif ($cmd eq 'add'){
				push(@{$file_map{$meta_files[0]}[1]}, @args);
			}
		}

		$options{tabs_to_add} = 2;
		if (!$options{verbose}){
			p(\%options, 0, "Relevant metadata files found:\n\t" . join("\n\t", keys(%file_map)) . "\n");
		} elsif ($options{verbose} > 0){
			p(\%options, 1, "Relevant metadata files found:\n");
			while (my ($meta_file, $info) = each(%file_map)){
				my $args = $info->[1];
				if ($cmd =~ /add|remove/){
					p(\%options, 1, "\t$meta_file:\n\t\t" . join("\n\t\t", @$args) . "\n");
				} else {
					p(\%options, 1, "\t$meta_file\n");
				}
			}
		}
		$options{'base-dir'} ||= '.';
		while (my ($meta_file, $info) = each(%file_map)){
			my ($meta, $args) = @$info;
			$ret = $commands{$cmd}(\%options, $meta, $meta_file, $args);
		}
	}

	return clean(\%options, $ret);
}

sub clean {
	my ($options, $ret) = @_;
	if ($options->{'close-log-file'}){
		close($options->{'log-file'});
	}
	return $ret;
}

sub p {
	my $options = shift;
	my $level = shift;
	if ($options->{verbose} >= $level){
		print {$options->{'log-file'}} join('', @_);
	}
}

sub read_meta {
	my ($options, $file) = @_;
	my %meta;
	open(my $fh, "<", $file) || die("Could not open meta file: $!");
	while (my $line = <$fh>){
		chomp($line);
		if ($line !~ /^\s*$|^#/){
			if ($line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+([a-f0-9]+)$/){
				$meta{$1} = {'local' => $1, 'remote' => $2, 'id' => $3, 'hash' => $4};
			} else {
				return;
			}
		}
	}
	close($fh);
	return \%meta;
}

sub write_meta {
	my ($options, $meta, $filename) = @_;
	my ($need_to_close, $fh) = (1);
	if ($options->{'out-file'} && $options->{'out-file'} eq '-'){
		$fh = $options->{'log-file'};
		$need_to_close = 0;
	} elsif ($options->{'out-file'}){
		if ($options->{debug} && $options->{'out-file'} eq $filename){
			print "--debug used, can't overwrite original file\n";
			return 1;
		}
		open($fh, ">", $options->{'out-file'}) || die("Could not open meta file: $!");
	} elsif ($options->{debug}){
		open($fh, ">", "$filename.tmp") || die("Could not open meta file: $!");
	} else {
		open($fh, ">", $filename) || die("Could not open meta file: $!");
	}

	print {$fh} "# local-path remote-path git-rev sha1-hash\n";
	for (sort(keys(%$meta))){
		my $e = $meta->{$_};
		print {$fh} "$e->{local} $e->{remote} $e->{id} $e->{hash}\n";
	}
	if ($need_to_close){
		close($fh);
	}
	return;
}

sub _remove {
	my ($options, $meta, $meta_filename, $files) = @_;
	my (@new_files, %remote_locations);
	for (@$files){
		delete($meta->{File::Spec->abs2rel($_, $options->{'base-dir'})});
	}
	return write_meta($options, $meta, $meta_filename);
}

sub _add {
	my ($options, $meta, $meta_filename, $files) = @_;

	for (@$files){
		my ($local, $remote) = ($_);
		if ($local =~ /^(.+?)=(.+?)$/){
			($local, $remote) = ($1, $2);
		}
		if (-d $local){
			p($options, 0, "$local is a directory, skipping\n");
			next;
		}
		unless (-f $local && -r $local){
			p($options, 0, "Can't read file $_\n");
			return 1;
		}
		if (-s $local == 0 && !$options->{force}){
			p($options, 0, "$local has a size of 0, use --force to continue adding\n");
			return 1;
		}
		my $canonical = File::Spec->abs2rel($local, $options->{'base-dir'});
		unless ($remote){
			$remote = $meta->{$canonical}{remote};
			unless ($remote){
				if ($options->{'remote-default-to-basename'}){
					$remote = basename($canonical);
				} elsif ($options->{'remote-default-from-local'}){
					$remote = $canonical;
				} elsif ($options->{'remote-default-from'}){
					$remote = File::Spec->abs2rel($local, $options->{'remote-default-from'});
				} elsif ($options->{'remote-default-from-base'}){
					$remote = File::Spec->abs2rel($local, $options->{'cmake-base'});
				}
			}
		}
		unless ($remote){
			p($options, 0, "No remote path known for $canonical\n");
			return 1;
		}
		my $sha = Digest::SHA->new('SHA1');
		$sha->addfile($local);
		my $digest = $sha->hexdigest;

		if (!$meta->{$canonical}{hash} || $meta->{$canonical}{hash} ne $digest || $options->{force}){
			$meta->{$canonical} = {
				'id' => $options->{id},
				'local' => $canonical, 'remote' => $remote,
				'hash' => $digest,
			};
		} else {
			p($options, 1, "Hash unchanged for $canonical, file skipped.  Use --force to update it anyway.\n");
		}
	}
	return write_meta($options, $meta, $meta_filename);
}

sub _push {
	my ($options, $meta, $meta_filename, $files) = @_;
	my $id = $options->{id};
	if ($id){
		my @dests = @{$options->{'rsync-path'}};
		my $found = 0;
		for my $dest (@dests){
			my ($dest_path, $dest_host) = ($dest, $dest);
			my $remote_is_local = 1;
			if ($dest_path =~ /:/){
				$dest_path =~ s/^.*?://;
				$dest_host =~ s/:.*?$//;
				$remote_is_local = 0;
			}
			my $src = File::Spec->abs2rel($options->{'base-dir'});
			$dest =~ s"([^/])$"$1/"g;
			while (my ($local, $file) = each(%$meta)){
				if ($file->{id} eq $id){
					$found++;
					p($options, 1, "Copying $src/$file->{local} => $dest$id/$file->{remote}\n");
					my ($remote_dir, @cmd) = (dirname("$dest_path/$id/$file->{remote}"));
					if ($remote_is_local){
						make_path($remote_dir, {mode => 0775});
						p($options, 2, "mkdir $remote_dir\n");
						@cmd = ("rsync", '-haq', "$src/$file->{local}", "$dest$id/$file->{remote}");
					} else {
						@cmd = ("rsync", '-haq', '--rsync-path', "mkdir -p -m 775 '$remote_dir' && rsync", "$src/$file->{local}", "$dest$id/$file->{remote}");
					}
					p($options, 2, join(' ', @cmd) . "\n");
					system(@cmd);
					if ($? != 0){
						p($options, 0, "rsync to $dest failed\n");
					} 
				}
			}
		}
		unless ($found){
			p($options, 0, "No files found for ID $id\n");
		}

	} else {
		p($options, 0, "No ID found to be synced, use --id\n");
	}
	return;
}
sub _init {
	my ($options, $meta, $meta_filename, $files) = @_;
	my $dest = shift(@$files) || '.';

	if (-f $meta_filename){
		$meta_filename = dirname($meta_filename);
	}
	$meta_filename = abs_path($meta_filename);
	$dest =~ s"/$"";

	open(my $f, '>', "$dest/cmake_install.cmake") || die("Couldn't create cmake_install.cmake file");
	print {$f} "# Install script for directory: $meta_filename\n";
	print {$f} "# This is fake, generated by CData to point to who created it\n";
	close($f);

	# create the metadata file if it does not exist
	$meta_filename = "$meta_filename/CDataList.txt";
	if(! -f $meta_filename) {
		open($f, '>', "$meta_filename") || die("Couldn't create CDataList.txt file");
		print {$f} "# local-path remote-path git-rev sha1-hash\n";
		close($f);
	}

	return;
}

sub _list {
	my ($options, $meta, $meta_filename, $files) = @_;

	for (sort(keys(%$meta))){
		my $e = $meta->{$_};
		p($options, 0, "\t"x$options->{tabs_to_add} . "$e->{local} $e->{id}/$e->{remote} sha1:$e->{hash}\n");
	}

	return;
}

sub _clean {
	my ($options, $meta, $meta_filename, $files) = @_;

	my $src = File::Spec->abs2rel($options->{'base-dir'});

	my %should_exist;
	for (sort(keys(%$meta))){
		my $e = $meta->{$_};
		$should_exist{"$src/$e->{local}"} = 1;
	}

	my @should_delete;
	find({
		no_chdir => 1,
		wanted => sub {
			my $f = $File::Find::name;
			if ($f !~ /cmake_install.cmake$/ && -f $f && !$should_exist{$f}){
				push(@should_delete, $f);
			}
		}
	}, $src);

	for (@should_delete){
		p($options, 0, "Deleting $_\n");
		unless ($options->{debug}){
			unlink($_);
		}
	}
	return;
}

sub _download {
	my ($options, $meta, $meta_filename, $files) = @_;

	my $src = File::Spec->abs2rel($options->{'base-dir'});
	my $dest = $options->{'rsync-path'}[0];
	$dest =~ s"([^/])$"$1/"g;
	my @failed;

	my $i=0;
	for (sort(keys(%$meta))){
		my $e = $meta->{$_};
		#p($options, 0, "\t"x$options->{tabs_to_add} . "$e->{local} $e->{id}/$e->{remote} sha1:$e->{hash}\n");
		my @cmd = ("rsync", '-haq', "$dest$e->{id}/$e->{remote}", "$src/$e->{local}");
		make_path(dirname("$src/$e->{local}"), {mode => 0775});
		p($options, 1, "Downloading $src/$e->{local}...\n");
		if ($options->{debug}){
			my @cmd = ("rsync", '-haq', "$dest$e->{id}/$e->{remote}", "$src/$e->{local}");
			p($options, 0, join(" ", @cmd) . "\n");
		} else {
			system(@cmd);
			if ($? != 0){
				p($options, 0, "Download to $src/$e->{local} failed\n");
				push(@failed, "$src/$e->{local}");
				unless ($options->{'continue-on-failure'}){
					p($options, 0, "Use --continue-on-failure to keep going anyway\n");
					last;
				}
			}
			$i++;
			if ($i % 10 == 0){
				p($options, 1, "Sleeping for a sec to allow Ctrl-C...\n");
				sleep(1);
			}
		}
	}

	if ($options->{'continue-on-failure'} && @failed){
		p($options, 0, "The following files failed to download:\n\t");
		p($options, 0, join("\n\t", @failed) . "\n");
	}

	return;
}

unless (caller){
	exit(cdata(@ARGV) || 0);
}

__END__

=head1 NAME

cdata - CData list file manipulator

=head1 SYNOPSIS

cdata [options] add <file>[=<remote path>] ...

cdata [options] remove|rm <file> ...

cdata [options] push

cdata [options] list

cdata [options] init <CDataList.txt> [<destination>]

cdata [options] download [--continue-on-failure]

cdata [options] clean

=head1 OPTIONS

=over 8

=item B<--help>,B<-h>,B<-?>

Prints a brief help message and exits

=item B<--meta-file>,B<-m>

Meta-file to update; if not found, searches in first line of cmake_install.cmake or $PWD if not found

=item B<--create>,B<-c>

Create metadata file if it doesn't exist

=item B<--remote-default-to-basename>

Remote paths default to <commit>/<basename> (flattened structure)

=item B<--remote-default-from-local>

Remote paths default to <commit>/<local-path-from-meta-file>

=item B<--base-dir>

What to consider the local directory

=item B<--remote-default-from>=DIR

Remote paths default to <commit>/<local-path-from-DIR>

=item B<--remote-default-from-base>

Remote paths default to <commit>/<local-path-from-CMake-base>

=item B<--cmake-base>

Force CMake base, defaults to climbing directory tree until no more CMakeList.txt files are found

=item B<--rsync-path>

Path (optionally, with host) to rsync added files to

=item B<--id>=ID

Don't autodetect current Git commit

=item B<--log-file>,B<-l>

Change status message location, defaults to STDOUT

=item B<--out-file>,B<--out>,B<-o>

Change output file; set to - to output to B<--log-file>

=item B<--verbose>,B<-v>,B<--quiet>,B<-q>

Adjust run verbosity

=item B<--debug>,B<-d>

Almost useless, prevents overwriting original file (sometimes)

=item B<--force>,B<-f>

Allow weird configurations

=back

=head1 DESCRIPTION

This program is used to modify CData metadata files, whose format is described below.  The
file is meant to chronicle which data files have changed throughout Git revisions for the
purposes of downloading the data files just before running a test and no sooner than required
and only if the actual test using them is run.  There is an accompanying CMake module that reads
them and generates the download commands (or can run them at CMake configure time) and, optionally,
makes a CMake script to be used as the test instead of the binary to be run (thus downloading at
test time).

=head2 FILE FORMAT

local-path remote-path git-rev sha1-hash

=head2 NOTES

Spaces in file paths are not supported at this time.

=head1 WORKFLOW

The files to be added require a Git revision before they can be added.  So, first commit the file
changes, use this to modify the CData files, commit the CData files, then push both commits at once.
This isn't strictly necessary anymore; if bypassed, files will be tagged with an "incorrect" Git
revision, but this revision is no longer useful, it might as well be a random UUID.

=cut
