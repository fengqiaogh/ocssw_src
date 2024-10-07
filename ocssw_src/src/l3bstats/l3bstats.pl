#!/usr/bin/perl -w

# This script reads a level-3 daily bin file, computes some global statistics
# from the data in that file, and writes the statistics out as a series of
# binary files -- one output file per data set in the input file.  The
# twenty-one four-byte fields in each output file are as follows.

#  1.  (unsigned integer) number of bins containing data on this day
#  2.  ( floating-point ) minimum of the means
#  3.  ( floating-point ) maximum of the means
#  4.  ( floating-point ) mean of the means
#  5.  ( floating-point ) median of the means
#  6.  ( floating-point ) standard deviation of the means
#  7.  ( floating-point ) minimum of the standard deviations
#  8.  ( floating-point ) maximum of the standard deviations
#  9.  ( floating-point ) mean of the standard deviations
# 10.  ( floating-point ) median of the standard deviations
# 11.  ( floating-point ) standard deviation of the standard deviations
# 12.  (unsigned integer) minimum of the number of observations
# 13.  (unsigned integer) maximum of the number of observations
# 14.  ( floating-point ) mean of the number of observations
# 15.  ( floating-point ) median of the number of observations
# 16.  ( floating-point ) standard deviation of the number of observations
# 17.  (unsigned integer) minimum of the number of scenes
# 18.  (unsigned integer) maximum of the number of scenes
# 19.  ( floating-point ) mean of the number of scenes
# 20.  ( floating-point ) median of the number of scenes
# 21.  ( floating-point ) standard deviation of the number of scenes
#
# All values are written in littleendian order.
#
# The unit string given in the bin file is appended to the above 21 values.

# Norman Kuring		   Feb-2015	Initial development
# Norman Kuring		11-May-2015	Add standard deviation of each
#					parameter to the output.
# Norman Kuring		13-Feb-2017	Retrieve units strings and place in
#					output files.
# Norman Kuring		 8-Jun-2017	Ignore qual_l3 dataset in SST files and
#					the vestigial qualityDim dataset.
# Norman Kuring         25-Oct-2017     Modified for use by John & Tommy's
#                                       data system.

use strict;
use POSIX;
use vars qw($opt_v $opt_x);
use Getopt::Std;
my $optok = getopts('vx:');

my $pwd = `pwd`;
chop $pwd;

my $infile = shift or die <<"--";
Usage:
$0 [-v] [-x h5dump_path] L3_bin_file
--
my $basename = $infile;
$basename =~ s/.*\///;
$basename =~ s/\.nc$//;
$infile =~ /^\// or substr($infile,0,0) = "$pwd/";

my $h5dump = $opt_x || 'h5dump';

# Extract the BinList dataset.
doit("
$h5dump -d /level-3_binned_data/BinList -b -o BinList $infile >/dev/null
");

# Read the bin numbers, weights, etc. into arrays.
# I am assuming that BinList records are structured as follows.
#     DATATYPE "binListType" H5T_COMPOUND {
#        H5T_STD_U32LE "bin_num";
#        H5T_STD_I16LE "nobs";
#        H5T_STD_I16LE "nscenes";
#        H5T_IEEE_F32LE "weights";
#        H5T_IEEE_F32LE "time_rec";
#     }
my (@bin_num,@nobs,@nscenes,@weights);
open F,'<','BinList' or die "Could not open BinList ($!)";
my $rec;
while(read F,$rec,16){
  my ($bin_num,$nobs,$nscenes,$weights) = unpack "L<s<s<f",$rec;
  push @bin_num,$bin_num;
  push @nobs   ,$nobs;
  push @nscenes,$nscenes;
  push @weights,$weights;
}
close F;

# Get a list of the geophysical datasets in the input file.
# I expect these datasets to look like this.
#     DATATYPE "binDataType" H5T_COMPOUND {
#        H5T_IEEE_F32LE "sum";
#        H5T_IEEE_F32LE "sum_squared";
#     }

my @datasets;
doit("
$h5dump -n $infile | grep dataset | grep -v '/[Bb]in' | grep -v 'qual_l3' | grep -v 'qualityDim'> datasets
");
open F,'<','datasets' or die "Could not open list of datasets ($!)";
while(<F>){
  my $ds = (split ' ')[-1];
  push @datasets,$ds;
}
close F;

# Get the unit strings for the products from the 'units' attribute.
$_ = `$h5dump -a units $infile`;
my ($units) = /DATA\s*\{\s*\(0\):\s*"([^"]+)/;
my %units = map {split /:/,$_,2} split /,/,$units;

# Unpack the datasets.
my (%sum,%sumsq);
foreach my $ds (@datasets){
  doit("
  $h5dump -d $ds -b -o L3data $infile >/dev/null
  ");

  open F,'<','L3data' or die "Could not open L3data ($!)";
  my $rec;
  while(read F,$rec,8){
    my ($sum,$sumsq) = unpack "ff",$rec;
    push @{$sum{  $ds}},$sum;
    push @{$sumsq{$ds}},$sumsq;
  }
  close F;
}

my $numbins = @bin_num;

foreach my $ds (@datasets){

  my $ds_short = (split /\//, $ds)[-1];
  my $outfile = "$pwd/$basename.$ds_short.cbs"; # cbs=cumulative bin statistics

  # Gather up the means and standard deviations for all the filled bins.
  my (@means,@stdevs);
  for(my $i=0; $i<$numbins; $i++){

    # This should not happen, but I'm checking anyway in case
    # of corrupted data.
    if($nobs[$i] < 1 or $weights[$i] <= 0){
      warn "
      bin_num: $bin_num[$i]
      nobs:    $nobs[$i]
      nscenes: $nscenes[$i]
      weight:  $weights[$i]
      ";
      next;
    }

    my $mean = $sum{$ds}[$i]/$weights[$i];
    push @means, $mean;
    my $variance = $sumsq{$ds}[$i]/$weights[$i] - $mean*$mean;

    # Rounding errors in the 32-bit floats sometimes cause
    # the variance to go negative, so I fix that here.
    $variance = 0 if $variance < 0;

    push @stdevs, sqrt($variance);
  }

  ####################################################################
  # Compute minima, maxima, means, standard deviations, and medians  #
  # of the mean, standard deviation, number of observations, and     #
  # number of scenes per bin.                                        #
  ####################################################################
  my ( $mmean,$mmedian,$mstd,$mmin,$mmax,
       $smean,$smedian,$sstd,$smin,$smax,
       $nmean,$nmedian,$nstd,$nmin,$nmax,
       $omean,$omedian,$ostd,$omin,$omax ) = (0) x 12;
  # First, sort the data.
  @means  = sort {$a<=>$b}  @means;
  @stdevs = sort {$a<=>$b}  @stdevs;
  @nobs   = sort {$a<=>$b}  @nobs;
  @nscenes= sort {$a<=>$b}  @nscenes;
  # The first and last elements of the sorted data are the extrema.
  ($mmin,$mmax) =   @means[0,-1];
  ($smin,$smax) =  @stdevs[0,-1];
  ($nmin,$nmax) =    @nobs[0,-1];
  ($omin,$omax) = @nscenes[0,-1];
  # Sum up all the data...
  grep { $mmean += $_; 0 }   @means;
  grep { $smean += $_; 0 }  @stdevs;
  grep { $nmean += $_; 0 }    @nobs;
  grep { $omean += $_; 0 } @nscenes;
  # ... and divide by the number of measurements to get the means.
  $mmean /= @means;
  $smean /= @stdevs;
  $nmean /= @nobs;
  $omean /= @nscenes;
  # Sum up the squared deviations from the mean, ...
  grep { $mstd += ($_ - $mmean)*($_ - $mmean); 0 }   @means;
  grep { $sstd += ($_ - $smean)*($_ - $smean); 0 }   @stdevs;
  grep { $nstd += ($_ - $nmean)*($_ - $nmean); 0 }   @nobs;
  grep { $ostd += ($_ - $omean)*($_ - $omean); 0 }   @nscenes;
  # ... divide by the number of measurements, and take the square root
  # to get the standard deviation.
  $mstd = sqrt($mstd/@means  );
  $sstd = sqrt($sstd/@stdevs );
  $nstd = sqrt($nstd/@nobs   );
  $ostd = sqrt($ostd/@nscenes);
  # The middle value of the sorted data is the median.
  my $mid = @means/2;
  $mmedian = @means%2   ? $means[int($mid)]
                        : ($means[$mid] + $means[$mid-1])/2;
  $smedian = @stdevs%2  ? $stdevs[int($mid)] 
                        : ($stdevs[$mid] + $stdevs[$mid-1])/2;
  $nmedian = @nobs%2    ? $nobs[int($mid)]
                        : ($nobs[$mid] + $nobs[$mid-1])/2;
  $omedian = @nscenes%2 ? $nscenes[int($mid)]
                        : ($nscenes[$mid] + $nscenes[$mid-1])/2;

  my $output = pack "LffffffffffLLfffLLfff",
  $numbins,
  $mmin,$mmax,$mmean,$mmedian,$mstd,
  $smin,$smax,$smean,$smedian,$sstd,
  $nmin,$nmax,$nmean,$nmedian,$nstd,
  $omin,$omax,$omean,$omedian,$ostd;

  open F,'>',$outfile or die "Could not create $outfile ($!)";
  print F $output,$units{$ds_short};
  close F;
}



sub doit{
  my $cmd = join ' ', split ' ', shift;
  $opt_v and print STDERR "[",scalar localtime,"] $cmd\n";
  system($cmd) == 0 or die "Error executing $cmd ($!)";
}
