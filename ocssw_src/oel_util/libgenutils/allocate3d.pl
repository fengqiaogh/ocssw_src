#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';

use FindBin;
use Data::Dumper;

my @types = (
	#['char', '_char'],
	#['unsigned char', '_uchar'],
	#['signed char', '_schar'],
	['short', '_short'],
	['int', '_int'],
	#['long', '_long'],
	['float', '_float'],
	['double', '_double'],
	#['unsigned short', '_ushort'],
	#['unsigned int', '_uint'],
	#['long long', '_longlong'],
	#['unsigned long long', '_ulonglong'],
);

sub generate_allocate3d {
	#my @args = @_;
	
	my %templates = read_templates();

	my $c_file = "$FindBin::RealDir/allocate3d.c";
	my $h_file = "$FindBin::RealDir/allocate3d.h";

	open(my $c_fh, '>', $c_file) || die("Couldn't open $c_file: $!");
	open(my $h_fh, '>', $h_file) || die("Couldn't open $h_file: $!");
	print {$c_fh} $templates{c_header};
	print {$h_fh} $templates{h_header};
	for (@types){
		my ($type, $suffix) = @$_;
		print {$c_fh} replace_variables($templates{c_functions}, type => $type, suffix => $suffix);
		print {$h_fh} replace_variables($templates{h_functions}, type => $type, suffix => $suffix);
	}
	print {$h_fh} $templates{h_footer};
	close($h_fh);
	close($c_fh);

	return;
}

sub replace_variables {
	my ($string, %variables) = @_;
	while (my ($variable, $value) = each(%variables)){
		$string =~ s/\$\Q$variable\E\$/$value/mg;
	}
	return $string;
}

sub read_templates {
	my %ret;
	my @data = <DATA>;
	my ($filename, $buffer);
	while (my $line = shift(@data)){
		if ($line =~ /^#####\s*(.*?)$/){
			if ($filename){
				$ret{$filename} = $buffer;
			}
			$filename = $1;
			$buffer = '';
		} else {
			$buffer .= $line;
		}
	}
	if ($filename){
		$ret{$filename} = $buffer;
	}
	return %ret;
}


unless (caller){
	exit(generate_allocate3d(@ARGV) || 0);
}

__DATA__
##### c_header
/** @file allocate3d.c
	@brief Utility functions for allocating and freeing three-dimensional arrays of various types.

	This file was created by allocate3d.pl and should not be edited manually.
*/

#include "allocate3d.h"

#include <stdio.h>
#include <stdlib.h>


##### c_functions
$type$ ***allocate3d$suffix$(size_t nz, size_t ny, size_t nx){
    $type$ *x_ptr = ($type$*) malloc(nz * ny * nx * sizeof($type$));
    if (x_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of data block failed.n", __FILE__, __LINE__);
        return NULL;
    }
    $type$ **y_ptr = ($type$**) malloc(nz * ny * sizeof($type$*));
    if (y_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of y array failed.n", __FILE__, __LINE__);
        return NULL;
    }
    $type$ ***z_ptr = ($type$***) malloc(nz * sizeof($type$**));
    if (z_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of z array failed.n", __FILE__, __LINE__);
        return NULL;
    }
    for(size_t z=0; z<nz; z++) {
        for(size_t y=0; y<ny; y++) {
            y_ptr[y] = x_ptr;
            x_ptr += nx;
        }
        z_ptr[z] = y_ptr;
        y_ptr += ny;
    }
    return z_ptr;
}
void free3d$suffix$($type$ ***p) {
    free(p[0][0]);
    free(p[0]);
    free(p);
}

##### h_header
/** @file allocate3d.h
	@brief Utility functions for allocating and freeing three-dimensional arrays of various types.

	This file was created by allocate3d.pl and should not be edited manually.
*/
#ifndef OEL_UTIL_LIBGENUTILS_ALLOCATE3D_H_
#define OEL_UTIL_LIBGENUTILS_ALLOCATE3D_H_

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

##### h_functions
/** @brief Allocate a three-dimensional array of type $type$ of a given size.

	The order of the parameters are (now) the same as when the data is being accessed.
        E.g., this is valid (ignoring the printf type specifier):

	$type$ ***array = allocate3d$suffix$(7, 5, 1);
	printf("%d\n", array[6][4][0]);
	
	@param[in] nz slowest incrimenting dimension of array in memory.
	@param[in] ny middle dimension of array.
	@param[in] nx fastest incrimenting dimension of array.

	@return A malloc'd array or NULL if any malloc fails.
*/
$type$ ***allocate3d$suffix$(size_t nz, size_t ny, size_t nx);
/** @brief Free a three-dimensional array created by allocate3d$suffix$.

	@param[in] p Pointer to array created by allocate3d$suffix$.
*/
void free3d$suffix$($type$ ***p);

##### h_footer
#ifdef __cplusplus
}
#endif

#endif /* OEL_UTIL_LIBGENUTILS_ALLOCATE3D_H_ */
