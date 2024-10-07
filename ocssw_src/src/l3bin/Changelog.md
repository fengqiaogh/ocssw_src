
# l3bin Changelog

## [Unreleased]

## 5.13 2018-10-31
 - fixed core dump if no valid pixels found to bin.

## 5.12 also 2018-09-01
 - forgot to bump version number
 - combined 32 and 64 bit version

## 5.12 2018-04-09
 - use new sensor info functions
 - implemented OBPG standard code formatting
 - fixed bug where core dumped when displaying help
 - added infile as an alias for ifile to make l3bin backward compatable
 - fixed buffer over run in l3bin
 - fixed compiler warnings and cleaned up code
 - Fixed prototype error

### Source Changed
  * Changelog.md
  * CMakeLists.txt
  * l3bin64.cpp
  * l3bin.cpp
  * l3bin_input.c
  * l3bin_input.h
  * l3binmerge64.cpp
  * l3binmerge.cpp
  * l3binmerge.mk
  * l3bin.mk
  * l3rebin_meris64.cpp
  * l3rebin_meris.cpp
  * l3rebin_meris.mk

## 5.11 2017-11-02
normalized GCMD keywords and DOI code

## 4.32 2016-10-31
Fix bug in hdf5_bin::create (Aquarius): Set n_data_prod to 0.

## 4.31 2016-08-06
Add composite_scheme parameter for land composite binning.

## 4.30 2016-07-12
Add support for land composite binning

## 4.20 2015-09-03
Add support for time_rec

## 4.10 2015-08-03
Add support for L2C SMAP binfile input

## 4.03 2015-04-28
In the case of reduce_fac != 1 for quality field allow
  "i" and "n_write" to differ by 2 or less due to adjacent
  map rows of integralized sinesodial projection being of different lengths. 

## 4.02 - 2015-03-12
### Added
  * Changelog.md
