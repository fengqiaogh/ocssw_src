# l2bin Changelog

This program reads L2 files output from the l2gen program and produces a L3 file
Changes for v4.0 and prior are reported in the code as comments...maybe...

## 6.2.0 - 2020-01-07
 - added default equater crossing time to 12
 - deleted old command line options related to oformat (only make netCDF now)

### Sources changed

l2bin.cpp
l2bin_input.cpp
l2bin_input.h


## 6.1.0 - 2019-11-15

 - added area weighting.

## 6.0.0 - 2019-06-18

 - changed to C++
 - now using the c++ API for netCDF
 - can no longer write HDF4 bin files

## 5.0.0 - 2019-02-18

 - fixed the resolution help message in l2bin
 - cleaned up a few things in l2bin like unused function forward declaration
 - deleted the old bin2ll routines
 - added sentinel 3B and made subsensors for OLCI

### Source Changed
 * l2bin_input.c
 * l2bin.c

## 4.5.0 - 2018-04-09
 - fixed l2bin subsampling lat/lon arrays and fixed 0 bins filled exit status
 - only keep 1 L2 file open at a time to fix the HDF5 compression leak
 - fixed order in which suite defaults files are read
 - use new sensor info functions
 - normalized how DOIs and GCMD Keywords are looked up.
 - added check that a pixel exists within the rowgroup - spurious bad inputs were causing segfaults...
 - implemented OBPG standard code formatting
 - fixed binner code for failing bounding box
 - fixed another l2bin east/west metadata bug
 - fixed case when only 1 pixel matched conditions for geospatial min/max metadata
 - removed roxml form l2bin
 - fixed a few memory over runs in l2bin using valgrind
 - changed compression init for l2bin
 - fixed compiler warnings and cleanup code
 - changed latbin to a double to increase precision.
 - Changed l3bprod to use validMin/Max from product.xml instead of displayMin/Max
 - Changed displayMin/Max to validMin/Max for l3bprod Min/Max Added clo (command line option) library for l2bin to use.
 - Adding clo to l2bin
 - Updates to add the min from the product.xml file as the default

### Source Changed
  * cdata.h
  * Changelog.md
  * CMakeLists.txt
  * dataday.c
  * get_l3bprod_index.c
  * l2bin64.c
  * l2bin.c
  * l2bin_input.c
  * l2bin_input.h
  * l2bin.mk
  * l2prod.h
  * l3bprod_struc.h
  * swfnav.h

## 4.4.2 - 2018-09-04
  * Removed support for HEALPIX

## 4.4.1 - 2017-11-02
  * bug fixes
  * normalized GCMD keywords and DOI code

## 4.4.0 2016-11-29
  * Added clo library to replace command line processing functions

## 4.3.0 2016-11-28
  * Added in xml code to read product min/max value from product.xml file.  Changed version to match l2bin64

## 4.1.2 2016-08-04
  * Add composite_scheme parameter for land composite binning.

## 4.1.1 2016-07-02
  * Add support for land composite binning.

## 4.1.0 2016-06-16
  * Allocate and clear dolat/dolon arrays in dataday code to avoid time bombs.

## 4.0.9 2015-01-06
  * Fixed issue introduced with 4.0.5 with default sday/eday which resulted 
    in a change in behavior
  * NOTE: if this program is still in use in 2038...may need a tweak :)
    Modified test for "regional" prodtype to be case insensitive
 
## 4.0.8 2015-11-30
  * Fixed opening group id for multiple files.
 
## 4.0.7 2015-10-26
  * Fix reading of year,day,time values for time_rec

## 4.0.6 2015-10-26
  * Accumulate time_rec sums at double precision, store at float
    nature of crossing the pole. 

## 4.0.5 2015-10-01
  * Added logic to handle orbit based files that cross the dateline by 
    nature of crossing the pole. 

## 4.0.4 2015-09-15
  * Add support for time_rec field in the binlist records

## 4.0.3 - 2015-04-11
### Source Changed
modified code
  * dataday.c - added more logic to ensure more than one valid scan of data exists, otherwise return a status 1
  * l2bin.c - added logic to output a qcfail file when fileuse is set

## 4.0.2 - 2015-04-03
### Source Changed
modified code
  * dataday.c - added more logic to handle bad navigation, modified daynight_outlines to return a status
  * dataday.h - daynight_outlines prototype update
  * l2bin.c - added logic to capture status from daynight_outlines, and better handle small scenes; modified memory allocation for dorn array (day or night) 

## 4.0.1 - 2015-04-01
### Source Changed
modified code
  * dataday.c - added logic to handle bad navigation

## 4.0 - 2015-03-10
### Added
  * Changelog.md

