# l3mapgen Changelog
This program reads L3 files and produces mapped output in a variety of projections
It is the replacement for the "smigen" program


## 2.3.0 - 2019-11-25

 - output 1D arrays for lat/lon for PlatteCarree instead of full 2D arrays


## 2.2.0 - 2019-07-29

 - generalized (using regex) the EPSG GTiff tag set-up for l3mapgen
 -  updated GIBS polar projections to a fixed EPSG:3413 and 3031 for north and south pole respectively
 - fixed a l3mapgen bug when a GIBS image crosses the date line
 - added transparency to l3mapgen tiff and png output.
 - make sure width and heaight ar at least 1 pixel
 - added "gibs" and "alaska" projection shortcuts; removed alt_projection and alt_thresh_lat options - unnecessary with "gibs" projection
 - can now output an image file for each product
 - fixed geospatial_lat_resolution and geospatial_lon_resolution to be in units of degrees.  Get the default resolution from the input file.
 - fixed l3mapgen deletes
 - fixed a missmatched new/free in l3mapgen.  Changed to delete.
 - const => constexpr
 - set l3mapgen output to fill value if proj4 returns a NAN
 - fixed a conversion problem in l3mapgen when the lat/lon is INF
 - fixed compression/chunking parameters
 - turned off lat/lon arrays for raw projection
 - add option to NOT write latitude, longitude arrays to NetCDF output
 -  write full latitude, longitude arrays to NetCDF output
 -  hdf4 output also uses new method OutFile::ProductStuff::calcOutputLineVals.
 -   NetCDF output is handled entirely with the c++ API.  QualityData (qual_sst) had been written as byte with _FillValue = -1b; is now written as ubyte with _FillValue = 255UB.
 -  add width parameter
 - fixed raw output for cm < 0.

### Source Changed
 * CMakeLists.txt
 * l3mapgen.h
 * OutFile.cpp
 * OutFile.h
 * l3mapgen.cpp
 * l3mapgen_input.cpp

## 2.1.0 - 2018-04-09
 - fixed l3mapgen bug when writing netCDF output file with no pixels
 - search share/common for suite defaults files first.
 - added subsensor search for default files in l3mapgen
 - l3mapgen exit status=110 when no mappable pixels.
 - made l3mapgen overwrite output file
 - use NetCDF-C++ API for most l3mapgen metadata.
 - use new sensor info functions

### Source Changed
  * CMakeLists.txt
  * l3mapgen_input.cpp
  * l3mapgen.cpp
  * l3mapgen64.cpp
  * OutFile.cpp
  * OutFile64.cpp
  * OutFile.h
  * OutFile64.h

## 2.0.0 - 2017-11-02
 - normalized GCMD keywords and DOI code
 - added suites
 - can do quality processing with multiple output products
 - write projection parameters when creating GeoTIFF files
 - added alternate projection if file above a threashold latitude
 - only use the fudge factor if pixel would have been empty
 - added land mask option

## 1.0.1 - 2015-04-28
fixed the filled pixels percentage output

### Source Changed
  * OutFile.cpp
  * l3mapgen.cpp

## 1.0 - 2015-03-10
- First version of the program.

### Added
  * Changelog.md

