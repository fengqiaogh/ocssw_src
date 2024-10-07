
# libnetcdfutils Changelog

## <VERSION STILL IN LIMBO> - 2018-04-10

 - target_include_directories, again.
 - squashed another memory leak
 - implemented OBPG standard code formatting.
 - get everything ready for shared libs
 - increased the chunk cache size to speed up full CONUS 300m mapping
 - added chunck cache for nc_griduils
 - Made nc_gridutils safe for c++
 - changed the NetCDF compression int function.
 - got rid of some warnings

### Source changed
  * cdl_utils.cpp
  * cdl_utils.h
  * CMakeLists.txt
  * createNCDF.c
  * nc4utils.c
  * nc4utils.h
  * nc_gridutils.c
  * nc_gridutils.h
  * nc_init_compress.c

## <VERSION STILL IN LIMBO> - 2015-05-08

Pulled the NetCDF functions into its own library. 

## <VERSION STILL IN LIMBO> - 2015-04-11
### Source Changed
  * ncdf_utils.c - added chunking and deflate options to writeQuality_nc function

## <VERSION> - 2015-03-12
### Added
  * Changelog.md
