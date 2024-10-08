
# libl2 Changelog

## <VERSION STILL IN LIMBO> - 2023-11-24

- scan line variables are now optional
- searches for latitude and longitude in a netCDF file in any NC group
- support for 3D product expansion has been added.
- added utilities for calculation of geospatial bounds if they are missing
- added utilities for attribute and variable search in a netCDF file
- reading l2 flags is optional
- l2flags values are read from the ncFile
- 
### Source Changed
* Changelog.md
* CMakeLists.txt
* readL2scan.c
* readL2scan.h
* setupflags.c
* setupflags.h

### New files added
* get_geospatial.cpp
* get_geospatial.hpp
* expand3D.cpp
* expand3D.hpp
* find_variable.hpp
* find_variable.h
* find_variable.c

## <VERSION STILL IN LIMBO> - 2018-04-09

 - removed "check and fail" for readL2scan when looking for orbit number
 - use new sensor info functions
 - fixed a bug in readL2 that was using a NULL after reading the units metadata.
 - fixed reading of netCDF string attributes fixed the generic L1B reader
 - implemented OBPG standard code formatting
 - get ready for shared libs
 - moved lots of memory from stack to heap
 - fixed reading sensor name metadata inro a null pointer.

### Source Changed
  * Changelog.md
  * CMakeLists.txt
  * get_product_table.h
  * readL2scan.c
  * readL2scan.h
  * setupflags.c
  * setupflags.h

## <VERSION STILL IN LIMBO> - 2015-05-07
moved the productInfo routines to the new piutils library

## <VERSION STILL IN LIMBO> - 2015-04-28
added a function to assemble the full product name from a ProductInfo structure.

### Source Changed
modified code
  * include/productInfo.h
  * productInfo.c

## <VERSION> - 2015-03-12
### Added
  * Changelog.md
