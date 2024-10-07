
# libpiutils Changelog

## <VERSION STILL IN LIMBO> - 2018-04-10

 - added sentinel 2B support
 - fixed memory leak
 - fixed the seawifs gac/lac dir name
 - reorg of subsensor defaults files
 - new setup for VIIRS and MODIS subsensors
 - simplified getDOI and added getGCMDKeywords
 - implemented OBPG standard code formatting.
 - get ready for shared libs
 - removed roxml
 - productInfo now uses pugi xml instead of roxml
 - added titleFormat to product structure
 - switched from roxml to pugixml for XML parsing of the DOI file.
 - added Sentinel-2A MSI, L5TM and L7ETM sensor definitions
 - got rid of dependency on boost.
 - added sensor SGLI

### Source Changed

  * Changelog.md
  * CMakeLists.txt
  * getDOI.cpp
  * getGCMDKeywords.cpp
  * piutils.h
  * productInfo.c
  * productInfo.cpp
  * productInfo.h
  * rdsensorinfo.c
  * sensorDefs.h
  * sensorInfo.c
  * sensorInfo.h

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
