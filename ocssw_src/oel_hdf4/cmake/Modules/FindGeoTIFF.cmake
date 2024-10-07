# GeoTIFF_FOUND - true if library and headers were found
# GeoTIFF_INCLUDE_DIRS - include directories
# GeoTIFF_LIBRARIES - library directories

find_path(GeoTIFF_INCLUDE_DIR geotiff.h
  HINTS $ENV{LIB3_DIR}/include
  PATH_SUFFIXES geotiff
  )

find_library(GeoTIFF_LIBRARY NAMES geotiff
  HINTS $ENV{LIB3_DIR}/lib
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GeoTIFF DEFAULT_MSG GeoTIFF_LIBRARY GeoTIFF_INCLUDE_DIR)
if(GeoTIFF_FOUND)
  set(GeoTIFF_LIBRARIES ${GeoTIFF_LIBRARY})
  set(GeoTIFF_INCLUDE_DIRS ${GeoTIFF_INCLUDE_DIR})

  find_package(TIFF)
  if(TIFF_FOUND)
    list(APPEND GeoTIFF_LIBRARIES ${TIFF_LIBRARIES})
    list(APPEND GeoTIFF_INCLUDE_DIRS ${TIFF_INCLUDE_DIRS})
  else()
    set(GeoTIFF_FOUND FALSE)
    unset(GeoTIFF_LIBRARIES)
    unset(GeoTIFF_INCLUDE_DIRS)
  endif()
  
endif()

mark_as_advanced(GeoTIFF_INCLUDE_DIR GeoTIFF_LIBRARY)
