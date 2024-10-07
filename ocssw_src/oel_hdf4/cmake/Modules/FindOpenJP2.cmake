# OpenJP2_FOUND - true if library and headers were found
# OpenJP2_INCLUDE_DIRS - include directories
# OpenJP2_LIBRARIES - library directories

find_path(OpenJP2_INCLUDE_DIR openjpeg.h
  HINTS $ENV{LIB3_DIR}/include
  )

find_library(OpenJP2_LIBRARY NAMES openjp2 opj_compress opj_decompress opj_dump
  HINTS $ENV{LIB3_DIR}/lib
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenJP2 DEFAULT_MSG OpenJP2_LIBRARY OpenJP2_INCLUDE_DIR)
if(OpenJP2_FOUND)
  set(OpenJP2_LIBRARIES ${OpenJP2_LIBRARY})
  set(OpenJP2_INCLUDE_DIRS ${OpenJP2_INCLUDE_DIR})
endif()
  
mark_as_advanced(OpenJP2_INCLUDE_DIR OpenJP2_LIBRARY)
