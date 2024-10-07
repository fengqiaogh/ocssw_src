# Cminpack_FOUND - true if library and headers were found
# Cminpack_INCLUDE_DIRS - include directories
# Cminpack_LIBRARIES - library directories

find_path(cminpack_INCLUDE_DIR cminpack.h
  PATHS $ENV{LIB3_DIR}/include
  )

find_library(cminpack_LIBRARY NAMES cminpack
  PATHS $ENV{LIB3_DIR}/lib
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(cminpack DEFAULT_MSG cminpack_LIBRARY cminpack_INCLUDE_DIR)
if(cminpack_FOUND)
  set(cminpack_LIBRARIES ${cminpack_LIBRARY})
  set(cminpack_INCLUDE_DIRS ${cminpack_INCLUDE_DIR})
endif()
  
mark_as_advanced(cminpack_INCLUDE_DIR cminpack_LIBRARY)
