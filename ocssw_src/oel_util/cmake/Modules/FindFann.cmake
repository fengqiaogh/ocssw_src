# Fann_FOUND - true if library and headers were found
# Fann_INCLUDE_DIRS - include directories
# Fann_LIBRARIES - library directories

find_path(Fann_INCLUDE_DIR fann.h
  HINTS $ENV{LIB3_DIR}/include
  )

find_library(Fann_LIBRARY NAMES fann
  HINTS $ENV{LIB3_DIR}/lib
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Fann DEFAULT_MSG Fann_LIBRARY Fann_INCLUDE_DIR)
if(Fann_FOUND)
  set(Fann_LIBRARIES ${Fann_LIBRARY})
  set(Fann_INCLUDE_DIRS ${Fann_INCLUDE_DIR})
endif()
  
mark_as_advanced(Fann_INCLUDE_DIR Fann_LIBRARY)
