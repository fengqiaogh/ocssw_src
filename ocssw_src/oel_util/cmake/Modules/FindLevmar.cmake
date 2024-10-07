# Levmar_FOUND - true if library and headers were found
# Levmar_INCLUDE_DIRS - include directories
# Levmar_LIBRARIES - library directories

find_path(Levmar_INCLUDE_DIR levmar.h
  HINTS $ENV{LIB3_DIR}/include
  )

find_library(Levmar_LIBRARY NAMES levmar
  HINTS $ENV{LIB3_DIR}/lib
  )

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Levmar DEFAULT_MSG Levmar_LIBRARY Levmar_INCLUDE_DIR)
if(Levmar_FOUND)
  set(Levmar_LIBRARIES ${Levmar_LIBRARY})
  set(Levmar_INCLUDE_DIRS ${Levmar_INCLUDE_DIR})
endif()

mark_as_advanced(Levmar_INCLUDE_DIR Levmar_LIBRARY)
