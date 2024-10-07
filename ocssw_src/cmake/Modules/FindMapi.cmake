# Mapi_FOUND - true if library and headers were found
# Mapi_INCLUDE_DIRS - include directories
# Mapi_LIBRARIES - library directories

find_path(Mapi_INCLUDE_DIR mapi.h
  HINTS $ENV{LIB3_DIR}/include
  )

find_library(Mapi_LIBRARY NAMES mapi
  HINTS $ENV{LIB3_DIR}/lib
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Mapi DEFAULT_MSG Mapi_LIBRARY Mapi_INCLUDE_DIR)
if(Mapi_FOUND)
  set(Mapi_LIBRARIES ${Mapi_LIBRARY})
  set(Mapi_INCLUDE_DIRS ${Mapi_INCLUDE_DIR})
endif()

find_package(PGSTK REQUIRED)
list(APPEND Mapi_LIBRARIES ${PGSTK_LIBRARIES})
list(APPEND Mapi_INCLUDE_DIRS ${PGSTK_INCLUDE_DIRS})

mark_as_advanced(Mapi_INCLUDE_DIR Mapi_LIBRARY)
