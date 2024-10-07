# SDST_FOUND - true if library and headers were found
# SDST_INCLUDE_DIRS - include directories
# SDST_LIBRARIES - library directories

find_path(SDST_INCLUDE_DIR SDST_TK.h
  HINTS $ENV{LIB3_DIR}/include
  )

find_library(SDST_LIBRARY NAMES sdst
  HINTS $ENV{LIB3_DIR}/lib
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SDST DEFAULT_MSG SDST_LIBRARY SDST_INCLUDE_DIR)
if(SDST_FOUND)
  set(SDST_LIBRARIES ${SDST_LIBRARY})
  set(SDST_INCLUDE_DIRS ${SDST_INCLUDE_DIR})
endif()

find_package(PGSTK REQUIRED)
list(APPEND SDST_LIBRARIES ${PGSTK_LIBRARIES})
list(APPEND SDST_INCLUDE_DIRS ${PGSTK_INCLUDE_DIRS})

mark_as_advanced(SDST_INCLUDE_DIR SDST_LIBRARY)
