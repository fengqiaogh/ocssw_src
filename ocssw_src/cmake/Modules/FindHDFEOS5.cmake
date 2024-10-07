# HDFEOS5_FOUND - true if library and headers were found
# HDFEOS5_INCLUDE_DIRS - include directories
# HDFEOS5_LIBRARIES - library directories

#set_ocssw_policy()

find_path(HDFEOS5_INCLUDE_DIR HE5_HdfEosDef.h
  HINTS $ENV{LIB3_DIR}/EOS/include
  )

find_library(HDFEOS5_LIBRARY NAMES he5_hdfeos
  HINTS $ENV{LIB3_DIR}/EOS/lib/${EOS_ARCH}
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HDFEOS5 DEFAULT_MSG HDFEOS5_LIBRARY HDFEOS5_INCLUDE_DIR)
if(HDFEOS5_FOUND)
  set(HDFEOS5_LIBRARIES ${HDFEOS5_LIBRARY})
  set(HDFEOS5_INCLUDE_DIRS ${HDFEOS5_INCLUDE_DIR})

  find_package(HDF5 COMPONENTS C)
  if(HDF5_FOUND)
    list(APPEND HDFEOS5_LIBRARIES ${HDF5_LIBRARIES})
    list(APPEND HDFEOS5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
  else()
    set(HDFEOS5_FOUND FALSE)
    unset(HDFEOS5_LIBRARIES)
    unset(HDFEOS5_INCLUDE_DIRS)
  endif()
  
endif()

mark_as_advanced(HDFEOS5_INCLUDE_DIR HDFEOS5_LIBRARY)
