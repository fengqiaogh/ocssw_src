# HDFEOS_FOUND - true if library and headers were found
# HDFEOS_INCLUDE_DIRS - include directories
# HDFEOS_LIBRARIES - library directories

find_path(HDFEOS_INCLUDE_DIR HdfEosDef.h
  HINTS $ENV{LIB3_DIR}/EOS/include
  )

find_library(HDFEOS_LIBRARY NAMES hdfeos
  HINTS $ENV{LIB3_DIR}/EOS/lib/${EOS_ARCH}
  )

find_library(GCTP_LIBRARY NAMES Gctp
  HINTS $ENV{LIB3_DIR}/EOS/lib/${EOS_ARCH}
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HDFEOS DEFAULT_MSG HDFEOS_LIBRARY GCTP_LIBRARY HDFEOS_INCLUDE_DIR)
if(HDFEOS_FOUND)
  set(HDFEOS_LIBRARIES ${HDFEOS_LIBRARY} ${GCTP_LIBRARY})
  set(HDFEOS_INCLUDE_DIRS ${HDFEOS_INCLUDE_DIR})

  find_package(HDF4)
  if(HDF4_FOUND)
    list(APPEND HDFEOS_LIBRARIES ${HDF4_LIBRARIES})
    list(APPEND HDFEOS_INCLUDE_DIRS ${HDF4_INCLUDE_DIRS})
  else()
    set(HDFEOS_FOUND FALSE)
    unset(HDFEOS_LIBRARIES)
    unset(HDFEOS_INCLUDE_DIRS)
  endif()
  
endif()

mark_as_advanced(HDFEOS_INCLUDE_DIR HDFEOS_LIBRARY)
