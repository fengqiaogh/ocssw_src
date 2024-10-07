# - Find HDF4
# Find the native HDF4 includes and library
#
# Set HDF4_DIR to give a hint of where HDF4 is installed
#
#  HDF4_INCLUDE_DIR  - user modifiable choice of where netcdf headers are
#  HDF4_LIBRARY      - user modifiable choice of where netcdf libraries are
#
# This module returns these variables for the rest of the project to use.
#
#  HDF4_FOUND          - True if NetCDF found including required interfaces (see below)
#  HDF4_LIBRARIES      - All netcdf related libraries.
#  HDF4_INCLUDE_DIRS   - All directories to include.

#search starting from user editable cache var
if (HDF4_INCLUDE_DIR AND HDF4_LIBRARY)
  # Already in cache, be silent
  set (HDF4_FIND_QUIETLY TRUE)
endif ()

find_path (HDF4_INCLUDE_DIR hdf.h
  HINTS ${HDF4_DIR}/include $ENV{LIB3_DIR}/include)
mark_as_advanced (HDF4_INCLUDE_DIR)

find_library (HDF4_LIBRARY NAMES mfhdf
  HINTS ${HDF4_DIR}/lib $ENV{LIB3_DIR}/lib)
mark_as_advanced (HDF4_LIBRARY)

find_library (DF_LIBRARY NAMES hdf
  HINTS ${HDF4_DIR}/lib $ENV{LIB3_DIR}/lib)
mark_as_advanced (DF_LIBRARY)

# handle the QUIETLY and REQUIRED arguments and set HDF4_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (HDF4 DEFAULT_MSG HDF4_LIBRARY HDF4_INCLUDE_DIR)
if(HDF4_FOUND)
  set(HDF4_LIBRARIES ${HDF4_LIBRARY} ${DF_LIBRARY})
  set(HDF4_INCLUDE_DIRS ${HDF4_INCLUDE_DIR})

  find_package(JPEG)
  if(JPEG_FOUND)
    list(APPEND HDF4_LIBRARIES ${JPEG_LIBRARIES})
    list(APPEND HDF4_INCLUDE_DIRS ${JPEG_INCLUDE_DIR})

    find_package(ZLIB)
    if(ZLIB_FOUND)
      list(APPEND HDF4_LIBRARIES ${ZLIB_LIBRARIES})
      list(APPEND HDF4_INCLUDE_DIRS ${ZLIB_INCLUDE_DIR})
    else()
      set(HDF4_FOUND FALSE)
      unset(HDF4_LIBRARIES)
      unset(HDF4_INCLUDE_DIRS)
    endif()
    
  else()
    set(HDF4_FOUND FALSE)
    unset(HDF4_LIBRARIES)
    unset(HDF4_INCLUDE_DIRS)
  endif()
  
endif()

