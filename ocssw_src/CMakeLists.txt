cmake_minimum_required (VERSION 3.12)

# set the Mac OSX target.  Normally not needed, but this
# is the only way I got the upgraded Monterey Mac to work
set(CMAKE_OSX_SYSROOT macosx)
set(CMAKE_OSX_DEPLOYMENT_TARGET "11")

project (OCSSW)

# delete the anaconda dirs from the PATH.  That way find_library
# will use system libs instead of anaconda libs.  (ie. zlib)
set(GOOD_PATH "")
string(REPLACE ":" ";" PATH_LIST $ENV{PATH})
foreach(TMP_DIR ${PATH_LIST})
  if(NOT TMP_DIR MATCHES "conda")
    if(GOOD_PATH)
      set(GOOD_PATH "${GOOD_PATH}:${TMP_DIR}")
    else()
      set(GOOD_PATH "${TMP_DIR}")
    endif()
  endif()
endforeach(TMP_DIR)
set(ENV{PATH} ${GOOD_PATH})
unset(TMP_DIR)
unset(GOOD_PATH)

list(APPEND CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/oel_util/cmake/Modules
    ${CMAKE_CURRENT_SOURCE_DIR}/oel_hdf4/cmake/Modules
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

#include(SetOcsswPolicy)

set(BUILD_SHARED_LIBS ON)

# where to look for the external libraries
set(CMAKE_PREFIX_PATH $ENV{LIB3_DIR})

# set this to Release or Debug
# note: rerun cmake if you change the environment var OCSSW_DEBUG
if (NOT CMAKE_BUILD_TYPE)
    if ($ENV{OCSSW_DEBUG})
      set(CMAKE_BUILD_TYPE Debug)
    else($ENV{OCSSW_DEBUG})
      set(CMAKE_BUILD_TYPE Release)
    endif($ENV{OCSSW_DEBUG})
endif (NOT CMAKE_BUILD_TYPE)

# Finer-grained control...
OPTION(BUILD_HISTORICAL "Build the historical code (e.g. old seawifs binaries no longer needed)" OFF)
OPTION(BUILD_MISC "Build non-essential code" OFF)
OPTION(BUILD_ALL "Build all including  Aquarius, historical, misc and VIIRS L1 code" OFF)
OPTION(USE_MPI "Set definition to use MPI" OFF)
OPTION(BUILD_CLOUD "Build libcloud code" OFF)

if ($ENV{OCSSW_MPI})
    set(USE_MPI ON)
endif()

if(BUILD_ALL)
   set(BUILD_HISTORICAL ON)
   set(BUILD_MISC ON)
   set(BUILD_CLOUD ON)
endif(BUILD_ALL)

enable_language(C)
enable_language(CXX)
if(BUILD_HISTORICAL OR BUILD_CLOUD)
  enable_language(Fortran)
endif()

#include(PrintAllVariables)
#print_all_variables()

# if Doxygen exists, create the API documentation 
find_package(Doxygen)
if(DOXYGEN_FOUND AND IS_DIRECTORY ${OCSSW_SOURCE_DIR}/doc)
  configure_file(${OCSSW_SOURCE_DIR}/doc/Doxyfile
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY
  )
  add_custom_target(doc
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${OCSSW_SOURCE_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
endif()

# make the build
set(CMAKE_EXE_LINKER_FLAGS "-mcmodel=large")
set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set up C and CXX flags
set(CMAKE_C_FLAGS_RELEASE "-O3 -Wno-unused-result")
if (DEFINED ENV{GCC_TUNE})
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}  $ENV{GCC_TUNE}")
endif(DEFINED ENV{GCC_TUNE}) 
set(CMAKE_C_FLAGS_DEBUG "-ggdb -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
set(CMAKE_C_FLAGS "-DLINUX -Dlinux -Wno-parentheses")
set(CMAKE_CXX_FLAGS "-DLINUX -Dlinux -Wno-unused-local-typedefs -Wno-write-strings -Wno-parentheses")

set(CMAKE_EXPORT_COMPILE_COMMANDS on)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	if(CMAKE_C_COMPILER_VERSION VERSION_GREATER "6.0")
	  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-format-overflow")
	  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-format-overflow")
	endif()
endif()

# Set up FORTRAN flags
set(CMAKE_Fortran_FLAGS "-DLINUX -Dlinux -ffixed-line-length-none -fconvert=big-endian -fno-automatic -fno-range-check -cpp")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops")
set(CMAKE_Fortran_FLAGS_DEBUG "-ggdb -Wall -Wno-unused-variable -Wno-unused-label -Wno-unused-dummy-argument -Wno-conversion")

# set up OS-specific bits
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
  # standard Linux bits...
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_BSD_SOURCE -D_XOPEN_SOURCE=700 -D_DEFAULT_SOURCE")
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib;$ORIGIN/../opt/lib")
  set(EOS_ARCH linux)
  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
  # 64bit bits...
    set(EOS_ARCH linux64)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
    #set (CMAKE_SYSTEM_LIBRARY_PATH "/usr/lib/x86_64-linux-gnu")
    #set (CMAKE_SYSTEM_INCLUDE_PATH "/usr/include/x86_64-linux-gnu")
  endif(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
# Mac bits
  set(CMAKE_MACOSX_RPATH 1)
  set(CMAKE_INSTALL_RPATH "@executable_path/../lib;@executable_path/../opt/lib")
  set(CMAKE_C_FLAGS "-DMACINTOSH ${CMAKE_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "-DMACINTOSH  -Wno-enum-constexpr-conversion ${CMAKE_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "-DMACINTOSH ${CMAKE_Fortran_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "-Wl,-flat_namespace ${CMAKE_EXE_LINKER_FLAGS}")
  set(CMAKE_SHARED_LINKER_FLAGS "-Wl,-flat_namespace ${CMAKE_SHARED_LINKER_FLAGS}")
  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm64")
    set(EOS_ARCH macarm64)
  else()
    set(EOS_ARCH macintel64)
  endif()
else(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
# Miscellaneous - not sure what we are...
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib;$ORIGIN/../opt/lib")
  message("CMAKE_SYSTEM_NAME: " ${CMAKE_SYSTEM_NAME})
  message("System not recognized, using other")
  set(EOS_ARCH linux)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

# where, oh where are your FORTRAN mod files...
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran_modules)

# where we put the binaries
set(CMAKE_INSTALL_PREFIX $ENV{OCSSWROOT})

file(STRINGS ".manifest_tag" MANIFEST_TAG)
add_custom_target(update
          echo install_ocssw --tag ${MANIFEST_TAG} --install_dir ${CMAKE_SOURCE_DIR} --update
  COMMAND install_ocssw --tag ${MANIFEST_TAG} --install_dir ${CMAKE_SOURCE_DIR} --update
  COMMAND echo
  COMMAND echo update_luts --eval all
  COMMAND update_luts --eval all
  )

# set a default for lib3 search paths
set(HDF5_ROOT $ENV{LIB3_DIR})
set(HDF5_FIND_QUIETLY TRUE)

set(NETCDF_DIR $ENV{LIB3_DIR})

# make CGAL be quiet since I took care of this in the lib build
set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)

# add the submodule repos
add_subdirectory(oel_util)

add_subdirectory(oel_hdf4)

add_subdirectory(src)
add_subdirectory(test)

enable_testing()

