# - Try to find the PugiXML libraries
#  Once done this will define
#
#  PUGIXML_FOUND - system has PugiXML
#  PUGIXML_INCLUDE_DIRS - the PugiXML include directory
#  PUGIXML_LIBRARIES - the PugiXML library

set(_NAME PugiXML)
set(_HEADER pugixml.hpp)
set(_LIB pugixml)

string(TOUPPER ${_NAME} _NAME_UPPER)

find_package(PkgConfig)
if (PKG_CONFIG_FOUND)
    pkg_search_module(${_NAME_UPPER} QUIET ${_LIB})
endif (PKG_CONFIG_FOUND)

if (_HEADER)
  find_path(${_NAME_UPPER}_INCLUDE_DIR ${_HEADER} HINTS ${${_NAME_UPPER}_INCLUDEDIR} ${${_NAME_UPPER}_LIBRARY_INCLUDE_DIRS} $ENV{LIB3_DIR}/include)
endif (_HEADER)

if (_LIB)
    find_library(${_NAME_UPPER}_LIBRARY NAMES ${_LIB} HINTS ${${_NAME_UPPER}_LIBDIR} ${${_NAME_UPPER}_LIBRARY_DIRS} $ENV{LIB3_DIR}/lib)
endif (_LIB)

if (${_NAME_UPPER}_INCLUDE_DIR)
    file(READ "${${_NAME_UPPER}_INCLUDE_DIR}/${_HEADER}" H_CONTENTS LIMIT 1000)
    STRING(REGEX REPLACE ".*#.*define.*PUGIXML_VERSION +([0-9]+).*" "\\1" ${_NAME_UPPER}_VERSION "${H_CONTENTS}")
    STRING(REGEX REPLACE "([0-9]+)[0-9][0-9]" "\\1" ${_NAME_UPPER}_MAJOR_VERSION "${${_NAME_UPPER}_VERSION}")
    STRING(REGEX REPLACE "[0-9]+([0-9][0-9])" "\\1" ${_NAME_UPPER}_MINOR_VERSION "${${_NAME_UPPER}_VERSION}")
    set(${_NAME_UPPER}_MICRO_VERSION "0")
    unset(H_CONTENTS)
endif (${_NAME_UPPER}_INCLUDE_DIR)

set(_VALID_VERSION 1)

if (${_NAME_UPPER}_VERSION AND ${_NAME}_FIND_VERSION)
    if (${_NAME}_FIND_VERSION_EXACT)
        if (NOT ${_NAME_UPPER}_VERSION VERSION_EQUAL ${_NAME}_FIND_VERSION)
            set(_VALID_VERSION 0)
            set(_EQUAL_SIGN "=")
        endif (NOT ${_NAME_UPPER}_VERSION VERSION_EQUAL ${_NAME}_FIND_VERSION)
    else (${_NAME}_FIND_VERSION_EXACT)
        if (${_NAME_UPPER}_VERSION VERSION_LESS ${_NAME}_FIND_VERSION)
            set(_VALID_VERSION 0)
            set(_EQUAL_SIGN ">=")
        endif (${_NAME_UPPER}_VERSION VERSION_LESS ${_NAME}_FIND_VERSION)
    endif (${_NAME}_FIND_VERSION_EXACT)
endif (${_NAME_UPPER}_VERSION AND ${_NAME}_FIND_VERSION)

if (_VALID_VERSION)
    include(FindPackageHandleStandardArgs)
    if (CMAKE_MAJOR_VERSION LESS 3)
        find_package_handle_standard_args(${_NAME} DEFAULT_MSG ${_NAME_UPPER}_LIBRARY ${_NAME_UPPER}_INCLUDE_DIR)
    else (CMAKE_MAJOR_VERSION LESS 3)
        find_package_handle_standard_args(${_NAME}
            FOUND_VAR ${_NAME_UPPER}_FOUND
            REQUIRED_VARS ${_NAME_UPPER}_LIBRARY ${_NAME_UPPER}_INCLUDE_DIR
            VERSION_VAR ${_NAME_UPPER}_VERSION)
    endif (CMAKE_MAJOR_VERSION LESS 3)
else (_VALID_VERSION)
    set(_MESSAGE "Could NOT find ${_NAME}${_EQUAL_SIGN}${${_NAME}_FIND_VERSION} (found version \"${${_NAME_UPPER}_VERSION}\")")
    if (${_NAME}_FIND_REQUIRED)
        message(SEND_ERROR ${_MESSAGE})
    else (${_NAME}_FIND_REQUIRED)
        message(STATUS ${_MESSAGE})
    endif (${_NAME}_FIND_REQUIRED)
    unset(_MESSAGE)
    unset(_EQUAL_SIGN)
endif (_VALID_VERSION)
unset(_VALID_VERSION)

if (${_NAME_UPPER}_FOUND)
    set(${_NAME_UPPER}_INCLUDE_DIRS ${${_NAME_UPPER}_INCLUDE_DIR})
    set(${_NAME_UPPER}_LIBRARIES ${${_NAME_UPPER}_LIBRARY})
endif (${_NAME_UPPER}_FOUND)

mark_as_advanced(${_NAME_UPPER}_INCLUDE_DIR ${_NAME_UPPER}_LIBRARY)

unset(_NAME)
unset(_NAME_UPPER)
unset(_HEADER)
unset(_LIB)
