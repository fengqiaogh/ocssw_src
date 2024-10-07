# - Try to find the VALGRIND libraries
#  Once done this will define
#
#  VALGRIND_FOUND - system has check
#  VALGRIND_INCLUDE_DIRS - the check include directory
#  VALGRIND_BIN - valgrind command

set(_NAME Valgrind)
set(_HEADER valgrind/valgrind.h)
set(_BIN valgrind)
set(_LIB valgrind)

string(TOUPPER ${_NAME} _NAME_UPPER)

find_package(PkgConfig)
if (PKG_CONFIG_FOUND)
	pkg_search_module(${_NAME_UPPER} QUIET ${_LIB})
endif (PKG_CONFIG_FOUND)

if (_HEADER)
	find_path(${_NAME_UPPER}_INCLUDE_DIR ${_HEADER} HINTS ${${_NAME_UPPER}_INCLUDEDIR} ${${_NAME_UPPER}_LIBRARY_INCLUDE_DIRS})
	set(${_NAME_UPPER}_INCLUDE_DIRS ${${_NAME_UPPER}_INCLUDE_DIR})
endif (_HEADER)

if (_LIB)
	find_library(${_NAME_UPPER}_LIBRARY NAMES ${_LIB} HINTS ${${_NAME_UPPER}_LIBDIR} ${${_NAME_UPPER}_LIBRARY_DIRS})
	set(${_NAME_UPPER}_LIBRARIES ${${_NAME_UPPER}_LIBRARY})
endif (_LIB)

if (_BIN)
	find_program(${_NAME_UPPER}_BIN NAMES ${_BIN} HINTS ENV PATH)
endif (_BIN)

if (${_NAME_UPPER}_INCLUDE_DIR/valgrind/config.h)
    file(READ "${${_NAME_UPPER}_INCLUDE_DIR}/valgrind/config.h" H_CONTENTS)
    STRING(REGEX REPLACE ".*#define +VERSION +\"([0-9]+\\.[0-9]+\\.[0-9]+)\".*" "\\1" ${_NAME_UPPER}_VERSION "${H_CONTENTS}")
    STRING(REGEX REPLACE ".*([0-9]+)\\.[0-9]+\\.[0-9]+.*" "\\1" ${_NAME_UPPER}_MAJOR_VERSION "${${_NAME_UPPER}_VERSION}")
    STRING(REGEX REPLACE ".*[0-9]+\\.([0-9]+)\\.[0-9]+.*" "\\1" ${_NAME_UPPER}_MINOR_VERSION "${${_NAME_UPPER}_VERSION}")
    STRING(REGEX REPLACE ".*[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" ${_NAME_UPPER}_MICRO_VERSION "${${_NAME_UPPER}_VERSION}")
    unset(H_CONTENTS)
endif (${_NAME_UPPER}_INCLUDE_DIR/valgrind/config.h)

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
		find_package_handle_standard_args(${_NAME} DEFAULT_MSG ${_NAME_UPPER}_BIN ${_NAME_UPPER}_INCLUDE_DIR)
	else (CMAKE_MAJOR_VERSION LESS 3)
		find_package_handle_standard_args(${_NAME}
			FOUND_VAR ${_NAME_UPPER}_FOUND
			REQUIRED_VARS ${_NAME_UPPER}_BIN ${_NAME_UPPER}_INCLUDE_DIR
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

mark_as_advanced(${_NAME_UPPER}_INCLUDE_DIR ${_NAME_UPPER}_LIBRARY)

unset(_NAME)
unset(_NAME_UPPER)
unset(_HEADER)
unset(_LIB)
unset(_BIN)
