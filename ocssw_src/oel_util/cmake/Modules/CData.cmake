# CData
#
# Read CData manifest files:
#
# READ var [CData-file ...] [APPEND]
#
# ---------------------------------------------------------------------------------------------------------------------
# Generate a target or a test to be run at make- or make test-time:
#
# GENERATE filename
#	READ var [CData-file ...] [APPEND]
#	DOWNLOAD var ...
#	URL url
# 	TARGET name [DEPENDENTS dep ...] [ALWAYS_RUN]
# 	TEST name
# 	COMMAND command
#	WORKING_DIRECTORY dir (for TEST/TARGET)
#	DIR dir
#	DOWNLOAD_BASE dir
#	STRIP_PARTS #
#	NO_MKDIR
#	SAVE_DIR var
#
# Only URL, either READ or DOWNLOAD, and either TEST or TARGET are required. (These "either" arguments are not
#	mutually exclusive.)
#
# If WORKING_DIRECTORY is supplied, a test and/or target is added that ensures the directory exists before it's 
#	associated test/target is run (a command cannot be run from a directory that doesn't exist).  NO_MKDIR 
#	disables this sanity measure.
#
# COMMAND is passed unfiltered to execute_process(COMMAND ${COMMAND}). This accepts multiple arguments, allowing you
#	to pass arguments to said command.  A FATAL_ERROR is thrown if the command returns a non-zero.
#
# See below for documentation on how DIR, DOWNLOAD_BASE, and STRIP_PARTS interact.  The final directory is stored
#	in SAVE_DIR, if requested.  This is useful for setting the WORKING_DIRECTORY of dependent tests.
#
# ---------------------------------------------------------------------------------------------------------------------
# Download data at CMake-time:
#
# DOWNLOAD var ... [VERBOSE]
#	READ var [CData-file ...] [APPEND]
#	URL url
#	DIR dir
#	DOWNLOAD_BASE dir
#	STRIP_PARTS #
#	VERBOSE
#
# Only URL is required.
#
# See below for documentation on how DIR, DOWNLOAD_BASE, and STRIP_PARTS interact.  The final directory is stored
#	in SAVE_DIR, if requested.  This is useful for setting the WORKING_DIRECTORY of dependent tests.
#
# ---------------------------------------------------------------------------------------------------------------------
# With any command, if READ is supplied but no files are added after <var>, CDATA_DEFAULT_FILE or CDataList.txt is
#	searched for.  If found, those files are downloaded, as well.
# With the READ option, <var> is used as an output variable unless APPEND is used, in which it's both an input and an
#	output.
#
# With DOWNLOAD or GENERATE, defaults can be set for a few options:
#	CDATA_DEFAULT_URL
#	CDATA_DEFAULT_BASE (DOWNLOAD_BASE)
#	CDATA_DEFAULT_NO_MKDIR
#	CDATA_DEFAULT_STRIP_PARTS
# 
# if DIR is supplied to DOWNLOAD/GENERATE:
#	download relative to DIR
# else:
#	if DOWNLOAD_BASE (or CDATA_DEFAULT_BASE) is supplied:
#		find the current source directory path relative to CMake's base directory
#		if STRIP_PARTS (or CDATA_DEFAULT_STRIP_PARTS) is supplied:
#			strip off parts at the beginning of the directory accordingly
#		download relative to whatever's left
#	else:
#		download relative to CMAKE_CURRENT_BINARY_DIR
#
# =====================================================================================================================
# EXAMPLES
# =====================================================================================================================
#
# For generating tests:
# ---------------------------------------------------------------------------------------------------------------------
# In test/CMakeLists.txt:
# -----------------------
# include(CData)
#
# set(CDATA_DEFAULT_URL "http://cdatahost/ocssw_test/")
# set(CDATA_DEFAULT_STRIP_PARTS 1)
# set(CDATA_DEFAULT_BASE "${CMAKE_SOURCE_DIR}/testdata")
# #set(CDATA_DEFAULT_NO_MKDIR 1) # might be unsafe, but makes cleaner testing/targets
#
# -------------------------
# In test/*/CMakeLists.txt:
# -------------------------
# Optional, for files common to multiple tests:
# cdata(
#	READ common_files  "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
#	SAVE_DIR working_dir
# )
# cdata(
#	READ common_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt-test1"
#	GENERATE "${CMAKE_CURRENT_BINARY_DIR}/test1"
#	TEST test1
# )
#
# add_executable(test2 test.c)
# cdata(
#	DOWNLOAD common_files
#	GENERATE "${CMAKE_CURRENT_BINARY_DIR}/test2"
#	TEST test2
#	EXECUTE test2 WORKING_DIRECTORY working_dir
# )
#

include(Map)
include(CMakeParseArguments)

function(CData)
	set(options APPEND VERBOSE ALWAYS_RUN NO_MKDIR)
	set(oneValueArgs URL DIR GENERATE TARGET WORKING_DIRECTORY TEST DOWNLOAD_BASE SAVE_DIR STRIP_PARTS)
	set(multiValueArgs READ DOWNLOAD COMMAND DEPENDENTS)
	cmake_parse_arguments(CDATA "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
	
	if (NOT DEFINED CDATA_DOWNLOAD_BASE AND DEFINED CDATA_DEFAULT_BASE)
		set(CDATA_DOWNLOAD_BASE ${CDATA_DEFAULT_BASE})
	endif ()
	if (NOT DEFINED CDATA_STRIP_PARTS AND DEFINED CDATA_DEFAULT_STRIP_PARTS)
		set(CDATA_STRIP_PARTS ${CDATA_DEFAULT_STRIP_PARTS})
	endif ()
	if (NOT DEFINED CDATA_URL AND DEFINED CDATA_DEFAULT_URL)
		set(CDATA_URL ${CDATA_DEFAULT_URL})
	endif ()
	if (NOT DEFINED CDATA_NO_MKDIR AND DEFINED CDATA_DEFAULT_NO_MKDIR)
		set(CDATA_NO_MKDIR ${CDATA_DEFAULT_NO_MKDIR})
	endif ()
	if (NOT CDATA_COMMAND OR NOT CDATA_WORKING_DIRECTORY)
		set(CDATA_WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
		set(CDATA_NO_MKDIR 1)
	endif ()

	set(read_result "")
	if (CDATA_READ)
		list(GET CDATA_READ 0 ret)
		if (CDATA_APPEND)
			set(read_result ${${ret}})
		endif ()

		list(REMOVE_AT CDATA_READ 0)
		
		if (NOT CDATA_READ)
			if (CDATA_DEFAULT_FILE)
				list(APPEND CDATA_READ ${CDATA_DEFAULT_FILE})
			else ()
				list(APPEND CDATA_READ "CDataList.txt")
			endif ()
		endif ()

		foreach(list ${CDATA_READ})
			if (EXISTS ${list})
				file(READ ${list} list_contents)
				string(REGEX REPLACE ";" "\\\\;" list_contents "${list_contents}")
				string(REGEX REPLACE "\n" ";" list_contents "${list_contents}")
				
				foreach(line ${list_contents})
					# local-path remote-path git-rev sha1-hash
					if ("${line}" MATCHES "^#")
						# no-op
					elseif ("${line}" MATCHES "^(.+) (.+) (.+) (.+)$")
						map(SET read_result ${CMAKE_MATCH_1} "${CMAKE_MATCH_3}/${CMAKE_MATCH_2} ${CMAKE_MATCH_4}")
					else ()
						message(FATAL_ERROR "Invalid data file line: ${line}")
					endif ()
				endforeach(line)
			else ()
				message(FATAL_ERROR "${list} does not exist")
			endif (EXISTS ${list})
		endforeach(list)

		set(${ret} ${read_result} PARENT_SCOPE)
	endif ()

	if (CDATA_DOWNLOAD OR CDATA_GENERATE)
		if (NOT CDATA_URL)
			message(FATAL_ERROR "No URL given to download data files")
		endif ()
		if (NOT CDATA_DIR)
			if (CDATA_DOWNLOAD_BASE)
				file(RELATIVE_PATH rel ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

				if (CDATA_STRIP_PARTS)
					math(EXPR CDATA_STRIP_PARTS "${CDATA_STRIP_PARTS} - 1")
					foreach (p RANGE ${CDATA_STRIP_PARTS})
						string(REGEX REPLACE "^[^/]+/" "" rel "${rel}")
					endforeach()
				endif ()

				set(CDATA_DIR "${CDATA_DOWNLOAD_BASE}/${rel}")
			else ()
				set(CDATA_DIR ${CMAKE_CURRENT_BINARY_DIR})
			endif ()
		endif ()
		if (CDATA_SAVE_DIR)
			set(${CDATA_SAVE_DIR} ${CDATA_DIR} PARENT_SCOPE)
		endif ()
	endif ()

	if (CDATA_GENERATE)
		set(CDATA_GENERATE_ORIG "${CDATA_GENERATE}")
		set(CDATA_GENERATE "${CDATA_GENERATE}.cmake")
		file(WRITE ${CDATA_GENERATE} "# CMake file generated to download any required data files\n\n")

		file(APPEND ${CDATA_GENERATE} "if (NOT EXISTS \"${CDATA_DIR}/cmake_install.cmake\")\n")
		file(APPEND ${CDATA_GENERATE} "\tfile(WRITE \"${CDATA_DIR}/cmake_install.cmake\" \"# Install script for directory: ${CMAKE_CURRENT_LIST_DIR}\\n\")\n")
		file(APPEND ${CDATA_GENERATE} "\tfile(APPEND \"${CDATA_DIR}/cmake_install.cmake\" \"# This is fake, generated by CData to point to who created it\\n\")\n")
		file(APPEND ${CDATA_GENERATE} "endif ()\n")
		file(APPEND ${CDATA_GENERATE} "\n")
		
		foreach(data_files ${CDATA_DOWNLOAD} read_result)
			foreach(path ${${data_files}})
				if (path MATCHES "^([^=]+)=(.+) (.+)$")
					file(APPEND ${CDATA_GENERATE} "file(DOWNLOAD \"${CDATA_URL}${CMAKE_MATCH_2}\" \"${CDATA_DIR}/${CMAKE_MATCH_1}\" SHOW_PROGRESS EXPECTED_HASH SHA1=${CMAKE_MATCH_3})\n")
				else ()
					message(FATAL_ERROR "Malformed file map ${path} ${data_files}")
				endif()
			endforeach(path)
		endforeach(data_files)

		if (CDATA_COMMAND)
			file(APPEND ${CDATA_GENERATE} "\n")
			file(APPEND ${CDATA_GENERATE} "execute_process(COMMAND ${CDATA_COMMAND} RESULT_VARIABLE res)\n")
			file(APPEND ${CDATA_GENERATE} "if (\${res})\n")
			file(APPEND ${CDATA_GENERATE} "\tmessage(FATAL_ERROR \"execution failed: \${res}\")\n")
			file(APPEND ${CDATA_GENERATE} "else ()\n")
			file(APPEND ${CDATA_GENERATE} "\tmessage(\"execution succeeded: \${res}\")\n")
			file(APPEND ${CDATA_GENERATE} "endif ()\n")
		endif ()
		
		if (NOT CDATA_NO_MKDIR)
			file(WRITE "${CDATA_GENERATE_ORIG}-mkdir.cmake" "file(MAKE_DIRECTORY \"${CDATA_WORKING_DIRECTORY}\")\n")
		endif ()

		if (CDATA_TARGET)
			if (CDATA_ALWAYS_RUN)
				add_custom_target(${CDATA_TARGET} ALL COMMAND "${CMAKE_COMMAND}" ARGS -P ${CDATA_GENERATE} WORKING_DIRECTORY ${CDATA_WORKING_DIRECTORY})
			else ()
				add_custom_target(${CDATA_TARGET} COMMAND "${CMAKE_COMMAND}" ARGS -P ${CDATA_GENERATE} WORKING_DIRECTORY ${CDATA_WORKING_DIRECTORY})
			endif ()
			if (NOT CDATA_NO_MKDIR)
				add_custom_target("${CDATA_TARGET}-mkdir" COMMAND "${CMAKE_COMMAND}" ARGS -P "${CDATA_GENERATE_ORIG}-mkdir.cmake")
				add_dependencies(${CDATA_TARGET} "${CDATA_TARGET}-mkdir")
			endif ()
			if (CDATA_DEPENDENTS)
				foreach(dependent ${CDATA_DEPENDENTS})
					add_dependencies(${dependent} ${CDATA_TARGET})
				endforeach(dependent)
			endif ()
		endif ()
		if (CDATA_TEST)
			if (NOT CDATA_NO_MKDIR)
				add_test(NAME "${CDATA_TEST}-mkdir" COMMAND "${CMAKE_COMMAND}" ARGS -P "${CDATA_GENERATE_ORIG}-mkdir.cmake")
			endif ()
			add_test(NAME ${CDATA_TEST} COMMAND "${CMAKE_COMMAND}" ARGS -P ${CDATA_GENERATE} WORKING_DIRECTORY ${CDATA_WORKING_DIRECTORY})
		endif ()
	elseif (CDATA_DOWNLOAD)
		foreach(map ${CDATA_DOWNLOAD} read_result)
			foreach(path ${${map}})
				if (path MATCHES "^([^=]+)=(.+) (.+)$")
					if (CDATA_VERBOSE)
						message(STATUS "Downloading ${CDATA_URL}${CMAKE_MATCH_2} => ${CDATA_DIR}/${CMAKE_MATCH_1}")
					endif ()
					file(DOWNLOAD "${CDATA_URL}${CMAKE_MATCH_2}" "${CDATA_DIR}/${CMAKE_MATCH_1}" SHOW_PROGRESS EXPECTED_HASH SHA1=${CMAKE_MATCH_3})
				else ()
					message(FATAL_ERROR "Malformed file map")
				endif()
			endforeach(path)
		endforeach(map)
	endif ()
endfunction(CData)

