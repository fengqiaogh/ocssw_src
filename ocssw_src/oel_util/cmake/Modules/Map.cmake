# Associative array library
#
# map(SET <map> <key> <value> [IF_EXISTS | IF_NOT_EXISTS [REQUIRED]])
# map(APPEND <map> <key> <value>)
# map(REMOVE <map> <key>)
# map(GET <map> <key> <output variable> [REQUIRED])
# map(FIND <map> <key> <output variable>)
#
#  This is a very weak implementation.  Current problems:
#  - No semi-colons or equal signs allowed.
#    - Equal signs are my fault, semi-colons are CMake's list delim (but can be escaped).
#  - Using APPEND allows duplicates (intentionally). 
#  - Linear array implementation won't scale too well.
#  - CMake doesn't allow random array editing, requiring a deletion and append.
#
#  Author's note: CMake's preferred argument ordering is super weird...

function(map command map)
	if (command STREQUAL FIND)
		set(find_key ${ARGV2})
		set(output ${ARGV3})

		set(index 0)
		foreach(kv_entry IN LISTS ${map})
			string(REGEX REPLACE "^([^=]+)=.*" "\\1" key "${kv_entry}")
			if (key STREQUAL find_key)
				set(${output} ${index} PARENT_SCOPE)
				return()
			endif ()
			math(EXPR index ${index}+1)
		endforeach()
		set(${output} -1 PARENT_SCOPE)
	elseif (command STREQUAL APPEND)
		list(APPEND ${map} "${ARGV2}=${ARGV3}")
		set(${map} "${${map}}" PARENT_SCOPE)
	elseif (command STREQUAL SET)
		set(key ${ARGV2})
		set(value ${ARGV3})

		map(FIND ${map} ${key} index)
		if (NOT index STREQUAL -1)
			if (NOT ARGV4 STREQUAL IF_NOT_EXISTS)
				list(REMOVE_AT ${map} ${index})
				list(APPEND ${map} "${key}=${value}")
				set(${map} ${${map}} PARENT_SCOPE)
			elseif (ARGV5 STREQUAL REQUIRED)
				message(FATAL_ERROR "\"${key}\" already exists in map")
			endif ()
		elseif (NOT ARGV4 STREQUAL IF_EXISTS)
			list(APPEND ${map} "${key}=${value}")
			set(${map} ${${map}} PARENT_SCOPE)
		elseif (ARGV5 STREQUAL REQUIRED)
			message(FATAL_ERROR "\"${key}\" doesn't exist in map")
		endif ()
		# map(APPEND ${map} ${key} ${value}) # can't use this 'cause PARENT_SCOPE only goes up level
	elseif (command STREQUAL REMOVE)
		set(key ${ARGV2})

		map(FIND ${map} ${key} index)
		if (NOT index STREQUAL -1)
			list(REMOVE_AT ${map} ${index})
		endif ()
		set(${map} ${${map}} PARENT_SCOPE)
	elseif (command STREQUAL GET)
		set(key ${ARGV2})
		set(output ${ARGV3})

		map(FIND ${map} ${key} index)
		if (index STREQUAL -1)
			if (ARGV4 STREQUAL REQUIRED)
				message(FATAL_ERROR "\"${key}\" not found in map")
			else ()
				set(${output} "" PARENT_SCOPE)
			endif ()
		else ()
			list(GET ${map} ${index} kv_entry)
			string(REGEX REPLACE "^[^=]+=(.*)" "\\1" value ${kv_entry})
			set(${output} ${value} PARENT_SCOPE)
		endif ()
	else ()
		message(FATAL_ERROR "Invalid map command '${command}'")
	endif ()
endfunction()

