
# original from http://stackoverflow.com/a/7788165
macro(get_subdirs result curdir)
	file(GLOB children RELATIVE ${curdir} ${curdir}/*)
	set(dirlist "")
	foreach(child ${children})
		if(IS_DIRECTORY ${curdir}/${child})
			list(APPEND dirlist ${child})
		endif()
	endforeach(child)
	set(${result} ${dirlist})
endmacro(get_subdirs)



# get_subdirs(SUBDIRS ${CMAKE_CURRENT_LIST_DIR})
# foreach(subdir ${SUBDIRS})
#	add_subdirectory(${subdir})
# endforeach(subdir)

