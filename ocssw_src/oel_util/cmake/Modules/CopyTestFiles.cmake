
macro(copy_test_files)
	foreach(dir ${ARGN})
		file(COPY ${CMAKE_CURRENT_LIST_DIR}/${dir} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
	endforeach(dir)
endmacro(copy_test_files)
