
# http://stackoverflow.com/a/9328525/1236508
macro(print_all_variables)
	get_cmake_property(_variableNames VARIABLES)
	foreach (_variableName ${_variableNames})
	    message(STATUS "    ${_variableName}=${${_variableName}}")
	endforeach()
endmacro(print_all_variables)