
enable_testing()

# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../src/phash)
	
set(TESTS alloc)
foreach(t ${TESTS})
	add_memcheck_test(${t})
	memcheck_target_link_libraries(${t} ${CHECK_LIBRARIES} genutils)
endforeach(t)

# add_test(<name> <command> [<arg>...])
# add_test(NAME <name> COMMAND <command> [<arg>...] [CONFIGURATIONS <config>...] [WORKING_DIRECTORY <dir>])
# target_link_libraries(<target> [item1 [item2 [...]]] [[debug|optimized|general] <item>] ...)
