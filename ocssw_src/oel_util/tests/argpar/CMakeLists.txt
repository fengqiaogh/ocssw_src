
enable_testing()

# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../src/argpar)
find_package(Jansson REQUIRED)

copy_test_files(parfiles)

set(TESTS argpar)
foreach(t ${TESTS})
	add_memcheck_test(${t})
	memcheck_target_link_libraries(${t} ${CHECK_LIBRARIES} genutils ${JANSSON_LIBRARIES})
endforeach(t)

# these aren't real tests, just ones for visually inspecting output
add_executable(argpar-help argpar-help.c)
target_link_libraries(argpar-help genutils ${JANSSON_LIBRARIES})

add_executable(argpar-params argpar-params.c)
target_link_libraries(argpar-params genutils ${JANSSON_LIBRARIES})
