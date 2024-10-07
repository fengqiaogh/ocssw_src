
set(CHECK_MIN_VERSION 0.9.10)
find_package(Check ${CHECK_MIN_VERSION})
#find_package(Valgrind)

set(CMAKE_CXX_FLAGS "-g")

if (CHECK_FOUND AND NOT CHECK_MIN_VERSION VERSION_GREATER CHECK_VERSION)
    include(GetSubdirs)
    include(CopyTestFiles)

    # message(STATUS "Check ${CHECK_VERSION} found, building test suite")
    include_directories(SYSTEM ${CHECK_INCLUDE_DIRS})
    link_directories(${CHECK_LIBRARY_DIRS})

    if (VALGRIND_FOUND)
        add_definitions(-DHAVE_VALGRIND_H)
        set(VALGRIND_OPTIONS "--tool=memcheck --error-exitcode=1 --leak-check=full --show-leak-kinds=all --show-reachable=yes --trace-children=yes")
        set(memcheck_command "${VALGRIND_BIN} ${VALGRIND_OPTIONS}")
        separate_arguments(memcheck_command)
    endif (VALGRIND_FOUND)

    function(add_memcheck_test name)
        if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${name}.c")
            add_executable(${TEST_PREFIX}${name}-test ${name}.c ${TEST_SOURCES})
        elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc")
            add_executable(${TEST_PREFIX}${name}-test ${name}.cc ${TEST_SOURCES})
        elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${name}.cpp")
            add_executable(${TEST_PREFIX}${name}-test ${name}.cpp ${TEST_SOURCES})
        else ()
            message(FATAL_ERROR "No viable file found at ${CMAKE_CURRENT_SOURCE_DIR}/${name}.{c,cc,cpp}")
        endif ()
        add_test(NAME ${TEST_PREFIX}${name}-test COMMAND ${TEST_PREFIX}${name}-test ${ARGN} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        if (VALGRIND_FOUND)
            add_test(NAME ${TEST_PREFIX}${name}-test-memcheck COMMAND ${memcheck_command} ${CMAKE_CURRENT_BINARY_DIR}/${TEST_PREFIX}${name}-test ${ARGN} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        endif (VALGRIND_FOUND)
    endfunction(add_memcheck_test)

    function(set_memcheck_test_properties name)
        set_tests_properties(${TEST_PREFIX}${name}-test ${ARGN})
        if (VALGRIND_FOUND)
            set_tests_properties(${TEST_PREFIX}${name}-test-memcheck ${ARGN})
        endif (VALGRIND_FOUND)
    endfunction(set_memcheck_test_properties)

    function(memcheck_target_link_libraries name)
        target_link_libraries(${TEST_PREFIX}${name}-test ${ARGN})
    endfunction(memcheck_target_link_libraries)

    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

    enable_testing()

    get_subdirs(SUBDIRS ${CMAKE_CURRENT_LIST_DIR})
    foreach(subdir ${SUBDIRS})
        add_subdirectory(${subdir})
    endforeach(subdir)
else (CHECK_FOUND AND NOT CHECK_MIN_VERSION VERSION_GREATER CHECK_VERSION)
    message(STATUS "tests/ directory requires Check>=${CHECK_MIN_VERSION}, no tests will be built")
endif (CHECK_FOUND AND NOT CHECK_MIN_VERSION VERSION_GREATER CHECK_VERSION)


