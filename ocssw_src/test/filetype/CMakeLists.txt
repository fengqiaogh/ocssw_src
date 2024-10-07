enable_testing()

set(CHECK_MIN_VERSION 0.9.10)
#find_package(Check ${CHECK_MIN_VERSION})

#if (CHECK_FOUND AND NOT CHECK_MIN_VERSION VERSION_GREATER CHECK_VERSION)
if (FALSE)
    cdata(READ data_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
     GENERATE "${CMAKE_CURRENT_BINARY_DIR}/filetype-data"
     VERBOSE
     TEST filetype-data
     SAVE_DIR dataDir)

    include_directories(${CHECK_INCLUDE_DIRS})

    set(TESTS filetype)
    foreach(t ${TESTS})
        add_executable(${t}-test ${t}.c ${TEST_SOURCES})
        add_test(NAME ${t}-test COMMAND ${t}-test WORKING_DIRECTORY ${dataDir})
        target_link_libraries(${t}-test ${CHECK_LIBRARIES} piutils genutils filetype)
    endforeach(t)

    # add_test(<name> <command> [<arg>...])
    # add_test(NAME <name> COMMAND <command> [<arg>...] [CONFIGURATIONS <config>...] [WORKING_DIRECTORY <dir>])
else ()
    message(STATUS "tests/ directory requires Check>=${CHECK_MIN_VERSION}, no tests will be built")
endif ()
