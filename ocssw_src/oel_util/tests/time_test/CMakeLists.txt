
enable_testing()

# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../src/pqueue)

add_executable(time_test time_test.c)
target_link_libraries(time_test timeutils)

# add_test(<name> <command> [<arg>...])
# add_test(NAME <name> COMMAND <command> [<arg>...] [CONFIGURATIONS <config>...] [WORKING_DIRECTORY <dir>])
