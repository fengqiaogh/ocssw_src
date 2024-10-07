cmake_minimum_required(VERSION 3.12)

#
# Make a version file containing the current version from git.
#
set (VERSION_MAJOR 1)
set (VERSION_MINOR 1)
set (VERSION_PATCH 0)
include(GetOcsswVersion)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(dataday PUBLIC ".")




add_executable(get_dataday
main_dataday.cpp
)


target_link_libraries(get_dataday dataday)
install (TARGETS get_dataday DESTINATION bin)

