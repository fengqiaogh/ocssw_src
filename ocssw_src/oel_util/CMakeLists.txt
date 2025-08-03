cmake_minimum_required(VERSION 3.12)

project (OEL_UTIL)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")

add_subdirectory("libgenutils")
add_subdirectory("libtimeutils")
add_subdirectory("libpiutils")
add_subdirectory("libnetcdfutils")
add_subdirectory("libgring")

if(NOT ${CMAKE_PROJECT_NAME} MATCHES "VIIRS_L1")
    add_subdirectory("libcgal_interp")

    add_subdirectory("tests")
    enable_testing()
endif()
