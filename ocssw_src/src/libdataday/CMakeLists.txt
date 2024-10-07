cmake_minimum_required(VERSION 3.12)

#
# Make a version file containing the current version from git.
#
set (VERSION_MAJOR 1)
set (VERSION_MINOR 2)
set (VERSION_PATCH 0)
include(GetOcsswVersion)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
# set(Boost_DEBUG 1)
find_package(Boost REQUIRED)

find_package(NetCDF REQUIRED COMPONENTS CXX)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES})
list(APPEND LIBS timeutils)
list(APPEND LIBS piutils)
list(APPEND LIBS l2)

add_library(dataday
  get_dataday.cpp
)



target_link_libraries(dataday ${LIBS})
target_include_directories(dataday SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})

target_include_directories(dataday PUBLIC ".")
install (TARGETS dataday DESTINATION lib)


