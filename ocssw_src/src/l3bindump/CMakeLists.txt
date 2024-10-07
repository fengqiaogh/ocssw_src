cmake_minimum_required(VERSION 3.12)

#
# Make a version file containing the current version from git.
#
set (VERSION_MAJOR 3)
set (VERSION_MINOR 1)
set (VERSION_PATCH 1)
include(GetOcsswVersion)

add_executable(l3bindump
  l3bindump_input.c
  compareObj.cpp
  main_l3bindump.cpp
)

target_link_libraries(l3bindump
  bin
)

install (TARGETS l3bindump DESTINATION bin)

