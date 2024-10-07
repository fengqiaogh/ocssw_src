cmake_minimum_required(VERSION 3.12)

find_package(HDF4 REQUIRED)
include_directories(${HDF4_INCLUDE_DIRS})
set(LIBS ${HDF4_LIBRARIES})

add_executable(l2stat_chk
  l2stat_chk.c
)

add_definitions(-ansi)

target_link_libraries(l2stat_chk
  ${LIBS}
  hdf_dio
  m
  seawifs
)

install (TARGETS l2stat_chk DESTINATION bin)

