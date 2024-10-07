cmake_minimum_required(VERSION 3.12)

set(LIBS genutils)

find_package(HDF4 REQUIRED)
include_directories(${HDF4_INCLUDE_DIRS})
list(APPEND LIBS ${HDF4_LIBRARIES})

add_executable(o3qc
  o3qc.c
  anc_daymon.c
  day2mday.c
  ANCroutines.c
  pexit.c
)

add_definitions(-DHDF)

include_directories (
  inc
)

target_link_libraries(o3qc
  ${LIBS}
  hdf_dio
  m
)

install (TARGETS o3qc DESTINATION bin)

