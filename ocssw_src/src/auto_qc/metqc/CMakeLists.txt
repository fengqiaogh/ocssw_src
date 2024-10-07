cmake_minimum_required(VERSION 3.12)

set(LIBS genutils)

find_package(HDF4 REQUIRED)
include_directories(${HDF4_INCLUDE_DIRS})
list(APPEND LIBS ${HDF4_LIBRARIES})

add_executable(metqc
  metqc.c
  rd_size.c
  resize_2d.c
  ANCroutines.c
  pexit.c
)

target_link_libraries(metqc
  hdf_dio
  ${LIBS}
)

install (TARGETS metqc DESTINATION bin)

