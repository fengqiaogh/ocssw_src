cmake_minimum_required(VERSION 3.12)

add_executable(l3stat_chk
  l3stat_chk.c
  hist.c
)

#add_definitions(-ansi)

include_directories (
  inc
)

target_link_libraries(l3stat_chk
  hdf4utils
  hdf_dio
)

install (TARGETS l3stat_chk DESTINATION bin)

