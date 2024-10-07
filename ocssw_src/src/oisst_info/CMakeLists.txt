cmake_minimum_required(VERSION 3.12)

add_executable(oisst_info
  oisst_info.c
)

target_link_libraries(oisst_info
  genutils
)

install (TARGETS oisst_info DESTINATION bin)

