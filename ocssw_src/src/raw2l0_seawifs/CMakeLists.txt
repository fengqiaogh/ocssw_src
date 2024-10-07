cmake_minimum_required(VERSION 3.12)

add_executable(raw2l0_seawifs
  raw2L0.c
)

target_link_libraries(raw2l0_seawifs
  genutils
)

install (TARGETS raw2l0_seawifs DESTINATION bin)

