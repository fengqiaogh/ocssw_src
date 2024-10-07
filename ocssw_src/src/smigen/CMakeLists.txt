cmake_minimum_required(VERSION 3.12)

add_executable(smigen
  smigen.cpp
  smigen_input.c
  put_smi.cpp
)

add_executable(fixrgb
  fixrgb.c
)

target_link_libraries(smigen
  map
  bin
)

install (TARGETS smigen fixrgb DESTINATION bin)
