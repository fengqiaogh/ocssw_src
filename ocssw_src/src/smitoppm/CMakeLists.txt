cmake_minimum_required(VERSION 3.12)

add_executable(smitoppm
  smitoppm.c
)

target_link_libraries(smitoppm
  map
  bin
)

install (TARGETS smitoppm DESTINATION bin)

