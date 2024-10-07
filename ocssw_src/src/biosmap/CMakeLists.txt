cmake_minimum_required(VERSION 3.12)

add_executable(biosmap
  main_biosmap.c
)

target_link_libraries(biosmap
  map
  bin
  )

install (TARGETS biosmap DESTINATION bin)

