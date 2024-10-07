cmake_minimum_required(VERSION 3.12)

add_executable(pdsinfo
  pdsinfo.c
)

target_link_libraries(pdsinfo
  m
)

add_executable(pdsmerge
  pdsmerge.c
)

target_link_libraries(pdsmerge
  m
)

install (TARGETS pdsinfo pdsmerge DESTINATION bin)
