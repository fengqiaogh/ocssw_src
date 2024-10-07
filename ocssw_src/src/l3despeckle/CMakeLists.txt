cmake_minimum_required(VERSION 3.12)

set(LIBS bin++)

add_executable(l3despeckle
  l3despeckle.cpp
  l3despeckle_input.c
)

add_executable(apply_despeckle
  apply_despeckle.cpp
)

target_link_libraries(l3despeckle ${LIBS})

target_link_libraries(apply_despeckle ${LIBS})

install (TARGETS l3despeckle apply_despeckle DESTINATION bin)

