cmake_minimum_required(VERSION 3.12)

set(LIBS dfutils bin)

ADD_LIBRARY(map
  get_l3m.c
  put_l3m.c
)

target_link_libraries(map ${LIBS})
target_include_directories(map PUBLIC ".")

install (TARGETS map DESTINATION lib)
