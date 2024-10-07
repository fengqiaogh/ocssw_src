cmake_minimum_required(VERSION 3.12)

set(LIBS nav)

ADD_LIBRARY(czcs
  czcscaleg.f
  czsuba.f
  satang.f
  subsatxy.f
)

target_link_libraries(czcs ${LIBS})
target_include_directories(czcs PUBLIC ".")

install (TARGETS czcs DESTINATION lib)
