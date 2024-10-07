cmake_minimum_required(VERSION 3.12)

set(LIBS nav timeutils genutils)

find_package(PROJ REQUIRED CONFIG)
list(APPEND LIBS PROJ::proj)
  
ADD_LIBRARY(bilio
  jplaeriallib.c
  read_aviris.c
  read_prism.c
)

target_link_libraries(bilio ${LIBS})
target_include_directories(bilio PUBLIC ".")

install (TARGETS bilio DESTINATION lib)
