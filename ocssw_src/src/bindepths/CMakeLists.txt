cmake_minimum_required(VERSION 3.12)

set(LIBS genutils)

find_package(HDF4 REQUIRED)
include_directories(${HDF4_INCLUDE_DIRS})
list(APPEND LIBS ${HDF4_LIBRARIES})

add_executable(bindepths
  bindepths.c
)

target_link_libraries(bindepths
  ${LIBS}
)

install (TARGETS bindepths DESTINATION bin)

