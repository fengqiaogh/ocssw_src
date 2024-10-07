cmake_minimum_required(VERSION 3.12)

set(LIBS timeutils genutils)

find_package(HDF4 REQUIRED)
include_directories(${HDF4_INCLUDE_DIRS})
list(APPEND LIBS ${HDF4_LIBRARIES})

add_executable(l1asubset_modis
  l1asubset_modis.c
)

target_link_libraries(l1asubset_modis
  ${LIBS}
)

install (TARGETS l1asubset_modis DESTINATION bin)
