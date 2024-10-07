cmake_minimum_required(VERSION 3.12)

find_package(HDF4 REQUIRED)
include_directories(${HDF4_INCLUDE_DIRS})
set(LIBS ${HDF4_LIBRARIES})

add_executable(l1aextract_modis
  l1aextract_modis.c
)

target_link_libraries(l1aextract_modis
  ${LIBS}
)

install (TARGETS l1aextract_modis DESTINATION bin)

