cmake_minimum_required(VERSION 3.12)

find_package(HDF4 REQUIRED)
include_directories(${HDF4_INCLUDE_DIRS})
set(LIBS ${HDF4_LIBRARIES})
  
add_executable(geocheck_modis
  geocheck_modis.c
)

target_link_libraries(geocheck_modis
  ${LIBS}
  m
)

install (TARGETS geocheck_modis DESTINATION bin)

