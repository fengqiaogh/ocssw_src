cmake_minimum_required(VERSION 3.12)

find_package(HDF4 REQUIRED)
include_directories(${HDF4_INCLUDE_DIRS})
set(LIBS ${HDF4_LIBRARIES})
  
add_executable(met_reverse
  met_reverse.c
  hio_close.c
  hio_open.c
  hio_r_sds.c
  hio_i_sds.c
)

target_link_libraries(met_reverse ${LIBS})

install (TARGETS met_reverse DESTINATION bin)

