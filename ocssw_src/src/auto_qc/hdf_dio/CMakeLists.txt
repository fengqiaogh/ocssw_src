cmake_minimum_required(VERSION 3.12)

find_package(HDF4 REQUIRED)
include_directories(${HDF4_INCLUDE_DIRS})
set(LIBS ${HDF4_LIBRARIES})

add_library (hdf_dio STATIC
  l1io_close.c
  l1io_open.c
  l1io_read.c
  l2io_read.c
  read_sd.c
  read_sds.c
  open_hdf.c
  read_g_attr.c
  nav_read.c
)

target_link_libraries(hdf_dio ${LIBS} )
target_include_directories(hdf_dio PUBLIC ".")
