cmake_minimum_required(VERSION 3.12)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(LIBS dfutils)

find_package(NetCDF REQUIRED COMPONENTS C CXX)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES})

ADD_LIBRARY(bin
  read_l3b_meta.c
  get_beg_ext.c
  wr_vdata.c
  write_l3b_meta.c
  ncdfbin_utils.c

  bin_io.cpp
  bin_util.cpp
  L3Shape.cpp
  L3ShapeIsine.cpp
  L3ShapeSMI.cpp
  L3File.cpp
  L3FileSMI.cpp
)

target_link_libraries(bin ${LIBS})
target_include_directories(bin PUBLIC ".")

install (TARGETS bin DESTINATION lib)
