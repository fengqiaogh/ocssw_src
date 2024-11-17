cmake_minimum_required(VERSION 3.12)

find_package(HDF4 REQUIRED)
include_directories(${HDF4_INCLUDE_DIRS})
set(LIBS ${HDF4_LIBRARIES})

find_package(NetCDF REQUIRED COMPONENTS C CXX)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES} nav)

list(APPEND LIBS genutils)

ADD_LIBRARY(anc
  getanc.c
  HDFroutines.c
  dataintp.c
)

target_link_libraries(anc ${LIBS})
target_include_directories(anc PUBLIC ".")

install (TARGETS anc DESTINATION lib)
