
cmake_minimum_required(VERSION 3.12)

project (oel_util)

set(LIBS piutils timeutils)

find_package(NetCDF REQUIRED COMPONENTS C CXX)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES})
  
add_library (netcdfutils
  nc_gridutils.c
  nc_init_compress.c
  nc4utils.c
  cdl_utils.cpp
  createNCDF.c
  global_attrs.cpp
  scaled_nc_var.cpp
  scaled_nc_group.cpp
  scaled_nc_file.cpp
  geo_region.cpp
)

add_executable(scaledNcVarTest
  scaled_nc_test.cpp
)

target_link_libraries(scaledNcVarTest netcdfutils ${LIBS})

target_link_libraries(netcdfutils ${LIBS})
target_include_directories(netcdfutils PUBLIC ".")


install (TARGETS netcdfutils DESTINATION lib) 
