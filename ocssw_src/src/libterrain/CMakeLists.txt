cmake_minimum_required(VERSION 3.12)


set(LIBS netcdfutils)

find_package(NetCDF REQUIRED)
find_package(HDF4 REQUIRED)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES} ${HDF4_LIBRARIES})

ADD_LIBRARY(libterrain
  get_dem_height.c
  get_nc_height.c
)

target_link_libraries(libterrain ${LIBS})
target_include_directories(libterrain PUBLIC ".")

install (TARGETS libterrain DESTINATION lib)

