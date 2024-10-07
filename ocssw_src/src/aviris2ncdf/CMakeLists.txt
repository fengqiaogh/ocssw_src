cmake_minimum_required(VERSION 3.12)

find_package(NetCDF REQUIRED)
include_directories(${NETCDF_INCLUDE_DIRS})
set(LIBS ${NETCDF_LIBRARIES})

list(APPEND LIBS bilio netcdfutils)

add_executable(avirisbil2ncdf
  avirisbil2ncdf.cpp
)

target_link_libraries(avirisbil2ncdf ${LIBS})

install (TARGETS avirisbil2ncdf DESTINATION bin)
