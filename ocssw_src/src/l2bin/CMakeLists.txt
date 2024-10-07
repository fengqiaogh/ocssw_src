cmake_minimum_required(VERSION 3.12)

set(LIBS l2 nav bin dataday)

find_package(NetCDF REQUIRED COMPONENTS C CXX)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES})

add_executable(l2bin
  l2bin.cpp
  l2bin_input.cpp
  )

target_link_libraries(l2bin
  ${LIBS}
  ${LibBin_LIBRARIES}
)

install (TARGETS l2bin DESTINATION bin)
