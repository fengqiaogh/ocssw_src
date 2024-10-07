cmake_minimum_required(VERSION 3.12)

set(LIBS filetype l1 nav libterrain)

find_package(NetCDF REQUIRED COMPONENTS C CXX)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES})

find_package(GeographicLib REQUIRED)
include_directories(${GeographicLib_INCLUDE_DIRS})
list(APPEND LIBS ${GeographicLib_LIBRARIES})

add_executable(l1cgen
  l1cgen.cpp
  l1c.cpp
  l1c_input.cpp
  l1c_filehandle.cpp
  l1c_str.cpp
  hawkeye_methods.cpp
  l2_str.cpp
  l1c_bin.cpp
  l1c_latlongrid.cpp
  )

target_link_libraries(l1cgen
  ${LIBS}
  )



install (TARGETS l1cgen DESTINATION bin)
