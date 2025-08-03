cmake_minimum_required(VERSION 3.12)

set(LIBS genutils netcdfutils libterrain nav gring)

find_package(NetCDF REQUIRED COMPONENTS C CXX)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
list(APPEND LIBS ${GSL_LIBRARIES})

add_executable(l1bgen_oci
  device.cpp
  aggregate.cpp
  corrections.cpp
  calibrations.cpp
  geolocate_oci.cpp
  l1b_file.cpp
  dark_data.cpp
  l1bgen_oci.cpp
  processing_tracker.cpp
  l1b_options.cpp
)

target_link_libraries(l1bgen_oci
  ${LIBS}
)

install (TARGETS l1bgen_oci DESTINATION bin)
