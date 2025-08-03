cmake_minimum_required(VERSION 3.12)

set(LIBS netcdfutils)

find_package(NetCDF REQUIRED COMPONENTS CXX)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES})

add_executable(l1agen_oci
  l0stream.cpp
  process_science_data.cpp
  l1agen_oci.cpp
  l1afile_oci.cpp
  l0PacketUtil_oci.cpp
  l1afile_oci_manager.cpp
  dimension_shape.cpp
)

target_link_libraries(l1agen_oci
  ${LIBS}
)

install (TARGETS l1agen_oci DESTINATION bin)
