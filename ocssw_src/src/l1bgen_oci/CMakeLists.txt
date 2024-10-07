cmake_minimum_required(VERSION 3.12)

set(LIBS genutils netcdfutils libterrain)

find_package(NetCDF REQUIRED COMPONENTS C CXX)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
list(APPEND LIBS ${GSL_LIBRARIES})

add_executable(l1bgen_oci
  l1bgen_oci.cpp
  geolocate_oci.cpp
  common.cpp
)

#add_executable(geolocate_oci
#  geolocate_oci.cpp
#  common.cpp
#)

target_link_libraries(l1bgen_oci
  ${LIBS}
)

#target_link_libraries(geolocate_oci
#  ${LIBS}
#)

#install (TARGETS l1bgen_oci geolocate_oci DESTINATION bin)
install (TARGETS l1bgen_oci DESTINATION bin)








