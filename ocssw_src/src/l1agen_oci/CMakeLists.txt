cmake_minimum_required(VERSION 3.12)
find_package(NetCDF REQUIRED)
include_directories(${NETCDF_INCLUDE_DIRS})
set(LIBS ${NETCDF_LIBRARIES} netcdfutils)
find_package(NetCDF REQUIRED COMPONENTS CXX)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES})

add_executable(l1agen_oci
  l0stream.cpp
  l1agen_oci.cpp
  common.cpp
)

# Using a python script in ocssw/src/scripts/l0info_oci.py instead
#
#add_executable(l0info_oci
#  l0info_oci.cpp
#  common.cpp
#)

target_link_libraries(l1agen_oci
  ${LIBS}
)

# Using a python script in ocssw/src/scripts/l0info_oci.py instead
#
#target_link_libraries(l0info_oci
#  ${LIBS}
#)

install (TARGETS l1agen_oci DESTINATION bin)
