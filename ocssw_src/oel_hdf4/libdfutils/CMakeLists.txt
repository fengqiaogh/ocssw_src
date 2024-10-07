cmake_minimum_required(VERSION 3.12)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(LIBS hdf4utils netcdfutils hdf5utils)

find_package(NetCDF REQUIRED)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES})

add_library (dfutils
  wrapper.c
  scale_sds.c
)

target_link_libraries(dfutils ${LIBS})
target_include_directories(dfutils PUBLIC ".")

install (TARGETS dfutils DESTINATION lib)
