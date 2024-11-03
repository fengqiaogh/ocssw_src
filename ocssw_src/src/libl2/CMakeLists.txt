cmake_minimum_required(VERSION 3.12)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if (BUILD_HISTORICAL)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D GEONAV")
  set(seawifs_lib seawifs)
endif()

set(LIBS l1 ${seawifs_lib} dfutils)

#find_package(NetCDF REQUIRED)
#include_directories(${NETCDF_INCLUDE_DIRS})
#list(APPEND LIBS ${NETCDF_LIBRARIES})

ADD_LIBRARY(l2
  readL2scan.c
  setupflags.c
  get_product_table.c
  expand3D.cpp
  get_geospatial.cpp
  find_variable.c
  find_variable.cpp
)

target_link_libraries(l2 ${LIBS})
target_include_directories(l2 PUBLIC ".")

install (TARGETS l2 DESTINATION lib)
