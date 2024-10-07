cmake_minimum_required(VERSION 3.12)

find_package(NetCDF REQUIRED)
include_directories(${NETCDF_INCLUDE_DIRS})
set(LIBS ${NETCDF_LIBRARIES})
  
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSTDC_HEADERS -DHAVE_STRING_H")
add_executable(nccmp
  opt.c
  strlist.c
  xmalloc.c
  ncinfo.c
  nccmp_user_type.c
  nccmp.cpp
)

target_link_libraries(nccmp ${LIBS})

install (TARGETS nccmp DESTINATION bin)
