cmake_minimum_required(VERSION 3.12)

#set(LIBS "")
set(LIBS genutils timeutils cgal_interp)

find_package(NetCDF REQUIRED COMPONENTS C CXX)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
list(APPEND LIBS ${GSL_LIBRARIES})

add_executable(ancgen
  ancgen.cpp
  copy_var_utils.cpp
)

add_executable(extractSPEX
  extractSPEX.cpp
  copy_var_utils.cpp
)

target_link_libraries(ancgen ${LIBS})
target_link_libraries(extractSPEX ${LIBS})

install (TARGETS ancgen extractSPEX DESTINATION bin)

