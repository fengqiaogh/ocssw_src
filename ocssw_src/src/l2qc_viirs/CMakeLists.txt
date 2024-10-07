cmake_minimum_required(VERSION 3.12)

find_package(NetCDF REQUIRED)
include_directories(${NETCDF_INCLUDE_DIRS})
set(LIBS ${NETCDF_LIBRARIES})

list(APPEND LIBS genutils piutils)

add_executable(l2qc_viirs
    l2qc_viirs.cpp
    l2qc_viirs_input.cpp
)

target_link_libraries(l2qc_viirs ${LIBS})

install (TARGETS l2qc_viirs DESTINATION bin)

