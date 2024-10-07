cmake_minimum_required(VERSION 3.12)

#set_ocssw_policy()

set(LIBS "")

find_package(HDF5 REQUIRED COMPONENTS C HL)
include_directories(${HDF5_INCLUDE_DIRS})
list(APPEND LIBS ${HDF5_LIBRARIES})

find_package(NetCDF REQUIRED COMPONENTS CXX F77)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES})

list(APPEND LIBS piutils dfutils genutils)

list(APPEND LIBS meris)

list(APPEND LIBS h5io)

list(APPEND LIBS bilio)

add_library(filetype filetype.c filetypeXml.cpp)
target_link_libraries(filetype ${LIBS})
target_include_directories(filetype PUBLIC ".")

install (TARGETS filetype DESTINATION lib)
