cmake_minimum_required(VERSION 3.12)

#set_ocssw_policy()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(LIBS piutils)

find_package(HDF5 REQUIRED COMPONENTS C)
include_directories(${HDF5_INCLUDE_DIRS})
list(APPEND LIBS ${HDF5_LIBRARIES})

add_library (hdf5utils
  hdf5util.cpp
)

target_link_libraries(hdf5utils ${LIBS})
target_include_directories(hdf5utils PUBLIC ".")

install (TARGETS hdf5utils DESTINATION lib)
