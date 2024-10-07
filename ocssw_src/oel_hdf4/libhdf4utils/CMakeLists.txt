cmake_minimum_required(VERSION 3.12)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(LIBS piutils)

find_package(HDF4 REQUIRED)
include_directories(${HDF4_INCLUDE_DIRS})
list(APPEND LIBS ${HDF4_LIBRARIES})

add_library (hdf4utils
  hdf_utils.c
  hdf4_utils.c
)

target_link_libraries(hdf4utils ${LIBS})
target_include_directories(hdf4utils PUBLIC "." ${HDF4_INCLUDE_DIRS})

install (TARGETS hdf4utils DESTINATION lib)
