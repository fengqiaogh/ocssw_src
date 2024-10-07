cmake_minimum_required(VERSION 3.12)

find_package(HDF4 REQUIRED)
include_directories(${HDF4_INCLUDE_DIRS})
set(LIBS ${HDF4_LIBRARIES})

ADD_LIBRARY(h4metadata
  granule_extract_metadata.c
)

target_link_libraries(h4metadata ${LIBS})
target_include_directories(h4metadata PUBLIC ".")

install (TARGETS h4metadata DESTINATION lib)
