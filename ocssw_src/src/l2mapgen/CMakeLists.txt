cmake_minimum_required(VERSION 3.12)

set(LIBS l2 hdf5utils)

find_package(GeoTIFF REQUIRED)
include_directories(${GeoTIFF_INCLUDE_DIRS})
list(APPEND LIBS ${GeoTIFF_LIBRARIES})

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIRS})
list(APPEND LIBS ${PNG_LIBRARIES})

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    find_package(X11 REQUIRED)
    include_directories(${X11_INCLUDE_DIR})
    list(APPEND LIBS ${X11_LIBRARIES})
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

add_executable(l2mapgen
  main_l2mapgen.c
  l2mapgen_input.c
  mipolyutil.c
)

target_link_libraries(l2mapgen ${LIBS})

install (TARGETS l2mapgen DESTINATION bin)

