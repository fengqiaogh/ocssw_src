cmake_minimum_required(VERSION 3.12)

set(LIBS bilio netcdfutils)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
list(APPEND LIBS ${GSL_LIBRARIES})

find_package(PROJ REQUIRED CONFIG)
list(APPEND LIBS PROJ::proj)

add_executable(avirisbil2oci
  avirisbil2oci.cpp
)

target_link_libraries(avirisbil2oci ${LIBS})

install (TARGETS avirisbil2oci DESTINATION bin)
