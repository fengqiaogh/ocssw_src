cmake_minimum_required(VERSION 3.12)

set(LIBS netcdfutils hdf5utils)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
list(APPEND LIBS ${GSL_LIBRARIES})
  
find_package(PROJ REQUIRED CONFIG)
list(APPEND LIBS PROJ::proj)
  
add_executable(prismbil2oci
  prismbil2oci.cpp
)

target_link_libraries(prismbil2oci
  bilio
  ${LIBS}
)

install (TARGETS prismbil2oci DESTINATION bin)
