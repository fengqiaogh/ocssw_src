cmake_minimum_required(VERSION 3.12)

find_package(NetCDF REQUIRED)
include_directories(${NETCDF_INCLUDE_DIRS})
set(LIBS ${NETCDF_LIBRARIES} netcdfutils)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
list(APPEND LIBS ${GSL_LIBRARIES})

add_executable(outernetDecode
  outernetDecode.c
  fec_seahawk.c
)

add_executable(l0gen_hawkeye
  l0gen_hawkeye.cpp
  hawkeyeUtil.cpp
)

add_executable(l1agen_hawkeye
  l1agen_hawkeye.cpp
  hawkeyeUtil.cpp
  HawkeyeDecode.c
)

add_executable(geolocate_hawkeye
  geolocate_hawkeye.cpp
  hawkeyeUtil.cpp
)


target_link_libraries(l0gen_hawkeye
  ${LIBS}
)

target_link_libraries(l1agen_hawkeye
  ${LIBS}
)

target_link_libraries(geolocate_hawkeye
  ${LIBS}
)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_EXE_LINKER_FLAGS "-Wl,-stack_size -Wl,0x20000000")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

install (TARGETS outernetDecode l0gen_hawkeye l1agen_hawkeye geolocate_hawkeye DESTINATION bin)
