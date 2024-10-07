cmake_minimum_required(VERSION 3.12)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(LIBS seawifs hdf4utils)

ADD_LIBRARY(osmi
  get_cal_osmi.c
  get_cal_misc_osmi.c
  calibrate_l1a_osmi.c
  calday.c
  obliq.c
  sunpos.c
  ctogd.c
  gmha.c
  tconv.c
  ctotc.c
  julian.c
  time-utils.c
  eanom.c
  locate.c
)

target_link_libraries(osmi ${LIBS})
target_include_directories(osmi PUBLIC ".")

install (TARGETS osmi DESTINATION lib)
