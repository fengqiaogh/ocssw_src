cmake_minimum_required(VERSION 3.12)

set(LIBS nav hdf4utils timeutils)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
list(APPEND LIBS ${GSL_LIBRARIES})

ADD_LIBRARY(seawifs
  geonav.f
  geonav_lonlat.f
  get_cal_swf.c
  get_cal.c
  get_cal_misc.c
  calibrate_l1a.c
  l1a.c
  lac_st.c
  gac_st.c
  out_band_corr.c
  st_lt.c
  GetStationInfo.c
)

target_link_libraries(seawifs ${LIBS})
target_include_directories(seawifs PUBLIC ".")

install (TARGETS seawifs DESTINATION lib)
