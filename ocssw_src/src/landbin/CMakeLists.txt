cmake_minimum_required(VERSION 3.12)

set(LIBS hdf4utils genutils)

add_executable(landbin
  make_L3_v1.1.c
  calib_calibrate_l1a.c
  calib_get_cal.c
  calib_get_cal_misc.c
  chand.c
  compute_l1b.c
  csalbr.c
  get_attributes.c
  get_calib_sds.c
  get_navig_sds.c
  make_psurf.c
  proj_cproj.c
  proj_hamfor.c
  proj_haminv.c
  proj_molwfor.c
  proj_molwinv.c
  proj_report.c
  proj_robfor.c
  proj_robinv.c
  read_write.c
)

add_executable(landtimebin
  landtimebin.c
)

add_definitions(-DHDF41r3)

target_link_libraries(landbin ${LIBS})

target_link_libraries(landtimebin ${LIBS})

install (TARGETS landbin landtimebin DESTINATION bin)

