cmake_minimum_required(VERSION 3.12)

set(CMAKE_Fortran_FLAGS "-std=legacy ${CMAKE_Fortran_FLAGS}")

set(LIBS nav genutils)

ADD_LIBRARY(swfnav
  actrng.f
  attcmp.f
  attdet.f
  cmpaxm.f
  conv_ins.f
  conv_soh.f
  conv_tilt.f
  earcnst.f
  earcomp.f
  earth.f
  earthsundist.f
  eartol.f
  eaxis.f
  ecef.f
  eci2ecef.f
  ellxfm.f
  euler.f
  fit3t.f
  fitrng.f
  fndflg.f
  get_node.f
  getorbit.f
  initnav.f
  inline.f
  intpos.f
  jul2yds.f
  leapchk.f
  matvec.f
  navtlm.f
  orb_interp.f
  orbcomp.f
  proctim2.f
  read_analog.f
  read_discrete.f
  read_double.f
  read_float.f
  read_gps.f
  read_long.f
  read_short.f
  readctl.f
  readqc.f
  runfit3t.f
  sc_att.f
  scpar.f
  sctime.f
  suncnst.f
  suncomp.f
  suntol.f
  sunvec.f
  sunwgt.f
  swfnav.f
  tiltcomp.f
  vmag.f
  xpose.f
  yd2ymd.f
  yds2jul.f
  ydsadd.f
  ymd2day.f
)

target_link_libraries(swfnav ${LIBS})
target_include_directories(swfnav PUBLIC ".")

install (TARGETS swfnav DESTINATION lib)
