cmake_minimum_required(VERSION 3.12)

set(CMAKE_Fortran_FLAGS "-std=legacy ${CMAKE_Fortran_FLAGS}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${GSL_LIBRARIES})
  
list(APPEND LIBS genutils m)

ADD_LIBRARY(nav
  cdata.f
  crossp.f
  esdist.f
  filenv.f
  gha2000.f
  l_sun.f
  ocorient.f
  sunangs.f
  get_zenaz.c
  compute_alpha.c
  compute_dist.c
  get_ut1.f
  asap_rots.f
  asap_rot_int.f
  ymdhms2jul.f
  jd.f
  asaps.f
  add_elements.f
  get_elements.f
  put_elements.f
  pderiv.f
  fitgps.f
  vec2mean.f
  jdate.f
  matvec2.f
  jddate.f
  julian.f
  kepler.f
  kozsak2.f
  coord.f
  setsun.f
  setthd.f
  rk78.f
  rk78cn.f
  pout.f
  der.f
  mediqr.f
  invert.f
  eqnox.f
  delm.f
  dens76.f
  legend.f
  angles.f
  ephem.f
  xthird.f
  select.f
  anomly.f
)

target_link_libraries(nav ${LIBS})
target_include_directories(nav PUBLIC ".")

install (TARGETS nav DESTINATION lib)
