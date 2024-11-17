cmake_minimum_required(VERSION 3.12)

set(CMAKE_Fortran_FLAGS "-std=legacy ${CMAKE_Fortran_FLAGS}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${GSL_LIBRARIES})
  
list(APPEND LIBS genutils m)
# Required packages
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWINPATH")
endif()


if (BUILD_HISTORICAL)
  set(src_nav_cns cdata.f)
endif()

ADD_LIBRARY(nav
  ${src_nav_cns}
  get_zenaz.c
  compute_alpha.c
  compute_dist.c
  nav.c
  math_utils.c
  crossp.c
  sunangs.c
  sun2000.c
  gha2000.c
  ephparms.c
  nutate.c
  l_sun.c
  ocorient.c
  esdist.c
  getenv.cpp
  ftn_interface.c
)

target_link_libraries(nav ${LIBS})
target_include_directories(nav PUBLIC ".")

install (TARGETS nav DESTINATION lib)
