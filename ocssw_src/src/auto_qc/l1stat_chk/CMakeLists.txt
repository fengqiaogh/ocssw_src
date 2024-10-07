cmake_minimum_required(VERSION 3.12)

find_package(HDF4 REQUIRED)
include_directories(${HDF4_INCLUDE_DIRS})
set(LIBS ${HDF4_LIBRARIES})

add_executable(l1stat_chk
  l1stat_chk.c
  geovex.f
  ck_trng.c
  anal_noise.c
  rpt_noise.c
  stat_chk_utl.c
  chk_inst_ana.c
  chk_gainv.c
  chk_tdiv.c
)

add_definitions(-DPROTOTYPE)

target_link_libraries(l1stat_chk
  ${LIBS}
  seawifs
  hdf_dio
  )

install (TARGETS l1stat_chk DESTINATION bin)

