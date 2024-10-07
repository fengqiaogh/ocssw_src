cmake_minimum_required(VERSION 3.12)

set(LIBS l2)

add_executable(viirs_sim_sdr
  viirs_sim_sdr.c
  viirs_sim_input.c
  rd_sim_init.c
  rd_geo_init.c
  init_sdr.c
  wr_attr_seq.c
  day2mday.c
  time_utl.c
  rd_sdr_scan.c
  rd_geo_scan.c
  gen_const_rad_scn.c
  wr_sdr_scan.c
  wr_geo_scan.c
  wr_bnd_scan.c
  fin_sdr.c
  gen_sdr_fname.c
  rd_rhos_scan.c
  rhos_to_lt.c
  mod_artifact.c
  viirs_oxt.c
  scan_cvt.c
  bbt_2_rad.c
  viirs_cal.c
  viirs_decal.c
  vset_cal.c
  vir_xf_scan.c
  hio_open.c
  hio_r_sds.c
  hio_close.c
  bnd_ix_2_sen_info.c
  viirs_ext.c
  viirs_noise.c
  viirs_straylt.c
)

include_directories(../l2gen)
include_directories(../h5io)
list(APPEND LIBS h5io)

target_link_libraries(viirs_sim_sdr ${LIBS})

install(TARGETS viirs_sim_sdr DESTINATION bin)

