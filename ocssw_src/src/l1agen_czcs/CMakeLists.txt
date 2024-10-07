cmake_minimum_required(VERSION 3.12)

add_executable(l1agen_czcs
  l1czcs.c
  read_crtt.c
  czcs_l1_write.c
  hdr_2_gattr.c
  czcs_ctl_pt.c
  lladjust.f
  cz_ll_upd.c
  wrt_czcs_sla.c
  time_str.c
  wrt_czcs_qual.c
  cz_clean.c
  cz_dat_alloc.c
  cz_sd_set.c
  set_czcs_ctl_data.c
  create_sds.c
  set_dim_names.c
  fill_orb_dat.c
  time_utl.c
  day_to_ofile.c
  rd_smmr_orb.c
  lonlat.f
  hdfio.c
  asap_int2.f
)

add_executable(l1amerge_czcs
  czl1merge.c
  cz_dat_alloc.c
  cz_l1_read.c
  cz_mov_scn.c
  cztimqual.c
  fill_mstr.c
  hdfio.c
  olap_resolve.c
  czcs_l1_write.c
  wrt_czcs_qual.c
  read_file_list.c
  set_czcs_ctl_data.c
  wrt_czcs_sla.c
  create_sds.c
  set_dim_names.c
  cz_meta_adj.c
  cz_ll_upd.c
)

add_definitions(-D__DEC_MACHINE)

target_link_libraries(l1agen_czcs
  czcs
  seawifs
)

target_link_libraries(l1amerge_czcs
  czcs
  seawifs
)

install (TARGETS l1agen_czcs l1amerge_czcs DESTINATION bin)
