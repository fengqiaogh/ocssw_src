cmake_minimum_required(VERSION 3.12)

add_executable(fmt_check
  fmt_check.c
  fmt_read.c
  fmt_rd_attr.c
  fmt_rd_dim.c
  fmt_rd_sds.c
  fmt_rd_ras.c
  fmt_rd_l3vg.c
  get_attr.c
  attr_disp.c
  chk_sds.c
  hdf_ras_chk.c
  ck_v_l3.c
  l3_get_org.c
  s_parse.c
  get_line.c
  var_decode.c
  chk_str.c
)

add_definitions(-D__GLX_MOTIF -DSYSV -D_NO_PROTO -DNOCHKHDF -DFCK_8LEN=5000)

target_link_libraries(fmt_check
  hdf4utils
  hdf_dio
)

install (TARGETS fmt_check DESTINATION bin)

