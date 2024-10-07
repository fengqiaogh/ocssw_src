cmake_minimum_required(VERSION 3.12)

set(LIBS seawifs)

add_executable(navassess_seawifs
  main_navassess.f
  dmatch.f
  doublt.f
  getcat.f
  sortcl.f
  swapcl.f
  corect.f
  coline.f
  dunvec.f
  dangle.f
  id_drv.f
  triplt.f
  dirini.f
  identy.f
  mediqr.f
  select.f
  find_segs.f
  find_all_segs.f
  check_segs.f
  b28_lin.f
  find_islands.f
  get_l1a_open.f
  get_l1a_record.f
  get_l1a_close.f
)

target_link_libraries(navassess_seawifs ${LIBS})

set_target_properties(navassess_seawifs PROPERTIES LINKER_LANGUAGE Fortran)

install (TARGETS navassess_seawifs DESTINATION bin)

