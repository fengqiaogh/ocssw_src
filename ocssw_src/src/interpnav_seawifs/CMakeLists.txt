cmake_minimum_required(VERSION 3.12)

add_executable(interpnav_seawifs
  main_interpnav.f
  readctl.f
  inline.f
  get_xfm.f
  interp_att.f
  ellxfm.f
  euler.f
  xpose.f
  eaxis.f
  get_l1a_open.f
  get_l1a_data.f
  put_l1a_data.f
  put_l1a_metadata.f
  checkorb.f
  mfsdff.f
)

#### do not link with swfnav.
target_link_libraries(interpnav_seawifs seawifs)

set_target_properties(interpnav_seawifs PROPERTIES LINKER_LANGUAGE Fortran)

install (TARGETS interpnav_seawifs DESTINATION bin)

