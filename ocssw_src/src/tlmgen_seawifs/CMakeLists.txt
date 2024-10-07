cmake_minimum_required(VERSION 3.12)

add_executable(tlmgen_seawifs
  main_tlmgen.f
)

set_target_properties(tlmgen_seawifs PROPERTIES LINKER_LANGUAGE Fortran)

install (TARGETS tlmgen_seawifs DESTINATION bin)

