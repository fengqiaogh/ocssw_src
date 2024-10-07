cmake_minimum_required(VERSION 3.12)

add_executable(ql02raw_seawifs
  ql02raw.f
)

set_target_properties(ql02raw_seawifs PROPERTIES LINKER_LANGUAGE Fortran)

install (TARGETS ql02raw_seawifs DESTINATION bin)

