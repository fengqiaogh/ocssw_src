cmake_minimum_required(VERSION 3.12)

add_executable(l2extract
  main_l2extract.c
  l2extract_netcdf.c
  l2extract_input.cpp
  l2extract_3d_wv.cpp
)

include_directories (
  ${CMAKE_SOURCE_DIR}/src/l2gen
)

target_link_libraries(l2extract
  l2
)

install (TARGETS l2extract DESTINATION bin)

