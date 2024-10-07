cmake_minimum_required(VERSION 3.12)

add_executable(gpsextract_seawifs
  main_gpsextract.c
)

target_link_libraries(gpsextract_seawifs
  genutils 
  timeutils
)

install (TARGETS gpsextract_seawifs DESTINATION bin)

