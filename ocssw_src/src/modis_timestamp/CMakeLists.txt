cmake_minimum_required(VERSION 3.12)

add_executable(modis_timestamp
  modis_timestamp.c
)

install (TARGETS modis_timestamp DESTINATION bin)

