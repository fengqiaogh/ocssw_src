cmake_minimum_required(VERSION 3.12)

project(UTM)

add_executable(get_utm_zone
  get_utm_zone.cpp
  LatLong-UTMconversion.cpp
)

target_link_libraries(get_utm_zone
  m
)


install (TARGETS get_utm_zone DESTINATION bin)

