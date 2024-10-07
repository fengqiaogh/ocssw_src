cmake_minimum_required(VERSION 3.12)

add_executable(neb2raw_seawifs
  neb2raw.c
)

install (TARGETS neb2raw_seawifs DESTINATION bin)

