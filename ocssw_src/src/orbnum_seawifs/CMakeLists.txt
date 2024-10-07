cmake_minimum_required(VERSION 3.12)

set(CMAKE_Fortran_FLAGS "-std=legacy ${CMAKE_Fortran_FLAGS}")

set(LIBS seawifs)

add_executable(orbnum_seawifs
   main_orbnum.f
   getorbit.f
)

target_link_libraries(orbnum_seawifs ${LIBS})

set_target_properties(orbnum_seawifs PROPERTIES LINKER_LANGUAGE Fortran)

install (TARGETS orbnum_seawifs DESTINATION bin)

