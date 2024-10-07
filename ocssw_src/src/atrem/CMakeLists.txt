cmake_minimum_required(VERSION 3.12)

# This is based on the Makefile_gfortran_atrem_app_refl_plus_gas_removal_2013
# provided by Bo-Cai Gao

add_executable(atrem
  atrem_app_refl_plus_gas_removal_2013.f
  cubeio.f90
  bndprms.f
  solar_irr_PC.f
  tpvmr_init.f
  6sm1.f
)
add_executable(get_atrem_cor
  get_atrem_cor3.c
  cubeio.f90
  bndprms.f
  solar_irr_PC.f
  tpvmr_init.f
)

target_link_libraries(get_atrem_cor
  hdfutils
  fann
  )

set_target_properties(atrem PROPERTIES LINKER_LANGUAGE Fortran)

set (CMAKE_Fortran_FLAGS "-DLINUX -Dlinux -ffixed-line-length-none -fno-automatic -fno-range-check -cpp")
set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops")
set (CMAKE_Fortran_FLAGS_DEBUG "-g3 -gdwarf-2 -Wall -Wno-unused-variable -Wno-unused-label -Wno-unused-dummy-argument -Wno-conversion")


install (TARGETS atrem DESTINATION bin)
