cmake_minimum_required(VERSION 3.12)

#find_package(HDF4 REQUIRED)
#include_directories(${HDF4_INCLUDE_DIRS})
#set(LIBS ${HDF4_LIBRARIES})

find_package(NetCDF REQUIRED COMPONENTS C CXX)
include_directories(${NETCDF_INCLUDE_DIRS})
set(LIBS ${NETCDF_LIBRARIES})

list(APPEND LIBS l1)
list(APPEND LIBS genutils)

add_definitions(-DMODIS_OD -DMODIS_INST -DRETRIEVE -DCT_1KM -DUSE_GDAS -DGEOS5_SFC)
set (CMAKE_Fortran_FLAGS "-Wno-tabs -ffree-line-length-0 ${CMAKE_Fortran_FLAGS}")

ADD_LIBRARY(cloud
  ch_xfr.f90
  mod06_run_settings.f90
  names.f90
  cld_fchk.f90
  mng_ms.cpp
  dim_mgr.cpp
  ch_cld_sci.f90
  ncio.c
  Driver_MOD_PR06OD.f90
  modis_cloudstructure.f90
  GeneralAuxType.f90
  modis_sciencestructure.f90
  nonscience_parameters.f90
  science_parameters.f90
  global_model_grids.f90
  ct_core_arrays.f90
  core_arrays.f90
  libraryinterpolates.f90
  specific_other.f90
  specific_ancillary.f90
  modis_numerical_module.f90
  MOD06AlbedoEcoModule.f90
  modis_albedo.f90
  ancillary_module.f90
  libraryarrays.f90
  modis_io_module.f90
  interpolate_libraries.f90
  modis_frontend_module.f90
  get_retrieval_uncertainty.f90
  spline_module.f90
  planck_functions.f90
  retrieval_solution_logic.f90
  retrieval_prep_logic.f90
  corescience_module.f90
  multi_layer_clouds.f90
  clear_sky_restoral.f90
  FASCODE_routines.f90
  pfaast.f90
  rtm_support.f90
  retrieval_irw.f90
  atmospheric_correction.f90
  cloud_phase.f90
  general_science_module.f90
  modis_science_module.f90
  mod_pr06od.f90
  local_message_handler.f90
  MODIS_SMF_SETDYNAMICMSG_V2.1.f
  nnr_001.mydl_Tau.f
  nnr_001.mydo_Tau.f
  nnr_001.modo_Tau.f
  nnr_001.modl_Tau.f
  profile_management.c
  extem101_64.f
  remove_edge_scenes.f90
  ffnet.f
  set_quality_data_module.f90
  get_cld_tbl.c
  cld_tbl_names.f90
)

target_link_libraries(cloud ${LIBS})
target_include_directories(cloud PUBLIC ".")

install (TARGETS cloud DESTINATION lib)
