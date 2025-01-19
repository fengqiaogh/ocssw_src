cmake_minimum_required(VERSION 3.12)

#
# Make a version file containing the current version from git.
#
set (VERSION_MAJOR 9)
set (VERSION_MINOR 10)
set (VERSION_PATCH 0)
include(GetOcsswVersion)

#####version 3 of Atrem - no fortran libraries, just C code
# TODO: get rid of numerical.c
set (ATREM_C_CODE numerical.c get_atrem_corl1v3.c rdatreminfo.c)
set (ATREM_F_CODE )
#####version 2 of Atrem - mix with fortran libraries and C code
#set (ATREM_C_CODE numerical.c  get_atrem_corl1v2.c rdatreminfo.c)
#set (ATREM_F_CODE atrem_app_refl_plus_gas_removal_l2v2.f90 atrem_bndprms.f atrem_tpvmr_init.f)
#####version 1 of Atrem - need fortran libraries and code
#set (ATREM_C_CODE numerical.c get_atrem_corl1.c rdatreminfo.c)
#set (ATREM_F_CODE atrem_app_refl_plus_gas_removal_l2.f90 atrem_bndprms.f atrem_tpvmr_init.f)

set(LIBS l2 l1 nav libterrain)

find_package(Levmar REQUIRED)
include_directories(${Levmar_INCLUDE_DIRS})
list(APPEND LIBS ${Levmar_LIBRARIES})

find_package(LAPACK REQUIRED)
include_directories(${LAPACK_INCLUDE_DIRS})
list(APPEND LIBS ${LAPACK_LIBRARIES})

find_package(PugiXML REQUIRED)
include_directories(${PUGIXML_INCLUDE_DIRS})
list(APPEND LIBS ${PUGIXML_LIBRARIES})

find_package(Fann REQUIRED)
include_directories(${Fann_INCLUDE_DIRS})
list(APPEND LIBS ${Fann_LIBRARIES})

find_package(NetCDF REQUIRED COMPONENTS CXX F77)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES})

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIRS})


if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    find_package(X11 REQUIRED)
    include_directories(${X11_INCLUDE_DIR})
    list(APPEND LIBS ${X11_LIBRARIES})
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

list(APPEND LIBS timeutils)
list(APPEND LIBS h5io)
list(APPEND LIBS filetype)
list(APPEND LIBS bilio)

set(L2GEN_COMMON_FILES
  filter.c
  viirs_pxcvt.c
  setflags_l2.c
  wavelength_3d.cpp
  geo_region.cpp
# fortran files to get rid of
  lenstr.cpp
)

if (BUILD_CLOUD)
  set(cloud_src get_cmp.c)
  set(cloud_lib cloud)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DBUILD_CLOUD")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBUILD_CLOUD")
endif()
if(BUILD_HISTORICAL)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DBUILD_HISTORICAL")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBUILD_HISTORICAL")
endif(BUILD_HISTORICAL)

set(L2GEN_PRODUCT_FILES
  prodgen.c
  prodlist.c
  get_l2prod_index.cpp
  get_chl.c
  get_es.c
  get_habs.c
  get_Kd.c
  get_ndvi.c
  get_smoke.c
  get_depth_classification.c
  get_poc.c
  get_poc_hybrid.c
  get_par.c
  get_toa_refl.c
  get_tricho.c
  get_npp.c
  get_avw.c
  get_Cphyt.c
  get_Cpicophyt.c
  photic_depth.c
  cdom_morel.c
  calc_par.c
  cdom_mannino.c
  par_utils.c
  ipar.c
  fluorescence.c
  get_rhos.c
  calcite.c
  amoeba.c
  aph.c
  gsm.c
  giop.c
  mgiop.c
  carder.c
  las_iop.c
  pml.c
  get_pml.c
  pml_iop_calculate.c
  pml_iop_config.c
  pml_iop_tables.c
  qaa.c
  get_qaa.c
  flags_iop.c
  niwa_iop.c
  get_niwa_iop.c
  sst_cloud_mask.cpp # separate the cloud mask rom the NLSST algorithm
  sst.c
  sst_adt.cpp
  sst_cloud_mask_utils.cpp
  sst_dsdi.c
  mumm.c
  get_bsi.c
  swim.c
  get_pft_hirata.c
  get_pft_uitz.c
  get_psd_ksm.c
  get_nitrate.c
  ${cloud_src}
  get_ctht.c
  int_4d.c
  axis_interp.c
  ams_oe_inversion.c
  get_bpar.c
  get_uncertainty.c
  get_ndwi.c
  get_ndii.c
  get_cci.c
  get_ndsi.c
  get_hyper_vi.cpp
  vegetation_indices.c
  get_sdp.c
  get_gpig.c
  trustopt.c
#  myprod.c
)

set(L2GEN_ANCILLARY_FILES
  gas_trans.c
  seawater.c
  seawater_get.c
  water_spectra.c
  setanc.c
  acq_sfc_albedo.c
  b128_msk_get.c
  b128_wd_bit.c
  b128_box_num.c
  read_mask.c
  get_height.c
  ice_mask.c
  get_ice_frac.c
  sstref.c
  sssref.c
  elev.c
  anc_acq.c
  met_cvt.c
  ncio.c
  get_mld.cpp
  get_zno3.c
  smi_climatology.c
)
# # trying out avx
# set_source_files_properties(sst_cloud_mask.cpp PROPERTIES COMPILE_FLAGS "-g -mavx2")
add_executable(l2gen
  main_l2gen.c
  msl12_input.c
  getl1rec.c
  loadl1.c
  read_pixel_anc_file.cpp
  cpl1rec.c
  target_io.c
  bioOptBandShift.c
  l2_generic.c
  alloc_l2.c
  alloc_target.c
  alloc_vcal.c
  convl12.c
  convl21.c
  atmocor1.c
  atmocor1_land.c
  whitecaps.c
  rayleigh.c
  polcor_hawkeye.cpp
  polcor_oci.cpp
  polcor.c
  water_vapor.c
  airmass.c
  aerosol.c
  atmocor2.c
  glint.c
  get_rhown_nir.c
  nlw_outband.c
  brdf.c
  bin_climatology.c
  aer_io.c
  alloc_aer.c
  vcal.c
  dtranbrdf.cpp
  turbid.c
  init_l2.c
  convert_band.c
  fuzzy_func_v3.c
  gammln.c
  gcf.c
  gser.c
  nrutil.c
  sprsax.c
  sprsin.c
  covariance_inversion.c
  lubksb.c
  ludcmp.c
  owt.c
  get_owmc.c
  lonlat2pixline.cpp
  virtual_constellation.c
  bioOptBandShift.c
  read_l3bin.cpp
  ${ATREM_C_CODE}
  ${ATREM_F_CODE}
  raman.c
  ${L2GEN_ANCILLARY_FILES}
  ${L2GEN_PRODUCT_FILES}
  ${L2GEN_COMMON_FILES}
)

add_executable(l3gen
  main_l3gen.cpp
  msl12_input.c
  getl1rec.c
  loadl1.c
  read_pixel_anc_file.cpp
  cpl1rec.c
  target_io.c
  bioOptBandShift.c
  l2_generic.c
  alloc_l2.c
  alloc_target.c
  alloc_vcal.c
  convl12.c
  convl21.c
  atmocor1.c
  atmocor1_land.c
  whitecaps.c
  rayleigh.c
  polcor_hawkeye.cpp
  polcor_oci.cpp
  polcor.c
  water_vapor.c
  airmass.c
  aerosol.c
  atmocor2.c
  glint.c
  get_rhown_nir.c
  nlw_outband.c
  brdf.c
  bin_climatology.c
  aer_io.c
  alloc_aer.c
  vcal.c
  dtranbrdf.cpp
  turbid.c
  init_l2.c
  convert_band.c
  fuzzy_func_v3.c
  gammln.c
  gcf.c
  gser.c
  nrutil.c
  sprsax.c
  sprsin.c
  covariance_inversion.c
  lubksb.c
  ludcmp.c
  owt.c
  get_owmc.c
  virtual_constellation.c
  bioOptBandShift.c
  read_l3bin.cpp
  ${ATREM_C_CODE}
  ${ATREM_F_CODE}
  raman.c
  ${L2GEN_ANCILLARY_FILES}
  ${L2GEN_PRODUCT_FILES}
  ${L2GEN_COMMON_FILES}
)


add_executable(vcalmerge
  main_vcalmerge.c
  calfile_utils.c
  msl12_input.c
  l2_generic.c
  loadl1.c
  read_pixel_anc_file.cpp
  atmocor1.c
  atmocor1_land.c
  whitecaps.c
  rayleigh.c
  polcor_hawkeye.cpp
  polcor_oci.cpp
  polcor.c
  water_vapor.c
  bioOptBandShift.c
  convl12.c
  convl21.c
  alloc_l2.c
  getl1rec.c
  alloc_aer.c
  alloc_target.c
  init_l2.c
  alloc_vcal.c
  vcal.c
  virtual_constellation.c
  bioOptBandShift.c
  lubksb.c
  ludcmp.c
  covariance_inversion.c
  nrutil.c
  sprsax.c
  sprsin.c
  gammln.c
  gcf.c
  gser.c
  fuzzy_func_v3.c
  owt.c
  get_owmc.c
  lonlat2pixline.cpp
  atmocor2.c
  brdf.c
  cpl1rec.c
  convert_band.c
  flags_iop.c
  dtranbrdf.cpp
  turbid.c
  glint.c
  get_rhown_nir.c
  aerosol.c
  nlw_outband.c
  airmass.c
  bin_climatology.c
  ${ATREM_C_CODE}
  ${ATREM_F_CODE}
  raman.c
  ${L2GEN_ANCILLARY_FILES}
  ${L2GEN_PRODUCT_FILES}
  ${L2GEN_COMMON_FILES}
)


add_executable(l1det2det
  main_l1det2det.c
  calfile_utils.c
  msl12_input.c
  l2_generic.c
  loadl1.c
  read_pixel_anc_file.cpp
  atmocor1.c
  atmocor1_land.c
  whitecaps.c
  rayleigh.c
  polcor_hawkeye.cpp
  polcor_oci.cpp
  polcor.c
  water_vapor.c
  bioOptBandShift.c
  convl12.c
  convl21.c
  alloc_l2.c
  read9km_mask.c
  getl1rec.c
  alloc_aer.c
  alloc_target.c
  init_l2.c
  alloc_vcal.c
  vcal.c
  virtual_constellation.c
  bioOptBandShift.c
  lubksb.c
  ludcmp.c
  covariance_inversion.c
  nrutil.c
  sprsax.c
  sprsin.c
  gammln.c
  gcf.c
  gser.c
  fuzzy_func_v3.c
  owt.c
  get_owmc.c
  lonlat2pixline.cpp
  atmocor2.c
  brdf.c
  cpl1rec.c
  convert_band.c
  flags_iop.c
  dtranbrdf.cpp
  turbid.c
  glint.c
  get_rhown_nir.c
  aerosol.c
  nlw_outband.c
  airmass.c
  bin_climatology.c
  ${ATREM_C_CODE}
  ${ATREM_F_CODE}
  raman.c
  ${L2GEN_ANCILLARY_FILES}
  ${L2GEN_PRODUCT_FILES}
  ${L2GEN_COMMON_FILES}
  )

add_executable(l1brsgen
  main_l1brsgen.c
  l1_imgscale.c
  msl12_input.c
  loadl1.c
  read_pixel_anc_file.cpp
  atmocor1.c
  atmocor1_land.c
  whitecaps.c
  rayleigh.c
  polcor_hawkeye.cpp
  polcor_oci.cpp
  polcor.c
  water_vapor.c
  get_rhos.c
  ${ATREM_C_CODE}
  ${ATREM_F_CODE}
  ${L2GEN_ANCILLARY_FILES}
  ${L2GEN_COMMON_FILES}
  )

add_executable(l1mapgen
  main_l1mapgen.c
  l1_imgscale.c
  msl12_input.c
  loadl1.c
  read_pixel_anc_file.cpp
  atmocor1.c
  atmocor1_land.c
  whitecaps.c
  rayleigh.c
  polcor_hawkeye.cpp
  polcor_oci.cpp
  polcor.c
  water_vapor.c
  mipolyutil.c
  get_rhos.c
  ${ATREM_C_CODE}
  ${ATREM_F_CODE}
  ${L2GEN_ANCILLARY_FILES}
  ${L2GEN_COMMON_FILES}
  )

add_executable(l1bgen_generic
  main_l1bgen.c
  msl12_input.c
  ${L2GEN_COMMON_FILES}
)

add_executable(l1info
  main_l1info.cpp
  msl12_input.c
  gringHelper.cpp
  ${L2GEN_COMMON_FILES}
  )

add_executable(lonlat2pixline
  main_lonlat2pixline.c
  msl12_input.c
  lonlat2pixline.cpp
  ${L2GEN_COMMON_FILES}
  )

add_executable(l2binmatch
  main_l2binmatch.cpp
  l2binmatch_input.cpp
  get_l2prod_index.cpp
  calfile_utils.c
  read_mask.c
  b128_msk_get.c
  b128_wd_bit.c
  b128_box_num.c
 )

target_link_libraries(l2gen
  anc
  bin
  ${cloud_lib}
  ${LIBS}
  )

target_link_libraries(l3gen
  anc
  bin
  ${cloud_lib}
  ${LIBS}
  )

target_link_libraries(l1bgen_generic
  ${LIBS}
  )

target_link_libraries(l1brsgen
  anc
  imageutils
  ${LIBS}
  )

target_link_libraries(l1mapgen
  anc
  ${PNG_LIBRARIES}
  ${LIBS}
  )

target_link_libraries(l1info
  ${LIBS}
  )

target_link_libraries(lonlat2pixline
  ${LIBS}
  )

target_link_libraries(vcalmerge
  anc
  ${cloud_lib}
  ${LIBS}
  )

target_link_libraries(l2binmatch
  bin
  ${LIBS}
  )

target_link_libraries(l1det2det
  anc
  ${cloud_lib}
  ${LIBS}
  )

#set (CMAKE_C_FLAGS_RELEASE "-O1")
# Use -Wpadded -Wpacked to warn about align-double. Use  '#pragma GCC diagnostic ignored "-Wpadded"' in file if you don't want to see messages
# -malign-double needed for accessing fortran common blocks from C structs

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm64")
    #set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -malign-double -mpreferred-stack-boundary=8")
    #set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -malign-double -mpreferred-stack-boundary=8")
  else()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -malign-double -mpreferred-stack-boundary=8")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -malign-double -mpreferred-stack-boundary=8")
  endif()
endif()

install (TARGETS l2gen l3gen l1bgen_generic l1brsgen l1mapgen l1info lonlat2pixline vcalmerge l2binmatch l1det2det DESTINATION bin)
