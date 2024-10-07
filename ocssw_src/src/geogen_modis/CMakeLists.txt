cmake_minimum_required(VERSION 3.12)

find_package(HDF4 REQUIRED)
include_directories(${HDF4_INCLUDE_DIRS})
set(LIBS ${HDF4_LIBRARIES})

find_package(SDST REQUIRED)
include_directories(${SDST_INCLUDE_DIRS})
list(APPEND LIBS ${SDST_LIBRARIES})

find_package(Mapi REQUIRED)
include_directories(${Mapi_INCLUDE_DIRS})
list(APPEND LIBS ${Mapi_LIBRARIES})

list(APPEND LIBS h4metadata)

add_executable(geogen_modis
  GEO_location_main.c
  GEO_locate_one_granule.c
  GEO_locate_one_scan.c
  GEO_get_inst_mirr_normal.c
  GEO_get_sample_time.c
  GEO_get_view_vec.c
  GEO_interp_mirr_ang.c
  GEO_interp_mirr_enc.c
  GEO_mat_vec_mul3.c
  GEO_poly_fit.c
  GEO_vec_unit3.c
  GEO_vec_length3.c
  GEO_vec_mul3.c
  GEO_poly_coef1.c
  GEO_vec_prod3.c
  GEO_DEM.c
  GEO_get_geoid.c
  GEO_landsea_mask.c
  GEO_terrain_correct.c
  GEO_check_ea_headers.c
  GEO_earth_location.c
  GEO_ellip_position.c
  GEO_ephem_attit.c
  GEO_get_ephatt_inputs.c
  GEO_get_T_inst2ecr.c
  GEO_aggregate.c
  GEO_hires.c
  GEO_solar_and_lunar_vectors.c
  GEO_interp_ECR.c
  GEO_maneuver.c
  GEO_prepare_l1a_data.c
  GEO_prepare_mirr_data.c
  GEO_read_param_file.c
  GEO_initialize_product.c
  GEO_derived_products.c
  GEO_update_L1A_metadata.c
  GEO_get_GRing_points.c
  GEO_get_bounding_coords.c
  GEO_create_swath.c
  GEO_get_utcpole_metadata.c
  GEO_get_version_metadata.c
  GEO_read_L1Aspecific_metadata.c
  GEO_read_L1Ascan_metadata.c
  GEO_read_L1Apacket_data.c
  GEO_read_L1AECS_metadata.c
  GEO_read_L1Atemp_data.c
  GEO_validate_derived_products.c
  GEO_validate_earth_location.c
  GEO_abs_limit_check.c
  GEO_del_limit_check.c
  GEO_find_next_flag.c
  GEO_write_ECS_metadata.c
  GEO_write_geospecific_metadata.c
  GEO_write_one_scan.c
  GEO_write_scan_data.c
  GEO_write_granule_metadata.c
  GEO_write_parameters.c
  GEO_write_scan_metadata.c
  GEO_write_input_metadata.c
  imsl_d_lin_sol_gen.c
  imsl_d_spline_interp.c
  imsl_d_spline_value.c
  imsl_error.c
  )

add_definitions(-w -DNOIMSL '-DMAKEFILE_REVISION="1.1"' -DDIRECT_BROADCAST)

target_link_libraries(geogen_modis
  ${LIBS}
)

install (TARGETS geogen_modis DESTINATION bin)
