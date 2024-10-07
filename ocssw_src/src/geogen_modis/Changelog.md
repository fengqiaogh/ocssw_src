
# geogen_modis Changelog

## 6.2.1 - 2019-06-18

 - moved large static arrays off the stack

### Source Changed

 * CMakeLists.txt
 * GEO_landsea_mask.c
 * GEO_earth_location.c
 * GEO_read_L1Ascan_metadata.c
 * GEO_DEM.c
 * GEO_ephem_attit.c
 * GEO_global_arrays.h
 * GEO_locate_one_scan.c
 * GEO_get_geoid.c
 * GEO_interp_mirr_enc.c
 * GEO_write_input_metadata.c

## 6.2 - 2018-04-09

 - isolated MODIS extract metadata code
 - removed extra libs from modis CMakeLists.txt
 - Collection "6.1" of MODIS L1 processing code: l1agen_modis v6.0.6, geogen_modis v6.1.0, l1bgen_modisa v6.2.1, l1bgen_modist v6.2.2.
 - Cleaned up our changes to MODIS L1 code
 - Updated MODIS L1 code to be compatible with opt (aka lib3) r14227, and simplified changes wrt SDST delivery so as to make future integration easier.
 - moving MODIS L1 files out of subdirectory
 - Fixed geogen_modis scripts and code to run with shared libraries & new dir structure.

### Source Changed

  * Changelog.md
  * CMakeLists.txt
  * from_SDST/HISTORY.txt
  * from_SDST/MOD03.mcf
  * from_SDST/MOD_PR03.mk
  * from_SDST/MOD_PR03.pcf
  * from_SDST/MOD_PR03.pff.doc
  * from_SDST/MOD_PR03_pr.txt
  * from_SDST/MYD03.mcf
  * from_SDST/MYD_PR03.pcf
  * from_SDST/PACKINGLIST.txt
  * from_SDST/README.txt
  * GEO_abs_limit_check.c
  * GEO_aggregate.c
  * GEO_basic.h
  * GEO_check_ea_headers.c
  * GEO_create_swath.c
  * GEO_del_limit_check.c
  * GEO_DEM.c
  * GEO_DEM.h
  * GEO_derived_products.c
  * GEO_earth.h
  * GEO_earth_location.c
  * GEO_ellip_position.c
  * GEO_ephem_attit.c
  * GEO_find_next_flag.c
  * GEO_geo.h
  * GEO_get_bounding_coords.c
  * GEO_get_ephatt_inputs.c
  * GEO_get_geoid.c
  * GEO_get_GRing_points.c
  * GEO_get_inst_mirr_normal.c
  * GEO_get_sample_time.c
  * GEO_get_T_inst2ecr.c
  * GEO_get_utcpole_metadata.c
  * GEO_get_version_metadata.c
  * GEO_get_view_vec.c
  * GEO_global_arrays.h
  * GEO_hires.c
  * GEO_initialize_product.c
  * GEO_input.h
  * GEO_inst.h
  * GEO_interp_ECR.c
  * GEO_interp_mirr_ang.c
  * GEO_interp_mirr_enc.c
  * GEO_landsea_mask.c
  * GEO_locate_one_granule.c
  * GEO_locate_one_scan.c
  * GEO_location_main.c
  * GEO_main_func.h
  * GEO_main.h
  * GEO_maneuver.c
  * GEO_mat_vec_mul3.c
  * GEO_output.h
  * GEO_parameters.h
  * GEO_poly_coef1.c
  * GEO_poly_fit.c
  * GEO_prepare_l1a_data.c
  * GEO_prepare_mirr_data.c
  * GEO_product.h
  * GEO_read_L1AECS_metadata.c
  * GEO_read_L1Apacket_data.c
  * GEO_read_L1Ascan_metadata.c
  * GEO_read_L1Aspecific_metadata.c
  * GEO_read_L1Atemp_data.c
  * GEO_read_param_file.c
  * GEO_solar_and_lunar_vectors.c
  * GEO_terrain_correct.c
  * GEO_update_L1A_metadata.c
  * GEO_util.h
  * GEO_validate_derived_products.c
  * GEO_validate_earth_location.c
  * GEO_validation.h
  * GEO_vec_length3.c
  * GEO_vec_mul3.c
  * GEO_vec_prod3.c
  * GEO_vec_unit3.c
  * GEO_write_ECS_metadata.c
  * GEO_write_geospecific_metadata.c
  * GEO_write_granule_metadata.c
  * GEO_write_input_metadata.c
  * GEO_write_one_scan.c
  * GEO_write_parameters.c
  * GEO_write_scan_data.c
  * GEO_write_scan_metadata.c
  * HISTORY.txt
  * imsl_d_lin_sol_gen.c
  * imsl_d_spline_interp.c
  * imsl_d_spline_value.c
  * imsl_error.c
  * imsl_wrap.h
  * L1a_data.h
  * MOD03.mcf
  * MODIS_35251.t
  * MOD_PR03.mk
  * MOD_PR03.pcf
  * MOD_PR03.pff.doc
  * MOD_PR03_pr.txt
  * MYD03.mcf
  * MYD_PR03.pcf
  * PACKINGLIST.txt
  * PGS_MODIS_35251.h
  * pseudo_imsl.h
  * README.txt
  * version.h

## <VERSION> - 2015-03-12
### Added
  * Changelog.md
