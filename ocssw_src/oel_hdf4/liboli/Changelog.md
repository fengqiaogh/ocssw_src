
# liboli Changelog

## <VERSION> - 2015-03-12

 - get everything ready for shared libs
 - modified calls to liboli to eliminate dependence on l2gen
 - pass geometries around as double rather than scaled short int
 - adding L5/7 and geometry per band support
 - added ability to handle the new OLI extended meta data angles file.
 - had to increase the size of the variable that reads the file OLI and file name
  
### Source Changed

  * alberfor.c
  * alberinv.c
  * alconfor.c
  * alconinv.c
  * azimfor.c
  * aziminv.c
  * br_gctp.c
  * Changelog.md
  * CMakeLists.txt
  * cproj.c
  * cproj.h
  * emeta_api.c
  * emeta_api.h
  * emeta.c
  * emeta_exploit.c
  * emeta_exploit.h
  * emeta_geometry.c
  * emeta_geometry.h
  * emeta.h
  * eqconfor.c
  * eqconinv.c
  * equifor.c
  * equiinv.c
  * for_init.c
  * gctp.c
  * gctp_create_transformation.c
  * gctp_dms2degrees.c
  * gctp.h
  * gctp_print_message.c
  * gctp_report.c
  * gctp_transform.c
  * gctp_utility.c
  * geographic.c
  * get_l57tm_nom_angles.c
  * get_l8_angles.c
  * get_oli_nom_angles.c
  * gnomfor.c
  * gnominv.c
  * goodfor.c
  * goodinv.c
  * gvnspfor.c
  * gvnspinv.c
  * hamfor.c
  * haminv.c
  * ias_const.h
  * ias_geo_convert_deg2dms.c
  * ias_geo_convert_dms2deg.c
  * ias_geo_convert_geod2cart.c
  * ias_geo_find_deg.c
  * ias_geo_find_min.c
  * ias_geo_find_sec.c
  * ias_geo.h
  * ias_geo_projection_transformation.c
  * ias_logging.c
  * ias_logging.h
  * ias_math_compute_3dvec_dot.c
  * ias_math_compute_unit_vector.c
  * ias_math_compute_vector_length.c
  * ias_math.h
  * ias_math_interpolate_lagrange_3dvec.c
  * ias_odl_free_tree.c
  * ias_odl_get_field.c
  * ias_odl.h
  * ias_odl_read_tree.c
  * ias_structures.h
  * ias_types.h
  * imolwfor.c
  * imolwinv.c
  * inv_init.c
  * isinfor.c
  * isin.h
  * isininv.c
  * lablib3.c
  * lablib3.h
  * lamazfor.c
  * lamazinv.c
  * lambert_conformal_conic.c
  * local.h
  * merfor.c
  * merinv.c
  * millfor.c
  * millinv.c
  * molwfor.c
  * molwinv.c
  * mtl_geometry.c
  * mtl_geometry.h
  * mtl_grid.c
  * mtl_grid.h
  * obleqfor.c
  * obleqinv.c
  * oblique_mercator.c
  * orthfor.c
  * orthinv.c
  * paksz.c
  * polar_stereographic.c
  * polyconic.c
  * proj.h
  * robfor.c
  * robinv.c
  * sinfor.c
  * sininv.c
  * smeta_api.c
  * smeta_api.h
  * smeta_exploit.c
  * smeta_exploit.h
  * smeta_geometry.c
  * smeta_geometry.h
  * smeta.h
  * som.c
  * sphdz.c
  * spload.c
  * state_plane.c
  * state_plane_table.h
  * sterfor.c
  * sterinv.c
  * tm.c
  * toolbox.h
  * untfz.c
  * vandgfor.c
  * vandginv.c
  * wivfor.c
  * wivinv.c
  * wviifor.c
  * wviiinv.c

## <VERSION> - 2015-03-12
### Added
  * Changelog.md
