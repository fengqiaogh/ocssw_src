module ch_xfr
  !  this will contain all data needed for chimaera use for a line
  !  (and any surrounding lines) from l2gen
  !
  !  Also, any new data to transfer out
  !
  !  W. Robinson, SAIC, Nov 2018
  !
  ! l1b data
  real, allocatable, dimension( :,:,: ) :: c2_refl, c2_bnd_unc
  real, allocatable, dimension( : ) :: c2_spec_unc, c2_unc_sf
  real, allocatable, dimension( :,: ) :: c2_lat, c2_lon, c2_senz
  real, allocatable, dimension( :,: ) :: c2_sena, c2_solz, c2_sola, c2_relaz
  ! anc data
  real, allocatable, dimension( :,:,: ) :: c2_prof_mixr, c2_prof_t
  real, allocatable, dimension( :,:,: ) :: c2_prof_p, c2_prof_hgt
  real, allocatable, dimension( :,: ) :: c2_tsfc, c2_psfc, c2_tot_o3 
  real, allocatable, dimension( :,: ) :: c2_wind
  real, allocatable, dimension( :,: ) :: c2_ice_frac, c2_snow_frac
  real, allocatable, dimension( :,: ) :: c2_alt_o3, c2_alt_icec
  ! anc sfc albedo
  real, allocatable, dimension( :,:,: ) :: c2_sfc_albedo
  ! cloud top hgt, p, t transfer arrays
  real, allocatable, dimension( :,: ) :: c2_cth, c2_ctp, c2_ctt
  !
  integer*1, allocatable, dimension( :,: ) :: c2_sfc_lvl, c2_trop_lvl
  integer, allocatable, dimension( :,: ) :: c2_alt_snowice
  ! cloud mask data
  integer*1, allocatable, dimension( :,: ) :: c2_cld_det, c2_conf_cld 
  integer*1, allocatable, dimension( :,: ) :: c2_clr_66, c2_clr_95, c2_clr_99
  integer*1, allocatable, dimension( :,: ) :: c2_sno_sfc, c2_wtr_sfc
  integer*1, allocatable, dimension( :,: ) :: c2_coast_sfc, c2_desert_sfc
  integer*1, allocatable, dimension( :,: ) :: c2_lnd_sfc, c2_night, c2_glint
  integer*1, allocatable, dimension( :,: ) :: c2_ocean_no_snow, c2_ocean_snow 
  integer*1, allocatable, dimension( :,: ) :: c2_lnd_no_snow, c2_lnd_snow
  integer*1, allocatable, dimension( :,: ) :: c2_tst_h_138, c2_tst_vis_refl
  integer*1, allocatable, dimension( :,: ) :: c2_tst_vis_ratio, c2_vis_cld_250
  integer*1, allocatable, dimension( :,: ) :: c2_appl_hcld_138
  integer*1, allocatable, dimension( :,: ) :: c2_appl_vis_refl
  integer*1, allocatable, dimension( :,: ) :: c2_appl_vis_nir_ratio
  !
  !  band presence array
  integer*1, allocatable, dimension( : ) :: c2_cmp_there
  !
  ! dimensions, location in whole chunck
  integer :: c2_npix, c2_nlin, c2_nbnd, c2_nbnd_albedo, c2_nlvl_model, &
    c2_scan, c2_st_samp, c2_g_year, c2_g_day
  !
  !  c2_sensor_id, the identification of which sensor is being used
  integer :: c2_sensor_id
  integer :: OCIS_ID = 31  ! code for OCIS from sensorDefs.h
  integer :: OCI_ID = 30  ! code for OCIS from sensorDefs.h
  !
  !  the cm_from_l2 is a switch to say read from the l2 inputs instead of 
  !  the work file
  integer :: cm_from_l2
  character(len=500) :: c2_cloud_hgt_file
  !
  !
  !  And now for data to transfer out to l2gen -- diagnostic for the point refl
  !  relative to the table
  real, allocatable, dimension( :,:,: ) :: prd_out_refl_loc_2100
  real, allocatable, dimension( :,:,: ) :: prd_out_refl_loc_1600
  real, allocatable, dimension( :,:,: ) :: prd_out_refl_loc_1621
  real, allocatable, dimension( :,:,: ) :: prd_out_refl_loc_2200
  integer :: tbl_lo_abs=1, tbl_hi_abs=2, tbl_lo_abs_ice=3, tbl_hi_abs_ice=4, &
    refl_abs=5, tbl_lo_nabs=6, tbl_hi_nabs=7, tbl_lo_nabs_ice=8, &
    tbl_hi_nabs_ice=9, refl_nabs=10

  ! Also the save arrays, so the 3-line results are saved
  real, allocatable, dimension( :,:,: ) :: prd_out_refl_loc_2100_sav
  real, allocatable, dimension( :,:,: ) :: prd_out_refl_loc_1600_sav
  real, allocatable, dimension( :,:,: ) :: prd_out_refl_loc_2200_sav
  real, allocatable, dimension( :,:,: ) :: prd_out_refl_loc_1621_sav
end module ch_xfr
