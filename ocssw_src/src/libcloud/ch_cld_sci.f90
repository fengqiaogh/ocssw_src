subroutine ch_cld_sci( tdat, nv, ubdat, nubyte, i32dat, ni32, &
  sensor_id, cloud_hgt_file )
!
!  ch_cld_sci  is the entry to start the CHIMAERA processing for the 3-
!    line radiance (and other) input data for 1 line of data
!  
!  tdat, nv  Float array of all real data and the number of values
!  ubdat, nubyte  Unsigned byte values for processing and the number of values
!  i32dat, ni32  integer values for processing and the number of values
!  sensor_id  integer sensor identification number from l2gen
!  cloud_hgt_file  string name of the cloud height file
!
! (see below for the unpacking of these values into all inputs)
!
!  W. Robinson, SAIC, Dec 2018
!
use ch_xfr
implicit none

integer :: npix, nlin, nbnd, nbnd_albedo, nlvl_model, scan, st_samp, &
   g_year, g_day
common /dim_ctl/ npix, nlin, nbnd, nbnd_albedo, nlvl_model, scan, st_samp, &
   g_year, g_day
integer :: nv, i,j,k, na, off, nubyte, ni32, sensor_id, nbnd_sfc_albedo
real, dimension(nv) :: tdat
integer, dimension(3) :: dim_3
integer*1, dimension(nubyte) :: ubdat
integer, dimension(ni32) :: i32dat
character(len=500) :: cloud_hgt_file
real, parameter :: bad_float = -32767.

c2_sensor_id = sensor_id  !  a global sensor name
if( ( c2_sensor_id == OCI_ID ) .OR. ( c2_sensor_id == OCIS_ID ) ) then
  nbnd_sfc_albedo = 6
else
  nbnd_sfc_albedo = 5
endif

!  set up the l2gen values that will go into Chimaera
if( .not. allocated( c2_refl ) )  then
  allocate(c2_refl(npix, nlin, nbnd))
  allocate(c2_bnd_unc(npix, nbnd_albedo, nlin))
  allocate(c2_spec_unc(nbnd_albedo))
  allocate(c2_unc_sf(nbnd_albedo))
  allocate(c2_lat(npix, nlin))
  allocate(c2_lon(npix, nlin))
  allocate(c2_senz(npix, nlin))
  allocate(c2_sena(npix, nlin))
  allocate(c2_sola(npix, nlin))
  allocate(c2_solz(npix, nlin))
  allocate(c2_relaz(npix, nlin))
  allocate(c2_prof_mixr(npix,nlin,nlvl_model))
  allocate(c2_prof_t(npix,nlin,nlvl_model))
  allocate(c2_prof_p(npix,nlin,nlvl_model))
  allocate(c2_prof_hgt(npix,nlin,nlvl_model))
  allocate(c2_tsfc(npix,nlin))
  allocate(c2_psfc(npix,nlin))
  allocate(c2_wind(npix,nlin))
  allocate(c2_tot_o3(npix,nlin))
  allocate(c2_ice_frac(npix,nlin))
  allocate(c2_snow_frac(npix,nlin))
  allocate(c2_alt_o3(npix,nlin))
  allocate(c2_alt_icec(npix,nlin))
  allocate(c2_sfc_albedo(npix,nlin,nbnd_sfc_albedo))
  allocate(c2_sfc_lvl(npix,nlin))
  allocate(c2_trop_lvl(npix,nlin))
  allocate(c2_cld_det(npix,nlin))
  allocate(c2_conf_cld(npix,nlin))
  allocate(c2_clr_66(npix,nlin))
  allocate(c2_clr_95(npix,nlin))
  allocate(c2_clr_99(npix,nlin))
  allocate(c2_sno_sfc(npix,nlin))
  allocate(c2_wtr_sfc(npix,nlin))
  allocate(c2_coast_sfc(npix,nlin))
  allocate(c2_desert_sfc(npix,nlin))
  allocate(c2_lnd_sfc(npix,nlin))
  allocate(c2_night(npix,nlin))
  allocate(c2_glint(npix,nlin))
  allocate(c2_ocean_no_snow(npix,nlin))
  allocate(c2_ocean_snow(npix,nlin))
  allocate(c2_lnd_no_snow(npix,nlin))
  allocate(c2_lnd_snow(npix,nlin))
  allocate(c2_tst_h_138(npix,nlin))
  allocate(c2_tst_vis_refl(npix,nlin))
  allocate(c2_tst_vis_ratio(npix,nlin))
  allocate(c2_vis_cld_250(npix,nlin))
  allocate(c2_appl_hcld_138(npix,nlin))
  allocate(c2_appl_vis_refl(npix,nlin))
  allocate(c2_appl_vis_nir_ratio(npix,nlin))
  allocate(c2_alt_snowice(npix,nlin))
  allocate(c2_cmp_there(nbnd))
  ! for out point diagnostics
  allocate(prd_out_refl_loc_2100(npix,nlin,10))
  allocate(prd_out_refl_loc_1600(npix,nlin,10))
  allocate(prd_out_refl_loc_2200(npix,nlin,10))
  allocate(prd_out_refl_loc_1621(npix,nlin,10))
  !
  allocate(c2_cth(npix,nlin))
  allocate(c2_ctp(npix,nlin))
  allocate(c2_ctt(npix,nlin))
  endif
c2_npix = npix
c2_nbnd = nbnd
c2_nlin = nlin
c2_nbnd_albedo = nbnd_albedo
c2_nlvl_model = nlvl_model
c2_scan = scan
c2_st_samp = st_samp
c2_g_year = g_year
c2_g_day = g_day
!
!  Start with the L1b information, rads, geoloc, view angles
na = npix * nbnd * nlin
dim_3(1) = npix
dim_3(2) = nbnd
dim_3(3) = nlin
c2_refl = reshape( tdat(1:na), (/ npix, nbnd, nlin /) )
!
off = na
na = npix * nbnd_albedo * nlin
c2_bnd_unc = reshape( tdat(1+off:na+off), (/ npix, nbnd_albedo, nlin /) )
!
off = off + na
na = nbnd_albedo
c2_spec_unc = tdat(1+off:na+off)
!
off = off + na
c2_unc_sf = tdat(1+off:na+off)
!
off = off + na
na = npix * nlin
c2_lat = reshape( tdat(1+off:na+off), (/ npix, nlin /) )
!
off = off + na
c2_lon = reshape( tdat(1+off:na+off), (/ npix, nlin /) )
!
off = off + na
c2_senz = reshape( tdat(1+off:na+off), (/ npix, nlin /) )
where( c2_senz >= 90. ) c2_senz = bad_float
! 
off = off + na 
c2_sena = reshape( tdat(1+off:na+off), (/ npix, nlin /) )
!
off = off + na
c2_sola = reshape( tdat(1+off:na+off), (/ npix, nlin /) )
!
off = off + na
c2_solz = reshape( tdat(1+off:na+off), (/ npix, nlin /) )
where( c2_solz >= 90. ) c2_solz = bad_float
!
off = off + na
c2_relaz = reshape( tdat(1+off:na+off), (/ npix, nlin /) )
!
!  The met information, profiles, 2D fields
off = off + na
na = npix * nlin * nlvl_model
c2_prof_mixr = reshape( tdat(1+off:na+off), (/npix, nlin,nlvl_model/) )
!
off = off + na
c2_prof_t = reshape( tdat(1+off:na+off), (/npix, nlin,nlvl_model/) )
!
off = off + na
c2_prof_p = reshape( tdat(1+off:na+off), (/npix, nlin,nlvl_model/) )
!
off = off + na
c2_prof_hgt = reshape( tdat(1+off:na+off), (/npix, nlin,nlvl_model/) )
!
off = off + na
na = npix * nlin
c2_tsfc = reshape( tdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_psfc = reshape( tdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_wind = reshape( tdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_tot_o3 = reshape( tdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_ice_frac = reshape( tdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_snow_frac = reshape( tdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_alt_o3 = reshape( tdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_alt_icec = reshape( tdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_sfc_albedo = reshape( tdat(1+off: na*nbnd_sfc_albedo + off), &
  (/npix, nlin, nbnd_sfc_albedo/))
!
!  the new cloud top height, press, temp
off = off + na * nbnd_sfc_albedo
c2_cth = reshape( tdat(1+off:na+off), (/npix, nlin/) )
off = off + na
c2_ctp = reshape( tdat(1+off:na+off), (/npix, nlin/) )
off = off + na
c2_ctt = reshape( tdat(1+off:na+off), (/npix, nlin/) )
!
!  byte met data
na = npix * nlin
off = 0
c2_sfc_lvl = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_trop_lvl = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
!  the rest are the cloud mask byte data
off = off + na
c2_cld_det = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_conf_cld = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_clr_66 = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_clr_95 = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_clr_99 = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_sno_sfc = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_wtr_sfc = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_coast_sfc = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_desert_sfc = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_lnd_sfc = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_night = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_glint = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_ocean_no_snow = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_ocean_snow = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_lnd_no_snow = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_lnd_snow = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_tst_h_138 = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_tst_vis_refl = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_tst_vis_ratio = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_vis_cld_250 = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_appl_hcld_138 = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_appl_vis_refl = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_appl_vis_nir_ratio = reshape( ubdat(1+off:na+off), (/npix, nlin/) )
!
off = off + na
c2_cmp_there = ubdat( 1+off:nbnd+off )
!
!
!  do the integer data
na = npix * nlin
c2_alt_snowice = reshape( i32dat(1:na), (/npix, nlin/) )
!
!  set a global height file name
c2_cloud_hgt_file = cloud_hgt_file
!
! Here's what needs to be done:
! - Make local values of the dim sizes from the common values
! - Extract real arrays to correctly allocated science-familliar arrays
! - call modified scienceinterface with the right stuff
!i = 1
!k = 2
!  init the output diagnostics
! WDR this done in general_science_module already prd_out_refl_loc_2100 = -999.
! WDR rewind prd_out_refl_loc_1600 = -999.
!prd_out_refl_loc_1600 = -999.
!
!  we zero out the L1B data and met less the center line for test purposes
!  and to see what needs surrounding lines
! WDR out now test done call zero_surround( )
!
!  OK, segway into the chimaera code
call Driver_MOD_PR06OD( )
 return 
end subroutine ch_cld_sci

!
!  so zero_surround will just zero the values in surrounding L1b lines
!  (for now) mainly to test that no area data is needed in doing the 
!  cloud work except for the cloud stuff
!
subroutine zero_surround( )
use ch_xfr
integer :: clin

clin = c2_nlin / 2 + 1
!  reflectance
c2_refl(:,:,1:clin - 1) = 0.
c2_refl(:,:,clin + 1:c2_nlin ) = 0.
!  uncertainty
c2_bnd_unc(:,:,1:clin - 1) = 0. 
c2_bnd_unc(:,:,1:clin + 1:c2_nlin ) = 0.
!  lat  ** leave out - quite possibly, used for the anc (which is still
!          processed in the CHIMAERA code)
!c2_lat(:,1:clin - 1) = 0. 
!c2_lat(:,clin + 1:c2_nlin ) = 0.
!  lon
!c2_lon(:,1:clin - 1) = 0.
!c2_lon(:,clin + 1:c2_nlin ) = 0.
!  senz
!c2_senz(:,1:clin - 1) = 0.
!c2_senz(:,clin + 1:c2_nlin ) = 0.
!  sena
!c2_sena(:,1:clin - 1) = 0.
!c2_sena(:,clin + 1:c2_nlin ) = 0.
!  solz
!c2_solz(:,1:clin - 1) = 0.
!c2_solz(:,clin + 1:c2_nlin ) = 0.
!  sola
!c2_sola(:,1:clin - 1) = 0.
!c2_sola(:,clin + 1:c2_nlin ) = 0.
!  relaz
!c2_relaz(:,1:clin - 1) = 0.
!c2_relaz(:,clin + 1:c2_nlin ) = 0.
end subroutine zero_surround

subroutine get_cmp_prod( prod_num, prod_array, n_prd )
!
!  get_cmp_prod  will grab the requested product from the CHIMAERA final 
!    data arrays
!
!  prod_num  integer  I  The l2gen-defined product ID number
!  prod_array  real   O  array to carry out the line of data
!  n_prd      integer I  # of sub-products in this product, ie the product 
!                        has c2_npix pixels and n_prd sub-products
!
!  W. Robinson, SAIC, 8 Feb 2019
use ch_xfr, only: c2_nlin, c2_npix, prd_out_refl_loc_2100, &
    prd_out_refl_loc_1600,prd_out_refl_loc_1621, &
    prd_out_refl_loc_2200, &
    tbl_lo_abs, tbl_hi_abs, tbl_lo_abs_ice, &
    tbl_hi_abs_ice, refl_abs, tbl_lo_nabs, tbl_hi_nabs, &
    tbl_lo_nabs_ice, tbl_hi_nabs_ice, refl_nabs
use core_arrays, only: effective_radius_21_final, cloud_top_pressure, &
  liquid_water_path_16, liquid_water_path_1621, liquid_water_path, &
  optical_thickness_1621_final, effective_radius_1621_final, &
  optical_thickness_16_final, optical_thickness_final, &
  effective_radius_16_final, cloud_top_temperature, atm_corr_refl, &
  surface_albedo, failure_metric, failure_metric_16, failure_metric_1621, &
  cloud_top_height, &
  effective_radius_22_final, optical_thickness_22_final, &
  liquid_water_path_22, failure_metric_22
  integer :: prod_num, n_prd, iprd
  !real, dimension(c2_npix, n_prd ) :: prod_array
  real, dimension(n_prd,c2_npix) :: prod_array
!  the same catalog #s defined in l2prod.h
integer, parameter :: CAT_Cld_p = 469, CAT_Cld_t = 470, CAT_Cld_h = 493
integer, parameter :: CAT_CER_2100 = 445
integer, parameter :: CAT_CER_1600 = 446, CAT_COT_2100 = 447
integer, parameter :: CAT_COT_1600 = 448, CAT_CER_1621 = 449
integer, parameter :: CAT_COT_1621 = 450, CAT_CWP_2100 = 451
integer, parameter :: CAT_CWP_1621 = 452, CAT_CWP_1600 = 453
!  these are used in the byte call get_cmp_byte()
!integer, parameter :: CAT_Cld_Sfc_Type = 454, CAT_Cld_Phase_2100 = 455
!integer, parameter :: CAT_Cld_Non_Abs_Band = 456, CAT_Cld_Phase_1600 = 457
!integer, parameter :: CAT_Cld_Phase_1621 = 458
integer, parameter :: CAT_Cld_Top_Refl_650 = 459
integer, parameter :: CAT_Cld_Top_Refl_860 = 460, CAT_Cld_Top_Refl_1200 = 461
integer, parameter :: CAT_Cld_Top_Refl_1600 = 462, CAT_Cld_Top_Refl_2100 = 463
integer, parameter :: CAT_Surface_Albedo_650 = 464, CAT_Surface_Albedo_860 = 465
integer, parameter :: CAT_Surface_Albedo_1200 = 466
integer, parameter :: CAT_Surface_Albedo_1600 = 467
integer, parameter :: CAT_Surface_Albedo_2100 = 468

integer, parameter :: CAT_COT_fail_2100 = 471
integer, parameter :: CAT_COT_fail_1600 = 472
integer, parameter :: CAT_COT_fail_1621 = 473
integer, parameter :: CAT_CER_fail_2100 = 474
integer, parameter :: CAT_CER_fail_1600 = 475
integer, parameter :: CAT_CER_fail_1621 = 476
integer, parameter :: CAT_CMP_fail_pct_2100 = 477
integer, parameter :: CAT_CMP_fail_pct_1600 = 478
integer, parameter :: CAT_CMP_fail_pct_1621 = 479
integer, parameter :: CAT_refl_loc_1600 = 480
integer, parameter :: CAT_refl_loc_2100 = 481
integer, parameter :: CAT_refl_loc_1621 = 482
! recent additions for the 2.2 um band
integer, parameter :: CAT_CER_2200 = 483
integer, parameter :: CAT_COT_2200 = 484
integer, parameter :: CAT_CWP_2200 = 485
! byte integer, parameter :: CAT_Cld_Phase_2200 = 486
integer, parameter :: CAT_Cld_Top_Refl_2200 = 487
integer, parameter :: CAT_Surface_Albedo_2200 = 488
integer, parameter :: CAT_COT_fail_2200 = 489
integer, parameter :: CAT_CER_fail_2200 = 490
integer, parameter :: CAT_CMP_fail_pct_2200 = 491
integer, parameter :: CAT_refl_loc_2200 = 492

real, parameter :: bad_float = -32767.
!
!  for the product type, transfer the center line
! 
select case( prod_num )
  case( CAT_Cld_p )
    prod_array( 1, : ) = cloud_top_pressure( :, c2_nlin / 2 + 1 )
  case( CAT_Cld_t )
    prod_array( 1, : ) = cloud_top_temperature( :, c2_nlin / 2 + 1 )
  case( CAT_Cld_h )
    prod_array( 1, : ) = cloud_top_height( :, c2_nlin / 2 + 1 )
  case( CAT_CER_2100 )
    prod_array( 1, : ) = effective_radius_21_final( :, c2_nlin / 2 + 1 )
  case( CAT_CER_1600 )
    prod_array( 1, : ) = effective_radius_16_final( :, c2_nlin / 2 + 1 )
  case( CAT_CER_2200 )
    prod_array( 1, : ) = effective_radius_22_final( :, c2_nlin / 2 + 1 )
  case( CAT_COT_2100 )
    prod_array( 1, : ) = optical_thickness_final( :, c2_nlin / 2 + 1 )
  case( CAT_COT_1600 )
    prod_array( 1, : ) = optical_thickness_16_final( :, c2_nlin / 2 + 1 )
  case( CAT_COT_2200 )
    prod_array( 1, : ) = optical_thickness_22_final( :, c2_nlin / 2 + 1 )
  case( CAT_CER_1621 )
    prod_array( 1, : ) = effective_radius_1621_final( :, c2_nlin / 2 + 1 )
  case( CAT_COT_1621 )
    prod_array( 1, : ) = optical_thickness_1621_final( :, c2_nlin / 2 + 1 )
  case( CAT_CWP_2100 )
    prod_array( 1, : ) = liquid_water_path( :, c2_nlin / 2 + 1 )
  case( CAT_CWP_1621 )
    prod_array( 1, : ) = liquid_water_path_1621( :, c2_nlin / 2 + 1 )
  case( CAT_CWP_1600 )
    prod_array( 1, : ) = liquid_water_path_16( :, c2_nlin / 2 + 1 )
  case( CAT_CWP_2200 )
    prod_array( 1, : ) = liquid_water_path_22( :, c2_nlin / 2 + 1 )
  case( CAT_Cld_Top_Refl_650 )
    prod_array( 1, : ) = atm_corr_refl( 1, :, c2_nlin / 2 + 1 )
  case( CAT_Cld_Top_Refl_860 )
    prod_array( 1, : ) = atm_corr_refl( 2, :, c2_nlin / 2 + 1 )
  case( CAT_Cld_Top_Refl_1200 )
    prod_array( 1, : ) = atm_corr_refl( 3, :, c2_nlin / 2 + 1 )
  case( CAT_Cld_Top_Refl_1600 )
    prod_array( 1, : ) = atm_corr_refl( 4, :, c2_nlin / 2 + 1 )
  case( CAT_Cld_Top_Refl_2100 )
    prod_array( 1, : ) = atm_corr_refl( 5, :, c2_nlin / 2 + 1 )
  case( CAT_Cld_Top_Refl_2200 )
    prod_array( 1, : ) = atm_corr_refl( 6, :, c2_nlin / 2 + 1 )
  case( CAT_Surface_Albedo_650 )
    prod_array( 1, : ) = 0.001 * surface_albedo( :, c2_nlin / 2 + 1, 1 )
  case( CAT_Surface_Albedo_860 )
    prod_array( 1, : ) = 0.001 * surface_albedo( :, c2_nlin / 2 + 1, 2 )
  case( CAT_Surface_Albedo_1200 )
    prod_array( 1, : ) = 0.001 * surface_albedo( :, c2_nlin / 2 + 1, 3 )
  case( CAT_Surface_Albedo_1600 )
    prod_array( 1, : ) = 0.001 * surface_albedo( :, c2_nlin / 2 + 1, 4 )
  case( CAT_Surface_Albedo_2100 )
    prod_array( 1, : ) = 0.001 * surface_albedo( :, c2_nlin / 2 + 1, 5 )
  case( CAT_Surface_Albedo_2200 )
    prod_array( 1, : ) = 0.001 * surface_albedo( :, c2_nlin / 2 + 1, 6 )

  case( CAT_COT_fail_2100 )
    prod_array( 1, : ) = failure_metric( :, c2_nlin / 2 + 1 )%tau / 100.
    where ( prod_array < -90. ) prod_array = bad_float
  case( CAT_COT_fail_1600 )
    prod_array( 1, : ) = failure_metric_16( :, c2_nlin / 2 + 1 )%tau / 100.
    where ( prod_array < -90. ) prod_array = bad_float
  case( CAT_COT_fail_1621 )
    prod_array( 1, : ) = failure_metric_1621( :, c2_nlin / 2 + 1 )%tau / 100.
    where ( prod_array < -90. ) prod_array = bad_float
  case( CAT_COT_fail_2200 )
    prod_array( 1, : ) = failure_metric_22( :, c2_nlin / 2 + 1 )%tau / 100.
    where ( prod_array < -90. ) prod_array = bad_float
  case( CAT_CER_fail_2100 )
    prod_array( 1, : ) = failure_metric( :, c2_nlin / 2 + 1 )%re / 100.
    where ( prod_array < -90. ) prod_array = bad_float
  case( CAT_CER_fail_1600 )
    prod_array( 1, : ) = failure_metric_16( :, c2_nlin / 2 + 1 )%re / 100.
    where ( prod_array < -90. ) prod_array = bad_float
  case( CAT_CER_fail_1621 )
    prod_array( 1, : ) = failure_metric_1621( :, c2_nlin / 2 + 1 )%re / 100.
    where ( prod_array < -90. ) prod_array = bad_float
  case( CAT_CER_fail_2200 )
    prod_array( 1, : ) = failure_metric_22( :, c2_nlin / 2 + 1 )%re / 100.
    where ( prod_array < -90. ) prod_array = bad_float
  case( CAT_CMP_fail_pct_2100 )
    prod_array( 1, : ) = failure_metric( :, c2_nlin / 2 + 1 )%RSS / 100.
    where ( prod_array < -90. ) prod_array = bad_float
  case( CAT_CMP_fail_pct_1600 )
    prod_array( 1, : ) = failure_metric_16( :, c2_nlin / 2 + 1 )%RSS / 100.
    where ( prod_array < -90. ) prod_array = bad_float
  case( CAT_CMP_fail_pct_1621 )
    prod_array( 1, : ) = failure_metric_1621( :, c2_nlin / 2 + 1 )%RSS / 100.
    where ( prod_array < -90. ) prod_array = bad_float
  case( CAT_CMP_fail_pct_2200 )
    prod_array( 1, : ) = failure_metric_22( :, c2_nlin / 2 + 1 )%RSS / 100.
    where ( prod_array < -90. ) prod_array = bad_float
  case( CAT_refl_loc_1600 )
    do iprd=1,n_prd
      prod_array( iprd, :) = prd_out_refl_loc_1600( :, c2_nlin / 2 + 1, iprd )
    end do
  !
  case( CAT_refl_loc_2100 )
    do iprd=1,n_prd
      prod_array( iprd, :) = prd_out_refl_loc_2100( :, c2_nlin / 2 + 1, iprd )
    end do
  case( CAT_refl_loc_1621 )
    do iprd=1,n_prd
      prod_array( iprd, :) = prd_out_refl_loc_1621( :, c2_nlin / 2 + 1, iprd )
    end do
  case( CAT_refl_loc_2200 )
    do iprd=1,n_prd
      prod_array( iprd, :) = prd_out_refl_loc_2200( :, c2_nlin / 2 + 1, iprd )
    end do
  case default
    print*, "Improper product ID, case encountered, ID:", prod_num
    prod_array = bad_float
  end select
!
!  set bad value to obpg standard
where(prod_array < -900. )
  prod_array = bad_float
end where
return
end subroutine get_cmp_prod

subroutine get_cmp_byt( prod_num, bprod )
!
!  get_cmp_byt is for all the flags / indexes of byte size
!
!  prod_num  integer    I  The l2gen-defined product ID number
!  bprod     integer*1  O  array to carry out the line of data
!
!  W. Robinson, SAIC, 5 Aug 2019
use ch_xfr, only: c2_nlin, c2_npix
use core_arrays, only: cloudmask, processing_information
integer :: prod_num, lin_mid
integer*1, dimension(c2_npix) :: bprod
integer, parameter :: CAT_Cld_Sfc_Type = 454, CAT_Cld_Phase_2100 = 455
integer, parameter :: CAT_Cld_Non_Abs_Band = 456, CAT_Cld_Phase_1600 = 457
integer, parameter :: CAT_Cld_Phase_1621 = 458
integer, parameter :: CAT_Cld_Phase_2200 = 486
integer, parameter :: CAT_Cld_water_cloud = 440, CAT_Cld_ice_cloud = 441
integer*1 :: iand_comp = 7

!
!  for the product type, transfer the center line

lin_mid = c2_nlin / 2 + 1
prod_sel: select case( prod_num )
  case( CAT_Cld_Sfc_Type )
    bprod = 0
    where(cloudmask(:,lin_mid)%ocean_snow == 1 ) &
      bprod = 1
    where(cloudmask(:,lin_mid)%land_no_snow == 1 ) &
      bprod = 2
    where(cloudmask(:,lin_mid)%land_snow == 1 ) &
      bprod = 3
  case( CAT_Cld_Phase_2100 )
    bprod = processing_information(:,lin_mid)%path_and_outcome
    bprod = ibclr( bprod, 3 )
  case( CAT_Cld_Non_Abs_Band )
    bprod = processing_information(:,lin_mid)%band_used_for_optical_thickness
  case( CAT_Cld_Phase_1600 )
    bprod = processing_information(:,lin_mid)%path_and_outcome_16
    bprod = ibclr( bprod, 3 )
  case( CAT_Cld_Phase_1621 )
    bprod = processing_information(:,lin_mid)%path_and_outcome_1621
    bprod = ibclr( bprod, 3 )
  case( CAT_Cld_Phase_2200 )
    bprod = processing_information(:,lin_mid)%path_and_outcome_22
    bprod = ibclr( bprod, 3 )
  case( CAT_Cld_water_cloud )
    bprod = 0
    where(iand(processing_information(:,lin_mid)%path_and_outcome,iand_comp) == 2 )
      bprod = 1
    elsewhere(iand(processing_information(:,lin_mid)%path_and_outcome,iand_comp) == 4 )
      bprod = 1
    endwhere
  case( CAT_Cld_ice_cloud )
    bprod = 0
    where(iand(processing_information(:,lin_mid)%path_and_outcome,iand_comp) == 3 ) &
      bprod = 1
  case default
    print*, "Improper product ID, case encountered, ID:", prod_num
    return

  end select prod_sel
  return
end subroutine get_cmp_byt
