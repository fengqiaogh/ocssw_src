module corescience_module

  implicit none
  private

  public :: corescience

  contains


subroutine corescience(xpoint,ypoint,process,measurements, Tc_liquid, &
  Tc_ice, debug, na_band_used, nearest_liq, nearest_ice, RSS_liq, RSS_ice, &
  alt_ray_liq, alt_ray_ice, do_retrievals_liq, do_retrievals_ice, status)

!  Arguments
!  xpoint   IN  pixel location in radiance chunk array
!  ypoint  IN  line location in radiance chunk array
!  process    IN  structure of pixel's cloud status, type processflag
!  measurements  IN  vector of top of cloud reflectances
!  Tc_liquid  IN/OUT  cloud temperature for liquid cloud gotten from 
!      cloud_top_temperature (IN) and may be re-calculated (OUT) - only 
!      used for 3.7 um band
!  Tc_ice  IN/OUT  cloud temperature for ice cloud - same in/out as above
!      - only for 3.7 um band
!  debug  IN  debugging switch USED?
!  na_band_used  OUT non-absorbing band used, 650 (land), 860 (water), or
!                    1240 (snow)
!  nearest_liq   OUT  size 4 array of 1 - nearest used, 0 - not used
!        for use of absorbing band: 1 - 1.6 um, 2 - 2.1 um, 3 - 3.7 | 2.2 um,
!        4 - 1.6 / 2.1 um bands
!  nearest_ice   OUT  size 4 array of 1 - nearest used, 0 - not used
!  RSS_liq    IN?OUT  size 4 array of some root sum square uncertainty for
!         each of the abs bands - liquid --OR-- a distance away from 
!         the closest solution in the alternate solution subplot/logic
!  RSS_ice  IN?OUT  size 4 array of some root sum square uncertainty for
!         each of the abs bands - ice
!   --- WDR lets try understanding std algorithm before alternate one (ASL)
!  alt_ray_liq
!  alt_ray_ice
!  do_retrievals_liq
!  do_retrievals_ice
!  status
   
  use GeneralAuxType
  use core_arrays
  use libraryarrays
  use libraryinterpolates
  use modis_cloudstructure, only: processflag
  use modis_numerical_module, only: linearinterpolation
  use get_retrieval_uncertainty
  use science_parameters
  use nonscience_parameters
  use mod06_run_settings
  use retrieval_solution_logic
  use retrieval_prep_logic
   
  use specific_other
  !  WDR for output of the some diagnostics
  use ch_xfr, only: prd_out_refl_loc_2100,  &
    prd_out_refl_loc_2200, &
    prd_out_refl_loc_1600, prd_out_refl_loc_1621, &
    tbl_lo_abs, tbl_hi_abs, tbl_lo_abs_ice, &
    tbl_hi_abs_ice, refl_abs, tbl_lo_nabs, tbl_hi_nabs, tbl_lo_nabs_ice, &
    tbl_hi_nabs_ice, refl_nabs, c2_sensor_id, OCI_ID, OCIS_ID

  implicit none
 
  logical,           intent(in)       :: debug
  type(processflag), intent(in)       :: process
  real, dimension(:), intent(in) :: measurements
  logical, dimension(:), intent(in) :: do_retrievals_liq, do_retrievals_ice
  real, intent(inout) :: Tc_liquid, Tc_ice, alt_ray_liq, alt_ray_ice
  integer,           intent(in)       :: xpoint,ypoint
  integer,           intent(out)      :: status
  integer, intent(out) :: na_band_used

  integer, dimension(:), intent(inout) :: nearest_liq, nearest_ice
  real, dimension(:), intent(inout) :: RSS_liq, RSS_ice

  type(cloudphase) :: local_process
  integer   :: nonabsorbingband_index, absorbingband_index, maxradii, &
    lib_vnir_index
  real                 :: thermal_trans_1way, thermal_trans_2way
  real, allocatable    :: optical_thickness_vector(:), residual(:)
  real :: retrievalopticalthickness, retrievalopticalthickness16,  &
          retrievalopticalthickness37,  retrievalopticalthickness1621,&
          retrievalradius21,  retrievalradius16,  retrievalradius37,  &
          retrievalradius1621, retrievalradius22, retrievalopticalthickness22

  integer :: idx16, idx21, idx37, channel_37, idx11, i, idx_alb37, &
    idx_alb16, quality_in, idx22
  real :: curTc, newTc, ray_temp
  logical :: PRN, use_nearest
  integer :: iii, idx1621(2)

  status = 0

  thermal_trans_1way = thermal_correction_oneway(1)   ! only the 3.7um is used
  thermal_trans_2way = thermal_correction_twoway(1)   ! only the 3.7um is used
 
  !  local band naming - var names are mostly obvious, more || less
  call get_band_idx(idx16, idx21, idx37, channel_37, idx11, idx_alb37, &
    idx_alb16)
  !  for OCI, get the 2.2
  if( ( c2_sensor_id == OCI_ID ) .OR. ( c2_sensor_id == OCIS_ID ) ) &
    idx22 = band_0226 - 1
 
  nonabsorbingband_index = band_0065
  processing_information(xpoint, ypoint)%band_used_for_optical_thickness =  0
  !
  !  set non-absorbing band based on the surface type
  !
  if (process%ocean_surface .or. process%coastal_surface) then!{
    nonabsorbingband_index = band_0086
    endif!}
 
  if (process%land_surface .or. process%desert_surface) then!{      
    nonabsorbingband_index = band_0065
    endif!}
  if (process%snowice_surface) then!{
    ! 1.2 um
    nonabsorbingband_index = band_0124
    endif!}
 
  if (platform_name == "AVIRIS") then 
    if (nonabsorbingband_index == 4) &
      processing_information(xpoint, ypoint)%band_used_for_optical_thickness = 3
    if (nonabsorbingband_index == 3) &
      processing_information(xpoint, ypoint)%band_used_for_optical_thickness = 2
    if (nonabsorbingband_index == 1) &
      processing_information(xpoint, ypoint)%band_used_for_optical_thickness = 1  
  else
    if (nonabsorbingband_index > 3) then 
      processing_information(xpoint, ypoint)%band_used_for_optical_thickness &
        = nonabsorbingband_index-1   
    else 
      processing_information(xpoint, ypoint)%band_used_for_optical_thickness &
        = nonabsorbingband_index   
    endif
  endif
  
  na_band_used = nonabsorbingband_index
 
  local_process%watercloud = 0
  local_process%icecloud = 0

  ! initialize all retrievals to fillvalue   
  retrievalopticalthickness = fillvalue_real
  retrievalopticalthickness1621 = fillvalue_real
  retrievalopticalthickness16 = fillvalue_real
  retrievalopticalthickness37  = fillvalue_real
  retrievalradius16 = fillvalue_real
  retrievalradius21 = fillvalue_real
  retrievalradius1621 = fillvalue_real
  retrievalradius37 = fillvalue_real
  if( ( c2_sensor_id == OCI_ID ) .OR. ( c2_sensor_id == OCIS_ID ) )then
    retrievalradius22 = fillvalue_real
    retrievalopticalthickness22 = fillvalue_real
    endif

  nearest_liq = 0
  nearest_ice = 0
  
  local_process%watercloud = 1

  !  set surface QA

  !  if the measured surface reflectance in the non absorbing bands is 
  !  greater than the surface albedo , and we have a cloudy pixel then
  !  we can try a retrieval

  if (process%cloudobserved ) then!{

    !  check for .86um saturation
    !  if the .86um saturates (wisconsin reader sets band_meas to neg. 
    !  if saturation) or there can be bad noisy neg. reflectance on low end,
    !  try the other .65um band (and change the QA flag)
    if (nonabsorbingband_index == band_0086 .and. &
      measurements(nonabsorbingband_index) < 0.) then!{
      if (measurements(band_0065) > 0.) nonabsorbingband_index = band_0065
      na_band_used = nonabsorbingband_index
      processing_information(xpoint, ypoint)%band_used_for_optical_thickness &
        = nonabsorbingband_index
    endif!}

    nonabsorbingband_index = na_band_used

    !  WDR save the used reflectance - we still need the table ranges found
    prd_out_refl_loc_2100(xpoint,ypoint,refl_nabs) = &
      measurements(nonabsorbingband_index)
    prd_out_refl_loc_1600(xpoint,ypoint,refl_nabs) = &
      measurements(nonabsorbingband_index)
    prd_out_refl_loc_1621(xpoint,ypoint,refl_nabs) = &
      measurements(band_0163)
    prd_out_refl_loc_2200(xpoint,ypoint,refl_nabs) = &
      measurements(nonabsorbingband_index)
!
    prd_out_refl_loc_2100(xpoint,ypoint,refl_abs) = &
      measurements(band_0213)
    prd_out_refl_loc_1600(xpoint,ypoint,refl_abs) = &
      measurements(band_0163)
    prd_out_refl_loc_1621(xpoint,ypoint,refl_abs) = &
      measurements(band_0213)
    prd_out_refl_loc_2200(xpoint,ypoint,refl_abs) = &
      measurements(band_0226)

    ! this is for rotten snow albedoes
    if (albedo_real4( nonabsorbingband_index) < 0. .or. &
      albedo_real4( band_0163) < 0. .or. &
      albedo_real4( band_0213) < 0.) return

    ! if the shortwave saturates, try the other shortwave band 
    allocate(optical_thickness_vector(number_waterradii))
    allocate(residual(number_waterradii))
    residual = 0  ! WDR UIV
    
    allocate(reflibA(number_taus+1, number_waterradii), &
      reflibB(number_taus+1, number_waterradii))
    reflibA = 0 ! WDR-UIV
    reflibB = 0 ! WDR-UIV
    allocate(rad37lib(number_taus+1, number_waterradii))
    
    rad37lib    = -999.0

    lib_vnir_index = nonabsorbingband_index
    !
    !  For water clouds, for the non-absorbing reflectance, find the 
    !  tau (optical thickness) values in the table for every eff. radius
    !
    call vis_nonabsorbing_science(measurements(nonabsorbingband_index), &
      lib_vnir_index, albedo_real4(nonabsorbingband_index), &
      library_taus, water_radii, spherical_albedo_water, &
      int_fluxdownwater_sensor, int_fluxdownwater_solar, & 
      int_reflectance_water,  sensor_zenith_angle(xpoint,ypoint), &
      solar_zenith_angle(xpoint,ypoint), &
      relative_azimuth_angle(xpoint,ypoint),  &
      cloud_top_pressure(xpoint,ypoint), local_process, &
      optical_thickness_vector)

    ! an option to knock out the SWIR channels completely so that we only 
    !  do optical thickness retrieval assuming an effective radius
    !  This appears to handle the other COT band values too
#if NOSWIR
    ! Retrieval using assumed effective radius
    retrievalradius16 = 14.0
    call solve_retrieval_noSWIR(optical_thickness_vector, water_radii, &
      retrievalradius16, retrievalopticalthickness16)

    retrievalopticalthickness = retrievalopticalthickness16
    retrievalradius21 = retrievalradius16

    ! Retrieval using assumed effective radius, minus 1-sigma CER
    retrievalradius16 = 6.0
    call solve_retrieval_noSWIR(optical_thickness_vector, water_radii, &
      retrievalradius16, retrievalopticalthickness16)

    optical_thickness_37_final(xpoint,ypoint) = &
      abs(retrievalopticalthickness-retrievalopticalthickness16)

    ! Retrieval using assumed effective radius, plus 1-sigma CER
    retrievalradius16 = 22.0
    call solve_retrieval_noSWIR(optical_thickness_vector, water_radii, &
      retrievalradius16, retrievalopticalthickness16)

    optical_thickness_37_final(xpoint,ypoint) = &
      optical_thickness_37_final(xpoint,ypoint) + &
      abs(retrievalopticalthickness-retrievalopticalthickness16)

    ! Calculate mean COT error due to unknown CER
    optical_thickness_37_final(xpoint,ypoint) = &
      optical_thickness_37_final(xpoint,ypoint) / 2.0

    retrievalopticalthickness1621 = fillvalue_real
        retrievalradius1621 = fillvalue_real

#else
    !  below is with SWIR
    !
    !  For 1.6 um, continue and find the re, tau using the absorbing 
    !  band reflectance  for water clouds
    !
    if (do_retrievals_liq(1)) then 
      absorbingband_index = band_0163 ! 1.6um band 
      if(measurements(absorbingband_index) > 0. ) then!{ 
           
        call vis_absorbing_science(optical_thickness_vector, &
          measurements(absorbingband_index), idx16, &
          albedo_real4(idx_alb16), library_taus, water_radii, &
          spherical_albedo_water, int_fluxdownwater_sensor,  &
          int_fluxdownwater_solar, int_reflectance_water, &
          residual,maxradii, debug)
                      
        if (maxradii > 2) then!{
          call solveretrieval(residual(1:maxradii), &
            optical_thickness_vector(1:maxradii),   &
            water_radii(1:maxradii), retrievalradius16, &
            retrievalopticalthickness16, debug, use_nearest, quality_in)
          if (use_nearest) then 
            nearest_liq(1) = 1
            call solveretrieval_nearest(xpoint,ypoint, &
              measurements(nonabsorbingband_index), &
              measurements(absorbingband_index), &
              (/nonabsorbingband_index, absorbingband_index/), &
              water_radii, retrievalopticalthickness16, &
              retrievalradius16, RSS_liq(re16), &
              .true., ray_temp, quality_in )    
          endif
    
        else!}{
          retrievalradius16 = fillvalue_real
          retrievalopticalthickness16 = fillvalue_real
        endif!}

      else!}{
        retrievalradius16 = fillvalue_real
        retrievalopticalthickness16 = fillvalue_real
      endif!}
    else
      retrievalradius16 = fillvalue_real
      retrievalopticalthickness16 = fillvalue_real
    endif
    !
    !  WDR get the table range
    prd_out_refl_loc_1600(xpoint,ypoint,tbl_lo_nabs) = MINVAL(reflibA)
    prd_out_refl_loc_1600(xpoint,ypoint,tbl_hi_nabs) = MAXVAL(reflibA)
    prd_out_refl_loc_1600(xpoint,ypoint,tbl_lo_abs) = MINVAL(reflibB)
    prd_out_refl_loc_1600(xpoint,ypoint,tbl_hi_abs) = MAXVAL(reflibB)
    !
    !  For 2.1 um, continue and find the re, tau using the absorbing
    !  band reflectance for water clouds
    !

    if (do_retrievals_liq(2) ) then 
      absorbingband_index = band_0213 ! 2.1um

      call vis_absorbing_science(optical_thickness_vector, &
        measurements(absorbingband_index), idx21, &
        albedo_real4(absorbingband_index), library_taus, water_radii, &
        spherical_albedo_water, int_fluxdownwater_sensor, &
        int_fluxdownwater_solar, int_reflectance_water, &
        residual, maxradii, debug)

      if (maxradii > 2) then!{  
        call solveretrieval(residual(1:maxradii), &
          optical_thickness_vector(1:maxradii), water_radii(1:maxradii), &
          retrievalradius21, retrievalopticalthickness, &
          debug, use_nearest, quality_in)
                     
        if (use_nearest) then 
          nearest_liq(2) = 1
          call solveretrieval_nearest(xpoint,ypoint, &
            measurements(nonabsorbingband_index), &
            measurements(absorbingband_index), &
            (/nonabsorbingband_index, absorbingband_index/), &
            water_radii, retrievalopticalthickness, retrievalradius21, &
            RSS_liq(re21), .true., alt_ray_liq, quality_in)
        endif
    
      else!}{
        retrievalradius21 = fillvalue_real
        retrievalopticalthickness = fillvalue_real
      endif!}
    else
      retrievalradius21 = fillvalue_real
      retrievalopticalthickness = fillvalue_real
    endif
    !
    !  WDR get the table range
    prd_out_refl_loc_2100(xpoint,ypoint,tbl_lo_nabs) = MINVAL(reflibA)
    prd_out_refl_loc_2100(xpoint,ypoint,tbl_hi_nabs) = MAXVAL(reflibA)
    prd_out_refl_loc_2100(xpoint,ypoint,tbl_lo_abs) = MINVAL(reflibB)
    prd_out_refl_loc_2100(xpoint,ypoint,tbl_hi_abs) = MAXVAL(reflibB)
    !

#endif
    !  WDR best do 2.2 here - same as 1.6 and 2.1
    if( ( c2_sensor_id == OCI_ID ) .OR. ( c2_sensor_id == OCIS_ID ) )then
      if (do_retrievals_liq(3) ) then
        absorbingband_index = band_0226 ! 2.2um
  
        call vis_absorbing_science(optical_thickness_vector, &
          measurements(absorbingband_index), idx22, &
          albedo_real4(absorbingband_index), library_taus, water_radii, &
          spherical_albedo_water, int_fluxdownwater_sensor, &
          int_fluxdownwater_solar, int_reflectance_water, &
          residual, maxradii, debug)
  
        if (maxradii > 2) then!{
          call solveretrieval(residual(1:maxradii), &
            optical_thickness_vector(1:maxradii), water_radii(1:maxradii), &
            retrievalradius22, retrievalopticalthickness22, &
            debug, use_nearest, quality_in)
  
          if (use_nearest) then
            nearest_liq(3) = 1
            call solveretrieval_nearest(xpoint,ypoint, &
              measurements(nonabsorbingband_index), &
              measurements(absorbingband_index), &
              (/nonabsorbingband_index, absorbingband_index/), &
              water_radii, retrievalopticalthickness22, retrievalradius22, &
              RSS_liq(re22), .true., alt_ray_liq, quality_in)
          endif
  
        else!}{
          retrievalradius22 = fillvalue_real
          retrievalopticalthickness22 = fillvalue_real
        endif!}
      else
        retrievalradius22 = fillvalue_real
        retrievalopticalthickness22 = fillvalue_real
      endif
      !
      !  WDR get the table range
      prd_out_refl_loc_2200(xpoint,ypoint,tbl_lo_nabs) = MINVAL(reflibA)
      prd_out_refl_loc_2200(xpoint,ypoint,tbl_hi_nabs) = MAXVAL(reflibA)
      prd_out_refl_loc_2200(xpoint,ypoint,tbl_lo_abs) = MINVAL(reflibB)
      prd_out_refl_loc_2200(xpoint,ypoint,tbl_hi_abs) = MAXVAL(reflibB)
    else   ! END 2.2 um 
      !
      !
      !  For 3.7 um, continue and find the re, tau using the absorbing
      !  band reflectance for water clouds
      !
      if (do_retrievals_liq(3)) then 
        absorbingband_index = band_0370 ! 3.7um
  
        ! sep, 4 May: absorbingband_index - 1 being passed in for all libarary 
        !  indices because 
        !  band_0370 = 7 and library index corresponding to this band 
        !  is 6. In the libraries, band_0935 is index 7.
        ! wind, 7 Dec: absorbingband_index - 1 has to be passed 
        !in for albedo as well. 
  
        call nir_absorbing_science(platform_name, optical_thickness_vector, &
          measurements(absorbingband_index), idx37, albedo_real4( idx_alb37), &
          xpoint, ypoint, cloud_top_temperature(xpoint, ypoint), &
          thermal_trans_1way, thermal_trans_2way, library_taus, &
          water_radii, spherical_albedo_water, int_fluxdownwater_sensor, &
          int_fluxdownwater_solar, int_fluxupwater_sensor, &
          int_reflectance_water, int_cloud_emissivity_water, &
          int_surface_emissivity_water, residual, maxradii, channel_37, &
          emission_uncertainty_pw_liq, emission_uncertainty_Tc_liq, &
          sigma_R37_PW_liq,debug)
  
        if ( maxradii > 2 ) then!{
          call solveretrieval(residual(1:maxradii), &
            optical_thickness_vector(1:maxradii), &
            water_radii(1:maxradii), retrievalradius37, &
            retrievalopticalthickness37, debug, use_nearest, quality_in)
                       
          if (use_nearest) then 
            nearest_liq(3) = 1
            call solveretrieval_nearest(xpoint,ypoint, &
              measurements(nonabsorbingband_index), &
              measurements(absorbingband_index), &
              (/nonabsorbingband_index, absorbingband_index/), &
              water_radii, retrievalopticalthickness37, retrievalradius37, &
              RSS_liq(re37), .true., ray_temp, quality_in,     &
              CH37_IDX = idx37, CTopT = cloud_top_temperature(xpoint, ypoint), &
              CH37_NUM =channel_37 , platFormName= platform_name)
          endif
  
        else!}{
          retrievalradius37 = fillvalue_real
          retrievalopticalthickness37 = fillvalue_real
        endif!}
        !
        ! now we iterate the retrieval.  I believe the emission in 3.7 is
        ! enough to require this iteration
        !
#if !AMS_INST
  
        ! if the cloud too thick, don't bother, the iteration will not 
        !  accomplish absolutely anything
        if ( retrievalopticalthickness37 > 0. &
          .and. retrievalradius37 > 0. .and. &
          irw_temperature(xpoint, ypoint) > 0. .and. nearest_liq(3) == 0 ) then 
  
          ! if the cloud is too thick, then don't bother doing anything, but 
          !  still set the temperature
          if (retrievalopticalthickness37 >= 10.) then 
            Tc_liquid = cloud_top_temperature(xpoint, ypoint)
          else
            curTc = irw_temperature(xpoint, ypoint)
  
            do i=1, 5
              if (retrievalopticalthickness37 < 0. .or. &
                  retrievalradius37 < 0. ) then 
                retrievalradius37 = fillvalue_real
                retrievalopticalthickness37 = fillvalue_real
                Tc_liquid = curTc
                exit
              endif
          
              call calculate_new_Tc(platform_name, &
                irw_temperature(xpoint, ypoint), &
                surface_temperature(xpoint, ypoint), &
                1.- surface_emissivity_land(xpoint, ypoint, 2), idx11, &
                retrievalopticalthickness37, retrievalradius37, library_taus, &
                water_radii, spherical_albedo_water, int_fluxdownwater_sensor, &
                int_fluxupwater_sensor, int_cloud_emissivity_water, &
                int_surface_emissivity_water, newTc, .false.)
                    
              if (newTc < 0.) then 
                Tc_liquid = curTc       
                exit
              endif
  
              call nir_absorbing_science(platform_name, &
                optical_thickness_vector, &
                measurements(absorbingband_index), idx37, &
                albedo_real4( idx_alb37), &
                xpoint, ypoint, newTc, thermal_trans_1way, thermal_trans_2way, &
                library_taus, water_radii, spherical_albedo_water, &
                int_fluxdownwater_sensor, int_fluxdownwater_solar, &
                int_fluxupwater_sensor, int_reflectance_water, &
                int_cloud_emissivity_water, int_surface_emissivity_water, &
                residual, maxradii, channel_37, emission_uncertainty_pw_liq, &
                emission_uncertainty_Tc_liq, sigma_R37_PW_liq,debug)
        
              if ( maxradii > 2 ) then!{
                call solveretrieval(residual(1:maxradii), &
                  optical_thickness_vector(1:maxradii), water_radii(1:maxradii), &
                  retrievalradius37, retrievalopticalthickness37, &
                  debug, use_nearest, quality_in)
                if (use_nearest) then 
                  nearest_liq(3) = 1      
                  call solveretrieval_nearest(xpoint,ypoint, &
                    measurements(nonabsorbingband_index), &
                    measurements(absorbingband_index), &
                    (/nonabsorbingband_index, absorbingband_index/), &
                    water_radii, retrievalopticalthickness37, retrievalradius37, &
                    RSS_liq(re37), .true., ray_temp, quality_in, &
                    CH37_IDX = idx37, &
                    CTopT = cloud_top_temperature(xpoint, ypoint), &
                    CH37_NUM =channel_37 , platFormName= platform_name)
  
                endif
              else!}{
                retrievalradius37 = fillvalue_real
                retrievalopticalthickness37 = fillvalue_real
                exit
              endif!}
  
              if ( abs(curTc - newTc) < 0.01 .or. retrievalradius37 < 0.) then 
                Tc_liquid = newTc
                exit
              endif
  
              curTc = newTc
            end do
          endif
  
        endif
#endif    
      else
        retrievalradius37 = fillvalue_real
        retrievalopticalthickness37 = fillvalue_real
      endif  ! END 3.7
    endif

#if NOSWIR
    retrievalradius1621 = fillvalue_real
    retrievalopticalthickness1621 = fillvalue_real
#else
    !
    !  For 1.6 and 2.1 um, continue and find the re, tau using the absorbing
    !  band reflectance for water clouds
    !
    if (do_retrievals_liq(4)) then 
      retrievalradius1621 = fillvalue_real
      retrievalopticalthickness1621 = fillvalue_real
    
      if( (process%ocean_surface .or. process%snowice_surface) .and. &
        measurements(band_0163) > 0. ) then!{
        nonabsorbingband_index = band_0163
        absorbingband_index = band_0213
 
        !  get the COT for all radii that gives the non-abs refl at this point
        call vis_nonabsorbing_science(measurements(nonabsorbingband_index), &
          idx16, albedo_real4(nonabsorbingband_index), library_taus, &
          water_radii, spherical_albedo_water, int_fluxdownwater_sensor, &
          int_fluxdownwater_solar, int_reflectance_water, &
          sensor_zenith_angle(xpoint,ypoint), &
          solar_zenith_angle(xpoint,ypoint), &
          relative_azimuth_angle(xpoint,ypoint), &
          cloud_top_pressure(xpoint,ypoint), local_process, &
          optical_thickness_vector)
                
        call vis_absorbing_science(optical_thickness_vector, &
          measurements(absorbingband_index), idx21, &
          albedo_real4(absorbingband_index), library_taus, water_radii, &
          spherical_albedo_water, int_fluxdownwater_sensor, &
          int_fluxdownwater_solar, int_reflectance_water, &
          residual,maxradii, debug)
        if (maxradii > 2) then!{
          call solveretrieval(residual(1:maxradii), &
            optical_thickness_vector(1:maxradii), water_radii(1:maxradii), &
            retrievalradius1621, retrievalopticalthickness1621, debug, &
            use_nearest, quality_in)
          if (use_nearest) then 
            nearest_liq(4) = 1

            call solveretrieval_nearest(xpoint,ypoint, &
              measurements(nonabsorbingband_index), &
              measurements(absorbingband_index), &
              (/nonabsorbingband_index, absorbingband_index/), &
              water_radii, retrievalopticalthickness1621, &
              retrievalradius1621, RSS_liq(re1621), &
              .true., ray_temp, quality_in)
          endif
        endif!}

      endif!}
    else
      retrievalradius1621 = fillvalue_real
      retrievalopticalthickness1621 = fillvalue_real
    endif
    !
    !  WDR get the table range
    prd_out_refl_loc_1621(xpoint,ypoint,tbl_lo_nabs) = MINVAL(reflibA)
    prd_out_refl_loc_1621(xpoint,ypoint,tbl_hi_nabs) = MAXVAL(reflibA)
    prd_out_refl_loc_1621(xpoint,ypoint,tbl_lo_abs) = MINVAL(reflibB)
    prd_out_refl_loc_1621(xpoint,ypoint,tbl_hi_abs) = MAXVAL(reflibB)
    !
#endif
    !
    ! End setting of water path for water phase clouds
    !  Repeat entire process for ice clouds
    !
    deallocate(rad37lib)

    deallocate(optical_thickness_vector, residual, reflibA, reflibB)

    optical_thickness_liquid         = retrievalopticalthickness
    optical_thickness_1621_liquid    = retrievalopticalthickness1621
    optical_thickness_16_liquid    = retrievalopticalthickness16
    optical_thickness_37_liquid    = retrievalopticalthickness37
    effective_radius_16_liquid       = retrievalradius16
    effective_radius_21_liquid       = retrievalradius21
    effective_radius_1621_liquid     = retrievalradius1621
    effective_radius_37_liquid       = retrievalradius37
    optical_thickness_22_liquid    = retrievalopticalthickness22
    effective_radius_22_liquid       = retrievalradius22
  
    !  This was commented out in initiall chimaera code
    ! if (nearest_liq(1) == 1) effective_radius_16_liquid(xpoint, ypoint) = &
    !   fillvalue_real
    ! if (nearest_liq(2) == 1) then 
    !   effective_radius_21_liquid(xpoint, ypoint) = fillvalue_real
    !   optical_thickness_liquid(xpoint, ypoint) = fillvalue_real
    ! endif
    ! if (nearest_liq(3) == 1) effective_radius_37_liquid(xpoint, ypoint) = &
    !   fillvalue_real
    ! if (nearest_liq(4) == 1) then 
    !   effective_radius_1621_liquid(xpoint, ypoint) = fillvalue_real
    !   optical_thickness_1621_liquid(xpoint, ypoint) = fillvalue_real
    ! endif
    !
    !
    ! now reset everything so the ice cloud can be processed
    !
    local_process%watercloud = 0
    local_process%icecloud = 1

    nonabsorbingband_index = na_band_used
    absorbingband_index = band_0163
    
    allocate(optical_thickness_vector(number_iceradii))
    allocate(residual(number_iceradii))
    residual = 0 ! WDR UIV
    allocate(reflibA(number_taus+1, number_iceradii), &
      reflibB(number_taus+1, number_iceradii))
    allocate(rad37lib(number_taus+1, number_iceradii))
    ! WDR uiv sort of
    reflibA = 0
    reflibB = 0
    
    rad37lib    = -999.0

    ! initialize all retrievals to fillvalue   
    retrievalopticalthickness = fillvalue_real
    retrievalopticalthickness1621 = fillvalue_real
    retrievalopticalthickness16 = fillvalue_real
    retrievalopticalthickness37  = fillvalue_real
    retrievalradius16 = fillvalue_real
    retrievalradius21 = fillvalue_real
    retrievalradius1621 = fillvalue_real
    retrievalradius37 = fillvalue_real
    retrievalopticalthickness22 = fillvalue_real
    retrievalradius22 = fillvalue_real

    lib_vnir_index = nonabsorbingband_index
    !
    !  For ice clouds, find the tau (optical thickness) values in the
    !  table for the non-absorbing reflectance
    !
    call vis_nonabsorbing_science(measurements(nonabsorbingband_index), &
      lib_vnir_index, albedo_real4(nonabsorbingband_index), &
      library_taus, ice_radii, spherical_albedo_ice, &
      int_fluxdownice_sensor, int_fluxdownice_solar, int_reflectance_ice, &
      sensor_zenith_angle(xpoint,ypoint), &
      solar_zenith_angle(xpoint,ypoint),  &
      relative_azimuth_angle(xpoint,ypoint), &
      cloud_top_pressure(xpoint,ypoint),  &
      local_process, optical_thickness_vector)

    !  code for work without SWIR bands
#if NOSWIR 
    ! Retrieval using assumed effective radius
    retrievalradius16 = 30.0
    call solve_retrieval_noSWIR(optical_thickness_vector, ice_radii, &
      retrievalradius16, retrievalopticalthickness16)

    retrievalopticalthickness = retrievalopticalthickness16
    retrievalradius21 = retrievalradius16

    ! Retrieval using assumed effective radius, minus 1-sigma CER
    retrievalradius16 = 13.0
    call solve_retrieval_noSWIR(optical_thickness_vector, ice_radii, &
      retrievalradius16, retrievalopticalthickness16)

    optical_thickness_1621_final(xpoint,ypoint) = &
      abs(retrievalopticalthickness-retrievalopticalthickness16)

    ! Retrieval using assumed effective radius, plus 1-sigma CER
    retrievalradius16 = 37.0
    call solve_retrieval_noSWIR(optical_thickness_vector, ice_radii, &
      retrievalradius16, retrievalopticalthickness16)

    optical_thickness_1621_final(xpoint,ypoint) = &
      optical_thickness_1621_final(xpoint,ypoint) + &
      abs(retrievalopticalthickness-retrievalopticalthickness16)

    ! Calculate mean COT error due to unknown CER
    optical_thickness_1621_final(xpoint,ypoint) = &
      optical_thickness_1621_final(xpoint,ypoint) / 2.0

    retrievalopticalthickness1621 = fillvalue_real
        retrievalradius1621 = fillvalue_real
#else
    !
    !  For 1.6 um, continue and find the re, tau using the absorbing
    !  band reflectance for ice clouds
    !
    if (do_retrievals_ice(1)) then 
      if (measurements(absorbingband_index) > 0.) then!{ 
        call vis_absorbing_science(optical_thickness_vector,         &
          measurements(absorbingband_index), idx16, albedo_real4(idx_alb16), &
          library_taus, ice_radii, spherical_albedo_ice, &
          int_fluxdownice_sensor, int_fluxdownice_solar, &
          int_reflectance_ice, residual,maxradii, debug)
        if ( maxradii > 2) then!{
          call solveretrieval(residual(1:maxradii), &
            optical_thickness_vector(1:maxradii), &
            ice_radii(1:maxradii), retrievalradius16, &
            retrievalopticalthickness16, debug, use_nearest, quality_in)
          if (use_nearest) then 
            nearest_ice(1) = 1
      
            call solveretrieval_nearest(xpoint,ypoint, &
              measurements(nonabsorbingband_index), &
              measurements(absorbingband_index), &
              (/nonabsorbingband_index, absorbingband_index/), &
              ice_radii, retrievalopticalthickness16, retrievalradius16, &
              RSS_ice(re16), .false., ray_temp, quality_in)
          endif
        else!}{
          retrievalradius16 = fillvalue_real
          retrievalopticalthickness16 = fillvalue_real
        endif!}
           
      else!}{
        retrievalradius16 = fillvalue_real
        retrievalopticalthickness16 = fillvalue_real
      endif!}
    else
      retrievalradius16 = fillvalue_real
      retrievalopticalthickness16 = fillvalue_real
    endif
    !
    !  WDR get the table range
    prd_out_refl_loc_1600(xpoint,ypoint,tbl_lo_nabs_ice) = MINVAL(reflibA)
    prd_out_refl_loc_1600(xpoint,ypoint,tbl_hi_nabs_ice) = MAXVAL(reflibA)
    prd_out_refl_loc_1600(xpoint,ypoint,tbl_lo_abs_ice) = MINVAL(reflibB)
    prd_out_refl_loc_1600(xpoint,ypoint,tbl_hi_abs_ice) = MAXVAL(reflibB)
    !
    !
    !  For 2.1 um, continue and find the re, tau using the absorbing
    !  band reflectance
    !
    if (do_retrievals_ice(2)) then 
      absorbingband_index = band_0213
      call vis_absorbing_science(optical_thickness_vector, &
        measurements(absorbingband_index), idx21, &
        albedo_real4(absorbingband_index), library_taus, ice_radii, &
        spherical_albedo_ice, int_fluxdownice_sensor, int_fluxdownice_solar, &
        int_reflectance_ice, residual,maxradii, debug)

      if ( maxradii > 2) then!{
        call solveretrieval(residual(1:maxradii), &
          optical_thickness_vector(1:maxradii), ice_radii(1:maxradii), &
          retrievalradius21, retrievalopticalthickness, &
          debug, use_nearest, quality_in)

        if (use_nearest) then 
          nearest_ice(2) = 1
          call solveretrieval_nearest(xpoint,ypoint, &
            measurements(nonabsorbingband_index), &
            measurements(absorbingband_index), &
            (/nonabsorbingband_index, absorbingband_index/), ice_radii, &
            retrievalopticalthickness, retrievalradius21, RSS_ice(re21), &
            .false., alt_ray_ice, quality_in)
        endif
      else!}{
        retrievalradius21 = fillvalue_real
        retrievalopticalthickness = fillvalue_real
      endif!}
    else
      retrievalradius21 = fillvalue_real
      retrievalopticalthickness = fillvalue_real
    endif
    !
    !  WDR get the table range
    prd_out_refl_loc_2100(xpoint,ypoint,tbl_lo_nabs_ice) = MINVAL(reflibA)
    prd_out_refl_loc_2100(xpoint,ypoint,tbl_hi_nabs_ice) = MAXVAL(reflibA)
    prd_out_refl_loc_2100(xpoint,ypoint,tbl_lo_abs_ice) = MINVAL(reflibB)
    prd_out_refl_loc_2100(xpoint,ypoint,tbl_hi_abs_ice) = MAXVAL(reflibB)
    !
#endif        
    ! the 2.2 um
    if( ( c2_sensor_id == OCI_ID ) .OR. ( c2_sensor_id == OCIS_ID ) )then
      if (do_retrievals_ice(3)) then
        absorbingband_index = band_0226
        call vis_absorbing_science(optical_thickness_vector, &
          measurements(absorbingband_index), idx22, &
          albedo_real4(absorbingband_index), library_taus, ice_radii, &
          spherical_albedo_ice, int_fluxdownice_sensor, int_fluxdownice_solar, &
          int_reflectance_ice, residual,maxradii, debug)
  
        if ( maxradii > 2) then!{
          call solveretrieval(residual(1:maxradii), &
            optical_thickness_vector(1:maxradii), ice_radii(1:maxradii), &
            retrievalradius22, retrievalopticalthickness22, &
            debug, use_nearest, quality_in)
  
          if (use_nearest) then
            nearest_ice(3) = 1
            call solveretrieval_nearest(xpoint,ypoint, &
              measurements(nonabsorbingband_index), &
              measurements(absorbingband_index), &
              (/nonabsorbingband_index, absorbingband_index/), ice_radii, &
              retrievalopticalthickness22, retrievalradius22, RSS_ice(re22), &
              .false., alt_ray_ice, quality_in)
          endif
        else!}{
          retrievalradius22 = fillvalue_real
          retrievalopticalthickness = fillvalue_real
        endif!}
      else
        retrievalradius22 = fillvalue_real
        retrievalopticalthickness = fillvalue_real
      endif
      !
      !  WDR get the table range
      prd_out_refl_loc_2200(xpoint,ypoint,tbl_lo_nabs_ice) = MINVAL(reflibA)
      prd_out_refl_loc_2200(xpoint,ypoint,tbl_hi_nabs_ice) = MAXVAL(reflibA)
      prd_out_refl_loc_2200(xpoint,ypoint,tbl_lo_abs_ice) = MINVAL(reflibB)
      prd_out_refl_loc_2200(xpoint,ypoint,tbl_hi_abs_ice) = MAXVAL(reflibB)
    else
      !
      !  For 3.7 um, continue and find the re, tau using the absorbing
      !  band reflectance for ice clouds
      !
      if (do_retrievals_ice(3)) then 
        absorbingband_index = band_0370
  
        ! sep, 4 May: absorbingband_index - 1 being passed in for all 
        ! libarary indices because 
        !   band_0370 = 7 and library index corresponding to this band 
        !   is 6. In the 
        !   libraries, band_0935 is index 7.
        !
        ! wind, 7 Dec: absorbingband_index - 1 has to be passed in for 
        !  albedo as well. 
        !
        call nir_absorbing_science(platform_name, optical_thickness_vector, &
          measurements(absorbingband_index), idx37, albedo_real4( idx_alb37), &
          xpoint, ypoint, cloud_top_temperature(xpoint, ypoint), &
          thermal_trans_1way, thermal_trans_2way, library_taus, &
          ice_radii, spherical_albedo_ice, int_fluxdownice_sensor, &
          int_fluxdownice_solar, int_fluxupice_sensor, int_reflectance_ice, &
          int_cloud_emissivity_ice, int_surface_emissivity_ice, &
          residual, maxradii, channel_37, emission_uncertainty_pw_ice, &
          emission_uncertainty_Tc_ice, sigma_R37_PW_ice, debug)
      
        if ( maxradii > 2 ) then!{
          call solveretrieval(residual(1:maxradii), &
            optical_thickness_vector(1:maxradii), ice_radii(1:maxradii), &
            retrievalradius37, retrievalopticalthickness37, debug, use_nearest, &
            quality_in)
          if (use_nearest) then 
            nearest_ice(3) = 1
  
          call solveretrieval_nearest(xpoint,ypoint, &
            measurements(nonabsorbingband_index), &
            measurements(absorbingband_index), &
            (/nonabsorbingband_index, absorbingband_index/), &
            ice_radii, retrievalopticalthickness37, &
            retrievalradius37, RSS_ice(re37), .false., ray_temp, quality_in, &
            CH37_IDX = idx37, CTopT = cloud_top_temperature(xpoint, ypoint), &
            CH37_NUM =channel_37 , platFormName= platform_name)
          endif
  
        else!}{
          retrievalradius37 = fillvalue_real
          retrievalopticalthickness37 = fillvalue_real
        endif!}
                                       
        ! now we iterate the retrieval.  I believe the emission in 3.7 is
        ! enough to require this iteration
#if !AMS_INST  
        if ( retrievalopticalthickness37 > 0. &
          .and. retrievalradius37 > 0. .and. &
          irw_temperature(xpoint, ypoint) > 0. .and. nearest_ice(3) == 0 ) then 
  
          ! if the cloud is too thick, then don't bother doing anything, 
          ! 1 but still set the temperature
          if (retrievalopticalthickness37 > 10.) then 
            Tc_ice = cloud_top_temperature(xpoint, ypoint) 
          else
            curTc = irw_temperature(xpoint, ypoint)
      
            do i=1, 5
          
              if (retrievalopticalthickness37 < 0. .or. &
                retrievalradius37 < 0. ) then 
                retrievalradius37 = fillvalue_real
                retrievalopticalthickness37 = fillvalue_real
                Tc_ice = curTc
                exit
              endif
          
              call calculate_new_Tc(platform_name, &
                irw_temperature(xpoint, ypoint), &
                surface_temperature(xpoint, ypoint), &
                1.- surface_emissivity_land(xpoint, ypoint, 2), &
                idx11, retrievalopticalthickness37, retrievalradius37, &
                library_taus, ice_radii, spherical_albedo_ice,   &
                int_fluxdownice_sensor, int_fluxupice_sensor,   &
                int_cloud_emissivity_ice, int_surface_emissivity_ice, &
                newTc, .false.)
                    
              if (newTc < 0.) then 
                Tc_ice  = curTc       
                exit
              endif
                                                              
              call nir_absorbing_science(platform_name, &
                optical_thickness_vector, measurements(absorbingband_index), &
                idx37, albedo_real4( idx_alb37), xpoint, ypoint, &
                newTc, thermal_trans_1way, thermal_trans_2way, &
                library_taus, ice_radii, spherical_albedo_ice, &
                int_fluxdownice_sensor, int_fluxdownice_solar, &
                int_fluxupice_sensor, int_reflectance_ice, &
                int_cloud_emissivity_ice, int_surface_emissivity_ice, &
                residual, maxradii, channel_37, emission_uncertainty_pw_ice, &
                emission_uncertainty_Tc_ice, sigma_R37_PW_ice,debug)
      
              if ( maxradii > 2 ) then!{
                call solveretrieval(residual(1:maxradii), &
                  optical_thickness_vector(1:maxradii), &
                  ice_radii(1:maxradii), retrievalradius37, &
                  retrievalopticalthickness37, &
                  debug, use_nearest, quality_in)
                if (use_nearest) then 
                  nearest_ice(3) = 1
        
                  call solveretrieval_nearest(xpoint,ypoint, &
                    measurements(nonabsorbingband_index), &
                    measurements(absorbingband_index), &
                    (/nonabsorbingband_index, absorbingband_index/), &
                    ice_radii, retrievalopticalthickness37, retrievalradius37, &
                    RSS_ice(re37), .false., ray_temp, quality_in, &
                    CH37_IDX = idx37, &
                    CTopT = cloud_top_temperature(xpoint, ypoint), &
                    CH37_NUM =channel_37 , platFormName= platform_name)
                endif
  
              else!}{
                retrievalradius37 = fillvalue_real
                retrievalopticalthickness37 = fillvalue_real
              endif!}
  
              if ( abs(curTc - newTc) < 0.01 .or. retrievalradius37 < 0.) then 
                Tc_ice = newTc
                exit
              endif
              curTc = newTc
            end do
  
          endif
        endif
#endif    
  
      else 
        retrievalradius37 = fillvalue_real
        retrievalopticalthickness37 = fillvalue_real
      endif   ! end 3.7
    endif
    
#if NOSWIR 
    retrievalradius1621 = fillvalue_real
    retrievalopticalthickness1621 = fillvalue_real
#else
    !
    !  For 1.6, 2.1 um, continue and find the re, tau using the absorbing
    !  band reflectance for ice clouds
    !
    if (do_retrievals_ice(4)) then 
      retrievalradius1621 = fillvalue_real
      retrievalopticalthickness1621 = fillvalue_real
    
      if( (process%ocean_surface .or. process%snowice_surface) .and. &
        measurements(band_0163) > 0. ) then!{

        nonabsorbingband_index = band_0163
        absorbingband_index = band_0213
        idx1621(1) = idx16
        idx1621(2) = idx21
    
        ! the Nakajima-King space for VIIRS (and AHI) 1.6-2.2 um ice 
        !  retrieval appears to be reversed
        ! where the 2.2um channel contains the tau information and the 
        ! 1.6um channel contains the 
        ! re information. So we play along and reverse the channels. 
#ifdef VIIRS_INST
        if ( .not. MODIS_MODE ) then 
          nonabsorbingband_index = band_0213
          absorbingband_index = band_0163
          idx1621(1) = idx21
          idx1621(2) = idx16
        endif     
#endif
#if AHI_INST | AMS_INST
        nonabsorbingband_index = band_0213
        absorbingband_index = band_0163
        idx1621(1) = idx21
        idx1621(2) = idx16
#endif
#if AVIRIS_INST
        if (set_bands(absorbingband_index) == 14) then 
          nonabsorbingband_index = band_0213
          absorbingband_index = band_0163
          idx1621(1) = idx21
          idx1621(2) = idx16
        endif
#endif
        call vis_nonabsorbing_science(measurements(nonabsorbingband_index), &
          idx1621(1), albedo_real4(nonabsorbingband_index), library_taus, &
          ice_radii, spherical_albedo_ice, int_fluxdownice_sensor, &
          int_fluxdownice_solar, int_reflectance_ice, &
          sensor_zenith_angle(xpoint,ypoint), &
          solar_zenith_angle(xpoint,ypoint),  &
          relative_azimuth_angle(xpoint,ypoint), &
          cloud_top_pressure(xpoint,ypoint), &
          local_process, optical_thickness_vector)

        call vis_absorbing_science(optical_thickness_vector, &
          measurements(absorbingband_index), idx1621(2), &
          albedo_real4(absorbingband_index), library_taus, ice_radii, &
          spherical_albedo_ice, int_fluxdownice_sensor, &
          int_fluxdownice_solar, int_reflectance_ice, &
          residual,maxradii,debug)

        if (maxradii > 2 ) then!{
          call solveretrieval(residual(1:maxradii), &
            optical_thickness_vector(1:maxradii), ice_radii(1:maxradii), &
            retrievalradius1621, retrievalopticalthickness1621, &
            debug, use_nearest, quality_in)
          if (use_nearest) then 
            nearest_ice(4) = 1
            call solveretrieval_nearest(xpoint,ypoint, &
              measurements(nonabsorbingband_index), &
              measurements(absorbingband_index), &
              (/nonabsorbingband_index, absorbingband_index/), ice_radii, &
              retrievalopticalthickness1621, retrievalradius1621, &
              RSS_ice(re1621), .false., ray_temp, quality_in)
      
          endif
        endif!}
      endif!}
    else
      retrievalradius1621 = fillvalue_real
      retrievalopticalthickness1621 = fillvalue_real
    endif
    !
    !  WDR get the table range
    prd_out_refl_loc_1621(xpoint,ypoint,tbl_lo_nabs_ice) = MINVAL(reflibA)
    prd_out_refl_loc_1621(xpoint,ypoint,tbl_hi_nabs_ice) = MAXVAL(reflibA)
    prd_out_refl_loc_1621(xpoint,ypoint,tbl_lo_abs_ice) = MINVAL(reflibB)
    prd_out_refl_loc_1621(xpoint,ypoint,tbl_hi_abs_ice) = MAXVAL(reflibB)
    !
#endif
    !
    !  end of ice cloud processing
    !
    deallocate(rad37lib)    
    deallocate(optical_thickness_vector, residual, reflibA, reflibB)

  else!}{
    ! there is no cloud
    retrievalopticalthickness = fillvalue_real
    retrievalopticalthickness1621 = fillvalue_real
    retrievalopticalthickness16 = fillvalue_real
    retrievalopticalthickness37  = fillvalue_real
    retrievalradius16 = fillvalue_real
    retrievalradius1621 = fillvalue_real
    retrievalradius37 = fillvalue_real
    if( ( c2_sensor_id == OCI_ID ) .OR. ( c2_sensor_id == OCIS_ID ) ) then
      retrievalopticalthickness22 = fillvalue_real
      retrievalradius22 = fillvalue_real
    endif
    cloud_layer_flag(xpoint, ypoint) = 0
    status = 1
  endif!}

  ! so the statistics are computed properly. The retrieval is only 
  ! successful for statistical purposes if the main
  ! re retrieval is successful. 
  if (retrievalradius21 == fillvalue_real) then 
      status = 1
  endif

  !assign the retrievals to the relevant arrays
  optical_thickness_ice         = retrievalopticalthickness
  optical_thickness_1621_ice    = retrievalopticalthickness1621
  optical_thickness_16_ice    = retrievalopticalthickness16
  optical_thickness_37_ice    = retrievalopticalthickness37
  effective_radius_16_ice       = retrievalradius16
  effective_radius_21_ice       = retrievalradius21
  effective_radius_1621_ice     = retrievalradius1621
  effective_radius_37_ice       = retrievalradius37
  if( ( c2_sensor_id == OCI_ID ) .OR. ( c2_sensor_id == OCIS_ID ) ) then
    optical_thickness_22_ice    = retrievalopticalthickness22
    effective_radius_22_ice       = retrievalradius22
  endif

  ! reset the extra retrievals for now. 
  ! if (nearest_ice(1) == 1) effective_radius_16_ice(xpoint, ypoint) = &
  !   fillvalue_real
  ! if (nearest_ice(2) == 1) then 
  !   effective_radius_21_ice(xpoint, ypoint) = fillvalue_real
  !   optical_thickness_ice(xpoint, ypoint) = fillvalue_real
  ! endif
  ! if (nearest_ice(3) == 1) effective_radius_37_ice(xpoint, ypoint) = &
  !  fillvalue_real
  ! if (nearest_ice(4) == 1) then 
  !   effective_radius_1621_ice(xpoint, ypoint) = fillvalue_real
  !   optical_thickness_1621_ice(xpoint, ypoint) = fillvalue_real
  ! endif

end subroutine corescience

end module corescience_module
