module modis_science_module

implicit none

private

public:: scienceinterface



contains

subroutine scienceinterface(threshold_solar_zenith,      &
                            threshold_sensor_zenith,     &
                            threshold_relative_azimuth,  &
                            debug,                       &
                            status)

  use set_quality_data_module
   use modis_cloudstructure, only:cloudphase
   use core_arrays, only:cloudsummary, CSR_flag_array, cloudmask, &
      band_measurements, irw_temperature, model_info, &
      surface_emissivity_land, irw_temperature, cloud_top_pressure, &
      c2_model_info, surface_temperature, irw_temperature, &
      cloud_top_temperature, surface_albedo, tau_liquid, tau_ice, &
      re21_liquid, re21_ice, albedo_real4, processing_information, &
      atm_corr_refl, optical_thickness_final, optical_thickness_1621_final, &
      optical_thickness_16_final, optical_thickness_37_final, &
      effective_radius_16_final, effective_radius_21_final,  &
      effective_radius_37_final, effective_radius_1621_final, &
      cloud_layer_flag, ml_test_flag, optical_thickness_error, &
      effective_radius_21_error, liquid_water_path_error, &
      optical_thickness_16_error, effective_radius_16_error, &
      liquid_water_path_16_error, optical_thickness_16_error, &
      optical_thickness_37_error, effective_radius_37_error, &
      liquid_water_path_37_error, optical_thickness_1621_error, &
      effective_radius_1621_error, liquid_water_path_1621_error, &
      failure_metric, failure_metric_16, failure_metric_37, &
      failure_metric_1621, liquid_water_path, liquid_water_path_16, &
      liquid_water_path_37, liquid_water_path_1621, albedo_fac, &
      bprime_tc, bprime_ts, const_c, effective_radius_1621_ice, &
      effective_radius_1621_liquid, effective_radius_16_ice, &
      effective_radius_16_liquid, effective_radius_21_ice, &
      effective_radius_21_liquid, effective_radius_37_ice, &
      effective_radius_37_liquid, emission_uncertainty_pw_ice, &
      emission_uncertainty_pw_liq, emission_uncertainty_tc_ice, &
      emission_uncertainty_tc_liq, optical_thickness_1621_ice, &
      optical_thickness_1621_liquid, optical_thickness_16_ice, &
      optical_thickness_16_liquid, optical_thickness_37_ice, &
      optical_thickness_37_liquid, optical_thickness_ice, &
      optical_thickness_liquid, platform_name, precip_water_094, &
      sigma_r37_pw_ice, sigma_r37_pw_liq, tc_high_for_delta, &
      tc_low_for_delta, transprime_1way, transprime_2way, &
      solar_zenith_angle, latitude, cloud_height_method, &
      abovecloud_watervapor, cloud_phase_infrared, band_uncertainty, &
      transmittance_twoway, sensor_zenith_angle, longitude, &
      thermal_correction_oneway_low, thermal_correction_twoway_low, &
      cloud_mask_spi, meandelta_trans, relative_azimuth_angle, &
      thermal_correction_oneway_high, thermal_correction_twoway_high, &
      transmittance_stddev, scn_loop_st, scn_loop_en, &
      !  for the 2.2 um band
      optical_thickness_22_final,effective_radius_22_final, &
      optical_thickness_22_error, effective_radius_22_error, &
      liquid_water_path_22_error, optical_thickness_22_error, &
      failure_metric_22, liquid_water_path_22, &
      effective_radius_22_ice, effective_radius_22_liquid, &
      optical_thickness_22_ice, optical_thickness_22_liquid
      
   use libraryarrays, only:rayleigh_liq, rayleigh_ice, ice_radii, &
      library_taus, transmit_correction_table, water_radii
   use libraryinterpolates, only:int_surface_emissivity_water,  &
      int_reflectance_water, im_cloudy_count, im_ice_cloud_count, &
      im_successful_retrieval_count, im_undet_count, im_water_cloud_count, &
      iterationx, number_of_iterationsx, pi, pixx, pixy
   !use science_parameters
   use corescience_module, only:corescience
   use multi_layer_clouds, only:compute_multilayer_map
   use clear_sky_restoral, only:cloudiness_test
   use mod06_run_settings, only:set_albedo_bands, band_0065, band_0086, &
      band_0124, band_0163, band_0213, band_0370, band_1100, channel_37um, &
      set_number_of_bands, do_cox_munk, do_csr, force_ice, force_water, &
      re16, re1621, re21, re37, re22, set_bands, band_0935, band_0226
   use nonscience_parameters, only:fillvalue_int2, fillvalue_real, &
      max_tau_rtrieved
   use interpolate_libraries, only:scatangle, libraryinterpolate
   use global_model_grids, only:get_model_idx, get_model_idx_geos5
   use ancillary_module, only: get_above_cloud_properties, given_P_get_T
   use get_retrieval_uncertainty, only: init_half_radii, getuncertainties
   use retrieval_irw, only:retrieve_irw_temp
   use atmospheric_correction_module, only:atmospheric_correction
   use retrieval_prep_logic, only:color_first_time, cox_munk, d2r, delta_pc, &
      delta_ts, last_2way_angle, last_cox_munk, lastinterp_relative_azimuth, &
      lastinterp_scat_angle, lastinterp_scat_angle_ss, &
      lastinterp_sensor_zenith, lastinterp_solar_zenith, lastinterp_wind_speed,&
      retr_scale_factor, solar_constant_37, solar_zenith_threshold, swir_error,&
      unc_scale_factor, vnir_error, watervapor_error, uncertain_sf, &
      uncertain_sf, spec_uncertain, init_retrieval, cleanup_retrieval
   use rtm_support, only:rtm_cloud_prof, rtm_cloud_prof_high, &
      rtm_cloud_prof_low, rtm_rad_atm_clr, rtm_rad_atm_clr_high, &
      rtm_rad_atm_clr_low, rtm_trans_atm_clr, rtm_trans_atm_clr_high, &
      rtm_trans_atm_clr_low, init_rtm_vars, get_rtm_parameters, &
      get_rtm_parameters
   use cloud_phase, only:clouddecision
   use specific_other, only:set_cox_munk_albedo
   use modis_numerical_module, only: linearinterpolation, bisectionsearch
   use names, only: MY_UNIT_LUN
   use general_science_module, only: init_science_arrays, &
      set_drel, assign_retrieval_error, set_interp_controls, &
      set_water_path_answers, set_failure_answers, split_pcl, &
      capture_arrays
   use planck_functions, only: modis_planck
   use ch_xfr, only: c2_scan, c2_cmp_there, c2_sensor_id, OCI_ID, OCIS_ID
!WDR for the new routine
! out in this vers   use wdr_wr_ch_vars, only: wr_ch_vars
   
   implicit none

   logical, intent(in) :: debug
   real,    intent(in) :: threshold_solar_zenith,      &
                          threshold_sensor_zenith,     &
                          threshold_relative_azimuth

   integer, intent(out) :: status


   integer             :: xdimension, ydimension,i,j, k,jj,count_interpolations, &
                          retrievalcount
   integer             :: retrieval_failcount, library_failcount, cloud_failcount, atmoscorr_failcount
   integer             :: cloudstatus, retrievalstatus, librarystatus, atmoscorrstatus,cloudiness_degree_250m
   real                :: diff_solar_zenith,           &
                          diff_sensor_zenith,          &
                          diff_relative_azimuth
  
   real :: scattering_angle, dscat, diff_scat_angle, diff_scat_angle_ss
   real :: dsol, dsen, drel, diff_wind_speed
   real :: cur_wind_speed
   real :: cloud_top_temperature_water, cloud_top_temperature_ice, ctt_1km
 
   integer :: istart, iend, jstart, jend, myrad
  
   logical :: ldebug, sunglint_dust_test, lowvariability_confidence_test, put_back_cloud
   logical :: ir_cloudphase_1km_watercloud, ir_cloudphase_1km_icecloud
   integer :: model_i, model_j
   integer :: uncertain_start

   integer :: na_band_used, R1R2wavelengthIdx(2), absorbingband_index, uncertainty_nonabsorbing_1621
   real:: corr_meas(set_number_of_bands), temp_meas(set_number_of_bands), alb_meas(set_albedo_bands)
   character(10) :: phase
   real    :: cloud_reflectance(2),delta_reflectance(2), albedo_holder(2) , &
              delta_transmittance(2), tauRegimeThreshold(2), unc_reflectance(2)

   integer :: start_time, end_time, cmax, crate

   integer :: cnt_sza, cnt_vza, cnt_raz, cnt_cm_switch, cnt_wspeed, cnt_scat, cnt_wspeed_only
   logical :: wind_speed_only, interp_MS, interp_SS
   real :: irw_pressure
   real :: Tc_liquid, Tc_ice, Pc_liquid, Pc_ice, Pw_liquid, Pw_ice
   real :: rad_clr(2), bt_clr(2)
   integer :: ice_near(4), liq_near(4), nearest_used(4)
   real :: RSS_ice(4), RSS_liq(4), RSS_final(4)
   
   real :: unc_tau_real, unc_tau_1621_real, unc_tau16_real, unc_tau37_real
   real :: unc_re21_real, unc_re16_real, unc_re37_real, unc_re1621_real
   real :: unc_lwp21_real, unc_lwp16_real, unc_lwp37_real, unc_lwp1621_real
   real :: unc_re22_real, unc_lwp22_real, unc_tau22_real
   
   real :: emission_pw(20), emission_Tc(20), sigma_R37_pw(20)

   real :: alt_ray_liq, alt_ray_ice, temp_pres
   real :: aod550, irw_dummy
   
   integer :: re_idx_low, re_idx_hi 
   integer :: clin ! WDR for directed processing

   real, dimension(:,:), allocatable :: aod550_store
   
    logical :: vis1km_test    ! KGM 3-4-13
   type(cloudphase) :: ir_cloudphase

   logical :: do_retrievals_liq(4), do_retrievals_ice(4)
   logical :: finalize_liq(4), finalize_ice(4)
   logical :: set_near(4)
   real :: zzz_test ! WDR just for testing

   unc_re16_real = 0 ! WDR-UIV
   unc_tau16_real = 0 ! WDR-UIV
   unc_tau_1621_real = 0 ! WDR-UIV
   unc_re1621_real = 0 ! WDR-UIV
   unc_lwp1621_real = 0 ! WDR-UIV
   unc_re22_real = 0 
   unc_lwp22_real = 0
   unc_re21_real = 0 ! WDR-UIV
   unc_lwp21_real = 0 ! WDR-UIV
   unc_tau37_real = 0 ! WDR-UIV
   unc_re37_real = 0 ! WDR-UIV
   unc_lwp37_real = 0 ! WDR-UIV
   unc_tau_real = 0 ! WDR-UIV

   re_idx_low = 0  ! WDR-UIV
   re_idx_hi = 0  ! WDR-UIV

   status = 0
 
   xdimension = size(optical_thickness_final, 1)
   ydimension = size(optical_thickness_final, 2)
 

   lastinterp_solar_zenith     = fillvalue_real
   lastinterp_sensor_zenith    = fillvalue_real
   lastinterp_relative_azimuth = fillvalue_real
   lastinterp_scat_angle = fillvalue_real
   lastinterp_wind_speed = fillvalue_real
   lastinterp_scat_angle_ss = fillvalue_real
   count_interpolations   = 0
   retrievalcount = 0
   retrieval_failcount = 0
   library_failcount = 0
   cloud_failcount  = 0
   atmoscorr_failcount = 0
 
 
   cnt_sza = 0
   cnt_vza = 0
   cnt_raz = 0
   cnt_cm_switch = 0
   cnt_wspeed = 0
   cnt_scat = 0
   cnt_wspeed_only = 0
 
   call init_science_arrays

   call init_rtm_vars
   call init_half_radii

   COX_MUNK = .false.
   last_COX_MUNK = .false. 
   wind_speed_only = .false. 
   interp_MS = .false. 
   interp_SS = .false.

   color_first_time = .true.
   last_2way_angle = fillvalue_real
   

   call set_drel(threshold_relative_azimuth, drel)

   call init_retrieval(library_taus)

#ifndef RETRIEVE
   allocate(aod550_store(xdimension, ydimension))
   aod550_store = fillvalue_real
   cloudsummary(:,:)%cloudobserved = .true.
#endif
!  WDR export the data needed for a small test
! out in this vers  call wr_ch_vars

   do i = 1, xdimension
!  do i=40, 40

! do not process lines 1 and 1354 of the data because there are issues with cloud mask, particularly for the last line
! of the data
      if (iterationX == 1 .and. i==1 .or. iterationX == number_of_iterationsX .and. i==xdimension ) then 
         cloudsummary(i,:)%cloudobserved = .false.
!        aod550_store(i,:) = fillvalue_real
         do j=1, ydimension
            call assign_retrieval_error(i,j)
            retrieval_failcount = retrieval_failcount + 1
         end do
         cycle
      endif

   !WDR try only the center line do j=1, ydimension
   ! (this makes the center line a cloud edge somehow and missing data goes 
   !  there, so don't do for now
   ! clin = ydimension / 2 + 1
   ! do j=clin, clin
   ! WDR switch to scn_loop_st,scn_loop_en  do j=1,ydimension
   do j=scn_loop_st, scn_loop_en
!  do j=665, 665
! WDR to see all lines, pix: print*, __FILE__, __LINE__," Doing pixel ", i, " line ", j
! WDR set up the single-point model_info at chunk point i, j
!  (fills c2_model_info)
  call fill_c2_mdl( i, j )
  zzz_test = effective_radius_22_final(231,2)
   
#ifndef RETRIEVE

   if (surface_albedo(i,j,1) >= 300 .or. (cloudmask(i,j)%sunglint == 1 .and. &
      cloudmask(i,j)%water_surface) ) then 
      aod550_store(i,j) = fillvalue_real
      cycle
   endif


#endif

   
      pixX = i
      pixY = j

      cloudsummary(i,j)%cloudmask_determined = .false.
      cloudsummary(i,j)%cloudobserved = .false.
      cloudsummary(i,j)%watercloud = .false.
      cloudsummary(i,j)%icecloud = .false.
      cloudsummary(i,j)%unknowncloud = .false.
      CSR_flag_array(i,j) = 0

      liq_near = 0
      ice_near = 0

      emission_uncertainty_pw_ice = fillvalue_real
      emission_uncertainty_pw_liq = fillvalue_real
      emission_uncertainty_Tc_ice = fillvalue_real
      emission_uncertainty_Tc_liq = fillvalue_real
      sigma_R37_PW_ice = fillvalue_real
      sigma_R37_PW_liq = fillvalue_real

! can't retrieve, sun too low or no cloud top
       if (solar_zenith_angle(i,j) > solar_zenith_threshold &
           .or. solar_zenith_angle(i,j) < 0. &
           .or. cloud_top_pressure(i,j) < 0. &
           .or. cloud_top_temperature(i,j) < 0. ) then 
         retrievalstatus = 1
         call assign_retrieval_error(i,j)
         retrieval_failcount = retrieval_failcount + 1
         ! This should get clear areas have a 'no cloud' in phase
         if (cloudmask(i,j)%cloudmask_determined) &
           cloudsummary(i,j)%cloudmask_determined = .true.
         cycle 
      endif


      scattering_angle = ScatAngle(solar_zenith_angle(i,j),sensor_zenith_angle(i,j),relative_azimuth_angle(i,j))
      lowvariability_confidence_test = .false.
      na_band_used = 0


! we have to set cloudy/not cloudy and surface type outside the cloud phase call so the retrieval actually works. 
      if (cloudmask(i,j)%cloudmask_determined) cloudsummary(i,j)%cloudmask_determined = .true.

      if (cloudmask(i,j)%confident_cloudy .or. cloudmask(i,j)%probablyclear_66) &
         cloudsummary(i,j)%cloudobserved = .true.

      if (cloudmask(i,j)%snowice_surface) then
         if (cloudmask(i,j)%land_surface) cloudsummary(i,j)%land_surface = .true.
         if (cloudmask(i,j)%water_surface) cloudsummary(i,j)%ocean_surface = .true.
         if (cloudmask(i,j)%desert_surface) cloudsummary(i,j)%desert_surface = .true.
         if (cloudmask(i,j)%coastal_surface) cloudsummary(i,j)%coastal_surface = .true.
      endif
      if (cloudmask(i,j)%land_surface .or. cloudmask(i,j)%coastal_surface .or. cloudmask(i,j)%desert_surface) then
         if (cloudsummary(i,j)%ocean_surface) then
            cloudsummary(i,j)%coastal_surface= .true.
            cloudsummary(i,j)%ocean_surface = .false.
         endif
      endif


      corr_meas = fillvalue_real



! We have a cloud, now we can attempt retrieval
      if (cloudsummary(i,j)%cloudobserved .and. cloud_top_pressure(i,j) > 0.) then



         if (iterationX == 1) then
            temp_meas = band_measurements(i, :, j)
            corr_meas = band_measurements(i, :, j)
         else
            temp_meas = band_measurements(i+1, :, j)
            corr_meas = band_measurements(i+1, :, j)
         endif

! DO_COX_MUNK is a model control flag. Set/unset this in mod06_run_settings.f90
         if (cloudsummary(i,j)%ocean_surface .and. .not. cloudsummary(i,j)%snowice_surface &
                     .and. temp_meas(2) > 0. .and. DO_COX_MUNK) then 
            COX_MUNK = .true.
         else
            COX_MUNK = .false.
         endif


         const_C = pi / ( cos(solar_zenith_angle(i,j)*d2r) * solar_constant_37)
         

#ifdef GEOS5

#ifdef MCARS
      model_i = i
      model_j = j
#else
      call get_model_idx_geos5(geos5_istart, geos5_jstart, latitude(i,j), longitude(i,j), model_i, model_j)
#endif

#else
         call get_model_idx(latitude(i,j), longitude(i,j), model_i, model_j)
#endif
         ! WDR cur_wind_speed = model_info(model_i, model_j)%wind_speed
         cur_wind_speed = c2_model_info%wind_speed

         call set_interp_controls(i,j, scattering_angle, cur_wind_speed, drel,  &
                     threshold_solar_zenith,      &
                            threshold_sensor_zenith,  &
                     wind_speed_only, interp_SS, interp_MS )

!        do k=1, model_levels
         
!           print*, model_info(model_i, model_j)%pressure_profile(k), &
!              model_info(model_i, model_j)%temp_profile(k), &
!              model_info(model_i, model_j)%mixr_profile(k)
!           print*, c2_model_inf%pressure_profile(k), &
!              c2_model_info%temp_profile(k), &
!              c2_model_info%mixr_profile(k)
         
!        end do


#ifdef RETRIEVE
            if( librarystatus == 5 ) then  ! not sure if this is best way, but
                                     ! keep trying till librarystatus is not 5
                                     ! for being outside the table WDR
              interp_MS = .true.
              interp_SS = .true.
              endif
            call libraryinterpolate(solar_zenith_angle(i,j),    &
                                  sensor_zenith_angle(i,j),   &
                                  relative_azimuth_angle(i,j), &
                          scattering_angle, &
                          cur_wind_speed, &
                          wind_speed_only, interp_MS, interp_SS, &
                          debug, &
                                  librarystatus, i, j)
            if( librarystatus == 5 ) then
              call assign_retrieval_error(i,j)
              cycle
              endif
#endif

            ctt_1km = cloud_top_temperature(i,j)   !(save 1km ctt for to pass to clouddecision)

! do the IRW retrieval
#ifndef RETRIEVE

      ! WDRcloud_top_pressure(i,j) = model_info(model_i, model_j)%&
      !pressure_profile(model_info(model_i, model_j)%surface_level-6)
      !WDR cloud_top_temperature(i,j) = model_info(model_i, model_j)%&
      !temp_profile(model_info(model_i, model_j)%surface_level-6)
      cloud_top_pressure(i,j) = &
        c2_model_info%pressure_profile(c2_model_info%surface_level-6)
      cloud_top_temperature(i,j) = &
        c2_model_info%temp_profile(c2_model_info%surface_level-6)

#else

      irw_temperature(i,j) = fillvalue_real
      !
      ! WDR 15 jun 22 route around this, seems to be only for IR bands
      if( ( c2_sensor_id /= OCI_ID ) .and. ( c2_sensor_id /= OCIs_ID ) ) then
        if (COX_MUNK) then 
          call get_rtm_parameters(platform_name, &
            int_surface_emissivity_water(1,:,1), &
            sensor_zenith_angle(i,j), solar_zenith_angle(i,j),  model_i, &
            model_j, i, j)
        else
          call get_rtm_parameters(platform_name, &
            surface_emissivity_land(i,j,:),  &
            sensor_zenith_angle(i,j), solar_zenith_angle(i,j), model_i, &
            model_j, i, j)       
        endif
        endif

!     print*, i,j, model_i, model_j
!     do k=1, model_levels
!       print*, k, model_info(model_i,model_j)%pressure_profile(k), &
!         model_info(model_i,model_j)%temp_profile(k), &
!         model_info(model_i,model_j)%mixr_profile(k)
!     end do

      if (cloud_height_method(i,j) == 6) then 
        ! retrieve regular temperature
        call retrieve_irw_temp(i,j, temp_meas(band_1100), &
          model_i, model_j, rtm_rad_atm_clr, rtm_trans_atm_clr, &
          rtm_cloud_prof, irw_temperature(i,j), irw_pressure, irw_dummy)      
        cloud_top_temperature(i,j) = irw_temperature(i,j)
        cloud_top_pressure(i,j) = irw_pressure
        !
        ! now do the stuff for uncertainty
        call retrieve_irw_temp(i,j, temp_meas(band_1100), &
          model_i, model_j, rtm_rad_atm_clr_low, rtm_trans_atm_clr_low, &
          rtm_cloud_prof_low , Tc_low_for_delta, temp_pres, irw_dummy )
        call retrieve_irw_temp(i,j, temp_meas(band_1100), &
          model_i, model_j, rtm_rad_atm_clr_high, rtm_trans_atm_clr_high, &
          rtm_cloud_prof_high, Tc_high_for_delta, temp_pres, irw_dummy )
      !
      else if (cloud_height_method(i,j) > 0 .and. &
               cloud_height_method(i,j) < 6) then 
        irw_temperature(i,j) = fillvalue_real
        !  for uncertainty we need to find the temperature that fits a 
        !  delta_P of 50 mb
        !  WDR call given_P_get_T(cloud_top_pressure(i,j)-delta_Pc, &
        !   model_info(model_i, model_j), Tc_low_for_delta)
        call given_P_get_T(cloud_top_pressure(i,j)-delta_Pc, c2_model_info, &
          Tc_low_for_delta)
        ! WDR if (cloud_top_pressure(i,j) + delta_Pc > &
        !  model_info(model_i, model_j)%Ps) then 
        if (cloud_top_pressure(i,j) + delta_Pc > c2_model_info%Ps) then 
          Tc_high_for_delta = surface_temperature(i,j)
        else 
          ! WDR call given_P_get_T(cloud_top_pressure(i,j)+delta_Pc, &
          !  model_info(model_i, model_j), Tc_high_for_delta)
          call given_P_get_T(cloud_top_pressure(i,j)+delta_Pc, c2_model_info, &
            Tc_high_for_delta)
        endif
      endif

      ! WDR if (cloud_top_pressure(i,j) < 0. .or. cloud_top_pressure(i,j) &
      !  > model_info(model_i, model_j)%Ps) then
      if (cloud_top_pressure(i,j) < 0. .or. &
        cloud_top_pressure(i,j) > c2_model_info%Ps) then
        ! WDR cloud_top_pressure(i,j) = model_info(model_i, model_j)%Ps
        cloud_top_pressure(i,j) = c2_model_info%Ps
        cloud_top_temperature(i,j) = surface_temperature(i,j)
      endif
#endif

! get above-cloud water vapor          
         ! WDR call get_above_cloud_properties(model_info(model_i,model_j)%pressure_profile(:),&   
         !                       model_info(model_i,model_j)%mixr_profile(:),  &   
         !                       model_info(model_i,model_j)%surface_level, &
         call get_above_cloud_properties(c2_model_info%pressure_profile(:),&   
                                 c2_model_info%mixr_profile(:),  &   
                                 c2_model_info%surface_level, &
                                 cloud_top_pressure(i,j),   &
                                 abovecloud_watervapor(i,j), &
                                 status )
   
         Tc_liquid = fillvalue_real
         Tc_ice = fillvalue_real




! do atmospheric correction
         
         ! WDR call atmospheric_correction(i,j, iterationX, corr_meas, model_info(model_i,model_j), &
         call atmospheric_correction(i,j, iterationX, corr_meas, c2_model_info, &
                        debug, atmoscorrstatus)


! now we need to compute derivatives that we'll hang on to in the uncertainty calculations
         if( ( c2_sensor_id == OCI_ID ) .or. ( c2_sensor_id == OCIS_ID ) )then
           Bprime_Tc = 0
           Bprime_Ts = 0
         else
           Bprime_Tc = &
           ( modis_planck(platform_name, Tc_high_for_delta, channel_37um, 1) &
           - modis_planck(platform_name, Tc_low_for_delta, channel_37um, 1)) / &
           (abovecloud_watervapor(i,j)*(2.*watervapor_error)) 

           Bprime_Ts = &
             ( modis_planck(platform_name, surface_temperature(i,j) + &
             delta_Ts, channel_37um, 1) - &
             modis_planck(platform_name, surface_temperature(i,j) - &
             delta_Ts, channel_37um, 1)) / (2.0*delta_Ts) 
         endif
                           
         Transprime_1way = ( thermal_correction_oneway_high(1) - thermal_correction_oneway_low(1)) / &
                           (abovecloud_watervapor(i,j)*(2.*watervapor_error))
         Transprime_2way = ( thermal_correction_twoway_high(1) - thermal_correction_twoway_low(1)) / &
                           (abovecloud_watervapor(i,j)*(2.*watervapor_error))



         albedo_real4 = surface_albedo(i,j,:)*albedo_fac
         if (COX_MUNK) albedo_real4 = 0.
         
         RSS_ice = -999.
         RSS_liq = -999.
         rayleigh_liq = fillvalue_real
         rayleigh_ice = fillvalue_real
         alt_ray_liq = fillvalue_real
         alt_ray_ice = fillvalue_real
         liq_near = 0
         ice_near = 0

         optical_thickness_liquid = fillvalue_real
         optical_thickness_16_liquid = fillvalue_real
         optical_thickness_37_liquid = fillvalue_real
         optical_thickness_1621_liquid = fillvalue_real
         effective_radius_16_liquid = fillvalue_real
         effective_radius_21_liquid = fillvalue_real
         effective_radius_37_liquid = fillvalue_real
         effective_radius_1621_liquid = fillvalue_real

         optical_thickness_ice = fillvalue_real
         optical_thickness_16_ice = fillvalue_real
         optical_thickness_37_ice = fillvalue_real
         optical_thickness_1621_ice = fillvalue_real
         effective_radius_16_ice = fillvalue_real
         effective_radius_21_ice = fillvalue_real
         effective_radius_37_ice = fillvalue_real
         effective_radius_1621_ice = fillvalue_real

         optical_thickness_22_liquid = fillvalue_real
         effective_radius_22_liquid = fillvalue_real
         optical_thickness_22_ice = fillvalue_real
         effective_radius_22_ice = fillvalue_real
         
         nearest_used = 0
         RSS_final = fillvalue_real

#ifdef RETRIEVE
         
                     !        1.6     2.1    3.7|2.2   1.6-2.1 
         do_retrievals_liq = (/ .true., .true.,  .true., .true. /)
         do_retrievals_ice = (/ .true., .true.,  .true., .true. /)
! WDR the 3rd do_retrievals_... is for 3.7 for MODIS and 2.2 for OCI

         call corescience (i, j, cloudsummary(i,j), corr_meas, &
               Tc_liquid, Tc_ice, &
               debug, na_band_used, liq_near, ice_near, &
               RSS_liq, RSS_ice, alt_ray_liq, alt_ray_ice, &
               do_retrievals_liq, do_retrievals_ice, retrievalstatus)


         if (allocated(tau_liquid)) then 
         
            if (optical_thickness_liquid > 0.) then 
               tau_liquid(i,j) = nint(optical_thickness_liquid / retr_scale_factor)
            else 
               tau_liquid(i,j) = fillvalue_int2
            endif       

            if (optical_thickness_ice > 0.) then 
               tau_ice(i,j) = nint(optical_thickness_ice / retr_scale_factor)
            else 
               tau_ice(i,j) = fillvalue_int2
            endif       

            if (effective_radius_21_liquid > 0.) then 
               re21_liquid(i,j) = nint(effective_radius_21_liquid / retr_scale_factor)
            else 
               re21_liquid(i,j) = fillvalue_int2
            endif       

            if (effective_radius_21_ice > 0.) then 
               re21_ice(i,j) = nint(effective_radius_21_ice / retr_scale_factor)
            else 
               re21_ice(i,j) = fillvalue_int2
            endif       

         endif

            
#endif


         if (COX_MUNK) &
            call set_cox_munk_albedo (albedo_real4(:), int_reflectance_water(1,:,1))


   
         if (retrievalstatus == 0 ) retrievalcount = retrievalcount+1
         if (retrieval_failcount /= 0) retrieval_failcount = retrieval_failcount +1


! there was no cloud
      else
         ! failure before retrieval 
         retrievalstatus = 1
         call assign_retrieval_error(i,j)
         retrieval_failcount = retrieval_failcount + 1
      endif



! now we do cloud phase
      cloudsummary(i,j)%cloudmask_determined = .true.
      cloudsummary(i,j)%cloudobserved = .false.
      cloudsummary(i,j)%watercloud = .false.
      cloudsummary(i,j)%icecloud = .false.
      cloudsummary(i,j)%unknowncloud = .false.

!      set the baum phase according to the "Cloud_Phase_Infrared_1km" SDS  (to pass to cloud phase and multi-layer alg.)

        ir_cloudphase%icecloud       = 0
        ir_cloudphase%watercloud     = 0
        ir_cloudphase%unknowncloud   = 1
      if (cloud_phase_infrared(i,j) == 1) ir_cloudphase%watercloud = 1
      if (cloud_phase_infrared(i,j) == 2) ir_cloudphase%icecloud = 1
      if (cloud_phase_infrared(i,j) == 1 .or. cloud_phase_infrared(i,j) == 2) ir_cloudphase%unknowncloud = 0

      call clouddecision(platform_name,                       &
        cloudmask(i,j),                     &
        corr_meas,                          &
        RSS_liq,                            &
        RSS_ice,                            &
        optical_thickness_liquid,           & 
        optical_thickness_ice,              & 
        effective_radius_16_liquid,      & 
        effective_radius_21_liquid,      & 
        effective_radius_37_liquid,      & 
        effective_radius_16_ice,          & 
        effective_radius_21_ice,          & 
        effective_radius_37_ice,          &                         
        ctt_1km,                            &
        cloud_mask_SPI(2,i,j)*0.01,              &
        cloud_height_method(i,j),           &
        ir_cloudphase,                      &
        processing_information(i,j)%band_used_for_optical_thickness, &
        cloudsummary(i,j),                  &
        ldebug,                             &
        cloudstatus, i,j) 

!  Normally, FORCE_ICE, FORCE_WATER are false and the cloud decision 
!  above sets the cloudsummary water or ice state
! force phase to ice ( Yes, for some reason, they set the contrary state and 
!  set the right state for cloudsummary(i,j)%icecloud, watercloud)
      if (FORCE_ICE .and. cloudsummary(i,j)%cloudobserved) then 
         cloudsummary(i,j)%watercloud = .false.
         cloudsummary(i,j)%icecloud = .false.
         cloudsummary(i,j)%unknowncloud = .false.
         cloudsummary(i,j)%icecloud = .true.
      endif
! force phase to water
      if (FORCE_WATER .and. cloudsummary(i,j)%cloudobserved) then 
         cloudsummary(i,j)%watercloud = .false.
         cloudsummary(i,j)%icecloud = .false.
         cloudsummary(i,j)%unknowncloud = .false.
         cloudsummary(i,j)%watercloud = .true.
      endif

#ifndef RETRIEVE
      cloudsummary(i,j)%cloudobserved = .true.
#endif

      if (cloudsummary(i,j)%cloudobserved) then 
        ! the channels get set regardless of phase, however
        ! 0.65um gets overwritten if rayleigh is applied later
        atm_corr_refl(band_0065,i,j) = corr_meas(band_0065) 
        atm_corr_refl(band_0086,i,j) = corr_meas(band_0086) 
        atm_corr_refl(band_0124,i,j) = corr_meas(band_0124)
        atm_corr_refl(band_0163,i,j) = corr_meas(band_0163)
        atm_corr_refl(band_0213,i,j) = corr_meas(band_0213)
        if( ( c2_sensor_id == OCI_ID ) .or. ( c2_sensor_id == OCIS_ID ) ) then
          atm_corr_refl(band_0226-1,i,j) = corr_meas(band_0226)
        else
          atm_corr_refl(band_0370-1,i,j) = fillvalue_real
        endif



! set the answers based on final phase decision and do the remaining science here
! we need to compute water path, multilayer and uncertainty and set the tau_out_of_bounds QA bit


         if (cloudsummary(i,j)%watercloud .or. cloudsummary(i,j)%unknowncloud) then 
            
! set the liquid water answers
            optical_thickness_final(i,j) = optical_thickness_liquid
            optical_thickness_1621_final(i,j) = optical_thickness_1621_liquid
            optical_thickness_16_final(i,j) = optical_thickness_16_liquid
            optical_thickness_37_final(i,j) = optical_thickness_37_liquid
            effective_radius_16_final(i,j) = effective_radius_16_liquid
            effective_radius_21_final(i,j) = effective_radius_21_liquid
            effective_radius_37_final(i,j) =  effective_radius_37_liquid
            effective_radius_1621_final(i,j) = effective_radius_1621_liquid
            if (Tc_liquid > 0.) then 
               irw_temperature(i,j) = Tc_liquid
               cloud_top_temperature(i,j) = Tc_liquid
            endif
            if( ( c2_sensor_id == OCI_ID  ) .or. &
                ( c2_sensor_id == OCIS_ID ) ) then
              optical_thickness_22_final(i,j) = optical_thickness_22_liquid
              effective_radius_22_final(i,j) = effective_radius_22_liquid
              endif

            nearest_used = liq_near
            RSS_final = RSS_liq
            emission_pw = emission_uncertainty_pw_liq
            emission_Tc = emission_uncertainty_Tc_liq 
            sigma_R37_pw = sigma_R37_PW_liq

! set the rayleigh refl here
            if (.not. COX_MUNK) then
               if (alt_ray_liq > 0.) then 
                  atm_corr_refl(band_0065, i,j) = alt_ray_liq
               else
                  call bisectionsearch(water_radii, effective_radius_21_liquid, &
                              re_idx_low,re_idx_hi)
                  if (rayleigh_liq(re_idx_low) > 0. .and. rayleigh_liq(re_idx_hi) > 0.) then 
                     atm_corr_refl(band_0065, i,j) = &
                        linearinterpolation( (/water_radii(re_idx_low), water_radii(re_idx_hi) /), &
                              (/rayleigh_liq(re_idx_low), rayleigh_liq(re_idx_hi)/), &
                                       effective_radius_21_liquid)   
                  endif
               endif          
            endif
! set the tau out of bounds bit

            if (optical_thickness_final(i,j) > 150.) then!{
               optical_thickness_final(i,j) = 150.
               processing_information(i,j)%optical_thickness_outofbounds = 2
            else!}{
               processing_information(i,j)%optical_thickness_outofbounds = 0
            endif!}


            if (optical_thickness_16_final(i,j) > 150.) &
                  optical_thickness_16_final(i,j) = 150.

            if( ( c2_sensor_id == OCI_ID ) .or. &
                ( c2_sensor_id == OCIS_ID ) )then
              if( ( optical_thickness_22_final(i,j) > 150.) ) &
                optical_thickness_22_final(i,j) = 150.
            else
              if (optical_thickness_37_final(i,j) > 150.) &
                optical_thickness_37_final(i,j) = 150.
            endif

            if (optical_thickness_1621_final(i,j) > 150.) &
                  optical_thickness_1621_final(i,j) = 150.

            finalize_liq = .false.
            finalize_ice = .false. 
            if (optical_thickness_16_final(i,j) > 0. .and. &
                effective_radius_16_final(i,j) > 0.) &
                  finalize_liq(1) = .true.

            if (optical_thickness_final(i,j) > 0. .and. &
              effective_radius_21_final(i,j) > 0.) &
              finalize_liq(2) = .true.

            if( ( c2_sensor_id == OCI_ID ) .or. &
                ( c2_sensor_id == OCIS_ID ) )then
              if( optical_thickness_22_final(i,j) > 0. .and. &
              effective_radius_22_final(i,j) > 0.) &
              finalize_liq(3) = .true.
            else
              if (optical_thickness_37_final(i,j) > 0. .and. &
                effective_radius_37_final(i,j) > 0.) &
                finalize_liq(3) = .true.
            endif

            if (optical_thickness_1621_final(i,j) > 0. .and. &
              effective_radius_1621_final(i,j) > 0.) &
              finalize_liq(4) = .true.
            
         else if(cloudsummary(i,j)%icecloud) then 

! set the ice cloud answers      
            optical_thickness_final(i,j) = optical_thickness_ice
            optical_thickness_1621_final(i,j) = optical_thickness_1621_ice
            optical_thickness_16_final(i,j) = optical_thickness_16_ice
            effective_radius_16_final(i,j) = effective_radius_16_ice
            effective_radius_21_final(i,j) = effective_radius_21_ice
            effective_radius_1621_final(i,j) = effective_radius_1621_ice

            if( ( c2_sensor_id == OCI_ID ) .or. &
                ( c2_sensor_id == OCIS_ID ) )then
              optical_thickness_22_final(i,j) = optical_thickness_22_ice
              effective_radius_22_final(i,j) =  effective_radius_22_ice
            else
              optical_thickness_37_final(i,j) = optical_thickness_37_ice
              effective_radius_37_final(i,j) =  effective_radius_37_ice
            endif

            if (Tc_ice > 0.) then 
               irw_temperature(i,j) = Tc_ice
               cloud_top_temperature(i,j) = Tc_ice
            endif

            nearest_used = ice_near
            RSS_final = RSS_ice
            emission_pw = emission_uncertainty_pw_ice
            emission_Tc = emission_uncertainty_Tc_ice 
            sigma_R37_pw = sigma_R37_PW_ice
            
            if (.not. COX_MUNK) then
               if (alt_ray_ice > 0.) then 
                  atm_corr_refl(band_0065, i,j) = alt_ray_ice
               else
                  call bisectionsearch(ice_radii,effective_radius_21_ice, &
                           re_idx_low,re_idx_hi)
                  if (rayleigh_ice(re_idx_low) > 0. .and. rayleigh_ice(re_idx_hi) > 0.) then 
                     atm_corr_refl(band_0065, i,j) = &
                        linearinterpolation( (/ice_radii(re_idx_low), ice_radii(re_idx_hi) /), &
                           (/rayleigh_ice(re_idx_low), rayleigh_ice(re_idx_hi)/), &
                                       effective_radius_21_ice)   
                  endif
               endif          
            endif

! set the tau out of bounds bit
! the new setting indicates flagging of tau only if it's more than 150. All others are considered perfectly valid

            if (optical_thickness_final(i,j) > 150.) then!{
               optical_thickness_final(i,j) = 150.
               processing_information(i,j)%optical_thickness_outofbounds = 2
            else!}{
               processing_information(i,j)%optical_thickness_outofbounds = 0
            endif!}


            if (optical_thickness_16_final(i,j) > 150.) &
                  optical_thickness_16_final(i,j) = 150.

            if( ( c2_sensor_id == OCI_ID ) .or. &
                ( c2_sensor_id == OCIS_ID ) ) then
              if( optical_thickness_22_final(i,j) > 150. ) &
                optical_thickness_22_final(i,j) = 150.
            else
              if( optical_thickness_37_final(i,j) > 150. ) &
                optical_thickness_37_final(i,j) = 150.
            endif

            if (optical_thickness_1621_final(i,j) > 150.) &
                  optical_thickness_1621_final(i,j) = 150.


            finalize_liq = .false.
            finalize_ice = .false. 
            if (optical_thickness_16_final(i,j) > 0. .and. effective_radius_16_final(i,j) > 0.) &
                  finalize_ice(1) = .true.
            if (optical_thickness_final(i,j) > 0. .and. effective_radius_21_final(i,j) > 0.) &
                  finalize_ice(2) = .true.

            if( ( c2_sensor_id == OCI_ID ) .or. &
                ( c2_sensor_id == OCIS_ID ) )then
              if( optical_thickness_22_final(i,j) > 0. .and. &
                effective_radius_22_final(i,j) > 0.) &
                finalize_ice(3) = .true.
            else
              if (optical_thickness_37_final(i,j) > 0. .and. &
                effective_radius_37_final(i,j) > 0.) &
                finalize_ice(3) = .true.
            endif

            if (optical_thickness_1621_final(i,j) > 0. .and. effective_radius_1621_final(i,j) > 0.) &
                  finalize_ice(4) = .true.
         
         endif

         call set_water_path_answers(i,j, finalize_liq, finalize_ice)

         if (optical_thickness_final(i,j) > 0. .and. effective_radius_21_final(i,j) > 0.) then
            IM_successful_retrieval_count = IM_successful_retrieval_count + 1
         endif

! assign the failure metric here: 

! nearest_used and RSS_final
      set_near = .false.
      if (nearest_used(re16) == 1 .and. &
        (.not. (optical_thickness_16_final(i,j) == MAX_TAU_RTRIEVED &
        .and. effective_radius_16_final(i,j) /= fillvalue_real))) then 
         set_near(re16) = .true.
      endif

      if (nearest_used(re21) == 1 .and. &
        (.not. (optical_thickness_final(i,j) == MAX_TAU_RTRIEVED &
        .and. effective_radius_21_final(i,j) /= fillvalue_real))) then 
         set_near(re21) = .true.
      endif

      if( ( c2_sensor_id == OCI_ID ) .or. ( c2_sensor_id == OCIS_ID ) ) then
        if( nearest_used(re22) == 1 .and. &
          (.not. (optical_thickness_22_final(i,j) == MAX_TAU_RTRIEVED &
          .and. effective_radius_22_final(i,j) /= fillvalue_real))) then
          set_near(re22) = .true.
        endif
      else
        if (nearest_used(re37) == 1 .and. &
          (.not. (optical_thickness_37_final(i,j) == MAX_TAU_RTRIEVED &
          .and. effective_radius_37_final(i,j) /= fillvalue_real))) then 
           set_near(re37) = .true.
        endif
      endif

      if (nearest_used(re1621) == 1) then 
         set_near(re1621) = .true. 
      endif

      call set_failure_answers(i,j,RSS_final, set_near)

! compute multilayer
         if ( (optical_thickness_final(i,j) > 0.) .and. &
              ( c2_cmp_there(band_0935) == 1) .and. &
              ( c2_cmp_there(band_1100) == 1) ) then 
            
               call Compute_Multilayer_Map(platform_name,  &
                           transmit_correction_table,             &   
                           temp_meas,   &   
                           cloudsummary(i,j),           &   
                           ir_cloudphase,          &   
                           ! WDR model_info(model_i,model_j)%pressure_profile,&   
                           ! model_info(model_i,model_j)%mixr_profile(:),  &   
                           ! model_info(model_i,model_j)%temp_profile(:),   &
                           ! model_info(model_i,model_j)%surface_level, &
                           c2_model_info%pressure_profile,&   
                           c2_model_info%mixr_profile(:),  &   
                           c2_model_info%temp_profile(:),   &
                           c2_model_info%surface_level, &
                           cloud_top_pressure(i,j),     &   
                           abovecloud_watervapor(i,j),  &   
                           sensor_zenith_angle(i,j),                   &   
                           solar_zenith_angle(i,j),                    &  
                           relative_azimuth_angle(i,j), &
                           optical_thickness_final(i,j),&
                           optical_thickness_1621_final(i,j), &
                           i, j,           &   
                           cloud_layer_flag(i,j), ml_test_flag(i,j))              
         
         else
            cloud_layer_flag(i,j) = 0
            ml_test_flag(i,j) = 0
         endif





! *** ATTENTION **** 
! to do the 3.7um uncertainty, it is not enough to replace the re and tau with the 3.7um values. It is also 
! necessary to set the absorbingband_index to be 3.7um to feed the libraries in. In addition to that
! the albedo_holder and R1R2wavelengthIdx arrays MUST be fed with absorbingband_index-1 !!
! If you fail to do so, you will end up with a royal mess and not know why the numbers don't make sense. 
! -- G. Wind 7.5.2006


!  get retrieval uncertainty estimate

            if ( cloudsummary(i,j)%icecloud ) then!{
               phase = 'ice'
            else!}{
               phase = 'water'
            endif!}



!        if ((nearest_used(re21) == 0 .or. (nearest_used(re21) == 1 .and. optical_thickness_final(i,j) == MAX_TAU_RTRIEVED ))&
!         .and. (optical_thickness_final(i,j) .ge. 0.01) .and. (effective_radius_21_final(i,j) .ge. 0.01) .and. &
!           (cloudsummary(i,j)%icecloud .or. cloudsummary(i,j)%watercloud .or. cloudsummary(i,j)%unknowncloud)) then!{

       ! WDR - I'm keeping original 'if' decision logic above but 
       !  cleaning it up and adding a test for na_band_used
         if( ( nearest_used(re21) == 0 .or. ( nearest_used(re21) == 1 &
               .and. optical_thickness_final(i,j) == MAX_TAU_RTRIEVED ) ) &
             .and. ( optical_thickness_final(i,j) .ge. 0.01 ) &
             .and. ( effective_radius_21_final(i,j) .ge. 0.01 ) &
             .and. ( cloudsummary(i,j)%icecloud &
                     .or. cloudsummary(i,j)%watercloud &
                     .or. cloudsummary(i,j)%unknowncloud ) &
             .and. ( na_band_used > 0 ) ) then!{
 
 
            absorbingband_index = band_0213
 
!if ( na_band_used <= 0 )  then
!  print*, __FILE__, __LINE__," WDR BAD condition, na_band_used = ", &
!    na_band_used
!endif
            albedo_holder =  (/albedo_real4(na_band_used), &
                     albedo_real4(absorbingband_index)/)
            cloud_reflectance = (/corr_meas(na_band_used), &
                           corr_meas(absorbingband_index)/) 
            R1R2wavelengthIdx = (/na_band_used, absorbingband_index/)
         
            if (iterationX == 1) then 
               uncertain_start = i
            else 
               uncertain_start   = i+1
            endif

            unc_reflectance(1) = spec_uncertain(na_band_used) * &
                        exp (band_uncertainty(uncertain_start,na_band_used, j)*1.0 / uncertain_sf(na_band_used)) * 0.01
            unc_reflectance(2) = spec_uncertain(absorbingband_index) * &
                        exp (band_uncertainty(uncertain_start,absorbingband_index, j)*1.0 / uncertain_sf(absorbingband_index)) * 0.01


            if (set_bands(na_band_used) < set_bands(band_0124)) unc_reflectance(1) = max (VNIR_error, unc_reflectance(1))
            if (set_bands(na_band_used) >= set_bands(band_0124)) unc_reflectance(1) = max (SWIR_error, unc_reflectance(1))
            if (set_bands(absorbingband_index) >= set_bands(band_0124)) unc_reflectance(2) = max(SWIR_error, unc_reflectance(2))



! FIVE PERCENT       
!           unc_reflectance = 0.05
! UNC_REFL + 2%
!           unc_reflectance = sqrt ( unc_reflectance**2 + 0.02**2 )

            call getuncertainties(  optical_thickness_final(i,j), &
                           effective_radius_21_final(i,j),         &
                           liquid_water_path(i,j), &
                           phase,                     & 
                           R1R2wavelengthIdx,         &
                           unc_reflectance, &
                           albedo_holder,             &
                           transmittance_twoway(na_band_used), &
                           transmittance_twoway(absorbingband_index), &
                           meandelta_trans(na_band_used), &
                           meandelta_trans(absorbingband_index), &
                           transmittance_stddev(na_band_used), &
                           transmittance_stddev(absorbingband_index), &
                           emission_pw, emission_Tc, sigma_R37_pw, &
                           unc_tau_real ,    &
                           unc_re21_real,   &
                           unc_lwp21_real, i, j)

            optical_thickness_error(i, j) = nint(unc_tau_real / unc_scale_factor)
            effective_radius_21_error(i, j) = nint(unc_re21_real / unc_scale_factor)
            liquid_water_path_error(i,j)      = nint(unc_lwp21_real / unc_scale_factor)
                     

         else!}{
            optical_thickness_error(i, j) = fillvalue_int2
            effective_radius_21_error(i, j) = fillvalue_int2
            liquid_water_path_error(i,j)      = fillvalue_int2
         endif!}

         if ( unc_tau_real .lt. epsilon(unc_tau_real) .or.  &
            unc_re21_real .lt. epsilon(unc_re21_real)  .or. &
            unc_lwp21_real .lt. epsilon(unc_lwp21_real)) then!{ 
         
            optical_thickness_error(i, j) = fillvalue_int2
            effective_radius_21_error(i, j) = fillvalue_int2
            liquid_water_path_error(i,j)       = fillvalue_int2
         endif!}




! get uncertainty estimate for 1.6um retrieval
         if ((nearest_used(re16) == 0 .or. (nearest_used(re16) == 1 .and. optical_thickness_16_final(i,j) == MAX_TAU_RTRIEVED ))&
            .and. (optical_thickness_16_final(i,j) .ge. 0.01) .and. (effective_radius_16_final(i,j) .ge. 0.01) .and. &
            (cloudsummary(i,j)%icecloud .or. cloudsummary(i,j)%watercloud .or. cloudsummary(i,j)%unknowncloud) &
            .and. ( na_band_used > 0 ) ) then!{
 

 
            absorbingband_index = band_0163

            albedo_holder =  (/albedo_real4(na_band_used), &
                     albedo_real4(absorbingband_index)/)
            cloud_reflectance = (/corr_meas(na_band_used), &
                           corr_meas(absorbingband_index)/) 
            R1R2wavelengthIdx = (/na_band_used, absorbingband_index/)
         
            if (iterationX == 1) then 
               uncertain_start = i
            else 
               uncertain_start   = i+1
            endif

            unc_reflectance(1) = spec_uncertain(na_band_used) * &
                        exp (band_uncertainty(uncertain_start,na_band_used, j)*1.0 / uncertain_sf(na_band_used)) * 0.01
            unc_reflectance(2) = spec_uncertain(absorbingband_index) * &
                        exp (band_uncertainty(uncertain_start,absorbingband_index, j)*1.0 / uncertain_sf(absorbingband_index)) * 0.01
   
            if (set_bands(na_band_used) < set_bands(band_0124)) unc_reflectance(1) = max (VNIR_error, unc_reflectance(1))
            if (set_bands(na_band_used) >= set_bands(band_0124)) unc_reflectance(1) = max (SWIR_error, unc_reflectance(1))
            if (set_bands(absorbingband_index) >= set_bands(band_0124)) unc_reflectance(2) = max(SWIR_error, unc_reflectance(2))

            call getuncertainties(  optical_thickness_16_final(i,j), &
                           effective_radius_16_final(i,j),         &
                           liquid_water_path_16(i,j), &
                           phase,                     & 
                           R1R2wavelengthIdx,         &
                           unc_reflectance, &
                           albedo_holder,             &
                           transmittance_twoway(na_band_used), &
                           transmittance_twoway(absorbingband_index), &
                           meandelta_trans(na_band_used), &
                           meandelta_trans(absorbingband_index), &
                           transmittance_stddev(na_band_used), &
                           transmittance_stddev(absorbingband_index), &
                           emission_pw, emission_Tc, sigma_R37_pw, &                         
                           unc_tau16_real ,    &
                           unc_re16_real,   &
                           unc_lwp16_real, i, j)
            
            optical_thickness_16_error(i, j) = nint(unc_tau16_real / unc_scale_factor)
            effective_radius_16_error(i, j) = nint(unc_re16_real / unc_scale_factor)
            liquid_water_path_16_error(i, j) = nint(unc_lwp16_real / unc_scale_factor)

         else!}{
            optical_thickness_16_error(i, j) = fillvalue_int2
            effective_radius_16_error(i, j) = fillvalue_int2
            liquid_water_path_16_error(i, j) = fillvalue_int2
         endif!}

         if ( unc_tau16_real .lt. epsilon(unc_tau16_real) .or.  &
            unc_re16_real .lt. epsilon(unc_re16_real) .or. &
            unc_lwp16_real .lt. epsilon(unc_lwp16_real) ) then!{ 
            optical_thickness_16_error(i, j) = fillvalue_int2
            effective_radius_16_error(i, j) = fillvalue_int2
            liquid_water_path_16_error(i, j) = fillvalue_int2
         endif!}

    ! WDR uncertainty for 2.2? we have no re uncertainty in so why bother

   if( ( c2_sensor_id /= OCI_ID ) .and. ( c2_sensor_id /= OCIS_ID ) )then
! get uncertainty estimate for 3.7um retrieval
! this is the initial part, without any emission uncertainty
! we are not using the transmittance table to do the atmospheric correction here, so 
! for the moment uncertainty due to PW table for 3.7um is set to be 0.0 
         if ((nearest_used(re37) == 0 .or. (nearest_used(re37) == 1 .and. optical_thickness_37_final(i,j) == MAX_TAU_RTRIEVED ))&
            .and. (optical_thickness_37_final(i,j) .ge. 0.01) .and. (effective_radius_37_final(i,j) .ge. 0.01) .and. &
            (cloudsummary(i,j)%icecloud .or. cloudsummary(i,j)%watercloud .or. cloudsummary(i,j)%unknowncloud) &
            .and. ( na_band_used > 0 ) ) then!{
            

            absorbingband_index = band_0370
 

            albedo_holder =  (/albedo_real4(na_band_used), &
                     albedo_real4(absorbingband_index-1)/)
            cloud_reflectance = (/corr_meas(na_band_used), &
                           corr_meas(absorbingband_index)/) 
            R1R2wavelengthIdx = (/na_band_used, absorbingband_index-1/)
         
            if (iterationX == 1) then 
               uncertain_start = i
            else 
               uncertain_start   = i+1
            endif
         
            unc_reflectance(1) = spec_uncertain(na_band_used) * &
                        exp (band_uncertainty(uncertain_start,na_band_used, j)*1.0 / uncertain_sf(na_band_used)) * 0.01
            unc_reflectance(2) = spec_uncertain(absorbingband_index-1) * &
                        exp (band_uncertainty(uncertain_start,absorbingband_index-1, j)*1.0 / &
                           uncertain_sf(absorbingband_index-1)) * 0.01

            if (set_bands(na_band_used) < set_bands(band_0124)) unc_reflectance(1) = max (VNIR_error, unc_reflectance(1))
            if (set_bands(na_band_used) >= set_bands(band_0124)) unc_reflectance(1) = max (SWIR_error, unc_reflectance(1))
            if (set_bands(absorbingband_index) >= set_bands(band_0124) .and. set_bands(absorbingband_index) <= set_bands(band_0213)) &
                                             unc_reflectance(2) = max(SWIR_error, unc_reflectance(2))

            call getuncertainties(  optical_thickness_37_final(i,j), &
                           effective_radius_37_final(i,j),         &
                           liquid_water_path_37(i,j), &
                           phase,                     & 
                           R1R2wavelengthIdx,         &
                           unc_reflectance, &
                           albedo_holder,             &
                           transmittance_twoway(na_band_used), &
                           transmittance_twoway(absorbingband_index), &
                           meandelta_trans(na_band_used), &
                           meandelta_trans(absorbingband_index), &
                           transmittance_stddev(na_band_used), &
                           transmittance_stddev(absorbingband_index), &
                           emission_pw, emission_Tc, sigma_R37_pw,&
                           unc_tau37_real ,    &
                           unc_re37_real,   &
                           unc_lwp37_real, i, j)

            optical_thickness_37_error(i, j) = nint(unc_tau37_real / unc_scale_factor)
            effective_radius_37_error(i, j) = nint(unc_re37_real / unc_scale_factor)
            liquid_water_path_37_error(i, j) = nint(unc_lwp37_real / unc_scale_factor)

         else!}{
            optical_thickness_37_error(i, j) = fillvalue_int2
            effective_radius_37_error(i, j) = fillvalue_int2
            liquid_water_path_37_error(i, j) = fillvalue_int2
         endif!}

         if ( unc_tau37_real .lt. epsilon(unc_tau37_real) .or.  &
             unc_re37_real .lt. epsilon(unc_re37_real) .or. &
            unc_lwp37_real .lt. epsilon(unc_lwp37_real) ) then!{ 
            optical_thickness_37_error(i, j) = fillvalue_int2
            effective_radius_37_error(i, j) = fillvalue_int2
            liquid_water_path_37_error(i, j) = fillvalue_int2
         endif!}
       endif !  WDR guard this in the OCI case


!  get 1621 retrieval uncertainty estimate
         if ((nearest_used(re1621) == 0 .or. (nearest_used(re1621) == 1 .and. optical_thickness_1621_final(i,j) == MAX_TAU_RTRIEVED ))&
            .and.  (optical_thickness_1621_final(i,j) .ge. 0.01) .and. (effective_radius_1621_final(i,j) .ge. 0.01) .and. &
            (cloudsummary(i,j)%icecloud .or. cloudsummary(i,j)%watercloud .or. cloudsummary(i,j)%unknowncloud) &
            .and. ( na_band_used > 0 ) ) then!{
      

            uncertainty_nonabsorbing_1621 = band_0163
            absorbingband_index = band_0213


            albedo_holder =  (/albedo_real4(uncertainty_nonabsorbing_1621), &
                        albedo_real4(absorbingband_index)/)
            cloud_reflectance = (/corr_meas(uncertainty_nonabsorbing_1621), &
                     corr_meas(absorbingband_index)/) 
            R1R2wavelengthIdx = (/uncertainty_nonabsorbing_1621, absorbingband_index/)
         
         
            if (iterationX == 1) then 
               uncertain_start = i
            else 
               uncertain_start   = i+1
            endif

            unc_reflectance(1) = spec_uncertain(uncertainty_nonabsorbing_1621) * &
                        exp (band_uncertainty(uncertain_start,uncertainty_nonabsorbing_1621, j)*1.0 / &
                        uncertain_sf(uncertainty_nonabsorbing_1621)) * 0.01
            unc_reflectance(2) = spec_uncertain(absorbingband_index) * &
                        exp (band_uncertainty(uncertain_start,absorbingband_index, j)*1.0 / uncertain_sf(absorbingband_index)) * 0.01
         
            if (set_bands(uncertainty_nonabsorbing_1621) < set_bands(band_0124)) unc_reflectance(1) = max (VNIR_error, unc_reflectance(1))
            if (set_bands(uncertainty_nonabsorbing_1621) >= set_bands(band_0124)) unc_reflectance(1) = max (SWIR_error, unc_reflectance(1))
            if (set_bands(absorbingband_index) >= set_bands(band_0124)) unc_reflectance(2) = max(SWIR_error, unc_reflectance(2))

            call getuncertainties(optical_thickness_1621_final(i,j), &
                           effective_radius_1621_final(i,j),         &
                           liquid_water_path_1621(i,j), &
                           phase,                     & 
                           R1R2wavelengthIdx,         &
                           unc_reflectance, &
                           albedo_holder,             &
                           transmittance_twoway(uncertainty_nonabsorbing_1621), &
                           transmittance_twoway(absorbingband_index), &
                           meandelta_trans(uncertainty_nonabsorbing_1621), &
                           meandelta_trans(absorbingband_index), &
                           transmittance_stddev(uncertainty_nonabsorbing_1621), &
                           transmittance_stddev(absorbingband_index), &
                           emission_pw, emission_Tc, sigma_R37_pw,&
                           unc_tau_1621_real ,    &
                           unc_re1621_real,   &
                           unc_lwp1621_real, i, j)

         

            optical_thickness_1621_error(i, j) = nint(unc_tau_1621_real / unc_scale_factor)
            effective_radius_1621_error(i, j) = nint(unc_re1621_real / unc_scale_factor)
            liquid_water_path_1621_error(i,j) = nint(unc_lwp1621_real / unc_scale_factor)

         else!}{
            optical_thickness_1621_error(i, j) = fillvalue_int2
            effective_radius_1621_error(i, j) = fillvalue_int2
            liquid_water_path_1621_error(i,j) = fillvalue_int2
         endif!}

         if ( unc_tau_1621_real .lt. epsilon(unc_tau_1621_real) .or.  &
            unc_re1621_real .lt. epsilon(unc_re1621_real)  .or. &
            unc_lwp1621_real .lt. epsilon(unc_lwp1621_real)) then!{ 

            optical_thickness_1621_error(i, j) = fillvalue_int2
            effective_radius_1621_error(i, j) = fillvalue_int2
            liquid_water_path_1621_error(i,j) = fillvalue_int2
         endif!}

      ! WDR HA! the 2.2 uncertainty would go here.  Seeing as 
      ! we could naver make uncertainty (no refl uncertainty) I'll
      ! leave this as a placeholder for when it can be done

         if (DO_CSR) then 
   
! let us remember that band measurements are being overscanned
               if (iterationX == 1) then 
                  if (i==1) then 
                     istart = i
                     iend = i+1
                  endif
               else
                  if (i==1) then 
                     istart = i
                     iend = i+2
                  endif
               endif
   
               if (i > 1 .and. i < xdimension) then 
                  istart = i-1
                  iend = i+1
               endif
               
               if (iterationX < number_of_iterationsX) then 
                   if (i == xdimension) then 
                     istart = i-1
                     iend = i+1
                   endif
               else 
                  if (i == xdimension) then 
                     istart = i-1
                     iend = i
                  endif
               endif
               
               
               if (j == 1) then 
                  jstart = j
                  jend = j+1
               endif
               if (j >=2 .and. j <= (ydimension-1)) then 
                  jstart = j-1
                  jend = j+1
               endif
               if (j == ydimension) then 
                  jstart = j-1
                  jend = j
               endif

               ! Check if 1km visible reflectance threshold or VIS/NIR ratio tests are applied and cloudy. Clear sky restoral
               ! Part V (CSR=3) test will only be applied over ocean if vis1km_test = .true. (i.e., either one of visible
               ! reflectance or VIS/NIR ratio tests are applied and cloudy).
               ! KGM 3-4-13
               vis1km_test = .false.
               if ((cloudmask(i,j)%applied_visiblereflectance==1 .and. cloudmask(i,j)%test_visiblereflectance==1) .or. &
                  (cloudmask(i,j)%applied_visnirratio==1 .and. cloudmask(i,j)%test_visnirratio==1) ) vis1km_test = .true.

               call cloudiness_test (cloudmask(i,j),           &
                                  cloudsummary(i,j),        &
                                  temp_meas, &
                                  band_measurements(istart:iend,:,jstart:jend), &
                                  sunglint_dust_test,       &
                                  lowvariability_confidence_test, &
                                  CSR_flag_array(i,j), latitude(i,j), &
                                  cloud_height_method(i,j), vis1km_test)   ! KGM 3-4-13 GW 3.28.13
 
               if (CSR_flag_array(i,j) == 2 .and. cloudsummary(i,j)%ocean_surface) then 
#ifdef RETRIEVE               
                  call compute_aod(i, j, scattering_angle, corr_meas, cur_wind_speed, aod550)

                  ! if the aerosol optical depth is too much then it's probably a cloud
                  ! and we can keep the retrieval, however we will mark it as a potentially
                  ! problematic cloud. 
                  if (aod550 > 0.95) then ! aod550 is a ln(aod+0.01) quantity
                     CSR_flag_array(i,j) = 0
                  endif
#endif
               endif

         endif

#ifndef RETRIEVE        



         if (surface_albedo(i,j,1) < 300 .and. .not. cloudmask(i,j)%desert_surface .and. &
               .not. cloudsummary(i,j)%snowice_surface ) then 
         

         
            call compute_aod(i, j, scattering_angle, corr_meas, cur_wind_speed, aod550)
            aod550_store(i,j) = (exp(aod550) - 0.01)  ! that's so it can be stored in a good way in an RE sds
            if (aod550_store(i,j) < 0.) aod550_store(i,j) = 0.01

         else
            aod550_store(i,j) = fillvalue_real
         endif

#endif       
         
      else
         retrievalstatus = 1
         call assign_retrieval_error(i,j)
         retrieval_failcount = retrieval_failcount + 1
      endif

   enddo
enddo
!  end of loop over the data blob

! WDR capture the working arrays before modification
   call capture_arrays

   call cleanup_retrieval

!  print*, "NUM_INTERP:", count_interpolations

!  print*, "num_sza: ", cnt_sza
!  print*, "num_vza: ", cnt_vza
!  print*, "num_raz: ", cnt_raz
!  print*, "num_scat: ", cnt_scat
!  print*, "num_cm_switch: ", cnt_cm_switch
!  print*, "num_wspeed: ", cnt_wspeed
!  print*, "** wind speed only: ", cnt_wspeed_only

! Now that we've done the clear sky restoral, we remove the edges of the clouds (ED)
! WDR temp to set so no remove done
   if (DO_CSR) then !{

      call remove_edge_scenes(cloudsummary, &
                           CSR_flag_array, &
                           xdimension, ydimension,  &
                           status)


!   reset cloudsummary variables for pixels "cleared" by CSR or ED, this 
!   step necessary
!   so that "cleared" pixels in cloud optical properties SDS get set to 
!   clear, and  
!   so that pertainent QA can be properly identified when set   GTA 6/7/05

!     endif

      where(CSR_flag_array == 2) 
         cloudsummary%cloudobserved = .false.
         cloudsummary%watercloud = .false.
         cloudsummary%icecloud = .false.
         cloudsummary%unknowncloud = .false.
         
         optical_thickness_final = fillvalue_real
         optical_thickness_16_final = fillvalue_real
         optical_thickness_37_final = fillvalue_real
         optical_thickness_1621_final = fillvalue_real
         effective_radius_16_final = fillvalue_real
         effective_radius_21_final = fillvalue_real
         effective_radius_37_final = fillvalue_real
         effective_radius_1621_final = fillvalue_real
         liquid_water_path = fillvalue_real
         liquid_water_path_16 = fillvalue_real
         liquid_water_path_37 = fillvalue_real
         liquid_water_path_1621 = fillvalue_real
         optical_thickness_error = fillvalue_int2
         effective_radius_21_error = fillvalue_int2
         effective_radius_16_error = fillvalue_int2
         effective_radius_37_error = fillvalue_int2
         liquid_water_path_error = fillvalue_int2
         liquid_water_path_16_error = fillvalue_int2
         liquid_water_path_37_error = fillvalue_int2
         cloud_layer_flag = 0
         ml_test_flag = 0
         optical_thickness_1621_error = fillvalue_int2
         optical_thickness_16_error = fillvalue_int2
         optical_thickness_37_error = fillvalue_int2
         effective_radius_1621_error = fillvalue_int2
         liquid_water_path_1621_error = fillvalue_int2
         
         failure_metric%tau = fillvalue_int2  
         failure_metric%re = fillvalue_int2  
         failure_metric%RSS = fillvalue_int2  

         failure_metric_16%tau = fillvalue_int2  
         failure_metric_16%re = fillvalue_int2  
         failure_metric_16%RSS = fillvalue_int2  

         failure_metric_1621%tau = fillvalue_int2  
         failure_metric_1621%re = fillvalue_int2  
         failure_metric_1621%RSS = fillvalue_int2  

         failure_metric_37%tau = fillvalue_int2  
         failure_metric_37%re = fillvalue_int2  
         failure_metric_37%RSS = fillvalue_int2  
         
         irw_temperature = fillvalue_real
         precip_water_094 = fillvalue_real
       end where

       if( ( c2_sensor_id == OCI_ID ) .or. ( c2_sensor_id == OCIS_ID ) )then
         where(CSR_flag_array == 2) 
           optical_thickness_22_final = fillvalue_real
           effective_radius_22_final = fillvalue_real
           liquid_water_path_22 = fillvalue_real
           effective_radius_22_error = fillvalue_int2
           liquid_water_path_22_error = fillvalue_int2
           optical_thickness_22_error = fillvalue_int2
           failure_metric_22%tau = fillvalue_int2
           failure_metric_22%re = fillvalue_int2
           failure_metric_22%RSS = fillvalue_int2
         end where
       endif

! Now split off the PCL retrievals
! if it's an edge pixel or a 250m variable pixels then 
! split off the retrieval to the PCL storage and get rid of the value in 
! main retrieval arrays
      call split_PCL(xdimension, ydimension)
                 
   endif !}




! now we need to compute the inventory metadata that relates to the cloudiness percentage
! water cloud percentage and ice cloud percentage. We need to aggregate the little suckers

   do i=1, xdimension
      do j = 1, ydimension

! first of all we need to make sure that we have a cloud
         if (solar_zenith_angle(i,j) <= solar_zenith_threshold .and. &
              cloudsummary(i,j)%cloudobserved) then 
            IM_cloudy_count = IM_cloudy_count + 1
            
! now that we're sure, we can count the ice and water cloud pixels
            if (cloudsummary(i,j)%watercloud) then 
               IM_water_cloud_count = IM_water_cloud_count + 1
            endif
            if (cloudsummary(i,j)%icecloud) then 
               IM_ice_cloud_count = IM_ice_cloud_count + 1
            endif

            if (cloudsummary(i,j)%unknowncloud) then 
               IM_undet_count = IM_undet_count + 1
            endif

          endif

      end do
   end do


!  optical_thickness_final = abovecloud_watervapor
!  effective_radius_16_final = cloud_top_pressure

#ifndef RETRIEVE
   effective_radius_21_final(:,:) = aod550_store(:,:)
   deallocate(aod550_store)
#endif
   
!  print*, abovecloud_watervapor(14,884), cloud_top_pressure(14,884), surface_temperature(14,884)

! WDR complete the processing_information setup with call to set_quality_data
  call set_quality_data( xdimension, ydimension )
 
!  print*, optical_thickness_final(19, 1992)
!  print*,     effective_radius_16_final (19, 1992)
!  print*,     effective_radius_21_final (19, 1992)
!  print*,     effective_radius_37_final(19, 1992)
   
 end subroutine scienceinterface

subroutine compute_aod(x, y, scat_ang, corr_meas, ws, aod550)

   use core_arrays
   use mod06_run_settings
   use science_parameters, only: d2r
   use nonscience_parameters
   use modis_numerical_module, only: linearinterpolation
   
   integer, intent(in) :: x, y
   real, intent(in) :: scat_ang, ws
   real, intent(inout) :: aod550
   real, dimension(:), intent(in) :: corr_meas

   external ffnet_terra, ffnet_aqua, ffnet_aqua_land, ffnet_terra_land

   real*8 :: input(15), output(1), input_land(14)
   real*8  :: ga, sza, vza, saz, vaz, raz, sca
   real :: check_val, temp_16
   integer :: i
   
   sza = solar_zenith_angle(x,y)
   vza = sensor_zenith_angle(x,y)
   saz = solar_azimuth_angle(x,y)
   vaz = sensor_azimuth_angle(x,y)
   raz = relative_azimuth_angle(x,y)
   
   ga = cos(sza*d2r) * cos(vza*d2r) + sin(sza*d2r) * sin(vza*d2r) *  cos(raz*d2r)
   
   sca = cos(scat_ang*d2r)*1.0d0
   saz = cos(saz*d2r)*1.0d0
   sza = cos(sza*d2r)*1.0d0
   vaz = cos(vaz*d2r)*1.0d0
   vza = cos(vza*d2r)*1.0d0
            
!tbsS'InputNames ocean'
!S'mRef470'
!aS'mRef550'
!aS'mRef660'
!aS'mRef870'
!aS'mRef1200'
!aS'mRef1600'
!aS'mRef2100'
!aS'ScatteringAngle'
!aS'GlintAngle'
!aS'SolarAzimuth'
!aS'SolarZenith'
!aS'SensorAzimuth'
!aS'SensorZenith'
!aS'cloud'
!aS'wind'
   
   
   if (cloudsummary(x,y)%ocean_surface) then 

      input = (/  corr_meas(band_0047)*1.0d0, &
               corr_meas(band_0055)*1.0d0, &
               corr_meas(band_0065)*1.0d0, &
               corr_meas(band_0086)*1.0d0, &
               corr_meas(band_0124)*1.0d0, &
               corr_meas(band_0163)*1.0d0, &
               corr_meas(band_0213)*1.0d0, &
               sca, ga, saz, sza, vaz, vza, 0.0d0, &
               ws*1.0d0 /)
   
      if ( corr_meas(band_0163) < 0.) then 
   
         temp_16 = linearinterpolation ( (/1.24, 2.13/), (/corr_meas(band_0124), corr_meas(band_0213)/), 1.63)
         input(6) = temp_16*1.0d0
         
      endif

      if (platform_name(1:4) == "Aqua") call ffnet_aqua(input, output)
      if (platform_name(1:5) == "Terra") call ffnet_terra(input, output)
   else


      input_land = (/   corr_meas(band_0055)*1.0d0, &
               corr_meas(band_0047)*1.0d0, &
               corr_meas(band_0065)*1.0d0, &
               corr_meas(band_0086)*1.0d0, &
               corr_meas(band_0124)*1.0d0, &
               corr_meas(band_0163)*1.0d0, &
               corr_meas(band_0213)*1.0d0, &
               sca, saz, sza, vaz, vza, 0.0d0, albedo_real4(band_0065)*1.0d0 /) ! this is really a 550nm albedo
                                                                ! look at ancillary module and see
   
      if ( corr_meas(band_0163) < 0.) then 
   
         temp_16 = linearinterpolation ( (/1.24, 2.13/), (/corr_meas(band_0124), corr_meas(band_0213)/), 1.63)
         input(6) = temp_16*1.0d0
         
      endif

      if (platform_name(1:4) == "Aqua") call ffnet_aqua_land(input_land, output)
      if (platform_name(1:5) == "Terra") call ffnet_terra_land(input_land, output)
   endif 

   aod550 = output(1) ! ln(aod+0.01)

end subroutine compute_aod

subroutine fill_c2_mdl( i, j )
  !
  !  WDR fill_c2_mdl will make a 'ancillary_type' structure like model_info
  !  has but for just 1 point of interest
  !
  !  i, j   pixel, line, coordinates in the current granule chunk
  !
  !  WDR 27 Nov 2018
  !
  use core_arrays, only: c2_model_info
  use ch_xfr, only: c2_prof_mixr, c2_prof_t, c2_prof_p, c2_prof_hgt, &
    c2_tsfc, c2_psfc, c2_wind, c2_tot_o3, c2_ice_frac, c2_snow_frac, &
    c2_alt_o3, c2_alt_icec, c2_sfc_lvl, c2_trop_lvl, c2_alt_snowice

  integer, intent(in) :: i, j

  c2_model_info%mixr_profile = c2_prof_mixr( i, j, : )
  c2_model_info%temp_profile = c2_prof_t( i, j, : )
  c2_model_info%height_profile = c2_prof_hgt( i, j, : )
  c2_model_info%pressure_profile = c2_prof_p( i, j, : )
  ! c2_model_info%o3_profile - not used, I think
  c2_model_info%o3_profile = 0
  c2_model_info%Ts = c2_tsfc( i, j )
  c2_model_info%Ps = c2_psfc( i, j )
  c2_model_info%wind_speed = c2_wind( i, j )
  c2_model_info%col_o3 = c2_tot_o3( i, j )
  c2_model_info%seaice_fraction = c2_ice_frac( i, j )
  c2_model_info%snow_fraction = c2_snow_frac( i, j )
  c2_model_info%surface_level = c2_sfc_lvl( i, j )
  c2_model_info%trop_level = c2_trop_lvl( i, j )
  ! c2_model_info%LSM = Not sure if used again
  c2_model_info%LSM = 0

end subroutine fill_c2_mdl

end module modis_science_module

