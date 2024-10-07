module core_arrays

#ifdef CT_CODE


   real,   allocatable, dimension (:,:)   :: solar_zenith_angle

#else


   use GeneralAuxType
   use modis_cloudstructure
   use modis_sciencestructure
   use mod06_run_settings, only: set_albedo_bands, set_number_of_bands

   use science_parameters
   use global_model_grids

   implicit none

   ! core algorithm arrays
   
	real :: optical_thickness_liquid, optical_thickness_ice
	real :: optical_thickness_16_liquid, optical_thickness_16_ice
	real :: optical_thickness_37_liquid, optical_thickness_37_ice
	real :: optical_thickness_1621_liquid, optical_thickness_1621_ice
	real :: effective_radius_21_liquid, effective_radius_21_ice
	real :: effective_radius_22_liquid, effective_radius_22_ice
   real :: optical_thickness_22_ice, optical_thickness_22_liquid
	real :: effective_radius_16_liquid, effective_radius_16_ice
	real :: effective_radius_1621_liquid, effective_radius_1621_ice
	real :: effective_radius_37_liquid, effective_radius_37_ice

	integer*2, dimension(:,:), allocatable :: tau_liquid, tau_ice, re21_liquid
   integer*2, dimension(:,:), allocatable ::  re21_ice, re22_liquid, &
      re22_ice

   real(single),   allocatable, dimension (:,:)   :: optical_thickness_final
   real(single),   allocatable, dimension (:,:)   :: optical_thickness_22_final
   real(single),   allocatable, dimension (:,:)   :: optical_thickness_16_final
   real(single),   allocatable, dimension (:,:)   :: optical_thickness_37_final
   real(single),   allocatable, dimension (:,:)   :: optical_thickness_1621_final
   real(single),   allocatable, dimension (:,:)   :: effective_radius_16_final
   real(single),   allocatable, dimension (:,:)   :: effective_radius_21_final
   real(single),   allocatable, dimension (:,:)   :: effective_radius_22_final
   real(single),   allocatable, dimension (:,:)   :: effective_radius_37_final
   real(single),   allocatable, dimension (:,:)   :: effective_radius_1621_final

   integer*2,   allocatable, dimension (:,:)   :: optical_thickness_final_PCL
   integer*2,   allocatable, dimension (:,:)   :: optical_thickness_22_final_PCL
   integer*2,   allocatable, dimension (:,:)   :: optical_thickness_16_final_PCL
   integer*2,   allocatable, dimension (:,:)   :: optical_thickness_37_final_PCL
   integer*2,   allocatable, dimension (:,:)   :: optical_thickness_1621_final_PCL
   integer*2,   allocatable, dimension (:,:)   :: effective_radius_16_final_PCL
   integer*2,   allocatable, dimension (:,:)   :: effective_radius_21_final_PCL
   integer*2,   allocatable, dimension (:,:)   :: effective_radius_22_final_PCL
   integer*2,   allocatable, dimension (:,:)   :: effective_radius_37_final_PCL
   integer*2,   allocatable, dimension (:,:)   :: effective_radius_1621_final_PCL
   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_PCL
   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_22_PCL
   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_16_PCL
   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_37_PCL
   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_1621_PCL

   real(single),   allocatable, dimension (:,:)   :: liquid_water_path
   real(single),   allocatable, dimension (:,:)   :: liquid_water_path_22
   real(single),   allocatable, dimension (:,:)   :: liquid_water_path_16
   real(single),   allocatable, dimension (:,:)   :: liquid_water_path_37
   real(single),   allocatable, dimension (:,:)   :: liquid_water_path_1621

! we can get away with integer*2 here because these arrays are for output storage only
! they are a final answer that is not used in any subsequent calculation at any time
   integer*2,   allocatable, dimension (:,:)   :: optical_thickness_error
   integer*2,   allocatable, dimension (:,:)   :: optical_thickness_22_error
   integer*2,   allocatable, dimension (:,:)   :: optical_thickness_16_error
   integer*2,   allocatable, dimension (:,:)   :: optical_thickness_37_error
   
   integer*2,   allocatable, dimension (:,:)   :: effective_radius_21_error
   integer*2,   allocatable, dimension (:,:)   :: effective_radius_22_error
   integer*2,   allocatable, dimension (:,:)   :: effective_radius_16_error
   integer*2,   allocatable, dimension (:,:)   :: effective_radius_37_error

   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_error
   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_22_error
   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_16_error
   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_37_error
   
   integer*2,   allocatable, dimension (:,:)   :: optical_thickness_1621_error
   integer*2,   allocatable, dimension (:,:)   :: effective_radius_1621_error
   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_1621_error

   integer(integer_onebyte),      allocatable, dimension (:,:)   :: cloud_layer_flag, ml_test_flag
   integer(integer_onebyte),      allocatable, dimension (:,:)   :: CSR_flag_array

#if MAS_INST || EMAS_INST
   real(single), allocatable, dimension(:,:) :: spatial_variability
   integer(integer_onebyte), allocatable, dimension(:,:) :: restoral_pos   
#endif

   type(processflag)  , allocatable,dimension (:,:)   :: cloudsummary

   !water vapor band cirrus retrieval arrays
   real(single),   allocatable, dimension (:,:)   :: optical_thickness_138
   real(single),   allocatable, dimension (:,:)   :: effective_radius_138
   real(single),   allocatable, dimension (:,:)   :: optical_thickness_138_error
   real(single),   allocatable, dimension (:,:)   :: effective_radius_138_error
   real(single),   allocatable, dimension (:,:,:) :: uncertainty_138_transmittance
   real(single),   allocatable, dimension (:,:,:) :: uncertainty_138_calibration
   real(single),   allocatable, dimension (:,:,:) :: uncertainty_138_effvar


   !input data arrays (measurement)
   real,   allocatable, dimension (:,:)   :: sensor_zenith_angle
   real,   allocatable, dimension (:,:)   :: solar_zenith_angle
   real,   allocatable, dimension (:,:)   :: sensor_azimuth_angle
   real,   allocatable, dimension (:,:)   :: solar_azimuth_angle

   real(single),   allocatable, dimension (:,:)   :: relative_azimuth_angle
   real(single),   allocatable, dimension (:,:,:) :: band_measurements
#ifdef SIM
   real, allocatable, dimension(:,:,:) :: band_uncertainty
   real, allocatable, dimension(:,:) :: scatter_angle
#else
   integer*1, allocatable, dimension(:,:,:) :: band_uncertainty
#endif
   type(cloudmask_type),allocatable,dimension (:,:)   :: cloudmask

   real    :: meandelta_trans(set_number_of_bands)
   real :: mean_delta_ozone, ozone_transmittance
   real    :: transmittance_twoway(set_number_of_bands)
   real    :: transmittance_stddev(set_number_of_bands)
   real    :: thermal_correction_twoway(set_number_of_bands)
   real    :: thermal_correction_oneway(set_number_of_bands)

   ! ancillary arrays
   real, parameter :: albedo_fac = 0.001
#ifdef SEVIRI_INST
 real,   allocatable, dimension (:,:,:) :: surface_albedo
#else
 integer*2,   allocatable, dimension (:,:,:) :: surface_albedo
#endif

   real :: albedo_real4(set_albedo_bands)
   
   real, dimension(:,:), allocatable :: flight_time
   integer*2, dimension(:,:,:), allocatable :: cloud_mask_SPI


   real(single),   allocatable, dimension (:,:)   :: latitude
   real(single),   allocatable, dimension (:,:)   :: longitude

   real(single),   allocatable, dimension (:,:)   :: surface_temperature
   real(single),   allocatable, dimension (:,:,:)   :: surface_emissivity_land
   
   real(single),   allocatable, dimension (:,:)   :: cloud_top_temperature, cloud_top_temperature_B
   real(single),   allocatable, dimension (:,:)   :: cloud_top_pressure
   real(single),   allocatable, dimension (:,:)   :: cloud_top_height
   real(single),   allocatable, dimension (:,:)   :: cloud_effective_emissivity
   
   real(single),   allocatable, dimension (:,:)   :: abovecloud_watervapor
   real(single),   allocatable, dimension (:,:)   :: column_ozone
	real, dimension(:,:), allocatable :: precip_water_094

   real, dimension(:,:,:), allocatable :: clear_sky_rad
   real, dimension(:,:,:), allocatable :: clear_sky_btemp


	integer*1, dimension(:,:), allocatable :: cloud_height_method
	integer*1, dimension(:,:), allocatable :: cloud_phase_infrared
	real, dimension(:,:), allocatable :: irw_temperature
	
	type(failed_type), dimension(:,:), allocatable :: failure_metric
	type(failed_type), dimension(:,:), allocatable :: failure_metric_22
	type(failed_type), dimension(:,:), allocatable :: failure_metric_16
	type(failed_type), dimension(:,:), allocatable :: failure_metric_37
	type(failed_type), dimension(:,:), allocatable :: failure_metric_1621
	
	real, dimension(:,:,:), allocatable :: atm_corr_refl
		
	   !Other
   character*15                                  :: platform_name
   
   type(ancillary_type), allocatable, dimension(:,:)  :: model_info
   ! WDR add a similar structure for just 1 point
   ! WDR MAY NOT NEED type(ancillary_type), allocatable :: c2_model_info
   type(ancillary_type) :: c2_model_info
   type(qualityanalysis), allocatable, dimension(:,:) :: processing_information
   type(qualityanalysis), allocatable, dimension(:,:) :: processing_information_sav

	integer*1, dimension(:,:), allocatable :: seviri_cloudphase
	integer*1, dimension(:,:), allocatable :: ahi_cloudphase
	integer*1, dimension(:,:), allocatable :: viirs_cloudphase
	
   integer(integer_fourbyte),allocatable, dimension (:)     :: bands
   
	type(stat_type) :: Statistics_1km
   
    real :: thermal_correction_oneway_low(2), thermal_correction_oneway_high(2)
    real :: thermal_correction_twoway_low(2), thermal_correction_twoway_high(2)
    real :: emission_uncertainty_pw_ice(20), emission_uncertainty_pw_liq(20)
    real :: emission_uncertainty_Tc_ice(20), emission_uncertainty_Tc_liq(20)
	real :: sigma_R37_PW_ice(20), sigma_R37_PW_liq(20)


    real :: Tc_low_for_delta, Tc_high_for_delta
   
   	real :: Bprime_Tc, Bprime_Ts, Transprime_1way, Transprime_2way
   	real :: const_C
   
   	real, dimension(:,:), allocatable :: snow_cover
  ! WDR we will have controls for init_science_arrays and modis_science_module
  ! to control processing or transferring of data for the 3 lines of data
  integer :: scn_loop_st, scn_loop_en, xfr_num, xfr_from(2), xfr_to(2)

! WDR - declarations for the capture or _sav arrays to preserve point-only 
! (not area-influenced) values of many arrays
   integer*2, allocatable, dimension (:,:) :: optical_thickness_final_PCL_sav
   integer*2, allocatable, dimension (:,:) :: optical_thickness_22_final_PCL_sav
   integer*2, allocatable, dimension (:,:) :: optical_thickness_1621_final_PCL_sav
   integer*2, allocatable, dimension (:,:) :: effective_radius_21_final_PCL_sav
   integer*2, allocatable, dimension (:,:) :: effective_radius_22_final_PCL_sav
   integer*2,   allocatable, dimension (:,:)   :: effective_radius_1621_final_PCL_sav
   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_PCL_sav
   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_22_PCL_sav
   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_1621_PCL_sav
   integer*2,   allocatable, dimension (:,:)   :: optical_thickness_37_final_PCL_sav
   integer*2,   allocatable, dimension (:,:)   :: optical_thickness_16_final_PCL_sav
   integer*2,   allocatable, dimension (:,:)   :: effective_radius_16_final_PCL_sav
   integer*2,   allocatable, dimension (:,:)   :: effective_radius_37_final_PCL_sav
   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_16_PCL_sav
   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_37_PCL_sav
   real, dimension(:,:), allocatable :: precip_water_094_sav
   real, dimension(:,:), allocatable :: irw_temperature_sav
   type(failed_type), dimension(:,:), allocatable :: failure_metric_37_sav
   type(failed_type), dimension(:,:), allocatable :: failure_metric_1621_sav
   type(failed_type), dimension(:,:), allocatable :: failure_metric_16_sav
   type(failed_type), dimension(:,:), allocatable :: failure_metric_sav
   type(failed_type), dimension(:,:), allocatable :: failure_metric_22_sav
   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_1621_error_sav
   integer*2,   allocatable, dimension (:,:)   :: effective_radius_1621_error_sav
   integer*2,   allocatable, dimension (:,:)   :: optical_thickness_37_error_sav
   integer*2,   allocatable, dimension (:,:)   :: optical_thickness_16_error_sav
   integer*2,   allocatable, dimension (:,:)   :: optical_thickness_1621_error_sav
   integer(integer_onebyte),      allocatable, dimension (:,:)   :: cloud_layer_flag_sav, ml_test_flag_sav
   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_37_error_sav
   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_16_error_sav
   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_error_sav
   integer*2,   allocatable, dimension (:,:)   :: liquid_water_path_22_error_sav
   integer*2,   allocatable, dimension (:,:)   :: effective_radius_37_error_sav
   integer*2,   allocatable, dimension (:,:)   :: effective_radius_16_error_sav
   integer*2,   allocatable, dimension (:,:)   :: effective_radius_21_error_sav
   integer*2,   allocatable, dimension (:,:)   :: effective_radius_22_error_sav
   integer*2,   allocatable, dimension (:,:)   :: optical_thickness_error_sav
   integer*2,   allocatable, dimension (:,:)   :: optical_thickness_22_error_sav
   real(single),   allocatable, dimension (:,:)   :: liquid_water_path_1621_sav
   real(single),   allocatable, dimension (:,:)   :: liquid_water_path_37_sav
   real(single),   allocatable, dimension (:,:)   :: liquid_water_path_16_sav
   real(single),   allocatable, dimension (:,:)   :: liquid_water_path_sav
   real(single),   allocatable, dimension (:,:)   :: liquid_water_path_22_sav
   real(single),   allocatable, dimension (:,:)   :: effective_radius_1621_final_sav
   real(single),   allocatable, dimension (:,:)   :: effective_radius_37_final_sav
   real(single),   allocatable, dimension (:,:)   :: effective_radius_21_final_sav
   real(single),   allocatable, dimension (:,:)   :: effective_radius_22_final_sav
   real(single),   allocatable, dimension (:,:)   :: effective_radius_16_final_sav
   real(single),   allocatable, dimension (:,:)   :: optical_thickness_1621_final_sav
   real(single),   allocatable, dimension (:,:)   :: optical_thickness_37_final_sav
   real(single),   allocatable, dimension (:,:)   :: optical_thickness_22_final_sav
   real(single),   allocatable, dimension (:,:)   :: optical_thickness_16_final_sav
   real(single),   allocatable, dimension (:,:)   :: optical_thickness_final_sav
   real(single),   allocatable, dimension (:,:)   :: optical_thickness_21_final_sav
   type(processflag)  , allocatable,dimension (:,:)   :: cloudsummary_sav
   integer(integer_onebyte),      allocatable, dimension (:,:)   :: CSR_flag_array_sav
   real, dimension(:,:,:), allocatable :: atm_corr_refl_sav
   integer*2, dimension(:,:), allocatable :: tau_liquid_sav, tau_ice_sav, re21_liquid_sav, re21_ice_sav

! WDR end of save arrays
#endif
   
end module core_arrays

