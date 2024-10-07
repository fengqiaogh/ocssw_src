module modis_sciencestructure

  type qualityanalysis

#if VIIRS_INST | AHI_INST | AMS_INST | ASTER_INST | AVIRIS_INST | EPIC_INST
	
	integer(1) :: bowtie_pixel
	integer(1) :: bad_radiance_65
	integer(1) :: bad_radiance_86
	integer(1) :: bad_radiance_12
	integer(1) :: bad_radiance_16
	integer(1) :: bad_radiance_21
	integer(1) :: bad_radiance_37
	integer(1) :: bad_radiance_11

	integer(1) :: spectral_VNSWIR_21
	integer(1) :: spectral_VNSWIR_1621
	integer(1) :: spectral_VNSWIR_16
	integer(1) :: spectral_VNSWIR_37

#endif

#if EPIC_INST
	integer(1) :: icetau_outofbounds
#endif


!   product quality and retrieval processing QA flags   
    integer(1) :: optical_thickness_GC
    integer(1) :: optical_thickness_outofbounds
    integer(1) :: effective_radius_GC
    integer(1) :: water_path_GC
    integer(1) :: rayleigh_correction
! since they have to go together as it is, let them sit together in one byte instead of two
    integer(1) :: path_and_outcome
    integer*1 :: path_and_outcome_PCL

    integer*1 :: path_and_outcome_16
    integer*1 :: path_and_outcome_16_PCL

    integer*1 :: path_and_outcome_37
    integer*1 :: path_and_outcome_37_PCL
    
    integer(1) :: path_and_outcome_1621
    integer*1 :: path_and_outcome_1621_PCL

    integer*1 :: path_and_outcome_22
    integer*1 :: path_and_outcome_22_PCL
    
    integer(1) :: band_used_for_optical_thickness
    integer(1) :: optical_thickness_1621_GC
    integer(1) :: effective_radius_1621_GC
    integer(1) :: water_path_1621_GC
    integer(1) :: multi_layer_cloud
    integer(1) :: CSR_flag
	integer(1) :: ml_test_mark
	
#if SEVIRI_INST | AHI_INST
	integer :: Tc_override
	integer :: day_flag
#endif	
	
  end type qualityanalysis
  
  type stat_type
  
	real :: retrieval_fraction 
	real :: land_fraction 
	real :: water_fraction 
	real :: snow_fraction 
	real :: cloud_fraction 
	real :: water_cloud_fraction 
  	real :: ice_cloud_fraction 
	real*8 :: mean_liquid_tau 
	real*8 :: mean_ice_tau
	real*8 :: mean_liquid_re21 
	real*8 :: mean_ice_re21
	real*8 :: ctp_liquid 
	real*8 :: ctp_ice
	real*8 :: ctp_undetermined
	real*8 :: ctt_liquid 
	real*8 :: ctt_ice
	real*8 :: ctt_undetermined
  
  end type stat_type
  
  type failed_type
  
  	integer*2 :: tau
  	integer*2 :: re
  	integer*2 :: RSS
  
  end type failed_type
  
  
  
end module modis_sciencestructure
