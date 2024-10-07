 module names
 
#if !CRTM_CT & !ACAERO 
  use mod06_run_settings
#endif 
 
  implicit none
  
  integer :: MYYEAR, MYDAY, MYTIME, MYMONTH, MYDATE, MYMSG
  real ::  FRAC_TIME
  character(len = 2) :: MYMISSION

  character(len = 1000) :: &
#if CRTM_CT | ACAERO
  Alevel1b_name, &
#else
  Alevel1b_name(set_number_of_bands),   &
#endif
                           Acloudmask_name,             &
                           Amod06_name,            &
                           ACTH_name, &
                           Ageolocation_name,           &
                           Agdas_name,                  &
                           Agdas_name2,                  &
						   Aozone_name, &
                           Ancepice_name,               &
                           Anise_name, ASHIS_name, Aangle_name, ALST_name               
	character(len = 1000) :: Awater_library,              &
                           Aice_library, ADEM_name, Amod06source_name
! GEOS-5 FP-IT product names
	character(len = 1000) :: Ageos2d_name1, Ageos2d_name2, & 
							Ageos3d_name1, Ageos3d_name2, &
							Ageos_ocn_name, Ageos_lnd_name					   
						   
	character(len =1000) :: Alibnames_ice(3), Alibnames_ice_sdev(3)					   
	character(len =1000) :: Alibnames_water(3), Alibnames_water_sdev(3)					   
	character(len = 1000) :: Aphase_library,				&
                           Atransmittance_library,      &
                           Aecosystem_data_name,        &
                           Asnowicealbedo_data_name
	character(len=1000) :: A138_library,			&
						   A138_slope_library,		&
						   A138_slope_stdev_library
	character(len=500) :: ACT_lib_path, Aemissivity_name, ACSRBias_name, Atmp_dirname

	character(len=500) :: MY_TEXT_FILE
	integer, parameter :: MY_UNIT_LUN = 20

! MOD_PRCRTMCT names	
	character(len=500) :: Aforward_uncertainty_name, Agas_name, ACRTM_path, AIce_path, ALiquid_path, ACRTM_output_name

! ACAERO names
	character(len=300) :: name_model, name_out
	character(len=3) :: product, ancillary_time
														    
 end module names
 
 
