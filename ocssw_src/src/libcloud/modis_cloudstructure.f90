module modis_cloudstructure

! These desginations can only really be understood in conjunction with 
! the following document:
! http://modis-atmos.gsfc.nasa.gov/MOD35_L2/format.html
! Commented field names are from this document

! mag Sept 2003

  type cloudmask_type

!    top level flags
     logical    ::  cloudmask_determined
     logical    ::  confident_cloudy
     logical    ::  probablyclear_66
     logical    ::  probablyclear_95
     logical    ::  probablyclear_99
    
!    processing path flags 
     integer*1    ::  night
     integer*1    ::  sunglint
     logical    ::  snowice_surface
	 
     logical    ::  water_surface
     logical    ::  coastal_surface
     logical    ::  desert_surface 
     logical    ::  land_surface

	 integer*1 :: ocean_no_snow
	 integer*1 :: ocean_snow
	 integer*1 :: land_no_snow
	 integer*1 :: land_snow

#ifndef SEVIRI_INST

#if EMAS_INST || AMS_INST || MAS_INST

	logical :: test_high_67
 	logical :: test_high_138
  	logical :: test_irdiff
	logical :: test_39_11
	logical :: test_visiblereflectance

 	logical :: applied_highcloud67
 	logical :: applied_highcloud138
  	logical :: applied_irtempdifference
  	logical :: applied_39_11
	logical :: applied_visiblereflectance

#else
!    1km cloud test flags
     integer*1    ::  test_high_138                 ! High Cloud Flag (1.38 micron Test) 
	 integer*1    :: test_visiblereflectance
     integer*1    ::  test_visnirratio	! KGM 3-4-13, GW 3.28.13

!    250m cloud test flags
     integer*1    ::  visible_cloudtest_250m	
 
!    1km cloud test application flags
 
     integer*1    ::  applied_highcloud138
     integer*1    ::  applied_visiblereflectance
     integer*1    ::  applied_visnirratio		! KGM 3-4-13  GW 3.28.13
#endif

#else

	 integer*1 :: partly_cloudy

#endif


  end type cloudmask_type

  type processflag   ! WDR-UIV
     logical    ::  cloudmask_determined = .false.
     logical    ::  cloudobserved = .false.
     logical    ::  land_surface = .false.
     logical    ::  ocean_surface = .false.
     logical    ::  coastal_surface = .false.
     logical    ::  desert_surface = .false.
     logical    ::  snowice_surface = .false.
     logical    ::  watercloud = .false.
     logical    ::  icecloud = .false.
     logical    ::  unknowncloud = .false.
  end type processflag

  type cloudphase
     integer*1    ::  watercloud
     integer*1    ::  icecloud
     integer*1    ::  unknowncloud
#ifdef CT_CODE
	 integer*1  :: mlayer
#endif

  end type cloudphase



end module modis_cloudstructure
