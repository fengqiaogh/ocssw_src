module mod06_run_settings
  implicit none

  real(4), parameter :: set_threshold_solar_zenith     = 3.,  &
                             set_threshold_sensor_zenith    = 3.,  &
                             set_threshold_relative_azimuth = 3.

  real, parameter :: threshold_wind_speed = 1.0

  real, parameter :: dscat1 = 1.0, &
					 dscat2 = 2.0, &
					 dscat3 = 3.0


  integer, parameter :: meas_uncertain_threshold = 15


  ! integer, parameter :: set_tilesize = 100, set_number_of_bands = 15
  ! WDR increase tile/chunk size to very large
  integer, parameter :: set_tilesize = 50000, set_number_of_bands = 15
  integer, parameter :: set_albedo_bands = 6

  integer, parameter, dimension(set_number_of_bands) :: set_bands = (/1,2,5,6,7,19,20,26,29,31,8,32,34,3,4/)
  

  integer, parameter :: band_0065 = 1, &
                        band_0086 = 2, &
                        band_0124 = 3, &
                        band_0163 = 4, &
                        band_0213 = 5, &
                        band_0935 = 6, &
                        band_0370 = 7, &
                        band_0138 = 8, &
                        band_0850 = 9, &
                        band_1100 = 10, &
						band_0410 = 11, &
						band_1200 = 12, &
						band_1350 = 13, &
						band_0047 = 14, &
						band_0055 = 15
  ! WDR what was the 3.7 um band for modis will be the 2.26 um band for OCI
  integer, parameter :: band_0226 = 7

  integer, parameter :: channel_37um = set_bands(band_0370)
  integer, parameter :: channel_11um = set_bands(band_1100)
  integer, parameter :: channel_12um = set_bands(band_1200)

! cost function indices, #3 will be 2.2 um for OCI, 3.7 for modis
  integer, parameter :: re16 = 1, re21 = 2, re37= 3, re1621 = 4, re22 = 3
  integer, parameter :: set_cf_bands = 4

  integer, parameter, dimension(2) :: set_start  = (/    0 ,    0 /), &
                                      set_edge   = (/ 1354 , 2050 /), &
                                      set_stride = (/    1 ,    1 /)

! SET THIS VARIABLE TO .true. TO ENABLE Clear Sky Restoral AND TO .false. TO DISABLE IT
  logical, parameter :: DO_CSR = .true.

! SET THIS VARIABLE TO .true. TO ENABLE Collection 004 phase algorithm AND TO .false. TO DISABLE IT
  logical, parameter :: DO_C4PhaseTest = .false.

! SET THIS VARIABLE TO .true. TO ENABLE forcing of cloud phase for the entire granule to ice phase
  logical, parameter :: FORCE_ICE = .false.

! SET THIS VARIABLE TO .true. TO ENABLE forcing of cloud phase for the entire granule to water phase
  logical, parameter :: FORCE_WATER = .false.

! SET THIS VARIABLE TO .true. TO ENABLE the Cox-Munk BRDF model for ocean surface
  logical, parameter :: DO_COX_MUNK = .true.


end module mod06_run_settings
