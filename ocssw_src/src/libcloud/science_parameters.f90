module science_parameters

  implicit none

! A solar zenith threshold of 81.36 degrees (cosine(.15)) is used. 
! Chosen to be a value within the current upper limit of the 
! forward reflectance library

  real, parameter :: solar_zenith_threshold = 81.36
  real, parameter :: sensor_zenith_threshold = 81.36
  logical :: COX_MUNK, last_COX_MUNK
  
! used in solution_re as the designation for the retrieval failure status, as distinguished from fill and valid values
  real, parameter         :: INVALID_ATTEMPTED_BUT_FAILED_RE = -1.

! used the in the calculation of the (effective) liquid water path
  real, parameter :: liquid_water_density = 1.0, ice_water_density = .93

! estimated relative uncertainties used in uncertainty computation
  real, parameter :: watervapor_error = .2    ! used to get transmittance uncertainty
  real, parameter :: albedo_error = .15
  real, parameter :: measurement_error = .05
  real, parameter :: wind_speed_error = .2
  
  real, parameter :: unc_scale_factor = 0.01
  real, parameter :: retr_scale_factor = 0.01
  real, parameter :: MAX_UNCERTAINTY = 200.0
  real, parameter :: VNIR_error = 0.02
  real, parameter :: SWIR_error = 0.03
  
  
  real, parameter :: delta_Ts = 1.0 !K
  real, parameter :: delta_Pc = 50. !mb
  
  real :: solar_constant_37
  integer, parameter :: model_levels = 101
  
  integer :: grid_xsize
  integer :: grid_ysize
  
  integer :: geos5_istart, geos5_jstart, geos5_iend, geos5_jend
  
  real :: max_solar_zenith, min_solar_zenith
  real :: max_sensor_zenith, min_sensor_zenith
  real :: max_rel_azimuth, min_rel_azimuth

  real :: lastinterp_scat_angle, lastinterp_solar_zenith,     &
                          lastinterp_sensor_zenith,    &
                          lastinterp_relative_azimuth, &
						  lastinterp_wind_speed
  real :: lastinterp_scat_angle_ss	! KGM 11-1-13
	integer*1 :: no_valid_data

	logical :: color_first_time
	real :: last_2way_angle

	real :: spec_uncertain(6), uncertain_sf(6)
	
	
	real, parameter :: d2r = 0.017453292519943295
	

end module science_parameters
