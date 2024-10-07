module libraryarrays

  use GeneralAuxType
  !library descriptor parameters
  integer                   :: number_wavelengths, number_waterradii, number_iceradii
  integer                   :: number_solarzenith, number_relazimuth, number_sensorazimuth
  integer                   :: number_taus, number_sensorzenith, number_phase_angles_water
  integer                   :: number_fluxsolarzenith, number_fluxsensorzenith, number_phase_angles_ice
  integer					:: number_wind_speed
  
  ! library bidirectional reflectance
  ! added the wind speed dimension for cox-munk, index 1 is the lambertian library
  ! WSDR for dimension manager, no need for these MS arrays
  !real(single), allocatable, dimension(:,:,:,:,:,:,:)   :: library_reflectance_water
  !real(single), allocatable, dimension(:,:,:,:,:,:,:)   :: library_reflectance_ice

  !real(single), allocatable, dimension(:,:,:,:,:,:,:)   :: library_reflectance_water_sdev
  !real(single), allocatable, dimension(:,:,:,:,:,:,:)   :: library_reflectance_ice_sdev

  ! water vapor band cirrus retrieval library
  real(single), allocatable, dimension(:,:,:,:,:,:)   :: library_reflectance_138, library_reflectance_138_sdev
  real(single), allocatable, dimension(:,:,:,:,:,:)   :: library_slope_w03, library_slope_stdv_w03
  real(single), allocatable, dimension(:,:,:,:,:,:)   :: library_slope_w07, library_slope_stdv_w07
  real(single), allocatable, dimension(:,:,:,:,:,:)   :: library_slope_w15, library_slope_stdv_w15


  ! flux parameters

  !  extinction  		ce             extinction coefficient
  !  singlescattering		w              single scattering albedo
  !  asymmetry 			g              asymetry factor
  !  spherical_albedo		sfr            spherical albedo
  !  flux_up			fr             upwards flux
  !  flux_down			ft             downwards flux
  ! 


  real, allocatable,dimension(:,:,:,:)          :: transmit_correction_table, transmit_stddev_table

  real(single), allocatable, dimension(:)       :: library_taus
  real(single), allocatable, dimension(:)       :: library_solar_zenith, library_relative_azimuth
  real(single), allocatable, dimension(:)       :: library_sensor_zenith
  real(single), allocatable, dimension(:)       :: library_fluxsolarzenith, library_fluxsensorzenith
  
  real(single), allocatable, dimension(:)       :: rayleigh_tau, aerosol_tau, aerosol_asym, aerosol_ssa

  real(single), allocatable, dimension(:)       :: library_138_wavelengths, library_138_taus
  real(single), allocatable, dimension(:)       :: library_slope_wavelengths, library_slope_wspeed
  real(single), allocatable, dimension(:)       :: library_138_solar_zenith, library_138_relative_azimuth, &
												   library_138_sensor_zenith
  real(single), allocatable, dimension(:)       :: ice_radii_138
  real(single), allocatable, dimension(:,:)     :: asymmetry_ice_138

  real(single), allocatable, dimension(:)       :: water_radii, phase_angles_water
  real(single), allocatable, dimension(:,:)     :: truncation_factor_water, phase_fun_norm_constant_water,&
													extinction_water, asymmetry_water
  real(single), allocatable, dimension(:,:)     :: singlescattering_water
  real(single), allocatable, dimension(:,:,:)   :: spherical_albedo_water, phase_funcs_water 
  real(single), allocatable, dimension(:,:,:,:) :: flux_up_water_solar, flux_down_water_solar
  real(single), allocatable, dimension(:,:,:,:) :: flux_up_water_sensor, flux_down_water_sensor

  
  real(single), allocatable, dimension(:)       :: ice_radii, phase_angles_ice
  real(single), allocatable, dimension(:,:)     :: singlescattering_ice, extinction_ice
  real(single), allocatable, dimension(:,:)     :: asymmetry_ice
  real(single), allocatable, dimension(:,:)     :: truncation_factor_ice, phase_fun_norm_constant_ice
  real(single), allocatable, dimension(:,:,:)   :: spherical_albedo_ice, phase_funcs_ice
  real(single), allocatable, dimension(:,:,:,:) :: flux_up_ice_solar, flux_down_ice_solar  
  real(single), allocatable, dimension(:,:,:,:) :: flux_up_ice_sensor, flux_down_ice_sensor
  real(single), allocatable, dimension(:,:,:,:,:) :: cloud_emissivity_water, surface_emissivity_water !4D plus wind speed
  real(single), allocatable, dimension(:,:,:,:,:) :: cloud_emissivity_ice, surface_emissivity_ice !4D plus wind speed

  real(single), allocatable, dimension(:,:,:,:,:) :: cloud_emissivity_water_sdev, surface_emissivity_water_sdev !4D plus wind speed
  real(single), allocatable, dimension(:,:,:,:,:) :: cloud_emissivity_ice_sdev, surface_emissivity_ice_sdev !4D plus wind speed

! this is for the alternate retrieval so we can have same exact code for Cox-Munk and Lambertian
   real, dimension(:,:), allocatable :: reflibA, reflibB, rad37lib
   real, dimension(:), allocatable :: rayleigh_liq, rayleigh_ice
   

  
end module libraryarrays
