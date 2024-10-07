module libraryinterpolates
  use GeneralAuxType
  
  real(single), allocatable, dimension(:,:,:) :: int_reflectance_water,int_reflectance_ice
  real(single), allocatable, dimension(:,:,:) :: int_reflectance_water_sdev,int_reflectance_ice_sdev
  
  real(single), allocatable, dimension(:,:,:,:) :: int_reflectance_water_wspeed, int_reflectance_ice_wspeed
  real(single), allocatable, dimension(:,:,:,:) :: int_refl_water_sdev_wspeed, int_refl_ice_sdev_wspeed
 
! sep, 6 May: The retrieval nomenclature, consistent in the three routines 
! vis_absorbing_science, nir_absorbing_science, and vis_nonabsorbing_science is
! as follows (note that there is now consistency in that '0' => solar, '1'=> sensor):
!
! int_fluxup_sensor           -> fri1
! int_fluxdown_solar          -> fti0
! int_fluxdown_sensor         -> fti1
! ---------------------


  real(single), allocatable, dimension(:,:,:) :: int_fluxupwater_sensor, int_fluxdownwater_solar, int_fluxdownwater_sensor 
  real(single), allocatable, dimension(:,:,:) :: int_fluxupice_sensor, int_fluxdownice_solar, int_fluxdownice_sensor

  real(single), allocatable, dimension(:,:,:) :: int_cloud_emissivity_water, int_surface_emissivity_water
  real(single), allocatable, dimension(:,:,:) :: int_cloud_emissivity_ice, int_surface_emissivity_ice
  real(single), allocatable, dimension(:,:,:) :: int_cloud_emissivity_water_sdev, int_surface_emissivity_water_sdev
  real(single), allocatable, dimension(:,:,:) :: int_cloud_emissivity_ice_sdev, int_surface_emissivity_ice_sdev

  real(single), allocatable, dimension(:,:,:,:) :: int_cloud_emis_water_wspeed, int_surface_emis_water_wspeed
  real(single), allocatable, dimension(:,:,:,:) :: int_cloud_emis_ice_wspeed, int_surface_emis_ice_wspeed
  real(single), allocatable, dimension(:,:,:,:) :: int_cloud_emis_water_sdev_wspeed, int_surface_emis_water_sdev_wspeed
  real(single), allocatable, dimension(:,:,:,:) :: int_cloud_emis_ice_sdev_wspeed, int_surface_emis_ice_sdev_wspeed




end module libraryinterpolates
