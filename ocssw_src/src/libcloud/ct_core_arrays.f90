module ct_core_arrays

   use modis_cloudstructure

   use science_parameters
   use global_model_grids

   implicit none

	type seviri_cm_type 
	
     integer*1    ::  cloudmask_determined
     integer*1    ::  confident_cloudy
    
!    processing path flags 
	 integer*1    ::  water_surface
     integer*1    ::  land_surface
     integer*1    :: desert_surface
     integer*1    :: snowice_surface
     
     integer*1    :: mlayer
	
	end type seviri_cm_type


   type(cloudphase),allocatable,dimension (:,:)   :: ir_cloudphase

   !input data arrays (measurement)
   real,   allocatable, dimension (:,:)   :: sensor_zenith_angle
   real,   allocatable, dimension (:,:)   :: scattering_angle

   real,   allocatable, dimension (:,:,:) :: band_measurements
#ifdef SEVIRI_INST
   type(seviri_cm_type),allocatable,dimension (:,:)   :: cloudmask
   integer*2,   allocatable, dimension (:,:,:)   :: surface_emissivity_land
   real, allocatable, dimension(:,:) :: cloud_emissivity
#else
   type(cloudmask_type),allocatable,dimension (:,:)   :: cloudmask
   real,   allocatable, dimension (:,:,:)   :: surface_emissivity_land
#endif

   real,   allocatable, dimension (:,:)   :: latitude
   real,   allocatable, dimension (:,:)   :: longitude
   	real, dimension(:,:), allocatable :: snow_cover

   real,   allocatable, dimension (:,:)   :: surface_temperature
#ifdef MODIS_INST
   real, allocatable, dimension(:,:,:) :: CSRBias
#else 
	real :: CSRBias(4)
#endif
   
   real,   allocatable, dimension (:,:)   :: cloud_top_temperature
   real,   allocatable, dimension (:,:)   :: cloud_top_pressure
   
#ifdef AMS_INST
   real,   allocatable, dimension (:,:)   :: abovecloud_watervapor   
   real,   allocatable, dimension (:,:)   :: abovecloud_watervapor_alt   
   real, dimension(:,:,:,:), allocatable :: transmit_correction_table

#endif	   

    integer*1, dimension(:,:), allocatable :: cloud_height_method
	integer*1, dimension(:,:), allocatable :: overshoot_top
	
   real, allocatable, dimension(:,:) :: emiss11, emiss85, emiss12, emiss13

   !Other
   character*12                                   :: platform_name
   
   type(ancillary_type), allocatable, dimension(:,:)  :: model_info


	real, dimension(:,:), allocatable :: beta_8511_store, beta_1211_store, beta_1311_store

    integer*1, dimension(:,:), allocatable :: seviri_cloudphase
   real,   allocatable, dimension (:,:)   :: cloud_top_height

   integer*2, allocatable, dimension(:,:) :: seviri_dem
   
   real, dimension(:,:,:), allocatable :: clear_sky_rad
   real, dimension(:,:,:), allocatable :: clear_sky_btemp
   
end module ct_core_arrays

