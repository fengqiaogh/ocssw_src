module modis_frontend_module

  implicit none
 
  include 'netcdf.inc'
 
  private

  integer   :: index_solarzenith, index_sensorzenith,&
              index_relativeazimuth,index_solarzenith_flux

  real, dimension(:), allocatable  :: solarzenith_all, sensorzenith_all, &
                                     relativeazimuth_all, solarzenith_flux_all


  public :: initialize_run, check_datasources, readlibraries_base, &
          readlibraries_extra, &
          allocate_arrays, init_qualitydata, deallocate_cleanup, allocate_model

  contains
 
  subroutine initialize_run( start, edge, stride,        &
                           tilesize,                   &
                           threshold_solar_zenith,     &
                           threshold_sensor_zenith,    &
                           threshold_relative_azimuth, &
                           status )
  !  WDR 
  !  tilesize  O  the size of the data block to process - in case of MODIS,
  !       AVHRR, it is edge (1) the # pixels?
  ! 
  use core_arrays
  use mod06_run_settings
  use nonscience_parameters
  use names
  use specific_ancillary
  use specific_other
  use ch_xfr, only : c2_sensor_id, OCI_ID, OCIS_ID

  implicit none


  integer,  intent (out)         :: status, tilesize 

  real,     intent(out)  :: threshold_solar_zenith,      &
                            threshold_sensor_zenith,     &
                            threshold_relative_azimuth

  integer,  dimension (2),  intent (out) :: start,edge,stride 

  integer   :: number_of_bands, index, file_version, localstatus 

  !  function definitions
  ! WDR integer get_platform_name
  ! WDR external get_platform_name

  status = 0

  ! WDR only set the name to Aqua
  platform_name = "Aqua"
  ! WDR 13 jun 22, use c2_sensor_id to set to OCI if so
  ! temp to test incrimentaly if( c2_sensor_id == OCI_ID ) platform_name = "OCI"
  !localstatus = get_platform_name(Alevel1b_name(1), platform_name)

  call get_channels
 
  ! WDR remove these entirely
  !if (platform_name == "RSP" .or. platform_name == "SSFR") then 
  !	call get_data_dims(Alevel1b_name(1), start, stride, edge)
  !else
  !	call get_data_dims(Amod06_name, start, stride, edge)
  !endif

  threshold_solar_zenith = set_threshold_solar_zenith
  threshold_sensor_zenith = set_threshold_sensor_zenith
  threshold_relative_azimuth = set_threshold_relative_azimuth

  !	print*, size(Alevel1b_name)
  
#ifdef VIIRS_INST
  if (MODIS_MODE) then 
    tilesize = set_tilesize_modis
	 number_of_bands = size(set_bands_modis)
	 allocate(bands(number_of_bands))
    bands(:) = set_bands_modis(:)
	 channel_37um = set_bands_modis(band_0370)
	 channel_11um = set_bands_modis(band_1100)
	 channel_12um = set_bands_modis(band_1200)
  else 
    tilesize = set_tilesize
	 number_of_bands = size(set_bands)
	 allocate(bands(number_of_bands))
	 bands(:) = set_bands(:)
	 channel_37um = set_bands(band_0370)
	 channel_11um = set_bands(band_1100)
	 channel_12um = set_bands(band_1200)
  endif
#else

	if (platform_name == "RSP" .or. platform_name == "SSFR") then 
	  tilesize = edge(1)
	else
     tilesize = set_tilesize
	endif
	number_of_bands = size(set_bands)
   ! WDR re-entrant
	if( .not. allocated(bands) ) allocate(bands(number_of_bands))
	bands(:) = set_bands(:)

#endif	

 end subroutine initialize_run

 subroutine check_datasources(status)

    use nonscience_parameters
	use names
   use ch_xfr, only : cm_from_l2
    implicit none
	
    integer, intent(out)       :: status
	integer :: i, allstatus

    status = 0
    allstatus = 0

! WDR remove this check
!    allstatus = checkfile(Alevel1b_name(1))
!    if (allstatus /= success) then
!      status = failure
!      call local_message_handler ('Error noted: level 1B file. Check input.' ,status, 'check_datasources')
!	  return
!    endif
!
! WDR only use work file if instructed
! WDR out forever, I think
!    if ( cm_from_l2 .EQ. 0 ) THEN
!      allstatus = checkfile(Acloudmask_name)
!      if (allstatus /= success) then
!        status = failure
!        call local_message_handler ('Error noted: cloud mask file. Check input.',  status,'check_datasources')
!    	  return
!      endif
!   END IF
! WDR remove this check
!    allstatus = checkfile(Ageolocation_name)
!    if (allstatus /= success) then
!      status = failure
!      call local_message_handler ('Error noted: geolocation file. Check input.' ,status,'check_datasources')
!	  return
!    endif

#ifdef GEOS5
#ifndef SSFR_INST
	allstatus = checkfile(Ageos3d_name1)
    if (allstatus /= success) then
      status = failure
      call local_message_handler ('Error noted: first GEOS5 3D file. Check input.' ,status, 'check_datasources')
	  return
    endif
	allstatus = checkfile(Ageos3d_name2)
    if (allstatus /= success) then
      status = failure
      call local_message_handler ('Error noted: second GEOS5 3D file. Check input.' ,status, 'check_datasources')
	  return
    endif

	allstatus = checkfile(Ageos2d_name1)
    if (allstatus /= success) then
      status = failure
      call local_message_handler ('Error noted: first GEOS5 2D file. Check input.' ,status, 'check_datasources')
	  return
    endif
	allstatus = checkfile(Ageos2d_name2)
    if (allstatus /= success) then
      status = failure
      call local_message_handler ('Error noted: second GEOS5 2D file. Check input.' ,status, 'check_datasources')
	  return
    endif
#endif	
#else    
! WDR remove this check
!    allstatus = checkfile(Agdas_name)    
!    if (allstatus /= success) then
!      status = failure
!      call local_message_handler ('Error noted: first GDAS file. Check input.' , status, 'check_datasources')
!	  return
!    endif
! WDR remove this check
!    allstatus = checkfile(Agdas_name2)    
!    if (allstatus /= success) then
!      status = failure
!      call local_message_handler ('Error noted: second GDAS file. Check input.' ,status,  'check_datasources')
!	  return
!    endif
#endif

! WDR remove this check
!    allstatus = checkfile(Ancepice_name)    
!    if (allstatus /= success) then
!      status = failure
!      call local_message_handler ('Error noted: NCEP Ice file. Check input.' ,status, 'check_datasources')
!	  return
!    endif

#ifdef USE_TOAST	
    allstatus = checkfile(Aozone_name)    
    if (allstatus /= success) then
      status = failure
      call local_message_handler ('Error noted: TOAST ozone file. Check input.' ,status, 'check_datasources')
	  return
    endif
#endif
	
! WDR remove this check
!    allstatus = checkfile(Anise_name)    
!    if (allstatus /= success) then
!      status = failure
!      call local_message_handler ('Error noted: NISE  file. Check input.' , status, 'check_datasources')
!	  return
!    endif

! WDR remove this check
!    allstatus = checkfile(Amod06_name)    
!    if (allstatus /= success) then
!      status = failure
!      call local_message_handler ('Error noted: MOD06_L2 file. Check input.' , status, 'check_datasources')
!	  return
!    endif
    allstatus = checkfile(Awater_library)    
    if (allstatus /= success) then
      status = failure
      call local_message_handler ('Error noted: Lambertian water library file. Check input.' ,  status, 'check_datasources')
	  return
    endif
    
    allstatus = checkfile(Aice_library)    
    if (allstatus /= success) then
      status = failure
      call local_message_handler ('Error noted: Lambertian ice library file. Check input.' ,status, 'check_datasources')
	  return
    endif
    allstatus = checkfile(Atransmittance_library)    
    if (allstatus /= success) then
      status = failure
      call local_message_handler ('Error noted: transmittance file. Check input.' , status, 'check_datasources')
	  return
    endif
    allstatus = checkfile(Aecosystem_data_name)    
    if (allstatus /= success) then
      status = failure
      call local_message_handler ('Error noted: IGBP Ecosystem file. Check input.' ,status, 'check_datasources')
	  return
    endif
    allstatus = checkfile(Asnowicealbedo_data_name)    
    if (allstatus /= success) then
      status = failure
      call local_message_handler ('Error noted: Static SnowIce Albedo file. Check input.' ,status, 'check_datasources')
	  return
    endif
    allstatus = checkfile(Aphase_library)    
    if (allstatus /= success) then
      status = failure
      call local_message_handler ('Error noted: phase library file. Check input.' ,status, 'check_datasources')
	  return
    endif
	
	do i=1, 3
		allstatus = checkfile(Alibnames_water(i))    
		if (allstatus /= success) then
			status = failure
			call local_message_handler ('Error noted: Cox-Munk water library file. Check input.' , status, 'check_datasources')
	  return
		endif
		allstatus = checkfile(Alibnames_water_sdev(i))    
		if (allstatus /= success) then
			status = failure
			call local_message_handler ('Error noted: Cox-Munk water sdev library file. Check input.' ,status, 'check_datasources')
	  return
		endif
		allstatus = checkfile(Alibnames_ice(i))    
		if (allstatus /= success) then
			status = failure
			call local_message_handler ('Error noted: Cox-Munk ice library file. Check input.' ,status, 'check_datasources')
	  return
		endif
		allstatus = checkfile(Alibnames_ice_sdev(i))    
		if (allstatus /= success) then
			status = failure
			call local_message_handler ('Error noted: Cox-Munk ice sdev library file. Check input.' ,status, 'check_datasources')
	  return
		endif

	end do

    allstatus = checkfile(Aemissivity_name)    
    if (allstatus /= success) then
      status = failure
      call local_message_handler ('Error noted: emissivity library file. Check input.' , status,  'check_datasources')
	  return
    endif
	
	
!	if (Alevel1b_name(2) /= "none" ) then 
!		allstatus = checkfile(Alevel1b_name(2))
!		if (allstatus /= success) then
!			status = failure
!			call local_message_handler ('Error noted: additional level 1B file. Check input.' , status, 'check_datasources')
!	  return
!		endif
!	endif
	   
 end subroutine check_datasources
 
 integer function checkfile(name)

  implicit none
  
  character(*), intent(in)  :: name
  logical                   :: exist
  character(len=11)         :: readability

  checkfile = 0

  inquire( file = name,      &
           exist= exist,   &
           read = readability)
           
  if (.not. exist) then
         checkfile = 1
  endif
  if (readability == 'NO') then
         checkfile = 2
  endif
 end function checkfile



! this subroutine returns the needed array bounds for a resized array
 subroutine find_bounds(array, min_val, max_val, min_bnd, max_bnd) 
 
	implicit none
	
	real, dimension(:), intent(in) :: array
	real, intent(in) :: min_val, max_val
	integer, intent(inout) :: min_bnd, max_bnd
	
	integer :: i, N
	
	N = size(array)
	
	do i=1, N
		if (array(i) >= min_val) then 
			min_bnd = i-1
			exit
		endif
	end do
 
	if (min_bnd < 1) min_bnd = 1
	if (min_bnd > N) min_bnd = N
 
	do i=1, N
		if (array(i) >= max_val) then 
			max_bnd = i
			exit
		endif
	end do

	if (max_bnd >= N) max_bnd = N
	if (max_bnd <= 1) max_bnd = 1

 end subroutine find_bounds
 
!-----------------------------------------------------------------------
!f90 readlibraries
!
!Description:
!
! Read bidirectional reflectance and flux library values for all
! angles, as defined in a library description file.  Library arrays
! are allocated as required to read all library data.
!
!input parameters:
!
!output parameters:
!
!revision history:
!
! v2.1  July 2002 arrays are now allocated to only include the angles
!               needed.
!
! v2.0  June 2002 Changed the read routines to work with the HDF
!         version of the libraries. -- Gala. wind@climate.gsfc.nasa.gov
!
! v1.0  November 2001 Initial work mag gray@climate.gsfc.nasa.gov
!
!team-unique header:
!
! Cloud Retrieval Group, NASA GSFC, Greenbelt, Maryland, USA 
!
!references and credits:
!
! Mark Gray
! gray@climate.gsfc.nasa.gov
! EmergentIT 
! Code 913, NASA GSFC 
! Greenbelt, MD 20771
!
! Gala Wind
! wind@climate.gsfc.nasa.gov
! L-3 Comm Analytics
! all same, all same
!
!design note:
!
!end
!----------------------------------------------------------------------

  subroutine readlibraries_base(debug,status)
  use GeneralAuxType
  use libraryarrays
  use libraryinterpolates
  use science_parameters
  use nonscience_parameters
  use interpolate_libraries 
  use names, only: Atransmittance_library, Awater_library, Aice_library, &
    Alibnames_ice, Aphase_library
  use ch_xfr, only: c2_sensor_id, OCI_ID, OCIS_ID
  implicit none  

  logical, intent(in)        :: debug
  integer,intent (out)       :: status

  integer                    :: file_id, var_id, var_index,i, j
  integer                    :: start(6), edge(6)
  character(MAX_SDS_NAME_LEN)   :: dummy_name
  integer                    :: dummy_type, dummy_num_attr, dummy_rank
  integer  :: dim_sizes(6),dim_sizes4d(4)
  integer  :: idim, len, dimids(6)

  ! WDR I havd replaced the more extensive sf... HDF 4 I/O with calls to
  !  fortran netcdf I/O (nf_...).  The advantage is that they can also 
  !  read hdf 4 table files that were used originally
  status = 0

  ! WDR try to read this only if the allocated array is not allocated
  if( allocated( transmit_correction_table ) ) then
    status = 0
    return
  endif

!  read transmittance correction library

   !start = 0
   start = 1   ! for netcdf

   status = nf_open( Atransmittance_library, NF_NOWRITE, file_id )
   call cld_fchk( status, __FILE__, __LINE__ )

   status = nf_inq_varid( file_id, "Transmittance", var_id )
   call cld_fchk( status, __FILE__, __LINE__ )

   status = nf_inq_vardimid( file_id, var_id, dimids )
   call cld_fchk( status, __FILE__, __LINE__ )
   do idim = 1, 4
     status = nf_inq_dimlen(file_id, dimids(idim), len )
     call cld_fchk( status, __FILE__, __LINE__ )
     dim_sizes(idim) = len
     end do
  ! WDR - if the OCI transmittance is read, it should have 8 bands
  ! otherwise (Aqua) it has 7
  if( ( c2_sensor_id == OCI_ID ) .or. ( c2_sensor_id == OCIS_ID ) ) then
    if( dim_sizes(4) /= 8 ) then
      print*, "Transmittance table has wrong # bands for OCI, ", &
      __FILE__, __LINE__
      stop 10
    endif
  else
    if( dim_sizes(4) /= 7 ) then
      print*, "Transmittance table has wrong # bands for Aqua, ", &
      __FILE__, __LINE__
      stop 10
    endif
  endif
  
  edge(1:4) = dim_sizes(1:4)

  allocate(transmit_correction_table(edge(1), edge(2), edge(3), edge(4)))
  allocate(transmit_stddev_table(edge(1), edge(2), edge(3), edge(4)))

  status = nf_get_vara_real( file_id, var_id, start(1:4), edge(1:4), &
      transmit_correction_table )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_varid( file_id, "Standard_Deviation", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:4), edge(1:4), &
    transmit_stddev_table )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_close( file_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  ! we need to find angle range for the granule, so we can only read a 
  ! limited amount from the library
  ! just enough to get the granule processed. 
  
  ! from the find_granule_angle_range routine we know just how wide the 
  ! angle space we need. 
  ! it's fixed for sensor zenith and relative azimuth, but not for solar 
  ! angle. However for MAS<TER>
  ! the relative azimuth space isn't fixed either and the sensor zenith 
  ! limits are different. 


  number_wind_speed = 3

  ! start with the ice library first
  status = nf_open( Aice_library, NF_NOWRITE, file_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  ! get the number of wavelengths, we never use them, so we don't need 
  ! to actually read them. 
  status = nf_inq_varid( file_id, "Wavelengths", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_vardimid( file_id, var_id, dimids )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_dimlen(file_id, dimids(1), len )
  call cld_fchk( status, __FILE__, __LINE__ )

  number_wavelengths = len

  start = 1

  ! get the radii
  status = nf_inq_varid( file_id, "ParticleRadius", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_vardimid( file_id, var_id, dimids )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_dimlen(file_id, dimids(1), len )
  call cld_fchk( status, __FILE__, __LINE__ )

  number_iceradii = len
  allocate (ice_radii(number_iceradii))
  edge(1) = number_iceradii

  status = nf_get_vara_real( file_id, var_id, start(1:1), edge(1:1), &
    ice_radii)
  call cld_fchk( status, __FILE__, __LINE__ )

  ! get the taus
  status = nf_inq_varid( file_id, "OpticalThickness", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_vardimid( file_id, var_id, dimids )
  call cld_fchk( status, __FILE__, __LINE__ )
  status = nf_inq_dimlen(file_id, dimids(1), len )
  call cld_fchk( status, __FILE__, __LINE__ )

  number_taus = len
  !number_taus = dim_sizes(1)
  allocate (library_taus(number_taus))
  edge(1) = number_taus
  status = nf_get_vara_real( file_id, var_id, start(1:1), edge(1:1), &
    library_taus )
  call cld_fchk( status, __FILE__, __LINE__ )

  ! get the full 1D angle arrays, they are small, we'll hang on to them 
  ! and do resizing after we read them
  ! in the read_libraries_extra() subroutine

  ! now get the relative azimuth
  status = nf_inq_varid( file_id, "ReflectanceRelativeAzimuth", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_vardimid( file_id, var_id, dimids )
  call cld_fchk( status, __FILE__, __LINE__ )
  status = nf_inq_dimlen(file_id, dimids(1), len )
  call cld_fchk( status, __FILE__, __LINE__ )

  index_relativeazimuth = len

  allocate (relativeazimuth_all(index_relativeazimuth))
  edge(1) = index_relativeazimuth
  status = nf_get_vara_real( file_id, var_id, start(1:1), edge(1:1), &
    relativeazimuth_all )
  call cld_fchk( status, __FILE__, __LINE__ )

  ! now get the solar zenith 
  status = nf_inq_varid( file_id, "ReflectanceSolarZenith", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_vardimid( file_id, var_id, dimids )
  call cld_fchk( status, __FILE__, __LINE__ )
  status = nf_inq_dimlen(file_id, dimids(1), len )
  call cld_fchk( status, __FILE__, __LINE__ )
  index_solarzenith = len

  allocate (solarzenith_all(index_solarzenith))
  edge(1) = index_solarzenith

  status = nf_get_vara_real( file_id, var_id, start(1:1), edge(1:1), &
    solarzenith_all )
  call cld_fchk( status, __FILE__, __LINE__ )
	
  ! now get the sensor zenith 
  status = nf_inq_varid( file_id, "ReflectanceSensorZenith", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_vardimid( file_id, var_id, dimids )
  call cld_fchk( status, __FILE__, __LINE__ )
  status = nf_inq_dimlen(file_id, dimids(1), len )
  call cld_fchk( status, __FILE__, __LINE__ )

  index_sensorzenith = len
  allocate (sensorzenith_all(index_sensorzenith))
  edge(1) = index_sensorzenith
  status = nf_get_vara_real( file_id, var_id, start(1:1), edge(1:1), &
    sensorzenith_all )
  call cld_fchk( status, __FILE__, __LINE__ )

  ! now get the flux angle 
  status = nf_inq_varid( file_id, "FluxSolarZenith", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_vardimid( file_id, var_id, dimids );
  call cld_fchk( status, __FILE__, __LINE__ )
  status = nf_inq_dimlen(file_id, dimids(1), len )
  call cld_fchk( status, __FILE__, __LINE__ )
  index_solarzenith_flux = len

  allocate (solarzenith_flux_all(index_solarzenith))
  edge(1) = index_solarzenith_flux
  status = nf_get_vara_real( file_id, var_id,start(1:1), edge(1:1), &
    solarzenith_flux_all )
  call cld_fchk( status, __FILE__, __LINE__ )

  ! Read the 2D arrays 

  start = 1
  edge(1) = number_wavelengths
  edge(2) = number_iceradii

  allocate(extinction_ice(number_wavelengths,number_iceradii))
  allocate(singlescattering_ice(number_wavelengths,number_iceradii))
  allocate(asymmetry_ice(number_wavelengths,number_iceradii))
  allocate(truncation_factor_ice(number_wavelengths,number_iceradii))
  allocate(phase_fun_norm_constant_ice(number_wavelengths,number_iceradii))

  status = nf_inq_varid( file_id, "ExtinctionCoefficient", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )
  
  status = nf_get_vara_real( file_id, var_id, start(1:2), edge(1:2), &
    extinction_ice )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_varid( file_id, "SingleScatterAlbedo", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:2), edge(1:2), &
    singlescattering_ice )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_varid( file_id, "Phase: AsymmetryParameter", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:2), edge(1:2), &
    asymmetry_ice )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_varid( file_id, "TruncationFactor", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:2), edge(1:2), &
    truncation_factor_ice )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_varid( file_id, "PhaseFuncNormConstant", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

   status = nf_get_vara_real( file_id, var_id, start(1:2), edge(1:2), &
     phase_fun_norm_constant_ice )
   call cld_fchk( status, __FILE__, __LINE__ )

  ! Read the 3D arrays
	
  start = 1
  edge(1) = number_taus
  edge(2) = number_wavelengths
  edge(3) = number_iceradii
   
  allocate(spherical_albedo_ice(number_taus,number_wavelengths,number_iceradii))
  status = nf_inq_varid( file_id, "SphericalAlbedo", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:3), edge(1:3), &
    spherical_albedo_ice )
  call cld_fchk( status, __FILE__, __LINE__ )

  ! Read the 4D arrays
  ! Now we finally read the 6D reflectance array

  status = nf_close( file_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  !  Now that we're done with ice, get started with water...
  start = 1

  !  Read all the water stuff there is to read.
  status = nf_open( Awater_library, NF_NOWRITE, file_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  ! get the radii
  status = nf_inq_varid( file_id, "ParticleRadius", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_vardimid( file_id, var_id, dimids )
  call cld_fchk( status, __FILE__, __LINE__ )
  status = nf_inq_dimlen(file_id, dimids(1), len )
  call cld_fchk( status, __FILE__, __LINE__ )

  number_waterradii = len

  allocate (water_radii(number_waterradii))
  edge(1) = number_waterradii
  status = nf_get_vara_real( file_id, var_id, start(1:1), edge(1:1), &
    water_radii )
  call cld_fchk( status, __FILE__, __LINE__ )

  ! Read the 2D arrays 

  start = 1
  edge(1) = number_wavelengths
  edge(2) = number_waterradii

  allocate(extinction_water(number_wavelengths,number_waterradii))
  allocate(singlescattering_water(number_wavelengths,number_waterradii))
  allocate(asymmetry_water(number_wavelengths, number_waterradii))
  allocate(truncation_factor_water(number_wavelengths,number_waterradii))
  allocate(phase_fun_norm_constant_water(number_wavelengths,number_waterradii))

  status = nf_inq_varid( file_id, "ExtinctionCoefficient", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:2), edge(1:2), &
    extinction_water )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_varid( file_id, "SingleScatterAlbedo", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:2), edge(1:2), &
    singlescattering_water )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_varid( file_id, "Phase: AsymmetryParameter", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:2), edge(1:2), &
    asymmetry_water )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_varid( file_id, "TruncationFactor", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:2), edge(1:2), &
    truncation_factor_water )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_varid( file_id, "PhaseFuncNormConstant", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:2), edge(1:2), &
    phase_fun_norm_constant_water )
  call cld_fchk( status, __FILE__, __LINE__ )

  ! Read the 3D arrays

	
  start = 1
  edge(1) = number_taus
  edge(2) = number_wavelengths
  edge(3) = number_waterradii
  
  allocate(spherical_albedo_water(number_taus,number_wavelengths, &
    number_waterradii))
  status = nf_inq_varid( file_id, "SphericalAlbedo", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:3), edge(1:3), &
    spherical_albedo_water )
  call cld_fchk( status, __FILE__, __LINE__ )
 
  status = nf_close( file_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  allocate(int_reflectance_water(number_taus+1, number_wavelengths,number_waterradii))
  allocate(int_reflectance_ice(number_taus+1, number_wavelengths,number_iceradii))

  allocate(int_reflectance_water_sdev(number_taus+1, number_wavelengths,number_waterradii))
  allocate(int_reflectance_ice_sdev(number_taus+1, number_wavelengths,number_iceradii))

  allocate(int_refl_water_sdev_wspeed(number_taus+1, number_wavelengths,number_waterradii, number_wind_speed))
  allocate(int_refl_ice_sdev_wspeed(number_taus+1, number_wavelengths,number_iceradii, number_wind_speed))

  allocate(int_reflectance_water_wspeed(number_taus+1, number_wavelengths,number_waterradii, number_wind_speed))
  allocate(int_reflectance_ice_wspeed(number_taus+1, number_wavelengths,number_iceradii, number_wind_speed))


  allocate(int_fluxupwater_sensor(number_taus, number_wavelengths,number_waterradii))
  allocate(int_fluxdownwater_solar(number_taus, number_wavelengths,number_waterradii))
  allocate(int_fluxdownwater_sensor(number_taus, number_wavelengths,number_waterradii))
  allocate(int_fluxupice_sensor(number_taus, number_wavelengths,number_iceradii))
  allocate(int_fluxdownice_solar(number_taus, number_wavelengths,number_iceradii))
  allocate(int_fluxdownice_sensor(number_taus, number_wavelengths,number_iceradii))

  ! WDR re-entrant
  if( .not. allocated(int_cloud_emissivity_water) ) &
    allocate(int_cloud_emissivity_water(number_taus+1, 2, number_waterradii))
  allocate(int_cloud_emissivity_ice(number_taus+1, 2, number_iceradii))
  if( .not. allocated(int_surface_emissivity_water) ) &
    allocate(int_surface_emissivity_water(number_taus+1, 2, number_waterradii))
  allocate(int_surface_emissivity_ice(number_taus+1, 2, number_iceradii))

  if( .not. allocated(int_cloud_emissivity_water_sdev) ) &
    allocate(int_cloud_emissivity_water_sdev(number_taus+1, 2, number_waterradii))
  allocate(int_cloud_emissivity_ice_sdev(number_taus+1, 2, number_iceradii))
  if( .not. allocated(int_surface_emissivity_water_sdev) ) &
    allocate(int_surface_emissivity_water_sdev(number_taus+1, 2, number_waterradii))
  allocate(int_surface_emissivity_ice_sdev(number_taus+1, 2, number_iceradii))

  allocate(int_cloud_emis_water_wspeed(number_taus+1, 2, number_waterradii, number_wind_speed))
  allocate(int_cloud_emis_ice_wspeed(number_taus+1, 2, number_iceradii, number_wind_speed))
  allocate(int_surface_emis_water_wspeed(number_taus+1, 2, number_waterradii, number_wind_speed))
  allocate(int_surface_emis_ice_wspeed(number_taus+1, 2, number_iceradii, number_wind_speed))

  allocate(int_cloud_emis_water_sdev_wspeed(number_taus+1, 2, number_waterradii, number_wind_speed))
  allocate(int_cloud_emis_ice_sdev_wspeed(number_taus+1, 2, number_iceradii, number_wind_speed))
  allocate(int_surface_emis_water_sdev_wspeed(number_taus+1, 2, number_waterradii, number_wind_speed))
  allocate(int_surface_emis_ice_sdev_wspeed(number_taus+1, 2, number_iceradii, number_wind_speed))


  allocate(rayleigh_liq(number_waterradii))
  allocate(rayleigh_ice(number_iceradii))
	
  ! Now for the phase function information

  status = nf_open( Aphase_library, NF_NOWRITE, file_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  start = 1
   
  ! get the ice phase func info
  status = nf_inq_varid( file_id, "ScatAnglesIce", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_vardimid( file_id, var_id, dimids )
  call cld_fchk( status, __FILE__, __LINE__ )
  status = nf_inq_dimlen(file_id, dimids(1), len )
  call cld_fchk( status, __FILE__, __LINE__ )

  number_phase_angles_ice = len
  allocate (phase_angles_ice(number_phase_angles_ice))
  edge(1) = number_phase_angles_ice

  status = nf_get_vara_real( file_id, var_id, start(1:1), edge(1:1), &
    phase_angles_ice )
  call cld_fchk( status, __FILE__, __LINE__ )
   
  status = nf_inq_varid( file_id, "ScatAnglesWater", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_vardimid( file_id, var_id, dimids )
  call cld_fchk( status, __FILE__, __LINE__ )
  status = nf_inq_dimlen(file_id, dimids(1), len )
  call cld_fchk( status, __FILE__, __LINE__ )

  number_phase_angles_water = len
  allocate (phase_angles_water(number_phase_angles_water))
  edge(1) = number_phase_angles_water
  status = nf_get_vara_real( file_id, var_id, start(1:1), edge(1:1), &
    phase_angles_water )
  call cld_fchk( status, __FILE__, __LINE__ )
		
  allocate(phase_funcs_water(number_phase_angles_water, number_wavelengths, number_waterradii))
  allocate(phase_funcs_ice(number_phase_angles_ice, number_wavelengths, number_iceradii))

  start = 1
  edge(1) = number_phase_angles_water
  edge(2) = number_wavelengths
  edge(3) = number_waterradii
	
  status = nf_inq_varid( file_id, "WaterPhaseFuncVals", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:3), edge(1:3), &
    phase_funcs_water )
  call cld_fchk( status, __FILE__, __LINE__ )
	
  edge(1) = number_phase_angles_ice
  edge(2) = number_wavelengths
  edge(3) = number_iceradii
	
  status = nf_inq_varid( file_id, "IcePhaseFuncVals", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:3), edge(1:3), &
    phase_funcs_ice )
  call cld_fchk( status, __FILE__, __LINE__ )
	
  status = nf_close( file_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  ! now for the aerosol and rayleigh property information

  allocate(rayleigh_tau(number_wavelengths))
  allocate(aerosol_tau(number_wavelengths))
  allocate(aerosol_asym(number_wavelengths))
  allocate(aerosol_ssa(number_wavelengths))

  status = nf_open( Alibnames_ice(1), NF_NOWRITE, file_id )
  call cld_fchk( status, __FILE__, __LINE__ )
	
  start = 1
  edge(1) = number_wavelengths
	
  status = nf_inq_varid( file_id, "RayleighOpticalThickness", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:1), edge(1:1), &
    rayleigh_tau )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_varid( file_id, "AerosolOpticalThickness", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:1), edge(1:1), &
    aerosol_tau )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_inq_varid( file_id, "AerosolAsymParameter", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )
  
  status = nf_get_vara_real( file_id, var_id, start(1:1), edge(1:1), &
    aerosol_asym )

  status = nf_inq_varid( file_id, "AerosolSSAlbedo", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:1), edge(1:1), &
    aerosol_ssa )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_close( file_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  call lib_init

  !	print*, water_radii
  !	print*, ice_radii

  end subroutine readlibraries_base
 
 
  subroutine readlibraries_extra(debug,status)
   
  use GeneralAuxType
  use libraryarrays
  use libraryinterpolates
  use science_parameters
  use nonscience_parameters
  use mod06_run_settings
  use names, only: Awater_library, Aice_library, Alibnames_ice, Alibnames_water, Alibnames_ice_sdev, Alibnames_water_sdev
  use ch_xfr, only : c2_sensor_id, OCI_ID, OCIS_ID

  !  WDR  try the setup for the call to mng_msl_init
  use iso_c_binding
  implicit none  
  interface
    function mng_ms_init(fwl, fil, fww1, fww1sd, fiw1, fiw1sd, &
       fww2, fww2sd, fiw2, fiw2sd, fww3, fww3sd, fiw3, fiw3sd )
    import
      character(kind=c_char), intent(in) :: fwl(*)
      character(kind=c_char), intent(in) :: fil(*)
      character(kind=c_char), intent(in) :: fww1(*)
      character(kind=c_char), intent(in) :: fww1sd(*)
      character(kind=c_char), intent(in) :: fiw1(*)
      character(kind=c_char), intent(in) :: fiw1sd(*)
      character(kind=c_char), intent(in) :: fww2(*)
      character(kind=c_char), intent(in) :: fww2sd(*)
      character(kind=c_char), intent(in) :: fiw2(*)
      character(kind=c_char), intent(in) :: fiw2sd(*)
      character(kind=c_char), intent(in) :: fww3(*)
      character(kind=c_char), intent(in) :: fww3sd(*)
      character(kind=c_char), intent(in) :: fiw3(*)
      character(kind=c_char), intent(in) :: fiw3sd(*)
    end function mng_ms_init
  end interface

  logical, intent(in)    :: debug
  integer,intent (out)   :: status

  integer                :: file_id, var_id, var_index,i, j
  integer					:: start(6), edge(6)
  character(MAX_SDS_NAME_LEN)  :: dummy_name
  integer                :: dummy_type, dummy_num_attr, dummy_rank
  integer                :: dim_sizes(6),dim_sizes4d(4)

  real, dimension(:,:,:,:), allocatable :: temp_emis

  integer :: max_bnd, min_bnd
  integer :: sensor_min, sensor_max, solar_min, solar_max, sza_min, &
    vza_min, ra_min
  real, dimension(:,:,:,:), allocatable :: temp_flux, temp_flux2

  integer :: start_time, end_time, crate, cmax

  real :: deg_to_rad, cos_vza_min, cos_vza_max, cos_sza_min, cos_sza_max
  integer :: i1, i2, i3, i4, i5, i6, nws

  ! WDR controls for limited reading of the ice... tables
  integer :: read_tables, ct_tbl_rd = 0
  real :: last_min_senz = -99999.
  real :: last_min_solz = -99999.
  real :: last_min_relaz = -99999.
  real :: last_max_senz = -99999.
  real :: last_max_solz = -99999.
  real :: last_max_relaz = -99999.
  real :: rng_margin_senz = 3.
  real :: rng_margin_solz = 5.  ! rng_margin_relaz not needed - take full 0-180
  ! WDR guard to set up the large array manager only once
  integer :: guard_mng_ms = 0;

  status = 0

  deg_to_rad = d2r
	
  ! WDR init the large refl table reading
  if( guard_mng_ms == 0 ) then
    status = mng_ms_init( trim(Awater_library), trim(Aice_library), &
      trim(Alibnames_water(1)), trim(Alibnames_water_sdev(1)), &
      trim(Alibnames_ice(1)), trim(Alibnames_ice_sdev(1)), &
      trim(Alibnames_water(2)), trim(Alibnames_water_sdev(2)), &
      trim(Alibnames_ice(2)), trim(Alibnames_ice_sdev(2)), &
      trim(Alibnames_water(3)), trim(Alibnames_water_sdev(3)), &
      trim(Alibnames_ice(3)), trim(Alibnames_ice_sdev(3) ) )
    guard_mng_ms = 1
  endif

  ! we need to find angle range for the granule, so we can only read a 
  ! limited amount from the library just enough to get the granule processed. 

  ! from the find_granule_angle_range routine we know just how wide the 
  ! angle space we need.  it's fixed for sensor zenith and relative azimuth, 
  ! but not for solar angle. However for MAS<TER> the relative azimuth space 
  ! isn't fixed either and the sensor zenith limits are different. 
  
  !!!!!!!!!!! WDR this logic will set the range needed to read the tables 
  ! just 1 time for test granule and also has refresh logic for test purposes
  !
  !  here we hard-code the ranges for init test
  ! TEMP
  min_rel_azimuth = 0.
  max_rel_azimuth = 180.
  !min_solar_zenith = 50.  ! WDR these are out for table range with margin
  !max_solar_zenith = 80.
  !min_sensor_zenith = 0.
  !max_sensor_zenith = 66.
  !
  !  here is the logic to check to see if a table read is needed
  if( last_min_senz < -1000. ) then
    read_tables = 1
  else
    if( ( min_rel_azimuth < last_min_relaz ) .or. &
        ( max_rel_azimuth > last_max_relaz ) .or. &
        ( min_sensor_zenith < last_min_senz ) .or. &
        ( max_sensor_zenith > last_max_senz ) .or. &
        ( min_solar_zenith < last_min_solz ) .or. &
        ( max_solar_zenith > last_max_solz ) ) then
     read_tables = 1
     endif
  endif
  
  if( read_tables /= 1 ) then
    return
  endif
  print*, "Reading trans tables for this line"
  
  last_min_relaz = min_rel_azimuth
  min_rel_azimuth = last_min_relaz
  last_max_relaz = max_rel_azimuth
  max_rel_azimuth = last_max_relaz
  last_min_senz = min_sensor_zenith - rng_margin_senz
  min_sensor_zenith = last_min_senz
  last_max_senz = max_sensor_zenith + rng_margin_senz
  max_sensor_zenith = last_max_senz
  last_min_solz = min_solar_zenith - rng_margin_solz
  min_solar_zenith = last_min_solz
  last_max_solz = max_solar_zenith + rng_margin_solz
  max_solar_zenith = last_max_solz
  read_tables = 0
  ct_tbl_rd = ct_tbl_rd + 1
  
  print*, "Table update # ", ct_tbl_rd
  print*, "relaz range: ", min_rel_azimuth, max_rel_azimuth
  print*, "solz range: ", min_solar_zenith, max_solar_zenith
  print*, "senz range: ", min_sensor_zenith, max_sensor_zenith
  
  ! start with the ice library first
  status = nf_open( Aice_library, NF_NOWRITE, file_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  start = 1

  ! now get the relative azimuth and resize it as needed
  if (allocated(library_relative_azimuth)) deallocate(library_relative_azimuth)	
  ! resize the array	
  call find_bounds(relativeazimuth_all, min_rel_azimuth, max_rel_azimuth, &
    min_bnd, max_bnd)
  number_relazimuth = max_bnd - min_bnd + 1

  ra_min = min_bnd
  allocate(library_relative_azimuth(number_relazimuth))
  library_relative_azimuth(1:number_relazimuth) = &
    relativeazimuth_all(min_bnd : max_bnd)

  ! now get the solar zenith and resize it as needed
  if (allocated(library_solar_zenith)) deallocate(library_solar_zenith)
  ! resize the array	
  cos_sza_max = cos(max_solar_zenith*deg_to_rad)
  cos_sza_min = cos(min_solar_zenith*deg_to_rad)

  call find_bounds(solarzenith_all, cos_sza_max, cos_sza_min, min_bnd, max_bnd)
  number_solarzenith = max_bnd - min_bnd + 1
		
  sza_min = min_bnd
  allocate(library_solar_zenith(number_solarzenith))
  library_solar_zenith(1:number_solarzenith) = &
    solarzenith_all(min_bnd : max_bnd)
	
  ! now get the sensor zenith and resize it as needed
  if (allocated(library_sensor_zenith)) deallocate(library_sensor_zenith)
  ! resize the array	
  cos_vza_max = cos(max_sensor_zenith*deg_to_rad)
  cos_vza_min = cos(min_sensor_zenith*deg_to_rad)

  call find_bounds(sensorzenith_all, cos_vza_max, cos_vza_min, min_bnd, max_bnd)

  number_sensorzenith = max_bnd - min_bnd + 1
  vza_min = min_bnd
  allocate(library_sensor_zenith(number_sensorzenith))
  library_sensor_zenith(1:number_sensorzenith) = &
    sensorzenith_all(min_bnd : max_bnd)

  ! now get the flux angle and resize it as needed
  if (allocated(library_fluxsensorzenith)) deallocate(library_fluxsensorzenith)
  if (allocated(library_fluxsolarzenith)) deallocate(library_fluxsolarzenith)
  ! resize the array once for solar and once for sensor
  call find_bounds(solarzenith_flux_all, cos_vza_max, cos_vza_min, &
    min_bnd, max_bnd)
  number_fluxsensorzenith = max_bnd - min_bnd + 1
  sensor_min = min_bnd
  sensor_max = max_bnd
  allocate(library_fluxsensorzenith(number_fluxsensorzenith))
  library_fluxsensorzenith(1:number_fluxsensorzenith) = &
    solarzenith_flux_all(min_bnd : max_bnd)

  call find_bounds(solarzenith_flux_all, cos_sza_max, cos_sza_min, &
    min_bnd, max_bnd)
  number_fluxsolarzenith = max_bnd - min_bnd + 1
  solar_min = min_bnd
  solar_max = max_bnd
  allocate(library_fluxsolarzenith(number_fluxsolarzenith))
  library_fluxsolarzenith(1:number_fluxsolarzenith) = &
    solarzenith_flux_all(min_bnd : max_bnd)

  ! Read the 3D arrays

  ! Read the 4D arrays

  start = 1

  allocate(temp_flux2 (number_iceradii, number_wavelengths, number_taus, &
    index_solarzenith_flux))

  if (allocated(flux_up_ice_solar)) deallocate(flux_up_ice_solar)
  if (allocated(flux_up_ice_sensor)) deallocate(flux_up_ice_sensor)
  if (allocated(flux_down_ice_solar)) deallocate(flux_down_ice_solar)
  if (allocated(flux_down_ice_sensor)) deallocate(flux_down_ice_sensor)

  allocate(flux_up_ice_solar (number_fluxsolarzenith, number_taus,number_wavelengths,number_iceradii))
  allocate(flux_down_ice_solar(number_fluxsolarzenith, number_taus,number_wavelengths,number_iceradii)) 
  allocate(flux_up_ice_sensor(number_fluxsensorzenith, number_taus,number_wavelengths,number_iceradii))
  allocate(flux_down_ice_sensor(number_fluxsensorzenith, number_taus,number_wavelengths,number_iceradii)) 

  edge(1) = number_iceradii
  edge(2) = number_wavelengths
  edge(3) = number_taus
  edge(4) = index_solarzenith_flux

	
	
  status = nf_inq_varid( file_id, "ReflectedFlux", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:4), edge(1:4), &
    temp_flux2 )
  call cld_fchk( status, __FILE__, __LINE__ )

  allocate(temp_flux (index_solarzenith_flux, number_taus,number_wavelengths, &
    number_iceradii))
  do i1=1, index_solarzenith_flux
    do i2=1, number_taus
      do i3=1, number_wavelengths
        do i4=1, number_iceradii
			
          temp_flux(i1, i2, i3, i4) = temp_flux2(i4, i3, i2, i1)
				
        end do
      end do
    end do
  end do

	
  flux_up_ice_solar(:,:,:,:) = temp_flux(solar_min:solar_max, :,:,:)
  flux_up_ice_sensor(:,:,:,:) = temp_flux(sensor_min:sensor_max, :,:,:)

  status = nf_inq_varid( file_id, "TransmittedFlux", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:4), edge(1:4), &
    temp_flux2 )

  do i1=1, index_solarzenith_flux
    do i2=1, number_taus
      do i3=1, number_wavelengths
        do i4=1, number_iceradii
				
          temp_flux(i1, i2, i3, i4) = temp_flux2(i4, i3, i2, i1)
				
        end do
      end do
    end do
  end do


  deallocate(temp_flux2)

  flux_down_ice_solar(:,:,:,:) = temp_flux(solar_min:solar_max, :,:,:)
  flux_down_ice_sensor(:,:,:,:) = temp_flux(sensor_min:sensor_max, :,:,:)
	
  deallocate(temp_flux)

  ! Now we finally read the 6D reflectance array
  ! but we don't read the 11um reflectance, as there sure isn't any. 
  ! we do however need the 11um flux to do emissivity over land surfaces.

  status = nf_close( file_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  if (allocated(cloud_emissivity_ice)) deallocate(cloud_emissivity_ice)
  if (allocated(surface_emissivity_ice)) deallocate(surface_emissivity_ice)
	
  allocate(cloud_emissivity_ice( number_sensorzenith, number_taus+1, 2, &
    number_iceradii, number_wind_speed))
  allocate(surface_emissivity_ice(number_sensorzenith, number_taus+1, 2, &
    number_iceradii, number_wind_speed))


  if (allocated(cloud_emissivity_ice_sdev)) &
    deallocate(cloud_emissivity_ice_sdev)
  if (allocated(surface_emissivity_ice_sdev)) &
    deallocate(surface_emissivity_ice_sdev)
	

  allocate(cloud_emissivity_ice_sdev(number_sensorzenith, number_taus+1, 2, &
	 number_iceradii, number_wind_speed))
  allocate(surface_emissivity_ice_sdev(number_sensorzenith, number_taus+1, 2, &
	 number_iceradii, number_wind_speed))

  ! wind speed 3
									
  allocate(temp_emis( number_iceradii, 2, number_taus + 1, &
    number_sensorzenith))

  if (DO_COX_MUNK) then ! DO_COX_MUNK

  do i=1, 3

    status = nf_open( Alibnames_ice(i), NF_NOWRITE, file_id )
    call cld_fchk( status, __FILE__, __LINE__ )
	
    edge(4) = number_sensorzenith
    edge(3) = number_taus+1
    edge(2) = 2
    edge(1) = number_iceradii
	
    start = 1
    start(4) = vza_min
	
    if ( ( c2_sensor_id /= OCI_ID ) .and. ( c2_sensor_id /= OCIS_ID ) ) then 
      ! WDR any emiss is not used
      status = nf_inq_varid( file_id, "CloudEmissivity", var_id )
      call cld_fchk( status, __FILE__, __LINE__ )

      status = nf_get_vara_real( file_id, var_id, start(1:4), edge(1:4), &
        temp_emis )
      call cld_fchk( status, __FILE__, __LINE__ )
      end if

    ! cloud_emissivity_ice(:,:,:,:,i) = temp_emis(:,:,:,:)

    do i1=1, number_sensorzenith
      do i4=1, number_taus+1
        do i5=1, 2
          do i6=1, number_iceradii
						
            cloud_emissivity_ice(i1, i4, i5, i6, i) = temp_emis(i6, i5, i4, i1)
						
          end do
        end do
      end do
    end do

    if ( ( c2_sensor_id /= OCI_ID ) .and. ( c2_sensor_id /= OCIS_ID ) ) then 
      ! WDR any emiss is not used
      status = nf_inq_varid( file_id, "SurfaceEmissivity", var_id )
      call cld_fchk( status, __FILE__, __LINE__ )

      status = nf_get_vara_real( file_id, var_id, start(1:4), edge(1:4), &
        temp_emis )
      end if

    !	surface_emissivity_ice(:,:,:,:,i) = temp_emis(:,:,:,:)

    do i1=1, number_sensorzenith
      do i4=1, number_taus+1
        do i5=1, 2
          do i6=1, number_iceradii
						
            surface_emissivity_ice(i1, i4, i5, i6, i) = temp_emis(i6, i5, i4, i1)
						
          end do
        end do
      end do
    end do


		
    status = nf_close( file_id )
    call cld_fchk( status, __FILE__, __LINE__ )

    status = nf_open( Alibnames_ice_sdev(i), NF_NOWRITE, file_id )
    call cld_fchk( status, __FILE__, __LINE__ )
	
    edge(4) = number_sensorzenith
    edge(3) = number_taus+1
    edge(2) = 2
    edge(1) = number_iceradii
	
    start = 1
    start(4) = vza_min
	
    !  The cloud folks used the wrong SD dataset namd, use 
    !  StdDevCloudEmissivity instead
    if ( ( c2_sensor_id /= OCI_ID ) .and. ( c2_sensor_id /= OCIS_ID ) ) then 
      ! WDR any emiss is not used
      status = nf_inq_varid( file_id, "StdDevCloudEmissivity", var_id )
      if( status .ne. 0 ) then
        status = nf_inq_varid( file_id, "CloudEmissivity", var_id )
        call cld_fchk( status, __FILE__, __LINE__ )
      endif

      status = nf_get_vara_real( file_id, var_id, start(1:4), edge(1:4), &
        temp_emis )
    end if

    !	cloud_emissivity_ice_sdev(:,:,:,:,i) = temp_emis(:,:,:,:)

    do i1=1, number_sensorzenith
      do i4=1, number_taus+1
        do i5=1, 2
          do i6=1, number_iceradii
						
            cloud_emissivity_ice_sdev(i1, i4, i5, i6, i) = &
              temp_emis(i6, i5, i4, i1)
						
          end do
        end do
      end do
    end do
	
    !  Same error here, use StdDevSurfaceEmissivity
    if ( ( c2_sensor_id /= OCI_ID ) .and. ( c2_sensor_id /= OCIS_ID ) ) then 
      ! WDR any emiss is not used
      status = nf_inq_varid( file_id, "StdDevSurfaceEmissivity", var_id )
      if( status .ne. 0 ) then
        status = nf_inq_varid( file_id, "SurfaceEmissivity", var_id )
        call cld_fchk( status, __FILE__, __LINE__ )
      end if

      status = nf_get_vara_real( file_id, var_id, start(1:4), edge(1:4), &
        temp_emis )
    endif

    !	surface_emissivity_ice_sdev(:,:,:,:, i) = temp_emis(:,:,:,:)

	 do i1=1, number_sensorzenith
      do i4=1, number_taus+1
        do i5=1, 2
           do i6=1, number_iceradii
						
            surface_emissivity_ice_sdev(i1, i4, i5, i6, i) = &
              temp_emis(i6, i5, i4, i1)
						
          end do
        end do
      end do
	 end do

		
    status = nf_close( file_id )
    call cld_fchk( status, __FILE__, __LINE__ )

  end do

  endif ! DO_COX_MUNK

  deallocate(temp_emis)

  !  Now that we're done with ice, get started with water...

  !  Read all the water stuff there is to read.
  status = nf_open( Awater_library, NF_NOWRITE, file_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  ! Read the 3D arrays

  ! Read the 4D arrays

  start = 1

  allocate(temp_flux2 (number_waterradii, number_wavelengths, number_taus, &
    index_solarzenith_flux))

  if (allocated(flux_up_water_solar)) deallocate(flux_up_water_solar)
  if (allocated(flux_up_water_sensor)) deallocate(flux_up_water_sensor)
  if (allocated(flux_down_water_solar)) deallocate(flux_down_water_solar)
  if (allocated(flux_down_water_sensor)) deallocate(flux_down_water_sensor)


  allocate(flux_up_water_solar (number_fluxsolarzenith, number_taus,number_wavelengths,number_waterradii))
  allocate(flux_down_water_solar(number_fluxsolarzenith, number_taus,number_wavelengths,number_waterradii)) 
  allocate(flux_up_water_sensor(number_fluxsensorzenith, number_taus,number_wavelengths,number_waterradii))
  allocate(flux_down_water_sensor(number_fluxsensorzenith, number_taus,number_wavelengths,number_waterradii)) 

  edge(1) = number_waterradii
  edge(2) = number_wavelengths
  edge(3) = number_taus
  edge(4) = index_solarzenith_flux
	
  allocate(temp_flux (index_solarzenith_flux, number_taus,number_wavelengths,number_waterradii))

  status = nf_inq_varid( file_id, "ReflectedFlux", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:4),edge(1:4), &
    temp_flux2 )
  call cld_fchk( status, __FILE__, __LINE__ )

  do i1=1, index_solarzenith_flux
    do i2=1, number_taus
      do i3=1, number_wavelengths
        do i4=1, number_waterradii
				
          temp_flux(i1, i2, i3, i4) = temp_flux2(i4, i3, i2, i1)
				
        end do
      end do
    end do
  end do
	
  flux_up_water_solar(:,:,:,:) = temp_flux(solar_min:solar_max, :,:,:)
  flux_up_water_sensor(:,:,:,:) = temp_flux(sensor_min:sensor_max, :,:,:)

  status = nf_inq_varid( file_id, "TransmittedFlux", var_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  status = nf_get_vara_real( file_id, var_id, start(1:4), edge(1:4), &
    temp_flux2 )
  call cld_fchk( status, __FILE__, __LINE__ )

  do i1=1, index_solarzenith_flux
    do i2=1, number_taus
      do i3=1, number_wavelengths
        do i4=1, number_waterradii
				
          temp_flux(i1, i2, i3, i4) = temp_flux2(i4, i3, i2, i1)

        end do
      end do
    end do
  end do

  deallocate(temp_flux2)

  flux_down_water_solar(:,:,:,:) = temp_flux(solar_min:solar_max, :,:,:)
  flux_down_water_sensor(:,:,:,:) = temp_flux(sensor_min:sensor_max, :,:,:)

  deallocate(temp_flux)

  ! Now we finally read the 6D reflectance array
  ! now, this is done on-demand by the dim_mgr
  status = nf_close( file_id )
  call cld_fchk( status, __FILE__, __LINE__ )

  if (allocated(cloud_emissivity_water)) deallocate(cloud_emissivity_water)
  if (allocated(surface_emissivity_water)) deallocate(surface_emissivity_water)
	
  allocate(cloud_emissivity_water(number_sensorzenith, number_taus+1, 2, &
			 number_waterradii, number_wind_speed))
  allocate(surface_emissivity_water(number_sensorzenith, number_taus+1, 2, &
			 number_waterradii, number_wind_speed))

  if (allocated(cloud_emissivity_water_sdev)) &
    deallocate(cloud_emissivity_water_sdev)
  if (allocated(surface_emissivity_water_sdev)) &
    deallocate(surface_emissivity_water_sdev)
  allocate(cloud_emissivity_water_sdev(number_sensorzenith, number_taus+1, 2, &
	 number_waterradii, number_wind_speed))
  allocate(surface_emissivity_water_sdev(number_sensorzenith, &
    number_taus+1, 2, number_waterradii, number_wind_speed))


  ! wind speed 3

  allocate(temp_emis( number_waterradii, 2, number_taus + 1, &
    number_sensorzenith))
	
  if (DO_COX_MUNK) then ! DO_COX_MUNK
  do i=1, 3
		
    status = nf_open( Alibnames_water(i), NF_NOWRITE, file_id )
    call cld_fchk( status, __FILE__, __LINE__ )
	
    edge(4) = number_sensorzenith
    edge(3) = number_taus+1
    edge(2) = 2
    edge(1) = number_waterradii
	
    start = 1
    start(4) = vza_min
	
    if ( ( c2_sensor_id /= OCI_ID ) .and. ( c2_sensor_id /= OCIS_ID ) ) then 
      ! WDR any emiss is not used
      status = nf_inq_varid( file_id, "CloudEmissivity", var_id )
      call cld_fchk( status, __FILE__, __LINE__ )

      status = nf_get_vara_real( file_id, var_id, start(1:4), edge(1:4), &
        temp_emis )
      call cld_fchk( status, __FILE__, __LINE__ )
    end if
	
    !	cloud_emissivity_water(:,:,:,:, i) = temp_emis(:,:,:,:)

    do i1=1, number_sensorzenith
      do i4=1, number_taus+1
        do i5=1, 2
          do i6=1, number_waterradii
						
           cloud_emissivity_water(i1, i4, i5, i6, i) = temp_emis(i6, i5, i4, i1)
						
          end do
        end do
      end do
    end do

    if ( ( c2_sensor_id /= OCI_ID ) .and. ( c2_sensor_id /= OCIS_ID ) ) then 
      ! WDR any emiss is not used	
      status = nf_inq_varid( file_id, "SurfaceEmissivity", var_id )
      call cld_fchk( status, __FILE__, __LINE__ )

      status = nf_get_vara_real( file_id, var_id, start(1:4), edge(1:4), &
        temp_emis )
      call cld_fchk( status, __FILE__, __LINE__ )
    end if

    !	surface_emissivity_water(:,:,:,:, i) = temp_emis(:,:,:,:)

    do i1=1, number_sensorzenith
      do i4=1, number_taus+1
        do i5=1, 2
          do i6=1, number_waterradii
						
            surface_emissivity_water(i1, i4, i5, i6, i) = &
              temp_emis(i6, i5, i4, i1)
						
          end do
        end do
      end do
    end do
	
    status = nf_close( file_id )
    call cld_fchk( status, __FILE__, __LINE__ )

    status = nf_open( Alibnames_water_sdev(i), NF_NOWRITE, file_id )
    call cld_fchk( status, __FILE__, __LINE__ )
	
	 edge(4) = number_sensorzenith
	 edge(3) = number_taus+1
	 edge(2) = 2
	 edge(1) = number_waterradii
	
	 start = 1
    start(4) = vza_min
	
    !  WDR the cloud group put the SD data in SurfaceEmissivity, but the
    !  netcdf files need this as StdDevSurfaceEmissivity, so look for that
    !  first
    if ( ( c2_sensor_id /= OCI_ID ) .and. ( c2_sensor_id /= OCIS_ID ) ) then 
      ! WDR any emiss is not used
      status = nf_inq_varid( file_id, "StdDevCloudEmissivity", var_id )
      if( status .ne. 0 ) then
        status = nf_inq_varid( file_id, "CloudEmissivity", var_id )
        call cld_fchk( status, __FILE__, __LINE__ )
      end if

      status = nf_get_vara_real( file_id, var_id, start(1:4), edge(1:4), &
        temp_emis )
      call cld_fchk( status, __FILE__, __LINE__ )
    endif
	
    !	cloud_emissivity_water_sdev(:,:,:,:,i) = temp_emis(:,:,:,:)

    do i1=1, number_sensorzenith
      do i4=1, number_taus+1
        do i5=1, 2
          do i6=1, number_waterradii
						
            cloud_emissivity_water_sdev(i1, i4, i5, i6, i) = &
              temp_emis(i6, i5, i4, i1)
						
          end do
        end do
      end do
    end do
	
    !  WDR the cloud group put the SD data in SurfaceEmissivity, but the 
    !  netcdf files need this as StdDevSurfaceEmissivity, so look for that 
    !  first

    if ( ( c2_sensor_id /= OCI_ID ) .and. ( c2_sensor_id /= OCIS_ID ) ) then 
      ! WDR any emiss is not used
      status = nf_inq_varid( file_id, "StdDevSurfaceEmissivity,", var_id )
      if( status .ne. 0 ) then
        status = nf_inq_varid( file_id, "SurfaceEmissivity", var_id )
        call cld_fchk( status, __FILE__, __LINE__ )
      end if

      status = nf_get_vara_real( file_id, var_id, start(1:4), edge(1:4), &
        temp_emis )
      call cld_fchk( status, __FILE__, __LINE__ )
    endif

    ! surface_emissivity_water_sdev(:,:,:,:,i) = temp_emis(:,:,:,:)

    do i1=1, number_sensorzenith
      do i4=1, number_taus+1
        do i5=1, 2
          do i6=1, number_waterradii
						
            surface_emissivity_water_sdev(i1, i4, i5, i6, i) = &
              temp_emis(i6, i5, i4, i1)
						
          end do
        end do
      end do
    end do
	
    status = nf_close( file_id )
    call cld_fchk( status, __FILE__, __LINE__ )

    ! call system_clock(end_time, crate, cmax)
    ! print*, "time elapsed: ", end_time - start_time

  end do
  endif ! DO_COX_MUNK


  deallocate(temp_emis)


  !	print*, "lamb: ", library_reflectance_water( 1, 1, 1, :, 4, 8, 1)
  !	print*, "ws3: ",library_reflectance_water( 1, 1, 1, :, 4, 8, 2)
  !	print*, "ws7: ",library_reflectance_water( 1, 1, 1, :, 4, 8, 3)
  !	print*, "ws15: ",library_reflectance_water( 1, 1, 1, :, 4, 8, 4)
  
  !	print*, "lamb: ", library_reflectance_ice( 1, 1, 1, :, 4, 6, 1)
  !	print*, "ws3: ",library_reflectance_ice( 1, 1, 1, :, 4, 6, 2)
  !	print*, "ws7: ",library_reflectance_ice( 1, 1, 1, :, 4, 6, 3)
  !	print*, "ws15: ",library_reflectance_ice( 1, 1, 1, :, 4, 6, 4)

  end subroutine readlibraries_extra
 
  subroutine allocate_model(st_iterX, st_iterY, grid_wid, grid_ht)
 
 	use core_arrays, only: model_info, c2_model_info
 	use science_parameters, only: model_levels
 	use GeneralAuxType
 
	integer, intent(in) :: st_iterX, st_iterY, grid_wid, grid_ht
	logical :: do_it
	integer :: i, j

#ifdef GEOS5
	do_it = .true.
#else
	do_it = .false.
   if (iterationX == st_iterX .and. iterationY == st_iterY) do_it = .true.
#endif

	if (do_it) then 
       !  WDR
       !  no need to allocate c2_model_info except the contained arrays
       ! WDR re-entrant change
       if( .not. allocated( c2_model_info%mixr_profile ) ) then
         allocate( c2_model_info%mixr_profile(model_levels))
         allocate( c2_model_info%temp_profile(model_levels))
         allocate( c2_model_info%height_profile(model_levels))
         allocate( c2_model_info%pressure_profile(model_levels))
         allocate( c2_model_info%o3_profile(model_levels))
         endif
      if( .not. allocated( model_info  ) ) then
    	   allocate(model_info(grid_wid, grid_ht))
    		do i=1, grid_wid
    			do j=1, grid_ht
    		
    				allocate(model_info(i,j)%mixr_profile(model_levels))
    				allocate(model_info(i,j)%temp_profile(model_levels))
    				allocate(model_info(i,j)%height_profile(model_levels))
    				allocate(model_info(i,j)%pressure_profile(model_levels))
    		
    			end do
    		end do
      end if
   end if
   
 end subroutine allocate_model
 
 
 

 subroutine allocate_arrays ( edge, meas_edge, st_iterX, st_iterY, status ) 

   use GeneralAuxType
   use core_arrays
   use nonscience_parameters
   use libraryarrays
   use mod06_run_settings
   use science_parameters
   use specific_other
   use ch_xfr, only: OCI_ID, OCIS_ID, c2_sensor_id
   
   implicit none
   
   integer, dimension(:), intent(in)     ::  edge, meas_edge
   integer,               intent (out)   :: status
   integer, intent(in) :: st_iterX, st_iterY

   integer                               :: checkvariable
   integer                               :: xdimension, meas_xdimension, ydimension, number_of_bands

   integer                               :: model_layers
   integer :: i, j
   
   logical       :: allocation_status, array_size_change
   
   status = 0
   model_layers = model_levels
   number_of_bands = size(bands)


   meas_xdimension = meas_edge(1)
   xdimension = edge(1)
   ydimension =  edge(2)
 
!  Use optical_thickness array as a size/allocation test to save time.
   allocation_status = allocated(optical_thickness_16_final) 

   array_size_change = .false.
   
   if (allocation_status) then
      if (size(optical_thickness_16_final,1) /= xdimension .or. &
          size(optical_thickness_16_final,2) /= ydimension .or. &
		  size(band_measurements, 1) /= meas_xdimension) then
         array_size_change  = .true.
      endif
   endif

#ifdef GEOS5
	if (allocated(model_info)) then 

		do i=1, grid_xsize
			do j=1, grid_ysize
		
				deallocate(model_info(i,j)%mixr_profile)
				deallocate(model_info(i,j)%temp_profile)
				deallocate(model_info(i,j)%height_profile)
				deallocate(model_info(i,j)%pressure_profile)
		
			end do
		end do

	    deallocate(model_info)
	endif
#endif

#ifndef GEOS5
	call allocate_model(st_iterX, st_iterY, grid_xsize, grid_ysize)
#endif

!  Core retrieval arrays   
   if(array_size_change) then

		deallocate(snow_cover)
		deallocate(cloud_height_method)
		deallocate(irw_temperature)
		deallocate(optical_thickness_final, stat = checkvariable)

		deallocate(optical_thickness_1621_final, stat = checkvariable)
		deallocate(effective_radius_21_final, stat = checkvariable)
		deallocate(effective_radius_1621_final, stat = checkvariable)
		deallocate(liquid_water_path_1621, stat = checkvariable)
		deallocate(liquid_water_path, stat = checkvariable)

		deallocate(optical_thickness_final_PCL, stat = checkvariable)
		deallocate(optical_thickness_1621_final_PCL, stat = checkvariable)
		deallocate(effective_radius_21_final_PCL, stat = checkvariable)
		deallocate(effective_radius_1621_final_PCL, stat = checkvariable)
		deallocate(liquid_water_path_1621_PCL, stat = checkvariable)
		deallocate(liquid_water_path_PCL, stat = checkvariable)

		deallocate(optical_thickness_error, stat = checkvariable)
		deallocate(effective_radius_21_error, stat = checkvariable)

		deallocate(liquid_water_path_error, stat = checkvariable)
		deallocate(optical_thickness_1621_error, stat = checkvariable)
		deallocate(effective_radius_1621_error, stat = checkvariable)
		deallocate(liquid_water_path_1621_error, stat = checkvariable)

      !  WDR only do 2.2 um band if OCI
      if( ( c2_sensor_id == OCI_ID ) .or. ( c2_sensor_id == OCIS_ID ) ) then
        deallocate(optical_thickness_22_final, stat = checkvariable)
        deallocate(effective_radius_22_final, stat = checkvariable)
        deallocate(liquid_water_path_22, stat = checkvariable)
        deallocate(optical_thickness_22_final_PCL, stat = checkvariable)
        deallocate(effective_radius_22_final_PCL, stat = checkvariable)
        deallocate(liquid_water_path_22_PCL, stat = checkvariable)
        deallocate(optical_thickness_22_error, stat = checkvariable)
        deallocate(effective_radius_22_error, stat = checkvariable)
        deallocate(liquid_water_path_22_error, stat = checkvariable)
        deallocate(failure_metric_22, stat = checkvariable)
        endif
		
		deallocate(cloud_mask_SPI, stat = checkvariable)
	
		deallocate(failure_metric, stat = checkvariable)
		deallocate(failure_metric_1621, stat = checkvariable)

		deallocate(atm_corr_refl, stat = checkvariable)

		deallocate(optical_thickness_16_final, stat = checkvariable)
		deallocate(optical_thickness_37_final, stat = checkvariable)
		deallocate(effective_radius_16_final, stat = checkvariable)
		deallocate(effective_radius_37_final, stat = checkvariable)
		deallocate(liquid_water_path_16, stat = checkvariable)
		deallocate(liquid_water_path_37, stat = checkvariable)

		deallocate(optical_thickness_16_final_PCL, stat = checkvariable)
		deallocate(optical_thickness_37_final_PCL, stat = checkvariable)
		deallocate(effective_radius_16_final_PCL, stat = checkvariable)
		deallocate(effective_radius_37_final_PCL, stat = checkvariable)
		deallocate(liquid_water_path_16_PCL, stat = checkvariable)
		deallocate(liquid_water_path_37_PCL, stat = checkvariable)

		deallocate(optical_thickness_16_error, stat = checkvariable)
		deallocate(optical_thickness_37_error, stat = checkvariable)
		deallocate(effective_radius_16_error, stat = checkvariable)
		deallocate(effective_radius_37_error, stat = checkvariable)

		deallocate(liquid_water_path_16_error, stat = checkvariable)
		deallocate(liquid_water_path_37_error, stat = checkvariable)
		deallocate(cloud_layer_flag, stat = checkvariable)
		deallocate(ml_test_flag, stat = checkvariable)
		deallocate(CSR_flag_array, stat = checkvariable)
		deallocate(precip_water_094, stat = checkvariable)

!******************
!		deallocate(clear_sky_btemp)
!		deallocate(clear_sky_rad)

		deallocate(latitude, stat = checkvariable)
		deallocate(longitude, stat = checkvariable)
		deallocate(sensor_zenith_angle, stat = checkvariable)
		deallocate(solar_zenith_angle, stat = checkvariable)
		deallocate(relative_azimuth_angle, stat = checkvariable)
		deallocate(band_measurements, stat = checkvariable)
		deallocate(band_uncertainty, stat = checkvariable)
	    deallocate(sensor_azimuth_angle, stat = checkvariable)
   		deallocate(solar_azimuth_angle, stat = checkvariable)
	
		deallocate(surface_albedo, stat = checkvariable)
		deallocate(surface_temperature, stat = checkvariable)
		deallocate(surface_emissivity_land, stat = checkvariable)
		deallocate(cloud_top_temperature, stat = checkvariable)
		deallocate(cloud_top_height, stat = checkvariable)
		deallocate(cloud_top_pressure, stat = checkvariable)
	    deallocate(cloud_phase_infrared, stat = checkvariable)
		deallocate(abovecloud_watervapor, stat = checkvariable)
		deallocate(column_ozone, stat = checkvariable)
		deallocate(cloudmask, stat = checkvariable)
		
		deallocate(cloud_mask_SPI, stat = checkvariable)

		deallocate(failure_metric_16, stat = checkvariable)
		deallocate(failure_metric_37, stat = checkvariable)	
		
		deallocate(cloudsummary, stat = checkvariable)


		call deallocate_extra
		deallocate(processing_information, stat = checkvariable)

	endif


	if (.not. allocated(optical_thickness_final)) then 

		
		allocate (optical_thickness_final(xdimension,ydimension), stat = checkvariable)	 
		allocate (optical_thickness_1621_final(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_21_final(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_1621_final(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_1621(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path(xdimension,ydimension), stat = checkvariable)

		allocate (optical_thickness_final_PCL(xdimension,ydimension), stat = checkvariable)	 
		allocate (optical_thickness_1621_final_PCL(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_21_final_PCL(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_1621_final_PCL(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_1621_PCL(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_PCL(xdimension,ydimension), stat = checkvariable)

		allocate (optical_thickness_error(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_21_error(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_error(xdimension,ydimension), stat = checkvariable)
		allocate (optical_thickness_1621_error(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_1621_error(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_1621_error(xdimension,ydimension), stat = checkvariable)

		allocate(failure_metric(xdimension, ydimension), stat = checkvariable)
		allocate(failure_metric_1621(xdimension, ydimension), stat = checkvariable)
      ! WDR only do 2.2 um for OCI
      if( ( c2_sensor_id == OCI_ID ) .or. ( c2_sensor_id == OCIS_ID ) ) then
        allocate (optical_thickness_22_final(xdimension,ydimension), &
          stat = checkvariable)
        allocate (effective_radius_22_final(xdimension,ydimension), &
          stat = checkvariable)
        allocate (liquid_water_path_22(xdimension,ydimension), &
          stat = checkvariable)
        allocate (optical_thickness_22_final_PCL(xdimension,ydimension), &
          stat = checkvariable)	 
        allocate (effective_radius_22_final_PCL(xdimension,ydimension), &
          stat = checkvariable)
        allocate (liquid_water_path_22_PCL(xdimension,ydimension), &
          stat = checkvariable)
        allocate (optical_thickness_22_error(xdimension,ydimension), &
          stat = checkvariable)
        allocate (effective_radius_22_error(xdimension,ydimension), &
          stat = checkvariable)
        allocate (liquid_water_path_22_error(xdimension,ydimension), &
          stat = checkvariable)
        allocate(failure_metric_22(xdimension, ydimension), &
          stat = checkvariable)
        endif

 		allocate(atm_corr_refl(set_albedo_bands, xdimension, ydimension), stat=checkvariable)
		allocate (cloud_mask_SPI(2,xdimension,ydimension), stat = checkvariable)
      cloud_mask_SPI = 0  ! WDR to set to a known


		allocate (optical_thickness_16_final(xdimension,ydimension), stat = checkvariable)	 
		allocate (optical_thickness_37_final(xdimension,ydimension), stat = checkvariable)	 
		allocate (effective_radius_16_final(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_37_final(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_16(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_37(xdimension,ydimension), stat = checkvariable)

		allocate (optical_thickness_16_final_PCL(xdimension,ydimension), stat = checkvariable)	 
		allocate (optical_thickness_37_final_PCL(xdimension,ydimension), stat = checkvariable)	 
		allocate (effective_radius_16_final_PCL(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_37_final_PCL(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_16_PCL(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_37_PCL(xdimension,ydimension), stat = checkvariable)
				
		allocate (optical_thickness_16_error(xdimension,ydimension), stat = checkvariable)
		allocate (optical_thickness_37_error(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_16_error(xdimension,ydimension), stat = checkvariable)
		allocate (effective_radius_37_error(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_16_error(xdimension,ydimension), stat = checkvariable)
		allocate (liquid_water_path_37_error(xdimension,ydimension), stat = checkvariable)
		allocate (cloud_layer_flag(xdimension,ydimension), stat = checkvariable)
		allocate (ml_test_flag(xdimension,ydimension), stat = checkvariable)
		allocate (CSR_flag_array(xdimension,ydimension), stat = checkvariable)
		allocate(precip_water_094(xdimension,ydimension), stat = checkvariable)
		allocate(snow_cover(xdimension,ydimension), stat = checkvariable)
 
		allocate(failure_metric_16(xdimension, ydimension), stat = checkvariable)
		allocate(failure_metric_37(xdimension, ydimension), stat = checkvariable)
 
 		 
!  Core input data arrays

		allocate (latitude(xdimension,ydimension), stat = checkvariable)
		allocate (longitude(xdimension,ydimension), stat = checkvariable)
		allocate (sensor_zenith_angle(xdimension,ydimension), stat = checkvariable)
		allocate (solar_zenith_angle(xdimension,ydimension), stat = checkvariable)
		allocate (sensor_azimuth_angle(xdimension,ydimension), stat = checkvariable)
		allocate (solar_azimuth_angle(xdimension,ydimension), stat = checkvariable)
		allocate (relative_azimuth_angle(xdimension,ydimension), stat = checkvariable)

		allocate (band_measurements(meas_xdimension,number_of_bands, ydimension), stat = checkvariable)
		allocate (band_uncertainty(meas_xdimension,set_albedo_bands, ydimension), stat = checkvariable)


!  Ancillary data array

		allocate (surface_albedo(xdimension,ydimension,set_albedo_bands), stat = checkvariable)
		allocate (surface_temperature(xdimension,ydimension), stat = checkvariable)
		allocate (cloud_height_method(xdimension,ydimension), stat = checkvariable)
		allocate (cloud_phase_infrared(xdimension,ydimension), stat = checkvariable)
		allocate (irw_temperature(xdimension,ydimension), stat = checkvariable)
		allocate (surface_emissivity_land(xdimension,ydimension, 2), stat = checkvariable)
		allocate (cloud_top_temperature(xdimension,ydimension), stat = checkvariable)
		allocate (cloud_top_height(xdimension,ydimension), stat = checkvariable)
		allocate (cloud_top_pressure(xdimension,ydimension), stat = checkvariable)
		allocate (abovecloud_watervapor(xdimension,ydimension), stat = checkvariable)
		allocate (column_ozone(xdimension,ydimension), stat = checkvariable)
      column_ozone = 0 ! WDR-UIV
		allocate (cloudmask(xdimension,ydimension), stat = checkvariable)
		allocate (cloud_mask_SPI(2,xdimension,ydimension), stat = checkvariable)
      cloud_mask_SPI = 0  ! WDR to set to a known

!******************
!		allocate (clear_sky_btemp(xdimension, ydimension, 2))
!		allocate (clear_sky_rad(xdimension, ydimension, 2))
		

!  QA and Processing arrays

		allocate (cloudsummary(xdimension,ydimension), stat = checkvariable)

!  atmospheric correction

!		allocate (transmittance_stddev(xdimension,ydimension,7), stat = checkvariable)
!		allocate (transmittance_twoway(xdimension,ydimension,7), stat = checkvariable)
!		allocate (meandelta_trans(xdimension,ydimension,7), stat = checkvariable)
!		allocate (thermal_correction_twoway(xdimension,ydimension,2), stat = checkvariable)
!		allocate (thermal_correction_oneway(xdimension,ydimension,2), stat = checkvariable)
!		allocate (thermal_correction_bands(2), stat = checkvariable)
  
		call allocate_extra(xdimension, ydimension)

		allocate (processing_information(xdimension,ydimension), stat = checkvariable)
		

    endif
  
end subroutine allocate_arrays


integer function findpoint( vector, value)
   use GeneralAuxType

   implicit none

   real(single) ,intent(in) :: vector(:)
   real(single) ,intent(in) :: value
   integer                  ::temp(1)

   real(single) ,allocatable   :: localvector(:)

   allocate(localvector(size(vector)))

   localvector = vector
   call realsingle_s_where_equal(localvector,value)

   temp = maxloc(localvector)
   findpoint = temp(1)
   deallocate(localvector)

end function findpoint


subroutine init_qualitydata
   use GeneralAuxType
   use core_arrays, only: processing_information

!   product quality and retrieval processing QA flags
    processing_information(:,:)%optical_thickness_GC = 0_integer_onebyte
    processing_information(:,:)%optical_thickness_outofbounds = 0_integer_onebyte
    processing_information(:,:)%effective_radius_GC = 0_integer_onebyte
    processing_information(:,:)%water_path_GC = 0_integer_onebyte
    processing_information(:,:)%rayleigh_correction = 0_integer_onebyte
    processing_information(:,:)%path_and_outcome = 0_integer_onebyte
    processing_information(:,:)%path_and_outcome_1621 = 0_integer_onebyte

    processing_information(:,:)%path_and_outcome_16 = 0_integer_onebyte
    processing_information(:,:)%path_and_outcome_16_PCL = 0_integer_onebyte

    processing_information(:,:)%path_and_outcome_37 = 0_integer_onebyte
    processing_information(:,:)%path_and_outcome_37_PCL = 0_integer_onebyte

    processing_information(:,:)%path_and_outcome_22 = 0_integer_onebyte
    processing_information(:,:)%path_and_outcome_22_PCL = 0_integer_onebyte

    processing_information(:,:)%path_and_outcome_PCL = 0_integer_onebyte
    processing_information(:,:)%path_and_outcome_1621_PCL = 0_integer_onebyte

    processing_information(:,:)%band_used_for_optical_thickness = 0_integer_onebyte
    processing_information(:,:)%optical_thickness_1621_GC = 0_integer_onebyte
    processing_information(:,:)%effective_radius_1621_GC = 0_integer_onebyte
    processing_information(:,:)%water_path_1621_GC= 0_integer_onebyte
    processing_information(:,:)%multi_layer_cloud = 0_integer_onebyte
    processing_information(:,:)%CSR_flag = 0_integer_onebyte
    processing_information(:,:)%ml_test_mark = 0_integer_onebyte

#ifdef SEVIRI_INST
    processing_information(:,:)%Tc_override = 0_integer_onebyte
#endif

end subroutine init_qualitydata

subroutine deallocate_cleanup(status)

	use core_arrays
	use libraryarrays
	use libraryinterpolates
	use specific_other
	use interpolate_libraries
	use science_parameters
   use ch_xfr, only: OCI_ID, OCIS_ID, c2_sensor_id

	integer, intent(inout) :: status
	integer :: checkvariable

	integer :: i, j

	status = 0

   if (allocated(model_info)) then 

		do i=1, grid_xsize
			do j=1, grid_ysize
					
				deallocate(model_info(i,j)%mixr_profile)
				deallocate(model_info(i,j)%temp_profile)
				deallocate(model_info(i,j)%height_profile)
				deallocate(model_info(i,j)%pressure_profile)
		
			end do
		end do

    	deallocate(model_info)

	endif

		deallocate(optical_thickness_final, stat = checkvariable)
		deallocate(optical_thickness_1621_final, stat = checkvariable)
		deallocate(effective_radius_21_final, stat = checkvariable)
		deallocate(effective_radius_1621_final, stat = checkvariable)
		deallocate(liquid_water_path_1621, stat = checkvariable)
		deallocate(liquid_water_path, stat = checkvariable)

		deallocate(optical_thickness_final_PCL, stat = checkvariable)
		deallocate(optical_thickness_1621_final_PCL, stat = checkvariable)
		deallocate(effective_radius_21_final_PCL, stat = checkvariable)
		deallocate(effective_radius_1621_final_PCL, stat = checkvariable)
		deallocate(liquid_water_path_1621_PCL, stat = checkvariable)
		deallocate(liquid_water_path_PCL, stat = checkvariable)

		deallocate(optical_thickness_error, stat = checkvariable)
		deallocate(effective_radius_21_error, stat = checkvariable)

		deallocate(liquid_water_path_error, stat = checkvariable)
		deallocate(optical_thickness_1621_error, stat = checkvariable)
		deallocate(effective_radius_1621_error, stat = checkvariable)
		deallocate(liquid_water_path_1621_error, stat = checkvariable)
		
		deallocate(cloud_mask_SPI, stat = checkvariable)
	
		deallocate(failure_metric, stat = checkvariable)
		deallocate(failure_metric_1621, stat = checkvariable)

		deallocate(atm_corr_refl, stat = checkvariable)
      ! WDR only for OCI
      if( ( c2_sensor_id == OCI_ID ) .or. ( c2_sensor_id == OCIS_ID ) )then
        deallocate(effective_radius_22_final, stat = checkvariable)
        deallocate(optical_thickness_22_final_PCL, stat = checkvariable)
        deallocate(effective_radius_22_final_PCL, stat = checkvariable)
        deallocate(liquid_water_path_22_PCL, stat = checkvariable)
        deallocate(optical_thickness_22_error, stat = checkvariable)
        deallocate(effective_radius_22_error, stat = checkvariable)
        deallocate(liquid_water_path_22_error, stat = checkvariable)
        deallocate(failure_metric_22, stat = checkvariable)
        endif

!  Core retrieval arrays   

		deallocate(optical_thickness_16_final, stat = checkvariable)
		deallocate(optical_thickness_37_final, stat = checkvariable)
		deallocate(effective_radius_16_final, stat = checkvariable)
		deallocate(effective_radius_37_final, stat = checkvariable)

     deallocate(liquid_water_path_16, stat = checkvariable)
     deallocate(liquid_water_path_37, stat = checkvariable)


		deallocate(optical_thickness_16_final_PCL, stat = checkvariable)
		deallocate(optical_thickness_37_final_PCL, stat = checkvariable)
		deallocate(effective_radius_16_final_PCL, stat = checkvariable)
		deallocate(effective_radius_37_final_PCL, stat = checkvariable)
		deallocate(liquid_water_path_16_PCL, stat = checkvariable)
		deallocate(liquid_water_path_37_PCL, stat = checkvariable)



     deallocate(optical_thickness_16_error, stat = checkvariable)
     deallocate(optical_thickness_37_error, stat = checkvariable)
     deallocate(effective_radius_16_error, stat = checkvariable)
     deallocate(effective_radius_37_error, stat = checkvariable)
	 deallocate(liquid_water_path_16_error, stat = checkvariable)
	 deallocate(liquid_water_path_37_error, stat = checkvariable)

     deallocate(cloud_layer_flag, stat = checkvariable)
     deallocate(ml_test_flag, stat = checkvariable)
     deallocate(CSR_flag_array, stat = checkvariable)
     deallocate(precip_water_094, stat = checkvariable)
     
	deallocate(failure_metric_16, stat = checkvariable)
	deallocate(failure_metric_37, stat = checkvariable)
     	
!  Core input data arrays
     deallocate(latitude, stat = checkvariable)
     deallocate(longitude, stat = checkvariable)
     deallocate(sensor_zenith_angle, stat = checkvariable)
     deallocate(solar_zenith_angle, stat = checkvariable)
     deallocate(relative_azimuth_angle, stat = checkvariable)
     deallocate(band_measurements, stat = checkvariable)
     deallocate(band_uncertainty, stat = checkvariable)
	    deallocate(sensor_azimuth_angle, stat = checkvariable)
   		deallocate(solar_azimuth_angle, stat = checkvariable)
	 
     deallocate(surface_albedo, stat = checkvariable)
     deallocate(surface_temperature, stat = checkvariable)
     deallocate(surface_emissivity_land, stat = checkvariable)
     deallocate(cloud_top_temperature, stat = checkvariable)
     deallocate(cloud_top_pressure, stat = checkvariable)
	 deallocate(cloud_phase_infrared, stat = checkvariable)
     deallocate(abovecloud_watervapor, stat = checkvariable)
     deallocate(column_ozone, stat = checkvariable)
     deallocate(cloudmask, stat = checkvariable)
	 deallocate(cloud_mask_SPI, stat = checkvariable)
	 deallocate(snow_cover, stat = checkvariable)
!     deallocate(model_height_profile, stat = checkvariable)
!     deallocate(model_temp_profile, stat = checkvariable)
!     deallocate(model_water_profile, stat = checkvariable)


!***********
!	deallocate(clear_sky_btemp)
!	deallocate(clear_sky_rad)

!  QA and Processing arrays
     deallocate(cloudsummary, stat = checkvariable)
     deallocate(processing_information, stat = checkvariable)

!  atmospheric correction

!     deallocate(transmittance_stddev, stat = checkvariable)
!     deallocate(transmittance_twoway, stat = checkvariable)
!     deallocate(meandelta_trans, stat = checkvariable)
!     deallocate(thermal_correction_twoway, stat = checkvariable)
!     deallocate(thermal_correction_oneway, stat = checkvariable)
!     deallocate(thermal_correction_bands, stat = checkvariable)

		deallocate(cloud_height_method)
		deallocate(irw_temperature)

! library arrays

     deallocate(transmit_correction_table)
     deallocate(transmit_stddev_table)
 
	 deallocate(ice_radii)
	 deallocate(library_taus)
	 deallocate(water_radii)

!  finally. Let's allocate the arrays we'll need later 
     deallocate(extinction_water, asymmetry_water, singlescattering_water, &
					truncation_factor_water, phase_fun_norm_constant_water)
     deallocate(spherical_albedo_water)
   
!  ice arrays, library, not interpolated


 	 deallocate(extinction_ice, asymmetry_ice, singlescattering_ice, &
					truncation_factor_ice, phase_fun_norm_constant_ice)
     deallocate(spherical_albedo_ice)
   

	 deallocate(rayleigh_tau, aerosol_tau, aerosol_asym, aerosol_ssa)
   
  
!  arrays for library interpolation
     deallocate(int_reflectance_water, int_reflectance_water_sdev)
	 deallocate( int_fluxupwater_sensor, int_fluxdownwater_solar,int_fluxdownwater_sensor)
     deallocate(int_reflectance_ice , int_reflectance_ice_sdev )
	 deallocate ( int_fluxupice_sensor, int_fluxdownice_solar,int_fluxdownice_sensor)
	
     deallocate(int_cloud_emis_water_wspeed, int_surface_emis_water_wspeed)
	 deallocate(int_cloud_emis_ice_wspeed, int_surface_emis_ice_wspeed)

   deallocate(int_reflectance_water_wspeed)
   deallocate(int_reflectance_ice_wspeed)
   deallocate(int_refl_water_sdev_wspeed)
   deallocate(int_refl_ice_sdev_wspeed)

	 deallocate(int_cloud_emissivity_ice, int_surface_emissivity_ice)
	 deallocate(int_cloud_emissivity_ice_sdev, int_surface_emissivity_ice_sdev)



     deallocate(int_cloud_emis_water_sdev_wspeed, int_surface_emis_water_sdev_wspeed)
	 deallocate(int_cloud_emis_ice_sdev_wspeed, int_surface_emis_ice_sdev_wspeed)

	deallocate(rayleigh_liq, rayleigh_ice)
	
! phase function information
	 deallocate(phase_angles_water)
	 deallocate(phase_angles_ice)
     deallocate(phase_funcs_water)
     deallocate(phase_funcs_ice)
	 
	 call deallocate_extra
	 call lib_clean

	if (no_valid_data == 0) then 

  	 deallocate(library_solar_zenith)
	 deallocate(library_sensor_zenith)
	 deallocate(library_relative_azimuth)
	 deallocate(library_fluxsolarzenith)
     deallocate(library_fluxsensorzenith)

     deallocate( solarzenith_all, sensorzenith_all, &
                    relativeazimuth_all, solarzenith_flux_all)
	
	 deallocate(cloud_emissivity_ice, surface_emissivity_ice)
	 deallocate(cloud_emissivity_ice_sdev, surface_emissivity_ice_sdev)

     deallocate(flux_up_ice_solar, flux_down_ice_solar)
     deallocate(flux_up_ice_sensor, flux_down_ice_sensor)

     deallocate(flux_up_water_solar, flux_down_water_solar)
     deallocate(flux_up_water_sensor, flux_down_water_sensor)
 
	 deallocate(cloud_emissivity_water, surface_emissivity_water)
	 deallocate(cloud_emissivity_water_sdev, surface_emissivity_water_sdev)

	endif

end subroutine deallocate_cleanup


end module modis_frontend_module
