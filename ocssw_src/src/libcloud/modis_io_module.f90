 module modis_io_module

 ! WDR use shared_io_module
 ! WDR no need use general_array_io

 implicit none

 
 
 private

 ! WDR public :: output_retrieval, get_modis_data_cube, write_statistics
 public :: get_modis_data_cube
 

 contains

!  subroutine write_statistics(filename)
!
!	use core_arrays
!
! include "hdf.f90"
! include "dffunc.f90"
! 
!
!	character(*), intent(in) :: filename
!
!	real :: data_to_write(17)
!	integer :: file_id, err_code, vref, vid
!	real, parameter :: fill_val = -999.
!	character(len=100) :: ln, un	
!	character(len=85) :: desc(19)
!	
!	data_to_write(1) = Statistics_1km%retrieval_fraction
!	data_to_write(2) = Statistics_1km%land_fraction 
!	data_to_write(3) = Statistics_1km%water_fraction 
!	data_to_write(4) = Statistics_1km%snow_fraction 
!	data_to_write(5) = Statistics_1km%cloud_fraction 
!	data_to_write(6) = Statistics_1km%water_cloud_fraction 
!  	data_to_write(7) = Statistics_1km%ice_cloud_fraction 
!	data_to_write(8) = Statistics_1km%mean_liquid_tau 
!	data_to_write(9) = Statistics_1km%mean_ice_tau
!	data_to_write(10) = Statistics_1km%mean_liquid_re21 
!	data_to_write(11) = Statistics_1km%mean_ice_re21
!	data_to_write(12) = Statistics_1km%ctp_liquid 
!	data_to_write(13) = Statistics_1km%ctp_ice
!	data_to_write(14) = Statistics_1km%ctp_undetermined
!	data_to_write(15) = Statistics_1km%ctt_liquid 
!	data_to_write(16) = Statistics_1km%ctt_ice
!	data_to_write(17) = Statistics_1km%ctt_undetermined
!	
!	file_id = hopen(filename, DFACC_WRITE, 0)
!	
!	err_code = vfstart(file_id)
!	
!	vref = vsffnd(file_id, "Statistics_1km")	
!	vid = vsfatch(file_id, vref, "w")
!
!	err_code = vsfwrt(vid, data_to_write, 17, 0)
!
!	err_code = vsfdtch(vid)
!	err_code = vfend(file_id)
!	err_code = hclose(file_id)
!
!! duplicate same information as an SDS
!! we can't use MOD_PR06CR because CR automatically turns any 1D SDS into VData
!! VData makes it difficult to see SDS attributes in common display tools
!! without the attributes the data is useless
!	file_id = sfstart(filename, DFACC_WRITE)
!	vid = sfselect(file_id, sfn2index(file_id, "Statistics_1km_sds"))
!	if (vid == -1) vid = sfcreate(file_id, "Statistics_1km_sds", DFNT_FLOAT, 1, (/ 17 /))
!	err_code = sfwdata(vid, (/ 0 /), (/ 1 /), (/ 17 /), data_to_write)
!
!	ln = "Granule Statistics for parameters at 1x1 resolution"
!	un = "see description attribute"
!	err_code = sfsattr(vid, "long_name", DFNT_CHAR, len(trim(ln)), ln )
!	err_code = sfsattr(vid, "units", DFNT_CHAR, len(trim(un)), un )
!	err_code = sfsattr(vid, "_FillValue", DFNT_FLOAT, 1, fill_val)
!
!	desc(1:19)(1:85) = ' '
!
!	desc(1) = " "
!	desc(2) = "Statistics_1km:"
!    desc(3) = "  1. Successful Retrieval Rate (%)"
!    desc(4) = "  2. Land Cover Fraction (%)"
!    desc(5) = "  3. Water Cover Fraction (%)"
!   	desc(6) = "  4. Snow Cover Fraction (%)"
!   	desc(7) = "  5. Cloud Cover Fraction (%)"
!   	desc(8) = "  6. Water Cloud Detected (%)"
!  	desc(9) = "  7. Ice Cloud Detected (%)"
!   	desc(10) = "  8. Mean of Water Cloud Optical Thickness"
!   	desc(11) = "  9. Mean of Ice Cloud Optical Thickness "
!   	desc(12) = "  10. Mean of Water Cloud Effective Particle Radius (microns)"
!   	desc(13) = "  11. Mean of Ice Cloud Effective Diameter (microns)"
!   	desc(14) = "  12. Mean Liquid Water Cloud Top Pressure (mb)"
!   	desc(15) = "  13. Mean Ice Cloud Top Pressure (mb)"
!   	desc(16) = "  14. Mean Undetermined Cloud Top Pressure (mb)"
!   	desc(17) = "  15. Mean Liquid Water Cloud Top Temperature (K)"
!   	desc(18) = "  16. Mean Ice Cloud Top Temperature (K) "
!   	desc(19) = "  17. Mean Undetermined Cloud Top Temperature (K)"
!	
!	desc(1:19)(85:85) = char(10)
!	
!	err_code = sfsattr(vid, "description", DFNT_CHAR, 19*85, desc)	
!
!	err_code = sfendacc(vid)
!	err_code = sfend(file_id)
!
!
!
!  end subroutine write_statistics



! subroutine output_retrieval(mapi_filedata, &
!                            filename,      &
!                            currentscanX, currentscanY, nscansX, nscansY,  &
!                            start, edge, stride, &
!                            status)
!   use GeneralAuxType
!   use core_arrays
!   use nonscience_parameters
!   use mod06_run_settings
!   use libraryarrays
!   use specific_other
!   implicit none
!
!   integer,      intent(in)          :: mapi_filedata(:), currentscanX, currentscanY, nscansX, nscansY
!   integer,      intent(in)          :: start(:), edge(:), stride(:)
!   character(*), intent (in)         :: filename
!
!   integer,      intent (out)        :: status
!
!   integer                              :: buffer_xsize, buffer_ysize, i,j, count
!   integer(integer_twobyte)             :: fillint_twobyte
!   integer(integer_twobyte),allocatable :: outputbuffer(:)
!   real(double)                         :: scale, add_offset 
!   integer :: localstart(3), localedge(3), localstride(3), xsize, ysize
!   integer(integer_onebyte),allocatable :: quality_assurance_1km(:,:,:)
!   real, allocatable :: retrieval_diff(:,:,:)
!	logical, parameter :: HDF4_OUTPUT = .true.    
!
!   integer*1, allocatable, dimension(:,:) :: cloud_phase_COP
!   
!   status = success
!   xsize = size(optical_thickness_final,1)
!   ysize = size(optical_thickness_final,2)
!
!	out_xsize = xsize
!	out_ysize = ysize
!
!   allocate(quality_assurance_1km(9, xsize,ysize ))
!
!   call set_quality_data(xsize, ysize)
!
!    call convert_binary_qa(quality_assurance_1km, status)
!   	call writeqaarray_toolkit(mapi_filedata, start, stride, edge, quality_assurance_1km,'Quality_Assurance_1km', status , HDF4_OUTPUT)
!
!   deallocate(quality_assurance_1km)
!
!
!   allocate(outputbuffer_twobyte(xsize,ysize))
!
!	call writefloatarray_toolkit(mapi_filedata, start, stride, edge, optical_thickness_final, 'Cloud_Optical_Thickness', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, optical_thickness_final_PCL, 'Cloud_Optical_Thickness_PCL', status, HDF4_OUTPUT)
!
!	call writefloatarray_toolkit(mapi_filedata, start, stride, edge, optical_thickness_1621_final, 'Cloud_Optical_Thickness_1621', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, optical_thickness_1621_final_PCL, 'Cloud_Optical_Thickness_1621_PCL', status, HDF4_OUTPUT)
!
!	call writefloatarray_toolkit(mapi_filedata, start, stride, edge, optical_thickness_16_final, 'Cloud_Optical_Thickness_16', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, optical_thickness_16_final_PCL, 'Cloud_Optical_Thickness_16_PCL', status, HDF4_OUTPUT)
!
!	call writefloatarray_toolkit(mapi_filedata, start, stride, edge, optical_thickness_37_final, 'Cloud_Optical_Thickness_37', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, optical_thickness_37_final_PCL, 'Cloud_Optical_Thickness_37_PCL', status, HDF4_OUTPUT)
!
!
!	call writefloatarray_toolkit(mapi_filedata, start, stride, edge, effective_radius_21_final, 'Cloud_Effective_Radius', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, effective_radius_21_final_PCL, 'Cloud_Effective_Radius_PCL', status, HDF4_OUTPUT)
!
!	call writefloatarray_toolkit(mapi_filedata, start, stride, edge, effective_radius_16_final, 'Cloud_Effective_Radius_16', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, effective_radius_16_final_PCL, 'Cloud_Effective_Radius_16_PCL', status, HDF4_OUTPUT)
!
!	call writefloatarray_toolkit(mapi_filedata, start, stride, edge, effective_radius_37_final, 'Cloud_Effective_Radius_37', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, effective_radius_37_final_PCL, 'Cloud_Effective_Radius_37_PCL', status, HDF4_OUTPUT)
!
!	call writefloatarray_toolkit(mapi_filedata, start, stride, edge, effective_radius_1621_final, 'Cloud_Effective_Radius_1621', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, effective_radius_1621_final_PCL, 'Cloud_Effective_Radius_1621_PCL', status, HDF4_OUTPUT)
!
!
!	call writefloatarray_toolkit(mapi_filedata, start, stride, edge, liquid_water_path, 'Cloud_Water_Path', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, liquid_water_path_PCL, 'Cloud_Water_Path_PCL', status, HDF4_OUTPUT)
!
!	call writefloatarray_toolkit(mapi_filedata, start, stride, edge, liquid_water_path_16, 'Cloud_Water_Path_16', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, liquid_water_path_16_PCL, 'Cloud_Water_Path_16_PCL', status, HDF4_OUTPUT)
!   
!	call writefloatarray_toolkit(mapi_filedata, start, stride, edge, liquid_water_path_37, 'Cloud_Water_Path_37', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, liquid_water_path_37_PCL, 'Cloud_Water_Path_37_PCL', status, HDF4_OUTPUT)
!
!	call writefloatarray_toolkit(mapi_filedata, start, stride, edge, liquid_water_path_1621, 'Cloud_Water_Path_1621', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, liquid_water_path_1621_PCL, 'Cloud_Water_Path_1621_PCL', status, HDF4_OUTPUT)
!   
!
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, optical_thickness_error, 'Cloud_Optical_Thickness_Uncertainty', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, optical_thickness_16_error, 'Cloud_Optical_Thickness_Uncertainty_16', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, optical_thickness_37_error, 'Cloud_Optical_Thickness_Uncertainty_37', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, optical_thickness_1621_error, 'Cloud_Optical_Thickness_Uncertainty_1621', status, HDF4_OUTPUT)
!
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, effective_radius_21_error, 'Cloud_Effective_Radius_Uncertainty', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, effective_radius_16_error, 'Cloud_Effective_Radius_Uncertainty_16', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, effective_radius_37_error, 'Cloud_Effective_Radius_Uncertainty_37', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, effective_radius_1621_error, 'Cloud_Effective_Radius_Uncertainty_1621', status, HDF4_OUTPUT)
!
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, liquid_water_path_error, 'Cloud_Water_Path_Uncertainty', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, liquid_water_path_16_error, 'Cloud_Water_Path_Uncertainty_16', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, liquid_water_path_37_error, 'Cloud_Water_Path_Uncertainty_37', status, HDF4_OUTPUT)
!   	call writeint2array_toolkit(mapi_filedata, start, stride, edge, liquid_water_path_1621_error, 'Cloud_Water_Path_Uncertainty_1621', status, HDF4_OUTPUT)
!
!
!	call writefloatarray_toolkit(mapi_filedata, start, stride, edge, precip_water_094, 'Above_Cloud_Water_Vapor_094', status, HDF4_OUTPUT)
!
!	call writefloatarray_toolkit(mapi_filedata, start, stride, edge, irw_temperature, 'IRW_Low_Cloud_Temperature_From_COP', status, HDF4_OUTPUT)
! 
!
!	if (allocated(tau_liquid)) then 
!
!  		call writeint2array_toolkit(mapi_filedata, start, stride, edge, tau_liquid, 'Cloud_Optical_Thickness_All_Liquid', status, HDF4_OUTPUT)
!  		call writeint2array_toolkit(mapi_filedata, start, stride, edge, tau_ice, 'Cloud_Optical_Thickness_All_Ice', status, HDF4_OUTPUT)
!
!  		call writeint2array_toolkit(mapi_filedata, start, stride, edge, re21_liquid, 'Cloud_Effective_Radius_All_Liquid', status, HDF4_OUTPUT)
!  		call writeint2array_toolkit(mapi_filedata, start, stride, edge, re21_ice, 'Cloud_Effective_Radius_All_Ice', status, HDF4_OUTPUT)
!
!	endif
!
!	deallocate(outputbuffer_twobyte)
!	
!
!   	call write_failed_array(mapi_filedata, start, stride, edge, failure_metric, 'Retrieval_Failure_Metric', status, HDF4_OUTPUT)
!   	call write_failed_array(mapi_filedata, start, stride, edge, failure_metric_16, 'Retrieval_Failure_Metric_16', status, HDF4_OUTPUT)
!   	call write_failed_array(mapi_filedata, start, stride, edge, failure_metric_37, 'Retrieval_Failure_Metric_37', status, HDF4_OUTPUT)
!   	call write_failed_array(mapi_filedata, start, stride, edge, failure_metric_1621, 'Retrieval_Failure_Metric_1621', status, HDF4_OUTPUT)
!
!    call writebytearray_toolkit(mapi_filedata, start, stride, edge, cloud_layer_flag, 'Cloud_Multi_Layer_Flag', status, HDF4_OUTPUT)
!
!
!	allocate(cloud_phase_COP(xsize,ysize))
!	cloud_phase_COP = processing_information%path_and_outcome
!	where(cloud_phase_COP > 4) 
!		cloud_phase_COP = cloud_phase_COP - 8
!	end where
!
!   	call writebytearray_toolkit(mapi_filedata, start, stride, edge, cloud_phase_COP, 'Cloud_Phase_Optical_Properties', status, HDF4_OUTPUT)
!
!
!
!	deallocate(cloud_phase_COP)
!
!   	call write3Dfloatarray(mapi_filedata, start, stride, edge, atm_corr_refl, 'Atm_Corr_Refl', status, HDF4_OUTPUT)
!
!
!	if (currentscanX == nscansX .and. currentscanY == nscansY) then 
!	
!		localstart = 0
!		localstride = 1
!		localedge(1) = number_wavelengths
!		localedge(2) = number_iceradii
!	
!		call write_float_array(mapi_filedata, "Extinction_Efficiency_Ice", &
!			localstart(1:2), localstride(1:2), localedge(1:2), extinction_ice, status)	
!		call write_float_array(mapi_filedata, "Single_Scatter_Albedo_Ice", &
!			localstart(1:2), localstride(1:2), localedge(1:2), singlescattering_ice, status)
!		call write_float_array(mapi_filedata, "Asymmetry_Parameter_Ice", &
!			localstart(1:2), localstride(1:2), localedge(1:2), asymmetry_ice, status)
!	
!		localedge(2) = number_waterradii
!
!		call write_float_array(mapi_filedata, "Extinction_Efficiency_Liq", &
!			localstart(1:2), localstride(1:2), localedge(1:2), extinction_water, status)	
!		call write_float_array(mapi_filedata, "Single_Scatter_Albedo_Liq", &
!			localstart(1:2), localstride(1:2), localedge(1:2), singlescattering_water, status)
!		call write_float_array(mapi_filedata, "Asymmetry_Parameter_Liq", &
!			localstart(1:2), localstride(1:2), localedge(1:2), asymmetry_water, status)
!	
!	
!	endif
!
! end subroutine output_retrieval

 
 
 subroutine get_modis_data_cube ( level1b_filedata,       &
                                 geolocation_filedata, &
                                 start, edge, stride, meas_start, meas_edge, scan_number, debug, status)

!  Core sds of interest in the MODIS level 1B file read to fill the following:
!   latitude (MOD03)
!   longitude (MOD03)
!   band_measurements
!   solar_zenith_angle (MOD03)
!   sensor_zenith_angle (MOD03)
!   solar_azimuth_angle (MOD03)
!   sensor_azimuth_angle (MOD03)


   use GeneralAuxType
   use core_arrays
! WDR out	use general_array_io
   use science_parameters
   use nonscience_parameters
   use mod06_run_settings
!WDR 	use modis_reader
!WDR new use
   use ch_xfr

   implicit none


   integer, dimension (2), intent (in)       :: start, edge, stride, meas_start, meas_edge
   integer, dimension(:), intent(in)         :: level1b_filedata, geolocation_filedata
   integer, intent(in)                       :: scan_number
  
   logical,                intent (in)       :: debug
   integer,                intent (out)      :: status

   integer                                   :: numberofbands ,checkvariable,  i, j, k

   integer                                   :: xdimension, meas_xdimension, ydimension
   real,dimension(:,:),allocatable           :: level1b_buffer, sza_temp
   integer(integer_onebyte),dimension(:,:), allocatable :: uncertainty_buffer
   logical                                   :: errorflag, useoffset
   real :: unc_spec, unc_scale
	integer :: unc_idx

   logical ::   Cal_type_is_refl


   status = success
   numberofbands = size(bands)

!  get level 1b data
 
   meas_xdimension = meas_edge(1)
   xdimension =   edge(1)
   ydimension =   edge(2)
 
!   allocate(level1b_buffer(meas_xdimension, ydimension))
!   allocate(uncertainty_buffer(meas_xdimension, ydimension))
!   allocate(sza_temp(meas_xdimension, ydimension))
   
   
   solar_constant_37 = 10.9295 ! average of Terra and Aqua values from table in Platnick and Fontenla (2008)
	
!WDR out with the reads
!   do i = 1, numberofbands
!     if(bands(i) < 20 .or. bands(i) == 26) then
!       Cal_type_is_refl = .true.
!     else
!       Cal_type_is_refl = .false.
!     endif
!
!     call read_L1B(level1b_filedata,     &
!                   bands(i), Cal_type_is_refl,       &
!                   meas_xdimension,ydimension,&
!				   meas_start(1), &
!                   level1b_buffer,       &
!                   uncertainty_buffer, unc_spec, unc_scale )
!     band_measurements(:,i,:) = level1b_buffer(:,:)
!
!	unc_idx = 0
!
!	if (bands(i) == 1) unc_idx = band_0065
!	if (bands(i) == 2) unc_idx = band_0086
!	if (bands(i) == 5) unc_idx = band_0124
!	if (bands(i) == 6) unc_idx = band_0163
!	if (bands(i) == 7) unc_idx = band_0213
!	if (bands(i) == channel_37um) unc_idx = band_0370 - 1
!	
!	if (unc_idx /= 0) then 
!	
!		band_uncertainty(:,unc_idx,:) = uncertainty_buffer(:,:)
!		spec_uncertain(unc_idx) = unc_spec
!		uncertain_sf(unc_idx) = unc_scale
!
!	endif
!
!   enddo

!WDR fill the arrays with the l2gen values
!print*, "WDR check shape of band_measurements ", shape(band_measurements)
!print*, "WDR check shape of c2_refl ", shape(c2_refl)
!print*, "WDR check shape of band_uncertainty ", shape(band_uncertainty)
!print*, "WDR check shape of c2_bnd_unc ", shape(c2_bnd_unc)
band_measurements = c2_refl
band_uncertainty = c2_bnd_unc
spec_uncertain = c2_spec_unc
uncertain_sf = c2_unc_sf

!   deallocate(level1b_buffer, uncertainty_buffer)

	no_valid_data = 0

!  get full resolution Latitude and Longitude arrays
!WDR out again
!	call read_float_array(geolocation_filedata, "Latitude", start, stride, edge, latitude, status)
!	call read_float_array(geolocation_filedata, "Longitude", start, stride, edge, longitude, status) 
!	
!	useoffset = .false.
!	call read_int_array(geolocation_filedata, "SensorZenith", start, stride, edge, useoffset, sensor_zenith_angle, status)
!	
!	call read_int_array(geolocation_filedata, "SensorAzimuth", start, stride, edge, useoffset, sensor_azimuth_angle, status) 
!							  
!							  
!	call read_int_array(geolocation_filedata, "SolarZenith", meas_start, stride, meas_edge, useoffset, sza_temp, status)
!
!	call read_int_array(geolocation_filedata, "SolarAzimuth", start, stride, edge, useoffset, solar_azimuth_angle, status) 


!WDR not needed here
!   do i = 1, numberofbands
!     if(bands(i) < 20 .or. bands(i) == 26) then
!        where (band_measurements(:,i,:) > 0.) 
!           band_measurements(:,i,:) = band_measurements(:,i,:) / cos(d2r*sza_temp)
!        end where
!     endif
!   enddo
   
! WDR insert angles`		
!	if (scan_number == 1) then 
!		solar_zenith_angle(:,:) = sza_temp(1:edge(1), :)
!	else
!		solar_zenith_angle(:,:) = sza_temp(2:edge(1)+1, :)
!	endif
  latitude = c2_lat
  longitude = c2_lon
  sensor_zenith_angle = c2_senz
  sensor_azimuth_angle = c2_sena
  solar_zenith_angle = c2_solz
  solar_azimuth_angle = c2_sola

!	deallocate(sza_temp)


!  calculate the relative azimuth
! WDR we already compute this and send it in so...
!   relative_azimuth_angle =  solar_azimuth_angle + 180. - sensor_azimuth_angle

!	do j = 1,  ydimension
!	do i = 1, xdimension
	 
!       if (relative_azimuth_angle(i,j) <= 0.) relative_azimuth_angle(i,j) = -relative_azimuth_angle(i,j) 
!       if (relative_azimuth_angle(i,j) > 180.) relative_azimuth_angle(i,j) = 360. - relative_azimuth_angle(i,j)  
!     enddo
!   enddo

!   relative_azimuth_angle = abs(relative_azimuth_angle)
!   where(relative_azimuth_angle > 180.) relative_azimuth_angle =360. - relative_azimuth_angle
 
!  WDR already prepared
  relative_azimuth_angle = c2_relaz 
! WDR I believe that just a (+) relaz is needed
  relative_azimuth_angle = abs(relative_azimuth_angle)
   
   max_rel_azimuth = 0.
   min_rel_azimuth = 999999.
   
   
   max_sensor_zenith = 0.
   min_sensor_zenith = 99999.
   
   max_solar_zenith = 0.
   min_solar_zenith = 99999.
   
   do j=1, ydimension
		do i=1, xdimension
		
		if (relative_azimuth_angle(i,j) < min_rel_azimuth .and. relative_azimuth_angle(i,j) >= 0.) &
			min_rel_azimuth = relative_azimuth_angle(i,j)
		if (relative_azimuth_angle(i,j) > max_rel_azimuth .and. relative_azimuth_angle(i,j) >= 0.) &
			max_rel_azimuth = relative_azimuth_angle(i,j)
						
		if (solar_zenith_angle(i,j) < min_solar_zenith .and. solar_zenith_angle(i,j) >= 0. ) &
			min_solar_zenith = solar_zenith_angle(i,j)
		if (solar_zenith_angle(i,j) > max_solar_zenith .and. solar_zenith_angle(i,j) >= 0.) &
			max_solar_zenith = solar_zenith_angle(i,j)

					

	end do
   end do


! the sensor zenith is constant along a column
	do j=1, ydimension
		
		do i=1, xdimension
			if (sensor_zenith_angle(i,j) < 0.) exit ! bad data, get a different line
			if (sensor_zenith_angle(i,j) < min_sensor_zenith .and. sensor_zenith_angle(i,j) >= 0.) &
				min_sensor_zenith = sensor_zenith_angle(i,j)
			if (sensor_zenith_angle(i,j) > max_sensor_zenith .and. sensor_zenith_angle(i,j) >= 0.) &
				max_sensor_zenith = sensor_zenith_angle(i,j)
		end do
		
!		if (i >= xdimension) then
!			exit ! we are done, a good line of data
!		endif

	end do
 !  in case of a partial or full fail
 if( ( min_sensor_zenith == 9999. ) .OR. ( max_sensor_zenith == 0. ) ) then
   min_sensor_zenith = 0.
   max_sensor_zenith = 89.9
 endif

 end subroutine get_modis_data_cube

 
 end module modis_io_module
