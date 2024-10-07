 module specific_other
 
 implicit none
 
 contains

  subroutine set_l1b_names(level1b_name)
  
   use names

   character(*), intent(inout) :: level1b_name(:)

      level1b_name(:) = "none"
 
	level1b_name(1) = Alevel1b_name(1)
      
  end subroutine set_l1b_names

  subroutine set_IRW_channel(IRW_channel)
	use mod06_run_settings
	
	integer, intent(inout) :: IRW_channel
	
	IRW_channel = 1
  
  end subroutine set_IRW_channel


! subroutine create_extra(file_id, modis_wid, modis_ht)
! 
! 	use core_arrays
! 	
! 	include "hdf.f90"
! 	include "dffunc.f90"
! 	
! 	integer, dimension(:), intent(in) :: file_id
! 	integer, intent(in) :: modis_wid, modis_ht
! 	
! 	integer :: var_id, err_code
! 	character(len=200) :: long_name, units
! 	integer*2 :: fillValue_integer2, range_integer2(2)
! 	real*8 :: add_offset, scale_factor
! 	
! 	fillValue_integer2 = -9999
! 	range_integer2 = (/0, 15000/)
! 	add_offset = 0.
! 	scale_factor = 0.01
! 	
! 	var_id = sfselect(file_id(1), sfn2index(file_id(1), "Cloud_Optical_Thickness_All_Liquid"))
! 	if (var_id == -1) then
!	 	var_id = sfcreate(file_id(1), "Cloud_Optical_Thickness_All_Liquid", DFNT_INT16, 2, (/modis_wid, modis_ht/))
!		long_name = "Optical Thickness from VNSWIR-2.1um retrieval, all liquid attempts"
!		err_code = sfsattr(var_id, "long_name", DFNT_CHAR, len(trim(long_name)), trim(long_name))
!		units = "none"
!		err_code = sfsattr(var_id, "units", DFNT_CHAR, len(trim(units)), trim(units))
!		err_code = sfsattr(var_id, "_FillValue", DFNT_INT16, 1, fillValue_integer2)
!		err_code = sfsattr(var_id, "valid_range", DFNT_INT16, 2, range_integer2)
!		err_code = sfsattr(var_id, "scale_factor", DFNT_DOUBLE, 1, scale_factor)
!		err_code = sfsattr(var_id, "add_offset", DFNT_DOUBLE, 1, add_offset)
!		err_code = sfendacc(var_id)
!	endif	 
!
! 	var_id = sfselect(file_id(1), sfn2index(file_id(1), "Cloud_Optical_Thickness_All_Ice"))
! 	if (var_id == -1) then
!	 	var_id = sfcreate(file_id(1), "Cloud_Optical_Thickness_All_Ice", DFNT_INT16, 2, (/modis_wid, modis_ht/))
!		long_name = "Optical Thickness from VNSWIR-2.1um retrieval, all ice attempts"
!		err_code = sfsattr(var_id, "long_name", DFNT_CHAR, len(trim(long_name)), trim(long_name))
!		units = "none"
!		err_code = sfsattr(var_id, "units", DFNT_CHAR, len(trim(units)), trim(units))
!		err_code = sfsattr(var_id, "_FillValue", DFNT_INT16, 1, fillValue_integer2)
!		err_code = sfsattr(var_id, "valid_range", DFNT_INT16, 2, range_integer2)
!		err_code = sfsattr(var_id, "scale_factor", DFNT_DOUBLE, 1, scale_factor)
!		err_code = sfsattr(var_id, "add_offset", DFNT_DOUBLE, 1, add_offset)
!		err_code = sfendacc(var_id)
!	endif	 
! 
!  	var_id = sfselect(file_id(1), sfn2index(file_id(1), "Cloud_Effective_Radius_All_Liquid"))
! 	if (var_id == -1) then
!	 	var_id = sfcreate(file_id(1), "Cloud_Effective_Radius_All_Liquid", DFNT_INT16, 2, (/modis_wid, modis_ht/))
!		long_name = "Effective Radius from VNSWIR-2.1um retrieval, all liquid attempts"
!		err_code = sfsattr(var_id, "long_name", DFNT_CHAR, len(trim(long_name)), trim(long_name))
!		units = "microns"
!		err_code = sfsattr(var_id, "units", DFNT_CHAR, len(trim(units)), trim(units))
!		err_code = sfsattr(var_id, "_FillValue", DFNT_INT16, 1, fillValue_integer2)
!		err_code = sfsattr(var_id, "valid_range", DFNT_INT16, 2, range_integer2)
!		err_code = sfsattr(var_id, "scale_factor", DFNT_DOUBLE, 1, scale_factor)
!		err_code = sfsattr(var_id, "add_offset", DFNT_DOUBLE, 1, add_offset)
!		err_code = sfendacc(var_id)
!	endif	 
!
! 	var_id = sfselect(file_id(1), sfn2index(file_id(1), "Cloud_Effective_Radius_All_Ice"))
! 	if (var_id == -1) then
!	 	var_id = sfcreate(file_id(1), "Cloud_Effective_Radius_All_Ice", DFNT_INT16, 2, (/modis_wid, modis_ht/))
!		long_name = "Effective Radius from VNSWIR-2.1um retrieval, all ice attempts"
!		err_code = sfsattr(var_id, "long_name", DFNT_CHAR, len(trim(long_name)), trim(long_name))
!		units = "microns"
!		err_code = sfsattr(var_id, "units", DFNT_CHAR, len(trim(units)), trim(units))
!		err_code = sfsattr(var_id, "_FillValue", DFNT_INT16, 1, fillValue_integer2)
!		err_code = sfsattr(var_id, "valid_range", DFNT_INT16, 2, range_integer2)
!		err_code = sfsattr(var_id, "scale_factor", DFNT_DOUBLE, 1, scale_factor)
!		err_code = sfsattr(var_id, "add_offset", DFNT_DOUBLE, 1, add_offset)
!		err_code = sfendacc(var_id)
!	endif	 
!
! 
! end subroutine create_extra

 subroutine set_esfc(os_flag_in, x, y, esfc, os_flag )
 
	use libraryinterpolates
	use science_parameters
	use core_arrays
 
	logical, intent(in) :: os_flag_in
	real, dimension(:), intent(inout) :: esfc
	logical, intent(inout) :: os_flag
	integer, intent(in):: x, y
	
	if (os_flag_in) then 
		os_flag = .true.
	else
		os_flag = .false.
	endif

	if (os_flag .and. COX_MUNK) then 
		esfc(1) = int_surface_emissivity_water(1,2,1)
		esfc(2) = int_surface_emissivity_water(1,1,1)
	else
		esfc(1) = surface_emissivity_land(x,y,2)
		esfc(2) = surface_emissivity_land(x,y,1)
	endif
 
 end subroutine set_esfc



! this subroutine is intentionally left blank
  subroutine set_cox_munk_albedo(albedo, lib_albedo)
 
	use mod06_run_settings
 
	real, dimension(:), intent(in) :: albedo
	real, dimension(:), intent(in) :: lib_albedo
	
 
 end subroutine set_cox_munk_albedo

  subroutine get_band_idx(idx16, idx21, idx37, channel_37, idx_11, idx_alb37, idx_alb16) 
	
	use mod06_run_settings
 
	integer, intent(inout) :: idx16, idx21, idx37, channel_37, idx_11, idx_alb37, idx_alb16
 
	idx16 = band_0163
	idx21 = band_0213
	idx37 = band_0370-1
	idx_11 = band_0370 !(it's 7 in MODIS)
	channel_37 = set_bands(band_0370)
	idx_alb37 = band_0370 - 1
	idx_alb16 = band_0163
 
 end subroutine get_band_idx

 ! this subroutine is intentionally left blank
 subroutine get_channels
	
	use mod06_run_settings
 
 end subroutine get_channels

 ! this subroutine is intentionally left blank
 subroutine allocate_extra(xdim, ydim)
	
	use core_arrays
	
	integer, intent(in) :: xdim, ydim
#if 0
 	allocate(tau_liquid(xdim, ydim))
 	allocate(re21_liquid(xdim, ydim))
 	allocate(tau_ice(xdim, ydim))
 	allocate(re21_ice(xdim, ydim))
#endif 
 
 end subroutine allocate_extra

 ! this subroutine is intentionally left blank
 subroutine deallocate_extra
	
	use core_arrays
#if 0
 	deallocate(tau_liquid, tau_ice, re21_liquid, re21_ice)
#endif 
 
 end subroutine deallocate_extra

 subroutine get_data_dims(filename, start, stride, edge)

	use mod06_run_settings
 
   integer, dimension(:), intent(inout) :: start, stride, edge
   character(len=*), intent(in) :: filename

   start = set_start
   stride = set_stride
   edge = set_edge


 end subroutine get_data_dims


 ! this subroutine is intentionally left blank
 subroutine set_process_time(file_id)

	integer, intent(in) :: file_id
 
 end subroutine set_process_time

 subroutine set_PH_desert(surface, R040, thres)
 
	logical, intent(in) :: surface
	real, intent(in) :: R040
	real, intent(inout) :: thres
	
	if (surface .and. R040 < 0.5) thres = 9999.
	

 end subroutine set_PH_desert

 logical function set_ice_ratio(ice_ratio)
 
	real, intent(in) :: ice_ratio
	
	if (ice_ratio < 1.3) then 
		set_ice_ratio = .true.
	else
		set_ice_ratio = .false.
	endif
 
 end function set_ice_ratio
 
 
! this subroutine is intentionally left blank
 subroutine set_albedo
 
 end subroutine set_albedo

! this subroutine is intentionally left blank
 subroutine set_processing_extra(file_id)
 
	integer, intent(in) :: file_id
	
 end subroutine set_processing_extra
 

! subroutine check_datasize(l1b_filedata, start, stride, edge, status)
!
!!   use GeneralAuxType
!   use nonscience_parameters
!
!   implicit none
!
!	include "hdf.f90"
!	include "dffunc.f90"
!
!   integer, dimension(:), intent(in)                   :: l1b_filedata
!   integer, dimension (2), intent(inout) :: start, edge, stride
!   integer, intent(out)                  :: status
!
!   integer      :: Scans_Per_Granule
!   character*40 :: att_N, dtype
!   integer      :: RTN
!
!   integer :: attr_id, file_id
!
!   status = success
!
!!  get number of scans
!   att_N = 'Number of Scans'
!   dtype = 'INTEGER*4'
!  
!	file_id = l1b_filedata(1)
!	attr_id = sffattr(file_id, att_N)
!	
!	RTN = sfrattr(file_id, attr_id, Scans_Per_Granule)
!	print*, Scans_Per_Granule
!  
!  
!   if (rtn /= 0) then
!     call MODIS_SMF_SETDYNAMICMSG(1, &
!                                  'MAPI function GMFIN for Number of Scans  failed', &
!                                  'check_datasize')
!     status = failure
!   endif
!
!   edge(2) = Scans_Per_Granule * 10
!
! end subroutine check_datasize
 
subroutine aggregate_statistics_1km

	use modis_sciencestructure
	use core_arrays

	implicit none
	
	integer :: i, j, wid, ht
	
	
	wid = size(optical_thickness_final, 1)
	ht = size(optical_thickness_final, 2)
	
	do j=1, ht
		do i=1, wid
		
			if (.not. cloudsummary(i,j)%ocean_surface) &
									Statistics_1km%land_fraction = Statistics_1km%land_fraction + 1
			if (cloudsummary(i,j)%snowice_surface) Statistics_1km%snow_fraction = Statistics_1km%snow_fraction + 1
			if (cloudsummary(i,j)%ocean_surface) Statistics_1km%water_fraction = Statistics_1km%water_fraction + 1
		
			
			if (cloudsummary(i,j)%watercloud) then
				if (optical_thickness_final(i,j) > 0.) &
					Statistics_1km%mean_liquid_tau = Statistics_1km%mean_liquid_tau + optical_thickness_final(i,j)
				if (effective_radius_21_final(i,j) > 0.) &
					Statistics_1km%mean_liquid_re21 = Statistics_1km%mean_liquid_re21 + effective_radius_21_final(i,j)
				if (cloud_top_pressure(i,j) > 0.) & 
					Statistics_1km%ctp_liquid = Statistics_1km%ctp_liquid + cloud_top_pressure(i,j)
				if (cloud_top_temperature(i,j) > 0.) &
					Statistics_1km%ctt_liquid = Statistics_1km%ctt_liquid + cloud_top_temperature(i,j)
			endif
			
			if (cloudsummary(i,j)%icecloud) then 
				if (optical_thickness_final(i,j) > 0.) &
					Statistics_1km%mean_ice_tau = Statistics_1km%mean_ice_tau + optical_thickness_final(i,j)
				if (effective_radius_21_final(i,j) > 0.) &
					Statistics_1km%mean_ice_re21 = Statistics_1km%mean_ice_re21 + effective_radius_21_final(i,j)
				if (cloud_top_pressure(i,j) > 0.) & 
					Statistics_1km%ctp_ice = Statistics_1km%ctp_ice + cloud_top_pressure(i,j)
				if (cloud_top_temperature(i,j) > 0.) &
					Statistics_1km%ctt_ice = Statistics_1km%ctt_ice + cloud_top_temperature(i,j)
			endif
			
			if (cloudsummary(i,j)%unknowncloud) then 
				if (cloud_top_pressure(i,j) > 0.) & 
					Statistics_1km%ctp_undetermined = Statistics_1km%ctp_undetermined + cloud_top_pressure(i,j)
				if (cloud_top_temperature(i,j) > 0.) &
					Statistics_1km%ctt_undetermined = Statistics_1km%ctt_undetermined + cloud_top_temperature(i,j)			
			endif
		
		end do
	end do
	

end subroutine aggregate_statistics_1km


!  subroutine openclose_files ( directive,           &
!                            l1b_filedata,        &
!                            cloudmask_filedata,  &
!                            geolocation_filedata,&
!                            mod06_filedata,      &
!                            status)
!
!  use nonscience_parameters
!  use mod06_run_settings
!  use names
!  use core_arrays, only: platform_name
!! WDR no need  use general_array_io, only : open_file, close_file
!  use ch_xfr, only : cm_from_l2
!  
!  
!  implicit none
!
!	include "hdf.f90"
!	include "dffunc.f90"
!
!  character(*), intent(in) :: directive
!  integer, dimension(:), intent(inout)    :: l1b_filedata, cloudmask_filedata,  &
!                               geolocation_filedata,mod06_filedata
!                               
!  integer, intent(out)      :: status
!
!  integer :: err_code, i, nbands
!
!  status = success 
!
!  nbands = 1
!
!  if (directive == 'open') then
!
!
!	l1b_filedata(:) = -1
!! WDR knock out the l1b opens
!!	do i=1, nbands
!!		if (trim(Alevel1b_name(i)) == "none") cycle 
!!		call open_file(Alevel1b_name(i), l1b_filedata(i), DFACC_READ)
!!  	end do
!!
!! WDR conditional use of work file cm data
!   if ( cm_from_l2 .EQ. 0 ) THEN
!! WDR out for l2gen xfr   	call open_file(Acloudmask_name, cloudmask_filedata(1), DFACC_READ)
!     end if
!! WDR don't open output (and input) file
!!	call open_file(Amod06_name, mod06_filedata(1), DFACC_WRITE)
!! WDR knock out the geo opens
!!	call open_file(Ageolocation_name, geolocation_filedata(1), DFACC_READ)
!
!    
!  else
!
!!	do i=1, nbands
!!		if (l1b_filedata(i) == -1) cycle
!!		call close_file(Alevel1b_name(i), l1b_filedata(i))
!!	end do
!	
!   if ( cm_from_l2 .EQ. 0 ) THEN
!! WDR out for l2gen xfr   	call close_file(Acloudmask_name, cloudmask_filedata(1))
!   endif
!! WDR	call close_file(Ageolocation_name, geolocation_filedata(1))
!! WDR	call close_file(Amod06_name, mod06_filedata(1))
!
!  endif
!
!
! end subroutine openclose_files



 subroutine convert_binary_qa( quality_assurance_1km, &
                             status)
 use core_arrays, only: processing_information, platform_name, cloudmask
 use nonscience_parameters
 use mod06_run_settings, only: set_bands, band_0213

 implicit none
 
 integer*1 , intent(out) ::  quality_assurance_1km(:,:,:)
 integer, intent(inout) :: status

 integer   :: i,j, blah, cm_wid, cm_ht
 
 quality_assurance_1km = 0

 cm_wid = size(processing_information, 1)
 cm_ht = size(processing_information, 2)


  do j= 1, cm_ht  
	do i = 1, cm_wid


!    Quality Assurance 1KM, byte 1 ----------------------------------------------------------------------------
		quality_assurance_1km(1,i,j) = processing_information(i,j)%optical_thickness_GC

			if (cloudmask(i,j)%ocean_no_snow == 1) then 
				! value of 00 for bits 4,3
			else if (cloudmask(i,j)%ocean_snow == 1) then 
	           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),3)
			else if (cloudmask(i,j)%land_no_snow == 1) then 
	           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),4)
			else if (cloudmask(i,j)%land_snow == 1) then 					
	           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),3)
	           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),4)
			endif

	 if (processing_information(i,j)%effective_radius_GC /= 0) then 
	 	quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),5)
	 	quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),6)
	 	quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),7)
	 endif



!    Quality Assurance 1KM, byte 2 ----------------------------------------------------------------------------

	 quality_assurance_1km(2,i,j) = ishft(processing_information(i,j)%path_and_outcome_1621, 3)

	 if (processing_information(i,j)%water_path_GC /= 0) then 
           quality_assurance_1km(2,i,j) =  ibset(quality_assurance_1km(2,i,j),0)
           quality_assurance_1km(2,i,j) =  ibset(quality_assurance_1km(2,i,j),1)
           quality_assurance_1km(2,i,j) =  ibset(quality_assurance_1km(2,i,j),2)
	 endif	 	

!    Quality Assurance 1KM, byte 3 ----------------------------------------------------------------------------

	  quality_assurance_1km(3,i,j) = processing_information(i,j)%path_and_outcome

     if(processing_information(i,j)%rayleigh_correction == 1) &
          quality_assurance_1km(3,i,j) = ibset(quality_assurance_1km(3,i,j),4) 

! atmospheric correction is always done. 
     quality_assurance_1km(3,i,j) = ibset(quality_assurance_1km(3,i,j),5) 

     if(processing_information(i,j)%band_used_for_optical_thickness ==1 ) then
          quality_assurance_1km(3,i,j) = ibset(quality_assurance_1km(3,i,j),6) 
     elseif(processing_information(i,j)%band_used_for_optical_thickness == 2 ) then
          quality_assurance_1km(3,i,j) = ibset(quality_assurance_1km(3,i,j),7)
     elseif(processing_information(i,j)%band_used_for_optical_thickness == 3 ) then
          quality_assurance_1km(3,i,j) = ibset(quality_assurance_1km(3,i,j),6)
          quality_assurance_1km(3,i,j) = ibset(quality_assurance_1km(3,i,j),7) 
     endif

!    Quality Assurance 1KM, byte 4 ----------------------------------------------------------------------------

	 quality_assurance_1km(4,i,j) = processing_information(i,j)%optical_thickness_1621_GC

	 if (processing_information(i,j)%effective_radius_1621_GC /= 0) then 
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),3)
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),4)
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),5)
	 endif
	 
     ! CSR QA added by G.Wind 4.7.05
     if (processing_information(i,j)%CSR_flag == 1) &
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),6)
     if (processing_information(i,j)%CSR_flag == 2) &
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),7)
     if (processing_information(i,j)%CSR_flag == 3) then
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),6)
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),7)
     endif
	

!    Quality Assurance 1KM, byte 5 ----------------------------------------------------------------------------

	 quality_assurance_1km(5,i,j) = processing_information(i,j)%water_path_1621_GC

		     if (processing_information(i,j)%multi_layer_cloud == 1) then
    		   quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),3)
    		 elseif(processing_information(i,j)%multi_layer_cloud == 2) then
    		   quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),4)
    		 elseif(processing_information(i,j)%multi_layer_cloud == 3) then
    		   quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),3)
      	 		quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),4)
     		elseif(processing_information(i,j)%multi_layer_cloud == 4) then
     		  quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),5)
    		 elseif(processing_information(i,j)%multi_layer_cloud == 5) then
       		quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),3)
      		 quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),5)
     		elseif(processing_information(i,j)%multi_layer_cloud == 6) then
       		quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),4)
      	 	quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),5)
     		elseif(processing_information(i,j)%multi_layer_cloud == 7) then
     		  quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),3)
     		  quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),4)
     		  quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),5)
     		endif 

    	 !This is a COPY of the retrieval_outcome bit for Level 3 compatibility
	     if(processing_information(i,j)%path_and_outcome > 4)  &
    	      quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),6)
	
	

!    Quality Assurance 1KM, byte 6 ----------------------------------------------------------------------------
		quality_assurance_1km(6,i,j) = processing_information(i,j)%ml_test_mark
!    Quality Assurance 1KM, byte 7 ----------------------------------------------------------------------------
		quality_assurance_1km(7,i,j) = processing_information(i,j)%path_and_outcome_16
		quality_assurance_1km(7,i,j) = ior(quality_assurance_1km(7,i,j), &
								ishft(processing_information(i,j)%path_and_outcome_16_PCL, 4))

!    Quality Assurance 1KM, byte 8 ----------------------------------------------------------------------------
		quality_assurance_1km(8,i,j) = processing_information(i,j)%path_and_outcome_37
		quality_assurance_1km(8,i,j) = ior(quality_assurance_1km(8,i,j), &
								ishft(processing_information(i,j)%path_and_outcome_37_PCL, 4))

!    Quality Assurance 1KM, byte 9 ----------------------------------------------------------------------------
		quality_assurance_1km(9,i,j) = processing_information(i,j)%path_and_outcome_1621_PCL
		quality_assurance_1km(9,i,j) = ior(quality_assurance_1km(9,i,j), &
								ishft(processing_information(i,j)%path_and_outcome_PCL, 4))


  enddo
 enddo   
  
 end subroutine convert_binary_qa
 

 
 
 end module specific_other
