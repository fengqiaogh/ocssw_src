 subroutine mod_pr06od( statistics, errorlevel)

!-----------------------------------------------------------------------
!f90 mod_pr06od
!
!Description:
!
!  retrieve cloud optical and microphysical properties from MODIS radiation
!  measurements
!
!input parameters:
!
!output parameters:
!
!revision history:
!
! v.1 July 2001 Initial work mag gray@climate.gsfc.nasa.gov
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
! Nakajima, T. and M.D. King, 1990: determination of the
!    optical thickness and effective particle radius of
!    clouds from reflected solar radiation measurements
!    part i: theory. j. atmos. sci., 47, 1878-1893.
!
!design note:
!
!end
!----------------------------------------------------------------------


   	use GeneralAuxType
   	use core_arrays
   	use nonscience_parameters
   	use global_model_grids,only:grids_are_read  
   	use ancillary_module
   	use specific_other
   	use modis_io_module
   	use modis_frontend_module
   	use modis_science_module
   	use libraryarrays
   	use general_science_module, only: init_science_arrays
   	use MOD06AlbedoEcoModule
   	use names
   	use science_parameters, only: no_valid_data
! WDR	use shared_io_module
! WDR need transfer information
   use ch_xfr

#if VIIRS_INST
   	use mod06_run_settings, only: IFF_yes, NASA_L1B
#endif

#if AHI_INST
   	use mod06_run_settings
#endif

#if !AMS_INST & !ASTER_INST & !AVIRIS_INST & !RSP_INST & !EPIC_INST & !SSFR_INST
   	use retrieval_irw
   	use FASCODE_routines
#endif

#ifdef EMAS_INST   
   	use modis_cirrus_module
#endif   
   
   	implicit none

#if HAVE_MPI
!	include 'mpif.h'
#endif	

!	include "hdf.f90"
!	include "dffunc.f90"

   	integer, parameter :: MODFILLEN = 10
   	real, intent(out)         :: statistics(7)
   	integer, intent(out)      :: errorlevel

   	integer, dimension (2)    :: start, edge, stride
   	integer, dimension (2)    :: localstart, localedge, localstride
!WDR added
      integer, dimension (2)    :: l2_start
   	integer, dimension(2) :: meas_start, meas_edge
 
   	integer                   :: status, tilesize
   	integer                   :: debug_start_iteration = 0, debug_stop_iteration = 0

   	real                      :: threshold_solar_zenith,      &
                                threshold_sensor_zenith,     &
                                threshold_relative_azimuth 

   	integer                   :: l1b_filedata(set_number_of_bands),        &
                                cloudmask_filedata(1),  &
                                geolocation_filedata(1), &
                                mod06_filedata(2), &
                                angle_filedata(1)

   	logical                   :: debug

	integer :: start_iteration

	integer :: file_id, var_id, err_code, mystart(2), mystride(2), myedge(2)

!	integer*1, dimension(:,:), allocatable :: eco, snow, sfc
!	real, dimension(:,:), allocatable :: albedo

	integer :: start_time, end_time, crate, cmax, total_start, total_end
	integer :: st1, et1, cr1, cm1, anc_filedata(2)

	integer :: start_iterationX, start_iterationY
	integer :: granule_offset(2)
	
	integer ::  mpi_stepX, mpi_stepY, mpi_index
	integer :: mpi_rc, mpi_err 
!WDR indicies
   integer :: iwdr, jwdr
! WDR need def here
   integer, parameter :: fail= -1

#if HAVE_MPI
	call MPI_INIT(mpi_err)
   	if (mpi_err /= MPI_SUCCESS) then
    	  print *,'Error starting MPI program. Terminating.'
    	  call MPI_ABORT(MPI_COMM_WORLD, mpi_rc, mpi_err)
   	end if

   	call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, mpi_err)
   	call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_nprocs, mpi_err)
#else
	mpi_rank = 0
	mpi_nprocs = 1
#endif	


!	call system_clock(total_start, crate, cmax)
!	print*, "system speed: ", crate, " cycles/sec"


	IM_cloudy_count = 0
	IM_water_cloud_count = 0
	IM_ice_cloud_count = 0
	IM_undet_count = 0
	
	Statistics_1km%retrieval_fraction = 0.
	Statistics_1km%land_fraction = 0.
	Statistics_1km%water_fraction = 0.
	Statistics_1km%snow_fraction = 0.
	Statistics_1km%cloud_fraction = 0.
	Statistics_1km%water_cloud_fraction = 0.
	Statistics_1km%ice_cloud_fraction = 0.
	Statistics_1km%mean_liquid_tau = 0.
	Statistics_1km%mean_ice_tau = 0.
	Statistics_1km%mean_liquid_re21 = 0.
	Statistics_1km%mean_ice_re21 = 0.
	Statistics_1km%ctp_liquid = 0.
	Statistics_1km%ctp_ice = 0.
	Statistics_1km%ctt_ice = 0.
	Statistics_1km%ctt_liquid = 0.
	Statistics_1km%ctp_undetermined = 0.
	Statistics_1km%ctt_undetermined = 0.
	

   	status = success
   	errorlevel = 0
   	statistics(:) = .0

!  get runtime setup parameters, array sizes, lat/long-pixel limits
    call initialize_run (   start, edge, stride,       &
                              tilesize,                  &
                              threshold_solar_zenith,    &
                              threshold_sensor_zenith,   &
                              threshold_relative_azimuth,&
                              status )
    if (status /= success) then
    	call local_message_handler('Problem reported in initialize_run. Check previous message/s', status,'mod_pr06od')
    endif

#ifndef SIM
!  check that all files exist and are openable
!  if not there is a complete algorithm failure
   	if( status  == success) then
    	call check_datasources(status)
      
     	if (status /= success) then
       		call local_message_handler('Problem reported in check_datasources. Check previous message/s', status,'mod_pr06od')
     	endif
   	endif
  
  
!  open core hdf files for reading and check
! WDR not needed
!   	if (status == success) then
!     	call openclose_files ('open',             &
!                           l1b_filedata,       &
!                           cloudmask_filedata, &
!                           geolocation_filedata,&  
!                           mod06_filedata,     &
!                           status)
!     	if (status /= success) then
!     		call local_message_handler('Problem reported in openclose_files . Check previous message/s', status,'mod_pr06od')
!     	endif
!   	endif
   
   
!  check the size of the array/s to be processed 
!   	if (status == success) then
! WDR knock out l1b io
!    	call check_datasize(l1b_filedata, start, stride, edge, status)
!     	if (status /= success) then
!       		call local_message_handler('Problem reported in check_datasize . Check previous message/s', status,'mod_pr06od')
!     	endif
!   	endif

! WDR try to insert the outside 'edge'
stride = (/1,1/)
edge(1) = c2_npix
edge(2) = c2_nlin

#else	
	start = 0
	stride = 1
	edge(1) = 90
	edge(2) = 10

	! WDR to not call ftn IO mod06_filedata(1) = sfstart(Amod06_name, DFACC_WRITE)
	l1b_filedata(1:4) = (/ 45, 46, 47, 48 /)	
	
#endif   
   
!  if the total number of lines to read is greater than the 
!  tilesize (number of lines to process/ iteration) then
!  we need to iterate.  

!     read libraries in preparation for later interpolation
    if (status == success ) then
    	call readlibraries_base(debug,status)
        if (status /= success) then
        	call local_message_handler('Problem reported in readlibraries',status,'mod_pr06od')
        endif
    endif
	

#if  !AMS_INST & !ASTER_INST & !AVIRIS_INST & !RSP_INST & !EPIC_INST & !SSFR_INST
	call init_FASCODE
	call init_irw
#endif

! for inventory metadata
   	total_number_of_pixels = edge(2)*edge(1)

	mystart = 0
	mystride = 1
	myedge(1) = edge(1)
	myedge(2) = edge(2)

	
!	allocate(eco(edge(1), edge(2)), snow(edge(1), edge(2)), sfc(edge(1), edge(2)), albedo(edge(1), edge(2)))
	
   	number_of_iterationsX = ceiling(real(edge(1)) / real(tilesize)) 
   	if (number_of_iterationsX == 0) number_of_iterationsX = 1
#if SEV_PR06OD || AHI_PR06OD || EPIC_OD
! because of GEO angle space, it's not possible to retrieve a single data stripe
! we have to use square blocks in order to keep the memory from going out of control 
   	number_of_iterationsY = ceiling(real(edge(2)) / real(tilesize))
   	if (number_of_iterationsY == 0) number_of_iterationsY = 1
#else
   	number_of_iterationsY = 1
#endif

   	grids_are_read = .false.
   	NISE_is_read = .false.
   	! WDR move to mod_pr06od.f90 snow_stats_are_read = .false.

   	debugPRN = .false.

   	anc_filedata(1) = mod06_filedata(1)
   	anc_filedata(2) = mod06_filedata(2) ! This may or may not be present, specific_ancillary module will figure that out. 

	if (platform_name == "RSP" .or. platform_name == "SSFR" .or. platform_name == "EPIC" ) & 
		anc_filedata(1) = cloudmask_filedata(1)

#if AHI_INST
	granule_offset = (/AHI_start_i, AHI_start_j/)
#else
	granule_offset = 0
#endif

!	print*, "NUM_ITER: X, Y, total: ", number_of_iterationsX, number_of_iterationsY, &
!											number_of_iterationsX * number_of_iterationsY

!	print*, edge
	
	mpi_stepX = mpi_nprocs
	mpi_stepY = 1
	start_iterationX = mpi_rank + 1 !1
	start_iterationY = 1


! stride is ALWAYS 1
	localstride = stride

!	open(unit=MY_UNIT_LUN, file = MY_TEXT_FILE)

!	print*, "proc: ", mpi_rank, "domain: ", start_iterationX, number_of_iterationsX, mpi_stepX

   	do iterationX = start_iterationX, number_of_iterationsX, mpi_stepX
		do iterationY = start_iterationY, number_of_iterationsY


!      		print*, "ITERATION X Y: ", iterationX, iterationY, "proc:", mpi_rank

			call system_clock(start_time, crate, cmax)
		
      		localstart  = start
      		localedge   = edge
      		
	  		meas_start = start
	  		meas_edge = edge


! we read one line more than we need on each side so that CSR can work correctly. 
!WDR out with this as we use the entire line available
!      		if (edge(1) > tilesize) then 
!      			if (iterationX == 1 ) then
!         			localedge(1) = tilesize
!		 			meas_edge(1) = tilesize+1
!	      		else  
!         			localstart(1)  = start(1) + (iterationX-1)*tilesize
!		 			meas_start(1) = start(1) + (iterationX-1)*tilesize-1
!         			localedge(1)   = tilesize
!		 			meas_edge(1) = tilesize+2
!      			endif
!			endif
!! this works because just above we set the meas_start(1) to be one less. 
!      		if ( edge(1) - iterationX*tilesize <= 0 ) then 
!				localedge(1) = edge(1) - (iterationX-1)*tilesize
!				meas_edge(1) = edge(1) - (iterationX-1)*tilesize+1
!	  		endif
!
!#if SEV_PR06OD || AHI_PR06OD || EPIC_OD
!      		if (edge(2) > tilesize) then 
!      			if (iterationY == 1) then
!         			localedge(2) = tilesize
!      			else 
!         			localstart(2)  = start(2) + (iterationY-1)*tilesize
!		 			localedge(2)   = tilesize
!      			endif
!			endif
!
!      		if ( edge(2) - iterationY*tilesize < 0 ) then 
!				localedge(2) = edge(2) - (iterationY-1)*tilesize
!	  		endif
!	  
!#endif
!
!			localstart = localstart + granule_offset

! these instruments do not currently do spatial variability, so there is no need to read the extra lines
!#if MAS_OD || SEV_PR06OD || MODIS_HKM || AMS_OD || ASTER_OD || AVIRIS_OD || RSP_OD || AHI_PR06OD || SIM || EPIC_OD || SSFR_OD
!			meas_start = localstart
!			meas_edge = localedge
!#endif




#ifndef GEOS5
! if we're using GDAS, we need to know the model size beforehand, it's constant anyhow. 
			grid_xsize = 360
			grid_ysize = 181
#endif

!print*, __FILE__, __LINE__," WDR pre allocate_arrays"
!print*, "WDR localedge: ", localedge
!print*, "WDR meas_edge: ", meas_edge

!     Setup and allocate all data arrays
      		if (status == success ) then
        		call allocate_arrays (localedge, &
							   meas_edge, &
							   start_iterationX, &
							   start_iterationY, &
                               status )
         		if (status /= success) then
           			call local_message_handler('Failure detected before allocate_modisarrays, check previous messages', &
                                      status, &
                                      'mod_pr06od')
         		endif
      		endif
      
      		call init_qualitydata
!	  		print*, "allocate: " , status

      		! get data cube to be processed and ancillary data arrays
! WDR make the change to start
! orig WDR l2_start = (/ c2_st_samp - 1, c2_scan /)
l2_start = (/ c2_st_samp, c2_scan - c2_nlin / 2 /)
!print*, __FILE__,__LINE__, " WDR l2_start is: ", l2_start
! WDR debug to see in/out sizes
!print*, __FILE__, __LINE__, " WDR debug localstart, localedge, localstride, meas_start, meas_edge"
!print*, localstart, localedge, localstride, meas_start, meas_edge
	 	     if (status == success ) then
    	     	call get_modis_data_cube (l1b_filedata, &
    	     	           geolocation_filedata, &
! WDR                      localstart,  &
                           l2_start,  &
                            localedge,   &
                            localstride, &
! WDR	            			meas_start, &
							l2_start, &
							 meas_edge, &
                            iterationX,   &
                            debug,       &
                              status )

        	 	if (status /= success) then
           			call local_message_handler('Failure detected before get_modis_data_cube check previous messages',&
                              status, &
        		           'mod_pr06od')
         			endif
      			endif

!  print*,  __FILE__, "  ", __LINE__," WDR after get_modis_data_cube"
!  print*, "and .86, 2.1 um rad vals"
!!  WDR for report the center line only
!  do jwdr = c2_nlin / 2 + 1, c2_nlin / 2 + 1
!    do iwdr = 1, c2_npix
!      write(*,"(i4,i4,f10.4,f10.4)"),iwdr,jwdr,band_measurements(iwdr, 2, jwdr), band_measurements(iwdr, 5, jwdr)
!      enddo
!    enddo

!					print*, "read data"	, status
! if it's an empty iteration, then don't even bother wasting time doing anything else
					if (no_valid_data == 1) then 
						print*, "no valid data"		  
! WDR 1may19 don't think this is needed		
!		 				call init_science_arrays
! WDR no need to call any output
!         				call output_retrieval( mod06_filedata, &
!                               Amod06_name,     &
!							   iterationX, iterationY, &
!							   number_of_iterationsX, number_of_iterationsY, &
!                               localstart - granule_offset, localedge, localstride, &
!                               status) 
!                               
                        cycle 

					endif




#ifdef GEOS5

#ifdef MCARS
			geos5_istart = localstart(1)
			geos5_jstart = 0
			geos5_iend = localstart(1) + localedge(1)-1
			geos5_jend = localstart(2) + localedge(2)-1
#else
	! GEOS-5 is NC4, so it makes sense to read it in pieces
			call find_geos5_bounds(geos5_istart, geos5_iend, geos5_jstart,  geos5_jend, latitude, longitude)
#endif

			grid_xsize = geos5_iend - geos5_istart + 1
			grid_ysize = geos5_jend - geos5_jstart + 1
			print*, geos5_istart, geos5_iend, geos5_jstart,  geos5_jend, grid_xsize, grid_ysize 
			call allocate_model(start_iterationX, start_iterationY, grid_xsize, grid_ysize)
#endif


!     get cloud decision info
      		if (status == success .and. no_valid_data == 0) then
!  print*, __FILE__, __LINE__," WDR l2_start for cloud processing is: ", l2_start
!
! WDR 12oct18 Why I had the l2_start reduction below...
!  if (l2_start(2) <= 0) then
!    l2_start(2) = 1
!    endif
   ! WDR change localstart to l2_start
!print*, __FILE__, __LINE__, " WDR localstart, localstride, localedge:", &
!  localstart, localstride, localedge
!
! WDR call this conditionally to read the CM from the CM file or
! take it from the l2 data
           if ( cm_from_l2 .EQ. 0 ) then
		    !	call modis_cloudprocessing(cloudmask_filedata, mod06_filedata, &
      ! WDR A lot of ftn I/O in this - will remove for now.  Just look in 
      !  earlier versions of mod_pr06od.f90
      print*, __FILE__, __LINE__," WDR Cloud mask option not available"
      status = fail
          !                          iterationX, l2_start, localstride, localedge, debug, status)		
         !		if (status /= success) then
         !		  	call local_message_handler('Failure detected before cloudprocessing. Check previous messages/s', &
          !                            status, &
          !                            'mod_pr06od')
         !		endif
            else
             ! WDR fill the cloudmask structure from transfer values
             call make_cld_msk_str()
             ! use the re-constituted modis_cloudprocessing
             ! out for now call modis_cloudprocessing()
             end if
      		endif
      		
!      		print*, "cloud mask read: ", status
      		
!     get ancillary data, right now that means albedos
      		if (status == success .and. no_valid_data == 0 ) then 
          ! WDR change localstart to l2_start
   ! print*, __FILE__, __LINE__, " WDR l2_start, localstride, localedge, localstart:", l2_start, localstride, localedge, localstart
        		 call set_ancillary(Agdas_name, Agdas_name2, Aozone_name,  &
                	            Ancepice_name,              &
                	            Anise_name,                 &
								anc_filedata, &
                    	        Aecosystem_data_name,       &
                        	    Asnowicealbedo_data_name,   &
								Aemissivity_name, &
                        	    l2_start, localstride, localedge, &
                        	    debug,                     &
                        	    status )
         		if (status /= success) then
         			call local_message_handler('Problem reported in set_ancillary. Check previous message/s', &
                                      status, &
                                      'mod_pr06od')
         		endif
      		endif            
! WDR print the line of result-read ancillary
!print*, __FILE__,__LINE__, " WDR p, l, cld P, T, Hgt meth, sfc T, IR phase"
! for center line controlled by wdr_wr_ch_vars values
!  do jwdr = c2_nlin / 2 + 1, c2_nlin / 2 + 1
!    do iwdr = 1, c2_npix
!      write(*,"(i4,i4,f10.4,f10.4,i5,f10.4,i5)"),iwdr,jwdr,cloud_top_pressure(iwdr, jwdr), cloud_top_temperature(iwdr, jwdr), cloud_height_method(iwdr, jwdr), surface_temperature(iwdr, jwdr), cloud_phase_infrared(iwdr, jwdr)
!      enddo
!    enddo

!			print*, "ancillary set", status

	      	if (status == 0  .and. no_valid_data == 0 ) then
	        	call readlibraries_extra(debug,status)
								
         		if (status /= success) then
         			call local_message_handler('Problem reported in readlibraries',status,'mod_pr06od')
         		endif
      		endif

!     process modis data for cloud optical thickness and cloud top effective
!     particle radius
!			print*, "libraries read", status


!	Thin cirrus retrieval using water vapor absorption band(s)
! This code massively segfaults when running old MAS.
#ifdef RETRIEVE_138
#if EMAS_INST
			call science_module_cirrus_twochannel(mod06_filedata,localstart,localedge,localstride,number_of_iterationsX,l1b_filedata)
#else
! For future implementation of a single channel 1.38µm cirrus retrieval
!	call science_module_cirrus(mod06_filedata,localstart,localedge,localstride)
#endif
#endif

      		if (status == success .and. no_valid_data == 0 ) then
         		call scienceinterface( threshold_solar_zenith,      &
                                threshold_sensor_zenith,     &
                                threshold_relative_azimuth,  &
                                debug,                       &
                                status)

         		if (status /= success) then
           			call local_message_handler('Problem reported in scienceinterface',status,'mod_pr06od')
         		endif
      		endif
!			print*, "done retrieval", status

!   eco(localstart(1)+1:localstart(1) + tilesize , &
!             localstart(2)+1:localstart(2) + localedge(2) ) = eco_out(:,:)
!     send output retrieval and qa data to a hdf file
!	print*, status, no_valid_data
      		if (status == success .and. no_valid_data == 0 ) then
! WDR debug to see data sizes
!print*, "WDR debug look at the output chunk sizes"
!print*, "localstart - granule_offset, localedge, localstride"
!print*, localstart - granule_offset, localedge, localstride
!
!WDR write out the basic radius and optical depth values for the box
!if ( iterationX == 1 ) then
!  print*,  __FILE__, "  ", __LINE__," WDR FINAL TEST REGION VALUES"
!  print*, "Cloud_Effective_Radius = effective_radius_21_final"
!  print*, "and Cloud_Optical_Thickness = optical_thickness_final"
!  print*, "and .86, 2.1 um rad vals"
!!print*, "OFF FOR NOW"
!!  WDR for report the center line only
!  do jwdr = c2_nlin / 2 + 1, c2_nlin / 2 + 1
!    !do iwdr = 1, c2_npix
!    do iwdr = 1,40
!      write(*,"(i4,i4,f10.4,f12.6,f10.4,f10.4)"),iwdr,jwdr,effective_radius_21_final(iwdr, jwdr), optical_thickness_final(iwdr, jwdr), band_measurements(iwdr, 2, jwdr), band_measurements(iwdr, 5, jwdr)
!      enddo
!    enddo
!! WDR same pattern as above but write the other rad/ thk
!  print*, "effective_radius and optical_thickness for 16, 37, 1621"
!  do jwdr = c2_nlin / 2 + 1, c2_nlin / 2 + 1
!    !do iwdr = 1, c2_npix
!    do iwdr = 1,40
!      write(*,"(i4,i4,f10.4,f12.6,f10.4,f12.6f10.4,f12.6)"),iwdr,jwdr,effective_radius_16_final(iwdr, jwdr),optical_thickness_16_final(iwdr, jwdr),effective_radius_37_final(iwdr, jwdr),optical_thickness_37_final(iwdr, jwdr),effective_radius_1621_final(iwdr, jwdr),optical_thickness_1621_final(iwdr, jwdr)
!      enddo
!    enddo
!  endif
!
! WDR knock out the last data write
!        		call output_retrieval(  mod06_filedata, &
!                               Amod06_name,     &
!							   iterationX, iterationY, &
!							   number_of_iterationsX, number_of_iterationsY, &
!                               localstart - granule_offset, localedge, localstride, &
!                               status)  
!         		if (status /= success) then
!           			call local_message_handler('Problem reported in output_retrieval',status,'mod_pr06od')
!         		endif
      		endif
!			print*, "output_retrieval", status
	  		call system_clock(end_time, crate, cmax)
	  
!	  		print*, "iteration_time: ", end_time - start_time
#if MODIS_INST & !MODIS_HKM & !SIM 
	  		call aggregate_statistics_1km
#endif

	 	enddo 
   	enddo

!print*, "WDR commenting out the deallocate_cleanup for TEST"

!    if (status == success) then
!		call deallocate_cleanup(status )
!        if (status /= success) then
!        	call local_message_handler('Problem reported in deallocate_cleanup', status,'mod_pr06od')
!        endif
!    endif
   
#ifndef HAVE_MPI
! WDR 		call set_processing_attributes(mod06_filedata) 
#endif

#ifndef SIM
!   	if (status == success) then
!   
!    	call openclose_files ('close',            &
!                           l1b_filedata,       &
!                           cloudmask_filedata, &
!                           geolocation_filedata,&
!                           mod06_filedata,     &
!                           status)
!     	if (status /= success) then
!       		call local_message_handler('Problem reported in openclose_files', status,'mod_pr06od')
!     	endif
!   	endif
#else

!	status = sfend(mod06_filedata(1))
!	print*, status

#endif

! here we do the final assignment of the Inventory metadata statistics
#ifndef MODIS_HKM
#if  MODIS_INST & !SIM 

	if (total_number_of_pixels > 0) then 
	   statistics(2) = real(IM_cloudy_count) / real(total_number_of_pixels) * 100.
	   Statistics_1km%land_fraction = Statistics_1km%land_fraction / real(total_number_of_pixels) * 100.
	   Statistics_1km%water_fraction = Statistics_1km%water_fraction / real(total_number_of_pixels) * 100.
	   Statistics_1km%snow_fraction = Statistics_1km%snow_fraction / real(total_number_of_pixels) * 100.
	else
	   statistics(2) = -9999.
	   Statistics_1km%land_fraction = -9999.
	   Statistics_1km%water_fraction = -9999.
	   Statistics_1km%snow_fraction = -9999.
	endif

	if (IM_cloudy_count > 0) then 
	   statistics(1) = real(IM_successful_retrieval_count) / real(IM_cloudy_count) * 100.
	   statistics(3) = real(IM_water_cloud_count) / real(IM_cloudy_count) * 100.
	   statistics(4) = real(IM_ice_cloud_count) / real(IM_cloudy_count) * 100.
	else
	   statistics(1) = -9999.
	   statistics(3:4) = -9999.
	endif

   
    Statistics_1km%retrieval_fraction = statistics(1)
	Statistics_1km%cloud_fraction = statistics(2)
	Statistics_1km%water_cloud_fraction = statistics(3)	
	Statistics_1km%ice_cloud_fraction = statistics(4)

	if (IM_water_cloud_count > 0) then 
		Statistics_1km%mean_liquid_tau = Statistics_1km%mean_liquid_tau / real(IM_water_cloud_count) 
		Statistics_1km%mean_liquid_re21 = Statistics_1km%mean_liquid_re21 / real(IM_water_cloud_count) 
		Statistics_1km%ctp_liquid = Statistics_1km%ctp_liquid / real(IM_water_cloud_count) 
		Statistics_1km%ctt_liquid = Statistics_1km%ctt_liquid / real(IM_water_cloud_count) 
	else
		Statistics_1km%mean_liquid_tau = -9999.
		Statistics_1km%mean_liquid_re21 = -9999.
		Statistics_1km%ctp_liquid = -9999.
		Statistics_1km%ctt_liquid = -9999.
	endif

	if (IM_ice_cloud_count > 0) then 
		Statistics_1km%mean_ice_tau = Statistics_1km%mean_ice_tau / real(IM_ice_cloud_count) 
		Statistics_1km%mean_ice_re21 = Statistics_1km%mean_ice_re21 / real(IM_ice_cloud_count) 
		Statistics_1km%ctp_ice = Statistics_1km%ctp_ice / real(IM_ice_cloud_count) 
		Statistics_1km%ctt_ice = Statistics_1km%ctt_ice / real(IM_ice_cloud_count) 
	else
		Statistics_1km%mean_ice_tau = -9999.
		Statistics_1km%mean_ice_re21 = -9999.
		Statistics_1km%ctp_ice = -9999.
		Statistics_1km%ctt_ice = -9999.
	endif
   
	if (IM_undet_count > 0) then 
		Statistics_1km%ctp_undetermined = Statistics_1km%ctp_undetermined / real(IM_undet_count) 
		Statistics_1km%ctt_undetermined = Statistics_1km%ctt_undetermined / real(IM_undet_count) 
    else
    	Statistics_1km%ctp_undetermined = -9999.
    	Statistics_1km%ctt_undetermined = -9999.
    endif
    
! WDR no outputs
!	call write_statistics(Amod06_name)
#endif
   
#endif
   
!    close (MY_UNIT_LUN)

   
   	call system_clock(total_end, crate, cmax)
!	print*, "Total Execution Time is: ", total_end - total_start, " cycles"
!    print*, "total execution time is: ", (total_end - total_start) / ( crate*1.0 ) / 60.0, " minutes"
   
#if HAVE_MPI
   call MPI_FINALIZE(mpi_err)
#endif


 end subroutine mod_pr06od

 subroutine make_cld_msk_str()
! WDR this will set up the cloudmask structure
 use core_arrays
 use ch_xfr
!
!  I think cloudmask is already set up, so just transfer the stuff into 
!  the structure
!
!  the logicals
 cloudmask(:,:)%cloudmask_determined = .false.
 cloudmask(:,:)%confident_cloudy = .false.
 cloudmask(:,:)%probablyclear_66 = .false.
 cloudmask(:,:)%probablyclear_95 = .false.
 cloudmask(:,:)%probablyclear_99 = .false.
 cloudmask(:,:)%snowice_surface = .false.
 cloudmask(:,:)%water_surface = .false.
 cloudmask(:,:)%coastal_surface = .false.
 cloudmask(:,:)%desert_surface = .false.
 cloudmask(:,:)%land_surface = .false.
 !
 do j = 1, c2_nlin
  do i = 1, c2_npix
   if ( c2_cld_det( i, j ) == 1 ) cloudmask(i,j)%cloudmask_determined = .true.
   if ( c2_conf_cld( i, j ) == 1 ) cloudmask(i,j)%confident_cloudy = .true.
   if ( c2_clr_66( i, j ) == 1 ) cloudmask(i,j)%probablyclear_66 = .true.
   if ( c2_clr_95( i, j ) == 1 ) cloudmask(i,j)%probablyclear_95 = .true.
   if ( c2_clr_99( i, j ) == 1 ) cloudmask(i,j)%probablyclear_99 = .true.
   if ( c2_sno_sfc( i, j ) == 1 ) cloudmask(i,j)%snowice_surface = .true.
   if ( c2_wtr_sfc( i, j ) == 1 ) cloudmask(i,j)%water_surface = .true.
   if ( c2_coast_sfc( i, j ) == 1 ) cloudmask(i,j)%coastal_surface = .true.
   if ( c2_desert_sfc( i, j ) == 1 ) cloudmask(i,j)%desert_surface = .true.
   if ( c2_lnd_sfc( i, j ) == 1 ) cloudmask(i,j)%land_surface = .true.
  end do
 end do
!
!  the byte
 cloudmask(:,:)%night = c2_night
 cloudmask(:,:)%sunglint = c2_glint
 cloudmask(:,:)%ocean_no_snow = c2_ocean_no_snow
 cloudmask(:,:)%ocean_snow = c2_ocean_snow
 cloudmask(:,:)%land_no_snow = c2_lnd_no_snow
 cloudmask(:,:)%land_snow = c2_lnd_snow
 cloudmask(:,:)%test_high_138 = c2_tst_h_138
 cloudmask(:,:)%test_visiblereflectance = c2_tst_vis_refl
 cloudmask(:,:)%test_visnirratio = c2_tst_vis_ratio
 cloudmask(:,:)%visible_cloudtest_250m = c2_vis_cld_250
 cloudmask(:,:)%applied_highcloud138 = c2_appl_hcld_138
 cloudmask(:,:)%applied_visiblereflectance = c2_appl_vis_refl
 cloudmask(:,:)%applied_visnirratio = c2_appl_vis_nir_ratio
!
! and end
 end subroutine make_cld_msk_str
