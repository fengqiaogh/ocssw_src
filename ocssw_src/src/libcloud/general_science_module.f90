module general_science_module

	implicit none

contains

	subroutine set_drel(threshold_relative_azimuth, drel)
		
		use mod06_run_settings
		use science_parameters
	
		real, intent(in) :: threshold_relative_azimuth
		real, intent(inout) :: drel
		

		if (min_solar_zenith > 0.98 .or. min_sensor_zenith > 0.98) then 
			drel = 30.
		else if (min_solar_zenith > 0.85 .or. min_sensor_zenith > 0.85) then 
			drel = 10.
		else
			drel = threshold_relative_azimuth
		endif	

	end subroutine set_drel

	subroutine set_interp_controls(i,j, scattering_angle, cur_wind_speed, drel, &
							threshold_solar_zenith,      &
                            threshold_sensor_zenith,     &
							wind_speed_only, interp_SS, interp_MS )

			use mod06_run_settings				
			use science_parameters
			use core_arrays
							
			real, intent(in) :: scattering_angle, cur_wind_speed, drel,threshold_solar_zenith,      &
                            threshold_sensor_zenith

			integer, intent(in) :: i, j

			logical, intent(inout) :: wind_speed_only, interp_SS, interp_MS

			real :: dsol, dsen, dscat
			real :: diff_scat_angle, diff_scat_angle_ss, diff_solar_zenith, &
					diff_sensor_zenith, diff_relative_azimuth, diff_wind_speed
			

	
			dsol = threshold_solar_zenith
			dsen = threshold_sensor_zenith

			diff_scat_angle = abs(scattering_angle - lastinterp_scat_angle)
			diff_scat_angle_ss = abs(scattering_angle - lastinterp_scat_angle_ss)
		
			if (scattering_angle < 60. .or. (scattering_angle > 120. .and. scattering_angle <= 160.)) then 
				dscat = 2.0
				dsen = 1.0
			else  if (scattering_angle >= 60 .and. scattering_angle <= 120.) then 
				dscat = dscat3
			else ! more than 160 degrees
				dscat = 1.0
				dsen = 0.5 
			endif
			
		
			diff_solar_zenith = abs(solar_zenith_angle(i,j) - lastinterp_solar_zenith)
			diff_sensor_zenith= abs(sensor_zenith_angle(i,j)- lastinterp_sensor_zenith)
			diff_relative_azimuth= abs(relative_azimuth_angle(i,j)-lastinterp_relative_azimuth)
		 
		
			diff_wind_speed = abs(cur_wind_speed - lastinterp_wind_speed)
		 
		 
			if (COX_MUNK .and. diff_wind_speed > threshold_wind_speed .and. &
				.not. (diff_solar_zenith > dsol) .and. &
				.not. (diff_sensor_zenith > dsen) .and. &
				.not. (diff_relative_azimuth > drel) .and. &
				.not. (diff_scat_angle > dscat) .and. &
				.not. (COX_MUNK .neqv. last_COX_MUNK)) then 
					wind_speed_only = .true.
			else
					! WDR temp change set .true. to do every pt
               wind_speed_only = .false. 
					!wind_speed_only = .true. 
			endif

! interpolate the libraries
			! WDR temp change set .true. to do every pt 
         interp_MS = .false.
			!interp_MS = .true.

			if ( diff_solar_zenith     > dsol  .or. &
				diff_sensor_zenith    > dsen .or. &
				diff_relative_azimuth > drel .or. &
				diff_scat_angle > dscat .or. &
				(COX_MUNK .neqv. last_COX_MUNK) .or. &
				(COX_MUNK .and. diff_wind_speed > threshold_wind_speed)) then 
				
					interp_MS = .true.

				if (diff_scat_angle > dscat) &
					lastinterp_scat_angle = scattering_angle
			
				if (diff_solar_zenith > dsol) &
					lastinterp_solar_zenith = solar_zenith_angle(i,j)
				if (diff_sensor_zenith > dsen) &
					lastinterp_sensor_zenith = sensor_zenith_angle(i,j)
				if (diff_relative_azimuth > drel) &
					lastinterp_relative_azimuth = relative_azimuth_angle(i,j)
			
				if (COX_MUNK .and. diff_wind_speed > threshold_wind_speed) &
					lastinterp_wind_speed = cur_wind_speed

				if (COX_MUNK .neqv. last_COX_MUNK) &
					last_COX_MUNK = COX_MUNK


			endif

         ! WDR temp change set .true. to do every pt
			interp_SS = .false.
			!interp_SS = .true.

			if ((diff_scat_angle_ss > 0.1) .or. (COX_MUNK .neqv. last_COX_MUNK) .or. &
               (COX_MUNK .and. diff_wind_speed > threshold_wind_speed)) then				
               	interp_SS = .true.
				lastinterp_scat_angle_ss = scattering_angle
			end if
	
	
	end subroutine set_interp_controls


	subroutine set_water_path_answers(i,j, finalize_liq, finalize_ice)

	    use libraryarrays
   		use libraryinterpolates
		use core_arrays
		use nonscience_parameters
		use retrieval_prep_logic
      use ch_xfr, only: c2_sensor_id, OCI_ID, OCIS_ID

		integer, intent(in) :: i, j
		logical, dimension(:), intent(in) :: finalize_liq, finalize_ice

		if (allocated(liquid_water_path)) &
			liquid_water_path(i,j) = fillvalue_real
		if (allocated(liquid_water_path_1621)) &
			liquid_water_path_1621(i,j) = fillvalue_real
		
		if (allocated(liquid_water_path_16)) &
			liquid_water_path_16(i,j) = fillvalue_real
		if (allocated(liquid_water_path_37)) &
			liquid_water_path_37(i,j) = fillvalue_real

		if (allocated(liquid_water_path_22)) &
			liquid_water_path_22(i,j) = fillvalue_real

		if (finalize_liq(1)) then 
			call compute_water_path(optical_thickness_16_final(i,j), &
								effective_radius_16_final(i,j), &
								liquid_water_density, &
								water_radii, &
								extinction_water(1, :), &	
								liquid_water_path_16(i,j))
		endif
		if (finalize_liq(2)) then 
			call compute_water_path(optical_thickness_final(i,j), &
								effective_radius_21_final(i,j), &
								liquid_water_density, &
								water_radii, &
								extinction_water(1, :), &
								liquid_water_path(i,j))			
		endif
		
      if( ( c2_sensor_id == OCI_ID ) .OR. ( c2_sensor_id == OCIS_ID ) )then
        if (finalize_liq(3)) then
            call compute_water_path(optical_thickness_22_final(i,j), &
                        effective_radius_22_final(i,j), &
                        liquid_water_density, &
                        water_radii, &
                        extinction_water(1, :), &
                        liquid_water_path_22(i,j))
         endif
      else
	   	if (finalize_liq(3)) then 
		   	call compute_water_path(optical_thickness_37_final(i,j), &
								effective_radius_37_final(i,j), &
								liquid_water_density, &
								water_radii, &
								extinction_water(1, :), &
								liquid_water_path_37(i,j))
    		endif
      endif
		
		if (finalize_liq(4)) then 
			call compute_water_path(optical_thickness_1621_final(i,j), &
								effective_radius_1621_final(i,j), &
								liquid_water_density, &
								water_radii, &
								extinction_water(1, :), &
								liquid_water_path_1621(i,j))	
		endif

		if (finalize_ice(1)) then 
			call compute_water_path(optical_thickness_16_final(i,j), &
								effective_radius_16_final(i,j), &
								ice_water_density, &
								ice_radii, &
								extinction_ice(1, :), &
								liquid_water_path_16(i,j))			
		endif

		if (finalize_ice(2)) then 
			call compute_water_path(optical_thickness_final(i,j), &
								effective_radius_21_final(i,j), &
								ice_water_density, &
								ice_radii, &
								extinction_ice(1, :), &
								liquid_water_path(i,j))
			
		endif
		
      if( ( c2_sensor_id == OCI_ID ) .OR. ( c2_sensor_id == OCIS_ID ) )then
        if (finalize_ice(3)) then
           call compute_water_path(optical_thickness_22_final(i,j), &
                        effective_radius_22_final(i,j), &
                        ice_water_density, &
                        ice_radii, &
                        extinction_ice(1, :), &
                        liquid_water_path_22(i,j))
        endif
      else
		  if (finalize_ice(3)) then 
			  call compute_water_path(optical_thickness_37_final(i,j), &
								effective_radius_37_final(i,j), &
								ice_water_density, &
								ice_radii, &
								extinction_ice(1, :), &
								liquid_water_path_37(i,j))
		  endif
		endif
		
		if (finalize_ice(4)) then 
			call compute_water_path(optical_thickness_1621_final(i,j), &
								effective_radius_1621_final(i,j), &
								ice_water_density, &
								ice_radii, &
								extinction_ice(1, :), &
								liquid_water_path_1621(i,j))		
		endif
				
	end subroutine set_water_path_answers


	subroutine set_failure_answers(i,j, RSS_final, set_near)

		use core_arrays
		use nonscience_parameters
		use mod06_run_settings
      use ch_xfr, only: c2_sensor_id, OCI_ID, OCIS_ID

		integer, intent(in) :: i,j
		real, dimension(:), intent(in) :: RSS_final
		logical, intent(in), dimension(:) :: set_near


		if (set_near(1)) then 
			failure_metric_16(i,j)%tau = nint(optical_thickness_16_final(i,j)/retr_scale_factor)
			if (optical_thickness_16_final(i,j) < 0.) failure_metric_16(i,j)%tau = fillvalue_int2
			failure_metric_16(i,j)%re = nint(effective_radius_16_final(i,j)/retr_scale_factor)
			if (effective_radius_16_final(i,j) < 0.) failure_metric_16(i,j)%re = fillvalue_int2
			failure_metric_16(i,j)%RSS = nint(RSS_final(re16)*100./retr_scale_factor)
			optical_thickness_16_final(i,j) = fillvalue_real
			effective_radius_16_final(i,j) = fillvalue_real
			liquid_water_path_16(i,j) = fillvalue_real					
		endif

		if (set_near(2)) then 
			failure_metric(i,j)%tau = nint(optical_thickness_final(i,j)/retr_scale_factor)
			if (optical_thickness_final(i,j) < 0) failure_metric(i,j)%tau = fillvalue_int2
			failure_metric(i,j)%re = nint(effective_radius_21_final(i,j)/retr_scale_factor)
			if (effective_radius_21_final(i,j) < 0.) failure_metric(i,j)%re = fillvalue_int2
			failure_metric(i,j)%RSS = nint(RSS_final(re21)*100./retr_scale_factor)
			effective_radius_21_final(i,j) = fillvalue_real
			optical_thickness_final(i,j) = fillvalue_real
			liquid_water_path(i,j) = fillvalue_real					
		endif
      if( ( c2_sensor_id == OCI_ID ) .OR. ( c2_sensor_id == OCIS_ID ) ) then
        if (set_near(3)) then
           failure_metric_22(i,j)%tau = nint(optical_thickness_22_final(i,j)/retr_scale_factor)
           if (optical_thickness_22_final(i,j) < 0) failure_metric_22(i,j)%tau = fillvalue_int2
           failure_metric_22(i,j)%re = nint(effective_radius_22_final(i,j)/retr_scale_factor)
           if (effective_radius_22_final(i,j) < 0.) failure_metric_22(i,j)%re = fillvalue_int2
           failure_metric_22(i,j)%RSS = nint(RSS_final(re22)*100./retr_scale_factor)
           effective_radius_22_final(i,j) = fillvalue_real
           optical_thickness_22_final(i,j) = fillvalue_real
           liquid_water_path_22(i,j) = fillvalue_real
        endif
      else
  		  if (set_near(3)) then 
  			 failure_metric_37(i,j)%tau = nint(optical_thickness_37_final(i,j)/retr_scale_factor)
  			 if (optical_thickness_37_final(i,j) < 0) failure_metric_37(i,j)%tau = fillvalue_int2
  			 failure_metric_37(i,j)%re = nint(effective_radius_37_final(i,j)/retr_scale_factor)
  			 if (effective_radius_37_final(i,j) < 0.) failure_metric_37(i,j)%re = fillvalue_int2
  			 failure_metric_37(i,j)%RSS = nint(RSS_final(re37)*100./retr_scale_factor)
  			 optical_thickness_37_final(i,j) = fillvalue_real
  			 effective_radius_37_final(i,j) = fillvalue_real
  			 liquid_water_path_37(i,j) = fillvalue_real					
  		  endif
      endif

		if (set_near(4)) then 
			failure_metric_1621(i,j)%tau = nint(optical_thickness_1621_final(i,j)/retr_scale_factor)
			if (optical_thickness_1621_final(i,j) < 0) failure_metric_1621(i,j)%tau = fillvalue_int2
			failure_metric_1621(i,j)%re = nint(effective_radius_1621_final(i,j)/retr_scale_factor)
			if (effective_radius_1621_final(i,j) < 0.) failure_metric_1621(i,j)%re = fillvalue_int2
			failure_metric_1621(i,j)%RSS = nint(RSS_final(re1621)*100./retr_scale_factor)
			optical_thickness_1621_final(i,j) = fillvalue_real
			effective_radius_1621_final(i,j) = fillvalue_real
			liquid_water_path_1621(i,j) = fillvalue_real					
		endif

	end subroutine set_failure_answers

	subroutine init_science_arrays
      ! W. Robinson, 1 may 2019 - set this to preserve lines processed 
      ! from a previous call and place them in the proper location 
      ! for the current scan
		use GeneralAuxType
	  	use nonscience_parameters
	  	use core_arrays
      use ch_xfr

      integer     :: scan_sav = -999 ! saved scan from previous call
      integer     :: iln

      ! WDR OK, this sets up the lines to just transfer array values
      ! to their new location relative to the new scan which points to 
      ! the center of the 3 line array (c2_scan).  The transfer is done 
      ! from xfr_from to xfr_to for a total of xfr_num lines.
      ! Otherwise, initialize here for lines scn_loop_st to scn_loop_en
      ! and process in modis_science_module.f90
      if( ( ( c2_scan - scan_sav ) .ge. 3 ) .or. &
          ( ( c2_scan - scan_sav ) .le. -3 ) ) then
        scn_loop_st = 1
        scn_loop_en = 3  ! for the 3-line standard
        !scn_loop_en = 5  ! for the 5-line test
        xfr_num = 0
      else if( ( c2_scan - scan_sav ) .eq. 2 ) then ! cur scan is 2 ahead
        scn_loop_st = 2
        scn_loop_en = 3
        xfr_num = 1
        xfr_from = (/ 3, 0 /)
        xfr_to = (/ 1, 0 /)
      else if( ( c2_scan - scan_sav ) .eq. 1 ) then ! cur scan is 1 ahead
        ! for no mandatory re-compute of center line (fastest)
        scn_loop_st = 3
        scn_loop_en = 3
        xfr_num = 2
        xfr_from = (/ 2, 3 /)
        xfr_to = (/ 1, 2 /)
        ! for mandatory re-compute of center line (does a re-compute of 
        !                                          output line)
        !scn_loop_st = 2
        !scn_loop_en = 3
        !xfr_num = 1
        !xfr_from = (/ 2, 0 /)
        !xfr_to = (/ 1, 0 /)
        ! for no efficiency at all (Do all 3 lines again)
        !scn_loop_st = 1
        !scn_loop_en = 3
        !xfr_num = 0
        ! test to go back to 3-lin use - well, it works
        !scn_loop_st = 1
        ! WDR for 5-lin test scn_loop_en = 3
        !scn_loop_en = 3
        !scn_loop_en = 5
        !xfr_num = 0
      else if( ( c2_scan - scan_sav ) .eq. 0 ) then ! cur scan is unchanged
        scn_loop_st = 0
        scn_loop_en = 0
        xfr_num = 0
      else if( ( c2_scan - scan_sav ) .eq. -1 ) then ! cur scan is 1 behind
        ! for no mandatory re-compute of center line
        scn_loop_st = 1
        scn_loop_en = 1
        xfr_num = 2
        xfr_from = (/ 2, 1 /)
        xfr_to = (/ 3, 2 /)
        !for mandatory re-compute of center line
        !scn_loop_st = 1
        !scn_loop_en = 2
        !xfr_num = 1
        !xfr_from = (/ 2, 0 /)
        !xfr_to = (/ 3, 0 /)
      else if( ( c2_scan - scan_sav ) .eq. -2 ) then ! cur scan is 2 behind
        scn_loop_st = 1
        scn_loop_en = 2
        xfr_num = 1
        xfr_from = (/ 1, 0 /)
        xfr_to = (/ 3, 0 /)
      endif
!print*, __FILE__, __LINE__
!print*, "scan_sav, c2_scan, scn_loop_st, scn_loop_en ", scan_sav, c2_scan, scn_loop_st, scn_loop_en
!print*, "xfr_num, xfr_from, xfr_to ", xfr_num, xfr_from, xfr_to
      scan_sav = c2_scan
  ! WDR We'll have a transfer phase (to put the lines in the right place)
  ! and a init phase to put fill in the open lines of the array
  !
  ! Transfer phase
  if( xfr_num > 0 ) then
    do iln = 1, xfr_num
		if (allocated(optical_thickness_final)) then 
        effective_radius_21_final(:, xfr_to(iln)) = effective_radius_21_final_sav(:, xfr_from(iln))
        optical_thickness_final(:, xfr_to(iln)) = optical_thickness_final_sav(:, xfr_from(iln))
        liquid_water_path(:, xfr_to(iln)) = liquid_water_path_sav(:, xfr_from(iln))
        effective_radius_21_final_PCL(:, xfr_to(iln)) = effective_radius_21_final_PCL_sav(:, xfr_from(iln))
        optical_thickness_final_PCL(:, xfr_to(iln)) = optical_thickness_final_PCL_sav(:, xfr_from(iln))
        liquid_water_path_PCL(:, xfr_to(iln)) = liquid_water_path_PCL_sav(:, xfr_from(iln))
        effective_radius_21_error(:, xfr_to(iln)) = effective_radius_21_error_sav(:, xfr_from(iln))
        optical_thickness_error(:, xfr_to(iln)) = optical_thickness_error_sav(:, xfr_from(iln))
        liquid_water_path_error(:, xfr_to(iln)) = liquid_water_path_error_sav(:, xfr_from(iln))
        
        effective_radius_1621_final(:, xfr_to(iln)) = effective_radius_1621_final_sav(:, xfr_from(iln))
        optical_thickness_1621_final(:, xfr_to(iln)) = optical_thickness_1621_final_sav(:, xfr_from(iln))
        liquid_water_path_1621(:, xfr_to(iln)) = liquid_water_path_1621_sav(:, xfr_from(iln))
        effective_radius_1621_final_PCL(:, xfr_to(iln)) = effective_radius_1621_final_PCL_sav(:, xfr_from(iln))
        optical_thickness_1621_final_PCL(:, xfr_to(iln)) = optical_thickness_1621_final_PCL_sav(:, xfr_from(iln))
        liquid_water_path_1621_PCL(:, xfr_to(iln)) = liquid_water_path_1621_PCL_sav(:, xfr_from(iln))
        effective_radius_1621_error(:, xfr_to(iln)) = effective_radius_1621_error_sav(:, xfr_from(iln))
        optical_thickness_1621_error(:, xfr_to(iln)) = optical_thickness_1621_error_sav(:, xfr_from(iln))
        liquid_water_path_1621_error(:, xfr_to(iln)) = liquid_water_path_1621_error_sav(:, xfr_from(iln))
        
        failure_metric(:,xfr_to(iln))%tau = failure_metric_sav(:,xfr_from(iln))%tau
        failure_metric(:,xfr_to(iln))%re = failure_metric_sav(:,xfr_from(iln))%re
        failure_metric(:,xfr_to(iln))%RSS = failure_metric_sav(:,xfr_from(iln))%RSS
        		
        failure_metric_1621(:,xfr_to(iln))%tau = failure_metric_1621_sav(:,xfr_from(iln))%tau
        failure_metric_1621(:,xfr_to(iln))%re = failure_metric_1621_sav(:,xfr_from(iln))%re
        failure_metric_1621(:,xfr_to(iln))%RSS = failure_metric_1621_sav(:,xfr_from(iln))%RSS
        
        atm_corr_refl(:,:, xfr_to(iln)) = atm_corr_refl_sav(:,:, xfr_from(iln))
        precip_water_094(:, xfr_to(iln)) = precip_water_094_sav(:, xfr_from(iln))
		endif
		
      if (allocated(optical_thickness_22_final)) then
        effective_radius_22_final(:, xfr_to(iln)) = effective_radius_22_final_sav(:, xfr_from(iln))
        optical_thickness_22_final(:, xfr_to(iln)) = optical_thickness_22_final_sav(:, xfr_from(iln))
        liquid_water_path_22(:, xfr_to(iln)) = liquid_water_path_22_sav(:, xfr_from(iln))
        effective_radius_22_final_PCL(:, xfr_to(iln)) = effective_radius_22_final_PCL_sav(:, xfr_from(iln))
        optical_thickness_22_final_PCL(:, xfr_to(iln)) = optical_thickness_22_final_PCL_sav(:, xfr_from(iln))
        liquid_water_path_22_PCL(:, xfr_to(iln)) = liquid_water_path_22_PCL_sav(:, xfr_from(iln))
        effective_radius_22_error(:, xfr_to(iln)) = effective_radius_22_error_sav(:, xfr_from(iln))
        optical_thickness_22_error(:, xfr_to(iln)) = optical_thickness_22_error_sav(:, xfr_from(iln))
        liquid_water_path_22_error(:, xfr_to(iln)) = liquid_water_path_22_error_sav(:, xfr_from(iln))
        failure_metric_22(:,xfr_to(iln))%tau = failure_metric_22_sav(:,xfr_from(iln))%tau
        failure_metric_22(:,xfr_to(iln))%re = failure_metric_22_sav(:,xfr_from(iln))%re
        failure_metric_22(:,xfr_to(iln))%RSS = failure_metric_22_sav(:,xfr_from(iln))%RSS
      endif
		if (allocated(seviri_cloudphase)) seviri_cloudphase(:, xfr_to(iln)) = seviri_cloudphase(:, xfr_from(iln))

  		optical_thickness_16_final(:, xfr_to(iln)) = optical_thickness_16_final_sav(:, xfr_from(iln))
  		optical_thickness_37_final(:, xfr_to(iln)) = optical_thickness_37_final_sav(:, xfr_from(iln))
  		effective_radius_16_final(:, xfr_to(iln)) = effective_radius_16_final_sav(:, xfr_from(iln))
  		effective_radius_37_final(:, xfr_to(iln)) = effective_radius_37_final_sav(:, xfr_from(iln))
  		liquid_water_path_16(:, xfr_to(iln)) = liquid_water_path_16_sav(:, xfr_from(iln))
  		liquid_water_path_37(:, xfr_to(iln)) = liquid_water_path_37_sav(:, xfr_from(iln))

  		optical_thickness_16_final_PCL(:, xfr_to(iln)) = optical_thickness_16_final_PCL_sav(:, xfr_from(iln))
  		optical_thickness_37_final_PCL(:, xfr_to(iln)) = optical_thickness_37_final_PCL_sav(:, xfr_from(iln))
  		effective_radius_16_final_PCL(:, xfr_to(iln)) = effective_radius_16_final_PCL_sav(:, xfr_from(iln))
  		effective_radius_37_final_PCL(:, xfr_to(iln)) = effective_radius_37_final_PCL_sav(:, xfr_from(iln))
  		liquid_water_path_16_PCL(:, xfr_to(iln)) = liquid_water_path_16_PCL_sav(:, xfr_from(iln))
  		liquid_water_path_37_PCL(:, xfr_to(iln)) = liquid_water_path_37_PCL_sav(:, xfr_from(iln))

  		optical_thickness_16_error(:, xfr_to(iln)) = optical_thickness_16_error_sav(:, xfr_from(iln))
  		optical_thickness_37_error(:, xfr_to(iln)) = optical_thickness_37_error_sav(:, xfr_from(iln))
  		effective_radius_16_error(:, xfr_to(iln)) = effective_radius_16_error_sav(:, xfr_from(iln))
  		effective_radius_37_error(:, xfr_to(iln)) = effective_radius_37_error_sav(:, xfr_from(iln))
  		liquid_water_path_16_error(:, xfr_to(iln)) = liquid_water_path_16_error_sav(:, xfr_from(iln))
  		liquid_water_path_37_error(:, xfr_to(iln)) = liquid_water_path_37_error_sav(:, xfr_from(iln))
 		cloud_layer_flag(:, xfr_to(iln)) = cloud_layer_flag_sav(:, xfr_from(iln))
  		ml_test_flag(:, xfr_to(iln)) = ml_test_flag_sav(:, xfr_from(iln))
  		CSR_flag_array(:, xfr_to(iln)) = CSR_flag_array_sav(:, xfr_from(iln))
  		irw_temperature(:, xfr_to(iln)) = irw_temperature_sav(:, xfr_from(iln))
  
	
      failure_metric_16(:,xfr_to(iln))%tau = failure_metric_16_sav(:,xfr_from(iln))%tau
      failure_metric_16(:,xfr_to(iln))%re = failure_metric_16_sav(:,xfr_from(iln))%re
      failure_metric_16(:,xfr_to(iln))%RSS = failure_metric_16_sav(:,xfr_from(iln))%RSS

      failure_metric_37(:,xfr_to(iln))%tau = failure_metric_37_sav(:,xfr_from(iln))%tau
      failure_metric_37(:,xfr_to(iln))%re = failure_metric_37_sav(:,xfr_from(iln))%re
      failure_metric_37(:,xfr_to(iln))%RSS = failure_metric_37_sav(:,xfr_from(iln))%RSS

		if (allocated(tau_liquid)) then 
			tau_liquid(:, xfr_to(iln)) = tau_liquid_sav(:, xfr_from(iln))
			tau_ice(:, xfr_to(iln)) = tau_ice_sav(:, xfr_from(iln))
			re21_liquid(:, xfr_to(iln)) = re21_liquid_sav(:, xfr_from(iln))
			re21_ice(:, xfr_to(iln)) = re21_ice_sav(:, xfr_from(iln))
		endif
      ! WDR this may at least need transfer and was not in the orig code:
      cloudsummary(:, xfr_to(iln)) = cloudsummary_sav(:, xfr_from(iln))
      ! WDR 31jul19 add in transfer of processing_information
      processing_information(:, xfr_to(iln)) = &
        processing_information_sav(:, xfr_from(iln))
      !
      !  WDR The point, table reflectance diagnostic arrays
      prd_out_refl_loc_2100( :, xfr_to(iln), : ) = &
        prd_out_refl_loc_2100_sav( :, xfr_from(iln), : )
      prd_out_refl_loc_1600( :, xfr_to(iln), : ) = &
        prd_out_refl_loc_1600_sav( :, xfr_from(iln), : )
      prd_out_refl_loc_2200( :, xfr_to(iln), : ) = &
        prd_out_refl_loc_2200_sav( :, xfr_from(iln), : )
      prd_out_refl_loc_1621( :, xfr_to(iln), : ) = &
        prd_out_refl_loc_1621_sav( :, xfr_from(iln), : )
     enddo
   endif  ! Transfer phase
   !
 	! Init phase
   if( scn_loop_st .ne. 0 ) then
      if (allocated(optical_thickness_final)) then

      optical_thickness_final(:, scn_loop_st:scn_loop_en) = fillvalue_real
      optical_thickness_1621_final(:, scn_loop_st:scn_loop_en) = fillvalue_real
      effective_radius_21_final(:, scn_loop_st:scn_loop_en) = fillvalue_real
      effective_radius_1621_final(:, scn_loop_st:scn_loop_en) = fillvalue_real
      liquid_water_path(:, scn_loop_st:scn_loop_en) = fillvalue_real
      liquid_water_path_1621(:, scn_loop_st:scn_loop_en) = fillvalue_real
      optical_thickness_final_PCL(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      optical_thickness_1621_final_PCL(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      effective_radius_21_final_PCL(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      effective_radius_1621_final_PCL(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      liquid_water_path_PCL(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      liquid_water_path_1621_PCL(:, scn_loop_st:scn_loop_en) = fillvalue_int2

      optical_thickness_error(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      effective_radius_21_error(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      liquid_water_path_error(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      optical_thickness_1621_error(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      effective_radius_1621_error(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      liquid_water_path_1621_error(:, scn_loop_st:scn_loop_en) = fillvalue_int2

      failure_metric(:,scn_loop_st:scn_loop_en)%tau = fillvalue_int2
      failure_metric(:,scn_loop_st:scn_loop_en)%re = fillvalue_int2
      failure_metric(:,scn_loop_st:scn_loop_en)%RSS = fillvalue_int2

      failure_metric_1621(:,scn_loop_st:scn_loop_en)%tau = fillvalue_int2
      failure_metric_1621(:,scn_loop_st:scn_loop_en)%re = fillvalue_int2
      failure_metric_1621(:,scn_loop_st:scn_loop_en)%RSS = fillvalue_int2

      atm_corr_refl(:,:, scn_loop_st:scn_loop_en) = fillvalue_real
      precip_water_094(:, scn_loop_st:scn_loop_en) = fillvalue_real

      endif
     
      if( allocated( optical_thickness_22_final ) ) then
        effective_radius_22_final(:, scn_loop_st:scn_loop_en ) = fillvalue_real
        optical_thickness_22_final(:, scn_loop_st:scn_loop_en ) = fillvalue_real
        liquid_water_path_22(:, scn_loop_st:scn_loop_en ) = fillvalue_real
        effective_radius_22_final_PCL(:, scn_loop_st:scn_loop_en ) = fillvalue_int2
        optical_thickness_22_final_PCL(:, scn_loop_st:scn_loop_en ) = fillvalue_int2
        liquid_water_path_22_PCL(:, scn_loop_st:scn_loop_en ) = fillvalue_int2
        effective_radius_22_error(:, scn_loop_st:scn_loop_en ) = fillvalue_real
        optical_thickness_22_error(:, scn_loop_st:scn_loop_en ) = fillvalue_real
        liquid_water_path_22_error(:, scn_loop_st:scn_loop_en ) = fillvalue_real
        prd_out_refl_loc_2200( :, scn_loop_st:scn_loop_en, : ) = fillvalue_real
        failure_metric_22(:,scn_loop_st:scn_loop_en)%tau = fillvalue_int2
        failure_metric_22(:,scn_loop_st:scn_loop_en)%re = fillvalue_int2
        failure_metric_22(:,scn_loop_st:scn_loop_en)%RSS = fillvalue_int2
      endif

      if (allocated(seviri_cloudphase)) seviri_cloudphase(:, scn_loop_st:scn_loop_en) = 0

      optical_thickness_16_final(:, scn_loop_st:scn_loop_en) = fillvalue_real
      optical_thickness_37_final(:, scn_loop_st:scn_loop_en) = fillvalue_real
      effective_radius_16_final(:, scn_loop_st:scn_loop_en) = fillvalue_real
      effective_radius_37_final(:, scn_loop_st:scn_loop_en) = fillvalue_real
      liquid_water_path_16(:, scn_loop_st:scn_loop_en) = fillvalue_real
      liquid_water_path_37(:, scn_loop_st:scn_loop_en) = fillvalue_real

      optical_thickness_16_final_PCL(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      optical_thickness_37_final_PCL(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      effective_radius_16_final_PCL(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      effective_radius_37_final_PCL(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      liquid_water_path_16_PCL(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      liquid_water_path_37_PCL(:, scn_loop_st:scn_loop_en) = fillvalue_int2

      optical_thickness_16_error(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      optical_thickness_37_error(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      effective_radius_16_error(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      effective_radius_37_error(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      liquid_water_path_16_error(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      liquid_water_path_37_error(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      cloud_layer_flag(:, scn_loop_st:scn_loop_en) = 0
      ml_test_flag(:, scn_loop_st:scn_loop_en) = 0
      CSR_flag_array(:, scn_loop_st:scn_loop_en) = 0
      irw_temperature(:, scn_loop_st:scn_loop_en) = fillvalue_real


      failure_metric_16(:,scn_loop_st:scn_loop_en)%tau = fillvalue_int2
      failure_metric_16(:,scn_loop_st:scn_loop_en)%re = fillvalue_int2
      failure_metric_16(:,scn_loop_st:scn_loop_en)%RSS = fillvalue_int2

      failure_metric_37(:,scn_loop_st:scn_loop_en)%tau = fillvalue_int2
      failure_metric_37(:,scn_loop_st:scn_loop_en)%re = fillvalue_int2
      failure_metric_37(:,scn_loop_st:scn_loop_en)%RSS = fillvalue_int2

      if (allocated(tau_liquid)) then
         tau_liquid(:, scn_loop_st:scn_loop_en) = fillvalue_int2
         tau_ice(:, scn_loop_st:scn_loop_en) = fillvalue_int2
         re21_liquid(:, scn_loop_st:scn_loop_en) = fillvalue_int2
         re21_ice(:, scn_loop_st:scn_loop_en) = fillvalue_int2
      endif
      !  also the refl_loc initialize
      prd_out_refl_loc_2100( :, scn_loop_st:scn_loop_en, : ) = fillvalue_real
      prd_out_refl_loc_1600( :, scn_loop_st:scn_loop_en, : ) = fillvalue_real
      prd_out_refl_loc_1621( :, scn_loop_st:scn_loop_en, : ) = fillvalue_real
  endif	 ! Init phase
	end subroutine init_science_arrays

  subroutine capture_arrays
    ! W. Robinson, 3 May 2019 new routine to save the states of the 3
    ! line science arrays after their derivation in modis_science_module.f90
    ! for use in the next line (really 3 lines) of data processing
    use GeneralAuxType
      use nonscience_parameters
      use core_arrays
      use ch_xfr

    integer   :: adims(3), np, nl, nb
    integer   :: checkvariable, sav_alloc = 0
    !
    ! if required, allocate the save arrays
    if( sav_alloc .eq. 0 ) then
      sav_alloc = 1
      adims = shape( atm_corr_refl )
      nb = adims(1)
      np = adims(2)
      nl = adims(3)

      allocate(optical_thickness_final_sav(np,nl), stat = checkvariable)
      allocate(optical_thickness_1621_final_sav(np,nl), stat = checkvariable)
      allocate(effective_radius_21_final_sav(np,nl), stat = checkvariable)
      allocate(effective_radius_1621_final_sav(np,nl), stat = checkvariable)
      allocate(liquid_water_path_sav(np,nl), stat = checkvariable)
      allocate(liquid_water_path_1621_sav(np,nl), stat = checkvariable)
      allocate(optical_thickness_final_PCL_sav(np,nl), stat = checkvariable)
      allocate(optical_thickness_1621_final_PCL_sav(np,nl), stat = checkvariable)
      allocate(effective_radius_21_final_PCL_sav(np,nl), stat = checkvariable)
      allocate(effective_radius_1621_final_PCL_sav(np,nl), stat = checkvariable)
      allocate(liquid_water_path_PCL_sav(np,nl), stat = checkvariable)
      allocate(liquid_water_path_1621_PCL_sav(np,nl), stat = checkvariable)
      allocate(optical_thickness_error_sav(np,nl), stat = checkvariable)
      allocate(effective_radius_21_error_sav(np,nl), stat = checkvariable)
      allocate(liquid_water_path_error_sav(np,nl), stat = checkvariable)
      allocate(optical_thickness_1621_error_sav(np,nl), stat = checkvariable)
      allocate(effective_radius_1621_error_sav(np,nl), stat = checkvariable)
      allocate(liquid_water_path_1621_error_sav(np,nl), stat = checkvariable)
      allocate(failure_metric_sav(np,nl), stat = checkvariable)
      allocate(failure_metric_1621_sav(np,nl), stat = checkvariable)
      allocate(atm_corr_refl_sav(nb,np,nl), stat = checkvariable)
      allocate(precip_water_094_sav(np,nl), stat = checkvariable)
      allocate(optical_thickness_16_final_sav(np,nl), stat = checkvariable)
      allocate(optical_thickness_37_final_sav(np,nl), stat = checkvariable)
      allocate(effective_radius_16_final_sav(np,nl), stat = checkvariable)
      allocate(effective_radius_37_final_sav(np,nl), stat = checkvariable)
      allocate(liquid_water_path_16_sav(np,nl), stat = checkvariable)
      allocate(liquid_water_path_37_sav(np,nl), stat = checkvariable)
      allocate(optical_thickness_16_final_PCL_sav(np,nl), stat = checkvariable)
      allocate(optical_thickness_37_final_PCL_sav(np,nl), stat = checkvariable)
      allocate(effective_radius_16_final_PCL_sav(np,nl), stat = checkvariable)
      allocate(effective_radius_37_final_PCL_sav(np,nl), stat = checkvariable)
      allocate(liquid_water_path_16_PCL_sav(np,nl), stat = checkvariable)
      allocate(liquid_water_path_37_PCL_sav(np,nl), stat = checkvariable)
      allocate(optical_thickness_16_error_sav(np,nl), stat = checkvariable)
      allocate(optical_thickness_37_error_sav(np,nl), stat = checkvariable)
      allocate(effective_radius_16_error_sav(np,nl), stat = checkvariable)
      allocate(effective_radius_37_error_sav(np,nl), stat = checkvariable)
      allocate(liquid_water_path_16_error_sav(np,nl), stat = checkvariable)
      allocate(liquid_water_path_37_error_sav(np,nl), stat = checkvariable)
      allocate(cloud_layer_flag_sav(np,nl), stat = checkvariable)
      allocate(ml_test_flag_sav(np,nl), stat = checkvariable)
      allocate(CSR_flag_array_sav(np,nl), stat = checkvariable)
      allocate(irw_temperature_sav(np,nl), stat = checkvariable)
      allocate(failure_metric_16_sav(np,nl), stat = checkvariable)
      allocate(failure_metric_37_sav(np,nl), stat = checkvariable)
      if (allocated(tau_liquid)) then
        allocate(tau_liquid_sav(np,nl), stat = checkvariable)
        allocate(tau_ice_sav(np,nl), stat = checkvariable)
        allocate(re21_liquid_sav(np,nl), stat = checkvariable)
        allocate(re21_ice_sav(np,nl), stat = checkvariable)
      endif
      allocate(cloudsummary_sav(np,nl), stat = checkvariable)
      allocate(processing_information_sav(np,nl), stat = checkvariable)
      ! WDR add the point / table diagnostics
      allocate( prd_out_refl_loc_2100_sav(np,nl,10), stat = checkvariable)
      allocate( prd_out_refl_loc_1600_sav(np,nl,10), stat = checkvariable)
      allocate( prd_out_refl_loc_1621_sav(np,nl,10), stat = checkvariable)
      if( ( c2_sensor_id == OCI_ID ) .OR. ( c2_sensor_id == OCIS_ID ) ) then
        allocate(optical_thickness_22_final_sav(np,nl), stat = checkvariable)
        allocate(effective_radius_22_final_sav(np,nl), stat = checkvariable)
        allocate(liquid_water_path_22_sav(np,nl), stat = checkvariable)
        allocate(optical_thickness_22_final_PCL_sav(np,nl), &
          stat = checkvariable)
        allocate(effective_radius_22_final_PCL_sav(np,nl), stat = checkvariable)
        allocate(liquid_water_path_22_PCL_sav(np,nl), stat = checkvariable)
        allocate(optical_thickness_22_error_sav(np,nl), stat = checkvariable)
        allocate(effective_radius_22_error_sav(np,nl), stat = checkvariable)
        allocate(liquid_water_path_22_error_sav(np,nl), stat = checkvariable)
        allocate(failure_metric_22_sav(np,nl), stat = checkvariable) 
        allocate( prd_out_refl_loc_2200_sav(np,nl,10), stat = checkvariable)
      endif
    endif
  ! and copy the data to the save arrays
  optical_thickness_final_sav = optical_thickness_final
  optical_thickness_1621_final_sav = optical_thickness_1621_final
  effective_radius_21_final_sav = effective_radius_21_final
  effective_radius_1621_final_sav = effective_radius_1621_final
  liquid_water_path_sav = liquid_water_path
  liquid_water_path_1621_sav = liquid_water_path_1621
  optical_thickness_final_PCL_sav = optical_thickness_final_PCL
  optical_thickness_1621_final_PCL_sav = optical_thickness_1621_final_PCL
  effective_radius_21_final_PCL_sav = effective_radius_21_final_PCL
  effective_radius_1621_final_PCL_sav = effective_radius_1621_final_PCL
  liquid_water_path_PCL_sav = liquid_water_path_PCL
  liquid_water_path_1621_PCL_sav = liquid_water_path_1621_PCL
  optical_thickness_error_sav = optical_thickness_error
  effective_radius_21_error_sav = effective_radius_21_error
  liquid_water_path_error_sav = liquid_water_path_error
  optical_thickness_1621_error_sav = optical_thickness_1621_error
  effective_radius_1621_error_sav = effective_radius_1621_error
  liquid_water_path_1621_error_sav = liquid_water_path_1621_error
  failure_metric_sav = failure_metric
  failure_metric_1621_sav = failure_metric_1621
  failure_metric_22_sav = failure_metric_22
  atm_corr_refl_sav = atm_corr_refl
  precip_water_094_sav = precip_water_094
  optical_thickness_16_final_sav = optical_thickness_16_final
  optical_thickness_37_final_sav = optical_thickness_37_final
  effective_radius_16_final_sav = effective_radius_16_final
  effective_radius_37_final_sav = effective_radius_37_final
  liquid_water_path_16_sav = liquid_water_path_16
  liquid_water_path_37_sav = liquid_water_path_37
  optical_thickness_16_final_PCL_sav = optical_thickness_16_final_PCL
  optical_thickness_37_final_PCL_sav = optical_thickness_37_final_PCL
  effective_radius_16_final_PCL_sav = effective_radius_16_final_PCL
  effective_radius_37_final_PCL_sav = effective_radius_37_final_PCL
  liquid_water_path_16_PCL_sav = liquid_water_path_16_PCL
  liquid_water_path_37_PCL_sav = liquid_water_path_37_PCL
  optical_thickness_16_error_sav = optical_thickness_16_error
  optical_thickness_37_error_sav = optical_thickness_37_error
  effective_radius_16_error_sav = effective_radius_16_error
  effective_radius_37_error_sav = effective_radius_37_error
  liquid_water_path_16_error_sav = liquid_water_path_16_error
  liquid_water_path_37_error_sav = liquid_water_path_37_error
  cloud_layer_flag_sav = cloud_layer_flag
  ml_test_flag_sav = ml_test_flag
  CSR_flag_array_sav = CSR_flag_array
  irw_temperature_sav = irw_temperature
  failure_metric_16_sav = failure_metric_16
  failure_metric_37_sav = failure_metric_37

  if( allocated( effective_radius_22_final ) ) then
    effective_radius_22_final_sav = effective_radius_22_final
    optical_thickness_22_final_sav = optical_thickness_22_final
    liquid_water_path_22_sav = liquid_water_path_22
    effective_radius_22_final_PCL_sav = effective_radius_22_final
    optical_thickness_22_final_PCL_sav = optical_thickness_22_final_PCL
    liquid_water_path_22_PCL_sav = liquid_water_path_22_PCL
    effective_radius_22_error_sav = effective_radius_22_error
    optical_thickness_22_error_sav = optical_thickness_22_error
    liquid_water_path_22_error_sav = liquid_water_path_22_error
    prd_out_refl_loc_2200_sav = prd_out_refl_loc_2200
  endif

  if (allocated(tau_liquid)) then
    tau_liquid_sav = tau_liquid
    tau_ice_sav = tau_ice
    re21_liquid_sav = re21_liquid
    re21_ice_sav = re21_ice
  endif
  cloudsummary_sav = cloudsummary
  processing_information_sav = processing_information
  ! table point refl diagnostic
  prd_out_refl_loc_2100_sav = prd_out_refl_loc_2100
  prd_out_refl_loc_1600_sav = prd_out_refl_loc_1600
  prd_out_refl_loc_1621_sav = prd_out_refl_loc_1621

  end subroutine capture_arrays

	subroutine assign_retrieval_error(xpoint, ypoint)
		use GeneralAuxType
	    use nonscience_parameters
	    use core_arrays
	    integer,           intent(in)       :: xpoint,ypoint


		if (allocated(optical_thickness_final)) then 
		
	    optical_thickness_final(xpoint, ypoint)     = fillvalue_real
	    optical_thickness_1621_final(xpoint, ypoint) = fillvalue_real
	    effective_radius_21_final(xpoint, ypoint)   = fillvalue_real
	    effective_radius_1621_final(xpoint, ypoint) = fillvalue_real
	    liquid_water_path(xpoint, ypoint) = fillvalue_real
	    liquid_water_path_1621(xpoint, ypoint) = fillvalue_real
	    optical_thickness_error(xpoint, ypoint) = fillvalue_int2
	    effective_radius_21_error(xpoint, ypoint) = fillvalue_int2
		liquid_water_path_error(xpoint, ypoint) = fillvalue_int2
   		optical_thickness_1621_error(xpoint, ypoint) = fillvalue_int2
   		effective_radius_1621_error(xpoint, ypoint) = fillvalue_int2
   		liquid_water_path_1621_error(xpoint, ypoint) = fillvalue_int2
		
		endif
		
	    optical_thickness_16_final(xpoint, ypoint)     = fillvalue_real
	    optical_thickness_37_final(xpoint, ypoint)     = fillvalue_real
	    effective_radius_16_final(xpoint, ypoint)   = fillvalue_real
	    effective_radius_37_final(xpoint, ypoint)   = fillvalue_real
	    liquid_water_path_16(xpoint, ypoint) = fillvalue_real
	    liquid_water_path_37(xpoint, ypoint) = fillvalue_real
	    optical_thickness_16_error(xpoint, ypoint) = fillvalue_int2
	    optical_thickness_37_error(xpoint, ypoint) = fillvalue_int2
	    effective_radius_16_error(xpoint, ypoint) = fillvalue_int2
	    effective_radius_37_error(xpoint, ypoint) = fillvalue_int2
   		liquid_water_path_16_error(xpoint, ypoint) = fillvalue_int2
   		liquid_water_path_37_error(xpoint, ypoint) = fillvalue_int2
   		cloud_layer_flag(xpoint, ypoint) = 0
   		ml_test_flag(xpoint, ypoint) = 0
   		CSR_flag_array(xpoint, ypoint) = 0 

      if( allocated( effective_radius_22_final ) ) then
        effective_radius_22_final(xpoint, ypoint) = fillvalue_real
        optical_thickness_22_final(xpoint, ypoint) = fillvalue_real
        liquid_water_path_22(xpoint, ypoint) = fillvalue_real
        effective_radius_22_final(xpoint, ypoint) = fillvalue_real
        effective_radius_22_error(xpoint, ypoint) = fillvalue_int2
        optical_thickness_22_error(xpoint, ypoint) = fillvalue_int2
        liquid_water_path_22_error(xpoint, ypoint) = fillvalue_int2
      endif
   		
		if (allocated(tau_liquid)) then 
			tau_liquid(xpoint, ypoint) = fillvalue_int2
			tau_ice(xpoint, ypoint) = fillvalue_int2
			re21_liquid(xpoint, ypoint) = fillvalue_int2
			re21_ice(xpoint, ypoint) = fillvalue_int2
		endif


	end subroutine assign_retrieval_error

	subroutine split_PCL(xdim, ydim)
	
		use core_arrays
	    use nonscience_parameters
		
		integer, intent(in) :: xdim, ydim
		
		integer :: i,j
	
	    do j=1, ydim
    		do i=1, xdim
    	
    			if (CSR_flag_array(i,j) == 1 .or. CSR_flag_array(i,j) == 3) then 
    		

					if (allocated(optical_thickness_final)) then 
						if (optical_thickness_final(i,j) > 0.) &
							optical_thickness_final_PCL(i,j) = nint(optical_thickness_final(i,j) / retr_scale_factor)
						optical_thickness_final(i,j) = fillvalue_real

						if (optical_thickness_1621_final(i,j) > 0.) &
							optical_thickness_1621_final_PCL(i,j) = nint(optical_thickness_1621_final(i,j) / retr_scale_factor)
						optical_thickness_1621_final(i,j) = fillvalue_real

						if (effective_radius_21_final(i,j) > 0.) &
							effective_radius_21_final_PCL(i,j) = nint(effective_radius_21_final(i,j) / retr_scale_factor)
						effective_radius_21_final(i,j) = fillvalue_real

						if (effective_radius_1621_final(i,j) > 0.) &
							effective_radius_1621_final_PCL(i,j) = nint(effective_radius_1621_final(i,j) / retr_scale_factor)
						effective_radius_1621_final(i,j) = fillvalue_real
						
						if (liquid_water_path(i,j) > 0.) &
							liquid_water_path_PCL(i,j) = nint(liquid_water_path(i,j))
						liquid_water_path(i,j) = fillvalue_real

						if (liquid_water_path_1621(i,j) > 0.) &
							liquid_water_path_1621_PCL(i,j) = nint(liquid_water_path_1621(i,j))
						liquid_water_path_1621(i,j) = fillvalue_real
    				endif

               if (allocated(optical_thickness_22_final)) then
                 if( optical_thickness_22_final(i,j) > 0.) &
                   optical_thickness_22_final_PCL(i,j) = nint(optical_thickness_22_final(i,j) / retr_scale_factor)
                 optical_thickness_22_final(i,j) = fillvalue_real

                 if( effective_radius_22_final(i,j) > 0.) &
                   effective_radius_22_final_PCL(i,j) = nint(effective_radius_22_final(i,j) / retr_scale_factor)
                 effective_radius_22_final(i,j) = fillvalue_real

                if( liquid_water_path_22(i,j) > 0.) &
                   liquid_water_path_22_PCL(i,j) = nint(liquid_water_path_22(i,j))
                 liquid_water_path_22(i,j) = fillvalue_real

               endif

					if (optical_thickness_37_final(i,j) > 0.) &
						optical_thickness_37_final_PCL(i,j) = nint(optical_thickness_37_final(i,j) / retr_scale_factor)
					optical_thickness_37_final(i,j) = fillvalue_real


					if (optical_thickness_16_final(i,j) > 0.) &
						optical_thickness_16_final_PCL(i,j) = nint(optical_thickness_16_final(i,j) / retr_scale_factor)
					optical_thickness_16_final(i,j) = fillvalue_real

					if (effective_radius_16_final(i,j) > 0.) &
						effective_radius_16_final_PCL(i,j) = nint(effective_radius_16_final(i,j) / retr_scale_factor)
					effective_radius_16_final(i,j) = fillvalue_real


					if (effective_radius_37_final(i,j) > 0.) &
						effective_radius_37_final_PCL(i,j) = nint(effective_radius_37_final(i,j) / retr_scale_factor)
					effective_radius_37_final(i,j) = fillvalue_real


					if (liquid_water_path_16(i,j) > 0.) &
						liquid_water_path_16_PCL(i,j) = nint(liquid_water_path_16(i,j))
					liquid_water_path_16(i,j) = fillvalue_real

					if (liquid_water_path_37(i,j) > 0.) &
						liquid_water_path_37_PCL(i,j) = nint(liquid_water_path_37(i,j))
					liquid_water_path_37(i,j) = fillvalue_real

    	
    			endif
    	
    		end do
   		end do
	
	
	
	end subroutine split_PCL




end module general_science_module

