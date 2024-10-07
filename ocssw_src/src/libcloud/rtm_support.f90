module rtm_support
!
!	Disclaimer: All this stuff is courtesy of UW-Madison. Author unknown.
!		There are no comments in the original code. 
!
! WDR  Note that for the transition to l2gen, we changed the 
!  model_info(...)%<val> to c2_model_info%<val> for l2gen anc values
!
!
! have to use MASTER channel numbers in modis_planck, modis_bright, but MAS channels in PFAAST

	use mod06_run_settings
	use names
	use science_parameters, only : model_levels, d2r 
	use planck_functions
   use ch_xfr, only : c2_sensor_id, OCI_ID, OCIS_ID
#if !CT_CODE
	use core_arrays, only : c2_model_info, platform_name
#else 
	use ct_core_arrays, only : c2_model_info, platform_name
#endif

!	use FASCODE_routines
	use pfaast	
#ifdef VIIRS_INST
	use pfaast_modis
#endif	
	
	implicit none

	private


	integer :: start_level

#ifdef CT_CODE
#if !MAS_INST & !EMAS_INST
	integer, parameter :: nchan = set_number_of_bands
#else
	integer, parameter :: nchan = set_number_of_bands-3
#endif
#else
	integer, parameter :: nchan = 2
#endif

	real :: rtm_trans_2way(nchan, model_levels)
	real :: rtm_trans_atm_clr(nchan, model_levels)
	real :: rtm_rad_atm_clr(nchan, model_levels)
	real :: rtm_cloud_prof(nchan, model_levels)
	
	real :: B_profile(nchan, model_levels)
	
#ifndef CT_CODE

	real :: rtm_trans_2way_mean(model_levels)
	
	real :: rtm_trans_2way_low(nchan, model_levels)
	real :: rtm_trans_atm_clr_low(nchan, model_levels)
	real :: rtm_rad_atm_clr_low(nchan, model_levels)
	real :: rtm_cloud_prof_low(nchan, model_levels)

	real :: rtm_trans_2way_high(nchan, model_levels)
	real :: rtm_trans_atm_clr_high(nchan, model_levels)
	real :: rtm_rad_atm_clr_high(nchan, model_levels)
	real :: rtm_cloud_prof_high(nchan, model_levels)
	

#endif	
	

	integer :: last_i, last_j
	real :: last_2way_angle
	real :: last_zenith
	
	

	public :: get_rtm_parameters, get_clear_toa_rad, rtm_cloud_prof, rtm_trans_2way, rtm_trans_atm_clr, rtm_rad_atm_clr
	public :: init_rtm_vars
#ifndef CT_CODE	
	public :: rtm_trans_2way_mean, rtm_rad_atm_clr_high, rtm_rad_atm_clr_low, &
				rtm_cloud_prof_high, rtm_cloud_prof_low
	public :: rtm_trans_2way_low, rtm_trans_atm_clr_low, rtm_trans_2way_high, rtm_trans_atm_clr_high
#endif

contains


  subroutine init_rtm_vars()
	
		last_zenith = -999.
		last_2way_angle = -999.
 
		last_i = -1
		last_j = -1
 
		start_level = 1
 
#if MAS_INST || EMAS_INST
		start_level = 36
#endif

  end subroutine init_rtm_vars

 subroutine get_rtm_parameters (platform, surface_emissivity, view_zenith, sun_zenith, i, j, x, y)
 
!  WDR - this looks to be entirely for 3.7 um and up.  Strictly speaking,
!  OCI should not need this, so try removing call to MODIS_fascode
!  True - OCI is unchanged, but Aqua does change - we'll have to do sensor
!  specific calling
	use mod06_run_settings
 
	real, intent(in) :: view_zenith, surface_emissivity(:), sun_zenith
	character(*), intent(in) :: platform
	integer, intent(in) :: i, j, x, y
	
	real :: ozone(model_levels)
	integer :: platform_id
	integer :: idet, iok, sfc_level
	real :: rtm_rad_clr, rtm_bt_clr
	integer :: channels(nchan)
	integer :: k, jj
	integer :: emis_idx, path_len, channel_use
	real :: angle_two_way, miu, miu0
	logical :: do_2way, newatm, newang, new_2way

	real :: temp_mixr(model_levels)
		
!	channels(1) = bands(band_0370)
#ifndef CT_CODE
	channels(1) = set_bands(band_1100)
	channels(2) = set_bands(band_0370)
#ifdef VIIRS_INST
	if (MODIS_MODE) then 
		channels(1) = set_bands_modis(band_1100)
		channels(2) = set_bands_modis(band_0370)
	endif
#endif	


	ozone = 0.

#else
   !  WDR, if I understand right, for modis, this is not used
	ozone = c2_model_info%o3_profile

#if !MAS_INST & !EMAS_INST
	channels = set_bands
#else
	channels = set_bands(1:set_number_of_bands-3)
#endif

#endif

	idet = 0

	
	! re-execute the transmittance calculations only if it makes sense to do so. 
	! There is more than one MODIS pixel that has the same vza and fits within same 1-degree box
	! so we will recycle the transmittance for points like that. 

	
	path_len = len(trim(ACT_lib_path))
	
		platform_id = -1
		if (platform(1:5) == "Terra" .or. platform(1:5) == "terra") then 
			platform_id = 0
		endif
		if (platform(1:4) == "Aqua" .or. platform(1:4) == "aqua" ) then 
			platform_id = 1
		endif
		
		miu0 = cos(d2r*sun_zenith)
		miu = cos(d2r*view_zenith)


		angle_two_way = acos ( miu*miu0 / (miu+miu0)) / d2r 
		
		newatm = .false.
		newang = .false.
		new_2way = .false.
		
		if (abs(view_zenith - last_zenith) > 0.005) newang = .true.
		if (angle_two_way /= last_2way_angle) new_2way = .true.
		if (last_i /= i .or. last_j /= j) newatm = .true.

		do_2way = .true.
		
#ifdef CT_CODE
		do_2way = .false.
		new_2way = .false. 
#endif
		if (newatm) then 
			do k=1, nchan
				do jj=start_level, model_levels
					B_profile(k,jj) = modis_planck(platform, c2_model_info%temp_profile(jj), channels(k), 1)
				end do
			end do
			
		end if
	
		do k=1, nchan
		

#ifdef MODIS_INST

         if( ( c2_sensor_id .ne. OCI_ID ) .and. &
           ( c2_sensor_id .ne. OCIS_ID ) ) then
			  call MODIS_fascode(trim(ACT_lib_path) // char(0), MYYEAR, &
             MYDAY, c2_model_info%temp_profile, c2_model_info%mixr_profile, &
             ozone, view_zenith, angle_two_way, platform_id, channels(k), &
             idet, rtm_trans_atm_clr(k,:), rtm_trans_2way(k,:), &
             newang, newatm, new_2way, do_2way, iok, x, y)			
           endif

#endif

#if MAS_INST || EMAS_INST
			if (platform == "MASTER") then 
#ifdef CT_CODE			
!				channel_use = set_bands_mas(k)
				channel_use = set_bands_master(k)
#else
				if (k==1) then 
!					channel_use = set_bands_mas(band_1100)
					channel_use = set_bands_master(k)
				else
					channel_use = channels(k)
				endif
#endif

			else
				channel_use = channels(k)
			endif
					
			call MAS_fascode(trim(ACT_lib_path), path_len, MYYEAR, MYDAY, c2_model_info%temp_profile(36:model_levels), &
					c2_model_info%mixr_profile(36:model_levels), ozone(36:model_levels), view_zenith, &
					angle_two_way, MYMISSION, channel_use, &
					 rtm_trans_atm_clr(k,36:model_levels), rtm_trans_2way(k, 36:model_levels), newang, newatm, new_2way, &
					do_2way, iok)
					
					
#endif

! need to get (and keep) information about which MSG satellite we have: 8 or 9
! this information must be passed to the code somehow by the production crew
#ifdef SEVIRI_INST	


			call SEVIRI_fascode(trim(ACT_lib_path), path_len, MYYEAR, MYDAY, c2_model_info%temp_profile, &
					c2_model_info%mixr_profile, ozone, view_zenith, angle_two_way, MYMSG, channels(k), &
					rtm_trans_atm_clr(k,:), rtm_trans_2way(k,:), newang, newatm, new_2way, do_2way, iok, x, y)

#endif

#ifdef AHI_INST	

		call AHI_fascode(trim(ACT_lib_path), path_len, MYYEAR, MYDAY, c2_model_info%temp_profile, &
					c2_model_info%mixr_profile, ozone, view_zenith, angle_two_way, channels(k), &
					rtm_trans_atm_clr(k,:), rtm_trans_2way(k,:), newang, newatm, new_2way, do_2way, iok, x, y)

#endif



#ifdef VIIRS_INST	
			if (MODIS_MODE) then 

				call MODIS_fascode(trim(ACT_lib_path) // char(0), MYYEAR, MYDAY, c2_model_info%temp_profile, &
					c2_model_info%mixr_profile, ozone, view_zenith, angle_two_way, &
					platform_id, channels(k), idet, rtm_trans_atm_clr(k,:), rtm_trans_2way(k,:), &
					newang, newatm, new_2way, do_2way, iok, x, y)			
			
			else
			

				call VIIRS_fascode(trim(ACT_lib_path), MYYEAR, MYDAY, c2_model_info%temp_profile, &
					c2_model_info%mixr_profile, ozone, view_zenith, angle_two_way, channels(k), &
					rtm_trans_atm_clr(k,:), rtm_trans_2way(k,:), newang, newatm, new_2way, do_2way, iok, x, y)
			endif		
#endif

		
		end do
		
#ifndef CT_CODE

#ifdef MODIS_INST
! we only need all this stuff for 3.7um channel, ever, so only run the 3.7um channel, don't waste time. 
		newatm = .true.
		temp_mixr = c2_model_info%mixr_profile*0.8
		do k=1, nchan
!		k = nchan
      if( ( c2_sensor_id .ne. OCI_ID ) .and. ( c2_sensor_id .ne. OCIS_ID ) ) then
		  call MODIS_fascode(trim(ACT_lib_path) // char(0), MYYEAR, MYDAY, &
          c2_model_info%temp_profile, temp_mixr, ozone, view_zenith, &
          angle_two_way, platform_id, channels(k), idet, &
          rtm_trans_atm_clr_low(k,:), rtm_trans_2way_low(k,:), &
			 newang, newatm, new_2way, do_2way, iok, x, y)					
        endif
		end do
		newatm = .true.
		temp_mixr = c2_model_info%mixr_profile*1.2
		do k=1, nchan
!		k = nchan
      if( ( c2_sensor_id .ne. OCI_ID ) .and. &
          ( c2_sensor_id .ne. OCIS_ID ) ) then
	  	  call MODIS_fascode(trim(ACT_lib_path) // char(0), MYYEAR, MYDAY, &
          c2_model_info%temp_profile, temp_mixr, ozone, view_zenith, &
          angle_two_way, platform_id, channels(k), idet, &
          rtm_trans_atm_clr_high(k,:), rtm_trans_2way_high(k,:), &
		    newang, newatm, new_2way, do_2way, iok, x, y)			
        endif
		end do

		rtm_trans_2way_mean(:) = ( abs(rtm_trans_2way_low(2,:) - rtm_trans_2way(2,:)) + &
											abs(rtm_trans_2way_high(2,:) - rtm_trans_2way(2,:)) ) / 2.

#endif

#if MAS_INST || EMAS_INST
! we only need all this stuff for 3.7um channel, ever, so only run the 3.7um channel, don't waste time. 
! the 3.7um channel is the same for MAS and MASTER
		newatm = .true.
		temp_mixr = c2_model_info%mixr_profile*0.8
		do k=1, nchan
!		k = nchan

			call MAS_fascode(trim(ACT_lib_path), path_len, MYYEAR, MYDAY, c2_model_info%temp_profile(36:model_levels), &
					temp_mixr(36:model_levels), ozone(36:model_levels), view_zenith, &
					angle_two_way, MYMISSION, channel_use, &
					 rtm_trans_atm_clr_low(k,36:model_levels), rtm_trans_2way_low(k, 36:model_levels), newang, newatm, new_2way, &
					do_2way, iok)
		end do

		newatm = .true.
		temp_mixr = c2_model_info%mixr_profile*1.2
		do k=1, nchan
!		k = nchan
			call MAS_fascode(trim(ACT_lib_path), path_len, MYYEAR, MYDAY, c2_model_info%temp_profile(36:model_levels), &
					temp_mixr(36:model_levels), ozone(36:model_levels), view_zenith, &
					angle_two_way, MYMISSION, channel_use, &
					 rtm_trans_atm_clr_high(k,36:model_levels), rtm_trans_2way_high(k, 36:model_levels), newang, newatm, new_2way, &
					do_2way, iok)
		end do

		rtm_trans_2way_mean(:) = ( abs(rtm_trans_2way_low(2,:) - rtm_trans_2way(2,:)) + &
											abs(rtm_trans_2way_high(2,:) - rtm_trans_2way(2,:)) ) / 2.

	
#endif



#ifdef SEVIRI_INST
! we only need all this stuff for 3.7um channel, ever, so only run the 3.7um channel, don't waste time. 
		newatm = .true.
		temp_mixr = c2_model_info%mixr_profile*0.8
		do k=1, nchan
!		k = nchan
			call SEVIRI_fascode(trim(ACT_lib_path), path_len, MYYEAR, MYDAY, c2_model_info%temp_profile, &
					temp_mixr, ozone, view_zenith, angle_two_way, MYMSG, channels(k), &
					rtm_trans_atm_clr_low(k,:), rtm_trans_2way_low(k,:), newang, newatm, new_2way, do_2way, iok, x, y)
		end do

		newatm = .true.
		temp_mixr = c2_model_info%mixr_profile*1.2
		do k=1, nchan
!		k = nchan
			call SEVIRI_fascode(trim(ACT_lib_path), path_len, MYYEAR, MYDAY, c2_model_info%temp_profile, &
					temp_mixr, ozone, view_zenith, angle_two_way, MYMSG, channels(k), &
					rtm_trans_atm_clr_high(k,:), rtm_trans_2way_high(k,:), newang, newatm, new_2way, do_2way, iok, x, y)
		end do
		
		rtm_trans_2way_mean(:) = ( abs(rtm_trans_2way_low(2,:) - rtm_trans_2way(2,:)) + &
											abs(rtm_trans_2way_high(2,:) - rtm_trans_2way(2,:)) ) / 2.


#endif

#ifdef AHI_INST
! we only need all this stuff for 3.7um channel, ever, so only run the 3.7um channel, don't waste time. 
		newatm = .true.
		temp_mixr = c2_model_info%mixr_profile*0.8
		do k=1, nchan
!		k = nchan
			call AHI_fascode(trim(ACT_lib_path), path_len, MYYEAR, MYDAY, c2_model_info%temp_profile, &
					temp_mixr, ozone, view_zenith, angle_two_way, channels(k), &
					rtm_trans_atm_clr_low(k,:), rtm_trans_2way_low(k,:), newang, newatm, new_2way, do_2way, iok, x, y)

		end do

		newatm = .true.
		temp_mixr = c2_model_info%mixr_profile*1.2
		do k=1, nchan
!		k = nchan

			call AHI_fascode(trim(ACT_lib_path), path_len, MYYEAR, MYDAY, c2_model_info%temp_profile, &
					temp_mixr, ozone, view_zenith, angle_two_way, channels(k), &
					rtm_trans_atm_clr_high(k,:), rtm_trans_2way_high(k,:), newang, newatm, new_2way, do_2way, iok, x, y)

		end do
		
		rtm_trans_2way_mean(:) = ( abs(rtm_trans_2way_low(2,:) - rtm_trans_2way(2,:)) + &
											abs(rtm_trans_2way_high(2,:) - rtm_trans_2way(2,:)) ) / 2.


#endif


#ifdef VIIRS_INST
! we only need all this stuff for 3.7um channel, ever, so only run the 3.7um channel, don't waste time. 

		newatm = .true.
		temp_mixr = c2_model_info%mixr_profile*0.8
		do k=1, nchan
!		k = nchan

		if (MODIS_MODE) then 
			call MODIS_fascode(trim(ACT_lib_path) // char(0), MYYEAR, MYDAY, c2_model_info%temp_profile, &
					temp_mixr, ozone, view_zenith, angle_two_way, &
					platform_id, channels(k), idet, rtm_trans_atm_clr_low(k,:), rtm_trans_2way_low(k,:), &
					newang, newatm, new_2way, do_2way, iok, x, y)					
		else
			call VIIRS_fascode(trim(ACT_lib_path), MYYEAR, MYDAY, c2_model_info%temp_profile, &
					temp_mixr, ozone, view_zenith, angle_two_way, channels(k), &
					rtm_trans_atm_clr_low(k,:), rtm_trans_2way_low(k,:), newang, newatm, new_2way, do_2way, iok, x, y)
		endif
		end do
		newatm = .true.
		temp_mixr = c2_model_info%mixr_profile*1.2
		do k=1, nchan
!		k = nchan

		if (MODIS_MODE) then 
			call MODIS_fascode(trim(ACT_lib_path) // char(0), MYYEAR, MYDAY, c2_model_info%temp_profile, &
					temp_mixr, ozone, view_zenith, angle_two_way, &
					platform_id, channels(k), idet, rtm_trans_atm_clr_high(k,:), rtm_trans_2way_high(k,:), &
					newang, newatm, new_2way, do_2way, iok, x, y)					
		else
			call VIIRS_fascode(trim(ACT_lib_path), MYYEAR, MYDAY, c2_model_info%temp_profile, &
					temp_mixr, ozone, view_zenith, angle_two_way, channels(k), &
					rtm_trans_atm_clr_high(k,:), rtm_trans_2way_high(k,:), newang, newatm, new_2way, do_2way, iok, x, y)
		endif
		end do

		rtm_trans_2way_mean(:) = ( abs(rtm_trans_2way_low(2,:) - rtm_trans_2way(2,:)) + &
											abs(rtm_trans_2way_high(2,:) - rtm_trans_2way(2,:)) ) / 2.
#endif



#endif					
					
				
	sfc_level = c2_model_info%surface_level

	last_zenith = view_zenith
	last_2way_angle = angle_two_way
	last_i = i
	last_j = j
		
	do k=1, nchan

#ifndef CT_CODE
	emis_idx = nchan - k + 1
#else
	emis_idx = k
#endif
	
		call clear_atm_rad(platform, B_profile(k,:), c2_model_info%Ts, sfc_level, &
						surface_emissivity(emis_idx), &
						rtm_trans_atm_clr(k,:), rtm_rad_atm_clr(k,:), rtm_cloud_prof(k, :), &
						channels(k), rtm_rad_clr, rtm_bt_clr)
																		
#ifndef CT_CODE

		call clear_atm_rad(platform, B_profile(k,:), c2_model_info%Ts, sfc_level, &
						surface_emissivity(emis_idx), &
						rtm_trans_atm_clr_low(k,:), rtm_rad_atm_clr_low(k,:), rtm_cloud_prof_low(k, :), &
						channels(k), rtm_rad_clr, rtm_bt_clr)

		call clear_atm_rad(platform, B_profile(k,:), c2_model_info%Ts, sfc_level, &
						surface_emissivity(emis_idx), &
						rtm_trans_atm_clr_high(k,:), rtm_rad_atm_clr_high(k,:), rtm_cloud_prof_high(k, :), &
						channels(k), rtm_rad_clr, rtm_bt_clr)
						

#endif						
						
					
	end do

 end subroutine get_rtm_parameters


 subroutine clear_atm_rad(platform, B_prof, Tsfc, sfc_level, esfc, &
		rt_trans_atm, rt_rad_atm_clr, rt_cloud_prof, channel, rt_rad_clr, rt_bt_clr)

	character(*), intent(in) :: platform
	real, intent(in) :: B_prof(:), rt_trans_atm(:) 
	real, intent(in) :: Tsfc, esfc
	integer, intent(in) :: sfc_level, channel
	real, intent(inout) :: rt_rad_clr, rt_bt_clr, rt_rad_atm_clr(:), rt_cloud_prof(:)
	
	integer :: i
	real :: B1, B2, dtrn
	
	rt_rad_atm_clr(1:start_level) = 0.
	
!	B1 = modis_planck(platform, temp_profile(start_level), channel, 1)
	B1 = B_prof(start_level)
	rt_rad_atm_clr(start_level) = 0.
	rt_cloud_prof(start_level) = B1*rt_trans_atm(start_level)
	
	do i=start_level+1, sfc_level
		
!		B2 = modis_planck(platform, temp_profile(i), channel, 1)
		B2 = B_prof(i)
		dtrn = -(rt_trans_atm(i) - rt_trans_atm(i-1))
		rt_rad_atm_clr(i) = rt_rad_atm_clr(i-1) + (B1 + B2) / 2.0  * dtrn
		rt_cloud_prof(i) = rt_rad_atm_clr(i) + B2*rt_trans_atm(i)
		B1 = B2

	end do
	
!	call clear_toa_rad(platform, rt_rad_atm_clr(sfc_level), rt_trans_atm(sfc_level), Tsfc, &
!						esfc, channel, rt_rad_clr, rt_bt_clr, .false.)
	

 end subroutine clear_atm_rad

 subroutine get_clear_toa_rad(platform,  Tsfc,  esfc, sfc_level, rad_clr, bt_clr, clear_rad_table, clear_trans_table, PRN)
			
	use mod06_run_settings
					
	character(*), intent(in) :: platform
	real, intent(in) :: esfc(:)
	integer, intent(in) :: sfc_level
	real, intent(inout) :: rad_clr(:), bt_clr(:)
	real, intent(in) :: clear_rad_table(:,:), clear_trans_table(:,:)
	integer :: channels(nchan)
	integer :: k
	real, intent(in) :: Tsfc
	logical, intent(in) :: PRN
	
!	channels(1) = bands(band_0370)
#ifndef CT_CODE
	channels(1) = set_bands(band_1100)
	channels(2) = set_bands(band_0370)
#ifdef VIIRS_INST
	if (MODIS_MODE) then 
		channels(1) = set_bands_modis(band_1100)
		channels(2) = set_bands_modis(band_0370)
	endif
#endif	


#else

#if !MAS_INST & !EMAS_INST
	channels = set_bands
#else
	channels = set_bands(1:set_number_of_bands-3)
#endif

#endif

! rtm_rad_atm_clr, rtm_trans_atm_clr
	do k=1, nchan
		call clear_toa_rad(platform,  clear_rad_table(k, sfc_level), clear_trans_table(k, sfc_level), Tsfc, &
			esfc(k), channels(k), rad_clr(k), bt_clr(k), PRN )
	end do
	
			
 end subroutine get_clear_toa_rad


 subroutine clear_toa_rad(platform, rad_atm, tau_atm, tsfc, &
		esfc, channel, rad_clr, bt_clr, PRN)
	
	character(*), intent(in) :: platform
	real, intent(in) :: rad_atm, tau_atm, tsfc, esfc
	real, intent(inout) :: rad_clr, bt_clr
	integer, intent(in) :: channel
	logical, intent(in) :: PRN

	real :: rad_temp

		
	rad_temp = modis_planck(platform, tsfc, channel, 1)

	rad_clr = rad_atm + esfc*rad_temp*tau_atm
		
	bt_clr = modis_bright(platform, rad_clr, channel, 1)


 end subroutine clear_toa_rad

end module rtm_support
