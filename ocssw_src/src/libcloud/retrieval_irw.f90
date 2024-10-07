module retrieval_irw

	use names, only : MYMONTH
	use science_parameters
#if MODIS_OD || MAS_OD || SEV_PR06OD || VIIRS_OD || AHI_PR06OD
	use core_arrays
#else
	use ct_core_arrays
#endif
	use mod06_run_settings
	use global_model_grids
	use rtm_support
	use modis_numerical_module, only : linearinterpolation
	use planck_functions

	implicit none

	private

	public :: retrieve_irw_temp, init_irw, retrieve_irw_temp_for_uncert

     real, parameter :: ma(60) = (/ &
                             2.9769800876, -0.0515871084,  0.0027409105,  0.0001135740,  0.00000113040 , &
                             3.3483238557,  0.1372575458,  0.0133258850,  0.0003043608,  0.00000218650 ,&
                             2.4060295675,  0.0372001609,  0.0096472724,  0.0002334206,  0.00000165450  ,&
                             2.6522386726,  0.0325728824,  0.0100892620,  0.0002601226,  0.00000198560  ,&
                             1.9578262599, -0.2112028966, -0.0057943564, -0.0001050464, -0.00000074313  ,&
                             2.7659753980, -0.1186500984,  0.0011626989,  0.0000936998,  0.00000101060 ,&
                             2.1106811602, -0.3073665907, -0.0090862456, -0.0000889596,  0.00000003552 ,&
                             3.0982173723, -0.1629588435, -0.0020384299,  0.0000286274,  0.00000060283 ,&
                             3.0760551826, -0.2043463270, -0.0053969994, -0.0000541329, -0.00000001768 ,&
                             3.6377215316, -0.0857783614,  0.0024313179,  0.0001495010,  0.00000170850 ,&
                             3.3206165420, -0.1411094234, -0.0026068389,  0.0000057937,  0.00000042113 , & 
                             3.0526632533, -0.1121521836, -0.0009912556,  0.0000179690,  0.00000027070   /)

      real, parameter::  mb(60) = (/&
                            2.9426577089, -0.0510674066,  0.0052420293,  0.0001097927, -0.00000372380 ,&
                            2.6499605646, -0.0105152229,  0.0042895903,  0.0000719741, -0.00000066735 ,&
                            2.3652046763,  0.0141129341,  0.0059242144, -0.0000158816, -0.00000265790 ,&
                            2.5433158163, -0.0046876415,  0.0059325140,  0.0000143938, -0.00000346360 ,&
                            2.4994027830, -0.0364706332,  0.0082001522,  0.0000843577, -0.00000768780 ,&
                            2.7641495752, -0.0728625243,  0.0088877822,  0.0001767765, -0.00001168390 ,&
                            3.1202042743, -0.1002374752,  0.0064054059,  0.0002620230, -0.00001078950 ,&
                            3.4331195144, -0.1021765880,  0.0010498850,  0.0001614861,  0.00000510150 ,&
                            3.4539389485, -0.1158261776,  0.0015449592,  0.0001711651,  0.00000248080 ,&
                            3.6013336912, -0.0775800028,  0.0041940388,  0.0000941307, -0.00000408720 ,&
                            3.1947419143, -0.1045316345,  0.0049986486,  0.0001910731, -0.00000505500 ,&
                            3.1276377012, -0.0707628268,  0.0055532926,  0.0001549833, -0.00000570980 /)

      real, parameter ::  mc(60) = (/ &
                            1.9009562748,  0.0236905223,  0.0086504022, -0.0002167013,  0.00000151230 ,&
                            2.4878735828, -0.0076514110,  0.0079443995, -0.0001773726,  0.00000114730 ,&
                            3.1251275103, -0.1214572133,  0.0146488407, -0.0003187508,  0.00000210290 ,&
                           13.3931706579, -1.2206947755,  0.0560380539, -0.0009873591,  0.00000598210 ,&
                            1.6432070460,  0.1151206937,  0.0033130967, -0.0001458434,  0.00000128610 ,&
                           -5.2366360253,  1.0105574562, -0.0355440449,  0.0005187964, -0.00000262410 ,&
                           -4.7396480830,  0.9625734101, -0.0355846807,  0.0005522497, -0.00000299860 ,&
                           -1.4424842734,  0.4769307320, -0.0139027010,  0.0001758823, -0.00000079846 ,&
                           -3.7140186247,  0.6720953861, -0.0210548327,  0.0002974491, -0.00000149380 ,&
                            8.2237401369, -0.5127532741,  0.0205285436, -0.0003015662,  0.00000157680 ,&
                           -0.4502046794,  0.2629679617, -0.0018419395, -0.0000368887,  0.00000048223 ,&
                            9.3930897423, -0.8836682302,  0.0460453172, -0.0008450362,  0.00000517810  /)

      real, parameter :: bs(24) = (/ &
                            - 3.8,  22.1 ,&
                            -21.5,  12.8 ,&
                            - 2.8,  10.7 ,&
                            -23.4,  29.4 ,&
                            -12.3,  14.9 ,&
                            - 7.0,  16.8 ,&
                            -10.5,  15.0 ,&
                            - 7.8,  19.5 ,&
                            - 8.6,  17.4 ,&
                            - 7.0,  27.0 ,&
                            - 9.2,  22.0 ,&
                            - 3.7,  19.0  /)
                        
	real, parameter :: PBOT = 700.
	real, parameter :: PTOP = 100.




	real :: month_coeffsa(5, 12), month_coeffsb(5, 12), month_coeffsc(5, 12), breakpts(2, 12)


contains

subroutine init_irw

	month_coeffsa = reshape(ma, (/5, 12/))
	month_coeffsb = reshape(mb, (/5, 12/))
	month_coeffsc = reshape(mc, (/5, 12/))
	breakpts = reshape(bs, (/2, 12 /))

end subroutine init_irw

subroutine retrieve_irw_temp_for_uncert(x, y, I11_meas, idx_i, idx_j, irw_temp_low, irw_temp_high )
#ifndef CT_CODE
	use rtm_support, only: rtm_rad_atm_clr_low, rtm_trans_atm_clr_low, rtm_cloud_prof_low, &
						   rtm_rad_atm_clr_high, rtm_trans_atm_clr_high, rtm_cloud_prof_high


	use specific_other, only: set_esfc, set_IRW_channel
	use mod06_run_settings
#endif

	integer, intent(in) :: x, y
	real, intent(inout) :: irw_temp_low, irw_temp_high
	integer, intent(in) :: idx_i, idx_j
	real, intent(in) :: I11_meas 

#ifndef CT_CODE
	real :: sfc_tmp
	integer :: idx_hi, idx_lo, idx_hi_mean, idx_lo_mean
	logical :: ocean_flag

	real :: esfc(2), clear_sky_radiance_low(2), clear_sky_bt_low(2), &
		clear_sky_radiance_high(2), clear_sky_bt_high(2)

	integer :: trop_lev
	integer :: sfc_lev
	real, dimension(model_levels) :: Iwin_low, Iwin_mean, Iwin_high
	integer :: i, icnt, num_work_lev
	real :: rad0, meas_temp
	real, dimension(:), allocatable :: pp, tt, zz, ww
	
	real ::  cpl, tpl, ppl, cph, tph, pph
	integer :: m, zph, zpl
	real :: dplat, c0, c1, c2, c3, c4, lapse_rate, dt, ctz
	real :: adjustment1,  irw_height_mean
	real :: irw_p_low, irw_p_high

	integer :: IRW_channel

	call set_IRW_channel(IRW_channel)
	
! get the temporally and spatially interpolated surface temperature
	sfc_tmp = surface_temperature(x,y)

! get the calculated 11um clear sky radiance and brightness temperature
	call set_esfc(cloudmask(x,y)%water_surface, x, y, esfc, ocean_flag)
	
	sfc_lev = c2_model_info%surface_level


	call get_clear_toa_rad(platform_name, sfc_tmp, esfc, &
						sfc_lev, clear_sky_radiance_low, clear_sky_bt_low, rtm_rad_atm_clr_low, rtm_trans_atm_clr_low, .false.)
	call get_clear_toa_rad(platform_name, sfc_tmp, esfc, &
						sfc_lev, clear_sky_radiance_high, clear_sky_bt_high, rtm_rad_atm_clr_high, rtm_trans_atm_clr_high, .false.)

!	print*, sfc_lev, clear_sky_radiance(band_1100), clear_sky_bt(band_1100)

! get the tropopause level

	trop_lev =  c2_model_info%trop_level

	Iwin_high = 0.
	Iwin_low = 0.
	icnt = 1

	do i=trop_lev, sfc_lev 
		Iwin_low(icnt) = rtm_cloud_prof_low(IRW_channel, i)
		Iwin_high(icnt) = rtm_cloud_prof_high(IRW_channel, i)
		icnt = icnt + 1
	end do
	num_work_lev = icnt-1

	rad0 = I11_meas
	

! get the profile piece that we will do work on
	allocate( tt(num_work_lev))

	icnt = 1
	do i=trop_lev, sfc_lev
		tt(icnt) = c2_model_info%temp_profile(i)
		icnt = icnt + 1
	end do

	
! if the observed BT is colder than minimum trop temperature, then put cloud at trop
	if (meas_temp < tt(1)) then 
		irw_temp_low = tt(1)
		irw_temp_high = tt(1)
! if the observed BT warmer than TSFC, put the cloud at surface
	else if (meas_temp > sfc_tmp) then 
		irw_temp_low = sfc_tmp
		irw_temp_high = sfc_tmp
	else 
! interpolate the temperature and pressure. This is a little more advanced than what UWisc is doing. They take closest level. 
		call locate_point_top_down(rad0, Iwin_low, num_work_lev, idx_lo, idx_hi)
		irw_temp_low = linearinterpolation( (/ Iwin_low(idx_lo), Iwin_low(idx_hi) /), (/ tt(idx_lo), tt(idx_hi) /), rad0)
		
		call locate_point_top_down(rad0, Iwin_high, num_work_lev, idx_lo, idx_hi)
		irw_temp_high = linearinterpolation( (/ Iwin_high(idx_lo), Iwin_high(idx_hi) /), (/ tt(idx_lo), tt(idx_hi) /), rad0)

	endif

	deallocate(tt)
#endif
end subroutine retrieve_irw_temp_for_uncert


subroutine retrieve_irw_temp(x, y, I11_meas, idx_i, idx_j, clear_rad_table, clear_trans_table, cloud_prof, irw_temp, irw_pressure, &
							irw_height) 
  ! WDR modified to replace model_info()%... with c2_model_info%...

	use specific_other, only: set_esfc, set_IRW_channel
	use mod06_run_settings

	integer, intent(in) :: x, y
	real, intent(inout) :: irw_temp, irw_pressure, irw_height
	integer, intent(in) :: idx_i, idx_j
	real, intent(in) :: I11_meas, clear_rad_table(:,:), clear_trans_table(:,:), cloud_prof(:,:)

	real :: sfc_tmp
	integer :: idx_hi, idx_lo, idx_hi_mean, idx_lo_mean
	logical :: ocean_flag
#ifndef CT_CODE 
	real :: esfc(2), clear_sky_radiance(2), clear_sky_bt(2)
#else
	real :: esfc(set_number_of_bands), &
				clear_sky_radiance(set_number_of_bands), clear_sky_bt(set_number_of_bands)
#endif
	integer :: trop_lev
	integer :: sfc_lev
	real, dimension(model_levels) :: Iwin, Iwin_mean
	integer :: i, icnt, num_work_lev
	real :: rad0, meas_temp
	real, dimension(:), allocatable :: pp, tt, zz, ww
	
	real ::  cpl, tpl, ppl, cph, tph, pph
	integer :: m, zph, zpl, ii
	real :: dplat, c0, c1, c2, c3, c4, lapse_rate, dt, ctz
	real :: adjustment1,  irw_height_mean

	integer :: IRW_channel

	call set_IRW_channel(IRW_channel)

! get the temporally and spatially interpolated surface temperature
	sfc_tmp = surface_temperature(x,y)

! get the info about NWP data	


! get the calculated 11um clear sky radiance and brightness temperature
	call set_esfc(cloudmask(x,y)%water_surface, x, y, esfc, ocean_flag)

	sfc_lev = c2_model_info%surface_level
	
	call get_clear_toa_rad(platform_name, sfc_tmp, esfc, &
						sfc_lev, clear_sky_radiance, clear_sky_bt, clear_rad_table, clear_trans_table, .false.)

! get the tropopause level

	trop_lev =  c2_model_info%trop_level

	Iwin = 0.
	icnt = 1

	do i=trop_lev, sfc_lev 
		Iwin(icnt) = cloud_prof(IRW_channel, i)
		icnt = icnt + 1
	end do
	num_work_lev = icnt-1

	rad0 = I11_meas
	meas_temp = modis_bright(platform_name, rad0, channel_11um, 1)

! get the profile piece that we will do work on
	allocate(pp(num_work_lev), zz(num_work_lev), tt(num_work_lev))

	icnt = 1
	do i=trop_lev, sfc_lev
		pp(icnt) = c2_model_info%pressure_profile(i)
		tt(icnt) = c2_model_info%temp_profile(i)
		zz(icnt) = c2_model_info%height_profile(i)
		icnt = icnt + 1
	end do

! if the observed BT is colder than minimum trop temperature, then put cloud at trop
	if (meas_temp < tt(1)) then 
		irw_temp = tt(1)
		irw_pressure = pp(1)
		irw_height = zz(1)

! if the observed BT warmer than TSFC, put the cloud at surface
	else if (meas_temp > sfc_tmp) then 
		irw_temp = sfc_tmp
		irw_pressure = c2_model_info%Ps
		irw_height = 0.
	else 
! interpolate the temperature and pressure. This is a little more advanced than what UWisc is doing. They take closest level. 
		call locate_point_top_down(rad0, Iwin, num_work_lev, idx_lo, idx_hi)
		irw_temp = linearinterpolation( (/ Iwin(idx_lo), Iwin(idx_hi) /), (/ tt(idx_lo), tt(idx_hi) /), rad0)
		irw_pressure = linearinterpolation( (/ Iwin(idx_lo), Iwin(idx_hi) /), (/ pp(idx_lo), pp(idx_hi) /), rad0)
		irw_height = linearinterpolation( (/ Iwin(idx_lo), Iwin(idx_hi) /), (/ zz(idx_lo), zz(idx_hi) /), rad0)

	endif


	if (irw_height < 0.) then 
		irw_height = 0.
		irw_pressure = c2_model_info%Ps
		irw_temp = sfc_tmp
	endif


! now we do the alternate IRW retrieval for extra-low clouds over ocean only. 

	if (irw_pressure > 600. .and. ocean_flag) then 
		
		m = MYMONTH
				
		dplat = latitude(x,y)

		if (dplat < breakpts(1, m)) then 
			c0 = month_coeffsa(1, m)
			c1 = month_coeffsa(2, m)
			c2 = month_coeffsa(3, m)
			c3 = month_coeffsa(4, m)
			c4 = month_coeffsa(5, m)
		else if (dplat >= breakpts(1, m) .and. dplat <= breakpts(2, m)) then 
			c0 = month_coeffsb(1, m)
			c1 = month_coeffsb(2, m)
			c2 = month_coeffsb(3, m)
			c3 = month_coeffsb(4, m)
			c4 = month_coeffsb(5, m)
		else 
			c0 = month_coeffsc(1, m)
			c1 = month_coeffsc(2, m)
			c2 = month_coeffsc(3, m)
			c3 = month_coeffsc(4, m)
			c4 = month_coeffsc(5, m)
		endif
	
		lapse_rate = c0 + c1*dplat + c2*dplat**2 + c3*dplat**3 + c4*dplat**4
		if (lapse_rate < 2.) lapse_rate = 2.
		if (lapse_rate > 10.) lapse_rate = 10.

		dt = clear_sky_bt(IRW_channel) - meas_temp
		if (dt > 0.) then 
			ctz = (dt / lapse_rate) * 1000.
		else
			ctz = 0.
		endif
		
		if (ctz == 0.) then 
			irw_temp = sfc_tmp
			irw_pressure = c2_model_info%Ps
		else 
			
			do i = sfc_lev, trop_lev, -1
			
				if (c2_model_info%height_profile(i) > ctz) then 
					zpl = i + 1
					if (zpl > sfc_lev) then 
						cpl = 0.
						zpl = i
						tpl = sfc_tmp
						ppl =  c2_model_info%Ps
					else 
						cpl = c2_model_info%height_profile(zpl)
						tpl = c2_model_info%temp_profile(zpl)
						ppl = c2_model_info%pressure_profile(zpl)				
					endif
			
					zph = i
					cph = c2_model_info%height_profile(zph)
					tph = c2_model_info%temp_profile(zph)
					pph = c2_model_info%pressure_profile(zph)				
					exit
				endif
				
			end do
		
			if (abs(ctz - cph) <= abs(ctz - cpl)) then 
				irw_height = cph
				irw_pressure = pph
			else
				irw_height = cpl
				irw_pressure = ppl
			endif
				
		
		
			if (irw_height < 0.) then 
				irw_height = 0.
				irw_pressure = c2_model_info%Ps
				irw_temp = sfc_tmp
			endif
		
		
		endif
		
	endif

#if CT_CODE && (MAS_INST || EMAS_INST)
	cloud_top_height(x,y) = irw_height
#endif

	deallocate(pp, zz, tt)

end subroutine retrieve_irw_temp



		SUBROUTINE locate_point_top_down(theta,thetaArray,ntheta,iAngLow, iAngHi)
	
		INTEGER,INTENT(IN)::ntheta
		REAL,INTENT(IN)::theta,thetaArray(1:ntheta)
		INTEGER,INTENT(OUT)::iAngLow, iAngHi
	
		INTEGER:: iAng
		integer :: i

		iAngLow = 1
		iAngHi = ntheta

		do i=1, ntheta
			if (theta < thetaArray(i)) then 
				iAngLow = i-1
				iAngHi = i
				exit
			endif
			
		end do
		
		if (iAngLow < 1) iAngLow = 1
		if (ianghi < 1) ianghi = 1
		if (iAngHi > ntheta) iAngHi = ntheta
		if (ianglow > ntheta) ianglow = ntheta
	
		END SUBROUTINE locate_point_top_down


end module retrieval_irw
