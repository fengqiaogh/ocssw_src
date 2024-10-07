#if SEV_PR06OD || VIIRS_OD

#else

 module multi_layer_clouds
!-------------------------------------------------------------------------------
! DESCRIPTION:
!        this module contains the subroutines needed for computation of 
!  the 0.94-based ancillary and the multi-layer map. 
!
!
! PROGRAMMER:
!            Gala Wind (L3 Comm GSI)
!            Climate and Radiation Branch
!            Code 913, NASA Goddard Space Flight Center
!            Greenbelt, Maryland 20771, U.S.A.
!
!            Mark A Gray (L3 Comm GSI)
!            Climate and Radiation Branch
!            Code 913, NASA Goddard Space Flight Center
!            Greenbelt, Maryland 20771, U.S.A.
!------------------------------------------------------------------------------
! REVISION:
! -------------------------------------------------------------------------------------------
!  6.7.04 Initial revision
!    7.04 corrected variable names for MODIS
!         implented logical cloud decision types
!         
!------------------------------------------------------------------------------
 use GeneralAuxType 
 use modis_cloudstructure
 use modis_numerical_module
 use mod06_run_settings
 use core_arrays
 use science_parameters
 use specific_other
 use planck_functions

 implicit none
 private
 
 ! there is only one public routine here
 public :: Compute_Multilayer_Map, Given_T_get_P, find_TPW, find_brightness_T
 
 
 ! All these variables are private to this module.
 INTEGER, PARAMETER :: NumberOfPressureHeights = 10
 REAL, DIMENSION(NumberOfPressureHeights), PARAMETER :: &
 heightScale = (/ 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000. /)

 INTEGER, PARAMETER :: NumberOfPrecipitableWater = 53
 REAL, DIMENSION(NumberOfPrecipitableWater), PARAMETER :: &
 pwScale = (/ 0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, &
                    0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6,1.8, &
                   2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, &
                   4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, &
                   6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, &
                   8.6 /)

 INTEGER, PARAMETER :: NumberOfMiu = 20 
 REAL, DIMENSION(NumberOfMiu), PARAMETER :: &
 miuscale = (/ 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, &
              0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, &
              0.95, 1.00 /) 


 integer, parameter :: nvza = 7
 integer, parameter :: nsza = 8
 integer, parameter :: nsct = 18


!	AVHRR cloud overlap function coefficients
 real*8, dimension(7,8) :: A_win_over, B_win_over, C_win_over, D_win_over, E_win_over
	
	! minimum 11-12um BTD for overlap detection
 real, dimension (7,8) :: MIN_win_over
	
	
	! VIIRS cloud overlap coefficients for water surface
 real*8, dimension(18) :: A_nir_over_water, B_nir_over_water, C_nir_over_water, D_nir_over_water, E_nir_over_water
	
	! VIIRS cloud overlap coefficients for grass surface
 real*8, dimension(18) :: A_nir_over_land, B_nir_over_land, C_nir_over_land, D_nir_over_land, E_nir_over_land


 
 contains



!=====================================
subroutine Compute_Multilayer_Map(platform_name, &
                                  BigTransTable,  &
								  measurements, &
                                  cloud_phase,    &
                                  Baum_phase,     &
                                  p_ncep,         &
                                  mixR_ncep,      &
                                  t_ncep,         &
								  surface, &
                                  MOD06_Pc,       &
                                  MOD06_PW,       &
                                  sat_zen,            &
                                  sol_zen,           &
								  rel_az, &
                                  tau, tau1621, xpoint, ypoint,           &
                                  mlayer, mlayer_test) 
!=====================================

!-----------------------------------------------------------------------
! !f90
!
! !Description:
!     This subroutine computes the multi-layer cloud value for one pixel
!
! !Input Parameters:
!     platform_name -- character(*) -- Terra or Aqua
!     BigTransTable -- real, dimension(:,:,:,:) -- entire transmittance table
!     measurements -- real, dimension(:) -- MODIS radiances/reflectances
!     cloud_phase   -- logical structure
!     Baum_phase    -- logical structure
!     p_ncep        -- real, dimension(:) -- NCEP pressure profile
!     mixR_ncep     -- real, dimension(:) -- NCEP mixing ratio profile
!     t_ncep        -- real, dimension(:) -- NCEP temperature profile
!     MOD06_Pc      -- real -- Cloud top pressure from MOD06CT
!     MOD06_PW      -- real -- above cloud precip. water derived from MOD06_Pc and NCEP by MOD06OD ancillary
!     sat_zen           -- real --  viewing angle
!     sol_zen          -- real -- solar zenith angle
!     rel_az		-- real -- relative azlimuth angle
!     tau           -- real(:,:) -- retrieved cloud optical depth 3x3 points
!     tau1621           -- real(:,:) -- retrieved cloud optical depth 1.6-2.1um 3x3 points
!     xpoint		-- integer -- x coordinate of pixel
!     ypoint        -- integer -- y coordinate of pixel
!
! !Output Parameters:
!     mlayer        -- integer*1 -- multi-layer cloud map for one pixel
!     mlayer_test   -- integer*1 -- multilayer qa
!
! !Revision History:
!     2004/06/02 wind: Gala Wind
!     Revision 1.0 Initial Revision
!     2009/03/25 wind: Gala Wind
!	  Revision 2.0 added Pavolonis-Heidinger and the delta-tau algorithms
!     2009/05/19 wind: Gala Wind
!     Revision 2.1 added the multilayer qa
!
!
! !Team-Unique Header:
!    Developed by the  Cloud Retrieval Group, NASA GSFC, Greenbelt, Maryland, USA.
!
! !References and Credits:
!      Written by:
!      Gala Wind
!   L3 GSI
!      Code 913, NASA/GSFC
!      Greenbelt, MD 20771
!      wind@climate.gsfc.nasa.gov
!
! !Design Notes:
! NONE
!
! !END
!
!-----------------------------------------------------------------------
  use ch_xfr, only : c2_cmp_there

! first of all we need to compute the quantities that usually come from the research retrieval 
! code. That would be the total column water and 0.94 pw
  character*(*),intent(in)        ::  platform_name
  real(single), dimension(:,:,:,:), intent(in) :: BigTransTable
  real, dimension(:), intent(in) :: measurements
  real, intent(in) :: MOD06_Pc, MOD06_PW, sat_zen, sol_zen, rel_az
  real, intent(in) :: tau, tau1621
  type(processflag), intent(in)                :: cloud_phase
  integer*1, intent(in) :: surface
  integer, intent(in) :: xpoint, ypoint
  type(cloudphase),intent(in)                  :: Baum_phase
  real(single), dimension(:), intent(in)       :: p_ncep, t_ncep, mixR_ncep

  integer(integer_onebyte), intent(out)        :: mlayer, mlayer_test
  real(single)                     :: refl_86, refl_94, IIR_data, refl_65, refl_12


  real :: T, P, PW, TotalPrecipWater, newT, temp, PW_at_900, miu, miu0
  logical :: FIRST
  real :: plant_ratio, pw_fraction, pw_fraction_900, ice_ratio
  integer(integer_onebyte) :: PW_test, PW_900_test, tau_test, cloud_phase_test
  integer :: ph_test, delta_tau
  integer :: DB_CO2
  real :: DB_TH
  integer :: mlayer_add
  logical :: ice_ratio_var
  

	integer :: i, n, j

  
	miu = cos(sat_zen * d2r)
	miu0 = cos(sol_zen * d2r)
  
   refl_65 = measurements(band_0065)
   refl_86 = measurements(band_0086)
   refl_94 = measurements(band_0935)
   refl_12 = measurements(band_0124)
   IIR_data = measurements(band_1100)
  
  ! Get the initial guess for brightness temperature
  T = modis_bright(platform_name,IIR_data,channel_11um,1)

  ! compute the initial pressure

  call Given_T_get_P(P, T, surface, p_ncep, t_ncep)

  ! find the initial amount of precipitable water based on the initial guess 
  ! for the cloud top pressure computed above

  call find_TPW(refl_86, refl_94, P, PW,  BigTransTable, miu, miu0, xpoint, ypoint)

    FIRST = .TRUE.
    call find_brightness_T(platform_name,PW, P, newT, TotalPrecipWater, p_ncep, t_ncep, mixR_ncep, surface, &
              IIR_data, BigTransTable, miu, FIRST)


    call Given_T_get_P(P, newT, surface, p_ncep, t_ncep)
    ! do another iteration 
    call find_TPW(refl_86, refl_94, P, PW,  BigTransTable, miu, miu0, xpoint, ypoint)
    call find_brightness_T(platform_name,PW, P, newT, temp, p_ncep, t_ncep, mixR_ncep, surface, &
              IIR_data, BigTransTable, miu, FIRST)

    T = newT

	if (tau >= 4.) precip_water_094(xpoint, ypoint) = PW


    ! now we're going to compute PW as if all the clouds were at 900 mb


	if (p_ncep(surface) >= 900.) then 
		call find_TPW(refl_86, refl_94, 900., PW_at_900,  BigTransTable, miu, miu0, xpoint, ypoint)
	else 
		call find_TPW(refl_86, refl_94, p_ncep(surface) - 100., PW_at_900, BigTransTable, miu, miu0, xpoint, ypoint)
	endif



    ! -----------------------------------------------------------------------
    !   NOW ANCILLARY PORTION IS FINISHED. STARTING LAYER DETECTION
    ! -----------------------------------------------------------------------

	if (PW > TotalPrecipWater) PW = TotalPrecipWater
	if (PW_at_900 > TotalPrecipWater) PW_at_900 = TotalPrecipWater

    plant_ratio = refl_86 / refl_65 ! vegetation ratio: 0.86/0.65 um
    pw_fraction = abs(MOD06_PW - PW) / TotalPrecipWater  ! precipitable water
	pw_fraction_900 = abs(MOD06_PW - PW_at_900) / TotalPrecipWater ! additional precip water ratio
    ice_ratio = refl_86 / refl_12 ! ice ratio 0.86/1.2 um
    
	ice_ratio_var = set_ice_ratio(ice_ratio)
	
	! there has to be enough precipitable water in the column for us to bother with retrieval. 
	
    if (pw_fraction > 0.08 .and. MOD06_Pc < 550.0 .and. plant_ratio < 1.25 .and. ice_ratio_var  ) then 
      PW_test = 1
    else 
      PW_test = 0
    endif 
    
    if (pw_fraction_900 > 0.08 .and. MOD06_Pc < 550.0 .and. plant_ratio < 1.25 .and. ice_ratio_var  ) then 
      PW_900_test = 1
    else 
      PW_900_test = 0
    endif 


! use the midpoint as that's what the original point would've been 
    if (tau < 4. .and. tau >= 0.) then 
      tau_test = 1
    else
      tau_test = 0
    endif
    
    ! cloud phase test: 0 if phases agree, 1 if SWIR-water and IR-ice or IR-mixed,
    ! 2 if SWIR-ice and IR-water
    cloud_phase_test = 0
    if (cloud_phase%watercloud .and. Baum_phase%watercloud == 1 ) cloud_phase_test = 0
    if (cloud_phase%icecloud .and. Baum_phase%icecloud == 1) cloud_phase_test = 0

    if (cloud_phase%watercloud .and. Baum_phase%icecloud == 1) cloud_phase_test = 1
    if (cloud_phase%icecloud .and. Baum_phase%watercloud == 1) cloud_phase_test = 2
    
	
	if (tau_test == 1) then 
		cloud_phase_test = 0
		PW_900_test = 0
		PW_test = 0
	endif
	

! the don't bother CO2 test for ice cloud layer too thick
	if ( abs(latitude(xpoint, ypoint)) < 30.) then 
		DB_TH = 1.7
	else
		DB_TH = 2.43
	endif
	
   if ( ( c2_cmp_there(band_1350) == 1 ) .and. &
        ( c2_cmp_there(band_1200) == 1 ) ) then
     if (measurements(band_1350) < DB_TH) then 
       DB_CO2 = 1 ! don't bother
     else
       DB_CO2 = 0 ! bother
     endif
   else
     DB_CO2 = 1 ! don't bother
   endif

        delta_tau = 0
        ph_test = 0

	if (DB_CO2 == 0 .and. tau_test == 0) then ! not too thick, not too thin, just right
	
! this is the Pavolonis-Heidinger algorithm	
		ph_test =  PH_multilayer (sol_zen, sat_zen, rel_az, measurements, cloud_phase, platform_name, latitude(xpoint, ypoint))
	
! this is the delta-tau method
		if (cloud_phase%icecloud) then 
				if (tau < 30. .and. tau1621 > 80.) then 
					delta_tau = 1
				else 
					delta_tau = 0
				endif
		endif
	else
		delta_tau = 0
		ph_test = 0
	endif
	
! calculate the confidence level of the multilayer answer. 
! confidence levels like this: phase MLC5 test 1 point, the other MLC5 tests 2 points, two MLC5 tests: 3 points, three MLC5 tests 4 points, delta-tau 1 point, 
! ph_test 3 points


    mlayer = 1 !to 1
	mlayer_test = 0 ! no tests done

	if (cloud_phase_test == 1) mlayer = mlayer + 1
	if (delta_tau == 1) mlayer = mlayer + 1
	if (PW_test == 1) mlayer = mlayer + 2
	if (PW_900_test == 1) mlayer = mlayer + 2
	if (ph_test == 1) mlayer = mlayer + 3

	if (cloud_phase_test == 1) mlayer_test = ibset(mlayer_test, 0)
	if (PW_test == 1) mlayer_test = ibset(mlayer_test, 1)
	if (PW_900_test == 1) mlayer_test = ibset(mlayer_test, 2)
	if (delta_tau == 1) mlayer_test = ibset(mlayer_test, 3)
	if (ph_test == 1) mlayer_test = ibset(mlayer_test, 4)
  
  end subroutine Compute_Multilayer_Map
 
 
 !=====================================
    subroutine find_TPW(refl_86, refl_94, P, PW, BigTransTable, miu1, miu0, xpoint, ypoint)
!=====================================
!-----------------------------------------------------------------------
! !f90
!
! !Description:
!     This subroutine computes 0.94um-based precipitable water above cloud for a given pixel
!
! !Input Parameters:
!     refl_86 -- real*4 -- 0.86um reflectance (band 2 in MODIS)
!     refl_94 -- real*4 -- 0.94um reflectance (band 9 in MODIS)
!     P -- real*4 -- cloud top pressure in mb
!     BigTransTable -- real*4, dimension(:,:,:,:) -- transmittance table
!     miu1 -- real*4 -- cos of viewing angle
!     miu0 -- real*4 -- cos of solar zenith angle
!
! !Output Parameters:
!     PW -- real*4 -- precipitable water amount in cm
!
! !Revision History:
!     2004/06/02 wind: Gala Wind
!     Revision 1.0 Initial Revision
!
! !Team-Unique Header:
!    Developed by the  Cloud Retrieval Group, NASA GSFC, Greenbelt, Maryland, USA.
!
! !References and Credits:
!      Written by:
!      Gala Wind
!    L3 GSI
!      Code 913, NASA/GSFC
!      Greenbelt, MD 20771
!      wind@climate.gsfc.nasa.gov
!
! !Design Notes:
!   NONE
!
! !END
!
!-----------------------------------------------------------------------

    use GeneralAuxType

    implicit none
            
    real, intent(in) :: P, miu1, miu0, refl_86, refl_94
    real, intent(inout) :: PW
    real(single), dimension(:,:,:,:), intent(in) :: BigTransTable
    integer, intent(in) :: xpoint, ypoint    
		
    integer :: k
        
    real :: tpw
    real, dimension(:), allocatable :: trans2way


    integer, dimension(4) :: mindex

    integer :: miu_index, pindex, pwindex, miu2index, first_pIndex, second_pIndex
  integer :: bandindexmapsw1(2)

    REAL :: miu, transValue



    real, dimension(53) :: pix086, pix094, ans_diff, ans_diff_abs
    real :: answer, trans_mid1, trans_mid2
    integer :: first_miuindex, second_miuindex
    real :: trans_p_1, trans_p_2 ! these are for the second pIndex
    real :: trans1, trans2
	
	real :: point1, point2, x1, x2
	integer :: idx1, idx2

  bandindexmapsw1(1) = 2
  bandindexmapsw1(2) = 3


   pIndex = NINT( (P - heightscale(1)) / 100.0 ) + 1
   if(pIndex < 1) then
     pIndex = 1
     first_pIndex = 1
     second_pIndex = 1
   elseif(pIndex > NumberOfPressureHeights) then
     pIndex = NumberOfPressureHeights
     first_pIndex = pIndex
     second_pIndex = pIndex
   endif

   if (heightscale(pIndex) >= P) then 
     first_pIndex = pIndex-1
     second_pIndex = pIndex
   else
     first_pIndex = pIndex
     second_pIndex = pIndex+1
   endif
 
   if (first_pIndex == 0 ) then 
     first_pIndex = 1
     second_pIndex = 2
   endif

   if (second_pIndex > NumberOfPressureHeights) &
		second_pIndex = NumberOfPressureHeights

   ! key: compute effective miu for a two way path:

    miu = (miu0*miu1)/(miu0+miu1)

!
!... find array indicies along p, pw and miu dimension:
!
  
        miu2Index = NINT( (miu - miuscale(1)) / 0.05 ) + 1
        if (miu2Index < 1) then
            miu2Index = 1
        elseif(miu2Index > NumberOfMiu) then
            miu2Index = NumberOfMiu
        endif
        

!
!... 2-way case (for shortwave channels upto 3.7 micron band):
!

   mindex(1) = pindex
   mindex(2) = 1
   mindex(3) = miu2Index
   mindex(4) = bandindexmapsw1(1)

        if (miuscale(miu2index) >= miu ) then
                first_miuindex = miu2index -1 
                second_miuindex = miu2index 
        else
                first_miuindex = miu2index
                second_miuindex = miu2index + 1
        endif

    do k=1, numberofprecipitablewater
!     for the first pressure
!      call interpolate(miuscale(first_miuindex), &
!                       bigtransTable(first_pIndex, k, first_miuindex, bandindexmapsw1(1)), &
!                       miuscale(second_miuindex), &
!                       bigtransTable(first_pIndex, k, second_miuindex, bandindexmapsw1(1)), &
!                       miu, trans_mid1)
       trans_mid1 = linearinterpolation ( (/miuscale(first_miuindex), miuscale(second_miuindex)/), &
                    (/bigtransTable(first_pIndex, k, first_miuindex, bandindexmapsw1(1)), &
                      bigtransTable(first_pIndex, k, second_miuindex, bandindexmapsw1(1))/), &
                    miu)
      
!      call interpolate(miuscale(first_miuindex), &
!                       bigtransTable(first_pIndex, k, first_miuindex, bandindexmapsw1(2)), &
!                       miuscale(second_miuindex), &
!                       bigtransTable(first_pIndex, k, second_miuindex, bandindexmapsw1(2)), &
!                       miu, trans_mid2)

      trans_mid2 = linearinterpolation( (/miuscale(first_miuindex), miuscale(second_miuindex)/), &
                   (/bigtransTable(first_pIndex, k, first_miuindex, bandindexmapsw1(2)), &
                    bigtransTable(first_pIndex, k, second_miuindex, bandindexmapsw1(2))/), &
                    miu) 
!     for the second pressure
!      call interpolate(miuscale(first_miuindex), &
!                       bigtransTable(second_pIndex, k, first_miuindex, bandindexmapsw1(1)), &
!                       miuscale(second_miuindex), &
!                       bigtransTable(second_pIndex, k, second_miuindex, bandindexmapsw1(1)), &
!                       miu, trans_p_1)

      trans_p_1 = linearinterpolation( (/miuscale(first_miuindex), miuscale(second_miuindex)/), &
                   (/ bigtransTable(second_pIndex, k, first_miuindex, bandindexmapsw1(1)), &
                      bigtransTable(second_pIndex, k, second_miuindex, bandindexmapsw1(1)) /), &
                      miu)
!      call interpolate(miuscale(first_miuindex), &
!                       bigtransTable(second_pIndex, k, first_miuindex, bandindexmapsw1(2)), &
!                       miuscale(second_miuindex), &
!                       bigtransTable(second_pIndex, k, second_miuindex, bandindexmapsw1(2)), &
!                       miu, trans_p_2)
      trans_p_2 = linearinterpolation( (/ miuscale(first_miuindex),  miuscale(second_miuindex) /), &
                     (/ bigtransTable(second_pIndex, k, first_miuindex, bandindexmapsw1(2)), &
                        bigtransTable(second_pIndex, k, second_miuindex, bandindexmapsw1(2)) /), &
                        miu)
 
!     final interpolation now along the pressure axis
      trans1 = linearinterpolation ((/heightscale(first_pindex), heightscale(second_pindex)/), &
                                    (/trans_mid1, trans_p_1/), &
                                    P)   
      trans2 = linearinterpolation ((/heightscale(first_pindex), heightscale(second_pindex)/), &
                                    (/trans_mid2, trans_p_2/), &  
                                    P)
      if (trans_mid1 < 0. .or. trans_p_1 < 0.) trans1 = -1.0
      if (trans_mid2 < 0. .or. trans_p_2 < 0.) trans2 = -1.0
    
      pix086(k) = refl_86 / trans1
      pix094(k) = refl_94 / trans2
    	
      if (trans1 < 0.0) pix086(k) = 20
      if (trans2 < 0.0) pix094(k) = 300
            
    end do

	do k=1, numberofprecipitablewater
		ans_diff(k) = pix086(k)-pix094(k)
		ans_diff_abs(k) = ans_diff(k)
		if (ans_diff(k) < 0) ans_diff_abs(k) = -ans_diff(k)
	end do



    answer = minval(ans_diff_abs)
    do k=1, numberofprecipitablewater
      if (real_s_equal(ans_diff_abs(k), answer)) then
        answer = pwscale(k)
        exit
      endif    
    end do

	if (ans_diff(k) < 0.) then 
	
		idx1 = k-1
		if (idx1 < 1) idx1 = 1
		idx2 = k
	
		if (idx1 == idx2) then 
			PW = pwscale(k)
			return
		endif
	
		point1 = ans_diff(idx1)
		point2 = ans_diff(idx2)
			
		x1 = pwscale(idx1)
		x2 = pwscale(idx2)
			
		if (point1 == -280.) then 
			PW = pwscale(idx2)
		else
			PW = linearinterpolation( (/ point1, point2 /), (/ x1, x2 /), 0.0)			
		endif
			
	else
		idx1 = k
		idx2 = k+1
		if (idx2 > numberofprecipitablewater) idx2 = numberofprecipitablewater

		if (idx1 == idx2) then 
			PW = pwscale(k)
			return
		endif

		point1 = ans_diff(idx1)
		point2 = ans_diff(idx2)
			
		x1 = pwscale(idx1)
		x2 = pwscale(idx2)
			
		if (point2 == -280.) then 
			PW = pwscale(idx1)
		else
			PW = linearinterpolation( (/ point1, point2 /), (/ x1, x2 /), 0.0)
		endif
		
	 endif

    end subroutine find_TPW


!=====================================
  subroutine Given_T_get_P(P,T, surface, p_ncep, t_ncep)
!=====================================

!-----------------------------------------------------------------------
! !f90
!
! !Description:
!     This subroutine gets NCEP cloud top pressure for a given 
!     cloud top temperature
!
! !Input Parameters:
!     T -- real*4 -- cloud top temperature in K
!     p_ncep -- real*4, dimension(:) -- NCEP pressure profile
!     t_ncep -- real*4, dimension(:) -- NCEP temperature profile
!
! !Output Parameters:
!     P -- real*4 -- cloud top pressure in mb
!
! !Revision History:
!     2004/06/02 wind: Gala Wind
!     Revision 1.0 Initial Revision
!
! !Team-Unique Header:
!    Developed by the  Cloud Retrieval Group, NASA GSFC, Greenbelt, Maryland, USA.
!
! !References and Credits:
!      Written by:
!      Gala Wind
!    L3 GSI
!      Code 913, NASA/GSFC
!      Greenbelt, MD 20771
!      wind@climate.gsfc.nasa.gov
!
! !Design Notes:
!   NONE
!
! !END
!
!-----------------------------------------------------------------------
    real, intent(in)  ::  T
    real, intent(inout) :: P
	integer*1, intent(in) :: surface
    real, dimension(:), intent(in) :: p_ncep, t_ncep
    
  integer :: k
  real :: factor
   
!  do k=1, model_levels
  do k=surface, 1, -1
    if ( T > t_ncep(k)) &
      exit
  end do
    
  if (k <= 1) then 
    P = p_ncep(1)
  elseif (k == surface) then 
    P = p_ncep(surface)
  else 

	P = linearinterpolation( (/ t_ncep(k), t_ncep(k+1) /), (/ p_ncep(k), p_ncep(k+1) /), T)

!    factor = (T - t_ncep(k)) / (t_ncep(k-1) - t_ncep(k))
!    P = p_ncep(k) + factor* (p_ncep(k-1) - p_ncep(k))

  endif

  if (P < 100.) P = 100. 


  end subroutine Given_T_get_P

 
 
!=====================================
    subroutine find_brightness_T(platform_name,PW, P, newT, TotalPrecipWater, p_ncep, t_ncep, mixR_ncep, surface, &
              IIR_data, BigTransTable, miu, FIRST)
!=====================================

!-----------------------------------------------------------------------
! !f90
!
! !Description:
!     This subroutine computes the cloud top temperature for a given pixel
!
! !Input Parameters:
!     P -- real*4 -- cloud top pressure in mb
!     PW -- real*4 -- precipitable water amount in cm
!     p_ncep -- real*4, dimension(:) -- NCEP pressure profile
!     t_ncep -- real*4, dimension(:) -- NCEP temperature profile
!     mixR_ncep -- real*4, dimension(:) -- NCEP mixing ratio profile
!     IIR_data -- real*4 -- 11um radiance
!     BigTransTable -- real*4, dimension(:,:,:,:) -- transmittance table
!     miu -- real*4 -- cosine of viewing angle
!     FIRST -- logical -- flag for integration of TPW
!
! !Output Parameters:
!     newT -- real*4 -- cloud top temperature in K
!     TotalPrecipWater -- real*4 -- integrated NCEP water vapor column amount
!
! !Revision History:
!     2004/06/02 wind: Gala Wind
!     Revision 1.0 Initial Revision
!
! !Team-Unique Header:
!    Developed by the  Cloud Retrieval Group, NASA GSFC, Greenbelt, Maryland, USA.
!
! !References and Credits:
!      Written by:
!      Gala Wind
!    L3 GSI
!      Code 913, NASA/GSFC
!      Greenbelt, MD 20771
!      wind@climate.gsfc.nasa.gov
!
! !Design Notes:
!   NONE
!
! !END
!
!-----------------------------------------------------------------------

  implicit NONE
  character*(*),intent(in)        ::  platform_name    
  real,    intent(in)                          :: PW, P, IIR_data, miu
  logical, intent(inout)                       :: FIRST
  integer*1, intent(in) :: surface
  real,    intent(inout)                       :: newT, TotalPrecipWater
  real(single), dimension(:), intent(in)       :: p_ncep, t_ncep, mixR_ncep
  real(single), dimension(:,:,:,:), intent(in) :: BigTransTable
        

  integer :: pIndex,  miu2index, pwIndex, first_miuindex, second_miuindex, k
  real ::  tpw, trans_mid1
  integer :: mindex(4), first_pIndex, second_pIndex
  real :: trans_p_1 ! these are for the second pIndex
  real :: trans1
  real :: total_PW, PW_layer, sum_T, slope,  mixR_lower, mixR_upper, T_lower, T_upper, T_layer
  real :: tmeanAbovecloud, tpwAboveCloud
  real :: B_11_crude
  real :: temp
  integer :: bandindexmaplw1(1)

  bandindexmaplw1(1) = 8

  pIndex = NINT( (P - heightscale(1)) / 100.0 ) + 1
  if(pIndex < 1) then
    pIndex = 1
    first_pIndex = 1
    second_pIndex = 1
  elseif(pIndex > NumberOfPressureHeights) then
    pIndex = NumberOfPressureHeights
    first_pIndex = pIndex
    second_pIndex = pIndex
  endif

  if (heightscale(pIndex) >= P) then 
    first_pIndex = pIndex-1
    second_pIndex = pIndex
  else
    first_pIndex = pIndex
    second_pIndex = pIndex+1
  endif

  if (first_pIndex == 0) then
    first_pIndex = 1
    second_pIndex = 2
  endif

   if (second_pIndex > NumberOfPressureHeights) &
		second_pIndex = NumberOfPressureHeights

!
!... find array indicies along p, pw and miu dimension:
!
   miu2Index = NINT( (miu - miuscale(1)) / 0.05 ) + 1
   if (miu2Index < 1) then
       miu2Index = 1
   elseif(miu2Index > NumberOfMiu) then
       miu2Index = NumberOfMiu
   endif


   if(pw < 0.0) then
       pwIndex = 1
   elseif(pw < 0.2) then
       pwIndex = NINT( (pw- pwScale(1)) / 0.02 ) + 1
   else
       pwIndex = NINT( (pw - pwScale(11)) / 0.20 ) + 11
       if(pwIndex > NumberOfPrecipitableWater) then
          pwIndex = NumberOfPrecipitableWater
       endif
   endif


   mindex(1) = pindex
   mindex(2) = pwIndex
   mindex(3) = miu2Index
   mindex(4) = bandindexmaplw1(1)

   if (miuscale(miu2index) >= miu ) then
     first_miuindex = miu2index -1 
     second_miuindex = miu2index 
   else
     first_miuindex = miu2index
     second_miuindex = miu2index + 1
   endif
!   call interpolate(miuscale(first_miuindex), bigtransTable(first_pIndex, pwindex, first_miuindex, bandindexmaplw1(1)), &
!        miuscale(second_miuindex), bigtransTable(first_pIndex, pwindex, second_miuindex, bandindexmaplw1(1)), &
!        miu, trans_mid1)

   trans_mid1 = linearinterpolation( (/miuscale(first_miuindex), miuscale(second_miuindex)/), &
                (/bigtransTable(first_pIndex, pwindex, first_miuindex, bandindexmaplw1(1)), & 
                   bigtransTable(first_pIndex, pwindex, second_miuindex, bandindexmaplw1(1))/), &
                miu) 


!   call interpolate(miuscale(first_miuindex), bigtransTable(second_pindex, pwindex, first_miuindex, bandindexmaplw1(1)), &
!        miuscale(second_miuindex), bigtransTable(second_pindex, pwindex, second_miuindex, bandindexmaplw1(1)), &
!        miu, trans_p_1)
   trans_p_1= linearinterpolation( (/miuscale(first_miuindex), miuscale(second_miuindex)/), &
               (/bigtransTable(second_pindex, pwindex, first_miuindex, bandindexmaplw1(1)),  &
                 bigtransTable(second_pindex, pwindex, second_miuindex, bandindexmaplw1(1))/), &
                miu)        


!    call interpolate(heightscale(first_pindex), trans_mid1, heightscale(second_pindex), trans_p_1, P, trans1)

    trans1 = linearinterpolation( (/heightscale(first_pindex), heightscale(second_pindex)/), &
                                  (/trans_mid1, trans_p_1/), &
                                   p)

    if (trans_mid1 >= 0. .and. trans_p_1 < 0.) trans1 = trans_p_1
    if (trans_mid1 < 0. .and. trans_p_1 >= 0.) trans1 = trans_mid1



! at this point we have the 1-way transmittance for the atmosphere at 11 microns and the 11 micron radiance
! stored in IIR_data that came from the read_mas_module.f90

! we can go ahead and find the amount of atmospheric correction for the 11 micron channel. 

 
      total_PW = 0.0
      sum_T  = 0.0


!     sep: start loop at NCEP's TOA level - 1
    
          
! add up the water profile in order to get the column amount later used for multi-layer cloud detection. 
      if (FIRST) then

!     do k=1, model_levels-1
      do k=surface, 2, -1
            PW_layer = (p_ncep(k)-p_ncep(k-1)) * (mixR_ncep(k-1) + mixR_ncep(k))*0.5/ 980.616 
            total_PW = total_PW+PW_layer
        end do
        TotalPrecipWater = total_PW
      endif
          
      
      PW_layer = 0.0
     total_PW = 0.0

!      k = model_levels - 1
	  k = 2
!      DO WHILE(p_ncep(k) < P .AND. k > 1)
      DO WHILE(p_ncep(k) < P .AND. k < surface)
      
      !     total precipitable water above cloud in g/cm2
      !   sep: dPW = rho*w*dz, dp=abs(rho*g*dz) or dz=dp/rho*g => dPW=w*dp/g
      !            {where w=layer mixing ratio, rho=density}
                        ! mixing ratio must be in g/kg, P in mb 1.22.04

         PW_layer = (p_ncep(k)-p_ncep(k-1)) * (mixR_ncep(k-1) + mixR_ncep(k))*0.5/ 980.616 

         T_layer  = (t_ncep(k)+t_ncep(k-1))*0.5
         total_PW = total_PW + PW_layer
         sum_T = sum_T + PW_layer*T_layer
                   
!         k = k - 1
         k = k + 1
      END DO

!	if (P > p_ncep(model_levels)) then 
!	   tpwAboveCloud = total_PW
!	   tmeanAbovecloud = sum_T

!	else


! all the (k+1) have now become (k-1)'s
      slope = (mixR_ncep(k) - mixR_ncep(k-1)) / (p_ncep(k) - p_ncep(k-1))
      mixR_upper = mixR_ncep(k-1)
      mixR_lower = slope * (P - p_ncep(k-1)) + mixR_upper
      PW_layer = (P - p_ncep(k-1))*(mixR_upper + mixR_lower)*0.5/980.616
      tpwAboveCloud = total_PW + PW_layer 

      slope = (t_ncep(k) - t_ncep(k-1)) / (p_ncep(k) - p_ncep(k-1))
      T_upper = t_ncep(k-1)
      T_lower = slope * (P - p_ncep(k-1)) + T_upper
      T_layer = (T_upper + T_lower)*0.5
      tmeanAbovecloud = (sum_T + PW_layer*T_layer)/tpwAboveCloud

!	endif

! now we have the mean temperature weighted by the the precip. water
! now we apply the approximation

    B_11_crude = modis_planck(platform_name, tmeanAboveCloud, channel_11um, 1) 
    if (trans1 < -0.0001) trans1 = 1.0
    B_11_crude = B_11_crude * (1. - trans1)   ! scale by 1. - 1waytransmittance

! got the amount of atmospheric correction


! now let's go ahead and subtract the correction from the measured radiance
    temp  = IIR_data - B_11_crude
! now we need to divide this mess by the transmittance that we computed earlier
    temp = temp / trans1

! now all this neat mess needs to be fed into the Planck function inversion
    newT = modis_bright(platform_name,temp,channel_11um,1)

    end subroutine find_brightness_T
	
!=====================================
!=====================================
	
! THESE ROUTINES BELONG TO THE PAVOLONIS-HEIDINGER ALGORITHM

!=====================================
!=====================================

	
!=====================================
	subroutine init_coeffs
!=====================================

!-----------------------------------------------------------------------
! !f90
!
! !Description:
!     This subroutine initializes the coefficient arrays for the Pavolonis-Heidinger algorithm
!
! !Input Parameters:
!		NONE
! !Output Parameters:
!		NONE
!
! !Revision History:
!     2009/03/20 wind: Gala Wind
!     Revision 1.0 Initial Revision
!
! !Team-Unique Header:
!    Developed by the  Cloud Retrieval Group, NASA GSFC, Greenbelt, Maryland, USA.
!
! !References and Credits:
!      Written by:
!      Gala Wind
!      SSAI
!      Code 613.2, NASA/GSFC
!      Greenbelt, MD 20771
!      Gala.Wind@nasa.gov
!
! !Design Notes:
!   NONE
!
! !END
!
!-----------------------------------------------------------------------
			real*8, dimension(56) :: temp
		
		
			temp =  (/0.70,0.70,0.70,0.70,0.75,0.80,0.80, &
                           0.70,0.70,0.70,0.70,0.75,0.80,0.80,&
                           0.70,0.70,0.70,0.70,0.75,0.80,0.80,&
                           0.70,0.70,0.70,0.70,0.75,0.80,0.80,&
                           0.70,0.70,0.70,0.70,0.75,0.80,0.80,&
                           0.70,0.70,0.70,0.70,0.75,0.90,0.90,&
                           0.75,0.75,0.75,0.80,0.80,0.90,0.90,&
                           0.75,0.75,0.75,0.80,0.80,0.90,0.90 /)

			MIN_win_over = reshape(temp, (/ 7,8 /))


			temp =  (/-5.09e+01,-3.83e+01,-4.20e+01,-5.52e+01,-4.80e+01,-3.16e+01,-3.16e+01 , &
                         -5.09e+01,-3.83e+01,-4.20e+01,-5.52e+01,-4.80e+01,-3.16e+01,-3.16e+01 , &
                         -6.02e+01,-4.36e+01,-4.01e+01,-4.14e+01,-4.67e+01,-6.36e+01,-6.36e+01 , &
                         -4.67e+01,-5.81e+01,-4.20e+01,-3.79e+01,-3.25e+01,-4.11e+01,-4.11e+01 , &
                         -6.25e+01,-5.49e+01,-5.42e+01,-4.99e+01,-5.63e+01,-7.13e+01,-7.13e+01 , &
                         -7.72e+01,-6.02e+01,-6.28e+01,-5.13e+01,-4.70e+01,-3.22e+01,-3.22e+01 , &
                         -1.00e+02,-1.12e+02,-8.96e+01,-7.65e+01,-7.04e+01,-8.40e+01,-8.40e+01 , &
                         -2.85e+02,-2.52e+02,-1.49e+02,-1.82e+02,-1.28e+02,-5.21e+01,-5.21e+01 /)
						 
			A_win_over =  reshape(temp, (/ 7,8 /))
			     
			temp =  (/ 8.58e+01, 6.05e+01, 6.69e+01, 9.35e+01, 7.91e+01, 4.24e+01, 4.24e+01 , &
                          8.58e+01, 6.05e+01, 6.69e+01, 9.35e+01, 7.91e+01, 4.24e+01, 4.24e+01 , &
                          1.05e+02, 7.15e+01, 6.50e+01, 6.78e+01, 7.92e+01, 1.07e+02, 1.07e+02 , &
                          7.88e+01, 1.03e+02, 7.16e+01, 6.47e+01, 5.46e+01, 6.82e+01, 6.82e+01 , &
                          1.11e+02, 9.79e+01, 1.01e+02, 9.37e+01, 1.08e+02, 1.33e+02, 1.33e+02 , &
                          1.41e+02, 1.08e+02, 1.20e+02, 1.03e+02, 9.63e+01, 7.01e+01, 7.01e+01 , &
                          1.78e+02, 2.08e+02, 1.70e+02, 1.53e+02, 1.43e+02, 1.87e+02, 1.87e+02 , &
						  5.08e+02, 4.55e+02, 2.75e+02, 3.69e+02, 2.66e+02, 1.40e+02, 1.40e+02 /)
			   

			B_win_over =  reshape(temp, (/ 7,8 /))

			temp =  (/-4.12e+01,-2.41e+01,-2.77e+01,-4.57e+01,-3.60e+01,-8.56e+00,-8.56e+00 , &
                         -4.12e+01,-2.41e+01,-2.77e+01,-4.57e+01,-3.60e+01,-8.56e+00,-8.56e+00 , &
                         -5.47e+01,-3.20e+01,-2.77e+01,-2.96e+01,-3.80e+01,-5.29e+01,-5.29e+01 , &
                         -3.76e+01,-5.51e+01,-3.41e+01,-3.01e+01,-2.40e+01,-3.02e+01,-3.02e+01 , &
                         -6.08e+01,-5.30e+01,-5.77e+01,-5.38e+01,-6.41e+01,-7.67e+01,-7.67e+01 , &
                         -8.14e+01,-5.96e+01,-7.29e+01,-6.59e+01,-6.21e+01,-4.88e+01,-4.88e+01 , &
                         -1.02e+02,-1.27e+02,-1.06e+02,-1.00e+02,-9.37e+01,-1.38e+02,-1.38e+02 , &
                         -3.09e+02,-2.80e+02,-1.69e+02,-2.56e+02,-1.86e+02,-1.23e+02,-1.23e+02 /)

			C_win_over =  reshape(temp, (/ 7,8 /))

			   
			temp =  (/ 9.36e-01,-3.25e+00,-2.60e+00, 1.71e+00,-7.43e-01,-8.27e+00,-8.27e+00 , &
                          9.36e-01,-3.25e+00,-2.60e+00, 1.71e+00,-7.43e-01,-8.27e+00,-8.27e+00 , &
                          4.48e+00,-1.20e+00,-2.31e+00,-1.98e+00, 1.48e-01, 2.65e+00, 2.65e+00 , &
                          3.65e-01, 4.94e+00,-2.40e-01,-1.20e+00,-2.60e+00,-2.09e+00,-2.09e+00 , &
                          6.62e+00, 4.96e+00, 6.72e+00, 5.76e+00, 8.27e+00, 9.71e+00, 9.71e+00 , &
                          1.21e+01, 6.67e+00, 1.12e+01, 1.05e+01, 9.62e+00, 7.24e+00, 7.24e+00 , &
                          1.67e+01, 2.45e+01, 1.99e+01, 1.99e+01, 1.81e+01, 3.37e+01, 3.37e+01 , &
                          6.92e+01, 6.26e+01, 3.53e+01, 6.57e+01, 4.58e+01, 3.56e+01, 3.56e+01 /)

			D_win_over =  reshape(temp, (/ 7,8 /))
			   
			temp =  (/ 2.94e+00, 3.14e+00, 3.15e+00, 3.03e+00, 3.27e+00, 3.77e+00, 3.77e+00 , &
                          2.94e+00, 3.14e+00, 3.15e+00, 3.03e+00, 3.27e+00, 3.77e+00, 3.77e+00 , &
                          2.76e+00, 3.04e+00, 3.14e+00, 3.20e+00, 3.23e+00, 3.25e+00, 3.25e+00 , &
                          2.95e+00, 2.75e+00, 3.03e+00, 3.15e+00, 3.34e+00, 3.48e+00, 3.48e+00 , &
                          2.62e+00, 2.71e+00, 2.65e+00, 2.80e+00, 2.80e+00, 2.97e+00, 2.97e+00 , &
                          2.26e+00, 2.59e+00, 2.33e+00, 2.43e+00, 2.62e+00, 3.01e+00, 3.01e+00 , &
                          1.94e+00, 1.29e+00, 1.65e+00, 1.65e+00, 1.88e+00, 6.49e-01, 6.49e-01 , &
                         -2.33e+00,-1.83e+00, 4.17e-01,-2.67e+00,-7.20e-01, 2.34e-01, 2.34e-01 /)

			E_win_over =  reshape(temp, (/ 7,8 /))

		
		
		

			A_nir_over_water = (/ 1.98e+01,-3.83e+00,-1.41e+01,-1.49e+01,-3.62e+00, &
                        -1.52e+01, 1.24e+01, 3.23e+00,-5.14e+00,-2.54e-01,&
                        -1.42e+00,-4.17e+00, 4.11e+00, 9.23e+00, 6.03e+00, &
                         3.11e+00, 3.11e+00, 3.11e+00 /)
			   
			B_nir_over_water = (/-1.44e+01, 3.48e+00, 1.27e+01, 1.37e+01, 4.81e+00,&
                         1.30e+01,-1.22e+01,-2.07e+00, 2.57e+00, 2.52e+00,&
                         4.71e+00, 4.44e-01,-3.24e+00,-1.10e+01,-6.34e+00,&
                        -3.65e+00,-3.65e+00,-3.65e+00 /)
			  
			C_nir_over_water = (/ 2.64e+00,-1.87e+00,-4.24e+00,-4.98e+00,-3.13e+00,&
                        -4.11e+00, 2.64e+00,-1.19e+00,-6.22e-01,-3.04e+00,&
                        -4.57e+00, 1.19e+00,-6.76e-01, 3.30e+00, 7.73e-01,&
                         9.53e-01, 9.53e-01, 9.53e-01 /)
			   
			D_nir_over_water = (/ 8.73e-01, 1.18e+00, 1.45e+00, 1.43e+00, 1.36e+00,&
                         1.18e+00, 6.15e-01, 1.09e+00, 6.95e-01, 1.26e+00,&
                         1.74e+00,-1.17e-01, 9.93e-01, 2.53e-01, 9.02e-01,&
                         2.00e-01, 2.00e-01, 2.00e-01 /)
			   
			E_nir_over_water = (/ 5.34e-02, 4.15e-02, 3.46e-02, 2.77e-02, 3.93e-02,&
                         6.08e-02, 7.28e-02, 5.06e-02, 7.56e-02, 3.34e-02,&
                         1.97e-02, 1.92e-01, 5.29e-02, 1.07e-01, 7.66e-02,&
                         3.03e-01, 3.03e-01, 3.03e-01 /)


		
			A_nir_over_land = (/ -1.94e+00, 9.09e-01, 2.51e+00,-1.01e+01, 1.01e+01, &
                        -7.71e-01, 1.31e+00, 5.42e+00,-2.80e-01, 7.54e-01,&
                         1.85e+00,-1.02e+01, 2.79e+00, 2.83e+00, 9.51e-01,&
                        -3.39e+00,-3.39e+00,-3.39e+00 /)
						
			B_nir_over_land =  (/ 2.55e+00,-9.93e-01,-2.95e+00, 7.17e+00,-9.65e+00, &
                        -2.59e+00,-2.52e+00,-7.57e+00,-1.93e+00,-1.93e+00,&
                        -2.33e+00, 7.87e+00,-3.29e+00,-3.75e+00,-5.92e-01,&
                         4.13e+00, 4.13e+00, 4.13e+00 /)

			C_nir_over_land = (/ -1.24e+00, 1.82e-01, 1.15e+00,-5.48e-01, 3.33e+00, &
                         3.19e+00, 1.93e+00, 3.63e+00, 2.45e+00, 1.87e+00, &
                         1.33e+00,-6.52e-01, 1.53e+00, 1.79e+00,-3.28e-02, &
                        -1.26e+00,-1.26e+00,-1.26e+00 /)
		
			D_nir_over_land = (/  2.33e-01,-4.16e-03,-1.74e-01,-5.56e-01,-5.56e-01, &
                        -1.13e+00,-7.03e-01,-7.78e-01,-1.01e+00,-7.98e-01, &
                        -5.12e-01,-7.32e-01,-4.81e-01,-4.57e-01,-1.50e-02, &
                        -1.05e-01,-1.05e-01,-1.05e-01 /)

			E_nir_over_land = (/  3.16e-01, 3.17e-01, 3.34e-01, 3.27e-01, 3.29e-01, &
                         3.43e-01, 3.12e-01, 3.19e-01, 3.46e-01, 3.25e-01, &
                         3.14e-01, 3.98e-01, 3.17e-01, 3.22e-01, 3.08e-01, &
                         4.77e-01, 4.77e-01, 4.77e-01 /)

	
		end subroutine init_coeffs

!=====================================
	function poly4(a, b, c, d, e, x)
!=====================================

!-----------------------------------------------------------------------
! !f90
!
! !Description:
!     This function evaluates a 4th degree polynomial
!
! !Input Parameters:
!		real*8 -- a, b, c, d, e -- coefficients of the poly
!       real -- x -- the abscissa to evaluate the poly at
! !Output Parameters:
!		y = ax^4 + bx^3 + cx^2 + dx + e
!
! !Revision History:
!     2009/03/20 wind: Gala Wind
!     Revision 1.0 Initial Revision
!
! !Team-Unique Header:
!    Developed by the  Cloud Retrieval Group, NASA GSFC, Greenbelt, Maryland, USA.
!
! !References and Credits:
!      Written by:
!      Gala Wind
!      SSAI
!      Code 613.2, NASA/GSFC
!      Greenbelt, MD 20771
!      Gala.Wind@nasa.gov
!
! !Design Notes:
!   NONE
!
! !END
!
!-----------------------------------------------------------------------
	
		real*8, intent(in) :: a, b, c, d, e
		real, intent(in) ::x
		real :: poly4
	
		poly4 =  a*x**4 + b*x**3 + c*x**2 + d*x + e
	
	end function poly4
	

!=====================================

	function PH_multilayer (sza, vza, rel_az, measurements, cloud_info, instrument, lat)

!=====================================

!-----------------------------------------------------------------------
! !f90
!
! !Description:
!     This function evaluates a 4th degree polynomial
!
! !Input Parameters:
!		real*4 -- sza -- solar zenith angle
!       real*4 -- vza -- sensor zenith angle
!       real*4 -- rel_az -- relative azimuth
!       real*4, dimension(:) -- measurements -- array of reflectances and radiances for one pixel
!       processflag -- cloud_info -- info about a pixel
!       character(*) -- instrument -- name of the platform (Terra or Aqua)
!       real*4 -- lat -- pixel latitude 
! !Output Parameters:
!		multilayer cloud result
!
! !Revision History:
!     2009/03/20 wind: Gala Wind
!     Revision 1.0 Initial Revision
!
! !Team-Unique Header:
!    Developed by the  Cloud Retrieval Group, NASA GSFC, Greenbelt, Maryland, USA.
!
! !References and Credits:
!      Written by:
!      Gala Wind
!      SSAI
!      Code 613.2, NASA/GSFC
!      Greenbelt, MD 20771
!      Gala.Wind@nasa.gov
!
! !Design Notes:
!
!   Pavolonis, M.J. Heidinger A.K. "Daytime Cloud Overlap Detection from AVHRR and VIIRS" 
!   Journal of Applied Meteorology, 2004, Vol.43, pp. 762-778
! 
!   the algorithm has been modified from the one described in this paper to include the 
!   improvements made to it since publication. 
!
! !END
!
!-----------------------------------------------------------------------
	use mod06_run_settings
	use specific_other


	real, intent(in) :: lat, sza, vza, rel_az
	real, dimension(:), intent(in) :: measurements
	type(processflag), intent(in) :: cloud_info
	character(len=*) :: instrument

	integer*1 :: PH_multilayer


	real :: R065, R040, R138, R163, I11, I12
	real :: dtor
	real :: BT11, BT12, BTdiff, miu0, scata
	integer :: index1, index2, index3
	real :: WIN_OVER_THRES, NIR_OVERLAP_THRES, MIN_REF26_OVER
	real, parameter :: BAD = 9999.
	logical :: OVER_AVHRR, OVER_VIIRS
	real :: MIN_REF26_OCEAN_LOW, MIN_REF26_OCEAN_HIGH, MIN_REF26_LAND_LOW, MIN_REF26_LAND_HIGH
	real :: REF26_WIN_CHECK_THRES, SNOW_REF6_THRES
	real :: REF26_WIN_CHECK_THRES_WATER, REF26_WIN_CHECK_THRES_LAND

	integer :: store1, store2, store3

	dtor = acos(-1.0)
	miu0 = cos (sza / 180. * dtor)

	scata = acos ( miu0 * cos(vza / 180. *dtor) + &
		sin( sza / 180. *dtor) * sin(vza / 180. * dtor) * sin (rel_az / 180. *dtor)) * 180. / dtor


	R065 = measurements(band_0065)
	R040 = measurements(band_0410)
	R138 = measurements(band_0138)
	R163 = measurements(band_0163)
	I11 = measurements(band_1100)
	I12 = measurements(band_1200)
	
	BT11 = modis_bright(instrument, I11, channel_11um, 1)
	BT12 = modis_bright(instrument, I12, channel_12um, 1)

	call init_coeffs

	MIN_REF26_OCEAN_HIGH = 0.1
	MIN_REF26_OCEAN_LOW = 0.025
	MIN_REF26_LAND_LOW = 0.027
	MIN_REF26_LAND_HIGH = 0.1
			
	REF26_WIN_CHECK_THRES_WATER = 0.08
	REF26_WIN_CHECK_THRES_LAND = 0.12

	store1 = vza / 10
	store2 = sza / 10
	store3 = scata / 10

	index1 = min(nvza-1, max(0, store1 ))  + 1
	index2 = min(nsza-1, max(0, store2 ) ) + 1 
	index3 = min(nsct-1, max(0, store3 ) ) + 1

	! This is how the BT11-BT12 threshold is defined. Don't ask me why.
	if (R065 >= 0.3 .and. R065 <= 0.6) then 
		WIN_OVER_THRES = max ( poly4(A_win_over(index1, index2), B_win_over(index1, index2), &
			C_win_over(index1, index2), D_win_over(index1, index2), E_win_over(index1, index2), &
			R065), MIN_win_over(index1, index2)) - 0.25
	else if (R065 > 0.6 .and. R065 <= 1.0) then 
		WIN_OVER_THRES = MIN_win_over(index1, index2) - 0.25
	else
		WIN_OVER_THRES = BAD
	endif


	! WARNING WARNING Pavolonis and Heidinger use the deep-blue band channel 8 to do additional testing
	! over desert. That is not mentioned anywhere in their paper. 
			
	call set_PH_desert(cloud_info%desert_surface, R040, WIN_OVER_THRES)
	

	! This is the 1.63um threshold definition
			
	if (R138 >= 0.4 .or. cloud_info%desert_surface) then 
		NIR_OVERLAP_THRES = BAD
		REF26_WIN_CHECK_THRES = REF26_WIN_CHECK_THRES_LAND
	else 
			
		if (cloud_info%ocean_surface) then 
		
			NIR_OVERLAP_THRES = poly4 (A_nir_over_water(index3), B_nir_over_water(index3), &
				C_nir_over_water(index3), D_nir_over_water(index3), E_nir_over_water(index3), &
				R138) + 0.03
					
			! Aqua correction also is not mentioned in the paper
			if (instrument == "Aqua") NIR_OVERLAP_THRES = NIR_OVERLAP_THRES - 0.09
			
			REF26_WIN_CHECK_THRES = REF26_WIN_CHECK_THRES_WATER
		
		else
		
			NIR_OVERLAP_THRES = max (poly4 (A_nir_over_land(index3), B_nir_over_land(index3), &
				C_nir_over_land(index3), D_nir_over_land(index3), E_nir_over_land(index3), R138), 0.25 ) + 0.05
			
			if (instrument == "Aqua") NIR_OVERLAP_THRES = NIR_OVERLAP_THRES - 0.09
			
			REF26_WIN_CHECK_THRES = REF26_WIN_CHECK_THRES_LAND
		endif 
			
	endif


			
	! THIS LOGIC IS NOT IN THE PAPER

	if (cloud_info%ocean_surface) then 
		if (lat >= 50. .or.  lat <= -50.) then 
			MIN_REF26_OVER = MIN_REF26_OCEAN_HIGH
		else
			MIN_REF26_OVER = MIN_REF26_OCEAN_LOW
		endif
	else if (cloud_info%desert_surface) then 
		if (lat >= 50. .or. lat <= -50.) then 
			MIN_REF26_OVER = MIN_REF26_LAND_HIGH
		else
			MIN_REF26_OVER = MIN_REF26_LAND_LOW
		endif
	else 
		if (lat >= 40. .or. lat <= -40.) then 
			MIN_REF26_OVER = MIN_REF26_LAND_HIGH
		else
			MIN_REF26_OVER = MIN_REF26_LAND_LOW
		endif
	endif

	if (lat >= 60. .or. lat <= -60) then 
		SNOW_REF6_THRES = 0.3
	else
		SNOW_REF6_THRES = 0.1
	endif

	BTdiff = BT11 - BT12


	PH_multilayer = 0

	if (R065 > 0.3 .and. R065 < 1.0 .and. BTdiff > WIN_OVER_THRES .and. BT11 < 270.) PH_multilayer = 1
				
	! this is the NIR reflectance method
	if (R138 < REF26_WIN_CHECK_THRES) then 
				
		if ( R138 > MIN_REF26_OVER .and. R163 > NIR_OVERLAP_THRES .and. R163/R065 < 1.0 .and. &
			BT11 < 280.0 .and. BTdiff > WIN_OVER_THRES .and. BT11 > 220.) PH_multilayer = 1
				
	else 
				
		if ( R138 > MIN_REF26_OVER .and. R163 > NIR_OVERLAP_THRES .and. R163/R065 < 1.0 .and. &
			BT11 < 280.0 .and. BT11 > 220.) PH_multilayer = 1
						
	endif
				
	! this is the IR split window technique
			
	if (BTdiff > WIN_OVER_THRES .and. BT11 < 270. .and. R163 > SNOW_REF6_THRES .and. BT11 > 220.) &
			PH_multilayer = 1




	end function PH_multilayer

	
	

 end module multi_layer_clouds

#endif
