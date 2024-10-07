module get_retrieval_uncertainty

 implicit none
 private

  logical :: print_info
  
  real :: new_water_radii(25), new_ice_radii(25)

 public ::  getUncertainties, getradiibounds, init_half_radii

 contains

subroutine init_half_radii

	use libraryarrays


	integer :: i, n
	
	n = number_waterradii
	new_water_radii(1) = water_radii(1)
	do i=2, n
		new_water_radii(i) = water_radii(i) - (water_radii(i) - water_radii(i-1)) / 2.0	
	end do
	new_water_radii(n+1) = water_radii(n)
	
	n = number_iceradii
	new_ice_radii(1) = ice_radii(1)
	do i=2, n
		new_ice_radii(i) = ice_radii(i) - (ice_radii(i)-ice_radii(i-1)) / 2.0
	end do
	new_ice_radii(n+1) = ice_radii(n)
	

end subroutine init_half_radii




subroutine getUncertainties  (thickness,               &
                              re,                      &
                              water_path,              &
                              phase,                   &
                              R1R2wavelengthIdx,       &
							  meas_error, &
                              surface_albedo,          &
                              transmittance_w1,        &
                              transmittance_w2,        &
                              delta_transmittance_w1,  &
                              delta_transmittance_w2,  &
                              transmittance_stddev_w1, &
                              transmittance_stddev_w2, &
                              emission_pw, &
                              emission_Tc, &
                              sigma_R37_pw, &
                              uTau,uRe,uWaterPath, xpoint, ypoint)
!-----------------------------------------------------------------------
! !f90
!
! !Description:
!     This subroutine computes the final uncertainties in thickness, re and water vapor
!
! !Input Parameters:
!     NONE
!
! !Output Parameters:
!     NONE
!
! !Revision History:
!     2004/06/01 bwind: Brad Wind
!     Revision 1.0 Initial Revision
!
! !Team-Unique Header:
!    Developed by the  Cloud Retrieval Group, NASA GSFC, Greenbelt, Maryland, USA.
!
! !References and Credits:
!    Written by:
!    Bradley P. Wind
!    L3 GSI
!    Code 913, NASA/GSFC
!    Greenbelt, MD 20771
!    bwind@climate.gsfc.nasa.gov
! 
!    Mark Gray
!    L3-GSI
!    Climate and Radiation Branch, Code 913
!    NASA/GSFC
!    Greenbelt MD 20771
!    gray@climate.gsfc.nasa.gov
!
! !Design Notes:
!   NONE
!
! !END
!
!-----------------------------------------------------------------------

! 
! S. Platnick, 2-23-05: 
!
!  Wrote variance calculations in terms of matrix formulation
!  Added covariance(tau,re) to water path uncertainty
!

  use GeneralAuxType
  use science_parameters
  use mod06_run_settings
  use ch_xfr, only: c2_sensor_id, OCI_ID, OCIS_ID
#if NOSWIR
  use core_arrays, only: mean_delta_ozone, ozone_transmittance, optical_thickness_37_final, &
    optical_thickness_1621_final
#else
  use core_arrays, only: mean_delta_ozone, ozone_transmittance
#endif
  use nonscience_parameters, only : fillvalue_real

  implicit none

	integer, intent(in) :: xpoint, ypoint
  real, intent(in) :: thickness, re, water_path,            &
                              surface_albedo(2),     &
                              transmittance_w1,transmittance_w2, &
                              delta_transmittance_w1,delta_transmittance_w2, &
                              transmittance_stddev_w1, transmittance_stddev_w2, meas_error(2), &
                              emission_pw(:), emission_Tc(:), sigma_R37_pw(:)

   
  integer, intent(in)       :: R1R2wavelengthIdx(2)
  character(10), intent(in) :: phase
  real, intent(out) :: uTau, uRe,uWaterPath

  real :: meanCloudTopReflW1, meanCloudTopReflW2
  real :: partialDerivTauWrtR1AtConstR2, partialDerivTauWrtR2AtConstR1 
  real :: partialDerivReWrtR2AtConstR1 , partialDerivReWrtR1AtConstR2
  real :: meanDeltaR1trans, meanDeltaR2trans
  real :: deltaRcloudtopmean(2), density, sdev_refl(2)

  real :: correlation_trans, correlation_WP
  real :: emis_pw_val, emis_Tc_val, sigma_R37_val

! declare matrices
  real, dimension(2,2) :: S, S_transpose, retrieval_covariance, dummy,  &
  								  refl_cov_total, refl_cov_sfc_albedo, &
  								  refl_cov_trans, refl_cov_trans_PW, refl_cov_trans_table, &
  								  refl_cov_calibration, refl_cov_cm_sdev, refl_cov_ws, &
  								  refl_cov_emission_pw, refl_cov_emission_Tc, refl_cov_trans_o3

	integer :: band37_index

#ifdef AIRBORNE
	band37_index = band_0370 + 3
#else

#ifdef AHI_INST	
	band37_index = band_0370
#else
	band37_index = band_0370 - 1
#endif

#endif
	
	
	print_info = .false.

! -----------------------------------------------------------------------
! Calculate change in cloud top reflectance due to sfc albedo uncertainty
! -----------------------------------------------------------------------
	sdev_refl = 0.

	if (COX_MUNK) then


		call getctrefl_windspeeddiff(re, thickness, phase,   &
							R1R2wavelengthIdx,  &  
							deltaRcloudtopmean, meanCloudTopReflW1, meanCloudTopReflW2)
	else
	
		call getctref_albedodiff(re, thickness, phase,     &
                          R1R2wavelengthIdx,  &
                          surface_albedo,     &
                          deltaRcloudtopmean, meanCloudTopReflW1, meanCloudTopReflW2)

	endif
	
	call get_refl_windvector(re, thickness, phase, R1R2wavelengthIdx, sdev_refl)
		
! -----------------------------------------------------------------------
! Calculate sensitivity derivatives and generate the sensitivity matrix S    
! -----------------------------------------------------------------------


	call sensitivityPartialDerivatives( R1R2wavelengthIdx, &
									re, thickness, phase,  &
									surface_albedo, &
									partialDerivTauWrtR1AtConstR2, &
									partialDerivTauWrtR2AtConstR1, &
									partialDerivReWrtR1AtConstR2,  &
									partialDerivReWrtR2AtConstR1 )

									
	S(1,1) = partialDerivTauWrtR1AtConstR2
	S(1,2) = partialDerivTauWrtR2AtConstR1  
	S(2,1) = partialDerivReWrtR1AtConstR2
	S(2,2) = partialDerivReWrtR2AtConstR1

	S_transpose = Transpose(S)
	

! ----------------------------------------------------------------------
! Calculate the uncertainties attributable to each source of uncertainty 
! ----------------------------------------------------------------------

! (1) Reflectance covariance matrix due to Surface albedo uncertainty

	refl_cov_cm_sdev(1,1) = sdev_refl(1)**2
#if NOSWIR
	refl_cov_cm_sdev(2,2) = 0.0
#else
	refl_cov_cm_sdev(2,2) = sdev_refl(2)**2
#endif
	refl_cov_cm_sdev(1,2) = 0.
	refl_cov_cm_sdev(2,1) = refl_cov_cm_sdev(1,2)

	if (COX_MUNK) then 
! wind speed uncertainty
		refl_cov_ws(1,1) = deltaRcloudtopmean(1)**2
#if NOSWIR
		refl_cov_ws(2,2) = 0.0
		
#else
		refl_cov_ws(2,2) = deltaRcloudtopmean(2)**2
#endif
		refl_cov_ws(1,2) = 0.
		refl_cov_ws(2,1) = refl_cov_ws(1,2)
! wind vector uncertainty		

		refl_cov_sfc_albedo = refl_cov_ws + refl_cov_cm_sdev

	else
		refl_cov_sfc_albedo(1,1) = deltaRcloudtopmean(1)**2
#if NOSWIR
		refl_cov_sfc_albedo(2,2) = 0.0
#else
		refl_cov_sfc_albedo(2,2) = deltaRcloudtopmean(2)**2
#endif
		refl_cov_sfc_albedo(1,2) = 0.
		refl_cov_sfc_albedo(2,1) = refl_cov_sfc_albedo(1,2)

		refl_cov_sfc_albedo = refl_cov_sfc_albedo + refl_cov_cm_sdev
	endif

! (2) Reflectance covariance matrix due to atmospheric transmittance uncertainty

!  2a. Calculate covariance matrix due to percipitable water (PW) first:

!  deltaR due to PW uncertainty only:
   meanDeltaR1trans =  meanCloudTopReflW1 * abs(delta_transmittance_w1)/transmittance_w1
#if NOSWIR
   meanDeltaR2trans =  0.0
#else
   meanDeltaR2trans =  meanCloudTopReflW2 * abs(delta_transmittance_w2)/transmittance_w2
#endif

   correlation_trans = 1.0
   refl_cov_trans_PW(1,1) = meanDeltaR1trans**2
   refl_cov_trans_PW(2,2) = meanDeltaR2trans**2

   if( ( c2_sensor_id /= OCI_ID ) .AND. ( c2_sensor_id /= OCIS_ID ) )then
     if (R1R2wavelengthIdx(2) == band37_index) then 
  
       ! the sigma_r37_val is already squared. 
  
  		 call get_emission_values(re,  sigma_R37_pw, phase, sigma_R37_val)
  		 refl_cov_trans_PW(1,2) = abs(meanDeltaR1trans * sqrt (sigma_R37_val + &
         refl_cov_trans_PW(2,2) ))*correlation_trans
  		 refl_cov_trans_PW(2,1) = refl_cov_trans_PW(1,2)
  		
  		 refl_cov_trans_PW(2,2) = refl_cov_trans_PW(2,2) + sigma_R37_val
  	  else
  	    refl_cov_trans_PW(1,2) = abs(meanDeltaR1trans*meanDeltaR2trans)* &
         correlation_trans
  	    refl_cov_trans_PW(2,1) = refl_cov_trans_PW(1,2)
  	
  	  endif
   endif

                
!  2b. Calculate covariance matrix due to table uncertainty next:

!  deltaR due to table uncertainty only:
   meanDeltaR1trans =  meanCloudTopReflW1 * transmittance_stddev_w1/transmittance_w1
#if NOSWIR
   meanDeltaR2trans =  0.0
#else
   meanDeltaR2trans =  meanCloudTopReflW2 * transmittance_stddev_w2/transmittance_w2
#endif

   refl_cov_trans_table(1,1) = meanDeltaR1trans**2
#if NOSWIR
   refl_cov_trans_table(2,2) = 0.0
#else
   refl_cov_trans_table(2,2) = meanDeltaR2trans**2
#endif
   refl_cov_trans_table(1,2) = 0.
   refl_cov_trans_table(2,1) = refl_cov_trans_table(1,2)

! 2c. uncertainty due to ozone

	refl_cov_trans_o3 = 0.
	if (R1R2wavelengthIdx(1) == band_0065 .and. mean_delta_ozone /= fillvalue_real &
			.and. ozone_transmittance /= fillvalue_real ) then 
		refl_cov_trans_o3(1,1) = (meanCloudTopReflW1 * abs(mean_delta_ozone)/ozone_transmittance)**2
	endif


!  2d. Combined covariance matrix:

   refl_cov_trans = refl_cov_trans_PW + refl_cov_trans_table + refl_cov_trans_o3

! (3) Reflectance covariance matrix due to calibration uncertainty

   refl_cov_calibration(1,1) = (meas_error(1) * meanCloudTopReflW1)**2
#if NOSWIR
   refl_cov_calibration(2,2) = 0.0
#else
   refl_cov_calibration(2,2) = (meas_error(2) * meanCloudTopReflW2)**2 
#endif
! delta F0 = 0.42 taken from difference between Fontenla and Thekaekara (band averaged Aqua and Terra)

   if( ( c2_sensor_id /= OCI_ID ) .and. ( c2_sensor_id /= OCIS_ID ) )then
     if (R1R2wavelengthIdx(2) == band37_index) then 
       refl_cov_calibration(2,2) = refl_cov_calibration(2,2) + &
         (( 0.42 * meanCloudTopReflW2 ) / solar_constant_37)**2
     endif
   endif
   
   refl_cov_calibration(1,2) = 0.
   refl_cov_calibration(2,1) = refl_cov_calibration(1,2)

! Total tau & re covariance matrix

   refl_cov_total = refl_cov_sfc_albedo + refl_cov_trans + refl_cov_calibration

! ------------------------------------------------------------------------------
! Calculate retrieval covariance matrix = S * retrieval_covariance * S_transpose
! ------------------------------------------------------------------------------

   dummy = MATMUL(refl_cov_total, S_transpose)
   retrieval_covariance = MATMUL(S, dummy)
#if NOSWIR
	! Use the standard deviation of all possible CER values (calculated from the LUT CER grid),
	! uniformly distributed (i.e., each is equally possible). It is important to have a somewhat
	! reasonable CER uncertainty for calculating the water path uncertainty.
	if ( ( phase == 'water') .and. ( c2_sensor_id /= OCI_ID ) .and &
     ( c2_sensor_id /= OCIS_ID ) ) then
		retrieval_covariance(1,1) = retrieval_covariance(1,1) + &
			(optical_thickness_37_final(xpoint,ypoint))**2	! Add uncertainty due to unknown CER
		retrieval_covariance(2,2) = 64.0	! stddev ~8.0, assuming 14µm CER
	else
		retrieval_covariance(1,1) = retrieval_covariance(1,1) + &
			(optical_thickness_1621_final(xpoint,ypoint))**2	! Add uncertainty due to unknown CER
		retrieval_covariance(2,2) = 289.0	! stddev ~17.0, assuming 30µm CER
	end if
#endif


! ---------------------
! Tau, re uncertainties
! ---------------------

   uTau = SQRT( retrieval_covariance(1,1) ) 
   uRe  = SQRT( retrieval_covariance(2,2) ) 


! ----------------------
! Water path uncertainty
! ----------------------


   if ((100.*utau/thickness) >= 200. .OR. (100.*ure/re) >= 200.  ) then
        uWaterPath = 200. 


   else

      if (phase == 'water') then
         density = liquid_water_density
      else
         density = ice_water_density
      endif
      
      uWaterPath = (re**2)*(uTau**2) + (thickness**2)*(uRe**2) + (uTau**2)*(uRe**2) + &
                    2*retrieval_covariance(1,2)*uTau*uRe + retrieval_covariance(1,2)**2
      if (uWaterPath < 0.) then
        uWaterPath = 200.
      else
         uWaterPath= sqrt ( density**2 * uWaterPath )
      end if
   endif


! ----------------------
! Relative uncertainties
! ----------------------

   utau = 100.*utau/thickness
   ure  = 100.*ure/re
   uwaterpath = 100.*uwaterpath/water_path

! ------------------------------------------------------
! Limit answer to a maximum of 200% relative uncertainty
! ------------------------------------------------------

!  sep: Was told by mag that max value should be 200 (or 255) 
!       because of the scaling value used in the filespec

   if (utau > 200.  ) utau = 200.
   if (ure > 200. ) ure = 200.
   if (uWaterPath > 200.) uWaterPath = 200.
   
   

end subroutine getUncertainties

subroutine get_emission_values(re, sigma37,  phase,  sigma37_val)

	use libraryarrays, only : water_radii, ice_radii, number_waterradii, number_iceradii
	use modis_numerical_module, only : bisectionsearch, linearinterpolation
	
	implicit none

	character*(*), intent(in) :: phase
	real, intent(in) :: sigma37(:), re
	real, intent(inout) :: sigma37_val
	
	integer :: i, n, ilow, ihigh

   ilow = 0 ! WDR-UIV
   ihigh = 0 ! WDR-UIV 

   ilow = 0 ! WDR-UIV
   ihigh = 0 ! WDR-UIV
	
	if (phase == 'water') then 
		call bisectionsearch(water_radii, re, ilow, ihigh)
		sigma37_val = linearinterpolation( (/water_radii(ilow), water_radii(ihigh) /), &
						(/ sigma37(ilow), sigma37(ihigh) /), re)
	else 
		call bisectionsearch(ice_radii, re, ilow, ihigh)	
		sigma37_val = linearinterpolation( (/ice_radii(ilow), ice_radii(ihigh) /), &
						(/ sigma37(ilow), sigma37(ihigh) /), re)
	endif


end subroutine get_emission_values



subroutine getctrefl_windspeeddiff(&
                              radius,      &
                             optical_thickness, &
                             phase,             &
                             wave_index,        &
                             reflectancediff, meancloudtopreflw1, meancloudtopreflw2)


  use GeneralAuxType
  use libraryinterpolates
  use libraryarrays
  use science_parameters
  use mod06_run_settings
  use interpolate_libraries, only: interpolate_wind_speed
  implicit none

  real, intent(in)  :: radius
  integer, intent(in)       :: wave_index(2)
  real, intent(in)  :: optical_thickness
  character(10),intent(in)  :: phase
  real, intent(out) :: reflectancediff(2), meancloudtopreflw1, meancloudtopreflw2

  integer :: radii_size , radius_upperbound, radius_lowerbound, radius_index, i,j, vectorsize
  real :: reflectance_lower(2), reflectance_middle(2), reflectance_upper(2)
  real, dimension(:), allocatable :: opticalthick_vector  !(6)
  
  real, dimension(:,:,:), allocatable :: temp_refl
  integer :: status
  
  integer :: dummy
  
  integer :: TOTAL_POINTS
  
     integer :: start_time, end_time, cmax, crate


  
  TOTAL_POINTS = size(library_taus) + 1
  
  
  allocate(opticalthick_vector(TOTAL_POINTS))
  
  
  reflectance_upper = 0.
  reflectance_lower = 0.
  reflectance_middle = 0.


  if (phase == 'water') then
    call getclosestradius(number_waterradii,  &
                   water_radii,radius,  &
                   i)

  else
    call getclosestradius(number_iceradii,  &
                        ice_radii,radius,  &
                        i)
  endif


    opticalthick_vector(1) = 0.
    opticalthick_vector(2:TOTAL_POINTS) = library_taus(1:(TOTAL_POINTS-1))


	if (phase == 'water') then 
		
		
		do j =1, 2
			call nonasymptotic_calcu_cox_munk(opticalthick_vector, &
											  optical_thickness, &
											  int_reflectance_water(:,wave_index(j),i), &
											  reflectance_middle(j))
		end do
	
		
		allocate(temp_refl(TOTAL_POINTS, number_wavelengths, number_waterradii))

		call interpolate_wind_speed(lastinterp_wind_speed*(1.0+wind_speed_error), &
										int_reflectance_water_wspeed, &
										temp_refl(:,:,:) )


		

		do j =1, 2
			call nonasymptotic_calcu_cox_munk(opticalthick_vector, &
											  optical_thickness, &
											  temp_refl(:,wave_index(j),i), &
											  reflectance_upper(j))
		end do

		call interpolate_wind_speed(lastinterp_wind_speed*(1.0-wind_speed_error), &
										int_reflectance_water_wspeed, &
										temp_refl(:,:,:) )


		do j =1, 2
			call nonasymptotic_calcu_cox_munk(opticalthick_vector, &
											  optical_thickness, &
											  temp_refl(:,wave_index(j),i), &
											  reflectance_lower(j))
		end do


		deallocate(temp_refl)

		do j=1, 2

			reflectancediff(j) = (abs(reflectance_lower(j) - reflectance_middle(j)) +  &
                         abs(reflectance_upper(j) - reflectance_middle(j)))/2.     
			if(j .eq. 1) then
				meancloudtopreflw1 = reflectance_middle(j)
			else
				meancloudtopreflw2 = reflectance_middle(j)
			endif


		end do

	else
	
		do j =1, 2
			call nonasymptotic_calcu_cox_munk(opticalthick_vector, &
											  optical_thickness, &
											  int_reflectance_ice(:,wave_index(j),i), &
											  reflectance_middle(j))
		end do

	
		allocate(temp_refl(TOTAL_POINTS, number_wavelengths, number_iceradii))

		call interpolate_wind_speed(lastinterp_wind_speed*(1.0+wind_speed_error), &
										int_reflectance_ice_wspeed, &
										temp_refl(:,:,:) )


		do j =1, 2
			call nonasymptotic_calcu_cox_munk(opticalthick_vector, &
											  optical_thickness, &
											  temp_refl(:,wave_index(j),i), &
											  reflectance_upper(j))
		end do

		call interpolate_wind_speed(lastinterp_wind_speed*(1.0-wind_speed_error), &
										int_reflectance_ice_wspeed, &
										temp_refl(:,:,:) )


		do j =1, 2
			call nonasymptotic_calcu_cox_munk(opticalthick_vector, &
											  optical_thickness, &
											  temp_refl(:,wave_index(j),i), &
											  reflectance_lower(j))
		end do


		deallocate(temp_refl)

		do j=1, 2

			reflectancediff(j) = (abs(reflectance_lower(j) - reflectance_middle(j)) +  &
                         abs(reflectance_upper(j) - reflectance_middle(j)))/2.     
			if(j .eq. 1) then
				meancloudtopreflw1 = reflectance_middle(j)
			else
				meancloudtopreflw2 = reflectance_middle(j)
			endif


		end do
	
	
	endif

   
   deallocate(opticalthick_vector)



end subroutine getctrefl_windspeeddiff

subroutine get_refl_windvector(radius, optical_thickness, phase, wave_index, sdev_refl)

	
  use GeneralAuxType
  use libraryinterpolates
  use libraryarrays
  use science_parameters
  use mod06_run_settings
  implicit none

  real, intent(in)  :: radius
  integer, intent(in)       :: wave_index(2)
  real, intent(in)  :: optical_thickness
  character(10),intent(in)  :: phase
  real, intent(out) :: sdev_refl(2)

  integer :: radii_size , radius_upperbound, radius_lowerbound, radius_index, i,j, vectorsize
  real, dimension(:), allocatable :: opticalthick_vector  !(6)
  
  integer :: status
  
  integer :: dummy
  
  integer :: TOTAL_POINTS
  
  TOTAL_POINTS = size(library_taus) + 1
  
  allocate(opticalthick_vector(TOTAL_POINTS))
  
  

  if (phase == 'water') then
    call getclosestradius(number_waterradii,  &
                   water_radii,radius,  &
                   i)

  else
    call getclosestradius(number_iceradii,  &
                        ice_radii,radius,  &
                        i)
  endif


    opticalthick_vector(1) = 0.
    opticalthick_vector(2:TOTAL_POINTS) = library_taus(1:(TOTAL_POINTS-1))


	if (phase == 'water') then 
		
		do j =1, 2
			call nonasymptotic_calcu_cox_munk(opticalthick_vector, &
											  optical_thickness, &
											  int_reflectance_water_sdev(:,wave_index(j),i), &
											  sdev_refl(j))
		end do

	else
	
		do j =1, 2
			call nonasymptotic_calcu_cox_munk(opticalthick_vector, &
											  optical_thickness, &
											  int_reflectance_ice_sdev(:,wave_index(j),i), &
											  sdev_refl(j))
		end do
	
	endif


	deallocate(opticalthick_vector)

end subroutine get_refl_windvector



subroutine getctref_albedodiff(&
                              radius,      &
                             optical_thickness, &
                             phase,             &
                             wave_index,        &
                             surface_albedo,    &
                             reflectancediff, meancloudtopreflw1, meancloudtopreflw2)
  use GeneralAuxType
  use libraryinterpolates
  use libraryarrays
  use science_parameters
  use mod06_run_settings
  implicit none

  real, intent(in)  :: radius
  integer, intent(in)       :: wave_index(2)
  real, intent(in)  :: optical_thickness, surface_albedo(2)
  character(10),intent(in)  :: phase
  real, intent(out) :: reflectancediff(2), meancloudtopreflw1, meancloudtopreflw2

  integer :: radii_size , radius_upperbound, radius_lowerbound, radius_index, i,j, vectorsize
  real :: reflectance_lower, reflectance_middle, reflectance_upper
  real, dimension(:), allocatable :: opticalthick_vector  !(6)
  integer :: dummy
  
  integer :: TOTAL_POINTS
  
  TOTAL_POINTS = size(library_taus) + 1
  
  
  allocate(opticalthick_vector(TOTAL_POINTS))
  
  
  reflectance_upper = 0.
  reflectance_lower = 0.
  reflectance_middle = 0.


  if (phase == 'water') then
    call getclosestradius(number_waterradii,  &
                   water_radii,radius,  &
                   i)

  else
    call getclosestradius(number_iceradii,  &
                        ice_radii,radius,  &
                        i)
						
  endif


    opticalthick_vector(1) = 0.
    opticalthick_vector(2:TOTAL_POINTS) = library_taus(1:(TOTAL_POINTS-1))

  if (phase == 'water') then
   do j = 1,2


      call nonasymptotic_albedovar(opticalthick_vector, &
                              optical_thickness,                         &
                              spherical_albedo_water(:,wave_index(j),i),&
                              int_reflectance_water(:,wave_index(j),i), &
                              int_fluxdownwater_sensor(:,wave_index(j),i), &
                              int_fluxdownwater_solar(:,wave_index(j),i), & 
                              surface_albedo(j),                              &
                              reflectance_lower, &
                              reflectance_middle, &
                              reflectance_upper)
 

		reflectancediff(j) = (abs(reflectance_lower - reflectance_middle) +  &
                         abs(reflectance_upper - reflectance_middle))/2.     
		if(j .eq. 1) then
			meancloudtopreflw1 = reflectance_middle
		else
			meancloudtopreflw2 = reflectance_middle
		endif
	enddo 
  else 
   do j = 1,2
 
	         call nonasymptotic_albedovar(opticalthick_vector, &
                              optical_thickness,                         &
                              spherical_albedo_ice(:,wave_index(j),i),&
                              int_reflectance_ice(:,wave_index(j),i), &
                              int_fluxdownice_sensor(:,wave_index(j),i), &
                              int_fluxdownice_solar(:,wave_index(j),i), &
                              surface_albedo(j),                              &
                              reflectance_lower, &
                              reflectance_middle, &
                              reflectance_upper)

			reflectancediff(j) = (abs(reflectance_lower - reflectance_middle) + &
                           abs(reflectance_upper - reflectance_middle))/2. 
     if(j .eq. 1) then
        meancloudtopreflw1 = reflectance_middle
     else
        meancloudtopreflw2 = reflectance_middle
     endif
    enddo
   endif
   
   
   deallocate(opticalthick_vector)

end subroutine getctref_albedodiff


subroutine getclosestradius(numberradii, radiidata, radius, index)

  integer, intent(in):: numberradii
  real, intent(in)   :: radiidata(numberradii), radius
  integer, intent(out):: index
  integer :: index1(1)

  index1 = MINLOC(abs(radiidata-radius))
  index=index1(1)
  
end subroutine getclosestradius



subroutine nonasymptotic_calcu_cox_munk(tau_vector,tc, rf, rf_calculated)
  use GeneralAuxType
  use modis_numerical_module
  use mod06_run_settings

  implicit none

  real, intent(in)    :: tau_vector(:),tc
  real, intent(in)    :: rf(:)
  real, intent(out)   :: rf_calculated


   integer  :: lowbound, highbound, taubounds(2)
   real     :: iftau(2), refvals(2)

   lowbound = 0 ! WDR-UIV
   highbound = 0 ! WDR-UIV

      call bisectionsearch(tau_vector, tc, lowbound, highbound)
      taubounds(1) = lowbound
      taubounds(2) = highbound
      iftau(1)    = tau_vector(taubounds(1))
      iftau(2)    = tau_vector(taubounds(2))

      refvals(1) = rf(lowbound)
      refvals(2) = rf(highbound)

      rf_calculated = linearinterpolation(iftau,refvals,tc)
	  
end subroutine nonasymptotic_calcu_cox_munk



subroutine nonasymptotic_calcu(tau_vector,tc, &
                            sfr, rf, ft1,ft0, &
                            albedo,rf_calculated)
  use GeneralAuxType
  use modis_numerical_module
  use mod06_run_settings
  use science_parameters

  implicit none

  real, intent(in)    :: tau_vector(:),tc,albedo
  real, intent(in)    :: rf(:), ft0(:), ft1(:), sfr(:)
  real, intent(out)   :: rf_calculated

  real, dimension(:), allocatable :: reflectance_vector

   integer  :: lowbound, highbound, taubounds(2)
   real     :: iftau(2), refvals(2)

  integer :: TOTAL_POINTS

  lowbound = 0 ! WDR-UIV
  highbound = 0 ! WDR-UIV
  
  TOTAL_POINTS = size(tau_vector)
 
  allocate(reflectance_vector(TOTAL_POINTS))

	reflectance_vector(1) =  albedo

	reflectance_vector(2:TOTAL_POINTS) = rf(1:(TOTAL_POINTS-1)) + &
		ft1(1:(TOTAL_POINTS-1)) * ft0(1:(TOTAL_POINTS-1)) * albedo /( 1. - albedo*sfr(1:(TOTAL_POINTS-1)))


      call bisectionsearch(tau_vector, tc, lowbound, highbound)
      taubounds(1) = lowbound
      taubounds(2) = highbound
      iftau(1)    = tau_vector(taubounds(1))
      iftau(2)    = tau_vector(taubounds(2))

      refvals(1) = reflectance_vector(lowbound)
      refvals(2) = reflectance_vector(highbound)

      rf_calculated = linearinterpolation(iftau,refvals,tc)

	deallocate(reflectance_vector)

end subroutine nonasymptotic_calcu

subroutine nonasymptotic_albedovar(tau_vector,tc, &
                            sfr, rf, ft1,ft0, &
                            albedo,rf_lower, rf_middle, rf_upper)
  use GeneralAuxType
  use science_parameters
  implicit none

  real, intent(in)    :: tau_vector(:),tc,albedo
  real, intent(in)    :: rf(:), ft0(:), ft1(:), sfr(:)
  real, intent(out)   :: rf_middle, rf_lower, rf_upper
  real :: albedoUpper


  call nonasymptotic_calcu(tau_vector,tc, &
                            sfr, rf, ft1,ft0, &
                            albedo,rf_middle)

  call nonasymptotic_calcu(tau_vector,tc, &
                            sfr, rf, ft1,ft0, &
                            (albedo * (1.0 - albedo_error)),rf_lower)

  albedoUpper = albedo * (1.0 + albedo_error)
  if(albedoUpper .gt. 0.99) albedoUpper = 0.99
  call nonasymptotic_calcu(tau_vector,tc, &
                            sfr, rf, ft1,ft0, &
                            albedoUpper,rf_upper)

end subroutine nonasymptotic_albedovar


  subroutine sensitivityPartialDerivatives( &
                                           R1R2wavelengthIdx,&
                                           re, tau,  &
                                           phase, surface_albedo, &
                                           partialDerivTauWrtR1AtConstR2,&
                                           partialDerivTauWrtR2AtConstR1,&
                                           partialDerivReWrtR1AtConstR2, &
                                           partialDerivReWrtR2AtConstR1)
!-----------------------------------------------------------------------
! !f90
!
! !Description:
!   This subroutine computes the partial derivatives of tau and re with 
!   respect to reflectance in one band as the other one is held constant
!
!
! !Input Parameters:
!
! !Output Parameters:
!   NONE
!
! !Revision History:
!   2004/05/28 bwind: Brad Wind
!   Revision 1.0 Initial Revision
!
! !Team-Unique Header:
!   Developed by the  Cloud Retrieval Group, NASA GSFC, Greenbelt, Maryland, USA.
!
! !References and Credits:
!   Written by:
!   Bradley P. Wind
!   L3 GSI
!   Code 913, NASA/GSFC
!   Greenbelt, MD 20771
!   bwind@climate.gsfc.nasa.gov
! 
!   Mark Gray
!   L3-GSI
!   Climate and Radiation Branch, Code 913
!   NASA/GSFC
!   Greenbelt MD 20771
!   gray@climate.gsfc.nasa.gov
!
! !Design Notes:
!   NONE
!
! !END
!
!-----------------------------------------------------------------------
  use GeneralAuxType
  use science_parameters
  use libraryarrays
  implicit none

  integer, intent(in)      :: R1R2wavelengthIdx(2)
  real, intent(in) :: re, tau
  real, intent(in) :: surface_albedo(2)
  character(10), intent(in):: phase
  real, intent(out) :: partialDerivTauWrtR1AtConstR2, &
                               partialDerivTauWrtR2AtConstR1, &
                               partialDerivReWrtR1AtConstR2,  &
                               partialDerivReWrtR2AtConstR1

  real :: partialDerivR1wrtTauAtConstRe, &
                  partialDerivR2wrtTauAtConstRe, &
                  partialDerivR1wrtReAtConstTau, &
                  partialDerivR2wrtReAtConstTau

  integer :: i
  integer :: wrtParam

  real :: partialDerivativesTau, partialDerivativesRe
  real :: upperLimit, lowerLimit,abscissaintervalstartingpoint
  real :: halfAbscissaInterval
  real, parameter :: halfAbscissaIntervalTau = 1., halfAbscissaIntervalRe = 1.
  real, parameter :: tauUpperLimit = 200.
  real, parameter :: tauLowerLimit = 0.

  integer, parameter :: tauParamPosition = 4
  integer, parameter :: reParamPosition = 5
  real :: small_number

! for tau, determine abscissa interval taking into account the 
! extrema and potentially rough transition point(s)

  halfAbscissaInterval = halfAbscissaIntervalTau
  abscissaIntervalStartingPoint = tau
  wrtParam = tauParamPosition


#if NOSWIR

	upperLimit=tauUpperLimit
	lowerLimit=tauLowerLimit

	i = 1
	call reflPartialDerivative( &
		upperLimit, lowerLimit, &
		halfAbscissaInterval, &
		abscissaIntervalStartingPoint, &
		re, tau,phase,  &
		surface_albedo(i), &
		wrtParam,                                &
		R1R2wavelengthIdx(i), &
		partialDerivativesTau)

	small_number = epsilon(partialDerivativesTau)

	if(abs(partialDerivativesTau) .lt. small_number) partialDerivativesTau = &
		small_number

	partialDerivTauWrtR1AtConstR2 = 1. / partialDerivativesTau
	partialDerivTauWrtR2AtConstR1 = 0.0
	partialDerivReWrtR1AtConstR2 = 0.0
	partialDerivReWrtR2AtConstR1 = 0.0

#else

! loop over wavelengths
  do i =1, 2

      upperLimit=tauUpperLimit
      lowerLimit=tauLowerLimit

!   Compute derivative for dR1/dT )re, and dR2/dT )re 
    call reflPartialDerivative( &
                               upperLimit, lowerLimit, &
                               halfAbscissaInterval, & 
                               abscissaIntervalStartingPoint, &
                               re, tau,phase,  &
                               surface_albedo(i), &
                               wrtParam,                                &
                               R1R2wavelengthIdx(i), &
                               partialDerivativesTau)
    if (i .eq. 1) then
!     dR1/dT | re
      partialDerivR1wrtTauAtConstRe = partialDerivativesTau
    else
!     dR2/dT | re
      partialDerivR2wrtTauAtConstRe = partialDerivativesTau
    endif
  end do

! for re, determine abscissa interval taking into account the extrema
  halfAbscissaInterval = halfAbscissaIntervalRe
  abscissaIntervalStartingPoint = re



! for re must determine which phase 
  if(phase == 'water') then
    upperLimit = water_radii(number_waterradii)
    lowerLimit = water_radii(1)
  else
    upperLimit = ice_radii(number_iceradii)
    lowerLimit = ice_radii(1)
  endif

! Compute derivative for dR1/dre )T, and dR2/dre)T
  wrtParam = reParamPosition

  do i = 1, 2
   call reflPartialDerivative(  &
                              upperLimit, lowerLimit,  &
                              halfAbscissaInterval,    &
                              abscissaIntervalStartingPoint, &
                              re, tau,phase,           &
                              surface_albedo(i),       &
                              wrtParam,                &
                              R1R2wavelengthIdx(i),    &
                              partialDerivativesRe)

    if (i == 1 ) then
!     dR1 / dre | T
      partialDerivR1wrtReAtConstTau = partialDerivativesRe
    else
!     dR2 / dre | T
      partialDerivR2wrtReAtConstTau = partialDerivativesRe
    endif

  enddo


! head-off divide by 0
	small_number = epsilon(partialDerivR1wrtTauAtConstRe)

  if(abs(partialDerivR1wrtTauAtConstRe) .lt. small_number) partialDerivR1wrtTauAtConstRe = &
     small_number
  if(abs(partialDerivR2wrtTauAtConstRe) .lt. small_number) partialDerivR2wrtTauAtConstRe = &
      small_number
  if(abs(partialDerivR1wrtReAtConstTau) .lt.  small_number) partialDerivR1wrtReAtConstTau = &
      small_number
  if(abs(partialDerivR2wrtReAtConstTau) .lt.  small_number) partialDerivR2wrtReAtConstTau = &
      small_number

! Compute the derivatives
  partialDerivTauWrtR1AtConstR2 = 1. / &
                                  (partialDerivR1wrtTauAtConstRe - &
                                   partialDerivR1wrtReAtConstTau * &
                                   (partialDerivR2wrtTauAtConstRe/partialDerivR2wrtReAtConstTau) )
    
  partialDerivTauWrtR2AtConstR1 = 1. / &
                                   (partialDerivR2wrtTauAtConstRe - &
                                    partialDerivR2wrtReAtConstTau * &
                                    (partialDerivR1wrtTauAtConstRe/partialDerivR1wrtReAtConstTau) )
   
  partialDerivReWrtR1AtConstR2 = 1. / &
                                  (partialDerivR1wrtReAtConstTau - &
                                   partialDerivR1wrtTauAtConstRe * &
                                   (partialDerivR2wrtReAtConstTau/partialDerivR2wrtTauAtConstRe) )
    
  partialDerivReWrtR2AtConstR1 = 1. / &
                                  (partialDerivR2wrtReAtConstTau - &
                                   partialDerivR2wrtTauAtConstRe * &
                                   (partialDerivR1wrtReAtConstTau/partialDerivR1wrtTauAtConstRe) )

#endif

if(debugPRN) then
print*, 'partialDerivTauWrtR1AtConstR2 = ', partialDerivTauWrtR1AtConstR2
print*, 'partialDerivTauWrtR2AtConstR1 = ', partialDerivTauWrtR2AtConstR1
print*, 'partialDerivReWrtR1AtConstR2 = ', partialDerivReWrtR1AtConstR2
print*, 'partialDerivReWrtR2AtConstR1 = ', partialDerivReWrtR2AtConstR1
endif

end subroutine sensitivityPartialDerivatives
  subroutine reflPartialDerivative(upperLimit, lowerLimit, &
                                   halfAbscissaInterval, abscissaIntervalStartingPoint, &
                                   re, tau, phase, &
                                   surface_albedo, &
                                   wrtParam, wavelengthIdx, &
                                   partialDerivative)
!-----------------------------------------------------------------------
! !f90
!
! !Description:
!     This subroutine computes a partial derivative of tau or re using the
!     central difference method. Special cases such as crossing to asymptotic regime
!     and being at libraries' ends are handled also. 
!
! !Input Parameters:
!     wrtParam -- integer -- flag controlling whether or not to compute the derivatives
!
! !Output Parameters:
!     partialDerivatives -- real*4, dimension(:) -- array of partial derivatives
!
! !Revision History:
!     2004/05/28 bwind: Brad Wind
!     Revision 1.0 Initial Revision
!
! !Team-Unique Header:
!    Developed by the  Cloud Retrieval Group, NASA GSFC, Greenbelt, Maryland, USA.
!
! !References and Credits:
!    Written by:
!    Bradley P. Wind
!    L3 GSI
!    Code 913, NASA/GSFC
!    Greenbelt, MD 20771
!    bwind@climate.gsfc.nasa.gov
!
!    Mark Gray
!    L3-GSI
!    Climate and Radiation Branch, Code 913
!    NASA/GSFC
!    Greenbelt MD 20771
!    gray@climate.gsfc.nasa.gov
!
! !Design Notes:
!      NONE
!
! !END
!
!-----------------------------------------------------------------------
  use GeneralAuxType
  use science_parameters
  implicit none

  integer, intent(in)       :: wrtParam, wavelengthIdx
  real, intent(in)  :: re, tau
  real, intent(in)  :: surface_albedo
  real, intent(in)  :: upperLimit, lowerLimit
  character(10),intent(in)  :: phase
  real, intent(out) :: partialDerivative

  real :: Rb, Ra
  real :: tmpTau, tmpRe
  
  real :: halfAbscissaInterval, abscissaIntervalStartingPoint
  real :: abscissaIntervalLower, abscissaIntervalUpper
  integer, parameter :: tauParamPosition = 4
  integer, parameter :: reParamPosition = 5

  integer :: idxJ
  call abscissaInterval(upperLimit, lowerLimit, &
                        halfAbscissaInterval, abscissaIntervalStartingPoint, &
                        abscissaIntervalLower, abscissaIntervalUpper)
    



  Rb = 0.
  Ra = 0.

! get reflectances
  select case (wrtParam)
    case(tauParamPosition)
!     in tau case
!     ...abscissaIntervalUpper/Lower are upper and lower optical thicknesses
      call getctref_constradius(      &
                                abscissaIntervalUpper,      &
                                abscissaIntervalLower,      &
                                re, phase,                  &
                                surface_albedo,             &
                                wavelengthIdx, Rb, Ra)
		partialDerivative = (Rb - Ra) 
		partialDerivative = partialDerivative / (abscissaIntervalUpper - abscissaIntervalLower)

    case(reParamPosition)
!     in re case
!     ...abscissaIntervalUpper/Lower are upper and lower effective radii

        call getCTref_consttau( &
                               abscissaIntervalUpper,      &
                               abscissaIntervalLower,      &
                               re, tau, phase,                  &
                               surface_albedo,             &
                               wavelengthIdx, partialDerivative)      


    case default
  end select


  ! perform the division

  end subroutine reflPartialDerivative


subroutine abscissainterval(upperLimit, lowerLimit, &
                            halfAbscissaInterval, abscissaIntervalStartingPoint, &
                            abscissaIntervalLower, abscissaIntervalUpper)

!-----------------------------------------------------------------------
! !f90
!
! !Description:
!   This subroutine figures out the interval for use in the
!   central difference method. Special cases such as crossing to asymptotic regime
!   and being at libraries' ends are handled also. 
!
! !Input Parameters:
!   NONE
!
! !Output Parameters:
!   NONE
!
! !Revision History:
!   2004/05/28 bwind: Brad Wind
!   Revision 1.0 Initial Revision
!
! !Team-Unique Header:
!   Developed by the  Cloud Retrieval Group, NASA GSFC, Greenbelt, Maryland, USA.
!
! !References and Credits:
!   Written by:
!   Bradley P. Wind
!   L3 GSI
!   Code 913, NASA/GSFC
!   Greenbelt, MD 20771
!   bwind@climate.gsfc.nasa.gov
!
!   Mark Gray
!   L3-GSI
!   Climate and Radiation Branch, Code 913
!   NASA/GSFC
!   Greenbelt MD 20771
!   gray@climate.gsfc.nasa.gov
!
! !Design Notes:
!   NONE
!
! !END
!
!-----------------------------------------------------------------------\
  use science_parameters

  implicit none

  real, intent(in) :: upperLimit, lowerLimit, &
                              halfAbscissaInterval, abscissaIntervalStartingPoint
  real, intent(out) :: abscissaIntervalLower, abscissaIntervalUpper

  ! Reminder: halfAbscissaIntervalTau = 1., halfAbscissaIntervalRe = 1.
   
  if((abscissaIntervalStartingPoint+halfAbscissaInterval) > upperLimit) then 
      abscissaIntervalUpper = abscissaIntervalStartingPoint
  else
     abscissaIntervalUpper = (abscissaIntervalStartingPoint+halfAbscissaInterval)
  end if
   
  if((abscissaIntervalStartingPoint-halfAbscissaInterval) < lowerLimit) then 
     abscissaIntervalLower = abscissaIntervalStartingPoint 
  else
     abscissaIntervalLower = (abscissaIntervalStartingPoint-halfAbscissaInterval)
  end if

end subroutine abscissainterval
subroutine getctref_constradius( &
                                optical_thickness_upper, &
                                optical_thickness_lower, &
                                radius,                  &
                                phase,                   &
                                surface_albedo,          &
                                wave_index,              &
                                reflectance_upper,       &
                                reflectance_lower)

  use libraryinterpolates
  use libraryarrays
  use modis_numerical_module
  use mod06_run_settings
  use science_parameters
  implicit none

  real, intent(in)  :: optical_thickness_upper, optical_thickness_lower
  integer, intent(in)       :: wave_index
  real, intent(in)  :: radius, surface_albedo
  character(10),intent(in)  :: phase
  real, intent(out) :: reflectance_upper, reflectance_lower
 
  integer :: radii_size , radius_upperbound, radius_lowerbound, radius_index, i,j, points
  real, dimension(:), allocatable :: opticalthick_vector!(6)
  real :: reflectance_vector_upper(3), reflectance_vector_lower(3), yy2(3)

	integer :: TOTAL_POINTS
	
	TOTAL_POINTS = size(library_taus) + 1


	allocate(opticalthick_vector(TOTAL_POINTS))


     opticalthick_vector(1) = 0.
      opticalthick_vector(2:TOTAL_POINTS) = library_taus(1:(TOTAL_POINTS-1))
 
 
	if (phase == 'water') then
       call getclosestradius(number_waterradii,  &
                             water_radii,radius,  &
                             i)

		if (COX_MUNK) then 

				call nonasymptotic_calcu_cox_munk(opticalthick_vector, &
											  optical_thickness_lower, &
											  int_reflectance_water(:,wave_index,i), &
											  reflectance_lower)

				call nonasymptotic_calcu_cox_munk(opticalthick_vector, &
											  optical_thickness_upper, &
											  int_reflectance_water(:,wave_index,i), &
											  reflectance_upper)
		
		
		else 
			call nonasymptotic_calcu(opticalthick_vector, optical_thickness_lower,  &
                                  spherical_albedo_water(:,wave_index,i), int_reflectance_water(:,wave_index,i), &
                                  int_fluxdownwater_sensor(:,wave_index,i), int_fluxdownwater_solar(:,wave_index,i), & 
                                  surface_albedo, reflectance_lower)
			call nonasymptotic_calcu(opticalthick_vector, optical_thickness_upper, &
                                  spherical_albedo_water(:,wave_index,i), int_reflectance_water(:,wave_index,i), &
                                  int_fluxdownwater_sensor(:,wave_index,i), int_fluxdownwater_solar(:,wave_index,i), & 
								  surface_albedo, reflectance_upper)
		endif

     else
       call getclosestradius(number_iceradii,  &
                             ice_radii,radius,  &
                             i)
							 
		 if (COX_MUNK) then 

				call nonasymptotic_calcu_cox_munk(opticalthick_vector, &
											  optical_thickness_lower, &
											  int_reflectance_ice(:,wave_index,i), &
											  reflectance_lower)

				call nonasymptotic_calcu_cox_munk(opticalthick_vector, &
											  optical_thickness_upper, &
											  int_reflectance_ice(:,wave_index,i), &
											  reflectance_upper)
										 
		 else 
			call nonasymptotic_calcu(opticalthick_vector, optical_thickness_lower,  &
                                  spherical_albedo_ice(:,wave_index,i), int_reflectance_ice(:,wave_index,i), &
                                  int_fluxdownice_sensor(:,wave_index,i), int_fluxdownice_solar(:,wave_index,i), & 
                                  surface_albedo, reflectance_lower)
			call nonasymptotic_calcu(opticalthick_vector, optical_thickness_upper,  &
                                  spherical_albedo_ice(:,wave_index,i), int_reflectance_ice(:,wave_index,i), &
                                  int_fluxdownice_sensor(:,wave_index,i), int_fluxdownice_solar(:,wave_index,i), & 
									surface_albedo, reflectance_upper)
		 endif
     endif

  
  deallocate(opticalthick_vector)
  
end subroutine getctref_constradius
subroutine getctref_consttau(&
                             radius_upper,      &
                             radius_lower,      &
                             radius, optical_thickness, &
                             phase,             &
                             surface_albedo,    &
                             wave_index,        &
                             partial_derivative)
  use GeneralAuxType
  use libraryinterpolates
  use libraryarrays
  use modis_numerical_module
  use mod06_run_settings
  use science_parameters
  implicit none

  real, intent(inout)  :: radius_upper, radius_lower
  integer, intent(in)       :: wave_index
  real, intent(in)  :: optical_thickness, radius, surface_albedo
  character(10),intent(in)  :: phase
  real, intent(out) :: partial_derivative
  
  real :: reflectance_upper, reflectance_lower

  integer :: radii_size , radius_upperbound, radius_lowerbound, radius_index, i,j, vectorsize
  real, dimension(:), allocatable :: opticalthick_vector, opticalthick_vector2!(6)
  real,dimension(20):: reflectance_vector
  
  integer :: dummy, il, jl, iu, ju, k
  real :: my_radius, myrefl(90), myres(90), myx(3), myy(3), half_rad(2)
  integer :: TOTAL_POINTS
  real :: mymid1, mymid2, reflectance_middle, delre
  real :: new_derivatives(25)
  integer :: istart, jstart
  
  TOTAL_POINTS = size(library_taus) + 1
  
  
  allocate(opticalthick_vector(TOTAL_POINTS), opticalthick_vector2(TOTAL_POINTS))
  

     opticalthick_vector(1) = 0.
      opticalthick_vector(2:TOTAL_POINTS) = library_taus(1:(TOTAL_POINTS-1))
     opticalthick_vector2 = opticalthick_vector
     if (phase == 'water') then

		call getradiibounds(number_waterradii, water_radii, radius, j, i)

!       radius_lower = water_radii(i)
!       radius_upper = water_radii(j)
!       call nonasymptotic_calcu(opticalthick_vector, optical_thickness,  &
!                                spherical_albedo_water(:,wave_index,i), int_reflectance_water(:,wave_index,i), &
!                                int_fluxdownwater_sensor(:,wave_index,i), int_fluxdownwater_solar(:,wave_index,i), & 
!                                surface_albedo, reflectance_lower)

!       call nonasymptotic_calcu(opticalthick_vector2, optical_thickness,  &
!                                spherical_albedo_water(:,wave_index,j), int_reflectance_water(:,wave_index,j), &
!                                int_fluxdownwater_sensor(:,wave_index,j), int_fluxdownwater_solar(:,wave_index,j), & 
!                                surface_albedo, reflectance_upper)

		new_derivatives = 0.
		
		i = i-1
		j = j+1
		if (i<1) i = 1
		if (j>number_waterradii) j = number_waterradii
		
		if (COX_MUNK) then 
		
			do il = i, j-1
				call nonasymptotic_calcu_cox_munk(opticalthick_vector, &
											  optical_thickness, &
											  int_reflectance_water(:,wave_index,il), &
											  reflectance_lower)

				call nonasymptotic_calcu_cox_munk(opticalthick_vector, &
											  optical_thickness, &
											  int_reflectance_water(:,wave_index,il+1), &
											  reflectance_upper)

				new_derivatives(il+1) = (reflectance_upper - reflectance_lower) / &
											(water_radii(il+1) - water_radii(il))
			end do
		
		else 
		
			do il = i, j-1
				call nonasymptotic_calcu(opticalthick_vector, optical_thickness,  &
								spherical_albedo_water(:,wave_index,il), int_reflectance_water(:,wave_index,il), &
                                int_fluxdownwater_sensor(:,wave_index,il), int_fluxdownwater_solar(:,wave_index,il), & 
                                surface_albedo, reflectance_lower)

				call nonasymptotic_calcu(opticalthick_vector2, optical_thickness,  &
                                spherical_albedo_water(:,wave_index,il+1), int_reflectance_water(:,wave_index,il+1), &
                                int_fluxdownwater_sensor(:,wave_index,il+1), int_fluxdownwater_solar(:,wave_index,il+1), & 
                                surface_albedo, reflectance_upper)
								
				new_derivatives(il+1) = (reflectance_upper - reflectance_lower) / &
											(water_radii(il+1) - water_radii(il))
																			
			end do 
		endif


		new_derivatives(1) = new_derivatives(2)
		new_derivatives(number_waterradii+1) = new_derivatives(number_waterradii)

		call getradiibounds(number_waterradii+1, new_water_radii, radius, iu, il)
			
		partial_derivative = linearinterpolation( (/ new_water_radii(il), new_water_radii(iu) /), &
										(/ new_derivatives(il), new_derivatives(iu) /), &
										   radius)


     else
	 
		call getradiibounds(number_iceradii, ice_radii, radius, j, i)


		new_derivatives = 0.
		
		i  = i-1
		j  = j+1
		if (i <1) i = 1
		if (j>number_iceradii) j = number_iceradii
		
		if (COX_MUNK) then 
			do il = i, j-1
				call nonasymptotic_calcu_cox_munk(opticalthick_vector, &
											  optical_thickness, &
											  int_reflectance_ice(:,wave_index,il), &
											  reflectance_lower)

				call nonasymptotic_calcu_cox_munk(opticalthick_vector, &
											  optical_thickness, &
											  int_reflectance_ice(:,wave_index,il+1), &
											  reflectance_upper)
								
				new_derivatives(il+1) = (reflectance_upper - reflectance_lower) / &
											(ice_radii(il+1) - ice_radii(il))
																			
			end do 
		
		else 
			do il = i, j-1
				call nonasymptotic_calcu(opticalthick_vector, optical_thickness,  &
                                spherical_albedo_ice(:,wave_index,il), int_reflectance_ice(:,wave_index,il), &
                                int_fluxdownice_sensor(:,wave_index,il), int_fluxdownice_solar(:,wave_index,il), & 
                                surface_albedo, reflectance_lower)

				call nonasymptotic_calcu(opticalthick_vector, optical_thickness ,  &
                                spherical_albedo_ice(:,wave_index,il+1), int_reflectance_ice(:,wave_index,il+1), &
                                int_fluxdownice_sensor(:,wave_index,il+1), int_fluxdownice_solar(:,wave_index,il+1), & 
                                surface_albedo, reflectance_upper)
								
				new_derivatives(il+1) = (reflectance_upper - reflectance_lower) / &
											(ice_radii(il+1) - ice_radii(il))
																			
			end do 
		endif

		new_derivatives(1) = new_derivatives(2)
		new_derivatives(number_iceradii+1) = new_derivatives(number_iceradii)


		call getradiibounds(number_iceradii+1, new_ice_radii, radius, iu, il)
					
		partial_derivative = linearinterpolation( (/ new_ice_radii(il), new_ice_radii(iu) /), &
										(/ new_derivatives(il), new_derivatives(iu) /), &
										   radius)

     endif


  deallocate(opticalthick_vector, opticalthick_vector2)


end subroutine getctref_consttau

subroutine getradiibounds(radiisize, radiidata, radius, upperbound, lowerbound)
  use GeneralAuxType
  implicit none

  integer,      intent(in)  :: radiisize
  real, intent(in)  :: radiidata(:), radius
  integer,      intent(out) :: upperbound, lowerbound

  integer :: i

  upperbound = -1
  lowerbound = -1
  do i = 1, radiisize - 1
    if (radius >= radiidata(i).and. radius <= radiidata(i+1)) then
      upperbound = i+1
      lowerbound = i
     exit
    endif
  enddo


end subroutine getradiibounds

end module get_retrieval_uncertainty 
