module atmospheric_correction_module

	use specific_ancillary

implicit none

private

public :: atmospheric_correction

!... scale in p, tpw and miu dimension:

INTEGER, PARAMETER :: NumberOfPressureHeights = 10
REAL, DIMENSION(NumberOfPressureHeights), PARAMETER :: &
heightScale = (/ 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000. /)


! changed by G. Wind 5.6.04 to match the new axis of transmittance table
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
miuScale = (/ 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, &
              0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, &
              0.95, 1.00 /) 


contains

#if AMS_INST | ASTER_INST | AVIRIS_INST | RSP_INST | EPIC_INST | SSFR_INST
subroutine fascode_vals(xpoint, ypoint, P, model, emission, trans2way, trans1way, trans2way_mean, &
						trans2way_low, trans2way_high, trans1way_low, trans1way_high) 
	use global_model_grids
	
	real, intent(in) :: P
	type(ancillary_type), intent(in) :: model
	real, intent(inout) :: emission, trans2way, trans1way
	real, intent(inout) :: trans2way_mean, trans1way_low, trans1way_high, trans2way_low, trans2way_high
	integer, intent(in) :: xpoint, ypoint

end subroutine fascode_vals

#else
subroutine fascode_vals(xpoint, ypoint, P, model, emission, trans2way, trans1way, trans2way_mean, &
						trans2way_low, trans2way_high, trans1way_low, trans1way_high) 
	
	use science_parameters, only: model_levels, delta_Pc
	use rtm_support, only: rtm_rad_atm_clr, rtm_trans_2way, rtm_trans_atm_clr, &
					rtm_trans_2way_mean, rtm_trans_2way_low, rtm_trans_atm_clr_low, &
					rtm_trans_2way_high, rtm_trans_atm_clr_high
	use core_arrays, only : cloud_height_method
	use global_model_grids

	real, intent(in) :: P
	type(ancillary_type), intent(in) :: model
	real, intent(inout) :: emission, trans2way, trans1way
	real, intent(inout) :: trans2way_mean, trans1way_low, trans1way_high, trans2way_low, trans2way_high
	integer, intent(in) :: xpoint, ypoint

	integer :: i, n, ilow, ihi, istart, iend
	
	real :: elow, ehi, t2low, t2hi, t1low, t1hi

	istart = model%trop_level
	iend = model%surface_level

	do i=istart, iend
		if (P <= model%pressure_profile(i)) exit
	end do
	
	if (i > iend) i = iend
	
	emission = rtm_rad_atm_clr(2, i)
	trans2way = rtm_trans_2way(2,i) ! 2 way correction for 3.7um, to be used later
	trans1way = rtm_trans_atm_clr(2,i)  ! 1 way correction for 3.7um, to be user later

	if (cloud_height_method(xpoint, ypoint) == 6) then 
		trans2way_mean = rtm_trans_2way_mean(i) ! 2 way correction for 3.7um, to be used later

		trans2way_low = rtm_trans_2way_low(2,i)
		trans2way_high = rtm_trans_2way_high(2,i)

		trans1way_low = rtm_trans_atm_clr_low(2,i)	
		trans1way_high = rtm_trans_atm_clr_high(2,i)	

	else

		do ilow=istart, iend
			if ((P-delta_Pc) <= model%pressure_profile(ilow)) exit
		end do
		do ihi=istart, iend
			if ((P+delta_Pc) <= model%pressure_profile(ihi)) exit
		end do
		
		t2low = rtm_trans_2way(2,ilow) ! 2 way correction for 3.7um, to be used later
		t1low = rtm_trans_atm_clr(2,ilow)  ! 1 way correction for 3.7um, to be user later
		
		t2hi = rtm_trans_2way(2,ihi) ! 2 way correction for 3.7um, to be used later
		t1hi = rtm_trans_atm_clr(2,ihi)  ! 1 way correction for 3.7um, to be user later

		trans2way_mean = ( abs(t2low - trans2way) + abs(t2hi - trans2way) ) / 2.0

		trans2way_low = t2low
		trans2way_high = t2hi

		trans1way_low = t1low
		trans1way_high = t1hi

		
	endif



end subroutine fascode_vals
#endif

subroutine atmospheric_correction(xpoint,ypoint, iteration, meas_out, model, debug,status)

   use GeneralAuxType
   use science_parameters
   use core_arrays
   use libraryarrays
   use libraryinterpolates
   use mod06_run_settings
   use nonscience_parameters
   use global_model_grids
   use specific_ancillary
   use ch_xfr, only: c2_cmp_there
   implicit none

   logical, intent(in)  :: debug
   integer, intent(in)  :: xpoint,ypoint, iteration
   integer, intent(out) :: status
   real, dimension(:), intent(inout) :: meas_out
   type(ancillary_type), intent(in) :: model

   real                 :: oneway_correction(2), twoway_correction(trans_nband), &
                           oneway_standarddev(2),twoway_standarddev(trans_nband), &
                           twoway_correction_lower(trans_nband), twoway_correction_upper(trans_nband)
   real                 :: cos_solarzenith, cos_viewingangle
   integer              :: errorlevel
   integer, parameter   :: oneway_correction_index = 2, twoway_correction_index =trans_nband
   real                 :: angle
	real :: emission37, trans2way37, trans1way37
	real :: trans2way37_mean, trans2way37_low, trans2way37_high, trans1way37_low, trans1way37_high
	real :: o3_trans_low, o3_trans_high, o3_trans_550, o3_trans_470
	
	real :: sdev_2way(trans_nband), sdev_1way(2), oneway_correction_lower(2), oneway_correction_upper(2)
	
   status = 0

	call setup_atm_corr


   cos_solarzenith = cos(d2r*solar_zenith_angle(xpoint,ypoint))
  
   cos_viewingangle = cos(d2r*sensor_zenith_angle(xpoint,ypoint))

   angle = 1./cos_solarzenith + 1./cos_viewingangle

!  currently the transmittance table breaks the MODIS band ordering by
!  choosing index the wavelength without regard to the SDS ordering implicit
!  in the data
!  bottom line; in the transmittance data .94 site between .86 and 1.2
!               in the measurements .94 sits betweem 2.1 and 3.7 (because .94 is in the 1km dataset)

   if ( cloud_top_pressure(xpoint,ypoint)     >= 0. .and.  &
        cloud_top_temperature(xpoint, ypoint) >= 0. ) then

     call gettransmittancedata(oneway_correction_index, &
                           twoway_correction_index, &
                           transmit_correction_table, &
                           transmit_stddev_table, &
                           cloud_top_pressure(xpoint,ypoint), &
                           abovecloud_watervapor(xpoint,ypoint), &
                           cos_solarzenith, &
                           cos_viewingangle, &
                           oneway_correction, &
                           twoway_correction, &
                           oneway_standarddev, &
                           twoway_standarddev, &
                           errorlevel)
	if (have_fascode) then 
     call gettransmittance_simple( &
                           transmit_correction_table, &
                           cloud_top_pressure(xpoint,ypoint), &
                           abovecloud_watervapor(xpoint,ypoint)*(1. - watervapor_error), &
                           cos_solarzenith, &
                           cos_viewingangle, &
                           twoway_correction_lower, &
                           errorlevel)
     call gettransmittance_simple( &
                           transmit_correction_table, &
                           cloud_top_pressure(xpoint,ypoint), &
                           abovecloud_watervapor(xpoint,ypoint)*(1. + watervapor_error), &
                           cos_solarzenith, &
                           cos_viewingangle, &
                           twoway_correction_upper, &
                           errorlevel)
	else
     call gettransmittancedata(oneway_correction_index, &
                           twoway_correction_index, &
                           transmit_correction_table, &
                           transmit_stddev_table, &
                           cloud_top_pressure(xpoint,ypoint), &
                           abovecloud_watervapor(xpoint,ypoint)*(1. - watervapor_error), &
                           cos_solarzenith, &
                           cos_viewingangle, &
                           oneway_correction_lower, &
                           twoway_correction_lower, &
                           sdev_1way, &
                           sdev_2way, &
                           errorlevel)
     call gettransmittancedata(oneway_correction_index, &
                           twoway_correction_index, &
                           transmit_correction_table, &
                           transmit_stddev_table, &
                           cloud_top_pressure(xpoint,ypoint), &
                           abovecloud_watervapor(xpoint,ypoint)*(1. + watervapor_error), &
                           cos_solarzenith, &
                           cos_viewingangle, &
                           oneway_correction_upper, &
                           twoway_correction_upper, &
                           sdev_1way, &
                           sdev_2way, &
                           errorlevel)	
	endif

 
! the values will never be fill, because there has been a change in the gettransmittance routines 
! used above. G.Wind 2.2.05


! we need to re-order the mean

     meandelta_trans(1:trans_nband) = ( abs(twoway_correction_lower - twoway_correction) + &
                                            abs(twoway_correction_upper - twoway_correction) ) /2.
     transmittance_twoway(1:trans_nband) = twoway_correction
     transmittance_stddev(1:trans_nband) = twoway_standarddev

! we are not using the table, so transmittance for 3.7um is from FASCODE and there is no standard deviation as it's not 
! statistical. 
	if (do_37) then
		if (have_fascode) then  
			call fascode_vals(xpoint, ypoint, cloud_top_pressure(xpoint, ypoint), model, emission37, trans2way37, trans1way37, &
						trans2way37_mean, trans2way37_low, trans2way37_high, trans1way37_low, trans1way37_high )
												
			meandelta_trans(band_0370) = trans2way37_mean
			transmittance_stddev(band_0370) = 0. 
		else
			emission37 = 0.
		endif

	endif

   else

!    these values are set to have a benign effect on subsequent computations
!    status is set an an early return performed.  It is expected that subsequest
!    retrieval computation on this scene will be terminated  
     oneway_correction = 1.
     twoway_correction = 1.

     transmittance_twoway(:) = 1.
     transmittance_stddev(:) = 0.
     meandelta_trans(:) = 0.

     status = 1
	 
	 emission37 = 0.
	 trans2way37 = 1.
	 trans1way37 = 1.
	 trans2way37_mean = 1.
	 
     return
   endif
!  Beer's law for ozone transmittance from M.D. King
   ozone_transmittance = exp (-1.0 * column_ozone(xpoint, ypoint)* ozone_absorp_coeff * angle)
   
   o3_trans_low = exp (-1.0 * 0.8*column_ozone(xpoint, ypoint)* ozone_absorp_coeff * angle)
   o3_trans_high = exp (-1.0 * 1.2*column_ozone(xpoint, ypoint)* ozone_absorp_coeff * angle)
   mean_delta_ozone = (abs(o3_trans_low - ozone_transmittance) + abs(o3_trans_high - ozone_transmittance))/2.0
	
! the only thing that channels 3 and 4 have absorption from is ozone
! used for discriminating very uniform stratocumulus from aerosol when CSR=2 for the GMAO FFNET code

   	if (have_band(1)) then
	    o3_trans_470 = exp (-1.0 * column_ozone(xpoint, ypoint)* o3_coef_470 * angle)
		if (meas_out(band_0047) /= fillvalue_real)&
   		    meas_out(band_0047) = meas_out(band_0047) / o3_trans_470
	endif
   	if (have_band(2)) then
	    o3_trans_550 = exp (-1.0 * column_ozone(xpoint, ypoint)* o3_coef_550 * angle)
		if (meas_out(band_0055) /= fillvalue_real)&
	   	    meas_out(band_0055) = meas_out(band_0055) / o3_trans_550
	endif

	if (have_band(3)) then
	   if (meas_out(band_0065) /= fillvalue_real) &
    	    meas_out(band_0065) = meas_out(band_0065)/ &
                          (ozone_transmittance*twoway_correction(band_0065)) ! .68
	   meandelta_trans(band_0065) = meandelta_trans(band_0065) * ozone_transmittance
	   transmittance_twoway(band_0065) = transmittance_twoway(band_0065) * ozone_transmittance
	   transmittance_stddev(band_0065 ) = transmittance_stddev(band_0065) * ozone_transmittance
	endif

	if (have_band(4)) then
	   if (meas_out(band_0086) /= fillvalue_real) &
	        meas_out(band_0086) = meas_out(band_0086)/&
                            twoway_correction(band_0086) ! .86
	endif
	! WDR orig, out: if (have_band(5)) then
	if (c2_cmp_there(band_0935) == 1 ) then
		if (meas_out(band_0935) /= fillvalue_real) &
	        meas_out(band_0935) = meas_out(band_0935)/&
                                 twoway_correction(band_0935) ! .94
	endif
	if (have_band(6)) then
	   if (meas_out(band_0124) /= fillvalue_real) &
	       meas_out(band_0124) = meas_out(band_0124)/&
	                      twoway_correction(band_0124) ! 1.2
	endif
	if (have_band(7)) then
	   if (meas_out(band_0163) /= fillvalue_real) &
    	    meas_out(band_0163) = meas_out(band_0163)/&
                                   twoway_correction(band_0163) ! 1.6
	endif
	if (have_band(8)) then
	   if (meas_out(band_0213) /= fillvalue_real) &
	        meas_out(band_0213) = meas_out(band_0213)/&
    		                       twoway_correction(band_0213) ! 2.1
	endif

   if (have_band(9)) then
      if (meas_out(band_0226) /= fillvalue_real) &
           meas_out(band_0226) = meas_out(band_0226)/&
                                twoway_correction(band_0226) ! 2.2
	endif

	
   thermal_correction_twoway(1:2) = 1.
   thermal_correction_oneway(1:2) = 1.

!  3.7 transmittance corrections cannot be applied until we have seperated the visible/thermal contribution to 
!  the total radiance
	if (do_37) then 
! subtract the emission due to atmosphere above cloud from measured radiance
#ifndef SIM_NORAD
		meas_out(band_0370) = meas_out(band_0370) - emission37
#endif

		if (have_fascode) then 
			transmittance_twoway(band_0370) = trans2way37
		else
#if !VIIRS_INST && !AHI_INST 
		    transmittance_twoway(band_0370) = transmittance_twoway(band_0370)
			trans2way37 = twoway_correction(band_0370)
			trans1way37 = oneway_correction (1)
			trans1way37_low = oneway_correction_lower(1)
			trans1way37_high = oneway_correction_upper(1)
			trans2way37_low = twoway_correction_lower(band_0370)
			trans2way37_high = twoway_correction_upper(band_0370)
#endif
		endif
	else
		trans2way37 = 1.
		trans1way37 = 1.
		trans1way37_low = 1.
		trans1way37_high = 1.
		trans2way37_low = 1.
		trans2way37_high = 1.
	endif
	
	  thermal_correction_twoway(1) = trans2way37 ! 2 way correction for 3.7um, to be used later
	  thermal_correction_oneway(1) = trans1way37  ! 1 way correction for 3.7um, to be user later

	  thermal_correction_oneway_low(1) = trans1way37_low 
	  thermal_correction_oneway_high(1) = trans1way37_high
	  thermal_correction_twoway_low(1) = trans2way37_low
	  thermal_correction_twoway_high(1) = trans2way37_high


end subroutine atmospheric_correction

SUBROUTINE GetTransmittanceData(KDIM_1WAY,   &
                                KDIM_2WAY,   &
                                BigTauTable, &
                                BigSdevTable,&
                                p,           &
                                tpw,         &
                                miu0,        &
                                miu1,        &
                                tau1way,     &
                                tau2way,     &
                                sdev1way,    &
                                sdev2way,    &
                                errorLevel)
!................................................................................
!
!................................................................................
! !F90
!
! !DESCRIPTION:     
!     given pressure, integrated precipitable water amount, miu0 and miu1, this
!     program returns: 
!     (1) one way transmittance for 3.7 and 11 um bands.
!     (2) two way transmittance for 0.67, 0.86, 1.2, 1.6, 2.1 and 3.7 um bands.
! 
!
! !INPUT PARAMETERS:
!
!       KDIM_1WAY = delcared band-dimension for 1-way tau array; >= 2
!       KDIM_2WAY = declared band-dimension for 2-way tau array; >= 6
!
!       p(XDIM, YDIM) = pressure height in mb.
!     tpw(XDIM, YDIM) = precipitable water (g/cm2) from level p to TOA.
!    miu0(XDIM, YDIM) = cosine of solar zenith angle.
!    miu1(XDIM, YDIM) = cosine of viewing zenith at level p.
!
!
! !OUTPUT PARAMETERS:
!
!      tau1way(XDIM, YDIM, 2) = contains 1-way spectral transmittance information.
!      tau2way(XDIM, YDIM, 6) = contains 2-way spectral transmittance information.
!      sdev1way(XDIM, YDIM, 2) = contains 1-way standard deviations.
!      sdev2way(XDIM, YDIM, 6) = contains 2-way standard deviations.
! 
!      *** Note ***: Fill value for tau1way and tau2way is -1.
!
!      errorLevel = return error status
!                   0 : none. 
!                   1 : second nearest neighbour maybe used (warning).
!                   2 : cannot find required files (fatal)
! 
!
!
!
! !REQUIRED EXTERNAL PROGRAMS: 
!                 netcdf.f90 and ModisAuxType.f90
!................................................................................
! !Revision History:
!
!
! Revision 2.5 2004/06/30    MAG
! fixed confusion between transmittance table and standard deviation table
!
! Revision 2.5 2004/05/06    wind
!   added standard deviation processing in order to facilitate uncertainty calculations
! Revision 2.0  2002/10/25           wind,mag
! add bilinear smoothing to fields
!
! Revision 1.5  2001/05/01 18:50:08  jyli
! read entire transmittance table once to speed up process.
!
! Revision 1.4  2001/04/18 21:56:35  jyli
! return both 1-way (for longwave) and 2-way (for shortwave) tau.
!
! Revision 1.3  2001/04/06 19:11:57  jyli
! first delivery version which includs MODIS run environment stuff.
!
! Revision 1.1  2001/04/06 18:46:11  jyli
! Initial revision
!
!................................................................................
! !Team-unique Header:
!             Cloud Retrieval Group, NASA/GSFC
!
! !PROGRAMMER:
!
!             Mark Gray (L3 Com)
!             Climate and Radiation Branch
!             NASA Goddard Space Flight Center
!             Greenbelt, Maryland 20771, U.S.A.
!
!             Gala Wind (L3 Com)
!             Climate and Radiation Branch
!             NASA Goddard Space Flight Center
!             Greenbelt, Maryland 20771, U.S.A.
!
!             Jason Li (Emergent-IT)
!             Climate and Radiation Branch
!             NASA Goddard Space Flight Center
!             Greenbelt, Maryland 20771, U.S.A.
!*******************************************************************************
! !END
use modis_numerical_module, only: linearinterpolation
use mod06_run_settings
use specific_ancillary
use ch_xfr, only : c2_cmp_there

IMPLICIT NONE


!... i/o variable delcarations:

INTEGER, intent(in) :: KDIM_1WAY, KDIM_2WAY
REAL , intent(in)            :: p, tpw, miu0, miu1
REAL, intent(out),DIMENSION(KDIM_1WAY) :: tau1way, sdev1way
REAL, intent(out),DIMENSION(KDIM_2WAY) :: tau2way, sdev2way
INTEGER, intent(out)                               :: errorLevel
REAL, intent(in), DIMENSION(:,:,:,:)   :: bigTauTable, bigSdevTable


INTEGER :: i, pIndex, miu1Index, miu2Index, pwIndex, set_fill
INTEGER, DIMENSION(4) :: start, counts, mindex

INTEGER           :: NumberOfPixels_1km, NumberOfLines_1km

REAL            :: miu, tauValue
real            :: tau11, tau12, tau21, tau22, tau_mid1, tau_mid2
real            :: sdev11, sdev12, sdev21, sdev22, sdev_mid1, sdev_mid2
Integer         :: first_pwIndex, first_miuindex, second_pwIndex, second_miuindex
real, dimension(8) :: temp_trans, temp_sdev
!................................................................................


! initialize values
temp_trans = -999.
temp_sdev = -999.
 

errorLevel = 0

   miu = (miu0*miu1)/(miu0+miu1)

!
!... find array indicies along p, pw and miu dimension:
!
   pIndex = NINT( (p - heightScale(1)) / 100.0 ) + 1
   set_fill = 0 
   if (p< 0. .and. tpw < 0.) set_fill = 1

   if(pIndex < 1) then
      pIndex = 1
      errorLevel = 1
   elseif(pIndex > NumberOfPressureHeights) then
      pIndex = NumberOfPressureHeights
      errorLevel = 1
   endif

   miu1Index = NINT( (miu1 - miuScale(1)) / 0.05 ) + 1
   if (miu1Index < 1) then
       miu1Index = 1
       errorLevel = 1
   elseif(miu1Index > NumberOfMiu) then
       miu1Index = NumberOfMiu
       errorLevel = 1
   endif

   miu2Index = NINT( (miu - miuScale(1)) / 0.05 ) + 1
   if (miu2Index < 1) then
       miu2Index = 1
       errorLevel = 1
   elseif(miu2Index > NumberOfMiu) then
       miu2Index = NumberOfMiu
       errorLevel = 1
   endif


! changed by G. Wind 5.6.04 to match the new bins in transmittance table
   if(tpw < 0.0) then
       pwIndex = 1
       errorLevel = 1
   elseif(tpw < 0.2) then
       pwIndex = NINT( (tpw- pwScale(1)) / 0.02 ) + 1
   else
       pwIndex = NINT( (tpw - pwScale(11)) / 0.20 ) + 11
       if(pwIndex > NumberOfPrecipitableWater) then
          pwIndex = NumberOfPrecipitableWater
          errorLevel = 1
       endif
   endif

   !... 
   !... 1-way case (for longwave 3.7 and 11 micron bands):
   !... 

   mindex(1) = pIndex
   mindex(2) = pwIndex
   mindex(3) = miu1Index
   mindex(4) = 1

   tauValue = bigTauTable(mindex(1),mindex(2),mindex(3),mindex(4))

! we need to step back until we find a value that is non-fill, the second neighbor will not be
! enough in many cases

   if (tauValue < 0.0) then  ! do I need 2nd nearest neighbour?
!       errorLevel = 1
!       mindex(2) = pwIndex - 1  
      
       do i=mindex(2), 1, -1
          tauValue =  bigTauTable(mindex(1),i,mindex(3),mindex(4))
          if (tauValue >= 0.) exit
       end do

! now we have a good value of trans, non-fill, so we set the index accordingly
       mindex(2) = i

   endif

!  grab the 4 points for bilinear interpolation: pw1,u1 pw2,u1 pw1,u2, pw2,u2

   if (pwscale(mindex(2)) >= tpw ) then
      second_pwIndex = mindex(2)
      if (mindex(2) == 1) then 
         first_pwIndex = mindex(2)
      else
         first_pwIndex = mindex(2) - 1
         if (bigtautable(mindex(1), first_pwIndex, mindex(3), mindex(4)) < 0.) &
            first_pwIndex = mindex(2)
      endif

   else
      first_pwIndex = mindex(2)
      if (mindex(2) == numberofprecipitablewater) then 
         second_pwIndex = mindex(2) 
      else 
         second_pwIndex = mindex(2) + 1
         if (bigTauTable(mindex(1),second_pwIndex,mindex(3),mindex(4)) < 0.0) &
            second_pwIndex = mindex(2)
      endif
   endif

   if (miuScale(miu1index) >= miu1 ) then
      first_miuindex = miu1index -1
      second_miuindex = miu1index
   else
      first_miuindex = miu1index
      second_miuindex = miu1index + 1
   endif

   if( c2_cmp_there(band_0370) == 1) then
     do i = 1, NumberOfLongWavelengths
        mindex(4) = bandIndexMapLW(i)
  
        tau11 = BigTauTable(pindex, first_pwIndex, first_miuindex, mindex(4))
        tau12 = BigTauTable(pindex, first_pwIndex, second_miuindex, mindex(4))
  
        tau21 = BigTauTable(pindex, second_pwIndex, first_miuindex, mindex(4))
        tau22 = BigTauTable(pindex, second_pwIndex, second_miuindex, mindex(4))
  
        sdev11 = BigSdevTable(pindex, first_pwIndex, first_miuindex, mindex(4))
        sdev12 = BigSdevTable(pindex, first_pwIndex, second_miuindex, mindex(4))
  
        sdev21 = BigSdevTable(pindex, second_pwIndex, first_miuindex, mindex(4))
        sdev22 = BigSdevTable(pindex, second_pwIndex, second_miuindex, mindex(4))
  
  
  !     first we interpolate along the PW lines
        tau_mid1 = linearinterpolation( (/pwscale(first_pwIndex), pwscale(second_pwIndex) /), &
                                        (/tau11, tau21 /), &
                                        tpw)
        tau_mid2 = linearinterpolation( (/pwscale(first_pwIndex), pwscale(second_pwIndex) /), &
                                        (/tau12, tau22 /), &
                                        tpw)
  
        sdev_mid1 =linearinterpolation( (/pwscale(first_pwIndex),  pwscale(second_pwIndex)/), &
                                        (/sdev11, sdev21/), &
                                        tpw)
        sdev_mid2 =linearinterpolation( (/pwscale(first_pwIndex),  pwscale(second_pwIndex)/), &
                                        (/sdev12, sdev22/), &
                                        tpw)
  
  
  !     now we interpolate along the MIU line through tau_mid1 and tau_mid2
  !     (bilinear interpolation)
        tau1way(i) = linearinterpolation( (/miuscale(first_miuindex), miuscale(second_miuindex)/), &
                                           (/tau_mid1, tau_mid2/), &
                                           miu1)
        sdev1way(i) = linearinterpolation( (/miuscale(first_miuindex), miuscale(second_miuindex)/), &
                                           (/sdev_mid1, sdev_mid2/), &
                                           miu1) 
        if (set_fill == 1 )then
           tau1way(i) = -999.
           sdev1way(i) = -999.
        endif
     enddo
  endif

   !
   !... 2-way case (for shortwave channels upto 3.7 micron band):
   !

   mindex(2) = pwIndex
   mindex(3) = miu2Index
   mindex(4) = 1

! the pw_index values had already been determined and they don't change with changing miu, 
! so pray tell why should we recompute them. This code had been deleted. G.Wind 2.2.05

   if (miuScale(miu2index) >= miu ) then
      first_miuindex = miu2index -1
      second_miuindex = miu2index
   else
      first_miuindex = miu2index
      second_miuindex = miu2index + 1
   endif

   do i = 1, NumberOfShortWavelengths
      mindex(4) = bandIndexMapSW(i)

      tau11 = BigTauTable(pindex, first_pwIndex, first_miuindex, mindex(4))
      tau12 = BigTauTable(pindex, first_pwIndex, second_miuindex, mindex(4))

      tau21 = BigTauTable(pindex, second_pwIndex, first_miuindex, mindex(4))
      tau22 = BigTauTable(pindex, second_pwIndex, second_miuindex, mindex(4))

      sdev11 = BigSdevTable(pindex, first_pwIndex, first_miuindex, mindex(4))
      sdev12 = BigSdevTable(pindex, first_pwIndex, second_miuindex, mindex(4))

      sdev21 = BigSdevTable(pindex, second_pwIndex, first_miuindex, mindex(4))
      sdev22 = BigSdevTable(pindex, second_pwIndex, second_miuindex, mindex(4))

!     first we interpolate along the PW lines
      tau_mid1 = linearinterpolation( (/pwscale(first_pwIndex), pwscale(second_pwIndex) /), &
                                      (/tau11, tau21 /), &
                                      tpw)
      tau_mid2 = linearinterpolation( (/pwscale(first_pwIndex), pwscale(second_pwIndex) /), &
                                      (/tau12, tau22 /), &
                                      tpw)
!     then we repeat for the standard deviations
      sdev_mid1 =linearinterpolation( (/pwscale(first_pwIndex),  pwscale(second_pwIndex)/), &
                                      (/sdev11, sdev21/), &
                                      tpw)
      sdev_mid2 =linearinterpolation( (/pwscale(first_pwIndex),  pwscale(second_pwIndex)/), &
                                      (/sdev12, sdev22/), &
                                      tpw)
!     now we interpolate along the MIU line through tau_mid1 and tau_mid2
!     nothing but rather standard bilinear interpolation
!     and that is our answer.
      temp_trans(i) = linearinterpolation( (/miuscale(first_miuindex),  miuscale(second_miuindex)/), &
                                        (/tau_mid1, tau_mid2/), &
                                        miu)      

!     and these are the final standard deviations
      temp_sdev(i) = linearinterpolation( (/miuscale(first_miuindex), miuscale(second_miuindex)/), &
                                         (/sdev_mid1, sdev_mid2/), &
                                         miu )
      if (set_fill == 1 )then
         temp_trans(i) = -999.
         temp_sdev(i) = -999.
      endif
         
   end do


! remap the bands, so they can be indexed by the main instrument band index
	call remap_bands(tau2way, temp_trans, sdev2way, temp_sdev)


9999 RETURN

END SUBROUTINE GetTransmittanceData
SUBROUTINE gettransmittance_simple( &
                                   BigTauTable, &
                                   p,           &
                                   tpw,         &
                                   miu0,        &
                                   miu1,        &
                                   tau2way,     &
                                   errorLevel)
!................................................................................
!
!................................................................................
! !F90
!
! !DESCRIPTION:     
!     given pressure, integrated precipitable water amount, miu0 and miu1, this
!     program returns: 
!     (1) one way transmittance for 3.7 and 11 um bands.
!     (2) two way transmittance for 0.67, 0.86, 1.2, 1.6, 2.1 and 3.7 um bands.
! 
!
! !INPUT PARAMETERS:
!
!       KDIM_2WAY(xdim, ydim, n) = declared band-dimension for 2-way tau array; >= n bands
!
!       p(XDIM, YDIM) = pressure height in mb.
!     tpw(XDIM, YDIM) = precipitable water (g/cm2) from level p to TOA.
!    miu0(XDIM, YDIM) = cosine of solar zenith angle.
!    miu1(XDIM, YDIM) = cosine of viewing zenith at level p.
!
!
! !OUTPUT PARAMETERS:
!
!      tau2way(XDIM, YDIM, n) = contains 2-way spectral transmittance information.
! 
!      *** Note ***: Fill value for tau1way and tau2way is -1.
!
!      errorLevel = return error status
!                   0 : none. 
!                   1 : second nearest neighbour maybe used (warning).
!                   2 : cannot find required files (fatal)
! 
!................................................................................
! !Revision History:
!................................................................................
! !Team-unique Header:
!             Cloud Retrieval Group, NASA/GSFC
!
! !PROGRAMMER:
!
!             Mark Gray (L3 Com)
!             Climate and Radiation Branch
!             NASA Goddard Space Flight Center
!             Greenbelt, Maryland 20771, U.S.A.
!
!             Gala Wind (L3 Com)
!             Climate and Radiation Branch
!             NASA Goddard Space Flight Center
!             Greenbelt, Maryland 20771, U.S.A.
!
!             Jason Li (Emergent-IT)
!             Climate and Radiation Branch
!             NASA Goddard Space Flight Center
!             Greenbelt, Maryland 20771, U.S.A.
!*******************************************************************************
! !END

use modis_numerical_module, only: linearinterpolation
use mod06_run_settings
use specific_ancillary

IMPLICIT NONE



!... i/o variable delcarations:

REAL , intent(in)            :: p, tpw, miu0, miu1
REAL, intent(out),DIMENSION(:) :: tau2way 
INTEGER, intent(out)                               :: errorLevel
REAL, intent(in), DIMENSION(:,:,:,:)   :: bigTauTable



INTEGER :: i, pIndex, miu1Index, miu2Index, pwIndex, set_fill
INTEGER, DIMENSION(4) :: start, counts, mindex

INTEGER           :: NumberOfPixels_1km, NumberOfLines_1km

REAL            :: miu, tauValue
real            :: tau11, tau12, tau21, tau22, tau_mid1, tau_mid2
real            :: sdev11, sdev12, sdev21, sdev22, sdev_mid1, sdev_mid2
Integer         :: first_pwIndex, first_miuindex, second_pwIndex, second_miuindex
real, dimension(8) :: temp_trans, temp_sdev, sdev2way
!................................................................................


! initialize values
temp_trans = -999.



errorLevel = 0

   miu = (miu0*miu1)/(miu0+miu1)

!
!... find array indicies along p, pw and miu dimension:
!
   pIndex = NINT( (p - heightScale(1)) / 100.0 ) + 1
   set_fill = 0 
   if (p< 0. .and. tpw < 0.) set_fill = 1

   if(pIndex < 1) then
      pIndex = 1
      errorLevel = 1
   elseif(pIndex > NumberOfPressureHeights) then
      pIndex = NumberOfPressureHeights
      errorLevel = 1
   endif

   miu1Index = NINT( (miu1 - miuScale(1)) / 0.05 ) + 1
   if (miu1Index < 1) then
       miu1Index = 1
       errorLevel = 1
   elseif(miu1Index > NumberOfMiu) then
       miu1Index = NumberOfMiu
       errorLevel = 1
   endif

   miu2Index = NINT( (miu - miuScale(1)) / 0.05 ) + 1
   if (miu2Index < 1) then
       miu2Index = 1
       errorLevel = 1
   elseif(miu2Index > NumberOfMiu) then
       miu2Index = NumberOfMiu
       errorLevel = 1
   endif


! changed by G. Wind 5.6.04 to match the new bins in transmittance table\
   pwIndex = 1
   if(tpw < 0.0) then
       errorLevel = 1
   elseif(tpw < 0.2) then
       pwIndex = NINT( (tpw- pwScale(1)) / 0.02 ) + 1
   else
       pwIndex = NINT( (tpw - pwScale(11)) / 0.20 ) + 11
       if(pwIndex > NumberOfPrecipitableWater) then
          pwIndex = NumberOfPrecipitableWater
          errorLevel = 1
       endif
   endif

   !
   !... 2-way case (for shortwave channels upto 3.7 micron band):
   !

   mindex(1) = pIndex
   mindex(2) = pwIndex
   mindex(3) = miu2Index
   mindex(4) = 1
      
   if (pwIndex  < 1) print*, "HERE" , mindex, tpw, p
      
   tauValue =  bigTauTable(mindex(1),mindex(2),mindex(3),mindex(4))
  
! we need to step back until we find a value that is non-fill, the second neighbor will not be
! enough in many cases

   if (tauValue < 0.0) then  ! do I need 2nd nearest neighbour?
!       errorLevel = 1
!       mindex(2) = pwIndex - 1  
      
       do i=mindex(2), 1, -1
          tauValue =  bigTauTable(mindex(1),i,mindex(3),mindex(4))
          if (tauValue >= 0.) exit
       end do

! now we have a good value of trans, non-fill, so we set the index accordingly
       mindex(2) = i

   endif


!  grab the 4 points for bilinear interpolation: pw1,u1 pw2,u1 pw1,u2, pw2,u2
   if (pwscale(mindex(2)) >= tpw ) then
      second_pwIndex = mindex(2)
      if (mindex(2) == 1) then 
         first_pwIndex = mindex(2)
      else
         first_pwIndex = mindex(2) - 1
         if (bigtautable(mindex(1), first_pwIndex, mindex(3), mindex(4)) < 0.) &
            first_pwIndex = mindex(2)
      endif

   else
      first_pwIndex = mindex(2)
      if (mindex(2) == numberofprecipitablewater) then 
         second_pwIndex = mindex(2) 
      else 
         second_pwIndex = mindex(2) + 1
         if (bigTauTable(mindex(1),second_pwIndex,mindex(3),mindex(4)) < 0.0) &
            second_pwIndex = mindex(2)
      endif
   endif


   if (miuScale(miu2index) >= miu ) then
      first_miuindex = miu2index -1
      second_miuindex = miu2index
   else
      first_miuindex = miu2index
      second_miuindex = miu2index + 1
   endif

   do i = 1, NumberOfShortWavelengths
      mindex(4) = bandIndexMapSW(i)

      tau11 = BigTauTable(pindex, first_pwIndex, first_miuindex, mindex(4))
      tau12 = BigTauTable(pindex, first_pwIndex, second_miuindex, mindex(4))

      tau21 = BigTauTable(pindex, second_pwIndex, first_miuindex, mindex(4))
      tau22 = BigTauTable(pindex, second_pwIndex, second_miuindex, mindex(4))

!     first we interpolate along the PW lines
      tau_mid1 = linearinterpolation( (/pwscale(first_pwIndex), pwscale(second_pwIndex) /), &
                                      (/tau11, tau21 /), &
                                      tpw)
      tau_mid2 = linearinterpolation( (/pwscale(first_pwIndex), pwscale(second_pwIndex) /), &
                                      (/tau12, tau22 /), &
                                      tpw)

!     now we interpolate along the MIU line through tau_mid1 and tau_mid2
!     nothing but rather standard bilinear interpolation
!     and that is our answer.
      temp_trans(i) = linearinterpolation( (/miuscale(first_miuindex),  miuscale(second_miuindex)/), &
                                        (/tau_mid1, tau_mid2/), &
                                        miu)

      if (set_fill == 1 )then
         temp_trans(i) = -999.
      endif

   end do


! remap the bands, so they can be indexed by the main instrument band index
	temp_sdev(1) = -10000. ! we don't do standard deviation in this subroutine
	call remap_bands(tau2way, temp_trans, sdev2way, temp_sdev)


9999 RETURN

END SUBROUTINE gettransmittance_simple




end module atmospheric_correction_module
