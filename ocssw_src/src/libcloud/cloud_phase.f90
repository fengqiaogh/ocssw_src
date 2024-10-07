module cloud_phase

 implicit none

contains

subroutine clouddecision(platform_name,                 &
                         cloudmask,                     &
                         measurements,                  &
                         RSSLiq,                        &
                         RSSIce,                        &
						 optical_thickness_liquid,      & 						 
						 optical_thickness_ice,         & 						 
						 effective_radius_16_liquid,	& 
						 effective_radius_21_liquid,	& 
						 effective_radius_37_liquid,	& 
						 effective_radius_16_ice,	    & 
						 effective_radius_21_ice,	    & 
						 effective_radius_37_ice,	    & 					 
                         cloud_top_temperature_1km,     &
                         cloud_mask_SPI,         	    &
                         cloudHeightMethod,             &
                         ir_phase,                      &
                         procflag_band_used_ot,         &
                         cloudsummary,                  &
                         debug,status, i, j)


  !  clouddecision will use cloud information values and the results
  !  from the CER/COT derivation for both ice and water phase clouds to
  !  decide the phase of the cloud at this pixel.  Essentially, a score is
  !  constructed from all the information: if > 0 the phase is water,
  !  < 0 ice phase and = 0 unknown phase
  !
  !  WDR added the use of c2_cmp_there to avoid certain tests (ie, if 3.7 um 
  !  data is not in the measurement suite)
  !
   ! input argument is the full cloudmask and level 1B reflectances
   ! output argument is the cloudsummary.  This reports on the underlying
   ! ecosystem and the cloud phase type

   use modis_cloudstructure
   use mod06_run_settings
   use GeneralAuxType 
   use ch_xfr, only: c2_cmp_there, c2_sensor_id, OCI_ID, OCIS_ID

   implicit none
   
   character*(*), intent(in)            :: platform_name
   logical, intent(in)                  :: debug 
   type(cloudmask_type),intent(in)      :: cloudmask
   type(processflag)  , intent(inout)     :: cloudsummary 
   type(cloudphase), intent(in)         :: ir_phase
   real, intent(in)                     :: measurements(:), cloud_top_temperature_1km, cloud_mask_SPI
   real, intent(in)                     :: RSSLiq(4), RSSIce(4)
   integer, intent(out)                 :: status
   integer, intent(in) :: i, j 
   integer*1, intent(in)                :: procflag_band_used_ot,cloudHeightMethod

   real                                 :: band_8_11_difference, xdimension,  &
                                           threshold_1, threshold_2, threshold_3, &
                                           threshold_4, bright_temp_11, bright_temp_85, &
                                           water_particle_threshold, ice_particle_threshold
   logical                              :: decision_made

   integer                              :: modis_c6_phase
   real                                 :: Re_threshold_01(3), Re_threshold_02(3)
   real                                 :: Tab_Water_Cloud_Effective_Radius(3), Tab_Ice_Cloud_Effective_Radius(3)

   real                                 :: optical_thickness_liquid,  optical_thickness_ice
   real                                 :: effective_radius_16_liquid, effective_radius_21_liquid 
   real                                 :: effective_radius_37_liquid, effective_radius_16_ice
   real                                 :: effective_radius_21_ice, effective_radius_37_ice

   logical                              :: ice_re_16_retrieval_failed, ice_re_21_retrieval_failed, ice_re_37_retrieval_failed
   logical                              :: liq_re_16_retrieval_failed, liq_re_21_retrieval_failed, liq_re_37_retrieval_failed
   logical                              :: Optical_Thickness_Ice_failed
   logical                              :: band06_reflectance_measurement_failed

   logical                              :: ASL_16, ASL_21, ASL_37
   real                                 :: min_liq_re, max_ice_re
   

   decision_made = .false. 
   status = 0 
   
   
   !initialize outputs, just to be sure 
   cloudsummary%cloudmask_determined = .true.
   cloudsummary%cloudobserved = .false.
   cloudsummary%watercloud = .false.
   cloudsummary%icecloud = .false.
   cloudsummary%unknowncloud = .false.

   !check if we have a cloud


   if (.not. cloudmask%cloudmask_determined) then
      cloudsummary%cloudmask_determined = .false.
      decision_made = .true.
      return
   endif
   
   if (cloudmask%night == 1) then
      decision_made = .true.
      return
   endif



!---------- MODIS C6 cloud phase Decision ----------!


   if (cloudmask%confident_cloudy .or. cloudmask%probablyclear_66 ) then

     cloudsummary%cloudobserved = .true.

     modis_c6_phase = 0

!---------- Re thresholds (MODIS C6 cloud phase) ----------!


     if ( cloud_mask_SPI  < 30.0 ) then 

       Re_threshold_01(1) = 30.0 
       Re_threshold_02(1) = 20.0
  
       Re_threshold_01(2) = 30.0
       Re_threshold_02(2) = 20.0

       Re_threshold_01(3) = 25.0
       Re_threshold_02(3) = 15.0

     else

	
       Re_threshold_01(1) = 90.0 
       Re_threshold_02(1) = 20.0
  
       Re_threshold_01(2) = 90.0
       Re_threshold_02(2) = 20.0
  
       Re_threshold_01(3) = 90.0
       Re_threshold_02(3) = 15.0
     endif


     min_liq_re = 4.0 
     max_ice_re = 60.0

!--------------------------------------------------!

     Tab_Water_Cloud_Effective_Radius(1) = effective_radius_16_liquid
     Tab_Water_Cloud_Effective_Radius(2) = effective_radius_21_liquid
     Tab_Water_Cloud_Effective_Radius(3) = effective_radius_37_liquid

     Tab_Ice_Cloud_Effective_Radius(1) = effective_radius_16_ice
     Tab_Ice_Cloud_Effective_Radius(2) = effective_radius_21_ice
     Tab_Ice_Cloud_Effective_Radius(3) = effective_radius_37_ice

!---------- Initialisation ----------!

     ice_re_16_retrieval_failed = .TRUE. 
     ice_re_21_retrieval_failed = .TRUE. 
     ice_re_37_retrieval_failed = .TRUE. 

     liq_re_16_retrieval_failed = .TRUE. 
     liq_re_21_retrieval_failed = .TRUE. 
     liq_re_37_retrieval_failed = .TRUE. 

     ASL_16 = .FALSE. 
     ASL_21 = .FALSE.
     ASL_37 = .FALSE.

     Optical_Thickness_Ice_failed = .TRUE.

     band06_reflectance_measurement_failed = .TRUE.

     if ( optical_thickness_ice > 0.0 ) Optical_Thickness_Ice_failed = .FALSE.

     if ( measurements(4) > 0.0 ) band06_reflectance_measurement_failed = .FALSE.

     if ( RSSLiq(1) > 0.0 .and. RSSIce(1) > 0.0 ) ASL_16 = .TRUE.
     if ( RSSLiq(2) > 0.0 .and. RSSIce(2) > 0.0 ) ASL_21 = .TRUE.
     if ( ( c2_cmp_there(band_0370) == 1 ) .and. RSSLiq(3) > 0.0 .and. &
          RSSIce(3) > 0.0 ) ASL_37 = .TRUE.
     if( ( c2_sensor_id == OCI_ID ) .OR. ( c2_sensor_id == OCIS_ID ) ) &
       ASL_37 = .false.


     if ( RSSIce(1) < 0 ) then
       if ( Tab_Ice_Cloud_Effective_Radius(1) > 0.0 ) ice_re_16_retrieval_failed = .FALSE.
     end if

     if (  RSSIce(2) < 0 ) then
       if ( Tab_Ice_Cloud_Effective_Radius(2) > 0.0 ) ice_re_21_retrieval_failed = .FALSE.
     end if

     if ( RSSIce(3) < 0 ) then
       if ( Tab_Ice_Cloud_Effective_Radius(3) > 0.0 ) ice_re_37_retrieval_failed = .FALSE.
     end if

     if ( RSSLiq(1) < 0 ) then
       if ( Tab_Water_Cloud_Effective_Radius(1) > 0.0 ) liq_re_16_retrieval_failed = .FALSE.
     end if

     if ( RSSLiq(2) < 0 ) then
       if ( Tab_Water_Cloud_Effective_Radius(2) > 0.0 ) liq_re_21_retrieval_failed = .FALSE.
     end if

     if ( RSSLiq(3) < 0 ) then
       if ( Tab_Water_Cloud_Effective_Radius(3) > 0.0 ) liq_re_37_retrieval_failed = .FALSE.
     end if


!---------- Cloud_Phase_Infrared_1km  (cloud phase test) ----------!

     if ( ir_phase%watercloud == 1 ) modis_c6_phase = modis_c6_phase + 1
     if ( ir_phase%icecloud == 1) modis_c6_phase = modis_c6_phase - 1 

!---------- Re 1.6 (cloud phase test) ----------!

     if (procflag_band_used_ot == 3 .and. optical_thickness_liquid < 3.0d0 ) then

       if ( measurements(band_0163) > 0.0d0 .and. measurements(band_0124) > 0.0d0  ) then

         if ( measurements(band_0163)/measurements(band_0124) < 0.45 ) modis_c6_phase = modis_c6_phase - 1 
         if ( measurements(band_0163)/measurements(band_0124) > 0.70 ) modis_c6_phase = modis_c6_phase + 1

       end if

     else
       if ( ASL_16 ) then

         if ( Tab_Water_Cloud_Effective_Radius(1) < min_liq_re + 0.01 ) modis_c6_phase = modis_c6_phase + 1
         if ( Tab_Ice_Cloud_Effective_Radius(1) > max_ice_re - 0.01 ) modis_c6_phase = modis_c6_phase - 1

       else

         if ( .not. ice_re_16_retrieval_failed .and. Tab_Ice_Cloud_Effective_Radius(1) .gt. Re_threshold_01(1) ) &
           modis_c6_phase = modis_c6_phase - 1
         if ( .not. ice_re_16_retrieval_failed .and. Tab_Ice_Cloud_Effective_Radius(1) < Re_threshold_02(1) ) &
           modis_c6_phase = modis_c6_phase + 1
         if ( ice_re_16_retrieval_failed .and. .not. liq_re_16_retrieval_failed ) modis_c6_phase = modis_c6_phase + 1

       end if
     end if

!---------- Re 2.1 (cloud phase test) ----------!

     if ( procflag_band_used_ot == 3 .and. optical_thickness_liquid < 3.0d0 ) then

       if( measurements(band_0213) > 0.0d0 .and. measurements(band_0124) > 0.0d0  )then

         if( measurements(band_0163) < 0.0d0 )then
           if ( measurements(band_0213)/measurements(band_0124) < 0.20 ) modis_c6_phase = modis_c6_phase - 2
           if ( measurements(band_0213)/measurements(band_0124) > 0.45 ) modis_c6_phase = modis_c6_phase + 2
         else
           if ( measurements(band_0213)/measurements(band_0124) < 0.20 ) modis_c6_phase = modis_c6_phase - 1
           if ( measurements(band_0213)/measurements(band_0124) > 0.45 ) modis_c6_phase = modis_c6_phase + 1 
         end if

       end if
 
     else
    
       if ( ASL_21 )then

	     if ( band06_reflectance_measurement_failed ) then
	  
	  	   if ( Tab_Water_Cloud_Effective_Radius(2) < min_liq_re + 0.01 ) modis_c6_phase = modis_c6_phase + 2
           if ( Tab_Ice_Cloud_Effective_Radius(2) > max_ice_re - 0.01 ) modis_c6_phase = modis_c6_phase - 2
      
         else
      
	       if ( Tab_Water_Cloud_Effective_Radius(2) < min_liq_re + 0.01 ) modis_c6_phase = modis_c6_phase + 1
	       if ( Tab_Ice_Cloud_Effective_Radius(2) > max_ice_re - 0.01 ) modis_c6_phase = modis_c6_phase - 1
	       
         end if

       else 

         if ( .not. ice_re_21_retrieval_failed .and. Tab_Ice_Cloud_Effective_Radius(2) > Re_threshold_01(2) ) then
           if ( band06_reflectance_measurement_failed ) modis_c6_phase = modis_c6_phase - 2 ! If no info from 1.6 increase the weight of 2.1 test
           if ( .not. band06_reflectance_measurement_failed ) modis_c6_phase = modis_c6_phase - 1
         end if

         if( .not. ice_re_21_retrieval_failed .and.  Tab_Ice_Cloud_Effective_Radius(2) < Re_threshold_02(2) ) then
           if ( band06_reflectance_measurement_failed ) modis_c6_phase = modis_c6_phase + 2
           if ( .not. band06_reflectance_measurement_failed ) modis_c6_phase = modis_c6_phase + 1
         end if

         if( ice_re_21_retrieval_failed .and. .not. liq_re_21_retrieval_failed ) then
           if ( band06_reflectance_measurement_failed ) modis_c6_phase = modis_c6_phase + 2
           if ( .not. band06_reflectance_measurement_failed ) modis_c6_phase = modis_c6_phase + 1
         end if

       end if
     end if

!---------- Re 3.7 (cloud phase test) ----------!

     if (procflag_band_used_ot /= 3 .or. (procflag_band_used_ot == 3 .and. optical_thickness_liquid >= 3.0d0 )) then

       if ( ASL_37 ) then

         if ( Tab_Water_Cloud_Effective_Radius(3) < min_liq_re + 0.01 ) modis_c6_phase = modis_c6_phase + 1
         if ( Tab_Ice_Cloud_Effective_Radius(3) > max_ice_re - 0.01 ) modis_c6_phase = modis_c6_phase - 1       
       else

         if ( .not. ice_re_37_retrieval_failed .and. Tab_Ice_Cloud_Effective_Radius(3) > Re_threshold_01(3) ) &
            modis_c6_phase = modis_c6_phase - 1
         if ( .not. ice_re_37_retrieval_failed .and. Tab_Ice_Cloud_Effective_Radius(3) < Re_threshold_02(3) ) &
           modis_c6_phase = modis_c6_phase + 1
         if ( ice_re_37_retrieval_failed .and. .not. liq_re_37_retrieval_failed ) modis_c6_phase = modis_c6_phase + 1

       end if
     end if

!---------- CTT (cloud phase test) ----------!

     if ( cloud_top_temperature_1km < 240.0 ) then
       if ( cloud_top_temperature_1km > 0.0 ) then
         modis_c6_phase = modis_c6_phase - 1
       end if
     end if

     if ( optical_thickness_liquid > 2.0d0 ) then	
	   if ( cloud_top_temperature_1km >= 270.0 ) then
	     modis_c6_phase = modis_c6_phase + 20
	   end if
	   if ( cloud_top_temperature_1km > 260.0 .and. cloud_top_temperature_1km < 270.0 ) then
	     modis_c6_phase = modis_c6_phase + 3
	   end if	
     else
	   if ( cloud_top_temperature_1km >= 270.0 ) then
	     modis_c6_phase = modis_c6_phase + 3
	   end if
	   if ( cloud_top_temperature_1km > 260.0 .and. cloud_top_temperature_1km < 270.0 ) then
	     modis_c6_phase = modis_c6_phase + 2
	   end if	

     end if	

 
!---------- 1.38 Test (cloud phase test) ----------!

     if ( .not. Optical_Thickness_Ice_failed .and. optical_thickness_ice < 2.0 ) then
       if ( cloudmask%test_high_138==1 .and. cloudmask%applied_highcloud138==1) modis_c6_phase = modis_c6_phase - 1
     end if


! Cold sanity check for undetermined and water cloud

   if (modis_c6_phase ==0 .or. modis_c6_phase == 1 .and. cloud_top_temperature_1km > 0.0) then 

       if (Optical_Thickness_Ice <= 40) then 

         if (ir_phase%icecloud == 1 .and. cloud_top_temperature_1km < 240.0 .and. cloudHeightMethod < 3) &
              modis_c6_phase  = modis_c6_phase - 20
       else
       
         if (ir_phase%icecloud == 1 .or. (cloud_top_temperature_1km < 240.0 .and. cloudHeightMethod < 3)) &
              modis_c6_phase  = modis_c6_phase - 20
       endif

     endif

!---------- MODIS C6 Cloud Phase Decision ----------!

     if ( modis_c6_phase > 0 ) cloudsummary%watercloud  = .true.
     if ( modis_c6_phase < 0 ) cloudsummary%icecloud  = .true.
     if ( modis_c6_phase ==  0 ) cloudsummary%unknowncloud  = .true. 

      
 
      

!---------- MODIS C6 Cloud Phase Decision Done ----------!

   end if
         
   decision_made = .true.  ! appears to be vestigial variable.

   return





end subroutine clouddecision


end module cloud_phase
