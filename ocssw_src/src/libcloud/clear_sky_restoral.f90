 module clear_sky_restoral
!-------------------------------------------------------------------------------
! DESCRIPTION:
!        this module contains the routines used to perform spatial variability
!        tests and 250m cloudiness evaluation.
!
!
! PROGRAMMER:
!            Gala Wind (L3 GSI)
!            Climate and Radiation Branch
!            Code 913, NASA Goddard Space Flight Center
!            Greenbelt, Maryland 20771, U.S.A.
!
!            Mark Gray (L3 GSI)
!            Climate and Radiation Branch
!            Code 913, NASA Goddard Space Flight Center
!            Greenbelt, Maryland 20771, U.S.A.
!            gray@climate.gsfc.nasa.gov
!------------------------------------------------------------------------------
! REVISION:
!            10.22.04  integrated with MOD_PR06OD, substantial style changes.
!            10.07.04  Initial revision
!            10.12.04  added more inputs and tests to be potentially included
!------------------------------------------------------------------------------

 use GeneralAuxType
 use modis_cloudstructure
 use mod06_run_settings

 implicit none
 
 private


 public :: cloudiness_test

 contains



!====================================================================
 subroutine cloudiness_test              &
                           ( cloudmask,        &
                             process_summary,  &
                             measurement,      &
                             reflectance_box,  &
                             not_cloud,        &
                             lowvariability_confidence_test, CSR_QA, latitude, &
							 chm, vis1km_test)		! KGM 3-4-13
!====================================================================

! NOTES TO Gala: 
!
! This subroutine name should be changed to "cloudiness_test", while the above routine could be something
! like "cloudiness_test_alternate"; similar changes in modis_science_module are also needed, of course. 
! This routine has not for a long time been equivalent to the original Riedi work, but includes contributions
! from all of us and so is a misnomer. Furthermore, I have added absolute variability calculations for testing
! purposes only (presently commented out) so this routine now essentially contains the "alternate" structure anyway.

!-----------------------------------------------------------------------
! !f90
!
! !Description:
!    This subroutine performs  Clear Sky Restoral (CSR) tests 
!    in an attempt to eliminate sunglint or heavy aerosol pixels falsely flagged as
!    cloudy by the cloud mask. Examples of this include bright ocean sunglint, 
!    smoke, and dust.
!
! !Input Parameters:
!     cloudmask       -- structure          -- all cloud mask tests
!     process_summary -- structure          -- process flags, surface type etc
!     measurements    -- 4 byte real vector -- reflectances/ radiances for all bands
!     reflectance_box -- 4 byte real array  -- all reflectances for 2x2 box
!
! !Output Parameters:
!     not_cloud     -- logical -- true indicates the presence of sunglint/dust or other non-cloud situation
!     lowvariability_confidence_test - Indicates went with a "cloudy" determination but are unsure so set confidence QA to a lower value
!     CSR_QA        --
!       CSR_flag_array=0 => cloudy, all CSR tests negative
!       CSR_flag_array=1 => edge test positive (determined in subroutine remove_edge_scenes.f90)
!       CSR_flag_array=2 => spatial, etc. tests positive
!       CSR_flag_array=3 => 250 m test positive
!
! !Revision History:
!     2009/04/28, G. Wind: 
!		 changed the way the 250m tests are handled by the code. Now everything will start from being all cloudy and cleared if clear
!        250m bit is encountered. This will help with any dead 250m detectors. 
!
!     2005/09/14, S. Plantick: 
!        An error was found in the implementation of the CSR flow chart. The 250 cloudiness test
!        was not being used for clouds brighter than the reflectance threshold. This was fixed, along with:
!        (1) Additional comments and edits of existing comments.
!        (2) Local variable name changes and additional variables to implement correct 250 test logic.
!        (3) For testing/comparison purposes - the calculation of an absolute spatial variability metric.
!
!     2005, May-June: several iterations with a number of authors.
!
!     2005/02/08 wind
!        implemented Riedi algorithm as per the flow chart from Dec.10.2004
!     2004/10/23 mag
!     integrated into MOD_PR06OD
!       - added logical structures, variable types used by MOD_PR06OD
!       - cleaned up variable naming, code formatting, style.
!
!     2004/10/07 wind: Gala Wind
!     Revision 1.0 Initial revision
!
! !Team-Unique Header:
!    Developed by the  Cloud Retrieval Group, NASA GSFC, Greenbelt, Maryland, USA.
!
! !References and Credits:
!    Written by:
!    Gala Wind
!      L3 GSI
!      Code 913, NASA/GSFC
!      Greenbelt, MD 20771
!      wind@climate.gsfc.nasa.gov
!
!    Mark Gray
!      L3 GSI
!      Code 913, NASA/GSFC
!      Greenbelt, MD 20771
!      gray@climate.gsfc.nasa.gov
!
!
! !END
!
!-----------------------------------------------------------------------

  type(processflag), intent(in)              :: process_summary
  type(cloudmask_type), intent(in)           :: cloudmask
  real(single), dimension(:), intent(in)     :: measurement
  real(single), dimension(:,:,:), intent(in) :: reflectance_box
  real(single), intent(in)                   :: latitude
  integer*1, intent(in)						 :: chm		! KGM 3-4-13
  integer*1, intent(inout)                   :: CSR_QA
  logical, intent(in)						 :: vis1km_test		! KGM 3-4-13
  logical, intent(out)                       :: lowvariability_confidence_test
  logical, intent(out)                       :: not_cloud

! local variables
  logical                                    :: cloudy_250
  integer(integer_fourbyte)                  :: cloudiness_degree 
  real                                       :: average, r21_r138, r124_r21, r086_r124, &
                                                standard_deviation_val, Rbox(3,3), th_Ri, th_Ri_QA, &
                                                dispersion_val, dispersion_Ri, dispersion_Ri_QA, &
                                                Rmeas, Rth_ocean, Rth_land, Rth_desert, Rth_snowice
  integer :: rbox_x, rbox_y

  Rbox = 0 ! WDR-UIV

! Initialize outputs for no restoral, and 250m cloudiness variables
  CSR_QA = 0 ; not_cloud = .false. ; lowvariability_confidence_test = .false.
  cloudiness_degree = 16 ; cloudy_250 = .true.




! Con't run this code if there is no cloud observed.
  if (.not. process_summary%cloudobserved) return

  r21_r138  = measurement(band_0213) - measurement(band_0138)
  r124_r21  = measurement(band_0124) - measurement(band_0213)
  r086_r124 = measurement(band_0086) - measurement(band_0124)


! For use in Part V, compute the degree of cloudiness from the 250m cloud mask bytes for this pixel
!  of course we need to check (for each test) that the test was actually run
!
! The cloudiness degree is a number between 0 and 16, representing the number of
!  cloudy pixels in the 250m box corresponding to a particular 1km pixel. The box is 4x4, so
!  if the 1km pixel is completely overcaset, the degree of cloudiness would be 16.

  If (process_summary%ocean_surface) then


    if (cloudmask%visible_cloudtest_250m <= 8) then
         cloudy_250 = .false.
    endif

  End if

! TO TURN OFF 250m TEST, UNCOMMENT THE FOLLOWING LINE:
! cloudy_250 = .true.

! ENTER PART I: PERFORM SPATIAL VARIABLILITY TESTS 

! original:
! Rth_ocean   = 0.35; Rth_land = 0.45; Rth_desert = NA;    Rth_snowice = 0.65

! threshold test 1: USED FOR COLLECTION 5 OPERATIONAL VERSION
  Rth_ocean   = 0.65; Rth_land = 0.65; Rth_desert = 0.65; Rth_snowice = 0.65

! threshold test 2:
! Rth_ocean   = 0.50; Rth_land = 0.50; Rth_desert = 0.55; Rth_snowice = 0.65


  If ( &
       (  process_summary%land_surface    .and. (measurement(band_0065) > Rth_land)    ) .OR. & 
       (  process_summary%desert_surface  .and. (measurement(band_0065) > Rth_desert)  ) .OR. &
       (  process_summary%snowice_surface .and. (measurement(band_0124) > Rth_snowice) ) .OR. &
       ( (process_summary%ocean_surface .or. process_summary%coastal_surface) &
                                     .and. (measurement(band_0086) > Rth_ocean) )   ) Then

     ! PART V (via PART I):
     if (.not. cloudy_250 .and. vis1km_test) then
        not_cloud = .true.
        CSR_QA = 3 ! restored by 250m test
     end If
     return
     
  End if

	rbox_x = size(reflectance_box, 1)
	rbox_y = size(reflectance_box, 3)
    
  Rmeas = 0 ! WDR-UIV
! Pick the reflectance box according to the surface type
  if (process_summary%land_surface .or. process_summary%desert_surface) then 
      Rbox(1:rbox_x, 1:rbox_y) = reflectance_box(:,band_0065,:)
      Rmeas = measurement(band_0065)
  endif
      
  if (process_summary%ocean_surface .or. process_summary%coastal_surface) then
      Rbox(1:rbox_x, 1:rbox_y) = reflectance_box(:,band_0086,:)
      Rmeas = measurement(band_0086)
  endif

  if (process_summary%snowice_surface) then
      Rbox(1:rbox_x, 1:rbox_y) = reflectance_box(:,band_0124,:)
      Rmeas = measurement(band_0124)
  endif

  average = array_mean(Rbox)
  standard_deviation_val = standard_deviation(Rbox(1:rbox_x, 1:rbox_y), average)

  if(average > 0.001) then
     dispersion_val = 100 * standard_deviation_val/average
  else
     dispersion_val = -1.
  end if


! *** For spatial varobility clear sky restoral based on box absolute sdev: TO BE USED FOR COLLECTION 5 OPERATIONAL VERSION

! 1km spatial variability (sdev) thresholds

!  early version of sdev thresholds:
!  if (process_summary%ocean_surface) then 
!     th_Ri = 0.004
!     th_Ri_QA = 0.003
!  else
!     th_Ri = 0.005
!     th_Ri_QA = 0.004
!  endif

! sdev threshold test 1:
!  if (process_summary%ocean_surface) then 
!     th_Ri = 0.005
!     th_Ri_QA = 0.004
!  else
!     th_Ri = 0.006
!     th_Ri_QA = 0.005
!  endif

! sdev threshold test 2: TO BE USED FOR COLLECTION 5 OPERATIONAL VERSION
  if (process_summary%ocean_surface) then 
     th_Ri    = 0.006
     th_Ri_QA = 0.005
  else
     th_Ri = 0.007
     th_Ri_QA = 0.006
  endif

  IF (Rmeas > 0. .and. standard_deviation_val < th_Ri) THEN ! candidate for restoral

! *** For spatial variability clear sky restoral based on box dispersion: FOR TESTING PURPOSES

! dispersion = 100*sdev/avg thresholds:
!  if (process_summary%ocean_surface) then 
!     dispersion_Ri    = 3.5
!     dispersion_Ri_QA = 3.0
!  else
!     dispersion_Ri    = 4.0
!     dispersion_Ri_QA = 3.5
!  endif

! IF (Rmeas > 0. .and. dispersion_val > 0. .and. dispersion_val < dispersion_Ri) THEN ! candidate for restoral


! ENTER PART II: ALTITUDE INDICATOR

! sep, 17May: previous version used IR phase, this version is using "not ice " from the decision tree phase 
!     If ( (.not. process_summary%icecloud) & !.and. measurement(band_0138) > -1.0 &
!          .and. measurement(band_0138) < 0.1) Then ! candidate for restoral
!	 If ( (.not. process_summary%icecloud) & !.and. measurement(band_0138) > -1.0 &
!		  .and. (measurement(band_0138) < 0.1) .and. (ctp > 350.0)) Then ! candidate for restoral
!		  .and. (measurement(band_0138) < 0.1) .and. (chm .eq. 6)) Then ! candidate for restoral
	 If ((measurement(band_0138) .lt. 0.1) .and. (chm .gt. 2)) Then ! candidate for restoral	! KGM 3-4-13

! ENTER PART III: SPECTRAL BEHAVIOR

        if (measurement(band_0138) > 0.015) then
           
           if ( (1.0 > r21_r138 .and. r21_r138 > r124_r21 .and. r124_r21 > r086_r124) .or. &
                (r21_r138 < r124_r21 .and. r124_r21 < r086_r124 .and. r086_r124 < 1.0)) then 
              
              !  Restore pixel and thereby pass through PART IV (branch 2 in flow chart)
              not_cloud = .true.
              CSR_QA = 2
           else
              !  PART V (via PART III): probaby thin cirrus
              if (.not. cloudy_250 .and. vis1km_test) then	! KGM 3-4-13
                 not_cloud = .true.
                 CSR_QA = 3 ! restored by 250m test
              end if
              return
           endif

        else
           !  Restore pixel and thereby pass through PART IV (branch 1 in flow chart)
           not_cloud = .TRUE.
           CSR_QA = 2
        endif

     Else

        !  PART V (via PART II): probably cirrus or high cloud
        if (.not. cloudy_250 .and. vis1km_test) then	! KGM 3-4-13
           not_cloud = .true.
           CSR_QA = 3 ! restored by 250m test
        end if
        return
        
     End If

  ELSE

     !  PART V (via PART I):
     if (.not. cloudy_250 .and. vis1km_test) then	! KGM 3-4-13
        not_cloud = .true.
        CSR_QA = 3 ! restored by 250m test
     end if
     return
     
  END IF

! ENTER PART IV: check variability again, if not_cloud is .true. but variability is still
!          fairly high then process as cloud but only lower the retrieval confidence QA

! *** For spatial varobility clear sky restoral based on box absolute sdev: TO BE USED FOR COLLECTION 5 OPERATIONAL VERSION
  if ( not_cloud .and. standard_deviation_val > th_Ri_QA ) then

! *** For spatial variability clear sky restoral based on box dispersion: FOR TESTING PURPOSES
! if ( not_cloud .and. dispersion_val > dispersion_Ri_QA ) then

         not_cloud = .FALSE.                      ! RESET PIXEL TO CLOUDY
         CSR_QA = 0
         lowvariability_confidence_test = .TRUE.  
                                                  
         ! PART V (via PART IV):
         if (.not. cloudy_250 .and. vis1km_test) then	! KGM 3-4-13
            not_cloud = .true.
            CSR_QA = 3 ! restored by 250m test
         end If

  endif

end subroutine cloudiness_test

function array_mean(array)
  real*4 :: array_mean, temp
  real*4, dimension(:,:), intent(in) :: array
  integer :: N1, N2, i, j, cnt

  N1 = size(array, 1)
  N2 = size(array, 2)

  temp = 0.
  cnt = 0
  do i=1, N1
     do j=1, N2
        if (array(i,j) > 0.) then 
           temp = temp + array(i,j)
           cnt = cnt + 1
        endif 
     end do
  end do

  if (cnt > 0) then 
     array_mean = temp / real(cnt)
  else
     array_mean = 0.
  endif 

end function array_mean

!=====================================
function standard_deviation(array, average)
!=====================================
!-----------------------------------------------------------------------
! !f90
!
! !Description:
!       This function computes standard deviation of a 2D array
!
! !Input Parameters:
!     array -- real*4, dimension(:,:) -- an array of data
!     average -- real*4 -- mean value of the array
!
! !Output Parameters:
!     standard_deviation -- real*4 -- the standard deviation value
!
! !Revision History:
!     2004/06/02 wind: Gala Wind
!     Revision 1.0 clean-up
!
! !Team-Unique Header:
!    Developed by the  Cloud Retrieval Group, NASA GSFC, Greenbelt, Maryland, USA.
!
! !References and Credits:
!      Written by:
!      Gala Wind
!      L3 GSI
!      Code 913, NASA/GSFC
!      Greenbelt, MD 20771
!      wind@climate.gsfc.nasa.gov
!
! !Design Notes:
!       NONE
!
! !END
!
!-----------------------------------------------------------------------
  real, intent(in) :: average
  real, dimension(:,:), intent(in) :: array
  real :: standard_deviation, temp
  integer :: i,j, cnt, N1, N2

  N1 = size(array, 1)
  N2 = size(array, 2)
  temp = 0.
  cnt = 0
  do i=1, N1
    do j=1, N2
      if ( array(i,j) > 0. ) then 
        temp = temp + (average - array(i,j))**2
        cnt = cnt + 1
      endif
    end do
  end do

  if ( (cnt-1) > 0) then 
     standard_deviation= sqrt ( temp / real(cnt - 1) )
  else
     standard_deviation = 0.
  endif

end function standard_deviation




end module clear_sky_restoral

