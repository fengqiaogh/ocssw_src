 module interpolate_libraries
 
	implicit none
	
	private
	
	public :: libraryinterpolate, ScatAngle, interpolate_wind_speed, lib_clean, lib_init ! GetReflForGivenWindSpeedIce, GetReflForGivenWindSpeedWater
	
!	REAL*8 					:: 	CONST_PI,	RAD_PER_DEG , DEG_PER_RAD 	
	real, parameter, dimension(3) :: wind_speed = (/3.0, 7.0, 15.0/)
	integer :: NUM_OF_TAUS, NUM_OF_WAVELENGTHS_RFTF
	
	real :: RLphase
	real, dimension(:), allocatable :: aeroPhase
	real, dimension(:,:), allocatable :: fOmega_water, fOmega_ice
	real, dimension(:,:,:), allocatable :: tauC_source_water, tauC_source_ice
	
	REAL,ALLOCATABLE::multiScatRefl_water_lamb(:,:,:)
	REAL,ALLOCATABLE::multiScatRefl_ice_lamb(:,:,:)

	REAL,ALLOCATABLE::multiScatRefl_water_CM(:,:,:,:)
	REAL,ALLOCATABLE::multiScatRefl_ice_CM(:,:,:,:)

! single-scattering arrays	
	real, allocatable :: fPhase_lamb(:,:), tauC(:,:)
	REAL, ALLOCATABLE :: gsquare(:), aero_denom(:)

	integer :: lib_call_counter

        REAL,ALLOCATABLE::phaseFunVals_water(:,:), phaseFunVals_ice(:,:), ssRefl_water(:,:,:), ssRefl_ice(:,:,:)
        real, dimension(:), allocatable :: LUT_TAU_VALUES

	real:: COS_SCAT, COS_VZA, COS_SZA


! ** extra SS stuff

	real:: prev_ratio 
	
	REAL, ALLOCATABLE ::  fPhase_liq(:,:),tauLR_liq(:,:),EXP0_liq(:,:),EXP1_liq(:,:)
	REAL, ALLOCATABLE ::  fPhase_ice(:,:),tauLR_ice(:,:),EXP0_ice(:,:),EXP1_ice(:,:)
	REAL, dimension(:), allocatable ::ocean_aero_Omega, ocean_aero_phase

! ** end extra SS stuff	


 contains
 
 
 subroutine lib_init
	
	use libraryarrays
 
	integer :: nwl, iwl, itau, ntau, i
	real, dimension(:,:), allocatable :: fQe_water, fQe_ice

	lib_call_counter = 0
	
	nwl = number_wavelengths
	ntau = number_taus
	allocate(aeroPhase(nwl))
 
	allocate(fQe_water(number_wavelengths, number_waterradii), fOmega_water(number_wavelengths, number_waterradii))
	allocate(fQe_ice(number_wavelengths, number_iceradii), fOmega_ice(number_wavelengths, number_iceradii))
 
! ** extra SS stuff	

	ALLOCATE(fPhase_liq(1:nwl,1:number_waterradii),tauLR_liq(1:nwl,1:number_waterradii), &
				EXP0_liq(1:nwl,1:number_waterradii),EXP1_liq(1:nwl,1:number_waterradii))
	ALLOCATE(fPhase_ice(1:nwl,1:number_iceradii),tauLR_ice(1:nwl,1:number_iceradii), &
				EXP0_ice(1:nwl,1:number_iceradii),EXP1_ice(1:nwl,1:number_iceradii))

	allocate(ocean_aero_Omega(nwl), ocean_aero_phase(nwl))
   ! WDR make re-entrant
	if( .not. allocated(aero_denom) ) allocate(aero_denom(nwl))

	do i=1, nwl
		if (aerosol_tau(i) > 0.) then 
			ocean_aero_Omega(i) = ((aerosol_tau(i) -0.1) * 1.0 + 0.1 * aerosol_ssa(i))/aerosol_tau(i)
		else
			ocean_aero_Omega(i) = 0. 
		endif
	end do

	aero_denom = (aerosol_tau -0.1) + 0.1 * aerosol_ssa 

	prev_ratio = -999.


! ** end extra SS stuff	


 
 	fQe_water(1,:) = 1.0
	do iwl = 2,nwl
		fQe_water(iwl,:) = extinction_water(iwl,:)/extinction_water(1,:)
	enddo
	
	fOmega_water = truncation_factor_water * singlescattering_water              !ssalbedo * ftrunc= (omega*f)

 	fQe_ice(1,:) = 1.0
	do iwl = 2,nwl
		fQe_ice(iwl,:) = extinction_ice(iwl,:)/extinction_ice(1,:)
	enddo
	
	fOmega_ice = truncation_factor_ice * singlescattering_ice              !ssalbedo * ftrunc= (omega*f)

	allocate(tauC_source_water(number_taus+1, number_wavelengths, number_waterradii))
	allocate(tauC_source_ice(number_taus+1, number_wavelengths, number_iceradii))

	tauC_source_water(1,:,:) = 0.0
        tauC_source_ice(1,:,:) = 0.0

	do itau = 2, number_taus+1
		tauC_source_water(itau,:,:) = (library_taus(itau-1) * fQe_water) * (1.0 - fOmega_water)   !Qe sacling, tau scaling, phase func truncation
		tauC_source_ice(itau,:,:) = (library_taus(itau-1) * fQe_ice) * (1.0 - fOmega_ice)   !Qe sacling, tau scaling, phase func truncation
	end do

	deallocate(fQe_water, fQe_ice)
 
!  ice_int(ntau+1,number_wavelengths,number_iceradii) ) 	

	allocate( multiScatRefl_water_lamb(ntau+1,number_wavelengths,number_waterradii), &
			  multiScatRefl_ice_lamb(ntau+1,number_wavelengths,number_iceradii), &
 			  multiScatRefl_water_CM(ntau+1,number_wavelengths,number_waterradii, 3), &
			  multiScatRefl_ice_CM(ntau+1,number_wavelengths,number_iceradii, 3) )
 		
 	allocate( fPhase_lamb(nwl, number_waterradii))
 			
 	allocate( gsquare(nwl))
		
 

	gsquare = aerosol_asym*aerosol_asym
	aeroPhase = 0.0

      NUM_OF_TAUS = number_taus + 1
      NUM_OF_WAVELENGTHS_RFTF = 2
     

      allocate(LUT_TAU_VALUES(NUM_OF_TAUS))
       LUT_TAU_VALUES(1) = 0.
       LUT_TAU_VALUES(2:NUM_OF_TAUS) = library_taus(1:number_taus)

   ALLOCATE(phaseFunVals_ice(number_wavelengths,number_iceradii) )
   !phaseFunVals_ice = 0.  ! WDR UIV again?
!        ALLOCATE(ssRefl_water(NUM_OF_TAUS,number_wavelengths,number_waterradii), &
!                 ssRefl_ice(NUM_OF_TAUS,number_wavelengths,number_iceradii), &
!                 phaseFunVals_water(number_wavelengths,number_waterradii),  &
!                 phaseFunVals_ice(number_wavelengths,number_iceradii))
   ALLOCATE(ssRefl_water(NUM_OF_TAUS,number_wavelengths,number_waterradii))
   ALLOCATE(phaseFunVals_water(number_wavelengths,number_waterradii) )
   !phaseFunVals_water = 999.  ! WDR UIV again?
   ALLOCATE(ssRefl_ice(NUM_OF_TAUS,number_wavelengths,number_iceradii) )

 
 end subroutine lib_init
 
 
 subroutine lib_clean
 
	deallocate(aeroPhase)
	deallocate(fOmega_water, fOmega_ice)
 
	deallocate(tauC_source_water, tauC_source_ice)
	
	deallocate(multiScatRefl_water_lamb, multiScatRefl_ice_lamb)
	deallocate(multiScatRefl_water_CM, multiScatRefl_ice_CM)
	
	deallocate(fPhase_lamb)
	deallocate(gsquare)
	deallocate(LUT_TAU_VALUES) 
        DEALLOCATE(phaseFunVals_water, phaseFunVals_ice, ssRefl_water, ssRefl_ice)

! ** extra SS stuff	
	deallocate(fPhase_liq, tauLR_liq, EXP0_liq, EXP1_liq)
	deallocate(fPhase_ice, tauLR_ice, EXP0_ice, EXP1_ice)
	deallocate(ocean_aero_Omega, ocean_aero_phase)
! ** end extra SS stuff	


 end subroutine lib_clean
 
 
 
 subroutine libraryinterpolate(local_solarzenith, &
                              local_sensorzenith, &
                              local_relativeazimuth, &
							  local_scatangle, &
							  local_wind_speed, &
							  wind_speed_only, interp_MS, interp_SS, &
                              debug, &
                              status, i, j)
  !
  !  so, this interpolates library values to the current conditions of
  !  1st 5 inputs
  !    What it interpolates?
   use core_arrays
   use GeneralAuxType
   use libraryarrays
   use libraryinterpolates
   use science_parameters
   
   implicit none

   real, intent(in)   :: local_solarzenith, local_sensorzenith, local_relativeazimuth, local_wind_speed
   real, intent(in) :: local_scatangle
   logical,intent(in)   :: debug
   integer, intent(in) ::  i, j
   integer,intent(out)  :: status
   logical, intent(in) :: wind_speed_only, interp_MS, interp_SS

	real :: dtheta
	integer :: iAngHi, iAngLow

   status = 0
   
   
!   	CONST_PI	=	4.d0 * ATAN(1.d0) !value of PI
!	RAD_PER_DEG = CONST_PI/180.d0	!radians per degree	
!	DEG_PER_RAD = 180.d0/CONST_PI		!degrees per radian

	COS_VZA = cos(local_sensorzenith * d2r)
	COS_SZA = cos(local_solarzenith * d2r)
	COS_SCAT  = cos(local_scatangle * d2r) 
!	NUM_OF_TAUS = number_taus + 1
!	NUM_OF_WAVELENGTHS_RFTF = 2

! For Lambertian libraries: land, also ocean when 0.86um saturates.
	if (.not. COX_MUNK) then 
	
        int_reflectance_ice = 0.
        int_reflectance_water = 0.

		call GetRefl(COS_SZA, COS_VZA, local_relativeazimuth, local_scatangle, &
					int_reflectance_water(:,:,:), int_reflectance_ice(:,:,:), interp_MS, interp_SS, status)
      if( status == 5 ) return  ! outside table

		if (interp_MS) then 

	        int_reflectance_water_sdev = 0.
    	    int_reflectance_ice_sdev = 0.

         call GetSdevReflLamb(COS_SZA, COS_VZA, local_relativeazimuth,&
                        int_reflectance_water_sdev(:,:,:), &
                        int_reflectance_ice_sdev(:,:,:), &
                        status)
         if( status == 5 ) return  ! out of table space
			call Setup_emissivity_flux(COS_SZA, library_fluxsolarzenith, iAngHi, iAngLow, dtheta)			

			call InterpolateFluxes(COS_SZA, library_fluxsolarzenith, flux_down_ice_solar, int_fluxdownice_solar, &
						iAngHi, iAngLow, dtheta, status)
			call InterpolateFluxes(COS_SZA, library_fluxsolarzenith, flux_down_water_solar, int_fluxdownwater_solar, &
						iAngHi, iAngLow, dtheta, status)
		
			call Setup_emissivity_flux(COS_VZA, library_fluxsensorzenith, iAngHi, iAngLow, dtheta)			

			call InterpolateFluxes(COS_VZA, library_fluxsensorzenith, flux_down_ice_sensor, int_fluxdownice_sensor, &
						iAngHi, iAngLow, dtheta, status)
			call InterpolateFluxes(COS_VZA, library_fluxsensorzenith, flux_up_ice_sensor, int_fluxupice_sensor, &
						iAngHi, iAngLow, dtheta, status)
	
			call InterpolateFluxes(COS_VZA, library_fluxsensorzenith, flux_down_water_sensor, int_fluxdownwater_sensor, &
						iAngHi, iAngLow, dtheta, status)
			call InterpolateFluxes(COS_VZA, library_fluxsensorzenith, flux_up_water_sensor, int_fluxupwater_sensor, &
						iAngHi, iAngLow, dtheta, status)

		endif

	else

        int_reflectance_ice = 0.
        int_reflectance_water = 0.

        call GetReflForGivenWindSpeed(COS_SZA, COS_VZA, local_relativeazimuth, &
          local_scatangle, COS_SCAT, local_wind_speed, &
          int_reflectance_water(:,:,:), int_reflectance_ice(:,:,:), &
          wind_speed_only, interp_MS, interp_SS, status)
      if( status == 5 ) return

		if (interp_MS) then 								
										
			int_reflectance_water_sdev = 0.
            int_reflectance_ice_sdev = 0.
							
        call GetSdevReflForGivenWindSpeed(COS_SZA, COS_VZA, &
           local_relativeazimuth, local_wind_speed, &
           int_refl_water_sdev_wspeed(:,:,:,:), &
           int_reflectance_water_sdev(:,:,:), &
           int_refl_ice_sdev_wspeed(:,:,:,:), &
           int_reflectance_ice_sdev(:,:,:), &          
           wind_speed_only, status)
         if( status == 5 ) return  !  point beyond table range

			call  Setup_emissivity_flux(COS_VZA, library_sensor_zenith, iAngHi, iAngLow, dtheta) 

			call InterpEmissForAGivenWSpeed(COS_VZA, local_wind_speed, surface_emissivity_water,&
							int_surface_emis_water_wspeed, int_surface_emissivity_water, &
							wind_speed_only, iAngHi, iAngLow, dtheta, status)
			call InterpEmissForAGivenWSpeed(COS_VZA, local_wind_speed, surface_emissivity_ice, &
							int_surface_emis_ice_wspeed, int_surface_emissivity_ice, &
							wind_speed_only, iAngHi, iAngLow, dtheta,  status)
			call InterpEmissForAGivenWSpeed(COS_VZA, local_wind_speed, cloud_emissivity_water, &
							int_cloud_emis_water_wspeed, int_cloud_emissivity_water, &
							wind_speed_only,  iAngHi, iAngLow, dtheta, status)
			call InterpEmissForAGivenWSpeed(COS_VZA, local_wind_speed, cloud_emissivity_ice, &
							int_cloud_emis_ice_wspeed, int_cloud_emissivity_ice, &
							wind_speed_only,  iAngHi, iAngLow, dtheta, status)

! These quantities are not used for Collection 6 calculations
			call InterpEmissForAGivenWSpeed(COS_VZA, local_wind_speed, surface_emissivity_water_sdev,&
							int_surface_emis_water_sdev_wspeed, int_surface_emissivity_water_sdev, &
							wind_speed_only,  iAngHi, iAngLow, dtheta, status)
			call InterpEmissForAGivenWSpeed(COS_VZA, local_wind_speed, surface_emissivity_ice_sdev, &
							int_surface_emis_ice_sdev_wspeed, int_surface_emissivity_ice_sdev, &
							wind_speed_only,  iAngHi, iAngLow, dtheta, status)

			call InterpEmissForAGivenWSpeed(COS_VZA, local_wind_speed, cloud_emissivity_water_sdev, &
						int_cloud_emis_water_sdev_wspeed, int_cloud_emissivity_water_sdev, &
						wind_speed_only,  iAngHi, iAngLow, dtheta, status)
			call InterpEmissForAGivenWSpeed(COS_VZA, local_wind_speed, cloud_emissivity_ice_sdev, &
						int_cloud_emis_ice_sdev_wspeed, int_cloud_emissivity_ice_sdev, &
						wind_speed_only,  iAngHi, iAngLow, dtheta, status)
		endif

	endif




end subroutine libraryinterpolate

REAL FUNCTION ScatAngle(solarAng,viewAng,relAzm)
!................................................................................
!
!................................................................................
! !F90
!
! !DESCRIPTION:     
!     This function computes the scattering angle for a given sun-satellite geometry! 
!
! !INPUT PARAMETERS:
!	solarAng : solar angle in degrees (scalar)
!	viewAng	 : view angle in degrees (scalar)
!	relAzm	 : rel azimuth (from the principle plane of the SUN) (scalar)
!
!
! !OUTPUT PARAMETERS:
!	scatAngle : scattering angle in degrees (scalar)
! 
!................................................................................
! !Revision History:
!		Revision 1.0  2009/06/03 GW
!		Initial integration of the standalone version. No changes from the standalone.
!		
!................................................................................
! !Team-unique Header:
!             Cloud Retrieval Group, NASA/GSFC
!
! !PROGRAMMER:
!
!             Nandana Amarasinghe (SSAI)
!             Climate and Radiation Branch
!             NASA Goddard Space Flight Center
!             Greenbelt, Maryland 20771, U.S.A.
!
!             Gala Wind (SSAI)
!             Climate and Radiation Branch
!             NASA Goddard Space Flight Center
!             Greenbelt, Maryland 20771, U.S.A.
!
!*******************************************************************************
! !END
	use science_parameters, only:d2r

	IMPLICIT NONE
	!CALLING ARGUMENTS
	REAL, INTENT(IN):: solarAng,viewAng,relAzm

	!LOCAL VARIABLES
	REAL::solarmu,viewmu !,pi,torad,todeg

!	pi		= acos(-1.d0 )
!	torad	= pi / 180.d0
!	todeg	= 180.d0 / pi

	solarmu = COS(solarAng * d2r)
	viewmu = COS(viewAng * d2r)

	scatAngle= -solarmu*viewmu + SQRT( ( 1.- solarmu**2 )* &
                      ( 1.- viewmu**2 ) )*COS( relAzm*d2r)


    IF(scatAngle < -1.0)scatAngle = -1.0 !guard against round off error 
    IF(scatAngle > 1.0)scatAngle = 1.0
    	
                     
	scatAngle = ACOS(scatAngle) / d2r !* todeg

END FUNCTION ScatAngle
 
 
SUBROUTINE GetPhaseFunctionValues(scatAngle, scatAngleArray, num_angles, phaseFunArray, phaseNormConst, & 
	                                                                 phaseFunVals, ierror)
!................................................................................
!
!................................................................................
! !F90
!
! !DESCRIPTION:     
!	This program interpolates the phase fun value array for a given scattering angle an returns
!	the corresponding values for particular lambda and re (see the inputs and outputs)!
! !INPUT PARAMETERS:
	!	scatAngle		:	scattering angle in degrees (scalar)
	!	scatAngleArray	:	1-D array of sacattering angles for which phase fun vals are precomputed
	!	phaseFunArray	:	3-D array of phase fun vals: 1st DIM = num of angles, 
	!						2nd DIM = num of wavelengths, 3rd DIM = num of REs
	!	phaseNormConst	:	Normalization factor for the phase function. If it is normalized to 1
	!						then it is OK. But, this number (= pmom(0)) from the dfit routines  has
	!						been saved and included in the library file in adition to d_fit coeffs.
!
!
! !OUTPUT PARAMETERS:
	!	phaseFunVals	:	2-D array of interpolated phase fun values at scatAngle
	!						1st dim = num of wavelengths, 2nd dim = num of rE's
	!	ierror			:	integer value
	!						ierror 	= 0  	success, 
	!								= 3		invalid data
!................................................................................
! !Revision History:
!		Revision 1.0  2009/06/03 GW
!		Initial integration of the standalone version. No changes from the standalone.
!		
!................................................................................
! !Team-unique Header:
!             Cloud Retrieval Group, NASA/GSFC
!
! !PROGRAMMER:
!
!             Nandana Amarasinghe (SSAI)
!             Climate and Radiation Branch
!             NASA Goddard Space Flight Center
!             Greenbelt, Maryland 20771, U.S.A.
!
!             Gala Wind (SSAI)
!             Climate and Radiation Branch
!             NASA Goddard Space Flight Center
!             Greenbelt, Maryland 20771, U.S.A.
!
!*******************************************************************************
! !END

 use libraryarrays, only: number_wavelengths
 use core_arrays, only: platform_name
 use ch_xfr, only : c2_sensor_id, OCI_ID, OCIS_ID
	IMPLICIT NONE
	!Arguments
	REAL,INTENT(IN)::scatAngle, scatAngleArray(:),phaseFunArray(:,:,:),phaseNormConst(:,:)
	REAL,INTENT(OUT)::phaseFunVals(:,:)
	integer, intent(in):: num_angles
	INTEGER,INTENT(OUT)::ierror

	!local variables
	INTEGER:: iAng, iAngLow, iAngHi, nwl
	REAL::dtheta
	CHARACTER(150)::tempstr

	!start execution
   !  WDR 13 Jun 2022 add the OCI platform_name in if()	
	if (platform_name == "AVIRIS" .or. platform_name == "RSP" .or. &
      platform_name == "AHI" .or. platform_name == "EPIC" .or. &
		platform_name == "SSFR" .or. c2_sensor_id == OCI_ID .or. &
      c2_sensor_id == OCIS_ID ) then 
		nwl = number_wavelengths
	else
		nwl = number_wavelengths-1
	endif
		

	ierror = 0 !start with no error
	iAngLow = 1
	iAngHi = num_angles
	IF(scatAngle > scatAngleArray(iAngHi))THEN
		iAngLow = iAngHi-1 !no need to do the search; do extrapolation
		
	ELSEIF(scatAngle < scatAngleArray(iAngLow))THEN
		iAngHi = iAngLow+1
	ELSE
		DO									!do a bisection search
			IF((iAngHi - iAngLow) == 1)EXIT
			iAng = (iAngHi + iAngLow)/2
			IF(scatAngleArray(iAng) < scatAngle)THEN
			iAngLow = iAng
			ELSE
				iAngHi = iAng
			ENDIF
		ENDDO
	ENDIF


	dtheta = scatAngleArray(iAngHi) - scatAngleArray(iAngLow)
	!print*,iAngHi,iAngLow,dtheta,scatAngleArray(iAngHi)-scatAngle,scatAngle - scatAngleArray(iAngLow)
!	IF(dtheta ==0)THEN
!		ierror = 3 			!invalid data read
!		RETURN
!	ENDIF
	
	phaseFunVals(1:nwl,:) = ( phaseFunArray(iAngLow,1:nwl,:)  * (scatAngleArray(iAngHi)-scatAngle) + &
                    phaseFunArray(iAngHi,1:nwl,:)  * (scatAngle - scatAngleArray(iAngLow) ))/dtheta 


    phaseFunVals(1:nwl,:) = phaseFunVals(1:nwl,:)/phaseNormConst(1:nwl,:)

	RETURN

END SUBROUTINE GetPhaseFunctionValues


subroutine get_aero_params(cos_scatAngle, aeroG)
	use libraryarrays, only: number_wavelengths
 use core_arrays, only: platform_name
 use ch_xfr, only : c2_sensor_id, OCI_ID, OCIS_ID
		
	real, intent(in) :: aeroG(:), cos_scatAngle
		
	INTEGER:: nwl
	
	if (platform_name == "AVIRIS" .or. platform_name == "RSP" .or. &
      platform_name == "AHI" .or. platform_name == "EPIC" .or. &
		platform_name == "SSFR" .or. c2_sensor_id == OCI_ID .or. &
      c2_sensor_id == OCIS_ID ) then 
		nwl = number_wavelengths
	else
		nwl = number_wavelengths-2
	endif
		
	RLphase = 0.75 * (1 + cos_scatAngle*cos_scatAngle)

! it is generally advisable to avoid a real power whenever possible in order to make the code run faster
! GW: 10.25.13
	aeroPhase(1:nwl) = (1.0 - gsquare(1:nwl))/ &
							(sqrt(1 + gsquare(1:nwl) &
							- 2.0 * aerog(1:nwl) * cos_scatAngle)**3)

end subroutine get_aero_params


SUBROUTINE InterpolateFluxes(solarOrViewAng, solarAngMuArray, inFluxArray, outFluxArray, iAngHi, iAngLow, dtheta, ierror)
!................................................................................
! !F90
!
! !DESCRIPTION:     
	!This program interpolates(linear) the RF or TF array for a given solar/oview angle and returns
	!the corresponding values for particular tau, lambda and re (see the inputs and outputs)
! !INPUT PARAMETERS:
	!	solarOrViewAng		:	input angle in degrees (scalar)
	!	solarAngleMuArray	:	1-D array of solar angles for which RF/TF vals are precomputed
	!	inFluxArray			:	4-D array of flux vals: 1st DIM = num of angles, 2nd DIM = ntaus
	!						 	3rdd DIM = num of wavelengths, 4th DIM = num of REs
!
!
! !OUTPUT PARAMETERS:
	!	outFluxArray		: 3-D array of interpolated RF or TF values at the input angle
	!						  1st dim = ntaus, 2nd dim = num of wavelengths, 3rdd dim = num of rE's
	!	ierror			:	integer value
	!						ierror 	= 0  	success, 
	!								= 3		invalid data
!................................................................................
! !Revision History:
!		Revision 1.0  2009/06/03 GW
!		Initial integration of the standalone version. No changes from the standalone.
!		
!................................................................................
! !Team-unique Header:
!             Cloud Retrieval Group, NASA/GSFC
!
! !PROGRAMMER:
!
!             Nandana Amarasinghe (SSAI)
!             Climate and Radiation Branch
!             NASA Goddard Space Flight Center
!             Greenbelt, Maryland 20771, U.S.A.
!
!             Gala Wind (SSAI)
!             Climate and Radiation Branch
!             NASA Goddard Space Flight Center
!             Greenbelt, Maryland 20771, U.S.A.
!
!*******************************************************************************
! !END
	IMPLICIT NONE
	!Arguments
	REAL,INTENT(IN):: solarOrViewAng, solarAngMuArray(:),inFluxArray(:,:,:,:)
	REAL,INTENT(OUT)::outFluxArray(:,:,:)
	INTEGER,INTENT(OUT)::ierror
	integer, intent(in) :: iAngLow, iAngHi
	real, intent(in) :: dtheta

	!local variables
	REAL*8::UMU

	real :: diff_one_way, diff_other_way

	!start esecution
	ierror = 0
	
	UMU  = solarOrViewAng 

	diff_one_way = solarAngMuArray(iAngHi)- UMU
	diff_other_way = UMU - solarAngMuArray(iAngLow)

	outFluxArray(:,:,:) = 0.
	outFluxArray(:,:,:) =  ( inFluxArray(iAngLow,:,:,:)  * diff_one_way + &
								inFluxArray(iAngHi,:,:,:)  * diff_other_way ) / dtheta

   RETURN

END SUBROUTINE InterpolateFluxes

	SUBROUTINE CalcSingScatPartOfReflectance(phaseFunVals_water, phaseFunVals_ice, &
							theta, theta0, ssRefl_water, ssRefl_ice, ierror)

!................................................................................
! !F90
!
! !DESCRIPTION:     
	!Calculates the single scattering part of the reflectance according to Nakajima/Tanaka method
	!implemented in DISORT Version2. We had to implement this way, because of the way single 
	!scattering part is takent out in DISORT. Please refer to the DISORT2 documentation
	!Formula Used here: refl = omega * PC * {1.0 - EXP[-tauc *(1/mu + 1/mu0)]}
	!                   refl = refl/(4.0 * (mu + mu0))
	!				    where tauc = tau * [Qe(lambda,re)/Qe(CH1,re)] * (1.0 - ftrunc * omega)
	!						   PC  = Ptrue / (1 - ftrunc * omega)
	!					also remember refl_lUT = (I_disort) * pi/( mu0* solarflux)
	!					so that multiplication is already taken care of here, no need to do it
	!					outside of this routine
! !INPUT PARAMETERS:
	!	PhaseFunVals : phase function calculated for all the wavelengths and effective radii
	!                  DIMS[1:nlambda,1:nradii] (look at the program getPhaseFunValues.f90)
	!   ssAlbedo     : single scattering albedo DIMS[1:nlambda,1:nradii]
	!	extCoeff	 : extinction coefficient  DIMS[1:nlambda,1:nradii]
	!	dfitF        : phase function truncation factor in dfit method, DIMS[1:nlambda,1:nradii]
	!   tauVals      : opticat thickness values in the LUT, DIMS[1:ntau]
	!   theta        : view angle in degrees
	!   theta0       : solar angle in degrees
!
!
! !OUTPUT PARAMETERS:
	!   ssRefl       : single scattering part of the reflectance , here no need to input the beam
	!                 strength, since we calculate the reflectance.
	!	ierror			:	integer value
	!						ierror 	= 0  	success, 
	!								= 2 	memory allocation error
	!								= 3		invalid data
!................................................................................
! !Revision History:
!		Revision 1.0  2009/06/03 GW
!		Initial integration of the standalone version. No changes from the standalone.
!		
!................................................................................
! !Team-unique Header:
!             Cloud Retrieval Group, NASA/GSFC
!
! !PROGRAMMER:
!
!             Nandana Amarasinghe (SSAI)
!             Climate and Radiation Branch
!             NASA Goddard Space Flight Center
!             Greenbelt, Maryland 20771, U.S.A.
!
!             Gala Wind (SSAI)
!             Climate and Radiation Branch
!             NASA Goddard Space Flight Center
!             Greenbelt, Maryland 20771, U.S.A.
!
!*******************************************************************************
! !END
	use libraryarrays
 use core_arrays, only: platform_name
 use ch_xfr, only : c2_sensor_id, OCI_ID, OCIS_ID

	REAL, INTENT(IN)  :: phaseFunVals_water(:,:), phaseFunVals_ice(:,:)
	REAL, INTENT(IN)  :: theta, theta0
	REAL, INTENT(OUT) :: ssRefl_water(:,:,:), ssRefl_ice(:,:,:)
	INtEGER,INTENT(OUT)::ierror
	
	
	INTEGER:: ntau, nwl, nre, itau,iwl, allocateStatus, ii,jj
	REAL ::  MUPLUSMU0,ratio, mytau
	
	ntau = number_taus !size(tauVals)

	if( platform_name == "AVIRIS" .or. platform_name == "RSP" .or. &
      platform_name == "AHI" .or. platform_name == "EPIC" .or. &
		platform_name == "SSFR" .or. c2_sensor_id == OCI_ID .or. &
      c2_sensor_id == OCIS_ID ) then 
		nwl = number_wavelengths + 1 ! we have to make sure not to skip the last channel
	else
		nwl = number_wavelengths
	endif
				

	MUPLUSMU0 = theta + theta0
	ratio = MUPLUSMU0/(theta * theta0)
	

	fPhase_lamb = phaseFunVals_water/ (1.0 - fOmega_water)   !add something here to safeguard div by zero
	
	ssRefl_water = 0.0
	ssRefl_ice = 0.0
	
	do itau = 1, ntau
	
		do ii=1, nwl-1
			do jj=1, number_waterradii

				mytau = tauC_source_water(itau+1,ii,jj)*ratio

				if (mytau == 0) then 
					EXP1_liq(ii,jj) = 1.
				else if (mytau > 10.) then 
					EXP1_liq(ii,jj) = 0.
				else  
					EXP1_liq(ii,jj) = exp(-mytau)

!					EXP1_liq(ii,jj) = 1. + (mytau ) + (mytau )**2/2. + (mytau )**3/6. + &
!                                                               (mytau )**4/24. + (mytau )**5/120. + mytau**6/720. 
!					EXP1_liq(ii,jj) = 1. / EXP1_liq(ii,jj)
				endif
	
			ssRefl_water(itau,ii,jj) = (singlescattering_water(ii,jj)* fPhase_lamb(ii,jj) * (1 - EXP1_liq(ii,jj))) / (4.0*MUPLUSMU0) 
			
			end do
		end do
	
	
	enddo
	
	fPhase_lamb(:, 1:number_iceradii) = phaseFunVals_ice/ (1.0 - fOmega_ice)   !add something here to safeguard div by zero
	
	do itau = 1, ntau


                do ii=1, nwl-1
                        do jj=1, number_iceradii


							mytau = tauC_source_ice(itau+1,ii,jj)*ratio

                                if (mytau == 0) then
                                        EXP1_ice(ii,jj) = 1.
                                else if (mytau > 10.) then
                                        EXP1_ice(ii,jj) = 0.
                                else 
										EXP1_ice(ii,jj) = exp(-mytau)
                                endif

						ssRefl_ice(itau,ii,jj) = (singlescattering_ice(ii,jj)* fPhase_lamb(ii,jj) * (1 - EXP1_ice(ii,jj))) / (4.0*MUPLUSMU0) 

                        end do
                end do


	enddo
	
	END SUBROUTINE CalcSingScatPartOfReflectance	


	SUBROUTINE Single_Scattering_Calcs_Ocean(phaseFunVals_liq, phaseFunVals_ice, ssAlbedo_liq, ssAlbedo_ice, &
												RLphase, aeroPhase, RLTau,aeroTau, aeroOmega,&
												  theta, theta0,ssRefl_liq, ssRefl_ice)

	use libraryarrays, only: number_wavelengths, number_waterradii, number_iceradii
 use core_arrays, only: platform_name
 use ch_xfr, only : c2_sensor_id, OCI_ID, OCIS_ID

	REAL, INTENT(IN)  :: phaseFunVals_liq(:,:), ssAlbedo_liq(:,:),  &
						 phaseFunVals_ice(:,:), ssAlbedo_ice(:,:), &
	                     RLTau(:), aeroTau(:), &
	                     aeroOmega(:), RLphase, aeroPhase(:)
	REAL, INTENT(IN)  :: theta, theta0
	REAL, INTENT(OUT) :: ssRefl_liq(:,:,:), ssRefl_ice(:,:,:)
	
	
	INTEGER:: ntau, nwl, itau, iwl, wave_limit
	REAL :: MUPLUSMU0,ratio
	
	real :: ocean_mul(20)
	logical:: changed_ratio
	integer :: ii, jj

	ntau = NUM_OF_TAUS
	
	if (platform_name == "AVIRIS" .or. platform_name == "RSP" .or. &
      platform_name == "AHI" .or. platform_name == "EPIC" .or. &
		platform_name == "SSFR" .or. c2_sensor_id == OCI_ID .or. &
      c2_sensor_id == OCIS_ID ) then
		nwl = number_wavelengths+1 ! capture the last channel, it is not to be skipped
	else
		nwl = number_wavelengths
	endif
								
	MUPLUSMU0 = theta + theta0
	ratio = MUPLUSMU0/(theta * theta0)

	fPhase_liq = phaseFunVals_liq/ (1.0 - fOmega_water)   !add something here to safeguard div by zero
	fPhase_ice = phaseFunVals_ice/ (1.0 - fOmega_ice)
	

#if AVIRIS_INST | RSP_INST | AHI_INST | EPIC_INST | SSFR_INST
	do iwl = 1,  number_wavelengths
#else
   if( ( c2_sensor_id == OCI_ID ) .or. ( c2_sensor_id == OCIS_ID ) ) then
     wave_limit = number_wavelengths - 3
   else
     wave_limit = number_wavelengths - 2
   endif
	do iwl = 1,  wave_limit
#endif
		ocean_aero_phase(iwl) = ( (aeroTau(iwl) -0.1) * RLphase  + 0.1 * aeroOmega(iwl) * aeroPhase(iwl) ) / aero_denom(iwl) 
		ocean_mul(iwl) = ocean_aero_phase(iwl) * ocean_aero_Omega(iwl)
	end do

	EXP1_liq = 0.
	EXP1_ice = 0.

	DO itau = 1, ntau

		tauLR_liq =  tauC_source_water(itau,:,:) * ratio
		tauLR_ice =  tauC_source_ice(itau,:,:) * ratio


		do ii=1, nwl-1
			do jj=1, number_waterradii

				if (tauLR_liq(ii,jj) == 0) then 
					EXP1_liq(ii,jj) = 1.
				else if (tauLR_liq(ii,jj) > 10.) then 
					EXP1_liq(ii,jj) = 0.
				else 
					EXP1_liq(ii,jj) = exp(-tauLR_liq(ii,jj))
				endif

			end do
		end do

                do ii=1, nwl-1
                        do jj=1, number_iceradii

                                if (tauLR_ice(ii,jj) == 0) then
                                        EXP1_ice(ii,jj) = 1.
                                else if (tauLR_ice(ii,jj) > 10.) then
                                        EXP1_ice(ii,jj) = 0.
                                else
										EXP1_ice(ii,jj) = exp(-tauLR_ice(ii,jj))
                                endif

                        end do
                end do

	
		ssRefl_liq(itau,:,:) = ssAlbedo_liq * fPhase_liq * (1. - EXP1_liq)
		ssRefl_ice(itau,:,:) = ssAlbedo_ice * fPhase_ice * (1. - EXP1_ice)
		
		EXP0_liq = EXP1_liq
		EXP0_ice = EXP1_ice

#if AVIRIS_INST | RSP_INST | AHI_INST | EPIC_INST | SSFR_INST
		DO iwl =1,  number_wavelengths
#else
     if( ( c2_sensor_id == OCI_ID ) .or. ( c2_sensor_id == OCIS_ID ) ) then
       wave_limit = number_wavelengths - 3
     else
       wave_limit = number_wavelengths - 2
     endif
     do iwl = 1,  wave_limit
#endif
                       tauLR_liq(iwl,:) =  (tauC_source_water(itau,iwl,:) + RLTau(iwl))*ratio       
                       tauLR_ice(iwl,:) =  (tauC_source_ice(itau,iwl,:) + RLTau(iwl))*ratio

			ii = iwl
                        do jj=1, number_waterradii

                                if (tauLR_liq(ii,jj) == 0) then
                                        EXP1_liq(ii,jj) = 1.
                                else if (tauLR_liq(ii,jj) > 10.) then
                                        EXP1_liq(ii,jj) = 0.
                                else
										EXP1_liq(ii,jj) = exp(-tauLR_liq(ii,jj))
                                endif

                        end do

			ii = iwl
                        do jj=1, number_iceradii

                                if (tauLR_ice(ii,jj) == 0) then
                                        EXP1_ice(ii,jj) = 1.
                                else if (tauLR_ice(ii,jj) > 10.) then
                                        EXP1_ice(ii,jj) = 0.
                                else
										EXP1_ice(ii,jj) = exp(-tauLR_ice(ii,jj))
                                endif

                        end do

			
			ssRefl_liq(itau,iwl,:) = ssRefl_liq(itau,iwl,:) +  & 
			                            RLphase * 1.0 * (EXP0_liq(iwl,:) - EXP1_liq(iwl,:))
			ssRefl_ice(itau,iwl,:) = ssRefl_ice(itau,iwl,:) +  & 
			                            RLphase * 1.0 * (EXP0_ice(iwl,:) - EXP1_ice(iwl,:))
			
			EXP0_liq(iwl,:) = EXP1_liq(iwl,:)
			EXP0_ice(iwl,:) = EXP1_ice(iwl,:)
			
                       tauLR_liq(iwl,:) =  (tauC_source_water(itau,iwl,:) + RLTau(iwl) + aeroTau(iwl)) * ratio
                       tauLR_ice(iwl,:) =  (tauC_source_ice(itau,iwl,:) + RLTau(iwl) + aeroTau(iwl)) * ratio
			
			!coming in aeroTau is the (rayleightau + 0.1) from DISORT run

                       ii = iwl
                        do jj=1, number_waterradii

                                if (tauLR_liq(ii,jj) == 0) then
                                        EXP1_liq(ii,jj) = 1.
                                else if (tauLR_liq(ii,jj) > 10.) then
                                        EXP1_liq(ii,jj) = 0.
                                else
										EXP1_liq(ii,jj) = exp(-tauLR_liq(ii,jj))
                                endif

                        end do

                        ii = iwl
                        do jj=1, number_iceradii
                       
                                if (tauLR_ice(ii,jj) == 0) then
                                        EXP1_ice(ii,jj) = 1.
                                else if (tauLR_ice(ii,jj) > 10.) then
                                        EXP1_ice(ii,jj) = 0.
                                else
										EXP1_ice(ii,jj) = exp(-tauLR_ice(ii,jj))
                                endif

                        end do



			ssRefl_liq(itau,iwl,:) = ssRefl_liq(itau,iwl,:) + ocean_mul(iwl) * &
			                                                (EXP0_liq(iwl,:) - EXP1_liq(iwl,:))
			ssRefl_ice(itau,iwl,:) = ssRefl_ice(itau,iwl,:) + ocean_mul(iwl) * &
			                                                (EXP0_ice(iwl,:) - EXP1_ice(iwl,:))
		ENDDO
				 
	ENDDO

	
	ssRefl_liq = ssRefl_liq/(4.0 * MUPLUSMU0)
	ssRefl_ice = ssRefl_ice/(4.0 * MUPLUSMU0)
		


	END SUBROUTINE Single_Scattering_Calcs_Ocean





	!----------------------SUBROUTINE CalcSingScatPartOfOceanReflectance---------------------------------
	!Calculates the single scattering part of the reflectance according to Nakajima/Tanaka method
	!implemented in DISORT Version2. We had to implement this way, because of the way single 
	!scattering part is takent out in DISORT. Please refer to the DISORT2 documentation
	!Formula Used here: refl = omega * PC * {1.0 - EXP[-tauc *(1/mu + 1/mu0)]}
	!                   refl = refl/(4.0 * (mu + mu0))
	!				    where tauc = tau * [Qe(lambda,re)/Qe(CH1,re)] * (1.0 - ftrunc * omega)
	!						   PC  = Ptrue / (1 - ftrunc * omega)
	!					also remember refl_lUT = (I_disort) * pi/( mu0* solarflux)
	!					so that multiplication is already taken care of here, no need to do it
	!					outside of this routine
	!INPUTS:
	!	PhaseFunVals : phase function calculated for all the wavelengths and effective radii
	!                  DIMS[1:nlambda,1:nradii] (look at the program getPhaseFunValues.f90)
	!	RLphase		 : Rayleigh phase func value (fn of the scattering angle)
	!	aeroPhse	 : mixed phase function of rayleigh layer and the aerosol layer
	!   ssAlbedo     : single scattering albedo DIMS[1:nlambda,1:nradii]
	!	extCoeff	 : extinction coefficient  DIMS[1:nlambda,1:nradii]
	!	dfitF        : phase function truncation factor in dfit method, DIMS[1:nlambda,1:nradii]
	!   tauVals      : opticat thickness values in the LUT, DIMS[1:ntau]
	!	RLTau		 : Rayleigh optical thickness (combined several layers upto 8km)
	!	aeroTau		 : aerosol optical thickness(0.1) + rayleigh layer(lowest)
	!	aeroOmega    : ss albedo of aerosol (not the mixed, mixing will be done in SS calc routine)
	!   theta        : view angle in degrees
	!   theta0       : solar angle in degrees
	!OUTPUTS:
	!   ssRefl       : single scattering part of the reflectance , here no need to input the beam
	!                 strength, since we calculate the reflectance.
	!	ierror			:	integer value
	!						ierror 	= 0  	success, 
	!								= 2 	memory allocation error
	!								= 3		invalid data
	!-----------------------------------------------------------------------------------------------	
	
	subroutine Setup_emissivity_flux(angle, angle_array, idx_hi, idx_lo, dtheta) 
	
		real, intent(in) :: angle, angle_array(:)
		real, intent(out) :: dtheta
		integer, intent(out) :: idx_hi, idx_lo
	
	
		real :: UMU
		integer :: iAng
	
		UMU  = angle !COS(solarOrViewAng * RAD_PER_DEG)
		idx_lo = 1
		idx_hi = size(angle_array)
		IF(UMU < angle_array(idx_lo))THEN
			idx_hi = idx_lo+1 !no need to do the search; do extrapolation
			if (idx_hi > size(angle_array)) idx_hi = idx_lo ! we don't have enough points
		ELSE
			DO									!do a bisection search
				IF((idx_hi - idx_lo) == 1)EXIT
				iAng = (idx_hi + idx_lo)/2
				IF(angle_array(iAng) < UMU)THEN
					idx_lo = iAng
				ELSE
				idx_hi = iAng
				ENDIF
			ENDDO
		ENDIF


	dtheta = angle_array(idx_hi) - angle_array(idx_lo)
		
	end subroutine Setup_emissivity_flux
	
	
	!--------------------SUBROUTINE InterpCloudEmissForAGivenWSpeed ---------------------------------
	!This program interpolates(linear) the cloud emissivity arrays for a given wind speed and
	!returns the corresponding values for particular tau, lambda and re(see the inputs and outputs)
	!INPUTS:
	!	iphase			 	: 	integer  1 =Water Cloud, 2 =Ice Cloud
	!	solarOrViewAng		:	input angle in degrees (scalar)
	!	wspeed				:	real scalar, wind speed in m/s
	!OUTPUTS:
	!	outFluxArray		: 3-D array ofinterpolated cloud emissivity values at the input ang
	!						  1st dim = ntaus, 2nd dim = num of wavelengths, 3rdd dim = num of rE's
	!	ierror				:	integer value
	!						 	ierror 	= 0  	success, 
	!								= 3		invalid data
	!-----------------------------------------------------------------------------------------------
	! 						
	SUBROUTINE InterpEmissForAGivenWSpeed(solarOrViewAng,wspeed, inEmissArray, outEmissArray_ws, outEmissArray, &
												wspeed_only, iAngHi, iAngLow, dtheta,ierror)

   use libraryarrays

	IMPLICIT NONE
	!Arguments
	REAL,INTENT(IN):: solarOrViewAng,wspeed, dtheta
	real, intent(in) :: inEmissArray(:,:,:,:,:)
	REAL,INTENT(OUT)::outEmissArray(:,:,:)
	real, intent(inout) :: outEmissArray_ws(:,:,:,:)
	INTEGER,INTENT(OUT)::ierror
	logical, intent(in) :: wspeed_only
	integer, intent(in) :: iAngHi, iAngLow


	real :: diff_one_way, diff_other_way


		diff_one_way = library_sensor_zenith(iAngHi)- solarOrViewAng
		diff_other_way = solarOrViewAng - library_sensor_zenith(iAngLow)

		outEmissArray_ws(:,:,:,1) = (inEmissArray(iAngLow,:,:,:,1)  * diff_one_way + &
									 inEmissArray(iAngHi,:,:,:,1) * diff_other_way ) / dtheta
		outEmissArray_ws(:,:,:,2) = (inEmissArray(iAngLow,:,:,:,2)  * diff_one_way + &
									 inEmissArray(iAngHi,:,:,:,2) * diff_other_way) / dtheta
		outEmissArray_ws(:,:,:,3) = (inEmissArray(iAngLow,:,:,:,3)  * diff_one_way + &
									 inEmissArray(iAngHi,:,:,:,3) * diff_other_way) / dtheta

	outEmissArray = 0. 
	call interpolate_wind_speed(wspeed, outEmissArray_ws, outEmissArray)

	RETURN

	END SUBROUTINE InterpEmissForAGivenWSpeed
	

!	SUBROUTINE TriLinearInterpolationRefl(solarAng,viewAng,azmAng, &
!												solarMuArray, viewMuArray, relAzmArray, &
!												multiScatReflLUT_water, multiScatRefl_water, &
!												multiScatReflLUT_ice, multiScatRefl_ice,ierror)
!!................................................................................
!! !F90
!!
!! !DESCRIPTION:     
!	!perform 3-D linear interpolation on multiple scattering part of the reflectance. Instead of 
!	!indexing, a search is done to find the the indices of the cube that enclose the interpolation
!	!point.
!! !INPUT PARAMETERS:
!	!	solarAng		:	solar zenith angle in degrees (scalar)
!	!	viewAng 		:	sensor zenith angle in degrees (scalar)
!	!	azmAng			:	rel azm ang in degrees (scalar)]
!	!	solarMuArray	:	array of solar zenith angles (from the LUT) 1-D array
!	!	viewMuArray		:	array of sensor zenith angles (from the LUT) 1D array
!	!	relAzmArray		:	array of rel azimuth angles (from the LUT) 1-D array
!	!	multiScatReflLUT:	multi scat part of the reflectance (from LUT) 6D array
!	!						DIMS[1:nview, 1:nsolar,1:nazm, 1:ntau, 1:nlambda, 1:nre]
!!
!!
!! !OUTPUT PARAMETERS:
!	!	multiScatRefl	:	interpolated MS reflectance for the input geometry 3D array
!	!						DIMS[1:ntau, 1:nlambda, 1:nre]
!	!	ierror			:	integer value
!	!						ierror 	= 0  	success, 
!	!								= 1 	I/O error
!	!								= 2 	memory allocation error
!	!								= 3		invalid data	
!!................................................................................
!! !Revision History:
!!		Revision 1.0  2009/06/03 GW
!!		Initial integration of the standalone version. No changes from the standalone.
!!		
!!................................................................................
!! !Team-unique Header:
!!             Cloud Retrieval Group, NASA/GSFC
!!
!! !PROGRAMMER:
!!
!!             Nandana Amarasinghe (SSAI)
!!             Climate and Radiation Branch
!!             NASA Goddard Space Flight Center
!!             Greenbelt, Maryland 20771, U.S.A.
!!
!!             Gala Wind (SSAI)
!!             Climate and Radiation Branch
!!             NASA Goddard Space Flight Center
!!             Greenbelt, Maryland 20771, U.S.A.
!!
!!*******************************************************************************
!! !END
!	use libraryarrays
!
!    IMPLICIT NONE
!
!	REAL,INTENT(IN):: solarAng, viewAng, azmAng
!	REAL,INTENT(IN):: solarMuArray(:), viewMuArray(:), relAzmArray(:)
!	REAL,INTENT(IN):: multiScatReflLUT_water(:,:,:,:,:,:)
!	REAL,INTENT(IN):: multiScatReflLUT_ice(:,:,:,:,:,:)
!
!	REAL,INTENT(OUT)::multiScatRefl_water(:,:,:)
!	REAL,INTENT(OUT)::multiScatRefl_ice(:,:,:)
!	INTEGER,INTENT(OUT)::ierror
!
!
!	REAL*8::prod11,prod12,prod21,prod22,vol1,vol2,vol3,vol4,vol5,vol6,vol7,vol8,volume
!	REAL*4::solarmu,viewmu
!	REAL*8::dmu0,dmu,dmu0_1,dmu0_2,dmu_1,dmu_2,dazm,dazm_1,dazm_2
!	INTEGER:: ishi,islow, ivhi,ivlow, iazhi,iazlow,n1,n2,n3
!
!	n1 = number_solarzenith !size(solarMuArray)
!	n2 = number_sensorzenith !size(viewMuArray)
!	n3 = number_relazimuth !size(relAzmArray)
!	
!	
!	ierror = 0   !no error
!	solarmu  = solarAng !COS(solarAng * RAD_PER_DEG)
!	viewmu   = viewAng !COS(viewAng  * RAD_PER_DEG)
!	
!	CALL BisectionSearch(solarmu,solarMuArray(1:n1),n1,islow,ishi)
!	
!	CALL BisectionSearch(viewmu,viewMuArray(1:n2),n2,ivlow,ivhi)
!	
!	CALL BisectionSearch(azmAng,relAzmArray(1:n3),n3,iazlow,iazhi)
!	
!	dmu = solarMuArray(ishi) - solarMuArray(islow)
!	dmu0  = viewMuArray(ivhi)  - viewMuArray(ivlow)
!	dazm = relAzmArray(iazhi) - relAzmArray(iazlow)
!	
!	dmu_1 = solarMu - solarMuArray(islow) ; dmu_2 = solarMuArray(ishi) - solarMu
!	dmu0_1  = viewMu  - viewMuArray(ivlow)  ; dmu0_2  = viewMuArray(ivhi)  - viewMu
!	dazm_1 = azmAng  -  relAzmArray(iazlow); dazm_2 = relAzmArray(iazhi) - azmAng
!		
!	volume = dmu0 * dmu * dazm
!	IF(volume == 0.0)THEN
!		ierror = 3		!invalid data 
!		RETURN
!	ENDIF
!	!print*,ivhi,ishi,iazlow
!	!print*,dmu0_1,dmu0_2,dmu_1,dmu_2,dazm_1,dazm_2
!
!	prod22 = dmu_2 * dmu0_2
!	prod12 = dmu_1 * dmu0_2
!	prod11 = dmu_1 * dmu0_1
!	prod21 = dmu_2 * dmu0_1
!
!	vol8 =  abs(prod22 * dazm_2)
!	vol7 =  abs(prod22 * dazm_1)
!	vol6 =  abs(prod12 * dazm_1)
!	vol5 =  abs(prod12 * dazm_2)
!
!	vol1 =  abs(prod11 * dazm_1)
!	vol2 =  abs(prod21 * dazm_1)
!	vol3 =  abs(prod21 * dazm_2)
!	vol4 =  abs(prod11 * dazm_2)
!	!volume = vol1 + vol2 + vol3 + vol4 + vol5 + vol6 + vol7 + vol8
!		
!	
!	multiScatRefl_water(:,:,:) = ( multiScatReflLUT_water(ivlow,islow,iazlow,:,:,:) * vol8  + &
!							 multiScatReflLUT_water(ivhi,islow,iazlow,:,:,:)  * vol3  + &
!                       		 multiScatReflLUT_water(ivhi,islow,iazhi,:,:,:)   * vol2  + &
!                       		 multiScatReflLUT_water(ivlow,islow,iazhi,:,:,:)  * vol7  + &
!					   		 multiScatReflLUT_water(ivlow,ishi,iazlow,:,:,:)  * vol5  + &
!					   		 multiScatReflLUT_water(ivhi,ishi,iazlow,:,:,:)   * vol4  + &
!                       		 multiScatReflLUT_water(ivhi,ishi,iazhi,:,:,:)    * vol1  + &
!                       		 multiScatReflLUT_water(ivlow,ishi,iazhi,:,:,:)   * vol6 )/volume
!
!	multiScatRefl_ice(:,:,:) = ( multiScatReflLUT_ice(ivlow,islow,iazlow,:,:,:) * vol8  + &
!							 multiScatReflLUT_ice(ivhi,islow,iazlow,:,:,:)  * vol3  + &
!                       		 multiScatReflLUT_ice(ivhi,islow,iazhi,:,:,:)   * vol2  + &
!                       		 multiScatReflLUT_ice(ivlow,islow,iazhi,:,:,:)  * vol7  + &
!					   		 multiScatReflLUT_ice(ivlow,ishi,iazlow,:,:,:)  * vol5  + &
!					   		 multiScatReflLUT_ice(ivhi,ishi,iazlow,:,:,:)   * vol4  + &
!                       		 multiScatReflLUT_ice(ivhi,ishi,iazhi,:,:,:)    * vol1  + &
!                       		 multiScatReflLUT_ice(ivlow,ishi,iazhi,:,:,:)   * vol6 )/volume
!
!
! 
!	RETURN
!                  
!
!	CONTAINS
!
!		SUBROUTINE BisectionSearch(theta,thetaArray,ntheta,iAngLow, iAngHi)
!	
!		INTEGER,INTENT(IN)::ntheta
!		REAL,INTENT(IN)::theta,thetaArray(1:ntheta)
!		INTEGER,INTENT(OUT)::iAngLow, iAngHi
!	
!		INTEGER:: iAng
!	
!
!		iAngLow = 1
!		iAngHi = ntheta
!
!		if (ntheta == 1) then 
!			iAngLow = 1
!			iAngHi = 1
!			return
!		endif
!
!		DO
!			IF((iAngHi - iAngLow) == 1)EXIT
!			iAng = (iAngHi + iAngLow)/2
!			IF(thetaArray(iAng) < theta)THEN
!				iAngLow = iAng
!			ELSE
!				iAngHi = iAng
!			ENDIF
!		ENDDO
!
!		if (iAngLow < 1) iAngLow = 1
!		if (ianghi < 1) ianghi = 1
!		if (iAngHi > ntheta) iAngHi = ntheta
!		if (ianglow > ntheta) ianglow = ntheta
!
!	
!		RETURN
!	
!		END SUBROUTINE BisectionSearch
!	
!	END SUBROUTINE TriLinearInterpolationRefl


	SUBROUTINE GetRefl(solarAng,viewAng,azmAng,in_scat, refl_water, refl_ice, interp_MS, interp_SS, ierror)
!................................................................................
! !F90
!
! !DESCRIPTION:     
	!This subroutine outouts the Reflectance table (R(tau,lambda, re)) for a given sun-satellite
	!geometry 
! !INPUT PARAMETERS:
	!	solarAng	:	solar zenith in degrees (scalar)
	!	viewAng		:	viewing angle in degrees (scalar)
	!	azmAng		:	relative azimuth between sat and solar planes
!
!
! !OUTPUT PARAMETERS:
	!	refl	:	reflectance as a function of tau,lambda,re : 3-D Array
	!					DIM-1 = ntau, DIM-2 = nlambda, DIM-3 = nre
	!	ierror			:	integer value
	!						ierror 	= 0  	success, 
	!								= 1 	I/O error
	!								= 2 	memory allocation error
	!								= 3		invalid data	
!................................................................................
! !Revision History:
!		Revision 1.0  2009/06/03 GW
!		Initial integration of the standalone version. Changed from the standalone
!       delivered version to not apply the albedo inside this routine. That is done in a
!       different place later in corescience. 
!		
!................................................................................
! !Team-unique Header:
!             Cloud Retrieval Group, NASA/GSFC
!
! !PROGRAMMER:
!
!             Nandana Amarasinghe (SSAI)
!             Climate and Radiation Branch
!             NASA Goddard Space Flight Center
!             Greenbelt, Maryland 20771, U.S.A.
!
!             Gala Wind (SSAI)
!             Climate and Radiation Branch
!             NASA Goddard Space Flight Center
!             Greenbelt, Maryland 20771, U.S.A.
!
!*******************************************************************************
! !END

	use libraryarrays
	use GeneralAuxType
   ! WDR to get the scan line
   use ch_xfr, only : c2_scan
	!
	REAL, INTENT(IN) :: solarAng,viewAng,azmAng, in_scat
	REAL,INTENT(OUT) :: refl_water(:,:,:), refl_ice(:,:,:)
	INTEGER,INTENT(OUT)::ierror
	logical, intent(in) :: interp_MS, interp_SS
	
	INTEGER::allocateStatus
	REAL :: scat_angle
   INTEGER::do_mgr

   !WDR TEMP set some loop controls for refl table output
   !INTEGER::ksen,klam,ktau,kre

   scat_angle = in_scat ! scatAngle(solarAng,viewAng,azmAng)

   if (interp_MS) then

   !interpolate multi scattering part of the reflection
     multiScatRefl_water_lamb = 0
     multiScatRefl_ice_lamb = 0
     call mng_ms_get_lib( solarAng, viewAng, azmAng, 0, 0, c2_scan, &
       multiScatRefl_water_lamb, multiScatRefl_ice_lamb, ierror )
     if( ierror /= 0 ) return
	endif                                    

	if (interp_SS) then

		!get the phase function values
		CALL GetPhaseFunctionValues(scat_angle, phase_angles_water, number_phase_angles_water, phase_funcs_water, &
		                             phase_fun_norm_constant_water, phaseFunVals_water,ierror)		                                    
		CALL GetPhaseFunctionValues(scat_angle, phase_angles_ice, number_phase_angles_ice, phase_funcs_ice, &
		                             phase_fun_norm_constant_ice, phaseFunVals_ice,ierror)		                                    


!		IF(ierror /= 0)RETURN      !if there is any error no point of continuing
	
		!calculate the single scattering part of the reflection
		CALL CalcSingScatPartOfReflectance(phaseFunVals_water, phaseFunVals_ice, &
					 viewAng, solarang, ssRefl_water, ssRefl_ice, ierror)

	endif


!	IF(ierror /= 0)RETURN      !if there is any error no point of continuing                        
								    

	refl_water = ssRefl_water + multiScatRefl_water_lamb  !total reflection	with asurf = 0.0	
	refl_ice = ssRefl_ice + multiScatRefl_ice_lamb  !total reflection	with asurf = 0.0	

	RETURN
	END SUBROUTINE GetRefl



	!--------------------------SUBROUTINE GetReflForGivenWindSpeedWater-----------------------------
	!This subroutine outputs the Reflectance table (R(tau,lambda, re)) for a given sun-satellite
	!geometry and a wind speed
	!INPUTS	:
	!	solarAng	:	solar zenith in degrees (scalar)
	!	viewAng		:	viewing angle in degrees (scalar)
	!	azmAng		:	relative azimuth between sat and solar planes
	!	wspeed		:	wind speed m/s
	!OUTPUTS :
	!	reflAsurf	:	reflectance as a function of tau,lambda,re : 3-D Array
	!					DIM-1 = ntau, DIM-2 = nlambda, DIM-3 = nre
	!	ierror			:	integer value
	!						ierror 	= 0  	success, 
	!								= 1 	I/O error
	!								= 2 	memory allocation error
	!								= 3		invalid data	
	!ROUTINES CAlLED :
	!	ReadPhaseFunctionDataHDF
	!	ReadWaterLibraryHDF	
	!	getAllPhaseFunctionValues
	!	TriLinearInterpolationReflSearch
	!	CalcSingScatPartOfReflectance`
	!-----------------------------------------------------------------------------------------------
	SUBROUTINE GetReflForGivenWindSpeed(solarAng,viewAng,azmAng, in_scat, cos_scat, wspeed, reflAsurf_water, &
											reflAsurf_ice, wind_speed_only, interp_MS, interp_SS, ierror)


	use libraryarrays
	use libraryinterpolates
   use ch_xfr, only : c2_scan
	!
	REAL, INTENT(IN) :: solarAng,viewAng,azmAng,wspeed, in_scat, cos_scat
	REAL,INTENT(OUT) :: reflAsurf_water(:,:,:), reflAsurf_ice(:,:,:)
	INTEGER,INTENT(OUT)::ierror
	logical, intent(in) :: wind_speed_only, interp_MS, interp_SS
	
	INTEGER::allocateStatus,iuhi,iulow
	REAL :: scat_angle,deltau
   INTEGER::do_mgr
	
	scat_angle = in_scat !scatAngle(solarAng,viewAng,azmAng)

	if (interp_SS) then

		!get the phase function values
		call get_aero_params(cos_scat, aerosol_asym)
								                                
		CALL GetPhaseFunctionValues(scat_angle, phase_angles_water, number_phase_angles_water, phase_funcs_water, &
		                             phase_fun_norm_constant_water, phaseFunVals_water,ierror)		

		CALL GetPhaseFunctionValues(scat_angle, phase_angles_ice, number_phase_angles_ice, phase_funcs_ice, &
		                             phase_fun_norm_constant_ice, phaseFunVals_ice,ierror)		                                    

		!calculate the single scattering part of the reflection
		call Single_Scattering_Calcs_Ocean(phaseFunVals_water, phaseFunVals_ice, singlescattering_water, singlescattering_ice, &
				RLPhase, aeroPhase, rayleigh_tau, aerosol_tau, aerosol_ssa, &
					solarAng,viewAng, ssRefl_water, ssRefl_ice)

	endif


									 
	if (interp_MS) then 

!     ! WDR wire up the sfc typ 1 for 3 m/s
       multiScatRefl_water_CM(:,:,:,1) = 0
       multiScatRefl_ice_CM(:,:,:,1) = 0
       call mng_ms_get_lib( solarAng, viewAng, azmAng, 1, 0, c2_scan, &
          multiScatRefl_water_CM(:,:,:,1), multiScatRefl_ice_CM(:,:,:,1), &
          ierror )
       if( ierror /= 0 ) return

!
!     ! WDR wire up the sfc typ 1 for 7 m/s
        multiScatRefl_water_CM(:,:,:,2) = 0
        multiScatRefl_ice_CM(:,:,:,2) = 0
        call mng_ms_get_lib( solarAng, viewAng, azmAng, 2, 0, c2_scan, &
          multiScatRefl_water_CM(:,:,:,2), multiScatRefl_ice_CM(:,:,:,2), &
          ierror )
       if( ierror /= 0 ) return

!
!     ! WDR wire up the sfc typ 1 for 15 m/s
       multiScatRefl_water_CM(:,:,:,3) = 0
       multiScatRefl_ice_CM(:,:,:,3) = 0
       call mng_ms_get_lib( solarAng, viewAng, azmAng, 3, 0, c2_scan, &
         multiScatRefl_water_CM(:,:,:,3), multiScatRefl_ice_CM(:,:,:,3), &
         ierror )
       if( ierror /= 0 ) return
	endif


	int_reflectance_water_wspeed(:,:,:,1) = multiScatRefl_water_CM(:,:,:,1) + ssRefl_water(:,:,:)					
	int_reflectance_water_wspeed(:,:,:,2) = multiScatRefl_water_CM(:,:,:,2) + ssRefl_water(:,:,:)
	int_reflectance_water_wspeed(:,:,:,3) = multiScatRefl_water_CM(:,:,:,3) + ssRefl_water(:,:,:)					

	int_reflectance_ice_wspeed(:,:,:,1) = multiScatRefl_ice_CM(:,:,:,1) + ssRefl_ice(:,:,:)					
	int_reflectance_ice_wspeed(:,:,:,2) = multiScatRefl_ice_CM(:,:,:,2) + ssRefl_ice(:,:,:)					
	int_reflectance_ice_wspeed(:,:,:,3) = multiScatRefl_ice_CM(:,:,:,3) + ssRefl_ice(:,:,:)					

															
																
	call interpolate_wind_speed(wspeed, int_reflectance_water_wspeed, reflAsurf_water)
	call interpolate_wind_speed(wspeed, int_reflectance_ice_wspeed, reflAsurf_ice)

	END SUBROUTINE GetReflForGivenWindSpeed


	
	subroutine interpolate_wind_speed(wspeed, data_in, data_out)
	
      use ch_xfr, only : c2_scan
		real, intent(in) :: wspeed
		real, dimension(:,:,:,:), intent(in) :: data_in
		real, dimension(:,:,:), intent(inout) :: data_out
		
		real :: deltau
		integer :: iuhi, iulow

		IF( (wspeed - WIND_SPEED(1)) <=  0.2)THEN
	
			data_out(:,:,:) = data_in(:,:,:,1)
	
		
		ELSEIF( (wspeed - WIND_SPEED(3)) >=  -0.2)THEN
		
			data_out(:,:,:) = data_in(:,:,:,3) 
			

		ELSEIF(ABS(wspeed - WIND_SPEED(2)) <= 0.2)THEN


			data_out(:,:,:) = data_in(:,:,:,2) 

		ELSE
		
			IF(wspeed - WIND_SPEED(2) < 0.0)THEN
				deltau = WIND_SPEED(2) - WIND_SPEED(1)
				iuhi  = 2
				iulow = 1

				data_out(:,:,:) = ( data_in(:,:,:,1)  * (WIND_SPEED(iuhi)- wspeed) + &
                    		data_in(:,:,:,2)  * (wspeed - WIND_SPEED(iulow) ))/deltau

				
			ELSE

			deltau = WIND_SPEED(3) - WIND_SPEED(2)
			iuhi = 3
			iulow = 2

				data_out(:,:,:) = ( data_in(:,:,:,2)  * (WIND_SPEED(iuhi)- wspeed) + &
                    		data_in(:,:,:,3)  * (wspeed - WIND_SPEED(iulow) ))/deltau



			ENDIF

	   
		ENDIF


	end subroutine interpolate_wind_speed
	
	
 
	!--------------------------SUBROUTINE GetReflForGivenWindSpeedIce-------------------------------
	!This subroutine outputs the Reflectance table (R(tau,lambda, re)) for a given sun-satellite
	!geometry and a wind speed(m/s)
	!INPUTS	:
	!	solarAng	:	solar zenith in degrees (scalar)
	!	viewAng		:	viewing angle in degrees (scalar)
	!	azmAng		:	relative azimuth between sat and solar planes
	!	wspeed		:	wind speed m/s
	!OUTPUTS :
	!	reflAsurf	:	reflectance as a function of tau,lambda,re : 3-D Array
	!					DIM-1 = ntau, DIM-2 = nlambda, DIM-3 = nre
	!	ierror			:	integer value
	!						ierror 	= 0  	success, 
	!								= 1 	I/O error
	!								= 2 	memory allocation error
	!								= 3		invalid data	
	!ROUTINES CAlLED :
	!	ReadPhaseFunctionDataHDF
	!	ReadIceLibraryHDF	
	!	getAllPhaseFunctionValues
	!	TriLinearInterpolationReflSearch
	!	CalcSingScatPartOfReflectance
	!-----------------------------------------------------------------------------------------------
	SUBROUTINE GetSdevReflForGivenWindSpeed(solarAng,viewAng,azmAng,wspeed, &
     inrefl_ws_water, reflAsurf_water, inrefl_ws_ice, reflAsurf_ice, &
     wind_speed_only, ierror)

	use libraryarrays
    use ch_xfr, only : c2_scan
	!
	REAL, INTENT(IN) :: solarAng,viewAng,azmAng,wspeed
	real, intent(inout) :: inrefl_ws_water(:,:,:,:), inrefl_ws_ice(:,:,:,:)
	REAL,INTENT(OUT) :: reflAsurf_water(:,:,:), reflAsurf_ice(:,:,:)
	INTEGER,INTENT(OUT) :: ierror
	logical, intent(in) :: wind_speed_only
	
	INTEGER :: allocateStatus,iuhi,iulow
	REAL :: deltau
	integer :: num_radii
   INTEGER::do_mgr
	
	if (.not. wind_speed_only) then 

!
!     ! WDR wire up the sfc typ 1 for 3 m/s
       inrefl_ws_water(:,:,:,1) = 0
       inrefl_ws_ice(:,:,:,1) = 0
       call mng_ms_get_lib( solarAng, viewAng, azmAng, 1, 1, c2_scan, &
         inrefl_ws_water(:,:,:,1), inrefl_ws_ice(:,:,:,1), ierror )
       if( ierror /= 0 ) return

!
!     ! WDR wire up the sfc typ 2 for 7 m/s
       inrefl_ws_water(:,:,:,2) = 0
       inrefl_ws_ice(:,:,:,2) = 0
       call mng_ms_get_lib( solarAng, viewAng, azmAng, 2, 1, c2_scan, &
         inrefl_ws_water(:,:,:,2), inrefl_ws_ice(:,:,:,2), ierror )
       if( ierror /= 0 ) return

!
!     ! WDR wire up the sfc typ 3 for 15 m/s
       inrefl_ws_water(:,:,:,3) = 0
       inrefl_ws_ice(:,:,:,3) = 0
       call mng_ms_get_lib( solarAng, viewAng, azmAng, 3, 1, c2_scan, &
         inrefl_ws_water(:,:,:,3), inrefl_ws_ice(:,:,:,3), ierror )
       if( ierror /= 0 ) return

     !deallocate( wtrsd_int, icesd_int )
	
	endif
	
	
	call interpolate_wind_speed(wspeed, inrefl_ws_water, reflAsurf_water)
	call interpolate_wind_speed(wspeed, inrefl_ws_ice, reflAsurf_ice)


	RETURN

	END SUBROUTINE GetSdevReflForGivenWindSpeed
	!-----------------------------------------------------------------------------------------------
   SUBROUTINE GetSdevReflLamb(solarAng,viewAng,azmAng, reflAsurf_water, &
     reflAsurf_ice, ierror)

	   									
	use libraryarrays
   ! WDR access the scan line
   use ch_xfr, only : c2_scan
	!
	REAL, INTENT(IN) :: solarAng,viewAng,azmAng
	REAL,INTENT(OUT) :: reflAsurf_water(:,:,:), reflAsurf_ice(:,:,:)
	INTEGER,INTENT(OUT) :: ierror
	
	INTEGER :: allocateStatus,iuhi,iulow
	REAL :: deltau
	integer :: num_radii
   INTEGER::do_mgr

     reflAsurf_water = 0
     reflAsurf_ice = 0
     call mng_ms_get_lib( solarAng, viewAng, azmAng, 0, 1, c2_scan, &
       reflAsurf_water, reflAsurf_ice, ierror )
     if( ierror /= 0 ) return

!   deallocate( wtrsd_int, icesd_int )

	END SUBROUTINE GetSdevReflLamb

 end module interpolate_libraries
