module planck_functions

! This module replaces a bunch of F77 code that has been causing us grief since time
! immemorial. -- Gala Wind 3.30.17

	implicit none

	private 
	
	public :: modis_planck, modis_bright
	
!  Fundamental constants required for the monochromatic
!  Planck function routines PLANCK_M, PLANC_M, BRIGHT_M, BRITE_M
!  Taken from the NIST Reference on Constants, Units, and Uncertainty

!  http://physics.nist.gov/cuu/Constants/

!  See also:

!    Mohr, P.J. and B.N. Taylor, "CODATA recommended values of the
!    fundamental physical constants: 1998", Reviews of Modern Physics,
!    Vol.72, No.2, 2000.

! ... Planck constant (Joule second)
      real(8), parameter :: h = 6.62606876d-34

! ... Speed of light in vacuum (meters per second)
      real(8), parameter :: c = 2.99792458d+08

! ... Boltzmann constant (Joules per Kelvin)      
      real(8), parameter :: k = 1.3806503d-23

! ... Derived constants      
 	  real(8), parameter :: c1 = 2.0d+0 * h * c * c
 	  real(8), parameter :: c2 = (h * c) / k
	
!     MODIS BANDS
!      20,  21,  22,  23,
!      24,  25,  27,  28,
!      29,  30,  31,  32,
!      33,  34,  35,  36,
! ... Effective central wavenumbers (inverse centimeters)
      real, parameter :: cwn_terra(16) = (/ 2.641775E+03, 2.505277E+03, 2.518028E+03, 2.465428E+03, &
     										2.235815E+03, 2.200346E+03, 1.477967E+03, 1.362737E+03, &
											1.173190E+03, 1.027715E+03, 9.080884E+02, 8.315399E+02, &
											7.483394E+02, 7.308963E+02, 7.188681E+02, 7.045367E+02 /)
											
! ... Temperature correction slopes (no units)
      real, parameter :: tcs_terra(16) = (/ 9.993411E-01, 9.998646E-01, 9.998585E-01, 9.998682E-01, &
											9.998820E-01, 9.998845E-01, 9.994878E-01, 9.994918E-01, &
											9.995496E-01, 9.997399E-01, 9.995607E-01, 9.997256E-01, &
											9.999161E-01, 9.999167E-01, 9.999192E-01, 9.999282E-01 /)
											
! ... Temperature correction intercepts (Kelvin)
      real, parameter :: tci_terra(16) = (/ 4.770522E-01, 9.264053E-02, 9.756834E-02, 8.928794E-02, &
											7.309468E-02, 7.061112E-02, 2.204659E-01, 2.045902E-01, &
											1.599076E-01, 8.249532E-02, 1.302885E-01, 7.181662E-02, &
											1.970605E-02, 1.912743E-02, 1.816222E-02, 1.579983E-02 /)


      real, parameter ::  cwn_aqua(16) = (/ 2.647409E+03, 2.511760E+03, 2.517908E+03, 2.462442E+03, &
											2.248296E+03, 2.209547E+03, 1.474262E+03, 1.361626E+03, &
											1.169626E+03, 1.028740E+03, 9.076813E+02, 8.308411E+02, &
											7.482978E+02, 7.307766E+02, 7.182094E+02, 7.035007E+02 /)
											
      real, parameter ::  tcs_aqua(16) = (/ 9.993363E-01, 9.998626E-01, 9.998627E-01, 9.998707E-01, &
											9.998737E-01, 9.998770E-01, 9.995694E-01, 9.994867E-01, &
											9.995270E-01, 9.997382E-01, 9.995270E-01, 9.997271E-01, &
											9.999173E-01, 9.999070E-01, 9.999198E-01, 9.999233E-01 /)
											
      real, parameter ::  tci_aqua(16) = (/ 4.818401E-01, 9.426663E-02, 9.458604E-02, 8.736613E-02, &
											7.873285E-02, 7.550804E-02, 1.848769E-01, 2.064384E-01, &
											1.674982E-01, 8.304364E-02, 1.343433E-01, 7.135051E-02, &
											1.948513E-02, 2.131043E-02, 1.804156E-02, 1.683156E-02 /)

!-----------------------------------------------------------------------
!     MAS BANDS
!       30, 31, 42, 45, 46, 48, 49
      real, parameter ::  cwn_mas(7) = (/ 2652.519, 2551.02, 1170.411, 912.99, 830.77, 751.16, 723.96 /)
      real, parameter ::  tcs_mas(7) = (/ 0.99937, 0.99946, 0.99937, 0.99955, 0.99967, 0.99978, 0.99977 /)
      real, parameter ::  tci_mas(7) = (/ 0.44728, 0.43909, 0.22513, 0.13007, 0.087251, 0.05295, 0.053137 /)
      
!-----------------------------------------------------------------------
!     MASTER BANDS
!       30, 43, 47, 49 plus SHIS MAS-equivalent bands 48 and 49
      real, parameter ::  cwn_master(6) = (/ 2670.227, 1155.4015, 941.17, 820.51, 751.16, 723.96 /)
      real, parameter ::  tcs_master(6) = (/ 0.99944, 0.99934, 0.99937, 0.99963, 0.99978, 0.99977 /)
      real, parameter ::  tci_master(6) = (/ 0.41362, 0.23203, 0.18665, 0.094882, 0.05295, 0.053137 /)
 
!-----------------------------------------------------------------------
!     SEVIRI BANDS
!       3.9, 7.3, 8.5, 11., 12. and 13.4 um     
      real, parameter ::  cwn_seviri(6) = (/ 2569.094, 1362.142, 1149.083, 930.659, 839.661, 746.27 /)
      real, parameter ::  tcs_seviri(6) = (/ 0.9959, 0.9991, 0.9996, 0.9983, 0.9988, 0.9981 /)
      real, parameter ::  tci_seviri(6) = (/ 3.471, 0.485, 0.181, 0.627, 0.397, 0.576 /)
      
!-----------------------------------------------------------------------
!     VIIRS BANDS
!       3., 7.3, 8.5, 11., 12. and 13.4 um
      real, parameter ::  cwn_viirs(5) = (/ 2708.3865, 2460.8193, 1166.1845, 935.10476, 845.79752 /)
      real, parameter ::  tcs_viirs(5) = (/ 0.999354, 0.999623,  0.999698, 0.998273, 0.998778 /)
      real, parameter ::  tci_viirs(5) = (/ 0.593537, 0.325879, 0.146942, 0.650338,  0.421701 /)
      
!-----------------------------------------------------------------------
!     AMS BANDS
!       3.7 and 10.2
      real, parameter ::  cwn_ams(2) = (/ 2656.04241, 979.9118 /)
      real, parameter ::  tcs_ams(2) = (/ 0.99942, 0.99778 /)
      real, parameter ::  tci_ams(2) = (/ 0.42667, 0.65987 /)
      
!-----------------------------------------------------------------------
!     EMAS BANDS
!       3.7, 6.7, 7.3, 8.5, 11, 12, 13.3, 13.6 and 13.9
!        26,  27, 28,  30, 33, 34,   36,   37,      38
      real, parameter ::  cwn_emas(9) = (/ 2683.8433, 1492.5373, 1372.4951, 1167.6786, 903.9956, &
                                           832.9168,  749.9062,  733.5681, 716.1784  /)
                                                                                        
      real, parameter ::  tcs_emas(9) = (/  0.99932, 0.99938, 0.99925, 0.99977, 0.99985, &
                                            0.99974, 0.99985, 0.99993, 0.99990 /)
                                                                                        
      real, parameter ::  tci_emas(9) = (/  0.50705, 0.27642, 0.32168, 0.082686, 0.041587, &
                                            0.069156, 0.036514, 0.017001, 0.022627 /)
          


!------------------------------------------------------------------------
!	ASTER bands
!       10, 11, 12, 13 and 14
      real, parameter ::  cwn_aster(5) = (/ 1208.1301, 1159.0481, 1102.0475, 939.39813, 886.77863 /)
      real, parameter ::  tcs_aster(5) = (/ 0.99967480, 0.99970657, 0.99965727, 0.99928176, 0.99941397 /)
      real, parameter ::  tci_aster(5) = (/ 0.12009826, 0.10460924, 0.11673445, 0.21520464, 0.16743699 /)
      
!------------------------------------------------------------------------
!	AHI bands
!       7-16
      real, parameter ::  cwn_ahi(10) = (/ 2575.7673, 1609.2411, 1442.0793, 1361.3868, &
      									   1164.4431, 1038.1084,  961.3326,  890.7408, &
      									   809.2418,   753.3690 /)
      									   
      real, parameter ::  tcs_ahi(10) = (/ 0.99936, 0.99649, 0.99927, 0.99986, &
      									   0.99963, 0.99971, 0.99971, 0.99938, &
      									   0.99908, 0.99975 /)
      									   
      real, parameter ::  tci_ahi(10) = (/ 0.45925, 1.62641, 0.30465, 0.05552, &
      									   0.13189, 0.09194, 0.08686, 0.17399, &
      									   0.23460, 0.05990 /)
	
contains

	  subroutine pick_a_band(platform_name, band, cwn, tcs, tci, cwn_array, tcs_array, tci_array)
	  
	  
	    character*(*), intent(in) :: platform_name
	    integer, intent(in) :: band
	    real, intent(inout) :: cwn, tcs, tci
		real, dimension(:), intent(in), optional :: cwn_array, tcs_array, tci_array	    

	    integer :: index  
	      
! if S-HIS data is present, we'll use defaults
! Aircraft CWN, TCS and TCI can change randomly so we have defaults, but 
! will probably never use them. 

		if (platform_name(1:6) == 'Master' .or. &
			platform_name(1:6) == 'master' .or. &
			platform_name(1:6) == 'MASTER') then

			if (band == 30) index = 1
			if (band == 43) index = 2
			if (band == 47) index = 3
			if (band == 49) index = 4
			if (band == 60) index = 5
			if (band == 61) index = 6

			if (present(cwn_array) .and. band < 60) then 
	        	cwn = cwn_array(band)
	        	tcs = tcs_array(band)
	        	tci = tci_array(band)
			else
	        	cwn = cwn_master(index)
	        	tcs = tcs_master(index)
	        	tci = tci_master(index)
			endif

		else if (platform_name(1:3) == 'Mas' .or. & 
        		platform_name(1:3) == 'mas' .or. &
              	platform_name(1:3) == 'MAS') then

			if (band == 30) index = 1
			if (band == 31) index = 2
			if (band == 42) index = 3
			if (band == 45) index = 4
			if (band == 46) index = 5
			if (band == 48 .or. band == 60) index = 6
			if (band == 49 .or. band == 61) index = 7

			if (present(cwn_array) .and. band < 60) then 
	        	cwn = cwn_array(band)
	        	tcs = tcs_array(band)
	        	tci = tci_array(band)
			else
	        	cwn = cwn_mas(index)
	        	tcs = tcs_mas(index)
	        	tci = tci_mas(index)
			endif

		else if (platform_name(1:6) == 'SEVIRI' .or. &
              	platform_name(1:6) == 'seviri' .or. &
              	platform_name(1:6) == 'Seviri') then

			if (band == 4) index = 1
			if (band == 6) index = 2
			if (band == 7) index = 3
			if (band == 9) index = 4
			if (band == 10) index = 5
			if (band == 11) index = 6

	        cwn = cwn_seviri(index)
	        tcs = tcs_seviri(index)
	        tci = tci_seviri(index)

		else if (platform_name(1:11) == 'NPP_:_VIIRS' .or. &
     			platform_name(1:3) == 'npp' .or.  &
     			platform_name(1:5) == 'VIIRS' ) then
			index = band - 11

	        cwn = cwn_viirs(index)
	        tcs = tcs_viirs(index)
    	    tci = tci_viirs(index)
	
		else if (platform_name(1:3) == 'AMS') then 
			if (band == 11) index = 1
			if (band == 12) index = 2

	      	cwn = cwn_ams(index)
	      	tcs = tcs_ams(index)
	      	tci = tci_ams(index)

		else if (platform_name(1:4) == 'EMAS') then 
		
			if (band == 26) index = 1
			if (band == 27) index = 2
			if (band == 28) index = 3
			if (band == 30) index = 4
			if (band == 33) index = 5
			if (band == 34) index = 6
			if (band == 36) index = 7
			if (band == 37) index = 8
			if (band == 38) index = 9

			if (present(cwn_array)) then 
	        	cwn = cwn_array(band)
	        	tcs = tcs_array(band)
	        	tci = tci_array(band)
			else
		      	cwn = cwn_emas(index)
		      	tcs = tcs_emas(index)
		      	tci = tci_emas(index)
			endif
			
		else if (platform_name(1:5) == 'ASTER') then 

			index = band - 9

	      	cwn = cwn_aster(index)
	      	tcs = tcs_aster(index)
	      	tci = tci_aster(index)

		else if (platform_name(1:3) == 'AHI') then 

			index = band-6

	      	cwn = cwn_ahi(index)
	      	tcs = tcs_ahi(index)
	      	tci = tci_ahi(index)      
		
		else 
! MODIS
! ... Get index into coefficient arrays
      		if (band .le. 25) then
        		index = band - 19
      		else
      		 	index = band - 20
      		endif
      		
      		if (platform_name(1:5) == 'Terra' .or. &
				platform_name(1:5) == 'terra' .or. &
				platform_name(1:5) == 'TERRA') then
        			cwn = cwn_terra(index)
        			tcs = tcs_terra(index)
        			tci = tci_terra(index)
        	else if (platform_name(1:4) == 'Aqua' .or. &
            	platform_name(1:4) == 'aqua' .or. &
            	platform_name(1:4) == 'AQUA') then
        			cwn = cwn_aqua(index)
       				tcs = tcs_aqua(index)
        			tci = tci_aqua(index)
        	else
        		print*, "Invalid platform: ", trim(platform_name)
        		stop
        	endif

		endif

	  
	  end subroutine pick_a_band



	REAL FUNCTION PLANCK_M(W, T)

!-----------------------------------------------------------------------
!!F77
!
!!DESCRIPTION:
!    Compute monochromatic Planck radiance given brightness temperature
!    (Radiance units: Watts per square meter per steradian per micron)
!
!!INPUT PARAMETERS:
!    W (REAL)           Wavelength (microns)
!    T (REAL)           Brightness temperature (Kelvin)
!
!!OUTPUT PARAMETERS:
!    PLANCK_M (REAL)    Monochromatic Planck radiance (Watts per
!                       square meter per steradian per micron)
!
!!REVISION HISTORY:
!
!!TEAM-UNIQUE HEADER:
!    Liam.Gumley@ssec.wisc.edu
!
!!END
!-----------------------------------------------------------------------

		real, intent(in) :: w, t
		real(8) :: ws

      	planck_m = -1.0
      
      	if (w <= 0.0 .or. t <= 0.0) return
                  
! ... Convert wavelength to meters
      	ws = w * 1.0d-6
! ... Compute Planck radiance
      	planck_m = sngl(1.0d-6 * (c1 / ws**5) / (exp(c2 / (ws * dble(t))) - 1.0d+0))

	END function planck_m
      
! same as planck_m, but the wavelength is in wavenumber instead of micron
    REAL FUNCTION PLANC_M(V, T)

      	real, intent(in) :: v, t
      	real(8) :: vs

      	planc_m = -1.0
      
      	if (v <= 0.0 .or. t <= 0.0) return
                  
! ... Convert wavenumber to inverse meters
      	vs = 1.0d+2 * v      
! ... Compute Planck radiance
      	planc_m = sngl(1.0d+5 * (c1 * vs**3) / (exp(c2 * vs / dble(t)) - 1.0d+0))
            
	END function planc_m


    REAL FUNCTION BRIGHT_M(W, R)

!-----------------------------------------------------------------------
!!F77
!
!!DESCRIPTION:
!    Compute brightness temperature given monochromatic Planck radiance
!    (Radiance units: Watts per square meter per steradian per micron)
!
!!INPUT PARAMETERS:
!    W (REAL)           Wavelength (microns)
!    R (REAL)           Monochromatic Planck radiance (Watts per
!                       square meter per steradian per micron)
!
!!OUTPUT PARAMETERS:
!    BRIGHT_M (REAL)    Brightness temperature (Kelvin)
!
!!REVISION HISTORY:
!
!!TEAM-UNIQUE HEADER:
!    Liam.Gumley@ssec.wisc.edu
!
!!END
!-----------------------------------------------------------------------

      	real, intent(in) :: w, r

      	real(8) :: ws

      	bright_m = -1.0
      
      	if (w <= 0.0 .or. r <= 0.0) return
                  
! ... Convert wavelength to meters
      	ws = w * 1.0d-6
! ... Compute brightness temperature
		bright_m = sngl(c2 / (ws * log(c1 / (1.0d+6 * dble(r) * ws**5) + 1.0d+0)))

	END function bright_m
     
! same as bright_m, only the wavelength is in wavenumber instead of microns      
	REAL FUNCTION BRITE_M(V, R)

    	real, intent(in) :: v, r

      	real(8) :: vs

      	brite_m = -1.0
      
      	if (v <= 0.0 .or. r <= 0.0) return
                  
! ... Convert wavenumber to inverse meters
      	vs = 1.0d+2 * v           
! ... Compute brightness temperature
      	brite_m = sngl(c2 * vs / log(c1 * vs**3 / (1.0d-5 * dble(r)) + 1.0d+0))
      
      END function brite_m



!-----------------------------------------------------------------------
!!F77
!
!!DESCRIPTION:
!    Compute brightness temperature for a MODIS infrared band
!    on Terra or Aqua.
!
!    Spectral responses for each IR detector were obtained from MCST:
!    ftp://mcstftp.gsfc.nasa.gov/pub/permanent/MCST in
!    directories PFM_L1B_LUT_4-30-99 and FM1_RSR_LUT_07-10-01
!
!    An average spectral response for each infrared band was
!    computed. The band-averaged spectral response data were used
!    to compute the effective central wavenumbers and temperature
!    correction coefficients included in this module.
!
!    NOTE: The plaform name ('Terra' or 'Aqua') is passed to this
!    function via the common block defined in 'platform_name.inc'.
!
!!INPUT PARAMETERS:
!    RAD (REAL)      Planck radiance (units are determined by UNITS)
!    BAND (LONG)     MODIS IR band number (20-25, 27-36)
!    UNITS (LONG)    Flag defining radiance units
!                    0 => milliWatts per square meter per
!                         steradian per inverse centimeter
!                    1 => Watts per square meter per
!                         steradian per micron
!
!!OUTPUT PARAMETERS:
!    MODIS_BRIGHT  Brightness temperature (Kelvin)
!                  Note that a value of -1.0 is returned if
!                  RAD .LE. 0.0, or BAND is not in range 20-25, 27-36.
!
!!REVISION HISTORY:
!    Liam.Gumley@ssec.wisc.edu
!
!!TEAM-UNIQUE HEADER:
!    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
!
!!END
!-----------------------------------------------------------------------


	REAL FUNCTION MODIS_BRIGHT(platform_name,RAD, BAND, UNITS, cwn_array, tcs_array, tci_array)


    	real, intent(in) :: rad
      	integer, intent(in) :: band, units
      	character*(*), intent(in) :: platform_name
		real, dimension(:), intent(in), optional :: cwn_array, tcs_array, tci_array

      	real :: cwn, tcs, tci

	    modis_bright = -1.0

		if (present(cwn_array)) then
			call pick_a_band(platform_name, band, cwn, tcs, tci, cwn_array=cwn_array, tcs_array=tcs_array, &	
															tci_array=tci_array)
		else
			call pick_a_band(platform_name, band, cwn, tcs, tci, cwn_array=cwn_array, tcs_array=tcs_array, &	
															tci_array=tci_array)
		endif	

! ... Compute brightness temperature
	    if (units == 1) then
! ...   Radiance units are
! ...   Watts per square meter per steradian per micron
    		modis_bright = (bright_m(1.0e+4 / cwn, rad) - tci) / tcs
      	else
! ...   Radiance units are
! ...   milliWatts per square meter per steradian per wavenumber
        	modis_bright = (brite_m(cwn, rad) - tci) / tcs
      	endif

		if (modis_bright == -1.) then 
			print*, "ERROR: Brightness temperature calculation: ", trim(platform_name), band, rad
		endif

	end function modis_bright

! same thing except it will return the planck radiance 

    REAL FUNCTION MODIS_PLANCK(platform_name, TEMP, BAND, UNITS)

    	real, intent(in) :: temp
      	integer, intent(in) :: band, units
      	character*(*), intent(in) :: platform_name

      	real :: cwn, tcs, tci

	    modis_planck = -1.0

		call pick_a_band(platform_name, band, cwn, tcs, tci)

! ... Compute Planck radiance

      	if (units == 1) then
! ...   Radiance units are
! ...   Watts per square meter per steradian per micron
      	  modis_planck = planck_m(1.0e+4 / cwn, temp * tcs + tci)
      	else
! ...   Radiance units are
! ...   milliWatts per square meter per steradian per wavenumber
      	  modis_planck = planc_m(cwn, temp * tcs + tci)
      	endif

		if (modis_planck == -1.) then 
			print*, "ERROR: Planck radiance calculation: ", trim(platform_name), band, temp
		endif
		
	end function modis_planck

end module planck_functions


