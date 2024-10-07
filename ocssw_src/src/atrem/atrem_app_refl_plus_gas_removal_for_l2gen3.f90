!------------------------------------------------------------------------------|
! The September 2012 version of atrem code, 'atrem_for_PRISM.f' (atrem for JPL |
!     PRISM Imaging Spectrometer, was obtained through modifying a 2011        |
!     version of 'atrem_f90_cubeio.f'.                                         |
! The March 2013 version of atrem code, 'atrem_for_PRISM_2013.f', was updated  |
!     from the September 2012 version of atrem code. The upgrades include:     |
!    a) replacing the MODTRAN 3.5 solar irradiance curve with a new solar curve|
!       constructed from the solar irradiance curve of Thuillier et al. 2004   |
!       ATLAS 3 (<=644.7 nm) and that of Kurucz 2005 (> 644.7 nm) built in     |
!       MODTRAN 5.2. This newly contructed solar curve is more suited for      |
!       modeling imaging imaging spectrometer data below 450 nm and with a     |
!       spectral resolution finer than 5 nm;                                   |
!    b) replacing high spectral resolution gas absorption tables with new      |
!       tables calculated with the HITRAN2008+ line database and a fast        |
!       line-by-line code originally developed by Bill Ridgway. After such     |
!       upgrades, the retrieved surface reflectances near 1.45 and 1.95 micron |
!       are improved significantly;                                            |
!    c) modifying volume mixing ratios of atmospheric carbon dioxide;          |
!    d) applying a scaling factor for improved modeling of atmospheric         |
!       oxygen band absorption centered near 1.265 micron.                     |
!------------------------------------------------------------------------------|
! The 2001 version of ATREM code uses the new 'cubeio.f90' for I/O operations. |
!     The output data cubes have no headers, unlike the previous generations   |
!     of ATREM codes containing 512 bytes of headers in output data cubes.     |
!------------------------------------------------------------------------------|
!********************************************************************************
!*									       *
!* Name:           atrem_app_refl_plus_gas_removal_2013         April,  2013    *
!*                     (ATmosphere REMoval Program)                	       *
!*									       *
!* Author: Bo-Cai Gao, Remote Sensing Division, Code 7230,                      *
!*                     Naval Research Laboratory, Washington, DC 20375 USA      *
!*									       *
!* Purpose: To derive surface reflectances from spectral imaging data collected *
!*             by the Airborne Visible Infrared Imaging Spectrometer (AVIRIS),  *
!*             HYDICE, HSI, PHILLS, HICO, NIS, PRISM, and possibly other        *
!*             airborne and spaceborne imaging spectrometers, and to derive a   *
!*             column water vapor image for each scene.                         *
!* 									       *
!* General Principles: In order to derive surface reflectances from image data  *
!*             cube, a thorough compensation for the atmospheric absorption and *
!*             scattering is required. The spatial and temporal variations of   *
!*             atmospheric water vapor amounts pose difficulties in removing    *
!*             water vapor absorption features in data cube using standard      *
!*             atmospheric models. In this algorithm, the amount of water vapor *
!*             on a pixel-by-pixel basis is derived from data cube themselves   *
!*             using the 0.94- and the 1.14-um water vapor bands and a          *
!*             three-channel ratioing technique. The derived water vapor        *
!*             values are then used for modeling water vapor absorption effects *
!*             in the entire 0.4-2.5 um region. The absorption effects of       *
!*             well mixed atmospheric gases, such as CO2, N2O, CH4, and O2,     *
!*             are modeled using standard atmospheric models. A line-by-line    *
!*             code developed by W. Ridgway at NASA Goddard Space Flight Center *
!*             is used to generate a database of gaseous absorption             *
!*             coefficients at high spectral resolution. This database and      *
!*             ozone cross sections and NO2 cross sections are used in          *
!*             calculating transmittances of eight atmospheric gases (H2O, CO2  *
!*             O3, N2O, CO, CH4, O2, and NO2). The scattering effect is modeled *
!*             using a modified version of 6S code (Vermote et al. 1996).       *
!*									       *
!* Algorithm:  1. An input parameter file is read in and global data are        *
!*		 initialized. 						       *
!*	      2. The solar zenith angle is derived based on the flight time    *
!*		 and latitude and longitude of the scene.		       *
!*	      3. A table consisting of transmittance spectra at the solar      *
!*		 and observational geometry but with varying amounts of water  *
!*	         vapor is created. 				               *
!*	         The transmittance spectra for water vapor (H2O), carbon       *
!*		 dioxide (CO2), ozone (O3), nitrous oxide (N2O), carbon        *
!*	         monoxide (CO), methane (CH4), oxygen (O2), and nitrogen       *
!*                dioxide (NO2) in the 0.4-2.5 micron region are calculated.    *
!*             4. Aerosol and molecular scattering terms are calculated using   *
!*		 a modified version of the 6S code.			       *
!*             5. The image data cube (2D spatial and 1D spectral) are accessed *
!*	         one spectral slice at a time from an image file in any storage*
!*	         order, Band SeQuential (BSQ), Band Interleaved by Pixel (BIP),*
!*                and Band Interleaved by Line (BIL). Each measured radiance    *
!*	         spectrum is divided by the solar irradiance curve above the   *
!*		 atmosphere to get apparent reflectances.  		       *
!*             6. The integrated water vapor amount for each apparent           *
!*		 reflectance spectrum is derived from the 0.94- and the 1.14-um*
!*	         water vapor absorption bands using a 3-channel ratioing       *
!*		 technique and a look-up table procedure. 		       *
!*	      7. The surface reflectances are derived using another look-up    *
!*		 table procedure. Routines written in C language are used to   *
!*                perform all the input/output operations.		       *
!*									       *
!*  Limitations:						         	       *
!*	  1.  This algorithm works best for areas where the surface elevation  *
!*             does not vary by more than 1 km.				       *
!*	  2.  The algorithm works well only for measurements taken on clear    *
!*             days. For hazy days with large aerosol optical depths, the       *
!*             algorithm does not work nearly as well.                          *
!*	  3.  The algorithm does not work well over dark areas such as rivers  *
!*	      and lakes. It does not perform corrections for "atmospheric      *
!*             adjacency" effects.					       *
!*	  4.  The algorithm assumes horizontal surfaces having Lambertian      *
!*             reflectances.                                       	       *
!*         5.  Additional work is needed to port this algorithm to different    *
!*             computing systems. This algorithm was developed on an SGI        *
!*             workstation and used FORTRAN statements to read several binary   *
!*             files. The binary file format on different computers can be      *
!*             different. The record lengths defined on an SGI machine are      *
!*             different from most other computers. Some of the Dec Alpha       *
!*             workstations use 8 bytes (64 bits) to store a C-pointer, instead *
!*             of 4 bytes (32 bits). All these factors need to be taken into    *
!*             considerations when porting this algorithm to other computers.   *
!*									       *
!*  User Input: 								       *
!*    1.  Name - Imaging spectrometer name, e.g., AVIRIS, HYDICE, HSI, PHYLLS.  *
!*    2.  Plane altitude in km, above the sea level.                            *
!*    3.  Date - The month, day, and year that the data were collected.         *
!*    4.  Time - The hour, minute, second (GMT) that the data were	       *
!*                 collected.                                                   *
!*    5.  Latitude - The mean latitude of the measurement area in degrees,      *
!*           minutes, seconds, and hemisphere.                          	       *
!*    6.  Longitude - The mean longitude of the measurement area in             *
!*           degrees, minutes, seconds, and hemisphere.                 	       *
!*    7.  View zenith angle - the angle  in degrees, minutes, and seconds       *
!*                            between the viewing direction and nadir           *
!*    8.  View azimuth angle in degree, minutes and seconds. The view azimuth   *
!*                            angle is defined as the angle between local North *
!*                            the view direction projected on the ground        *
!*                            (clockwise count). This angle can vary from 0     *
!*                            to 360 deg. Ths definition here is consistent with*
!*                            navigator's convention, but not consistent with   *
!*                            astronomer's convention. Here is an example: for  *
!*                            a downward looking instrument pointed off-nadir   *
!*                            to the West, the view azimuthal angle is 90 deg.  *
!*                            Conversely, if the instrument is pointed toward   *
!*                            East, the view azimuthal angle is 270 deg.        *
!*									       *
!*    9.  Resolution of input spectra in nm (~ 10 nm for AVIRIS) or "0"         *
!*        if the FWHM values (in micron) are present in the wavelength file     *
!*   10.  Filename of spectrometer wavelength table - 2 or 3 column ASCII data. *
!*           The first column is the channel number, the second column is       *
!*           the wavelength and the optional third column is the FWHM value for *
!*           each channel.   	      					       *
!*   11.  Channel ratio parameters - Center positions and widths of window and  *
!*           water vapor channels used in channel ratios for deriving   	       *
!*	    the amount of water vapor from imaging data. Each of the window    *
!*           and water vapor channels typically consists of several narrow      *
!*           imaging spectrometer channels.                        	       *
!*   12.  Atmospheric model - Temperature, pressure, water vapor volume mixing  *
!*           ratio profiles. The model can be either a predefined model or a    *
!*           user defined model.          				       *
!*   13.  Gases - An indication of which of the 8 gases (H2O, CO2, O3,          *
!*	    N2O, CO, CH4, O2, and NO2) should be included in the spectral      *
!* 	    calculations.                  		               	       *
!*   14. Ozone - Column amount of O3 - changes with season and                  *
!*	    latitude.  Typical value is .34 atm-cm.                    	       *
!*   15. Aerosol model number and visibility when measurements were taken.      *
!*   16. Aerosol optical depth at 550 nm (Optional, only if visibility V        *
!*           is set to 0)						       *
!*   17. Surface elevation - The average surface elevation of an imaging scene  *
!*           in km.                                                             *
!*   18. Filename of input image file 					       *
!*   19. Dimensions of input image.					       *
!*   20. Filename of output image file in same storage order as input file      *
!*   21. Resolution of output surface reflectance spectra in micron             *
!*   22. Scale factor for output reflectance values                             *
!*   23. Filename of output water vapor image.				       *
!*   24. Filename of output atmospheric transmittance look-up table             *
!*									       *
!*  Output:  								       *
!*     a. Output to user specified files: 				       *
!*        1. Surface reflectance cube data - retrieved from the image data cube *
!*        2. Water vapor image - an image of the derived column water vapor     *
!*           amounts at each pixel.					       *
!*        3. Transmittance Look-up table - consisting of channel ratio values   *
!*           for the .94- and 1.14-um channels corresponding to 60 water vapor  *
!*           amounts, and the associated atmospheric transmittance spectra.     *
!*     b. Output to standard output:                      	       	       *
!*        1. Debugging information - if the source is compiled with the         *
!*           debug flag set (-d_lines for the ULTRIX compiler),         	       *
!*           then debug information is written out.               	       *
!*        2. Error messages - if any of the user input is invalid, or there is  *
!*           an I/O error, a message is written out, and the program halts.     *
!*        3. Progress indicator - a message, which shows the number of the      *
!*           spectral slices that have been processed, to give the user an      *
!*           indication of the progress of the program execution.               *
!*									       *
!*  Special Considerations: 						       *
!*     a. Make sure that the imaging spectrometer's channel positions specified *
!*           in the input wavelength table are correct. The incorrect channel   *
!*           positions will introduce sharp features in derived surface         *
!*           reflectance spectra.                                               *
!*     b. The output reflectance cube file has the same size as the input       *
!*           data cube file. Make sure there is enough space in the file        *
!*           system for the output cube before running this program.	       *
!*									       *
!*  Change History:							       *
!*     ATREM 1.2 - uses full-width half-max values to smooth solar irradiance   *
!*                 curve and a new scaling factor for methane		       *
!*     ATREM 1.3 - used new solar irradiance curve			       *
!*		- supports 1992 AVIRIS data				       *
!*     ATREM 1.3.1 - allows 0 and above for output resolution                   *
!*     ATREM 2.0 - allows variable scale factors for each spectrometer          *
!*                 (code submitted by Roger Clark, USGS)                        *
!*               - user can input output scale factor for reflectance values    *
!*               - user can input radiance file header size in bytes            *
!*               - new solar irradiance curve (Neckel and Labs plus ATMOS)      *
!*                 used (Green and Gao, 1993).                                  *
!*     ATREM 3.0 - replaced the 5S code by the 6S code for modeling atmospheric *
!*                 scattering effects and for modeling measurements from        *
!*                 low-altitude aircrafts. ATREM 3.0, like previous versions    *
!*                 of ATREM, only models nadir viewing geometry.                *
!*               - changed algorithms for calculating atmospheric gaseous       *
!*                 transmittances for the upward surface-sensor path            *
!*		- users need to specify aircraft altitude (km) in input file   *
!*     ATREM 4.0 - replaced the band model by a line-by-line-based algorithm    *
!*                 for calculating atmospheric gaseous transmittances. This     *
!*                 allows the modeling of imaging spectrometer data at          *
!*                 spectral resolutions better than 10 nm.                      *
!*               - increased the buffer size to 1024x1024 in order to handle    *
!*                 images larger than the AVIRIS' size (614x512).               *
!*               - modified the algorithm to allow off-nadir pointing geometry. *
!*               - users need to specify view zenith angle and view azimuth     *
!*                 angle. The view azimuth angle is defined as the angle        *
!*                 between local North and the view direction projected on the  *
!*                 ground (clockwise count). This angle can vary from 0 to      *
!*                 360 deg.                                                     *
!*               - replaced the solar irradiance curve (Neckel and Labs plus    *
!*                 ATMOS) by the high spectral resolution solar irradiance      *
!*                 curve from MODTRAN 3.5 released in December of 1996.         *
!*     ATREM 4.1 - included atmospheric NO2 in transmittance calculations.      *
!*									       *
!*  Acknowledgments:							       *
!*           This work was partially supported over several years by grants     *
!*           from NASA Jet Propulsion Laboratory, California Institute of       *
!*           Technology, from NASA Headquarters, and from the Office of Naval   *
!*           Research. Special thanks goes to Kathy Heidebrecht at the Center   *
!*           for the Study of Earth from Space, University of Colorado at       *
!*           Boulder, Colorado, for supporting the development of the earlier   *
!*           versions of ATREM codes. Special thanks also goes to William L.    *
!*           Ridgway for providing the line-by-line gaseous absorption database *
!*           used in the present version of ATREM.                              *
!*                     							       *
!*  References: 								       *
!*        Gao, B.-C., K. Heidebrecht, and A. F. H. Goetz, Derivation of scaled  *
!*           surface reflectances from AVIRIS data, Remote Sens. Env., 44,      *
!*           165-178, 1993.						       *
!*        Gao, B.-C., and C. O. Davis, Development of an operational algorithm  *
!*           for removing atmospheric effects from HYDICE and HSI data,         *
!*           in SPIE'96 Conference Proceedings, Vol. 2819, 45-55, 1996.         *
!*        Gao, B.-C., and A. F. H. Goetz, Column atmospheric water vapor and    *
!*           vegetation liquid water retrievals from airborne imaging	       *
!*           spectrometer data, J. Geophys. Res., 95, 3549-3564, 1990.	       *
!*        Goetz, A. F. H., and M. Herring, The high resolution imaging	       *
!*           spectrometer (HIRIS) for Eos, IEEE Trans. Geosci. Remote Sens., 27,*
!*           136-144, 1989.						       *
!*        Goetz, A. F. H., G. Vane, J. Solomon, and B. N. Rock, Imaging	       *
!*           spectrometry for Earth remote sensing, Science, 228, 1147-1153,1985*
!*        Green, R. O., and B.-C. Gao, A Proposed Update to the Solar Irradiance*
!*           Spectrum used in LOWTRAN and MODTRAN, in Summaries of the Fourth   *
!*           Annual JPL Airborne Geoscience Workshop, October 25-29, (Editor,   *
!*           R. O. Green), JPL Publ. 93-26, Vol. 1, pp. 81-84, Jet Propul. Lab, *
!*           Pasadena, Calif., 1993.					       *
!*        Kneizys, F. X., E. P. Shettle, L. W. Abreu, J. H. Chetwynd, G. P.     *
!*           Anderson, W. O. Gallery, J. E. A. Selby, and S. A. Clough, Users   *
!*           guide to LOWTRAN7, AFGL-TR-8-0177, Air Force Geophys. Lab.,        *
!*           Bedford, Mass., 1988.					       *
!*        Iqbal, M., An Introduction To Solar Radiation, Academic, San Diego,   *
!*           Calif., 1983.						       *
!*        Malkmus, W., Random Lorentz band model with exponential-tailed S line *
!*           intensity distribution function, J. Opt. Soc. Am., 57, 323-329,1967*
!*        Press, W. H., B. P. Flannery, S. A. Teukolsky, and W. T.  Vetterling, *
!*           Numerical Recipes-The ART of Scientific Computing, Cambridge       *
!*           University Press, 1986.					       *
!*        Rothman, L. S., et al., The HITRAN 2008 molecular spectroscopic       *
!*           database, JQSRT, 110, 533-572, 2009.                               *
!*        Solomon, S., R. W. Portmann, R. W. Sanders, J. S. Daniel, W. Madsen,  *
!*           B. Bartram, and E. G. Dutton, On the role of nitrogen dioxide in   *
!*           the absorption of solar radiation, J. Geophys. Res., 104,          *
!*           12047-12058, 1999.                                                 *
!*        Tanre, D., C. Deroo, P. Duhaut, M. Herman, J. J. Morcrette, J. Perbos,*
!*           and P. Y. Deschamps, Description of a computer code to simulate    *
!*           the satellite signal in the solar spectrum: the 5S code, Int.      *
!*           J. Remote Sens., 11, 659-668, 1990.				       *
!*        Tanre, D., C. Deroo, P. Duhaut, M. Herman, J. J. Morcrette, J. Perbos,*
!*           and P. Y. Deschamps, Simulation of the satellite signal in the     *
!*           solar spectrum (5S), Users' Guide (U. S. T. De Lille, 59655        *
!*           Villeneu d'Ascq, France: Laboratoire d'Optique Atmospherique),     *
!*	    1986. 							       *
!*        Thuillier, G., et al., Solar irradiance reference spectra for two     *
!*           solar active levels, Adv. Space Res., 34, 256-261, 2004.           *
!*        Vane, G., R. O. Green, T. G. Chrien, H. T. Enmark, E. G. Hansen, and  *
!*           W. M. Porter, The Airborne Visible/Infrared Imaging Spectrometer,  *
!*           Remote Sens. Env., 44, 127-143, 1993.                              *
!*        Vane, G. (Ed), Airborne visible/infrared imaging spectrometer	       *
!*	    (AVIRIS), JPL Publ. 87-38, Jet Propul. Lab, Pasadena, Calif., 1987.*
!*        Vermote, E., D. Tanre, J. L. Deuze, M. Herman, and J. J. Morcrette,   *
!*           Second simulation of the satellite signal in the solar spectrum    *
!*           (6S), 6S User's Guide Version 1, NASA-GSFC, Greenbelt, Maryland,   *
!*           134 pages, 1994.                                                   *
!*									       *
!********************************************************************************
!
!
!********************************************************************************
!*									       *
!*  Name: GET_INPUT							       *
!*  Purpose: Reads and verifies all user provided input.			       *
!*  Parameters: none							       *
!*  Algorithm: Data is read in from standard input. The user is not prompted    *
!*             interactively, so data should be redirected from an input file.  *
!*             The data is validated after it is read in.  If the data is       *
!*             invalid, a message is printed and the program stops. 	       *
!*  Globals used:   TPVMR - two dimensional array containing 6 predefined       *
!*		          atmospheric models.  The models contain values for   *
!*			  the number of atmospheric layer boundaries, altitude,*
!*			  temperature, pressure, and water vapor volume mixing *
!*			  ratio.					       *
!*  Global output:  IH2OVP, ICO2, IO3, IN2O, ICO, ICH4, IO2 - set to 0 to       *
!*                         indicate that the gas should NOT be used in the      *
!*                         calculations, and set to 1 if the gas should be used *
!*		   H(), T(), P(), VMR(), NB, NL, MODEL - altitude (km),        *
!*			  temperature (K), pressure (atm), water vapor volume  *
!*			  mixing ratio (ppm), number of atmospheric layer      *
!*			  boundaries, number of atmospheric layers (NB-1),     *
!*                         and model number for the selected atmospheric model. *
!*		   VRTO3 - column amount of O3.				       *
!*		   WAVOBS() - wavelengths for all channels.                    *
!*		   FWHM() -   resolutions for each channel.                    *
!*                  NOBS - number of AVIRIS wavelengths.                        *
!*                  HSURF - the mean surface elevation of the imaging scene     *
!*                  DLT, DLT2 - resolution, in units of nm, of input            *
!*                         spectra and resolution of output surface reflectance *
!*                         spectra. If DLT2>DLT, output spectra are smoothed    *
!*                         using a gaussian function.   			       *
!*	           WNDOW1, WNDOW2, WP94C - center wavelength positions of two  *
!*			  broad window channels and one broad .94-um water     *
!*			  vapor channel.				       *
!*		   WNDOW3, WNDOW4, W1P14C - center wavelength positions of two *
!*			  broad window channels and one broad 1.14-um water    *
!*			  vapor channel.  				       *
!*		   NB1, NB2, NBP94, NB3, NB4, NB1P14 - number of individual    *
!*			  narrow channels that form the corresponding broad    *
!*			  window and water vapor absorption channels.	       *
!*		   IMN, IDY, IYR, IH, IM, IS - month, day, year, hour, minute, *
!*                         and second of measurement.			       *
!*		   XLATD, XLATM, XLATS, LATHEM - degrees, minutes, seconds, and*
!*                         hemisphere of latitude of measured area.             *
!*		   XLONGD, XLONGM, XLONGS, LNGHEM - degrees, minutes, seconds, *
!*                         and hemisphere of latitude of measured area.         *
!*                  NAME_INSTRU: Imaging spectrometer name, e.g., AVIRIS, HYDICE*
!*                  XPSS:  = HSURF, an interface for using 6S.                  *
!*                  XPPP:  plane height (km, above the sea level).              *
!*                  							       *
!*  Return Codes: none							       *
!*  Special Considerations: None						       *
!*									       *
!********************************************************************************
      
      SUBROUTINE GET_INPUT

      use cubeio

      INCLUDE 'COMMONS_INC.f'

      INTEGER              LUN_IN, LUN_OUT, LUN_VAP, I_RET
      COMMON /INOUT_UNITS/ LUN_IN, LUN_OUT, LUN_VAP


!C  Common variables
      DIMENSION H(25), T(25), P(25), VMR(25)
      DIMENSION WAVOBS(1024),FWHM(1024)
      DIMENSION TPVMR(7,81)
      CHARACTER*1 LATHEM, LNGHEM

      CHARACTER (LEN = 1000) :: FINAV,FOCUB,FOH2O
      INTEGER SORDER,HDREC

!C  Local variables
      INTEGER ANS
      CHARACTER (LEN = 1000) :: FINPWV,FTPVMR      
      LOGICAL GOOD_DATA
      CHARACTER (LEN = 1000) :: FOUT1
      INTEGER DIMS(4)
      INTEGER ST_ORDER
      
      COMMON /GETINPUT1/ IH2OVP,ICO2,IO3,IN2O,ICO,ICH4,IO2,INO2
      COMMON /GETINPUT3/ H,T,P,VMR,NB,NL,MODEL,IAER,V,TAER55,VRTO3,SNO2
      COMMON /GETINPUT4/ WAVOBS,FWHM
      COMMON /GETINPUT5/ NOBS,HSURF,DLT,DLT2
      COMMON /GETINPUT6/ WNDOW1,WNDOW2,WP94C,WNDOW3,WNDOW4,W1P14C
      COMMON /GETINPUT7/ NB1,NB2,NBP94,NB3,NB4,NB1P14
      COMMON /GETINPUT8/ IMN,IDY,IYR,IH,IM,IS
      COMMON /GETINPUT9/ XLATD,XLATM,XLATS,LATHEM
      COMMON /GETINPUT10/XLONGD,XLONGM,XLONGS,LNGHEM
      COMMON /GETINPUT11/HDREC,NSAMPS,NLINES,NBANDS,SORDER
      COMMON /GETINPUT12/SCALEF
      COMMON /TPVMR_INIT1/ TPVMR

!C Commons for use with the C programs for cube I/O
      COMMON /OUTCUBE/ FOCUB
      COMMON /INCUBE/ FINAV
      COMMON /OUTH2OVAP/ FOH2O

!C Parameters for names of imaging spectrometers
      CHARACTER (LEN = 80) :: NAME_INSTRU, NAMES(10)
      COMMON /GETINPUT13/ NAME_INSTRU, NAMES

      REAL XPSS, XPPP
      COMMON /GETINPUT14/ XPSS, XPPP

      REAL XVIEWD, XVIEWM, XVIEWS
      REAL XAZMUD, XAZMUM, XAZMUS
      COMMON /GETINPUT15/ XVIEWD,XVIEWM,XVIEWS, XAZMUD,XAZMUM,XAZMUS

      GOOD_DATA = .TRUE.	! initialize flag for good data

!C  Get the name of imaging spectrometer, e.g., AVIRIS, HYDICE
      READ(5,627) NAME_INSTRU
 627  FORMAT(A10)
!C
!C
!C***Temp code for assigning names of imaging spectrometers. Based on these
!C    names, different scale factors should be used for different instrument
!C    when converting measured radiances to standard radiance units.
!C*** The coding here may need to be moved to the file "COMMONS_INC"
      NAMES(1) = 'AVIRIS'
      NAMES(2) = 'HYDICE'
      NAMES(3) = 'HSI'
      NAMES(4) = 'TRWIS-III'
      NAMES(5) = 'PHYLLS'
      NAMES(6) = 'Hyperion'
      NAMES(7) = 'HICO'
      NAMES(8) = 'NIS'
      NAMES(9) = 'PRISM'
      NAMES(10)= 'OTHERS?'
!C***End of temp coding ------------

      print*, 'Instrument Name: ', NAME_INSTRU

!C Get the plane altitude (in km, above sea level). Must > HSURF.
      READ(*,*) XPPP

!C--      print*, 'XPPP = ', XPPP

!C  Get date and time of measurements.  Format: MM DD YYYY HH MM SS
      READ (5,*) IMN,IDY,IYR,IH,IM,IS
      IF((IMN.LT.1).OR.(IMN.GT.12).OR.(IDY.LT.1).OR.(IDY.GT.31).OR. &
        (IYR.LT.1987).OR.(IYR.GT.2020)) THEN
        WRITE(*,*)'Invalid date:',IMN,IDY,IYR
        WRITE(*,*)'Format is MM DD YYYY where 0<MM<12, 0<DD<32,  &
      YYYY>1986.'
        GOOD_DATA = .FALSE. 
      ENDIF
      IF((IH.LT.0).OR.(IH.GT.24).OR.(IM.LT.0).OR.(IM.GT.60).OR. &
        (IS.LT.0).OR.(IS.GT.60)) THEN
        WRITE(*,*)'Invalid time:',IH,IM,IS
        WRITE(*,*)'Format is HH MM SS where 0<HH<24, 0<MM<60, SS<60.'
        GOOD_DATA = .FALSE. 
      ENDIF

!C  Get latitude of measurements.  Format: degrees minutes seconds
      READ (5,*) XLATD,XLATM,XLATS
      IF((XLATD.GT.90).OR.(XLATM.GT.60).OR.(XLATS.GT.60)) THEN
        WRITE(*,*)'Invalid latitude:',XLATD,XLATM,XLATS
        WRITE(*,*)'Format: degrees minutes seconds'
        WRITE(*,*)'Valid values are degrees < 90, minutes < 60, &
      seconds<60.'
        GOOD_DATA = .FALSE. 
      ENDIF

!C Get hemisphere of latitude.  Format: "N" or "S"
      READ (5,127) LATHEM
  127 FORMAT(A1)
      IF((LATHEM.NE.'N').AND.(LATHEM.NE.'S')) THEN
        WRITE(*,*)'Invalid hemisphere value:',LATHEM
        WRITE(*,*)'Valid values are "N" or "S".'
        GOOD_DATA = .FALSE. 
      ENDIF

!C  Get longitude of measurements.  Format: degrees minutes seconds
      READ (5,*) XLONGD,XLONGM,XLONGS
      IF((XLONGD.GT.180).OR.(XLONGM.GT.60).OR.(XLONGS.GT.60)) THEN
        WRITE(*,*)'Invalid longitude:',XLONGD,XLONGM,XLONGS
        WRITE(*,*)'Format: degrees minutes seconds'
        WRITE(*,*)'Valid values are degrees < 180 minutes < 60, &
      seconds<60.'
        GOOD_DATA = .FALSE. 
      ENDIF

!C Get hemisphere of longitude.  Format: "E" or "W"
      READ (5,127) LNGHEM
      IF((LNGHEM.NE.'E').AND.(LNGHEM.NE.'W')) THEN
        WRITE(*,*)'Invalid hemisphere:',LNGHEM
        WRITE(*,*)'Valid values are "E" or "W".'
        GOOD_DATA = .FALSE. 
      ENDIF

!C  Get view zenith angle.  Format: degrees minutes seconds
      READ (5,*) XVIEWD,XVIEWM,XVIEWS
      IF((XVIEWD.GT.90).OR.(XVIEWM.GT.60).OR.(XVIEWS.GT.60)) THEN
        WRITE(*,*)'Invalid latitude:',XVIEWD,XVIEWM,XVIEWS
        WRITE(*,*)'Format: degrees minutes seconds'
        WRITE(*,*)'Valid values are degrees < 90, minutes < 60 & seconds<60.'
        GOOD_DATA = .FALSE. 
      ENDIF

!C  Get relative azimuth angle between the solar ray and view direction with
!C      both projected on the ground.  Format: degrees minutes seconds
      READ (5,*) XAZMUD,XAZMUM,XAZMUS
      IF((XAZMUD.GT.360.).OR.(XAZMUM.GT.60).OR.(XAZMUS.GT.60)) THEN
        WRITE(*,*)'Invalid view azimuth ang.:',XAZMUD,XAZMUM,XAZMUS
        WRITE(*,*)'Format: degrees minutes seconds'
        WRITE(*,*)'Valid values are degrees < 360, minutes < 60, & seconds<60.'
        GOOD_DATA = .FALSE. 
      ENDIF

!C Get spectral resolutions (FWHM) in micron. If the value is 0, that means
!C read the FWHM data from the wavelength file.  Otherwise, a non-zero value
!C will be used as a the constant FWHM value.
      READ(*,*)DLT
      IF (DLT.GT.0.0) THEN
        IF((DLT.LT.0.0005).OR.(DLT.GT.0.1)) THEN
          WRITE(*,*)'Invalid spectral resolution value:',DLT
          WRITE(*,*)'Valid values are 0.0005-0.1 micron.'
          GOOD_DATA = .FALSE. 
        ELSE
!C Initialize the FWHM array
!C***          DLTNEW=DLT/1000.
          DO 521 I=1,1024
            FWHM(I)=DLT
  521     CONTINUE
        ENDIF
      ENDIF

!C  Get imaging spectrometer wavelength file name
      READ(*,85) FINPWV
  85  FORMAT(A1000)
      OPEN(11,FILE=FINPWV,STATUS='OLD',IOSTAT=ISTAT)
      IF(ISTAT.NE.0) THEN
        WRITE(*,*)'wavelength file did not open successfully:',FINPWV
        GOOD_DATA = .FALSE. 
      ELSE

!C  Read wavelength data.  For 1991 and before, the file contains the channel
!C  number and the wavelength for each channel.  For 1992, the file also contains
!C  a third column with the resolution (FWHM) in micron for each channel.

        NOBS = 1024           !initially set to max allowable # of observations
        DO 510 I=1,NOBS
          IF (DLT.GT.0.) THEN
            READ(11,*,END=520) X,WAVOBS(I)
          ELSE
            READ(11,*,END=520) X,WAVOBS(I), FWHM(I)
          ENDIF
 510    CONTINUE
 520    NOBS=I-1              !set to actual # of observations in input file
        CLOSE(11)
      ENDIF

!C Determine if channel ratio parameters should be read.  Format: 0=no 1=yes
      READ(*,*) ANS
      IF((ANS.NE.0).AND.(ANS.NE.1)) THEN
	WRITE(*,*)'Invalid value to indicate whether the channel ratio parameters should be read:',ANS
        WRITE(*,*)'Valid values: 0=no 1=yes.'
        GOOD_DATA = .FALSE. 
      ELSE
        IF (ANS .EQ. 1) THEN
!C  Get ranges on curve for the two atmospheric windows surrounding the .94-um
!C    water vapor absorption feature and the center point of the .94-um water
!C    vapor absorption feature.  Enter:
!C         1. the midpoint of first window (0.6-2.5)
!C         2. number of points to average for first window (1-10)
!C         3. the midpoint of second window (0.6-2.5)
!C         4. number of points to average for second window (1-10)
!C         5. the midpoint of water vapor absorption feature (0.6-2.5)
!C         6. the number of points to average absorption feature (1-30)

  114     READ(*,*)WNDOW1, NB1, WNDOW2, NB2, WP94C, NBP94

          IF((WNDOW1.LT.0.6).OR.(WNDOW1.GT.2.5)) THEN
  	    WRITE(*,*)'Invalid wavelength position for first atmospheric window in the .94-um region:',WNDOW1
            WRITE(*,*)'Valid values are 0.6-2.5.'
            GOOD_DATA = .FALSE. 
          ENDIF

          IF((NB1.LT.1).OR.(NB1.GT.50)) THEN
   	    WRITE(*,*)'Invalid number of channels for first wavelength position in the .94-um region:',NB1
            WRITE(*,*)'Valid values are 1-50.'
            GOOD_DATA = .FALSE. 
          ENDIF

          IF((WNDOW2.LT.0.6).OR.(WNDOW2.GT.2.5)) THEN
  	    WRITE(*,*)'Invalid wavelength position for second atmospheric window in the .94-um region:',WNDOW2
            WRITE(*,*)'Valid values are 0.6-2.5.'
            GOOD_DATA = .FALSE. 
          ENDIF

          IF((NB2.LT.1).OR.(NB2.GT.50)) THEN
 	    WRITE(*,*)'Invalid number of channels for second wavelength position in the .94-um region:',NB2
            WRITE(*,*)'Valid values are 1-50.'
            GOOD_DATA = .FALSE. 
          ENDIF

          IF((WP94C.LE.WNDOW1).OR.(WP94C.GE.WNDOW2)) THEN
 	    WRITE(*,*)'Invalid water vapor wavelength position for the .94-um region:',WP94C
            WRITE(*,*)'Valid range is:',WNDOW1,' < value < ',WNDOW2,'.'
            GOOD_DATA = .FALSE. 
          ENDIF

          IF((NBP94.LT.1).OR.(NBP94.GT.90)) THEN
 	    WRITE(*,*)'Invalid number of channels for water vapor wavelength position in the .94-um region:',NBP94
            WRITE(*,*)'Valid values are 1-90.'
            GOOD_DATA = .FALSE. 
          ENDIF

!C  Get ranges on curve for the two atmospheric windows surrounding the 1.14-um
!C    water vapor absorption feature and the center point of the 1.14-um water
!C    vapor absorption feature.  Enter:
!C         1. the midpoint of third window (0.6-2.5)
!C         2. number of points to average for third window (1-10)
!C         3. the midpoint of fourth window (0.6-2.5)
!C         4. number of points to average for fourth window (1-10)
!C         5. the midpoint of 1.14-um absorption feature (0.6-2.5)
!C         6. the number of points to average for the absorption feature (1-30)


          READ(*,*)WNDOW3, NB3, WNDOW4, NB4, W1P14C, NB1P14

          IF((WNDOW3.LT.0.6).OR.(WNDOW3.GT.2.5)) THEN
  	    WRITE(*,*)'Invalid wavelength position for first atmospheric window in the 1.14-um region:',WNDOW3
            WRITE(*,*)'Valid values are 0.6-2.5'
            GOOD_DATA = .FALSE. 
          ENDIF

          IF((NB3.LT.1).OR.(NB3.GT.50)) THEN
 	    WRITE(*,*)'Invalid number of channels for first wavelength position in the 1.14-um region:',NB3
            WRITE(*,*)'Valid values are 1-50.'
            GOOD_DATA = .FALSE. 
          ENDIF

          IF((WNDOW4.LT.0.6).OR.(WNDOW4.GT.2.5)) THEN
  	    WRITE(*,*)'Invalid wavelength position for second atmospheric window in the 1.14-um region:',WNDOW4
            WRITE(*,*)'Valid values are 0.6-2.5'
            GOOD_DATA = .FALSE. 
          ENDIF

          IF((NB4.LT.1).OR.(NB4.GT.50)) THEN
 	    WRITE(*,*)'Invalid number of channels for second wavelength position in the 1.14-um region:',NB4
            WRITE(*,*)'Valid values are 1-50.'
            GOOD_DATA = .FALSE. 
          ENDIF

          IF((W1P14C.LE.WNDOW3).OR.(W1P14C.GE.WNDOW4)) THEN
	    WRITE(*,*)'Invalid water vapor wavelength position for 1.14-um:',W1P14C
            WRITE(*,*)'Valid range is:',WNDOW3,' < value <', WNDOW4,'.'
            GOOD_DATA = .FALSE. 
          ENDIF

          IF((NB1P14.LT.1).OR.(NB1P14.GT.110)) THEN
	    WRITE(*,*)'Invalid number of channels for water vapor wavelength position in the 1.14-um region:',NB1P14
            WRITE(*,*)'Valid values are 1-110.'
            GOOD_DATA = .FALSE. 
          ENDIF
        ELSE

!C Initialize default atmospheric window and water vapor absorption regions
!C     for 3-channel ratioing calculations.
          WNDOW1 = 0.865
          NB1    = 3
          WNDOW2 = 1.030
          NB2    = 3
          WP94C  = 0.940
          NBP94  = 5
          WNDOW3 = 1.050
          NB3    = 3
          WNDOW4 = 1.235
          NB4    = 3
          W1P14C = 1.1375
          NB1P14 = 7
        ENDIF
      ENDIF  

!C  Determine which atmospheric model to use
!C     1 = tropical
!C     2 = mid latitude summer
!C     3 = mid latitude winter
!C     4 = subarctic summer
!C     5 = subarctic winter
!C     6 = US standard 1962
!C     7 = user defined model
      READ(*,*)MODEL 
      IF((MODEL.LT.1).OR.(MODEL.GT.7)) THEN
	WRITE(*,*)'Invalid atmospheric model value:',MODEL
        WRITE(*,*)'Valid values are 1-7.'
        GOOD_DATA = .FALSE. 
        MODEL = 1		! init model to a valid value for use below
      ENDIF
      
      IF (MODEL.EQ.7) THEN
!C Get file name of atmospheric model information
        READ(*,101) FTPVMR
  101   FORMAT(A1000)
        OPEN(20,FILE=FTPVMR,STATUS='OLD',IOSTAT=ISTAT)
        IF(ISTAT.NE.0) THEN
          WRITE(*,*)'Atmospheric model file did not open  successfully:',FTPVMR
          GOOD_DATA = .FALSE. 
        ELSE
!C Read in number of boundaries, temperature (K), pressure (mb), and water
!C vapor volume mixing ratio (VMR, in units of parts per million (ppm)) profiles
          READ(20,*) NB
          IF((NB.LE.0).OR.(NB.GT.25)) THEN
	    WRITE(*,*)'Invalid number of boundaries:',NB
            WRITE(*,*)'Valid values are 1-25.'
            GOOD_DATA = .FALSE. 
          ELSE
            DO 300 I=1,NB
              READ(20,*) H(I),P(I),T(I),VMR(I)
  300       CONTINUE
            CLOSE(20,IOSTAT=ISTAT)
          ENDIF
        ENDIF
      ELSE
!C       Initialize NB, H, P, T, VMR from predefined atmospheric model array
        NB = TPVMR(MODEL,1)
        DO 120 I=1,NB
          H(I) = TPVMR(MODEL,2+(4*(I-1)))
          P(I) = TPVMR(MODEL,3+(4*(I-1)))
          T(I) = TPVMR(MODEL,4+(4*(I-1)))
          VMR(I) = TPVMR(MODEL,5+(4*(I-1)))
  120   CONTINUE
      ENDIF
      NL=NB-1

        DO I     = NB+1, 25
          H(I)   = 1000.
          P(I)   = 0.0
          T(I)   = 300.
          VMR(I) = 0.0
        END DO
!
!C Determine if various gases should be included in the calculations.  Format:
!C 1=yes or 0=no.  The order should be:
!C                     1. water vapor
!C                     2. carbon dioxide
!C                     3. ozone
!C                     4. nitrous oxide
!C                     5. carbon monoxide
!C                     6. methane
!C                     7. oxygen
!C                     8. nitrogen dioxide

      READ(*,*)IH2OVP, ICO2, IO3, IN2O, ICO, ICH4, IO2, INO2
      IF((IH2OVP.NE.0).AND.(IH2OVP.NE.1)) THEN
	WRITE(*,*)'Invalid selection for H2O vapor:',IH2OVP
        WRITE(*,*)'Valid values: 0=no, 1=yes.'
        GOOD_DATA = .FALSE. 
      ENDIF

      IF((ICO2.NE.0).AND.(ICO2.NE.1)) THEN
	WRITE(*,*)'Invalid selection for CO2:',ICO2
        WRITE(*,*)'Valid values: 0=no, 1=yes.'
        GOOD_DATA = .FALSE. 
      ENDIF

      IF((IO3.NE.0).AND.(IO3.NE.1)) THEN
	WRITE(*,*)'Invalid selection for O3:',IO3
        WRITE(*,*)'Valid values: 0=no, 1=yes.'
        GOOD_DATA = .FALSE. 
      ENDIF

      IF((IN2O.NE.0).AND.(IN2O.NE.1)) THEN
	WRITE(*,*)'Invalid selection for N2O:',IN2O
        WRITE(*,*)'Valid values: 0=no, 1=yes.'
        GOOD_DATA = .FALSE. 
      ENDIF

      IF((ICO.NE.0).AND.(ICO.NE.1)) THEN
	WRITE(*,*)'Invalid selection for CO:',ICO
        WRITE(*,*)'Valid values: 0=no, 1=yes.'
        GOOD_DATA = .FALSE. 
      ENDIF

      IF((ICH4.NE.0).AND.(ICH4.NE.1)) THEN
	WRITE(*,*)'Invalid selection for CH4:',ICH4
        WRITE(*,*)'Valid values: 0=no, 1=yes.'
        GOOD_DATA = .FALSE. 
      ENDIF

      IF((IO2.NE.0).AND.(IO2.NE.1)) THEN
	WRITE(*,*)'Invalid selection for O2:',IO2
        WRITE(*,*)'Valid values: 0=no, 1=yes.'
        GOOD_DATA = .FALSE. 
      ENDIF

      IF((INO2.NE.0).AND.(INO2.NE.1)) THEN
	WRITE(*,*)'Invalid selection for NO2:',INO2
        WRITE(*,*)'Valid values: 0=no, 1=yes.'
        GOOD_DATA = .FALSE. 
      ENDIF

!C Read in the column ozone amount and scaling factor for NO2 column amount.
!C    The typical ozone amount is: 0.28-0.55 (atm-cm). The built-in NO2
!C    column amount is 5.0E+15 molecules/cm^2.
      READ(*,*)VRTO3, SNO2
      IF((VRTO3.LT.0.1).OR.(VRTO3.GT.0.6)) THEN
	WRITE(*,*)'Invalid vertical column amount of ozone:',VRTO3
        WRITE(*,*)'Valid values are 0.1-0.6.'
        WRITE(*,*)'Reasonable values are 0.28-0.55.'
        GOOD_DATA = .FALSE. 
      ENDIF

!C Get aerosol model, visibility, or aerosol optical depth at 0.55 um (TAER55)
!C for SSSSS(). If V > 0, no TAER55 is needed in the input file. If V = 0,
!C TAER55 is needed.
!C IAER: 0=no aerosol, set V = -1
!C IAER: 1=continental aerosol, 2=maritime aerosol, 3=urban aerosol, set V>5 km.

      READ(*,*) IAER,V
      IF((IAER.LT.0).OR.(IAER.GT.6)) THEN
	WRITE(*,*)'Invalid aerosol model value:',IAER
        WRITE(*,*)'Valid values are 0-3.'
        GOOD_DATA = .FALSE. 
      ENDIF

      IF(IAER.EQ.0) THEN
        V = -1
        TAER55=0.0
       ELSE
        IF((V.LT.0).OR.(V.GT.300)) THEN
	  WRITE(*,*)'Invalid visibility value:',V
          WRITE(*,*)'Value must be greater than 0 and less than 300.'
          GOOD_DATA = .FALSE. 
       ENDIF
      ENDIF

      IF(V.EQ.0) THEN
        READ(*,*) TAER55
        IF((TAER55.GT.10.).OR.(TAER55.LE.0)) THEN
	  WRITE(*,*)'Invalid aerosol optical depth at 550 nm: ',TAER55
          WRITE(*,*)'Valid values are between 0 and 10.'
          GOOD_DATA = .FALSE. 
        ENDIF
      ENDIF

!C Get the mean surface elevation.  Range: 0.0 - H(NB-1) km
      READ(*,*)HSURF 
!C--
      XPSS = HSURF
!C--
      IF((HSURF.LT.0).OR.(HSURF.GT.H(NB-1))) THEN
	WRITE(*,*)'Invalid surface elevation value:',HSURF
        WRITE(*,*)'Value must be less than the maximum elevation in the atmospheric model.'
        GOOD_DATA = .FALSE. 
      ENDIF

!C Test to see if plane altitude is greater than bottom surface elevation
      IF(XPPP .LT. HSURF) THEN
        WRITE(*,*)'Invalid plane altitude:',XPPP
        WRITE(*,*)'Plane altitude must be greater than the bottom surface elevation.'
        GOOD_DATA = .FALSE.
      ENDIF

!C  Read in file name of data cube and get data dimensions in DIMS, ST_ORDER
!C  If no SIPS header on data, then default dimensions are header size = 0,
!C  number of samples = 614, number of lines = 512, number of bands = 224,
!C  storage order = BIL
      READ(*,85) FINAV
!C------      CALL OPENINFILE(IRET,DIMS,ST_ORDER)
      LUN_IN  = 47
!      CALL OPENINFILE(LUN_IN, I_RET)

!      IF(I_RET.NE.0) THEN
!        WRITE(*,*)'Input image file did not open
!     & successfully:',FINAV
!        GOOD_DATA = .FALSE.
!      ENDIF

!C Determine if cube dimensions should be read in.  If not, the dimensions
!C in the image file header will be used.
      READ(*,*) ANS
      IF((ANS.NE.0).AND.(ANS.NE.1)) THEN
	WRITE(*,*)'Invalid value to indicate whether the cube dimensions should be read:',ANS
        WRITE(*,*)'Valid values: 0=no 1=yes.'
        GOOD_DATA = .FALSE. 
      ELSE
        IF (ANS .EQ. 1) THEN
!C         Read the header size in bytes, samples, lines, bands, and
!C           storage order (0=BSQ, 1=BIP, 2=BIL)'
          READ(*,*) HDREC, NSAMPS, NLINES, NBANDS, SORDER
          IF ((HDREC.LT.0).OR.(NSAMPS.LE.0).OR. (NLINES.LE.0).OR.(NBANDS.LE.0).OR. &
             (SORDER.LE.0).OR.(SORDER.GT.2)) THEN
	    WRITE(*,*)'Invalid cube parameters:',HDREC,NSAMPS,NLINES, &
      NBANDS,SORDER
	    WRITE(*,*)'Values must be greater than or equal to 1 and the storage order must be less than 3.'
      WRITE(*,*)'Storage Order (0=BSQ, 1=BIP, 2=BIL)'
      WRITE(*,*)' 0 = BSQ is not supported in this version of ATREM'
            GOOD_DATA = .FALSE. 
          ENDIF
        ELSE
	  HDREC  = DIMS(1)
          NSAMPS = DIMS(2)
          NLINES = DIMS(3)
          NBANDS = DIMS(4)
	  SORDER = ST_ORDER
        ENDIF  
      ENDIF  

!C Get the output reflectance cube file name.  The output cube always has a 512
!C byte header record.
      READ(*,85) FOCUB
!C---      CALL OPENOUTFILE(IRET,512,NSAMPS,NLINES,NBANDS,SORDER)
      LUN_OUT = 48
!      CALL OPENOUTFILE(LUN_OUT, I_RET)

!      IF(I_RET.NE.0) THEN
!        WRITE(*,*)'Output image file did not open
!     & successfully:',FOCUB
!        GOOD_DATA = .FALSE.
!      ENDIF

!C Get the resolution of output spectra.  Format: 0-100 nm.
      READ(*,*)DLT2
      IF((DLT2.LT.0.).OR.(DLT2.GT.100.)) THEN
	WRITE(*,*)'Invalid output resolution value:',DLT2
        WRITE(*,*)'Valid values are 0-100 nm.'
        GOOD_DATA = .FALSE. 
      ENDIF

!C Get the scaling factor for output reflectance values.  Range: 1. to 32000.
      READ(*,*)SCALEF
      IF((SCALEF.LT.1.).OR.(SCALEF.GT.32000.)) THEN
	WRITE(*,*)'Invalid output resolution value:',SCALEF
        WRITE(*,*)'Valid values are 1-32000 nm.'
        GOOD_DATA = .FALSE. 
      ENDIF

!C Get the output water vapor file name.  The water vapor output has a 512-byte
!C header, is always one band, and is in BSQ format (BSQ=0).
      READ(*,85) FOH2O
!C---      CALL OPENVAPFILE(IRET,512,NSAMPS,NLINES,1,0)
      LUN_VAP = 49
!      CALL OPENVAPFILE(LUN_VAP, I_RET)

!      IF(I_RET.NE.0) THEN
!        WRITE(*,*)'Output water vapor file did not open
!     & successfully:',FOH2O
!        GOOD_DATA = .FALSE.
!      ENDIF

!C  Get the file name for output table consisting of atmospheric
!C  transmission spectra at the solar and observational geometry but
!C  with varying amounts of water vapor.

!COMMENTED out RJH
!      READ(*,85) FOUT1
!      OPEN(35,FILE=FOUT1,STATUS='UNKNOWN',IOSTAT=ISTAT)
!      IF(ISTAT.NE.0) THEN
!        WRITE(*,*)'Output file did not open successfully:',FOUT1
!        GOOD_DATA = .FALSE.
!      ENDIF
!
!      IF(GOOD_DATA .NEQV. .TRUE.) THEN
!        WRITE(*,*)'**** ERROR: Invalid input detected. ****'
!        STOP
!      ENDIF

!C--D     WRITE(*,*)'IH2OVP=',IH2OVP,' ICO2=',ICO2,' IO3=',IO3,' IN2O=',IN2O
!C--D     WRITE(*,*)'ICO=',ICO,' ICH4=',ICH4,' IO2=',IO2
!C--D     WRITE(*,*)'H         T        P        VMR'
!C--D     DO 998 I=1,25
!C--D 998   WRITE(*,*)H(I),' ',T(I),' ',P(I),' ',VMR(I)
!C--D     WRITE(*,*)'NB=',NB,' NL=',NL,' MODEL=',MODEL
!C--D     WRITE(*,*)' IAER=',IAER,' V=',V,' VRTO3=',VRTO3
!C--D     DO 997 I=1,NOBS
!C--D 997   WRITE(*,*)'I=',I,' WAVOBS(I)=',WAVOBS(I),' FWHM(I)=',FWHM(I)
!C--D     WRITE(*,*)'NOBS=',NOBS,' HSURF=',HSURF,' DLT=',DLT,' DLT2=',DLT2
!C--D     WRITE(*,*)'SCALEF=',SCALEF
!C--D     WRITE(*,*)'WNDOW1=',WNDOW1,' WNDOW2=',WNDOW2,' WP94C=',WP94C
!C--D     WRITE(*,*)'WNDOW3=',WNDOW3,' WNDOW4=',WNDOW4,' W1P14C=',W1P14C
!C--D     WRITE(*,*)'NB1=',NB1,' NB2=',NB2,' NBP94=',NBP94
!C--D     WRITE(*,*)'NB3=',NB3,' NB4=',NB4,' NB1P14=',NB1P14
!C--D     WRITE(*,*)'IMN=',IMN,' IDY=',IDY,' IYR=',IYR
!C--D     WRITE(*,*)'IH=',IH,' IM=',IM,' IS=',IS
!C--D     WRITE(*,*)'XLATD=',XLATD,' XLATM=',XLATM,' XLATS=',XLATS
!C--D     WRITE(*,147)LATHEM
!C--D 147 FORMAT('LATHEM=',A1)
!C--D     WRITE(*,*)'XLONGD=',XLONGD,' XLONGM=',XLONGM,' XLONGS=',XLONGS
!C--D     WRITE(*,148)LNGHEM
!C--D 148 FORMAT('LNGHEM=',A1)
!C--D     WRITE(*,*)'HDREC=',HDREC,' NSAMPS=',NSAMPS,' NLINES=',NLINES
!C--D     WRITE(*,*)'NBANDS=',NBANDS,' SORDER=',SORDER

      RETURN
      END

!********************************************************************************
!*            								       *
!*  Name: MODEL_ADJ							       *
!*  Purpose: resets the bottom boundary of the input model if the surface       *
!*           elevation is greater than 0, and calculate the column water vapor  *
!*           amount in the selected model.				       *
!*  Parameters: none							       *
!*  Algorithm: If the surface elevation > 0, the bottom layer temperature and   *
!*             water vapor volume mixing ratio are obtained through linear      *
!*             interpolation, while the bottom layer pressure is obtained       *
!*             through exponential interpolation.			       *
!*  Globals used:  H, T, P, VMR, NB, NL, HSURF - these values are adjusted if   *
!*                      HSURF > 0					       *
!*  Global output:  CLMVAP - Column water vapor amount in unit of cm.           *
!*                       Q - Number of molecules above the surface at one       *
!*                           atmosphere in units of molecules/cm**2	       *
!*  Return Codes: none							       *
!*  Special Considerations: none						       *
!*									       *
!********************************************************************************
      SUBROUTINE MODEL_ADJ

      INCLUDE 'COMMONS_INC.f'

!C  Common variables
      DIMENSION H(25), T(25), P(25), VMR(25)

      COMMON /GETINPUT3/ H,T,P,VMR,NB,NL,MODEL,IAER,V,TAER55,VRTO3,SNO2
      COMMON /GETINPUT5/ NOBS,HSURF,DLT,DLT2
      COMMON /MODEL_ADJ1/ CLMVAP,Q

      DIMENSION HP(25), TP(25), PP(25), VMRP(25)
      COMMON /MODEL_ADJ2/ HP, TP, PP, VMRP
      COMMON /MODEL_ADJ3/ K_PLANE, DVAP_PLANE, DVAP_LAYER, &
                         DP_PLANE, DP_LAYER, CLMVAPP
      COMMON /MODEL_ADJ4/ K_SURF

      REAL XPSS, XPPP
      COMMON /GETINPUT14/ XPSS, XPPP

!C     H = Layer boundary height, SL=Slant path length, DSLODH=DSL/DH

      RE=6380.0
!C     Q=# of molecules above the surface at one atmosphere (molecules/cm**2)
      Q=2.152E25

!C--        print*, 'original model atmosphere ...'
!C--      DO I = 1, 25
!C--        print*, 'H,T,P,VMR = ', H(I), T(I), P(I), VMR(I)
!C--      END DO

!C--- Convert the atmospheric pressure from "mb" to "atm.":
      DO I = 1, NB
         P(I) = P(I) / 1013.
      END DO
!C--- End of conversion

!C     Convert the VMR from the ppm unit in the model to absolute unit.
      DO 310 I=1,NB
  310   VMR(I)=VMR(I)*1.0E-06

!C     Do special processing if surface altitude is greater than 0
!C---      IF(HSURF.NE.0.0) THEN
!C       Reset H(I) to a smaller value if HSURF.EQ.H(I) to avoid possible
!C         problems in table searching using LOCATE
        DO 7455 I=1,NB
 7455     IF(HSURF.EQ.H(I)) HSURF=H(I)+0.0001

!C       Determine index of H() such that H(index) < HSURF < H(index+1)
        CALL LOCATE(H,NB,HSURF,K)
!C
!C K_SURF is an index relative to the original model atmosphere (0 - 100 km)
        K_SURF = K
!C---        print*, 'K_SURF = ', K_SURF
!C
        IF(K.EQ.0) THEN
          WRITE(6,5237) 
 5237     FORMAT(2X,'***WARNING: Surface elevation smaller then lowest boundary of the model atmosphere.')

          K_SURF = 1
          GOTO 5255
        ENDIF

        IF(K.GT.0) THEN
          DHK=H(K+1)-H(K)
          DHS=HSURF -H(K)

!C         linear interpolation for surface temperature (TSURF) and VMR (  VMRS)
          TSURF=T(K)+(DHS/DHK)*(T(K+1)-T(K))
          VMRS =VMR(K)+(DHS/DHK)*(VMR(K+1)-VMR(K))

!C         exponential interpolation for surface pressure (PSURF)
          PSURF=P(K)*EXP(-ALOG(P(K)/P(K+1))*DHS/DHK)
          H(1)=HSURF
          P(1)=PSURF
          T(1)=TSURF
          VMR(1)=VMRS

          NB=NB-K+1
          NL=NB-1
!C---          print*,'NL = ', NL, ' NB = ', NB, 'in Model_Adj'
!C
          DO 5240 I=2,NB
            H(I)=H(I+K-1)
            P(I)=P(I+K-1)
            T(I)=T(I+K-1)
            VMR(I)=VMR(I+K-1)
 5240     CONTINUE
!C  Zero out pressures and VMRS of top atmospheric layers.
          DO 5245 I=NB+1,25
            H(I)=1000.
            P(I)=0.0
            T(I)=300.
            VMR(I)=0.0
 5245     CONTINUE
        ENDIF
!C---      ENDIF

 5255 CONTINUE

      AMTVRT=0.0

      DO 350 I=1,NL
        DAMTVT=Q*(P(I)-P(I+1))*(VMR(I)+VMR(I+1))/2.0
        AMTVRT=AMTVRT+DAMTVT
 350  CONTINUE

      CLMVAP=AMTVRT/3.34E+22

      WRITE(*,*)'Column vapor amount in model atmosphere from ground'
      WRITE(*,*)'       to space = ', CLMVAP, ' cm'

!C---      WRITE(*,*)'In MODEL_ADJ, NB = ',NB,' NL = ',NL
!C---        print*, 'After adjusting for elevated surface ...'
!C---      DO I = 1, 25
!C---        print*, 'H,T,P,VMR = ', H(I), T(I), P(I), VMR(I)
!C---      END DO
!
!C
!C
!C Setting the upward atmospheric path's T, P, and VMR profiles:
!C
!C  1st duplicate the entire atmospheric profiles from the downward path
!C      to the upward path
!C
      DO I = 1, 25
         HP(I)    = H(I)
         TP(I)    = T(I)
         PP(I)    = P(I)
         VMRP(I)  = VMR(I)
      END DO


      HPLANE = XPPP
!C   Set the highest plane altitude to the upper bound of model atmosphere
!C---      IF(HPLANE.GT.100.0) HPLANE = 100. - 0.0001
      IF(HPLANE.GE.100.0) HPLANE = 100. - 0.0001
!C
!C  Do special processing if the plane height (HPLANE) is greater than HP(1)
      IF(HPLANE.GT.HP(1)) THEN
!C         Reset Plane altitude HPLANE (= XPPP) to a larger value if
!C             HPLANE.EQ.HP(I) to avoid possible problems in table
!C             searching using LOCATE
        DO 7456 I=1,25
 7456     IF(HPLANE.EQ.HP(I)) HPLANE=HP(I)-0.0001

!C       Determine index of HP() such that HP(index) < HPLANE < H(index+1)
        CALL LOCATE(HP,NB,HPLANE,KK)

        IF(KK.EQ.0) THEN
          WRITE(6,5239) 
 5239     FORMAT(2X,'***WARNING: Plane altitude less then lowest boundary of the model atmosphere.')
          GOTO 5256
        ENDIF

        IF(KK.GT.0) THEN
          DHKK = HP(KK+1) - HP(KK)
          DHSS = HPLANE   - HP(KK)

!C         linear interpolation for plane temperature (TPLANE) and VMR (  VMRSP)
          TPLANE = TP(KK)   + (DHSS/DHKK)*(TP(KK+1)-TP(KK))
          VMRSP  = VMRP(KK) + (DHSS/DHKK)*(VMRP(KK+1)-VMRP(KK))

!C         exponential interpolation for plane pressure (PPLANE)
          PPLANE     = PP(KK)*EXP(-ALOG(PP(KK)/PP(KK+1))*DHSS/DHKK)
          HP(KK+1)   = HPLANE
          PP(KK+1)   = PPLANE
          TP(KK+1)   = TPLANE
          VMRP(KK+1) = VMRSP

!C  Zero out pressures and VMRP of top atmospheric layers.
          IF(KK.LT.24) THEN
           DO I=KK+2,25
              HP(I)=1000.
              PP(I)=0.0
              TP(I)=300.
              VMRP(I)=0.0
           END DO 
          END IF

        ENDIF
      ENDIF

 5256 CONTINUE

      AMTVRTP=0.0

      DO 357 I=1,KK
        DAMTVTP=Q*(PP(I)-PP(I+1))*(VMRP(I)+VMRP(I+1))/2.0
        AMTVRTP=AMTVRTP+DAMTVTP
 357  CONTINUE

      CLMVAPP=AMTVRTP/3.34E+22

      WRITE(*,*)'Column vapor below plane (CLMVAPP) = ', &
      CLMVAPP, ' cm'

!C---      WRITE(*,*)'In MODEL_ADJ, NB = ',NB,' KK = ', KK
!C---        print*, 'After further adjusting for plane height ...'
!C---      DO I = 1, 25
!C---        print*, 'HP,TP,PP,VMRP = ', HP(I), TP(I), PP(I), VMRP(I)
!C---      END DO

!C--- Indices and parameters for the plane layer
      K_PLANE = KK

      DVAP_PLANE = Q*(PP(K_PLANE) - PP(K_PLANE+1))* &
       (VMRP(K_PLANE) + VMRP(K_PLANE+1))/2.0 / 3.34E+22

      DVAP_LAYER = Q*(P(K_PLANE) - P(K_PLANE+1))*  &
       (VMR(K_PLANE) + VMR(K_PLANE+1))/2.0 / 3.34E+22

      DP_PLANE = PP(K_PLANE) - PP(K_PLANE+1)
      DP_LAYER = P(K_PLANE)  - P(K_PLANE+1)

!C---      print*, 'K_PLANE, DVAP_PLANE, DVAP_LAYER = ',
!C---     &         K_PLANE, DVAP_PLANE, DVAP_LAYER
!C---      print*, 'DP_PLANE, DP_LAYER = ', DP_PLANE, DP_LAYER


      RETURN
      END
!********************************************************************************
!*            								       *
!*  Name: GEOMETRY							       *
!*  Purpose: Calculates the solar and the observational geometric factors.      *
!*  Parameters: none							       *
!*  Algorithm: The solar geometry was obtained based on the latitude, longitude,*
!*             GMT time using programs written by W. Mankin at National Center  *
!*             for Atmospheric Research in Boulder, Colorado. The	       *
!*             geometric factors for CO2, O3, N2O, CO, CH4, and O2 are based    *
!*             only on the solar and observational angles. Sixty artificial     *
!*             geometric factors for H2O are set up to produce a transmittance  *
!*             table for different atmospheric water vapor amounts.	       *
!*  Globals used:  VRTO3	- Column O3 amount in units of atm-cm		       *
!*      IMN,IDY,IYR,IH,IM,IS - time and date of data measurements	       *
!*      XLATD,XLATM,XLATS,LATHEM	- Latitude of measured area		       *
!*      XLONGD,XLONGM,XLONGS,LNGHEM - Longitude of measured area		       *
!*      CLMVAP - Column water vapor in unit of cm in the model atmosphere       *
!*  Global output:  							       *
!*      SOLZNI,SOLAZ,OBSZNI,OBSPHI,IDAY - Solar zenith angle, solar azimuth     *
!*            angle, observational zenith angle, observational azimuth angle,   *
!*            and the day number in the year 				       *
!*      GCO2,GO3,GN2O,GCO,GCH4,GO2,SSH2O,GGEOM,TOTLO3 - Geometric factors for   *
!*            different gases, and total O3 amount in the sun-surface ray path. *
!*            The geometric factor is defined as: if the vertical column amount *
!*            of the gas is equal 1, then GGAS is equal to the total amount of  *
!*            the gas in the combined Sun-surface-sensor ray path.	       *
!*  Return Codes: none							       *
!*  Special Considerations: none						       *
!*									       *
!********************************************************************************

      SUBROUTINE GEOMETRY

      INCLUDE 'COMMONS_INC.f'

      DIMENSION VAPVRT(60), VAP_SLANT(60)
      DIMENSION MD(12)
      DIMENSION SSH2O(60)
      DIMENSION H(25), T(25), P(25), VMR(25)
      CHARACTER*1 LATHEM,LNGHEM

!C The VAPVRT array contains 60 column water vapor values used for generating
!C  a table of transmittance spectra with different amount of water vapor.
!C  The values in VAPVRT is designed so that there is approximately 2% change
!C  in the .94-um H2O channel ratio for column water vapor in the range .4-6 cm.

      DATA VAPVRT/.00, .02, .06,  .11, .16,  .21,  .26,  .31,  .36, .40, &
                 .43, .46, .50,  .54, .58,  .62,  .66,  .70,  .75, .80, &
                 .86, .92, .98, 1.06,1.14, 1.22,  1.3,  1.4,  1.5, 1.6, &
                 1.7, 1.8, 1.9, 2.05, 2.2, 2.35, 2.55, 2.75, 2.95, 3.2, &
                 3.5, 3.8, 4.1,  4.4, 4.7,  5.0,  5.3,  5.6,  6.0, 6.4, &
                 7.0, 7.7, 8.5,  9.4,10.4, 11.6, 13.0, 15.0, 25.0, 50./

      DATA MD/0,31,59,90,120,151,181,212,243,273,304,334/

      COMMON /GETINPUT3/ H,T,P,VMR,NB,NL,MODEL,IAER,V,TAER55,VRTO3,SNO2

      COMMON /GETINPUT8/ IMN,IDY,IYR,IH,IM,IS
      COMMON /GETINPUT9/ XLATD,XLATM,XLATS,LATHEM
      COMMON /GETINPUT10/XLONGD,XLONGM,XLONGS,LNGHEM
      COMMON /MODEL_ADJ1/ CLMVAP,Q
      COMMON /GEOMETRY1/ SOLZNI,SOLAZ,OBSZNI,OBSPHI,IDAY
      COMMON /GEOMETRY2/ GCO2,GO3,GN2O,GCO,GCH4,GO2,SSH2O,TOTLO3,GGEOM

      COMMON /MODEL_ADJ3/ K_PLANE, DVAP_PLANE, DVAP_LAYER, &
                          DP_PLANE, DP_LAYER, CLMVAPP

      DIMENSION G_VAP(25), G_OTHER(25)
      COMMON /GEOMETRY3/ G_VAP, G_OTHER, G_VAP_EQUIV, VAP_SLANT_MDL

      REAL XPSS, XPPP
      COMMON /GETINPUT14/ XPSS, XPPP
 
      REAL XVIEWD, XVIEWM, XVIEWS
      REAL XAZMUD, XAZMUM, XAZMUS
      COMMON /GETINPUT15/ XVIEWD,XVIEWM,XVIEWS, XAZMUD,XAZMUM,XAZMUS

      HPLANE = XPPP
!C
!C VAP_SLANT is a new array for containing two-way total vapor amounts
!C     Here the number "2" can be changed to other numbers, e.g., 2.5,
!C     without any major effects in derived water vapor values from
!C     imaging spectrometer data.
!C

      DO I = 1, 60
        VAP_SLANT(I) = VAPVRT(I) * 2.0
      END DO

      XH=IH
      XM=IM
      XS=IS
      TT=6.28318*((XH)/24+XM/1440+XS/86400)

      XLAT = XLATD + XLATM/60.0 + XLATS/3600.0
      XLONG = XLONGD + XLONGM/60.0 + XLONGS/3600.0

!C northern hemisphere is positive, and western hemisphere is positive
!C (NOTE: this is not standard, but has been internally coded in this program)
      IF(LATHEM.NE.'N') XLAT = -XLAT
      IF(LNGHEM.NE.'W') XLONG = -XLONG

325   XLATR = XLAT/57.2958
      XLONGR = XLONG/57.2958

      CALL SUNCOR(IDY,IMN,IYR,TT,DEC,HAZ)
      CALL HAZEL(HAZ+TT-XLONGR,DEC,SOLAZ,EL,XLATR)
!C---Note: DEC, SOLAZ,and EL DEC are in radiant units

      IF (EL .LE. 0.) THEN
        WRITE(*,*)'ERROR: Sun is below the horizon!!!'
        WRITE(*,*)'Check input date, time, latitude and longitude.'
        STOP
      ENDIF

!C converting observational zenith angle and relative azimuth angle into degrees
      OBSZNI = XVIEWD + XVIEWM/60.0 + XVIEWS/3600.0
      OBSPHI = XAZMUD + XAZMUM/60.0 + XAZMUS/3600.0

      write(91,*) 'OBSZNI = ',OBSZNI, ' OBSPHI = ',OBSPHI, ' degrees'

      OBSZNI = OBSZNI / 57.2958
      OBSPHI = OBSPHI / 57.2958

      SOLZNI = 90.0/57.2958 - EL

      write(91,*) 'solzni=',solzni
      write(91,*) 'KPLANE=',k_plane

      GGEOM  = 1./COS(SOLZNI) + 1./COS(OBSZNI)

      GCO2   = GGEOM

      GO3    = GGEOM
      IF(HPLANE.LT.27.) GO3 = GGEOM - 1./COS(OBSZNI)

      GN2O   = GGEOM
      GCO    = GGEOM
      GCH4   = GGEOM
      GO2    = GGEOM

      write(91,*) 'GGEOM, GCO2, GO3, GN2O, GCO, GCH4, GO2 = '
      write(91,*)  GGEOM, GCO2, GO3, GN2O, GCO, GCH4, GO2

      TOTLO3 = GO3 * VRTO3

      WRITE(91,*) 'TOTLO3 = ', TOTLO3

!C Initialize newly created geometrical factors for each atmospheric
!C   layers (here G_VAP and G_OTHER are true geometrical factors that
!C   account for the actual Sun-surface-plane path lengths in the
!C   model atmosphere):
!C
!C---For layers below the plane layer---
      DO I = 1, K_PLANE - 1
         G_VAP(I)   = GGEOM
         G_OTHER(I) = GGEOM
      END DO
!C
!C---For layers above the plane layer
      DO I = K_PLANE + 1, 25
         G_VAP(I)   = GGEOM - 1./COS(OBSZNI)
         G_OTHER(I) = GGEOM - 1./COS(OBSZNI)
      END DO
!C
!C---Special treatment for the plane layer to take account the
!C     "shorter" upward path length
      G_VAP(K_PLANE)   = GGEOM - 1./COS(OBSZNI)  &
                         + DVAP_PLANE/DVAP_LAYER/COS(OBSZNI)
      G_OTHER(K_PLANE) = GGEOM - 1./COS(OBSZNI)  &
                         + DP_PLANE/DP_LAYER/COS(OBSZNI)

     write(91,*) ' G_VAP, G_OTHER, I ='
      DO I = 1, 25
         write(91,*) G_VAP(I), G_OTHER(I), I
      END DO

!C Calculate the water vapor SCALING factor relative to the total amount
!C    of water vapor in the model atmosphere in the L-shaped
!C    Sun-surface-plane ray path.
!C
      VAP_SLANT_MDL = CLMVAP/COS(SOLZNI) + CLMVAPP/COS(OBSZNI)
      write(91,*) 'VAP_SLANT_MDL =', VAP_SLANT_MDL, ' cm'
!C
!C The "equivalent" geometrical factor corresponding to the total
!C    slant vapor amount of VAP_SLANT_MDL':
!C
      G_VAP_EQUIV = VAP_SLANT_MDL / CLMVAP
       write(91,*) 'G_VAP_EQUIV = ', G_VAP_EQUIV

      DO 310 I=1,60
        SSH2O(I) = VAP_SLANT(I) / VAP_SLANT_MDL
        write(91,*) 'SSH2O(I), I = ', SSH2O(I), I
  310 CONTINUE

!C Calculate the number of days that have passed in this year.  Take leap year
!C into account.
      IDAY = MD(IMN) + IDY
      LPYR = IYR - (4 * (IYR/4))
      IF((LPYR.EQ.0).AND.(IDAY.GT.59).AND.(IMN.NE.2)) IDAY = IDAY + 1

!C--D     WRITE(*,*)'SOLZNI=',SOLZNI,' SOLAZ=',SOLAZ,' OBSZNI=',OBSZNI
!C--D     WRITE(*,*)'OBSPHI=',OBSPHI,' IDAY = ',IDAY
!C--D     WRITE(*,*)'GCO2=',GCO2,' GO3=',GO3,' GN2O=',GN2O,' GCO=',GCO
!C--D     WRITE(*,*)'GCH4=',GCH4,' GO2=',GO2,' TOTLO3=',TOTLO3
!C--D     WRITE(*,*)'GGEOM=',GGEOM
!C--D     DO 311 I=1,60
!C--D 311   WRITE(*,*)'I=',I,' SSH2O(I)',SSH2O(I)

      RETURN
      END

!********************************************************************************
!*									       *
!*  Name: INIT_SPECCAL							       *
!*  Purpose: initialize global data for spectrum calculations.		       *
!*  Parameters: none							       *
!*  Algorithm: initialize data.   					       *
!*  Globals used: AH2O, APH2O, BH2O, BPH2O, SODLT, SOGAM, O3CF - Band model     *
!*                             parameters for spectral calculations.            *
!*                WNDOW1, WNDOW2, WP94C, WNDOW3, WNDOW4, W1P14C - Center        *
!*                             positions of window and water vapor absorption   *
!*                             channels used in 3-channel ratio calculations.   *
!*                NB1, NB2, NBP94, NB3, NB4, NB1P14 - Number of narrow channels *
!*                             to form broader window and absorption channels.  *
!*  Global output:							       *
!*     IH2OLQ,RLQAMT,NGASTT,NH2O,VSTART,VEND - Flag for including liquid	       *
!*                     water, liquid water amount (cm), total number of gases   *
!*                     (typically 8), number of water vapor values, starting    *
!*                     and ending wavelengths in internal calculations.         *
!*     NO3PT,NCV,NCVHAF,NCVTOT,VMIN,ISTART,IEND	 - Number of O3 abs. coef.     *
!*                     points, parameters for gaussian function and spectral    *
!*                     calculations.					       *
!*     ISTCAL,IEDCAL,DP,PM,TM,VMRM - Parameters for spectral calculations       *
!*     IST1,IED1,IST2,IED2,ISTP94,IEDP94 - 3-channel ratioing parameters	for    *
!*                     the 0.94-um water vapor band.  			       *
!*     IST3,IED3,IST4,IED4,IST1P14,IED1P14 - 3-channel ratioing parameters for  *
!*                     the 1.14-um water vapor band.  			       *
!*     WT1,WT2,WT3,WT4,JA - Relative weights for the four window channels       *
!*                     used in channel-ratioing calculations. JA is a	       *
!*                     output parameter from a table searching routine.	       *
!*     NCV2,NCVHF2,NCVTT2,ISTRT2,IEND2,FINST2 - Parameters for smoothing	       *
!*                     output reflectance spectra.                 	       *
!*     NATOT,NBTOT,NCTOT,NDTOT - Number of channels for the four AVIRIS'	       *
!*                     grating spectrometers (A, B, C, and D).		       *
!*  Return Codes: None.							       *
!*  Special Considerations:  Some parameters may need to be fine-tuned.	       *
!*									       *
!********************************************************************************
!*  Notes about water vapor VMRS and related quantities:		 	       *
!*									       *
!*     VAPVRT(60)   - a table containing 60 column vapor values (in unit of cm) *
!*									       *
!*     VAP_SLANT(I) = VAPVRT(I) * 2.0, VAP_SLANT is a new table for containing  *
!*                    two-way total vapor amounts. Here the number "2" can be   *
!*                    changed to other numbers, e.g., 2.5, without major        *
!*                    effects on retrieved water vapor values.                  *
!*									       *
!*     G_VAP(I = 1,..., NL) = true vapor geometric factor for each layer in     *
!*                    the model atmosphere (after adjusting for the elevated    *
!*                    surface.                                                  *
!*									       *
!*     VMRM(I) = VMRM(I)*G_VAP(I). The VMRS are multiplied by the geometrical   *
!*                    factor. We can calculate the vapor transmittance on the   *
!*                    Sun-surface-sensor path by assuming a vertical path in    *
!*                    the model atmosphere with geometric-factor-adjusted VMRS. *
!*									       *
!*     CLMVAP  = vertical column amount from ground to space in model atmosphere*
!*     CLMVAPP = vertical column amount from ground to aircraft or satellite    *
!*                    sensor in model atmosphere                                *
!*     Q       = 2.152E25 = # of molecules above the surface at one atmosphere  *
!*                    (in unit of  molecules/cm**2)			       *
!*									       *
!*     VAP_SLANT_MDL= CLMVAP/COS(SOLZNI) + CLMVAPP/COS(OBSZNI) = total amount   *
!*                    of water vapor in the model atmosphere in the L-shaped    *
!*                    Sun-surface-plane ray path.                               *
!*									       *
!*     G_VAP_EQUIV  = VAP_SLANT_MDL / CLMVAP = the "equivalent" geometrical     *
!*                    factor corresponding to the total slant vapor amount      *
!*                    VAP_SLANT_MDL and the column vapor amount CLMVAP.         *
!*									       *
!*     SSH2O(I) (I = 1, ..., 60) - a pure scaling factor relative to the total  *
!*                    slant vapor amount of VAP_SLANT_MDL, and                  *
!*		     SSH2O(I) = VAP_SLANT(I) / VAP_SLANT_MDL		       *
!*									       *
!*     SH2O = one value of SSH2O(I). SH2O is used during generation of the      *
!*		     look-up table.       				       *
!*									       *
!*     VAPTT  = VAP_SLANT_MDL*SH2O, is the absolute total vapor amount on the   *
!*                    L-shaped path corresponding to a spectrum stored in the   *
!*                    look-up table.                               	       *
!*									       *
!*     CLMWVP = 0.5*(VAPTTA+VAPTTB)/G_VAP_EQUIV, is the retrieved column water  *
!*                    vapor amount from imaging spectrometer data.	       *
!********************************************************************************
      SUBROUTINE INIT_SPECCAL

      INCLUDE 'COMMONS_INC.f'

!C  Common variables
      DIMENSION H(25), T(25), P(25), VMR(25)
      DIMENSION SSH2O(60)
      DIMENSION WAVOBS(1024),FWHM(1024)
      DIMENSION DP(25), PM(25), TM(25), VMRM(25)
      DIMENSION FINST2(100)

!C  Local variables

      COMMON /GETINPUT1/ IH2OVP,ICO2,IO3,IN2O,ICO,ICH4,IO2,INO2
      COMMON /GETINPUT3/ H,T,P,VMR,NB,NL,MODEL,IAER,V,TAER55,VRTO3,SNO2
      COMMON /GETINPUT4/ WAVOBS,FWHM
      COMMON /GETINPUT5/ NOBS,HSURF,DLT,DLT2
      COMMON /GETINPUT6/ WNDOW1,WNDOW2,WP94C,WNDOW3,WNDOW4,W1P14C
      COMMON /GETINPUT7/ NB1,NB2,NBP94,NB3,NB4,NB1P14
      COMMON /GEOMETRY2/ GCO2,GO3,GN2O,GCO,GCH4,GO2,SSH2O,TOTLO3,GGEOM
      COMMON /MODEL_ADJ1/ CLMVAP,Q

      COMMON /INIT_SPECCAL3/ NH2O
      COMMON /INIT_SPECCAL5/ DP,PM,TM,VMRM
      COMMON /INIT_SPECCAL6/ IST1,IED1,IST2,IED2,ISTP94,IEDP94
      COMMON /INIT_SPECCAL7/ IST3,IED3,IST4,IED4,IST1P14,IED1P14
      COMMON /INIT_SPECCAL8/ WT1,WT2,WT3,WT4,JA

      COMMON /INIT_SPECCAL10/ NCV2,NCVHF2,NCVTT2,ISTRT2,IEND2,FINST2
      COMMON /INIT_SPECCAL11/ NATOT,NBTOT,NCTOT,NDTOT

      DIMENSION G_VAP(25), G_OTHER(25)
      COMMON /GEOMETRY3/ G_VAP, G_OTHER, G_VAP_EQUIV, VAP_SLANT_MDL
 
      DIMENSION O3CF(5001)
      COMMON /O3CF_INIT1/ O3CF
 
      DIMENSION TRAN_O3_STD(5001)
      COMMON /INIT_SPECCAL16/ TRAN_O3_STD
 
      DIMENSION RNO2CF(5001)
      COMMON /NO2CF_INIT1/ RNO2CF
 
      DIMENSION TRAN_NO2_STD(5001)
      COMMON /INIT_SPECCAL17/ TRAN_NO2_STD

      COMMON /MODEL_ADJ4/ K_SURF
      INTEGER start(2)/1,1/
      INTEGER cnt(2)/NP_HI,19/


      NH2O = 60         !Number of column water vapor values
!C
!C Initialization for high resolution spectral calculations -
!C First initialize arrays to smooth high resolution spectra (0.05 cm-1) to
!C    medium resolution spectra (0.2 nm):
!C
!C*** Note: The array WAVNO_HI is not used in actual computing, the array
!C*          should be removed from COMMONS_INC and the following DO LOOP
!C*          of this program, the purpose of keeping WAVNO_HI now is to
!C*          check array indices and make sure the correctness of the
!C*          indices ****************************************************
      DO I = 1, NP_HI
         WAVNO_HI(I)  = 3000. + FLOAT(I-1)*DWAVNO   ! Wavenumber (cm-1)of high 
                                                    !    resolution spectrum,
         TRAN_HI(I)   = 1.0                         !    3000 - 18000 cm-1
      END DO
!C---        print*,'WAVNO_HI ,1, 10000,NP_HI = ',WAVNO_HI(1),
!C---     &     WAVNO_HI(10000), WAVNO_HI(NP_HI)
!C
!C
      DO I = 1, NP_MED
         WAVLN_MED(I) = VSTART + FLOAT(I-1)*DWAVLN  ! Wavelength of medium 
                                                    ! resolution spectrum,
         TRAN_MED(I)  = 1.0                         ! FWHM=.2 nm, .56-3.1 um.
      END DO
!C-----
!C---        print*,'WAVLN_MED ,1, 10000,NP_MED = ',WAVLN_MED(1),
!C---     &     WAVLN_MED(10000), WAVLN_MED(NP_MED)
!C-----
!C
!C  NOTE:  WAVLN_STD starts from 0.3 um, instead of 0.56 um
      DO I = 1, NP_STD
         WAVLN_STD(I)  = 0.3   + FLOAT(I-1)*DWAVLN  ! Wavelength of medium 
                                                    ! resolution spectrum,
         TRAN_STD(I)   = 1.0                        ! FWHM=.2 nm, .3-3.1 um.
      END DO
!C-----
!C---        print*,'WAVLN_STD ,1, 10000,NP_STD = ',WAVLN_STD(1),
!C---     &     WAVLN_STD(10000), WAVLN_STD(NP_STD)
!C-----
!C
!C Note: The grids of WAVNO_HI do not match the grids of 10000./WAVLN_MED.
!C       INDEX_MED is a specially designed index for finding close matches
!C       between the two kinds of grids.
!C
      DO I = 1, NP_MED  
         INDEX_MED(I) = ( (10000./WAVLN_MED(I) - 3000.)/DWAVNO + 1.)
      END DO
!C-----
!C---        print*,'INDEX_MED ,1, 10000,NP_MED = ',INDEX_MED(1),
!C---     &     INDEX_MED(10000), INDEX_MED(NP_MED)
!C-----
!C
!C Note:     WAVLN_MED_INDEX(I) is very close to WAVLN_MED(I),
!C       and WAVLN_MED_INDEX(I) >= WAVLN_MED(I)
!C
      DO I = 1, NP_MED
         WAVLN_MED_INDEX(I) = 10000. /(FLOAT(INDEX_MED(I)-1)*DWAVNO &
                                     + 3000.)
      END DO
!C-----
!C---        print*,'WAVLN_MED_INDEX ,1, 10000,NP_MED = ',WAVLN_MED_INDEX(1),
!C---     &     WAVLN_MED_INDEX(10000), WAVLN_MED_INDEX(NP_MED)
!C-----


      DO I = 1, NP_MED
         FWHM_WAVNO(I) = 10000.*DLT_MED &
                         /(WAVLN_MED_INDEX(I)*WAVLN_MED_INDEX(I))
      END DO
!C-----
!C---        print*,'FWHM_WAVNO ,1, 10000,NP_MED = ',FWHM_WAVNO(1),
!C---     &     FWHM_WAVNO(10000), FWHM_WAVNO(NP_MED)
!C-----
!C

      DO I = 1, NP_MED
         NCVHF_WAVNO(I) = ( FACDLT * FWHM_WAVNO(I) / DWAVNO + 1.)
      END DO
!C-----
!C---        print*,'NCVHF_WAVNO ,1, 10000,NP_MED = ',NCVHF_WAVNO(1),
!C---     &     NCVHF_WAVNO(10000), NCVHF_WAVNO(NP_MED)
!C-----

!C Initialize arrays for smoothing medium resolution spectrum (DLT_MED = 0.2 nm,
!C         and point spacing DWAVLN = 0.0001 micron) to coarser spectral
!C         resolution data from imaging spectrometers.

      DO I = 1, NOBS
         NCVHF(I) = ( FACDLT * FWHM(I) / DWAVLN + 1.)
      END DO
!C
!C parameters and arrays to smooth output surface reflectance spectrum
      WAVCV2=FACDLT*DLT2

!C Find the largest value in the FWHM array, and use this value in calculation
!C    of indices for smoothing output reflectance spectra. This smoothing
!C    algorithm should work well with grating spectrometers having nearly
!C    constant spectral resolutions, but not so well for prism spectrometers
!C    having variable spectral resolution.
      DWVAVR = FWHM(1)

      DO I = 2, NOBS
         IF(DWVAVR.LT.FWHM(I)) DWVAVR = FWHM(I)
      END DO

      RNCV2=WAVCV2/DWVAVR
      NCV2=RNCV2
      NCVHF2=NCV2+1
      NCVTT2=2*NCV2+1

      CONS2=DLT2*SQRT(3.1415926/CONST1)

      IF (DLT2 .NE. 0.0) THEN
        SUMINS=0.0
        DO 585 I=NCVHF2,NCVTT2
          FINST2(I)=EXP(-CONST1*(FLOAT(I-NCVHF2)*DWVAVR/DLT2)**2.0)
          SUMINS=SUMINS+FINST2(I)
  585   CONTINUE

        DO 590 I=1,NCVHF2-1
          FINST2(I)=FINST2(NCVTT2-I+1)
          SUMINS=SUMINS+FINST2(I)
  590   CONTINUE

        SUMINS=SUMINS*DWVAVR

        DO 595 I=1,NCVTT2
          FINST2(I)=FINST2(I)*DWVAVR/SUMINS
  595   CONTINUE
      ENDIF

      ISTRT2=NCVHF2
      IEND2=NOBS-NCVHF2

!C  number of channels of the four AVIRIS spectrometers.  These are used
!C  in removing null AVIRIS radiance values in the overlap portions of two
!C  adjacent spectrometers.
      NCHNLA=32
      NCHNLB=64
      NCHNLC=64
      NCHNLD=64

      NATOT=NCHNLA
      NBTOT=NCHNLA+NCHNLB
      NCTOT=NCHNLA+NCHNLB+NCHNLC
      NDTOT=NCHNLA+NCHNLB+NCHNLC+NCHNLD

!C Resetting window wavelength positions and calculating weights for
!C  window and absorption channels used in 3-channel ratioing.
      IWNDW1=FINDMATCH(WAVOBS,NOBS,WNDOW1)
      IWNDW2=FINDMATCH(WAVOBS,NOBS,WNDOW2)

      WNDOW1=WAVOBS(IWNDW1)
      WNDOW2=WAVOBS(IWNDW2)

      JJ=MOD(NB1,2)
      IF(JJ.EQ.0) NB1=NB1+1
      KK=MOD(NB2,2)
      IF(KK.EQ.0) NB2=NB2+1
      NB1HAF=(NB1-1)/2
      NB2HAF=(NB2-1)/2

      IST1=IWNDW1-NB1HAF
      IED1=IWNDW1+NB1HAF
      IST2=IWNDW2-NB2HAF
      IED2=IWNDW2+NB2HAF

      IWP94C=FINDMATCH(WAVOBS,NOBS,WP94C)
      WP94C=WAVOBS(IWP94C)

      LL=MOD(NBP94,2)
      IF(LL.EQ.0) NBP94=NBP94+1
      NB3HAF=(NBP94-1)/2
      ISTP94=IWP94C-NB3HAF
      IEDP94=IWP94C+NB3HAF

      WT1=(WNDOW2-WP94C)/(WNDOW2-WNDOW1)
      WT2=(WP94C-WNDOW1)/(WNDOW2-WNDOW1)

      IWNDW4=FINDMATCH(WAVOBS,NOBS,WNDOW3)
      IWNDW5=FINDMATCH(WAVOBS,NOBS,WNDOW4)

      WNDOW3=WAVOBS(IWNDW4)
      WNDOW4=WAVOBS(IWNDW5)

      JJ=MOD(NB3,2)
      IF(JJ.EQ.0) NB3=NB3+1
      KK=MOD(NB4,2)
      IF(KK.EQ.0) NB4=NB4+1

      NB4HAF=(NB3-1)/2
      NB5HAF=(NB4-1)/2

      IST3=IWNDW4-NB4HAF
      IED3=IWNDW4+NB4HAF
      IST4=IWNDW5-NB5HAF
      IED4=IWNDW5+NB5HAF
      IW1P14C=FINDMATCH(WAVOBS,NOBS,W1P14C)

      W1P14C=WAVOBS(IW1P14C)
      LL=MOD(NB1P14,2)
      IF(LL.EQ.0) NB1P14=NB1P14+1
      NB6HAF=(NB1P14-1)/2
      IST1P14=IW1P14C-NB6HAF
      IED1P14=IW1P14C+NB6HAF

      WT3=(WNDOW4-W1P14C)/(WNDOW4-WNDOW3)
      WT4=(W1P14C-WNDOW3)/(WNDOW4-WNDOW3)

!C Initialization for searching water vapor table (TRNTBL)
      JA = 30
!C
!C
!C Calculate medium resolution O3 transmittances (0.2 nm, 2-way) with
!C     a point spacing of 0.1 nm between 0.3 and 0.8 micron.
!C
      IF(IO3.EQ.1) THEN
         DO I=1,NO3PT
            TRAN_O3_STD(I) = EXP(-TOTLO3*O3CF(I))
         END DO
      END IF

!C
!C If O3 is not intended to be included in total atmospheric gaseous
!C     transmittance calculations, assigning TRAN_O3_STD = 1.0:
!C
      IF(IO3.NE.1) THEN
         DO I=1,NO3PT
            TRAN_O3_STD(I) = 1.0
         END DO
      END IF
!C
!C Calculate medium resolution NO2 transmittances (0.2 nm, 2-way) with
!C     a point spacing of 0.1 nm between 0.3 and 0.8 micron.
!C
         NO2PT  = NO3PT
         VRTNO2 = 5.0E+15
         VRTNO2 = SNO2 * VRTNO2

         GNO2   = GO3
         TOTNO2 = GNO2 * VRTNO2

      IF(INO2.EQ.1) THEN
         DO I=1,NO2PT
            TRAN_NO2_STD(I) = EXP(-TOTNO2*RNO2CF(I))
         END DO
      END IF

!C---Temp Code:
!C---      OPEN(57,file='zzzzzz_tst_NO2_trn',status='unknown')
!C---         DO I = 1, NO2PT
!C---          write(57,*)WAVLN_STD(I), TRAN_NO2_STD(I)
!C---         END DO
!C---      CLOSE(57)
!C---End of Temp Code
!
!C---      print*,' TOTNO2 = ', TOTNO2
!C---      print*,'VRTNO2, SNO2, GNO2 = ', VRTNO2, SNO2, GNO2

!C
!C If NO2 is not intended to be included in total atmospheric gaseous
!C     transmittance calculations, assigning TRAN_NO2_STD = 1.0:
!C
      IF(INO2.NE.1) THEN
         DO I=1,NO2PT
            TRAN_NO2_STD(I) = 1.0
         END DO
      END IF


!C
!C Initial arrays for mean layer pressure and temperatures

      DO 320 I=1,NL
        DP(I)=P(I)-P(I+1)
        PM(I)=(P(I)+P(I+1))/2.0
        TM(I)=(T(I)+T(I+1))/2.0
  320 CONTINUE

!C
!C Calculate high resolution transmittances (0.05 cm-1) of CO2, N2O, CO,
!C     CH4, and O2 in the 0.56 - 3.1 micron range, and save values for
!C     calculating total atmospheric transmittances later.
!C     Because water vapor amounts are allowed to vary,
!C     the high resolution water vapor transmittances are calculated
!C     in subroutines TRAN_TABLE and TRANCAL. TRAN_TABLE provides variable
!C     water vapor amounts, and calls TRANCAL for the calculation of
!C     corresponding vapor transmittance spectrum.
!C
!C*** Note: The array WAVNO_HI is not used in actual computing, the array
!C*          should be removed from COMMONS_INC and the following DO LOOP
!C*          of this program, the purpose of keeping WAVNO_HI now is to
!C*          check array indices and make sure the correctness of the
!C*          indices.
!C*
!C Initialize the TRAN_HI_OTHERS array for high resolution spectrum:
        DO I = 1, NP_HI
           TRAN_HI_OTHERS(I) = 1.0
        END DO

       ncid = ncopn('abscf_gas.nc',NCNOWRIT,IRCODE)

!C---------------------------------------
!C For CO2 transmittance calculation -
        IF(ICO2.EQ.1) THEN
!C----        SCLCO2=1.0
!C----  On 2/7/2013 B.-C. Gao made the modification - Increased SCLCO2 i
!C----     from 1.0 to 1.1 to reflect the fact that the CO2 VMR reached the
!C----     2012 level of 390 ppmv.
!C----
          SCLCO2=1.1

          DO 322 I=1,NL
            VMRM(I)=SCLCO2*355.0*1.0E-06
!C Scale the VMRM by the two-way path geometrical factors. The geometric
!C           factors, G_OTHER, varies with atmospheric layer number for
!C           aircraft observational geometries.
            VMRM(I)= VMRM(I)*G_OTHER(I)
  322     CONTINUE
  
  
       
!        NRHID = NCVID (NCID, 'waveno_hi_co2', IRCODE)
!       CALL NCVGT (NCID, NRHID, START(1), CNT(1), WAVNO_HI, IRCODE)
!       if (Ircode .ne.0) then
!          write(*,*) 'Error reading abscf_gas.nc: wavno_hi: rcode=',ircode
!          stop 0
!       end if
       NRHID = NCVID (NCID, 'abscf_co2', IRCODE)
       CALL NCVGT (NCID, NRHID, START, CNT, ABSCF, IRCODE)
        if (ircode .ne.0) then
          write(*,*) 'Error reading abscf_gas.nc: abscf_co2: rcode=',ircode
          stop 0
       end if

!       OPEN(31,file='abscf_co2_PC',status='old', &
!         form='unformatted',access='direct',recl=4*300000)

!          READ(31,REC=1) (WAVNO_HI(I), I = 1, 300000)
!
       DO J = K_SURF, 19
!          READ(31,REC=J+1) (A(I), I = 1, 300000)
!
             TRAN_HI_OTHERS(:) = TRAN_HI_OTHERS(:) * EXP( - ABSCF(:,J)*Q* &
               DP(J-K_SURF+1)*VMRM(J-K_SURF+1) * 28.966 / &
               6.0225E+23 / 1.0E-06)
       END DO
!
!       CLOSE(31)
      ENDIF

!C--------------------------------------------
!C       For N2O transmittance calculation.

        IF(IN2O.EQ.1) THEN
          DO 324 I=1,NL
            VMRM(I)=0.3*1.0E-06
            VMRM(I)= VMRM(I)*G_OTHER(I)
  324     CONTINUE
  
!        NRHID = NCVID (NCID, 'waveno_hi_n2o', IRCODE)
!       CALL NCVGT (NCID, NRHID, START(1), CNT(1), WAVNO_HI, IRCODE)
!       if (Ircode .ne.0) then
!          write(*,*) 'Error reading abscf_gas.nc: wavno_hi: rcode=',ircode
!          stop 0
!       end if
       NRHID = NCVID (NCID, 'abscf_n2o', IRCODE)
       CALL NCVGT (NCID, NRHID, START, CNT, ABSCF, IRCODE)
        if (ircode .ne.0) then
          write(*,*) 'Error reading abscf_gas.nc: abscf_n2o: rcode=',ircode
          stop 0
       end if

!       OPEN(31,file='abscf_n2o_PC',status='old', &
!          form='unformatted',access='direct',recl=4*300000)

!          READ(31,REC=1) (WAVNO_HI(I), I = 1, 300000)
        
       DO J = K_SURF, 19
!          READ(31,REC=J+1) (A(I), I = 1, 300000)

             TRAN_HI_OTHERS(:) = TRAN_HI_OTHERS(:) * EXP( - ABSCF(:,J)*Q* &
                DP(J-K_SURF+1)*VMRM(J-K_SURF+1) * 28.966 / &
                6.0225E+23 / 1.0E-06)
       END DO

!       CLOSE(31)
      ENDIF

!C--------------------------------------------
!C       For CO transmittance calculation.
        IF(ICO.EQ.1) THEN
          DO 325 I=1,NL
            VMRM(I)=0.1*1.0E-06
            VMRM(I)= VMRM(I)*G_OTHER(I)
  325     CONTINUE
  
!        NRHID = NCVID (NCID, 'waveno_hi_co', IRCODE)
!       CALL NCVGT (NCID, NRHID, START(1), CNT(1), WAVNO_HI, IRCODE)
!       if (Ircode .ne.0) then
!          write(*,*) 'Error reading abscf_gas.nc: wavno_hi: rcode=',ircode
!          stop 0
!       end if
       NRHID = NCVID (NCID, 'abscf_co', IRCODE)
       CALL NCVGT (NCID, NRHID, START, CNT, ABSCF, IRCODE)
        if (ircode .ne.0) then
          write(*,*) 'Error reading abscf_gas.nc: abscf_co: rcode=',ircode
          stop 0
       end if

!       OPEN(31,file='abscf_co_PC',status='old', &
!         form='unformatted',access='direct',recl=4*300000)

 !         READ(31,REC=1) (WAVNO_HI(I), I = 1, 300000)

       DO J = K_SURF, 19
 !         READ(31,REC=J+1) (A(I), I = 1, 300000)

             TRAN_HI_OTHERS(:) = TRAN_HI_OTHERS(:) * EXP( - ABSCF(:,J)*Q* &
               DP(J-K_SURF+1)*VMRM(J-K_SURF+1) * 28.966 / &
               6.0225E+23 / 1.0E-06)
       END DO

 !      CLOSE(31)
      ENDIF

!C--------------------------------------------
!C       For CH4 transmittance calculation.
!C For assigning CH4 VMRM
!C NOTE: The scaling factor of 0.8 for the CH4 VMRS was obtained by comparing
!C       transmittance spectra calculated using our program, which is based on
!C       the Malkmus narrow band spectral model, with a ratioed spectrum
!C       provided by G. C. Toon at Jet Propulsion Laboratory (JPL). The JPL
!C       ratio spectrum was obtained by ratioing a high resolution (0.005 cm-1)
!C       solar spectrum measured at ground level against a solar spectrum
!C       measured at an altitude of approximately 35 km with the same Fourier
!C       Transform spectrometer. The high resolution ratio spectrum was
!C       degraded to a resolution of 10 nm during our derivation of the
!C       scaling factor for the CH4 VMRS.

        IF(ICH4.EQ.1) THEN
!C***          SCLCH4=0.8
          SCLCH4=1.0
          DO 326 I=1,NL
            VMRM(I)=SCLCH4*1.6*1.0E-06
            VMRM(I)= VMRM(I)*G_OTHER(I)
  326     CONTINUE

!        NRHID = NCVID (NCID, 'waveno_hi_ch4', IRCODE)
!       CALL NCVGT (NCID, NRHID, START(1), CNT(1), WAVNO_HI, IRCODE)
!       if (Ircode .ne.0) then
!          write(*,*) 'Error reading abscf_gas.nc: wavno_hi: rcode=',ircode
!          stop 0
!       end if
       NRHID = NCVID (NCID, 'abscf_ch4', IRCODE)
       CALL NCVGT (NCID, NRHID, START, CNT, ABSCF, IRCODE)
        if (ircode .ne.0) then
          write(*,*) 'Error reading abscf_gas.nc: abscf_ch4: rcode=',ircode
          stop 0
       end if

!       OPEN(31,file='abscf_ch4_PC',status='old', &
!          form='unformatted',access='direct',recl=4*300000)

!          READ(31,REC=1) (WAVNO_HI(I), I = 1, 300000)

       DO J = K_SURF, 19
!          READ(31,REC=J+1) (A(I), I = 1, 300000)

              TRAN_HI_OTHERS(:) = TRAN_HI_OTHERS(:) * EXP( - ABSCF(:,J)*Q* &
               DP(J-K_SURF+1)*VMRM(J-K_SURF+1) * 28.966 / &
               6.0225E+23 / 1.0E-06)
       END DO

!       CLOSE(31)
      ENDIF

!C--------------------------------------------
!C       For O2 transmittance calculation.
        IF(IO2.EQ.1) THEN

!C***Modified by Bo-Cai Gao on 2/7/2013 to increase O2 absorption
!C---  coefficients by the factor SCL_O2 for wavelengths > 1.2 micron
!C---  in order to model properly the atmospheric O2 band centered
!C---  near 1.265 micron.

          SCL_O2  = 2.60

          DO 327 I=1,NL
            VMRM(I)=0.21
            VMRM(I)= VMRM(I)*G_OTHER(I)
  327     CONTINUE

!        NRHID = NCVID (NCID, 'waveno_hi_o2', IRCODE)
!       CALL NCVGT (NCID, NRHID, START(1), CNT(1), WAVNO_HI, IRCODE)
!       if (Ircode .ne.0) then
!          write(*,*) 'Error reading abscf_gas.nc: wavno_hi: rcode=',ircode
!          stop 0
!       end if
       NRHID = NCVID (NCID, 'abscf_o2', IRCODE)
       CALL NCVGT (NCID, NRHID, START, CNT, ABSCF, IRCODE)
        if (ircode .ne.0) then
          write(*,*) 'Error reading abscf_gas.nc: abscf_o2: rcode=',ircode
          stop 0
       end if

!       OPEN(31,file='abscf_o2_PC',status='old', &
!         form='unformatted',access='direct',recl=4*300000)

!          READ(31,REC=1) (WAVNO_HI(I), I = 1, 300000)

       DO J = K_SURF, 19
!          READ(31,REC=J+1) (A(I), I = 1, 300000)

!C***Modified by Bo-Cai Gao on 2/7/2013 to increase O2 absorption
!C---  coefficients by the factor SCL_O2 for wavelengths > 1.2 micron
!C---  i& < 1.3333 micron in order to model properly the atmospheric
!C---  O2 band centered near 1.265 micron.

             TRAN_HI_OTHERS(1:9000) = TRAN_HI_OTHERS(1:9000) * EXP( - ABSCF(1:9000,J)*Q* &
               DP(J-K_SURF+1)*VMRM(J-K_SURF+1) * 28.966 / &
               6.0225E+23 / 1.0E-06)

             TRAN_HI_OTHERS(9001:106600) = TRAN_HI_OTHERS(9001:106600) * EXP( - ABSCF(9001:106600,J)*Q* &
             DP(J-K_SURF+1)*VMRM(J-K_SURF+1)*SCL_O2 *28.966/ &
               6.0225E+23 / 1.0E-06)
 
             TRAN_HI_OTHERS(106601:NP_HI) = TRAN_HI_OTHERS(106601:NP_HI) * EXP( - ABSCF(106601:NP_HI,J)*Q* &
               DP(J-K_SURF+1)*VMRM(J-K_SURF+1) * 28.966 / &
               6.0225E+23 / 1.0E-06)
       END DO

       CLOSE(31)
      ENDIF

       CALL NCCLOS(NCID, RCODE)

!C--------------------------------------------
!C End of calculation of high resolution transmittances for CO2, N2O, CO,
!C     CH4, and O2.
!C--------------------------------------------
!C Initial water vapor VMRs for repeated use in other subroutines and
!C    adjust layered water vapor VMRM with geometrical factors.
      DO I = 1, NL
         VMRM(I) = (VMR(I)+VMR(I+1))/2.0
!C Scale the VMRM by the two-way path geometrical factors. The geometric
!C           factors, G_VAP, varies with atmospheric layer number for
!C           aircraft observational geometries.
         VMRM(I) = VMRM(I)*G_VAP(I)
      END DO
!C--------------------------------------------
!C
!C--D     WRITE(*,*)'NPSHIF=',NPSHIF,' DWAVLN=',DWAVLN
!C--D     WRITE(*,*)'NO3PT=',NO3PT,' VMIN=',VMIN,' ISTART=',ISTART
!C--D     WRITE(*,*)'IH2OLQ=',IH2OLQ,' RLQAMT=',RLQAMT,' NGASTT=',NGASTT
!C--D     WRITE(*,*)'NH2O=',NH2O,' VSTART=',VSTART,' VEND=',VEND
!C--D     WRITE(*,*)'IEND=',IEND
!C--D     WRITE(*,*)'ISTCAL=',ISTCAL,' IEDCAL=',IEDCAL
!C--D     DO 545 I=1,NL
!C--D 545   WRITE(*,*)I,DP(I),PM(I),TM(I),VMRM(I)
!C--D     WRITE(*,*)'IST1=',IST1,' IED1=',IED1,' IST2=',IST2,' IED2=',IED2
!C--D     WRITE(*,*)'ISTP94=',ISTP94,' IEDP94=',IEDP94
!C--D     WRITE(*,*)'IST3=',IST3,' IED3=',IED3,' IST4=',IST4,' IED4=',IED4
!C--D     WRITE(*,*)'IST1P14=',IST1P14,' IED1P14=',IED1P14
!C--D     WRITE(*,*)'WT1=',WT1,' WT2=',WT2,' WT3=',WT3,' WT4=',WT4,' JA=',JA
!C--D     WRITE(*,*)'NCV2=',NCV2,' NCVHF2=',NCVHF2,' NCVTT2=',NCVTT2,
!C--D    &' ISTRT2=',ISTRT2,' IEND2=',IEND2
!C--D     DO 544 I=1,30
!C--D 544   WRITE(*,*)'I=',I,' FINST2(I)=',FINST2(I)
!C--D     WRITE(*,*)'NATOT=',NATOT,' NBTOT=',NBTOT
!C--D     WRITE(*,*)'NCTOT=',NCTOT,' NDTOT=',NDTOT
     
      RETURN
      END

!********************************************************************************
!*   									       *
!*  Name: TRAN_TABLE							       *
!*  Purpose: This subroutine generates a table consisting of 60 atmospheric     *
!*           transmittance spectra at the solar and observational               *
!*           geometry and with 60 column water vapor values. The table also     *
!*           includes the total amounts of column water vapor used in the       *
!*           calculations, and the 3-channel ratios calculated from the window  *
!*           and absorption channels in and around the 0.94- and 1.14-um water  *
!*           vapor bands.						       *
!*  Parameters: none							       *
!*  Algorithm: For each of the 60 water vapor amounts, calculate the 	       *
!*             atmospheric transmittance, and save 			       *
!*  Globals Used: NH2O   - number of column water vapor values		       *
!*		 VAPTT  - geometrically adjusted water vapor total	       *
!*		 R094   - channel ratio for .94 um region		       *
!*		 R114   - channel ratio for 1.14 um region 		       *
!*		 TRNCAL - atmospheric transmittance spectra		       *
!*  Global Output: VAPTOT() - array containing geometrically adjusted water     *
!*		  	     vapor values				       *
!*		  ROP94()  - array containing channel ratios for .94 um region *
!*		  R1P14()  - array containing channel ratios for 1.14 um region*
!*		  TRNTBL() - 2 dimensional array containing one transmittance  *
!*			     spectrum for each column water vapor amount       *
!*  Return Codes: none                                                          *
!*  Special Considerations: none                                                *
!*									       *
!********************************************************************************
      SUBROUTINE TRAN_TABLE

      INCLUDE 'COMMONS_INC.f'

!C  Common variables
      DIMENSION TRNCAL(1024)
      DIMENSION WAVOBS(1024),FWHM(1024)
      DIMENSION DP(25), PM(25), TM(25), VMRM(25)
      DIMENSION SSH2O(60)
      DIMENSION VAPTOT(60), R0P94(60), R1P14(60), TRNTBL(1024,60)
      
      DIMENSION O3CF(5001)
      COMMON /O3CF_INIT1/ O3CF

      COMMON /GETINPUT1/ IH2OVP,ICO2,IO3,IN2O,ICO,ICH4,IO2,INO2
      COMMON /GETINPUT4/ WAVOBS,FWHM
      COMMON /GETINPUT5/ NOBS,HSURF,DLT,DLT2
      COMMON /INIT_SPECCAL3/ NH2O
      COMMON /INIT_SPECCAL5/ DP,PM,TM,VMRM

      COMMON /TRAN_TABLE1/ SH2O,VAPTOT,R0P94,R1P14,TRNTBL
      COMMON /TRANCAL1/ TRNCAL,VAPTT
      COMMON /GEOMETRY2/ GCO2,GO3,GN2O,GCO,GCH4,GO2,SSH2O,TOTLO3,GGEOM
      COMMON /CHNLRATIO1/ R094,R114

!C For each water vapor amount, calculate the geometrically adjusted water
!C     vapor amount, the channel ratios, and the transmittance spectrum.
      VAPTT = 0.

      DO 530 I=1,NH2O
        SH2O=SSH2O(I)

        CALL TRANCAL

        VAPTOT(I)= VAPTT
        R0P94(I) = R094
        R1P14(I) = R114
        DO 540 J=1,NOBS
          TRNTBL(J,I)=TRNCAL(J)
 540    CONTINUE
 530  CONTINUE

!C Write above calculated values to a file.
      NJ=NH2O
      NK=FLOAT(NJ)/10.
      NJDIFF=NJ-NK*10
 
!RJH      DO 380 KK=1,NK
!        WRITE(35,77) (VAPTOT(I+(KK-1)*10),I=1,10)
!   77   FORMAT(7X,10(1X,E11.4))
!        WRITE(35,77) ( R0P94(I+(KK-1)*10),I=1,10)
!        WRITE(35,77) ( R1P14(I+(KK-1)*10),I=1,10)
!        DO 385 II=1,NOBS
!  385     WRITE(35,78) WAVOBS(II),(TRNTBL(II,I+(KK-1)*10),I=1,10)
!   78     FORMAT(1X,F6.4,10(1X,E11.4))
!  380 CONTINUE

!      WRITE(35,77) (VAPTOT(I+NK*10),I=1,NJDIFF)
!      WRITE(35,77) ( R0P94(I+NK*10),I=1,NJDIFF)
!      WRITE(35,77) ( R1P14(I+NK*10),I=1,NJDIFF)
!
!      IF(NJDIFF.GE.1) THEN
!        DO 386 II=1,NOBS
!  386     WRITE(35,78) WAVOBS(II),(TRNTBL(II,I+NK*10),I=1,NJDIFF)
!      ENDIF
!      CLOSE(35,IOSTAT=ISTAT)

      RETURN
      END

!********************************************************************************
!*									       *
!*  Name: TRANCAL							       *
!*  Purpose: This program calculates combined transmittances of H2O, CO2, O3,   *
!*      N2O, CO, CH4, and O2.                                                   *
!*  Parameters: none.							       *
!*  Algorithm: The calculations were based on the line-by-line absorption       *
!*      parameters supplied by William R. Ridgway of NASA/GSFC.		       *
!*  Global output:VAPTT  - geometrically adjusted water vapor total.	       *
!*		 R094   - channel ratio for 0.94 um region.		       *
!*		 R114   - channel ratio for 1.14 um region.		       *
!*		 TRNCAL - total transmittances of all gases that match the     *
!*                         resolutions of imaging spectrometers.		       *
!*  Return Codes: none.							       *
!*  Special Considerations: The high resolution (0.05 cm-1) line-by-line        *
!*      absorption parameters cover the 0.555 - 3.33 micron spectral range      *
!*      (3000 - 18000 cm-1). The medium resolution ozone absorption             *
!*      coefficients covers the 0.3-0.8 um spectral range. The line-by-line     *
!*      high resolution spectra were first smoothed to medium resolution        *
!*      spectra (resolution = 0.2 nm, wavelength spacing = 0.1 nm) covering     *
!*      the 0.56 - 3.1 micron spectral region. The medium resolution spectra    *
!*      of O3 and other gases are combined (in SUBROUTINE TRAN_SMOOTH) to form  *
!*      a single medium resolution spectrum from 0.3 to 3.1 micron. This        *
!*      combined spectrum (medium resolution) is then smoothed to lower         *
!*      resolutions to match the resolutions of imaging spectrometers. The      *
!*      smoothing is also done in SUBROUTINE TRAN_SMOOTH.                       *
!*									       *
!********************************************************************************

      SUBROUTINE TRANCAL

      INCLUDE 'COMMONS_INC.f'
      INCLUDE 'netcdf.inc'
!C  Common variables
      DIMENSION TRNCAL(1024)
      DIMENSION WAVOBS(1024),FWHM(1024)
      DIMENSION DP(25), PM(25), TM(25), VMRM(25)
      DIMENSION VAPTOT(60), R0P94(60), R1P14(60), TRNTBL(1024,60)

!C  Local variables
      DIMENSION TRANS(1024)
      DIMENSION TRNCV(1024)
      DIMENSION TRNSTD_SM(1050)
      INTEGER   INDEX
      
      COMMON /GETINPUT4/ WAVOBS,FWHM
      COMMON /GETINPUT5/ NOBS,HSURF,DLT,DLT2
      COMMON /MODEL_ADJ1/ CLMVAP,Q
      COMMON /INIT_SPECCAL5/ DP,PM,TM,VMRM
      COMMON /INIT_SPECCAL6/ IST1,IED1,IST2,IED2,ISTP94,IEDP94
      COMMON /INIT_SPECCAL7/ IST3,IED3,IST4,IED4,IST1P14,IED1P14
      COMMON /TRAN_TABLE1/ SH2O,VAPTOT,R0P94,R1P14,TRNTBL
      COMMON /TRANCAL1/ TRNCAL,VAPTT

      DIMENSION G_VAP(25), G_OTHER(25)
      COMMON /GEOMETRY3/ G_VAP, G_OTHER, G_VAP_EQUIV, VAP_SLANT_MDL

      COMMON /MODEL_ADJ4/ K_SURF
      INTEGER FIRST/1/
      INTEGER start(2)/1,1/
      INTEGER cnt(2)/NP_HI,19/
      
!      REAL ABSCF(NP_HI,19)
      SAVE FIRST

!C
!C Calculate water vapor transmittances with different scaling factors:
!C

       DO I = 1, NP_HI
          TRAN_HI(I) = 1.0
       END DO

      if (first.eq.1) then
       ncid = ncopn('abscf_gas.nc',NCNOWRIT,IRCODE)
!        NRHID = NCVID (NCID, 'waveno_hi_h2o', IRCODE)
!       CALL NCVGT (NCID, NRHID, START(1), CNT(1), WAVNO_HI, IRCODE)
!       if (Ircode .ne.0) then
!          write(*,*) 'Error reading abscf_gas.nc: wavno_hi: rcode=',ircode
!          stop 0
!       end if
       NRHID = NCVID (NCID, 'abscf_h2o', IRCODE)
       CALL NCVGT (NCID, NRHID, START, CNT, ABSCF, IRCODE)
        if (ircode .ne.0) then
          write(*,*) 'Error reading abscf_gas.nc: abscf_h2o: rcode=',ircode
          stop 0
       end if
       CALL NCCLOS(NCID, RCODE)
        
!       OPEN(31,file='abscf_h2o_PC',status='old', &
!         form='unformatted',access='direct',recl=4*300000)

!          READ(31,REC=1) (WAVNO_HI(I), I = 1, 300000)

!       DO J = 1, 19
!          READ(31,REC=J+1) (ABSCF(I,J), I = 1, 300000)

!       END DO

!       CLOSE(31)
       FIRST=0
      endif

       DO J = K_SURF, 19

             TRAN_HI(:) = TRAN_HI(:) * EXP( - ABSCF(:,J) * Q * &
               DP(J-K_SURF+1)*VMRM(J-K_SURF+1) * SH2O * 28.966 / &
               6.0225E+23 / 1.0E-06)

       END DO

!C Multiplying the high resolution water vapor transmittance with combined
!C    transmittances of CO2, N2O, CO, CH4, and O2:
!C
       TRAN_HI(:) = TRAN_HI(:)*TRAN_HI_OTHERS(:)
!       DO I = 1, NP_HI
!          TRAN_HI(I) = TRAN_HI(I) * TRAN_HI_OTHERS(I)
!       END DO

!C
!C Smooth the high resolution spectra to resolutions of measured spectrum
!C
      CALL TRAN_SMOOTH
!C
!C---      DO 470 I=1,NOBS
!C---        WRITE(*,*)'I=',I,' TRNCAL(I)=',TRNCAL(I)
!C---  470 CONTINUE

!C Calculate 3-channel ratio values, R094 and R114, from simulated spectra.
      CALL CHNLRATIO

!C Total amount of water vapor (in unit of cm) corresponding to the spectrum.
      VAPTT = VAP_SLANT_MDL * SH2O
!C--       print*, 'VAPTT= ', VAPTT,'VAP_SLANT_MDL= ',VAP_SLANT_MDL

 9999 RETURN
      END

!********************************************************************************
!*                                                                              *
!*  Name: TRAN_SMOOTH                                                           *
!*  Purpose: This program is to smooth the line-by-line high resolution         *
!*           spectrum to lower resolution spectrum that matches the resolutions *
!*           of imaging spectrometer data.                                      *
!*  Parameters: none.                                                           *
!*  Algorithm: The smoothing is done in two stages. The 1st stage is to smooth  *
!*             the high resolution spectrum to medium resolution spectrum at a  *
!*             constant FWHM (0.2 nm) and a constant wavelength interval        *
!*             (0.1 nm). The 2nd stage smoothing is to smooth the medium        *
!*             resolution spectrum to resolutions of input imaging spectrometer *
!*             data.                                                            *
!*  Globals used:  The global variables used are contained in the file          *
!*                         "COMMONS_INC"                                        *
!*  Global output: TRNCAL - total transmittances of all gases that match the    *
!*                          resolutions of imaging spectrometers.               *
!*  Return Codes:  none.                                                        *
!*                                                                              *
!********************************************************************************

      SUBROUTINE TRAN_SMOOTH

      INCLUDE 'COMMONS_INC.f'
 
      DIMENSION TRAN_O3_STD(5001)
      COMMON /INIT_SPECCAL16/ TRAN_O3_STD

      DIMENSION TRAN_NO2_STD(5001)
      COMMON /INIT_SPECCAL17/ TRAN_NO2_STD

      DIMENSION WAVOBS(1024),FWHM(1024)
      COMMON /GETINPUT4/ WAVOBS,FWHM

      DIMENSION TRNCAL(1024)
      COMMON /TRANCAL1/ TRNCAL,VAPTT

      COMMON /GETINPUT5/ NOBS,HSURF,DLT,DLT2

!C First stage of smoothing - smooth line-by-line high resolution spectrum with
!C     over 300,000 points (point spacing of 0.05 cm-1) to medium resolution
!C     spectrum (resolution of 0.2 nm and point spacing of 0.1 nm) with about
!C     25,000 points.
!C
!C     The smoothing of line-by-line spectrum is done in wavenumber domain. For
!C     a spectrum with a constant 0.2 nm resolution in wavelength domain, it has
!C     variable resolution in wavenumber domain. This effect is properly taken
!C     care of in the design of smoothing functions.
!C
!C     Because the high resolution spectrum is in wavenumber units (cm-1), while
!C     the medium resolution spectrum is in wavelength units, the two kinds of
!C     grids do not automatically match. In order to match the grids, arrays
!C     of INDEX_MED and TRAN_MED_INDEX are specially designed. The desired
!C     medium resolution spectrum, TRAN_MED, at constant 0.1 nm point spacing
!C     and 0.2 nm resolution is obtained through linear interpolation of
!C     TRAN_MED_INDEX array.
!C
      DO 466 J =1, NP_MED

             TRAN_MED_INDEX(J)  = 0.0
             NCVTOT_WAVNO = 2 * NCVHF_WAVNO(J) - 1
    
             SUMINS = 0.0

          DO 560 I = NCVHF_WAVNO(J), NCVTOT_WAVNO
             FINSTR_WAVNO(I) = &
                EXP( -CONST1*(FLOAT(I-NCVHF_WAVNO(J))*DWAVNO &
                /FWHM_WAVNO(J))**2.0)
             SUMINS = SUMINS + FINSTR_WAVNO(I)
 560      CONTINUE

          DO 565 I = 1, NCVHF_WAVNO(J)-1
             FINSTR_WAVNO(I) = FINSTR_WAVNO(NCVTOT_WAVNO-I+1)
             SUMINS = SUMINS + FINSTR_WAVNO(I)
 565      CONTINUE

          SUMINS = SUMINS * DWAVNO

          DO 570 I = 1, NCVTOT_WAVNO
             FINSTR_WAVNO(I) = FINSTR_WAVNO(I)*DWAVNO/SUMINS
 570      CONTINUE
  
!C*** !!!High resolution transmittances of CO2, N2O, CO, CH4, and O2 should
!C         also be calculated somewhere else (wavelength start = 0.56 micron).
!C***  Here assuming TCO2*TN2O*TCO*TCH4*TO2 is already calculated previously,
!C     i.e., TRANS(I) = TRAN_CO2(I)*TRAN_N2O(I)*TRAN_CO(I)*TRAN_CH4(I)*TRAN_O2(I)
!C      and TRAN_HI(I) = TRAN_HI_H2O(I) * TRANS(I), and TRAN_HI_H2O for varying
!C      water vapor amounts is calculated in this subroutine.

          DO 491 K = INDEX_MED(J)-(NCVHF_WAVNO(J)-1), &
                                 INDEX_MED(J)+NCVHF_WAVNO(J)-1
             TRAN_MED_INDEX(J) = TRAN_MED_INDEX(J) + TRAN_HI(K)*  &
                          FINSTR_WAVNO(K-INDEX_MED(J)+NCVHF_WAVNO(J))
 491      CONTINUE

 466  CONTINUE

!C
!C Linear interpolation to get TRAN_MED from TRAN_MED_INDEX:
!C     (Note that WAVLN_MED_INDEX(J) >= WAVLN_MED(J)    )
!C
         TRAN_MED(1)      = TRAN_MED_INDEX(1)  
         TRAN_MED(NP_MED) = TRAN_MED_INDEX(NP_MED)  

      DO J = 2, NP_MED-1
         IF(WAVLN_MED_INDEX(J).LE.WAVLN_MED(J)) THEN
           TRAN_MED(J) = TRAN_MED_INDEX(J)
         ELSE
           DLT  =  WAVLN_MED_INDEX(J) - WAVLN_MED_INDEX(J-1)
           FJM1 = (WAVLN_MED_INDEX(J) - WAVLN_MED(J))        /DLT
           FJ   = (WAVLN_MED(J)       - WAVLN_MED_INDEX(J-1))/DLT
           TRAN_MED(J) = FJM1*TRAN_MED_INDEX(J-1) + FJ*TRAN_MED_INDEX(J)
!C---
!C---           print*,j,fjm1,fj
!C---
         END IF
      END DO
!C
!C--- Here multiplying O3 and NO2 spectra and other spectrum at medium resolution:
!C
       DO I = 1, NP_STD
          TRAN_STD(I) = 1.
       END DO

       DO I = 1, NO3PT
          TRAN_STD(I) = TRAN_O3_STD(I) * TRAN_NO2_STD(I)
       END DO

       DO I = NPSHIF+1, NP_STD
          TRAN_STD(I) = TRAN_STD(I)*TRAN_MED(I-NPSHIF)
       END DO
 

!C The 2nd stage of smoothing - smooth the medium resolution spectrum (resolution
!C     of 0.2 nm and point spacing of 0.1 nm) with about 25,000 points to match
!C     the coarser and variable resolution spectrum from imaging spectrometers.
!C
!C Initialize some index parameters:
      IA = 1000
!C
      DO 1466 J =1, NOBS

             TRNCAL(J) = 0.0
             TRAN_IA     = 0.0
             TRAN_IAP1   = 0.0

             NCVTOT = 2 * NCVHF(J) - 1
!C---
!C---        print*,'J= ',j, 'NCVHF =', NCVHF(J), 'NCVTOT=',NCVTOT
    
!C Calculate instrumental response functions...
             SUMINS = 0.0

          DO 1560 I = NCVHF(J), NCVTOT
             FINSTR(I) = &
                EXP( -CONST1*(FLOAT(I-NCVHF(J))*DWAVLN &
                /FWHM(J))**2.0)
             SUMINS = SUMINS + FINSTR(I)
1560      CONTINUE

          DO 1565 I = 1, NCVHF(J)-1
             FINSTR(I) = FINSTR(NCVTOT-I+1)
             SUMINS = SUMINS + FINSTR(I)
1565      CONTINUE

          SUMINS = SUMINS * DWAVLN

          DO 1570 I = 1, NCVTOT
             FINSTR(I) = FINSTR(I)*DWAVLN/SUMINS
1570      CONTINUE
  
!C---          IF(j.eq.87) then
!C---           print*, ' J = 87 =', J
!C---           do i = 1, NCVTOT
!C---            print*,i, FINSTR(i)
!C---           end do
!C---          end if
 
!C  Index searching...
!C
          CALL HUNT(WAVLN_STD, NP_STD, WAVOBS(J), IA)
!C
!C---          print*,'J =', J, ' IA =', IA
!
!C  Smoothing...
!C
          DO 1491 K = IA-(NCVHF(J)-1), IA+NCVHF(J)-1
             TRAN_IA = TRAN_IA + TRAN_STD(K)* &
                          FINSTR(K-IA+NCVHF(J))

!C---             IF(J.eq.1) then
!C---          print*,'j =', j, 'K = ',K, ' IA= ',IA,
!C---     &       K-IA+NCVHF(J),
!C---     &       FINSTR(K-IA+NCVHF(J))
!C---             End IF

1491      CONTINUE

            IA_P1 = IA + 1
          DO 1492 K = IA_P1-(NCVHF(J)-1), IA_P1+NCVHF(J)-1
             TRAN_IAP1 = TRAN_IAP1 + TRAN_STD(K)* &
                          FINSTR(K-IA_P1+NCVHF(J))
1492      CONTINUE
!C
!C Linear interpolation to get TRNCAL from TRAN_IA and TRAN_IAP1:
!C
           DLT_IA  =  WAVLN_STD(IA_P1) - WAVLN_STD(IA)
           FIA     = (WAVLN_STD(IA_P1) - WAVOBS(J)) /DLT_IA
!C          FIA_P1  = (WAVOBS(J)     - WAVLN_STD(IA))/DLT_IA
           FIA_P1  = 1. - FIA
           TRNCAL(J) = FIA*TRAN_IA + FIA_P1*TRAN_IAP1
!C---
!C---       print*,'j=',j,'IA =',IA,'FIA =',FIA,'FIA_P1=',FIA_P1,DLT_IA
!C---
!C
1466  CONTINUE



      RETURN
      END
!********************************************************************************
!*            								       *
!*  Name: CHNLRATIO							       *
!*  Purpose: Calculate 3-channel ratios.					       *
!*  Parameters: none							       *
!*  Algorithm: The 0.94-um water vapor absorption channel is ratioed against    *
!*             the linear combination of two window channels near 0.86 and      *
!*             1.03 um to obtain one channel ratio for the 0.94-um band.	       *
!*             Similar calculation is done for the 1.14-um water vapor band.    *
!*  Globals used: NB1,NB2,NBP94,NB3,NB4,NB1P14 - number of points used in       *
!*                          channel ratios for both the .94- and 1.14-um regions*
!*                IST1,IED1,IST2,IED2,ISTP94,IEDP94 - 3-channel ratioing        *
!*                          parameters for the 0.94-um water vapor band	       *
!*                IST3,IED3,IST4,IED4,IST1P14,IED1P14 - 3-channel ratioing      *
!*                          parameters for the 1.14-um water vapor band.	       *
!*		 WT1,WT2,WT3,WT4,JA - Relative weights for the four window     *
!*                          channels used in channel-ratioing calculations. JA  *
!*			   is an output parameter from a table searching       *
!*			   routine.					       *
!*		 TRNCAL -  Atmospheric transmittance spectra.		       *
!*  Global output:R094,R114 - 3-channel ratio values for the 0.94- and 1.14-um  *
!*                          water vapor bands.				       *
!*  Return Codes: none							       *
!*  Special Considerations: none						       *
!*									       *
!********************************************************************************

      SUBROUTINE CHNLRATIO

!C  Common variables
      DIMENSION TRNCAL(1024)

      COMMON /GETINPUT7/ NB1,NB2,NBP94,NB3,NB4,NB1P14
      COMMON /INIT_SPECCAL6/ IST1,IED1,IST2,IED2,ISTP94,IEDP94
      COMMON /INIT_SPECCAL7/ IST3,IED3,IST4,IED4,IST1P14,IED1P14
      COMMON /INIT_SPECCAL8/ WT1,WT2,WT3,WT4,JA
      COMMON /TRANCAL1/ TRNCAL,VAPTT
      COMMON /CHNLRATIO1/ R094,R114

!C Calculate average of spectra over window and water vapor absorption regions.
      CONST1=0.0
      DO 560 I=IST1,IED1
        CONST1=CONST1+TRNCAL(I)
  560 CONTINUE
      CONST1=CONST1/FLOAT(NB1)

      CONST2=0.0
      DO 570 I=IST2,IED2
        CONST2=CONST2+TRNCAL(I)
  570 CONTINUE
      CONST2=CONST2/FLOAT(NB2)

      CONST3=0.0
      DO 575 I=ISTP94,IEDP94
        CONST3=CONST3+TRNCAL(I)
  575 CONTINUE
      CONST3=CONST3/FLOAT(NBP94)

      R094=CONST3/((WT1*CONST1) + (WT2*CONST2))

      CONST4=0.0
      DO 580 I=IST3,IED3 
        CONST4=CONST4+TRNCAL(I)
  580 CONTINUE
      CONST4=CONST4/FLOAT(NB3)

      CONST5=0.0
      DO 590 I=IST4,IED4
        CONST5=CONST5+TRNCAL(I)
  590 CONTINUE
      CONST5=CONST5/FLOAT(NB4)

      CONST6=0.0
      DO 595 I=IST1P14,IED1P14
        CONST6=CONST6+TRNCAL(I)
  595 CONTINUE
      CONST6=CONST6/FLOAT(NB1P14)

      R114=CONST6/((WT3*CONST4) + (WT4*CONST5))

      RETURN
      END

!********************************************************************************
!*									       *
!*  Name: PROCESS_CUBE							       *
!*  Purpose: Processes an input cube one spectral slice at a time to derive     *
!*	    surface reflectance for each spectrum and calculate the column     *
!*	    water vapor amount for each pixel.  The derived surface reflectance*
!*	    values are written to an output image file with the same dimensions*
!*	    as the input image, and the column water vapor amounts are written *
!*	    to a separate file as a single channel image.		       *
!*  Parameters: none							       *
!*  Algorithm: C programs are used to perform all of the file I/O.  A spectral  *
!*	      slice is read from the input file.  Each spectrum in the slice   *
!*	      is divided by the solar irradiance curve, then the water vapor   *
!*	      ratio is calculated.  The corresponding ratio and associated     *
!*	      transmission spectrum is used to derive the surface reflectance  *
!*	      values.							       *
!*									       *
!*	      The output is an image file in the same storage order format as  *
!*	      the input image file.  It has a 512-byte SIPS header.  The values*
!*	      are surface reflectance values multiplied by an input scale      *
!*             factor to make them integer*2 data.			       *
!*									       *
!*	      The output water vapor file is in BSQ format. It has a 512-byte  *
!*	      SIPS header.  It contains 1 channel with each pixel value        *
!*	      representing the column water vapor amount in units of cm        *
!*	      multiplied by a scale factor of 1000 to make the values          *
!*	      integer*2.						       *
!*	       								       *
!*  Globals used: NSAMPS,NLINES,NBANDS - Input image dimensions used in I/O and *
!*           cube processing.		       				       *
!*           FPIN - File pointer for input cube passed to rdslice()	       *
!*           FPOCUB - File pointer for output cube passed to wtslice()          *
!*           FPOH2O - File pointer for water vapor image passed to wtslice()    *
!*           SPATIAL,ST_SAMPLE,ST_ROW,SAMPLES_TODO,ROWS_TODO,ISLICE,DIMS - Cube *
!*             dimensions used in rdslice() and wtslice().		       *
!*           BUFFER - Holds the slice of data read by rdslice() and written by  *
!*             wtslice().						       *
!*           H2OBUF - Holds the water vapor image.			       *
!*  Global output: None.							       *
!*  Return Codes: None.							       *
!*  Special Considerations: None.					       *
!*									       *
!********************************************************************************
      SUBROUTINE PROCESS_CUBE

      use cubeio

      INTEGER              LUN_IN, LUN_OUT, LUN_VAP, I_RET
      COMMON /INOUT_UNITS/ LUN_IN, LUN_OUT, LUN_VAP

!C  Common variables
      DIMENSION YY(1024)       ! Observed radiances

!C  parameters for cube I/O
      CHARACTER (LEN = 1000) :: FINAV,FOCUB,FOH2O
      INTEGER*2 BUFFER(8388608)
      INTEGER*2 H2OBUF(8388608)
      INTEGER SORDER,HDREC

      INTEGER I_Sample, J_Line

!C Local variables
      REAL TBUF(1024)
      INTEGER OFFSET

      INTEGER J, ISAMP, IBAND
      REAL*4  SCALEF,CLMWVP

      DIMENSION YIRR(1024)
      COMMON /SOLAR_IRR1/YIRR !RJH
      DATA PI/3.1415926535897/ !RJH
!C Putting these large arrays in global memory forces memory allocation at
!C program initialization.  Otherwise, the program chokes on the first
!C call to PROCESS_CUBE.
!C---      COMMON /DUMMYGLOB/ BUFFER,H2OBUF

      COMMON /PROCCUBE1/ YY
      COMMON /GETINPUT11/HDREC,NSAMPS,NLINES,NBANDS,SORDER
      COMMON /GETINPUT12/SCALEF

!C Commons for use with the C programs for cube I/O
      COMMON /OUTCUBE/ FOCUB
      COMMON /INCUBE/ FINAV
      COMMON /OUTH2OVAP/ FOH2O
      COMMON /RJHDEBUG/ CONST1,CONST2,CONST3,CONST4,CONST5,CONST6, &
             R094CO, R114CO, JB
      COMMON /INIT_SPECCAL8/ WT1,WT2,WT3,WT4,JA
      COMMON /INIT_SPECCAL6/ IST1,IED1,IST2,IED2,ISTP94,IEDP94
       open(90,file='atrem_yy_pix512line2000_before.txt',form='formatted')
       open(91,file='atrem_yy_pix512line2000_after.txt',form='formatted')

!C Read input cube one slice at a time into the array BUFFER.
!C BUFFER contains a spectral slice of data.  The first NBANDS elements are
!C the first spectrum, the second NBANDS elements are the second spectrum, etc
!C Process each spectrum, then write to output file one slice at a time.
      DO 11 J=1,NLINES

!C Let us know the progress of program.
!        WRITE(*,*)'LINE=',J

!C-Note: BUFFER is now passed in through COMMON /DUMMYGLOB/ BUFFER,H2OBUF
!C    for reducing the required total computer memory
!C---        CALL RD_SPECTRAL_SLICE(J,HDREC,NSAMPS,NLINES,NBANDS,
!C---     &                         SORDER)
!        CALL RD_SLICE(LUN_IN,NSAMPS,NBANDS,SORDER,BUFFER)

!C For each spectrum in the slice, assign it to the YY array, and use it in
!C REFLDRV to create a surface reflectance spectrum.
        DO 46 ISAMP= 1,NSAMPS
          OFFSET = (ISAMP-1) * NBANDS
          DO 45 IBAND=1,NBANDS           ! fill YY array with this spectrum
            YY(IBAND) = BUFFER(OFFSET + IBAND)
   45     CONTINUE
    
!C Write out sample input spectrum.
        IF((ISAMP.EQ.512).AND.(J.EQ.2000)) THEN
!C            WRITE(90,*)'Islice=3, Isamp=2, Spectrum before REFLDRV'
           DO 49 I=1,NBANDS
  49         WRITE(90,*) I,ISAMP,J,PI*YY(I)/(50*YIRR(I))  !HICO CORRECTION

         ENDIF

          CALL App_Refl_With_Gas_Removal_DRV(CLMWVP)
!C         IF((ISAMP.EQ.512).AND.(J.EQ.2000)) THEN
!C           WRITE(91,*)'Islice=3, Isamp=2, Spectrum before REFLDRV'
           DO  I=1,NBANDS
            WRITE(91,'(3i5,9F9.4,4I5)') I,ISAMP,J,YY(I), &
           R094CO,R114CO,CONST1,CONST2,CONST3, &
           CONST4,CONST5,CONST6,JA,JB,IST2,IED2
           end do
!C         ENDIF
!C---          CALL REFLDRV(CLMWVP)

!C The YY array now contains true surface reflectance values. Multiply
!C by input scale factor to convert values to integer*2 format.
          DO 55 IBAND=1,NBANDS
            TBUF(IBAND)=SCALEF*YY(IBAND)
            BUFFER(OFFSET + IBAND) = NINT(TBUF(IBAND))
   55     CONTINUE

!C Save the water vapor amount (also scaled by 1000 to convert to integer*2).
!C Later, H2OBUF will be written to the water vapor output file.

          H2OBUF(ISAMP) = NINT(1000.*CLMWVP)

!C Write out sample processed spectrum.
!C--D         IF((ISAMP.EQ.2).AND.(J.EQ.3)) THEN
!C--D           WRITE(*,*)'Islice=3, Isamp=2, Spectrum after App_Refl_DRV'
!C--D           DO 48 I=1,NBANDS
!C--D             WRITE(*,*)'I=',I,' YY(I)=',YY(I)
!C--D  48         WRITE(*,*)'I=',I,' TBUF(I)=',TBUF(I)
!C--D         ENDIF

   46   CONTINUE

!C Write slice BUFFER to output file.  The output file is the same storage order
!C as the input file, and it always has a header that is one record (512 bytes).
!C-Note: BUFFER is now passed through COMMON /DUMMYGLOB/ BUFFER,H2OBUF
!C---        CALL WT_SPECTRAL_SLICE(FPOCUB,J,512,NSAMPS,NLINES,NBANDS,
!C---     &			       SORDER)
         
!        CALL WT_SLICE(LUN_OUT,NSAMPS,NBANDS,SORDER,BUFFER)
!
!        CALL WT_LINE(LUN_VAP,NSAMPS,H2OBUF)

   11 CONTINUE

!C Write the water vapor image to a file.  The file is one spatial image
!C with one header record (512 bytes).
!C-Note: H2OBUF is now passed through COMMON /DUMMYGLOB/ BUFFER,H2OBUF
!C---      CALL WT_SPATIAL_SLICE(FPOH2O,1,512,NSAMPS,NLINES,NBANDS,0)

!      CALL CLOSEINFILE(LUN_IN)
!      CALL CLOSEOUTFILE(LUN_OUT)
!      CALL CLOSEVAPFILE(LUN_VAP)

!        close(90)
!        close(91)
      RETURN
      END

!********************************************************************************
!*									       *
!*  Name: App_Refl_With_Gas_Removal_DRV 					       *
!*  Purpose: To derive column water vapor amount and apparent reflectance curve *
!*            with atmospheric gaseous absorption effect removed from an input  *
!*            spectrum using a 3-channel ratioing technique.                    *
!*  Parameters: none							       *
!*  Algorithm: the algorithm includes following procedures:		       *
!*     1. An input radiance spectrum is divided by a solar radiance curve above *
!*            the atmosphere to derive 'APPARENT REFLECTANCE SPECTRUM T(Lambda)'*
!*     2. Two 3-channel ratios:  T(0.94 um)/(WT1*T(0.86)+WT2*T(1.02))	       *
!*            and T(1.14 um)/(WT3*T(1.05)+WT4*T(1.23)), are calculated from     *
!*            the observed spectrum T(Lambda).				       *
!*     3. Based on the two ratios using look-up table procedures plus linear    *
!*            interpolation techniques, the column water vapor values and the   *
!*            atmospheric transmittance spectrum are derived.		       *
!*     4. The apparent reflectance spectrum is divided by the calculated        *
!*            transmittance spectrum to remove atmospheric gaseous absorption   *
!*            features.                                                         *
!*     5. The output apprent reflectance spectrum with atmospheric gaseous      *
!*            absorption effect is smoothed if desired.                         *
!*  Globals used: 							       *
!*  Global output:  							       *
!*  Return Codes: none							       *
!*  Special Considerations: None						       *
!* 									       *
!********************************************************************************

      SUBROUTINE App_Refl_With_Gas_Removal_DRV(CLMWVP)

!C  Common variables
      DIMENSION VAPTOT(60), R0P94(60), R1P14(60), TRNTBL(1024,60)
      DIMENSION YIRR(1024)
      DIMENSION YY(1024) 
      DIMENSION RR(1024) 
      DIMENSION ROTOT(1050), TTOT(1050), STOT(1050)
      DIMENSION FINST2(100)
      DIMENSION SSH2O(60)
      DIMENSION TRNCAL(1024)

!C Local variables
      DIMENSION TRNCV(1024)

      COMMON /GETINPUT5/ NOBS,HSURF,DLT,DLT2
      COMMON /GETINPUT7/ NB1,NB2,NBP94,NB3,NB4,NB1P14
      COMMON /GETINPUT8/ IMN,IDY,IYR,IH,IM,IS
      COMMON /GEOMETRY2/ GCO2,GO3,GN2O,GCO,GCH4,GO2,SSH2O,TOTLO3,GGEOM
      COMMON /INIT_SPECCAL3/ NH2O
      COMMON /INIT_SPECCAL6/ IST1,IED1,IST2,IED2,ISTP94,IEDP94
      COMMON /INIT_SPECCAL7/ IST3,IED3,IST4,IED4,IST1P14,IED1P14
      COMMON /INIT_SPECCAL8/ WT1,WT2,WT3,WT4,JA
      COMMON /INIT_SPECCAL10/ NCV2,NCVHF2,NCVTT2,ISTRT2,IEND2,FINST2
      COMMON /INIT_SPECCAL11/ NATOT,NBTOT,NCTOT,NDTOT
      COMMON /TRAN_TABLE1/ SH2O,VAPTOT,R0P94,R1P14,TRNTBL
      COMMON /TRANCAL1/ TRNCAL,VAPTT
      COMMON /SOLAR_IRR1/YIRR
      COMMON /PROCCUBE1/ YY
      COMMON /SIXS1/ ROTOT, TTOT, STOT
      COMMON /RJHDEBUG/ CONST1,CONST2,CONST3,CONST4,CONST5,CONST6, &
             R094CO, R114CO, JB
      DIMENSION G_VAP(25), G_OTHER(25)
      COMMON /GEOMETRY3/ G_VAP, G_OTHER, G_VAP_EQUIV, VAP_SLANT_MDL

!C Parameters for plane observations:
      CHARACTER (LEN = 80) :: NAME_INSTRU, NAMES(10)
      COMMON /GETINPUT13/ NAME_INSTRU, NAMES

!C Arrays related to look-up table:
!C     VAPTOT: TOTAL SUN-SURFACE-SENSOR PATH WATER VAPOR IN UNITS OF CM
!C     R0P94 : Tc(0.94 um)/(WT1*Tc(0.86)+WT2*Tc(1.02))
!C     R1P14 : Tc(1.14 um)/(WT3*Tc(1.05)+WT4*Tc(1.23))

      DIMENSION SPECA(1024),SPECB(1024),SPECAV(1024)

	integer*4 itmp1
	real*4 XFA
	real*4 XFB
	real*4 XFC
	real*4 XFD

      DATA PI,DTORAD /3.1415926,0.0174533/

!C
!C Loops for processing image data cube ...
      IF((NAME_INSTRU.EQ.NAMES(1)).AND.(IYR.LE.2009)) THEN

!C  Ratioing observed spectrum against solar irradiance curve.
!C  The radiance data from each year is scaled by a different constant.
!C    XFA = factor for A spectrometer
!C    XFB = factor for B spectrometer
!C    XFC = factor for C spectrometer
!C    XFD = factor for D spectrometer
!C
!C rewritten to handle each spectrometer differently
!C                                         9/14/1995 Roger N. Clark, U.S.G.S.

      IF(IYR.LE.1989) THEN
	XFA = 10.0
	XFB = 10.0
	XFC = 10.0
	XFD = 10.0
      ENDIF

      IF((IYR.GE.1990).AND.(IYR.LE.1991)) THEN
	XFA = 20.0
	XFB = 20.0
	XFC = 20.0
	XFD = 20.0
      ENDIF

      IF((IYR.GE.1992).AND.(IYR.LE.1994)) THEN
	XFA = 50.0
	XFB = 50.0
	XFC = 50.0
	XFD = 50.0
      ENDIF

      IF((IYR.GE.1995).AND.(IYR.LE.2009)) THEN
	XFA = 50.0
	XFB = 50.0
	XFC = 50.0
	XFD = 100.0
      ENDIF

!C   AVIRIS A Spectrometer (ch = 1-32)

      IF(NOBS.GE.32) THEN
	itmp1 = 32
      ELSE
	itmp1 = NOBS
      ENDIF

      DO 551 I=1,itmp1
          YY(I)=PI*YY(I)/(XFA*YIRR(I))
  551 CONTINUE

!C   AVIRIS B Spectrometer (ch = 33-96)

      IF(NOBS.GE.96) THEN
	itmp1 = 96
      ELSE if ((NOBS.GE.33).and.(NOBS.LT.96)) THEN
	itmp1 = NOBS
      ELSE
	go to 559
      ENDIF

      DO 552 I=33,itmp1
          YY(I)=PI*YY(I)/(XFB*YIRR(I))
  552 CONTINUE

!C   AVIRIS C Spectrometer (ch = 97-160)

      IF(NOBS.GE.160) THEN
	itmp1 = 160
      ELSE if ((NOBS.GE.97).and.(NOBS.LT.160)) THEN
	itmp1 = NOBS
      ELSE
	go to 559
      ENDIF

      DO 553 I=97,itmp1
          YY(I)=PI*YY(I)/(XFC*YIRR(I))
  553 CONTINUE

!C   AVIRIS D Spectrometer (ch = 161-224)

      IF(NOBS.GE.161) THEN
	itmp1 = NOBS
      ELSE
	go to 559
      ENDIF

      DO 554 I=161,itmp1
          YY(I)=PI*YY(I)/(XFD*YIRR(I))
  554 CONTINUE

!C branch point in case not all spectrometers are being done.
  559   CONTINUE

!C End of R. Clark contribution
      END IF
!C
!C Loops for processing AVIRIS data starting from year 2010---
      IF((NAME_INSTRU.EQ.NAMES(1)).AND.(IYR.GE.2010)) THEN

          F_AV_1 = 30.
          F_AV_2 = 60.
          F_AV_3 = 120.

        DO I = 1, 110
           YY(I) = PI*YY(I) / (F_AV_1 * YIRR(I))
        END DO

        DO I = 111, 160
           YY(I) = PI*YY(I) / (F_AV_2 * YIRR(I))
        END DO

        DO I = 161, 224
           YY(I) = PI*YY(I) / (F_AV_3 * YIRR(I))
        END DO

      END IF
 
!C Loops for processing HYDICE data ...
!C
      IF(NAME_INSTRU.EQ.NAMES(2)) THEN

          F_HYDICE = 75.

        DO I = 1, NOBS
           YY(I) = PI*YY(I) / (F_HYDICE * YIRR(I))
        END DO

      END IF 

 
!C Loops for processing HSI data ...
!C
      IF(NAME_INSTRU.EQ.NAMES(3)) THEN

          F_HSI = 100.

        DO I = 1, NOBS
           YY(I) = PI*YY(I) / (F_HSI * YIRR(I))
        END DO

      END IF 
 

!C Loops for processing TRWIS-III data ...
!C
      IF(NAME_INSTRU.EQ.NAMES(4)) THEN

          F_TRWIS = 100.

        DO I = 1, NOBS
           YY(I) = PI*YY(I) / (F_TRWIS * YIRR(I))
        END DO

      END IF 

!C Loops for processing EO-1 Hyperion data ...
!C
      IF(NAME_INSTRU.EQ.NAMES(6)) THEN

          F_Hyperion_A = 40.
          F_Hyperion_B = 80.

        DO I = 1, 70
           YY(I) = PI*YY(I) / (F_Hyperion_A * YIRR(I))
        END DO

        DO I = 71, NOBS 
           YY(I) = PI*YY(I) / (F_Hyperion_B * YIRR(I))
        END DO

      END IF

!C Loops for processing HICO data ...
!C
      IF(NAME_INSTRU.EQ.NAMES(7)) THEN

!C---          F_HICO = 50./1.32
          F_HICO = 50.

        DO I = 1, NOBS
           YY(I) = PI*YY(I) / (F_HICO * YIRR(I))
        END DO

      END IF

!C Loops for processing NIS (NEON Imaging Spectrometer) data ...
!C
      IF(NAME_INSTRU.EQ.NAMES(8)) THEN

          F_NIS = 100.

        DO I = 1, NOBS
           YY(I) = PI*YY(I) / (F_NIS * YIRR(I))
        END DO

      END IF 

!C Loops for processing PRISM Imaging Spectrometer data ...
!C
      IF(NAME_INSTRU.EQ.NAMES(9)) THEN

          F_PRISM = 100.

        DO I = 1, NOBS
           YY(I) = PI*YY(I) / (F_PRISM * YIRR(I))
        END DO

      END IF
 
!C
!C Calculating 3-channel ratios from an observed spectrum, using a
!C  look-up table procedure to derive the column amount of water vapor
!C  and to find the associated simulated transmittance spectrum.
      CONST1=0.0
      DO 560 I=IST1,IED1
        CONST1=CONST1+YY(I)
  560 CONTINUE
      CONST1=CONST1/FLOAT(NB1)

      CONST2=0.0
      DO 570 I=IST2,IED2
        CONST2=CONST2+YY(I)
  570 CONTINUE
      CONST2=CONST2/FLOAT(NB2)

      CONST3=0.0
      DO 575 I=ISTP94,IEDP94
        CONST3=CONST3+YY(I)
  575 CONTINUE
      CONST3=CONST3/FLOAT(NBP94)

      R094CO=CONST3/((WT1*CONST1) + (WT2*CONST2))
      R094C =R094CO

      IF(R094CO.GT.1.0) THEN
        CONST1=0.0
        DO 561 I=IST1,IED1
          CONST1=CONST1+1.0/YY(I)
  561   CONTINUE
        CONST1=CONST1/FLOAT(NB1)

        CONST2=0.0
        DO 571 I=IST2,IED2
          CONST2=CONST2+1.0/YY(I)
  571   CONTINUE
        CONST2=CONST2/FLOAT(NB2)

        CONST3=0.0
        DO 576 I=ISTP94,IEDP94
          CONST3=CONST3+1.0/YY(I)
  576   CONTINUE
        CONST3=CONST3/FLOAT(NBP94)

        R094C=CONST3/((WT1*CONST1) + (WT2*CONST2))
      ENDIF

      CONST4=0.0
      DO 580 I=IST3,IED3
        CONST4=CONST4+YY(I)
  580 CONTINUE
      CONST4=CONST4/FLOAT(NB3)

      CONST5=0.0
      DO 590 I=IST4,IED4
        CONST5=CONST5+YY(I)
  590 CONTINUE
      CONST5=CONST5/FLOAT(NB4)

      CONST6=0.0
      DO 595 I=IST1P14,IED1P14
        CONST6=CONST6+YY(I)
  595 CONTINUE
      CONST6=CONST6/FLOAT(NB1P14)

      R114CO=CONST6/((WT3*CONST4) + (WT4*CONST5))
      R114C =R114CO

      IF(R114CO.GT.1.0) THEN
        CONST4=0.0
        DO 581 I=IST3,IED3
          CONST4=CONST4+1.0/YY(I)
  581   CONTINUE
        CONST4=CONST4/FLOAT(NB3)

        CONST5=0.0
        DO 591 I=IST4,IED4
          CONST5=CONST5+1.0/YY(I)
  591   CONTINUE
        CONST5=CONST5/FLOAT(NB4)

        CONST6=0.0
        DO 596 I=IST1P14,IED1P14
          CONST6=CONST6+1.0/YY(I)
  596   CONTINUE
        CONST6=CONST6/FLOAT(NB1P14)

        R114C=CONST6/((WT3*CONST4) + (WT4*CONST5))
      ENDIF

      CALL HUNT(R0P94,NH2O,R094C,JA)
      IF (JA.GT.0.AND.JA.LT.60) THEN
        DLTA  =R0P94(JA+1)-R0P94(JA)
        FJA   =(R0P94(JA+1)-R094C)/DLTA
        FJAP1 =(R094C-R0P94(JA))/DLTA

        VAPTTA=FJA*VAPTOT(JA)+FJAP1*VAPTOT(JA+1)
        IF(R094CO.GT.1.) VAPTTA=-VAPTTA
      ELSE
        IF(JA.LE.0) VAPTTA = VAPTOT(JA+1)
        IF(JA.GE.60) VAPTTA = VAPTOT(JA)
      ENDIF

      IF(R094CO.LE.1.)THEN
        DO 900 I=1,NOBS
          IF (JA.GT.0.AND.JA.LT.60) THEN
            SPECA(I)=FJA*TRNTBL(I,JA)+FJAP1*TRNTBL(I,JA+1)
          ELSE
            IF(JA.LE.0) SPECA(I)=TRNTBL(I,JA+1)
            IF(JA.GE.60) SPECA(I)=TRNTBL(I,JA)
          ENDIF
  900   CONTINUE
      ENDIF

      IF(R094CO.GT.1.)THEN
        DO 901 I=1,NOBS
          IF (JA.GT.0.AND.JA.LT.60) THEN
            SPECA(I)=1.0/(FJA*TRNTBL(I,JA)+FJAP1*TRNTBL(I,JA+1))
          ELSE
            IF(JA.LE.0) SPECA(I)=1.0/TRNTBL(I,JA+1)
            IF(JA.GE.60) SPECA(I)=1.0/TRNTBL(I,JA)
          ENDIF
  901   CONTINUE
      ENDIF

      JB = JA
      CALL HUNT(R1P14,NH2O,R114C,JB)
      IF (JB.GT.0.AND.JB.LT.60) THEN
        DLTB  =R1P14(JB+1)-R1P14(JB)
        FJB   =(R1P14(JB+1)-R114C)/DLTB
        FJBP1 =(R114C-R1P14(JB))/DLTB
  
        VAPTTB=FJB*VAPTOT(JB)+FJBP1*VAPTOT(JB+1)
        IF(R114CO.GT.1.) VAPTTB=-VAPTTB
      ELSE
        IF(JB.LE.0) VAPTTB = VAPTOT(JB+1)
        IF(JB.GE.60) VAPTTB = VAPTOT(JB)
      ENDIF

      IF(R114CO.LE.1.)THEN
        DO 910 I=1,NOBS
          IF (JB.GT.0.AND.JB.LT.60) THEN
            SPECB(I)=FJB*TRNTBL(I,JB)+FJBP1*TRNTBL(I,JB+1)
          ELSE
            IF(JB.LE.0) SPECB(I)=TRNTBL(I,JB+1)
            IF(JB.GE.60) SPECB(I)=TRNTBL(I,JB)
          ENDIF
  910   CONTINUE
      ENDIF

      IF(R114CO.GT.1.)THEN
        DO 911 I=1,NOBS
          IF (JB.GT.0.AND.JB.LT.60) THEN
            SPECB(I)=1.0/(FJB*TRNTBL(I,JB)+FJBP1*TRNTBL(I,JB+1))
          ELSE
            IF(JB.LE.0) SPECB(I)=1.0/TRNTBL(I,JB+1)
            IF(JB.GE.60) SPECB(I)=1.0/TRNTBL(I,JB)
          ENDIF
  911   CONTINUE
      ENDIF

!C---      CLMWVP=0.5*(VAPTTA+VAPTTB)/GGEOM
      CLMWVP=0.5*(VAPTTA+VAPTTB)/G_VAP_EQUIV
!C
!C  Derivation of surface reflectances
      DO 915 I=1,NOBS
        SPECAV(I)=0.5*(SPECA(I)+SPECB(I))
!C---        TRNTMP=(YY(I)/SPECAV(I))-ROTOT(I)
!C---        YY(I)=TRNTMP/(TTOT(I)+(STOT(I)*TRNTMP))
        YY(I) = YY(I) / SPECAV(I)
  915 CONTINUE
!C
!C--       print*, 'i, rotot(i), ttot(i), stot(i)'
!C--      DO I = 1, NOBS
!C--       print*, i, rotot(i), ttot(i), stot(i)
!C--      END DO

!C
!C Smooth the derived surface reflectance spectra

      IF(DLT2.GT.DLT) THEN
!C
!C       First, replace radiances <= 0 near the spectral overlapping parts of the
!C       four AVIRIS spectrometers by radiances in the nearby AVIRIS' channels.
        DO 920 I=NATOT-2,NATOT+2
          IF(YY(I).LE.0.0) YY(I)=YY(I-1)
  920   CONTINUE

        DO 921 I=NBTOT-2,NBTOT+2
          IF(YY(I).LE.0.0) YY(I)=YY(I-1)
  921   CONTINUE

        DO 922 I=NCTOT-4,NCTOT+4
          IF(YY(I).LE.0.0) YY(I)=YY(I-1)
  922   CONTINUE

        DO 923 I=NDTOT-7,NOBS
          IF(YY(I).LE.0.0) YY(I)=YY(I-1)
  923   CONTINUE

        DO 3480 I=ISTRT2,IEND2
          TRNCV(I)=0.0
          II=I-NCV2
          DO 3490 J=I-NCV2,I+NCV2
            TRNCV(I)=TRNCV(I)+YY(J)*FINST2(J-II+1)
 3490     CONTINUE
 3480   CONTINUE
        DO 3495 I=ISTRT2,IEND2
          YY(I)=TRNCV(I)
 3495   CONTINUE
      ENDIF

      RETURN
      END

!********************************************************************************
!*            								       *
!*  Name: LOCATE								       *
!*  Purpose: given an array XX of length N, and given a value X, returns a value*
!*           J such that X is between XX(J) and XX(J+1).  XX must be monotonic, *
!*  Parameters: XX - monotonic array of values				       *
!*              N  - number of elements in XX				       *
!*              X  - value that will be matched to the XX array 		       *
!*              J  - index into the XX array where XX(J) <= X <= XX(J+1)	       *
!*  Algorithm:  bisectional table searching, copied from Numerical Recipes.     *
!*  Globals used: none							       *
!*  Global output: none							       *
!*  Return Codes: J=0 or J=N is returned to indicate that X is out of range     *
!*  Special Considerations: none						       *
!*									       *
!********************************************************************************
      SUBROUTINE locate(xx,n,x,j)
      INTEGER j,n
      REAL x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      if(x.eq.xx(1))then
        j=1
      else if(x.eq.xx(n))then
        j=n-1
      else
        j=jl
      endif
      return
      END

!********************************************************************************
!*            								       *
!*  Name: CUBSPLN 							       *
!*  Purpose: an interface for performing cubic spline interpolations.	       *
!*  Parameters: N - number of elements in XORGN, YORGN			       *
!*              XORGN - original x values				       *
!*              YORGN - original y values				       *
!*              XINT  - interpolated x values				       *
!*              YINT  - interpolated y values				       *
!*  Algorithm: Straight forward calculations				       *
!*  Globals used: NOBS - number of spectral points.			       *
!*  Global output: none 							       *
!*  Return Codes: none							       *
!*  Special Considerations: none						       *
!*									       *
!********************************************************************************

      SUBROUTINE CUBSPLN(N,XORGN,YORGN,XINT,YINT)

      DIMENSION XORGN(1050),YORGN(1050),Y2(1050)
      DIMENSION XINT(1024),YINT(1024)
      INTEGER N                               !number of elements in XORGN,YORGN
      
      COMMON /GETINPUT5/ NOBS,HSURF,DLT,DLT2
!C
!C  YP1 and YP2 are two parameters specifying the values of 2nd derivatives at
!C     XORGN(1) and XORGN(N) and were used in cubic spline interpolations. The
!C     setting of both YP1 and YP2 to 2.0E30 is equivalent to set the 2nd
!C     derivatives at the two boundaries as zero, according to Numeric Recipes.
      YP1=2.0E30
      YPN=2.0E30

      CALL SPLINE(XORGN,YORGN,N,YP1,YPN,Y2)

      DO 320 I=1,NOBS
        X=XINT(I)
        CALL SPLINT(XORGN,YORGN,Y2,N,X,Y)
        IF(Y.LT.0.0) Y=0.0
        YINT(I)=Y

  320 CONTINUE

      RETURN
      END
!
!********************************************************************************
!*     									       *
!*  Name: SPLINE								       *
!*  Purpose: program for cubic spline interpolation --- copyed from Numerical   *
!*           Recipes, p.88-89 						       *
!*  Parameters: X - x-values						       *
!*              Y - y-values						       *
!*              N - length of X						       *
!*              YP1 - 2nd derivative at X(1)				       *
!*              YPN - 2nd derivative at X(N)				       *
!*              Y2  - 2nd derivatives at all X positions			       *
!*  Algorithm: Given arrays X and Y of length N containing a tabulated function,*
!*      i.e.,Yj=f(Xj), with X1 < X2...<XN, and given values YP1 and YPN for     *
!*      the first derivative of the interpolating function at points 1 	       *
!*      and N, respectively, this routine returns an array Y2 of length N       *
!*      which contains the second derivatives of the interpolating function     *
!*      at the tabulated points Xj. If YP1 and/or YPN are equal to 1.0E30 or    *
!*      larger, the routine is signalled to set the corresponding boundary      *
!*      condition for a natural spline, with zero second derivative on          *
!*      that boundary.							       *
!*  Globals used: none							       *
!*  Global output: none							       *
!*  Return Codes: none							       *
!*  Special Considerations: none						       *
!*									       *
!********************************************************************************

      subroutine spline(x,y,n,yp1,ypn,y2)
      parameter (nmax=1050)
      integer n,i,k
      real x(n),y(n),y2(n),u(nmax)
      real yp1,ypn,sig,p,qn,un
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
           /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      end

!********************************************************************************
!*            								       *
!*  Name: SPLINT								       *
!*  Purpose: calculate a cubic-spline interpolated Y value, -- copied from      *
!*           Numerical Recipes.						       *
!*  Parameters: XA  - original X array					       *
!*              YA  - original Y array					       *
!*              Y2A - 2nd derivative array from SUBROUTINE SPLINE	       *
!*              N   - number of X elements				       *
!*              X   - an X position at which interpolation should be made       *
!*              Y   - interpolated y value at position X  		       *
!*  Algorithm: Given the arrays XA and YA of length N, which tabulate a function*
!*      (with the XAj's in order), and given the array Y2A, which is the output *
!*      from SPLINE above, and given a value of X, this routine returns a       *
!*      cubic-spline interpolated value Y.				       *
!*  Globals used: none 							       *
!*  Global output: none							       *
!*  Return Codes: none							       *
!*  Special Considerations: SPLINE is called only once to process an entire     *
!*      tabulated function in arrays Xi and Yi. Once this has been done, values *
!*      of the interpolated function for any value of X are obtained by calls   *
!*      (as many as desired) to a separate routine SPLINT (for cubic spline     *
!*      interpolation")							       *
!*									       *
!********************************************************************************

      subroutine splint(xa,ya,y2a,n,x,y)
      integer n,klo,khi,k
      real xa(n),ya(n),y2a(n)
      real x,y,h,a,b
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) stop 'bad xa input.'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+ &
           ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end

!********************************************************************************
!*            								       *
!*  Name: HUNT								       *
!*  Purpose: finds the element in array XX that is closest to value X.  Array AA*
!*	    must be monotonic, either increasing or decreasing.		       *
!*  Parameters: XX  - array to search					       *
!*              N - number of elements in the array			       *
!*              X - element to search for closest match			       *
!*	       JLO - index of the closest matching element		       *
!*  Algorithm: this subroutine was copied from Numerical Recipes		       *
!*  Globals used: none 							       *
!*  Global output: none							       *
!*  Return Codes: none							       *
!*  Special Considerations: none						       *
!*									       *
!********************************************************************************

      SUBROUTINE hunt(xx,n,x,jlo)
      INTEGER jlo,n
      REAL x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).ge.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)then
        if(x.eq.xx(n))jlo=n-1
        if(x.eq.xx(1))jlo=1
        return
      endif
      jm=(jhi+jlo)/2
      if(x.ge.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END

!********************************************************************************
!*            								       *
!*  Name: SUNCOR							      	       *
!*  Purpose: computes at some reference time TZ the declination of the sun and  *
!*           the hour angle-TZ.						       *
!*  Algorithm: At any time T, the local solar hour angle is HAZ+T(GMT)-XLONG,   *
!*           where XLONG is the longitude measured positive west of Greenwich.  *
!*           Both time and angles are in radians. To compute azimuth and        *
!*           elevation at latitude XLAT, longitude XLONG, and time T, call      *
!*           HAZEL(HAZ+T-XLONG,DEC,AZ,EL,XLAT). This routine was copied	       *
!*           from W. Mankin at NCAR, Boulder, CO.			       *
!*  Globals used: none 							       *
!*  Global output: none 							       *
!*  Return Codes: none							       *
!*  Special Considerations: none						       *
!*									       *
!********************************************************************************

      SUBROUTINE SUNCOR(IDAY,MONTH,IYR,TZ,DEC,HAZ)

      JD=JULIAN(IDAY,MONTH,IYR)
      FJD=0.5+TZ/6.283185307
      CALL SOLCOR(JD,FJD,RAS,DEC,GSDT,BZERO,PZERO,SOLONG)
      HAZ=GSDT-RAS-TZ

      RETURN
      END

!********************************************************************************
!*  Name:  SOLCOR							       *
!*  Purpose: A routine for solar geometry calculations --- copied from	       *
!*                  W. Mankin at NCAR.					       *
!*  Parameters: none							       *
!*  Algorithm:  This routine was supplied by W. Mankin at NCAR, Boulder, CO     *
!*  Globals used: 							       *
!*  Global output:  							       *
!*  Return Codes: none							       *
!*  Special Considerations: none						       *
!*									       *
!********************************************************************************
      SUBROUTINE SOLCOR(JD,FJD,RAS,DECS,GSDT,BZRO,P,SOLONG)

      PI=3.141592654
      D=(JD-2415020)+FJD
      IYR=D/365.25
      G=-.026601523+.01720196977*D-1.95E-15*D*D -2.*PI*IYR
      XLMS=4.881627938+.017202791266*D+3.95E-15*D*D-2.*PI*IYR
      OBL=.409319747-6.2179E-9*D
      ECC=.01675104-1.1444E-9*D
      F=D-365.25*IYR
      GSDT=1.739935476+2.*PI*F/365.25+1.342027E-4*D/365.25
      GSDT=GSDT+6.2831853*(FJD-0.5)
      XLTS=XLMS+2.*ECC*SIN(G)+1.25*ECC*ECC*SIN(2.*G)
      SNDC=SIN(XLTS)*SIN(OBL)
      DECS=ASIN(SNDC)
      CSRA=COS(XLTS)/COS(DECS)
      RAS=ACOS(CSRA)
      IF(SIN(XLTS).LT.0.) RAS=2.*PI-RAS
      OMEGA=1.297906+6.66992E-7*D
      THETAC=XLTS-OMEGA
      BZRO=ASIN(.126199*SIN(THETAC))
      P=-ATAN(COS(XLTS)*TAN(OBL))-ATAN(.127216*COS(THETAC))
      XLMM=ATAN2(.992005*SIN(THETAC),COS(THETAC))
      JDR=JD-2398220
      IROT=(JDR+FJD)/25.38
      FROT=(JDR+FJD)/25.38-IROT
      SOLONG=XLMM-2.*PI*FROT+PI-3.E-4
      IF(SOLONG.LT.0.) SOLONG=SOLONG+2.*PI
      IF(SOLONG.GT.(2.*PI)) SOLONG=SOLONG-2.*PI

      RETURN
      END


!********************************************************************************
!*     									       *
!*  Name: JULIAN								       *
!*  Purpose: computes the julian day 					       *
!*  Parameters: IDAY - day of the year					       *
!*              MONTH - day of the year					       *
!*  Algorithm: This routine was supplied by William Mankin at NCAR, Boulder, CO *
!*  Globals used: none							       *
!*  Global output:  none							       *
!*  Return Codes: none							       *
!*  Special Considerations: none						       *
!*									       *
!********************************************************************************

      FUNCTION JULIAN(IDAY,MONTH,IYR)

!C Local variables
      DIMENSION MD(12)

      DATA MD/0,31,59,90,120,151,181,212,243,273,304,334/

      IF(IYR.LT.100) IYR=IYR+1900
      JYR=IYR-1600
      I1=JYR/400
      I2=(JYR-400*I1)/100
      I3=(JYR-400*I1-100*I2)/4
      JULIAN=2305447+365*JYR+97*I1+24*I2+I3
      LEAP=0
      IF(MOD(JYR,4).EQ.0) LEAP=1
      IF(MOD(JYR,100).EQ.0) LEAP=0
      IF(MOD(JYR,400).EQ.0) LEAP=1

!C     LEAP=1 if iyr is a leap year
      JULIAN=JULIAN+MD(MONTH)+IDAY
      IF(MONTH.LE.2) JULIAN=JULIAN-LEAP

      RETURN
      END

!********************************************************************************
!*            								       *
!*  Name: HAZEL 								       *
!*  Purpose: Calculates azimuth and elevation				       *
!*  Parameters: none							       *
!*  Algorithm: This routine was supplied by William Mankin at NCAR, Boulder, CO *
!*  Globals used: None. 							       *
!*  Global output:   None.						       *
!*  Return Codes: None							       *
!*  Special Considerations: None.					       *
!*									       *
!********************************************************************************
      SUBROUTINE HAZEL (H,D,A,E,XLAT)

!C     H = HOUR ANGLE
!C     D = DECLINATION
!C     A = AZIMUTH
!C     E = ELEVATION
!C     XLAT = LATITUDE
!C     ALL ANGLES IN RADIANS

      PI = 3.14159265
      SNE = SIN(D)*SIN(XLAT)+COS(D)*COS(XLAT)*COS(H)
      E=ASIN(SNE)
      SNA = COS(D)*SIN(H)
      CSA=(SIN(XLAT)*COS(H)*COS(D)-SIN(D)*COS(XLAT))
      A=ATAN2(SNA,CSA)+PI
      END

!********************************************************************************
!*            								       *
!*  Name: FINDMATCH							       *
!*  Purpose: finds the closest match for ELEM in LIST			       *
!*  Parameters:  LIST - array of values to match ELEM to.  Elements is array    *
!*	          should increase in value.				       *
!*  Algorithm: linearly compare ELEM to each element.  When ELEM is smaller     *
!*             than the LIST(I), then it is assumed to be closest to LIST(I-1)  *
!*             or LIST(I).  The one that has the smallest absolute difference   *
!*             to ELEM is returned as the closest match.			       *
!*  Globals used: none							       *
!*  Global output: none							       *
!*  Return Codes: the closest matching element index			       *
!*  Special Considerations: none						       *
!*									       *
!********************************************************************************

      FUNCTION FINDMATCH(LIST,NOBS,ELEM)

      DIMENSION LIST(1024)
      INTEGER   NOBS
      REAL      ELEM,LIST

      DO 460 I=1,NOBS
        IF(LIST(I).GT.ELEM) GOTO 470
  460 CONTINUE
  470 CONTINUE
      DIFF1=ABS(LIST(I-1)-ELEM)
      DIFF2=ABS(LIST(I)-ELEM)
      IF (DIFF1.LT.DIFF2) THEN
        FINDMATCH=I-1
      ELSE
        FINDMATCH=I
      ENDIF

      RETURN
      END
!C--------1---------2---------3---------4---------5---------6---------7--
