!------------------------------------------------------------------------------|
! 2015 - Modifications for l2gen including speed enhancements by R. Healy      |
!        (richard.healy@nasa.gov)                                              |
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

!INTERFACE
!  SUBROUTINE ecdf(xcdf, ycdf, nbin, xs, nsamp) BIND(C)
!    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_FLOAT
!    IMPLICIT NONE
!    FLOAT(C_FLOAT) :: xcdf,ycdf,xs
!    INTEGER(C_INT) :: nsamp,nbin
!  END SUBROUTINE addnums
!END INTERFACE
      
      SUBROUTINE GET_INPUT

!      use cubeio

      INCLUDE 'COMMONS_INC_v2.f'

      INTEGER              LUN_IN, LUN_OUT, LUN_VAP, I_RET
      COMMON /INOUT_UNITS/ LUN_IN, LUN_OUT, LUN_VAP


!C  Common variables
      DIMENSION H(25), T(25), P(25), VMR(25)
      DIMENSION WAVOBS(NOBS_MAX),FWHM(NOBS_MAX)
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
      COMMON /GETINPUT5/ NOBS,IFULLCALC,HSURF,DLT,DLT2
      COMMON /GETINPUT6/ WNDOW1,WNDOW2,WP94C,WNDOW3,WNDOW4,W1P14C
      COMMON /GETINPUT7/ NB1,NB2,NBP94,NB3,NB4,NB1P14
      COMMON /GETINPUT8/ IMN,IDY,IYR,IH,IM,IS
      COMMON /GETINPUT9/ XLATD,XLATM,XLATS,LATHEM
      COMMON /GETINPUT10/XLONGD,XLONGM,XLONGS,LNGHEM
      COMMON /GETINPUT11/HDREC,NSAMPS,NLINES,NBANDS,SORDER
      COMMON /GETINPUT12/SCALEF
      COMMON /TPVMR_INIT1/ TPVMR
      COMMON /MODEL_ADJ1/ CLMVAP,Q

      COMMON /GEOMETRY1/ SOLZNI,SOLAZ,OBSZNI,OBSPHI,IDAY
      REAL   SSH2O(NH2O_MAX)
      COMMON /GEOMETRY2/ GCO2,GO3,GN2O,GCO,GCH4,GO2,SSH2O,TOTLO3,GGEOM
      REAL   G_VAP(25), G_OTHER(25)
      COMMON /GEOMETRY3/ G_VAP, G_OTHER, G_VAP_EQUIV
      COMMON /GEOMETRY4/VAP_SLANT_MDL
      REAL   MU,MU0,SSH2O_S(NH2O_MAX,2),VMRM_S(25,2)
      COMMON /GEOMETRY5/MU,MU0,SSH2O_S

      DIMENSION HP(25), TP(25), PP(25), VMRP(25)
      COMMON /MODEL_ADJ2/ HP, TP, PP, VMRP
      COMMON /MODEL_ADJ3/ K_PLANE, DVAP_PLANE, DVAP_LAYER, &
                         DP_PLANE, DP_LAYER, CLMVAPP
      COMMON /MODEL_ADJ4/ K_SURF

!C Commons for use with the C programs for cube I/O
      COMMON /OUTCUBE/ FOCUB
      COMMON /INCUBE/ FINAV
      COMMON /OUTH2OVAP/ FOH2O

!C Parameters for names of imaging spectrometers
      CHARACTER (LEN = 80) :: NAME_INSTRU, NAMES(10)
      COMMON /GETINPUT13/ NAME_INSTRU, NAMES

      DIMENSION TRAN_O3_STD(NO3PT)
      COMMON /INIT_SPECCAL16/ TRAN_O3_STD
      DIMENSION TRAN_NO2_STD(NO3PT)
      COMMON /INIT_SPECCAL17/ TRAN_NO2_STD
      DIMENSION O3CF(NO3PT)
      COMMON /O3CF_INIT1/ O3CF
      DIMENSION RNO2CF(NO3PT)
      COMMON /NO2CF_INIT1/ RNO2CF
      DIMENSION DP(25), PM(25), TM(25), VMRM(25)
      COMMON /INIT_SPECCAL5/ DP,PM,TM,VMRM

      REAL XPSS, XPPP
      COMMON /GETINPUT14/ XPSS, XPPP

      REAL XVIEWD, XVIEWM, XVIEWS
      REAL XAZMUD, XAZMUM, XAZMUS
      COMMON /GETINPUT15/ XVIEWD,XVIEWM,XVIEWS, XAZMUD,XAZMUM,XAZMUS

      GOOD_DATA = .TRUE.    ! initialize flag for good data

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

  85  FORMAT(A1000)

!C        ELSE

!C Initialize default atmospheric window and water vapor absorption regions
!C     for 3-channel ratioing calculations.
!C          WNDOW1 = 0.865
!C          NB1    = 3
!C          WNDOW2 = 1.030
!C          NB2    = 3
!C          WP94C  = 0.940
!C          NBP94  = 5
!C          WNDOW3 = 1.050
!C          NB3    = 3
!C          WNDOW4 = 1.235
!C          NB4    = 3
!C          W1P14C = 1.1375
!C          NB1P14 = 7
!C        ENDIF
!C      ENDIF

!C  Determine which atmospheric model to use
!C     1 = tropical
!C     2 = mid latitude summer
!C     3 = mid latitude winter
!C     4 = subarctic summer
!C     5 = subarctic winter
!C     6 = US standard 1962
!C     7 = user defined model
      MODEL = 1
!C       Initialize NB, H, P, T, VMR from predefined atmospheric model array
        NB = TPVMR(MODEL,1)
        DO 120 I=1,NB
          H(I) = TPVMR(MODEL,2+(4*(I-1)))
          !Convert the atmospheric pressure from "mb" to "atm.":
          P(I) = TPVMR(MODEL,3+(4*(I-1))) / 1013.
          T(I) = TPVMR(MODEL,4+(4*(I-1)))
          !Convert the VMR from the ppm unit in the model to absolute unit
          VMR(I) = TPVMR(MODEL,5+(4*(I-1)))*1.0E-06
  120   CONTINUE
      NL=NB-1

        DO I     = NB+1, 25
          H(I)   = 1000.
          P(I)   = 0.0
          T(I)   = 300.
          VMR(I) = 0.0
        END DO


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

!C      READ(*,*)IH2OVP, ICO2, IO3, IN2O, ICO, ICH4, IO2, INO2

    IH2OVP = 1
    ICO2   = 0
    IO3    = 0
    IN2O   = 0
    ICO    = 0
    ICH4   = 0
    IO2    = 0
    INO2   = 0

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

      INCLUDE 'COMMONS_INC_v2.f'

!C  Common variables
      DIMENSION H(25), T(25), P(25), VMR(25)

      COMMON /GETINPUT3/ H,T,P,VMR,NB,NL,MODEL,IAER,V,TAER55,VRTO3,SNO2
      COMMON /GETINPUT5/ NOBS,IFULLCALC,HSURF,DLT,DLT2
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
!      Moved this conversion to GET_INPUT and init_tpvmr - rjh 1/8/2016
!      DO I = 1, NB
!         P(I) = P(I) / 1013.
!      END DO
!C--- End of conversion

!C     Convert the VMR from the ppm unit in the model to absolute unit.
!      Moved this conversion to GET_INPUT and init_tpvmr - rjh 1/8/2016
!      DO 310 I=1,NB
!  310   VMR(I)=VMR(I)*1.0E-06

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
!C        print*, 'K_SURF = ', K_SURF
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
!            write(91,*) i,k,' RJH: P=',p(i)
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
!        write(91,*) i,' RJH: AMTVRT=',amtvrt,damtvt,p(i),p(i+1),vmr(i),vmr(i+1)
 350  CONTINUE

      CLMVAP=AMTVRT/3.34E+22

      WRITE(91,*)'Column vapor amount in model atmosphere from ground'
      WRITE(91,*)'       to space = ', CLMVAP, ' cm'

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

      WRITE(91,*)'Column vapor below plane (CLMVAPP) = ', &
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

!      write(91,*) 'K_SURF, K_PLANE, DVAP_PLANE, DVAP_LAYER = ', &
!              K_SURF, K_PLANE, DVAP_PLANE, DVAP_LAYER
!      write(91,*) 'DP_PLANE, DP_LAYER = ', DP_PLANE, DP_LAYER


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

      INCLUDE 'COMMONS_INC_v2.f'

      DIMENSION VAPVRT(NH2O_MAX), VAP_SLANT(NH2O_MAX)
      DIMENSION MD(12)
      DIMENSION SSH2O(NH2O_MAX)
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
      COMMON /GEOMETRY3/ G_VAP, G_OTHER, G_VAP_EQUIV
      COMMON /GEOMETRY4/VAP_SLANT_MDL
      REAL   MU,MU0, SSH2O_S(NH2O_MAX,2)
      COMMON /GEOMETRY5/MU,MU0,SSH2O_S
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

      DO I = 1, NH2O_MAX
        VAP_SLANT(I) = VAPVRT(I) * 2.0
      END DO

      OBSZNI = SENZN_L2
      OBSPHI = SENAZ_L2
      SOLZNI = SOLZN_L2

      OBSZNI = OBSZNI / RADEG
      OBSPHI = OBSPHI / RADEG

      SOLZNI = SOLZNI/RADEG

!      write(91,*) 'KPLANE=',k_plane

      MU0 = 1./COS(SOLZNI)
      MU  = 1./COS(OBSZNI)
      GGEOM  = MU0 + MU
      write(91,*) 'GGEOM =',GGEOM,' OBSZNI = ',OBSZNI, ' OBSPHI = ',OBSPHI, &
       'solzni=',solzni,' degrees :: MU0, MU = ',MU0, MU
      GCO2   = GGEOM

      GO3    = GGEOM
      IF(HPLANE.LT.27.) GO3 = GGEOM - 1./COS(OBSZNI)

      GN2O   = GGEOM
      GCO    = GGEOM
      GCH4   = GGEOM
      GO2    = GGEOM

!      write(91,*) 'GGEOM, GCO2, GO3, GN2O, GCO, GCH4, GO2 = '
!      write(91,*)  GGEOM, GCO2, GO3, GN2O, GCO, GCH4, GO2

      TOTLO3 = GO3 * VRTO3

!      WRITE(91,*) 'TOTLO3 = ', TOTLO3, VRTO3

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

!     write(91,*) ' G_VAP, G_OTHER, I ='
!      DO I = 1, 25
!         write(91,*) G_VAP(I), G_OTHER(I), I
!      END DO

!C Calculate the water vapor SCALING factor relative to the total amount
!C    of water vapor in the model atmosphere in the L-shaped
!C    Sun-surface-plane ray path.
!C
      VAP_SLANT_MDL = CLMVAP/COS(SOLZNI) + CLMVAPP/COS(OBSZNI)
      VAP_SOL = CLMVAPP*MU0
      VAP_SEN = CLMVAPP*MU
!      write(91,*) 'VAP_SLANT_MDL =', VAP_SLANT_MDL, ' cm ', SLANT_MDL
!C
!C The "equivalent" geometrical factor corresponding to the total
!C    slant vapor amount of VAP_SLANT_MDL':
!C
      G_VAP_EQUIV = VAP_SLANT_MDL / CLMVAP
       write(91,*) 'G_VAP_EQUIV = ', G_VAP_EQUIV, VAP_SLANT_MDL, CLMVAP
       write(91,*) 'VAP_SOL,VAP_SEN = ', VAP_SOL, VAP_SEN, WTRVPR

      DO 310 I=1,NH2O_MAX
        SSH2O(I) = VAP_SLANT(I) / VAP_SLANT_MDL
        write(91,*) 'SSH2O(I), I = ', SSH2O(I), I, VAP_SLANT_MDL, WTRVPR, ISPLITP
        if (WTRVPR.gt.0) THEN
            SSH2O_S(I,1) = VAPVRT(I) / VAP_SOL
            SSH2O_S(I,2) = VAPVRT(I) / VAP_SEN
            write(91,*) 'SSH2O_S(I,1) = ', I,SSH2O_S(I,1), VAPVRT(I), VAP_SOL, SOLZNI*RADEG
            write(91,*) 'SSH2O_S(I,2) = ', I,SSH2O_S(I,2), VAPVRT(I), VAP_SEN, OBSZNI*RADEG
        else
            if (ISPLITP.ne.0) then
                write(6,*) 'ATREM: Split paths is not working because WaterVapor is 0'
                stop 0
            end if
        end if
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

      INCLUDE 'COMMONS_INC_v2.f'

!C  Common variables
      DIMENSION H(25), T(25), P(25), VMR(25)
      DIMENSION SSH2O(NH2O_MAX)
      DIMENSION WAVOBS(NOBS_MAX),FWHM(NOBS_MAX)
      DIMENSION FINST2(100)
      DIMENSION SUMCF(NP_HI)

!C  Local variables

      COMMON /GETINPUT1/ IH2OVP,ICO2,IO3,IN2O,ICO,ICH4,IO2,INO2
      COMMON /GETINPUT3/ H,T,P,VMR,NB,NL,MODEL,IAER,V,TAER55,VRTO3,SNO2
      COMMON /GETINPUT4/ WAVOBS,FWHM
      COMMON /GETINPUT5/ NOBS,IFULLCALC,HSURF,DLT,DLT2
      COMMON /GETINPUT6/ WNDOW1,WNDOW2,WP94C,WNDOW3,WNDOW4,W1P14C
      COMMON /GETINPUT7/ NB1,NB2,NBP94,NB3,NB4,NB1P14
      COMMON /GEOMETRY2/ GCO2,GO3,GN2O,GCO,GCH4,GO2,SSH2O,TOTLO3,GGEOM

      COMMON /MODEL_ADJ1/ CLMVAP,Q

      COMMON /INIT_SPECCAL3/ NH2O
      DIMENSION DP(25), PM(25), TM(25), VMRM(25)
      COMMON /INIT_SPECCAL5/ DP,PM,TM,VMRM
      COMMON /INIT_SPECCAL6/ IST1,IED1,IST2,IED2,ISTP94,IEDP94
      COMMON /INIT_SPECCAL7/ IST3,IED3,IST4,IED4,IST1P14,IED1P14
      COMMON /INIT_SPECCAL8/ WT1,WT2,WT3,WT4,JA

      COMMON /INIT_SPECCAL10/ NCV2,NCVHF2,NCVTT2,ISTRT2,IEND2,FINST2
      COMMON /INIT_SPECCAL11/ NATOT,NBTOT,NCTOT,NDTOT

      REAL, ALLOCATABLE :: ABSCF_CO2(:,:),ABSCF_O2(:,:),ABSCF_N2O(:,:),ABSCF_CH4(:,:),ABSCF_CO(:,:)
      DIMENSION G_VAP(25), G_OTHER(25)
      COMMON /GEOMETRY3/ G_VAP, G_OTHER, G_VAP_EQUIV
      COMMON /GEOMETRY4/VAP_SLANT_MDL
 
      DIMENSION O3CF(NO3PT)
      COMMON /O3CF_INIT1/ O3CF
 
      DIMENSION TRAN_O3_STD(NO3PT)
      COMMON /INIT_SPECCAL16/ TRAN_O3_STD
 
      DIMENSION RNO2CF(NO3PT)
      COMMON /NO2CF_INIT1/ RNO2CF
 
      DIMENSION TRAN_NO2_STD(NO3PT)
      COMMON /INIT_SPECCAL17/ TRAN_NO2_STD

      COMMON /MODEL_ADJ4/ K_SURF
      INTEGER start(2)/1,1/
      INTEGER cnt(2)/NP_HI,19/
      CHARACTER(len=4096) :: filename
      INTEGER FINDMATCH
      external findMatch
      INTEGER FIRST/1/
      SAVE FIRST,ABSCF_CO2,ABSCF_O2,ABSCF_N2O,ABSCF_CH4,ABSCF_CO

      IF (FIRST.eq.1) THEN

!      CALL get_environment_variable("OCDATAROOT", datadir)
!      WRITE (*,*) 'OCDATAROOT='//TRIM(homedir)

      NH2O = NH2O_MAX         !Number of column water vapor values
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
!      DO I = 1, NP_HI
!         WAVNO_HI(I)  = 3000. + FLOAT(I-1)*DWAVNO   ! Wavenumber (cm-1)of high
!                                                    !    resolution spectrum,
!         TRAN_HI(I)   = 1.0                         !    3000 - 18000 cm-1
!      END DO

!      TRAN_HI(:,:)  = 1.0

!C---        print*,'WAVNO_HI ,1, 10000,NP_HI = ',WAVNO_HI(1),
!C---     &     WAVNO_HI(10000), WAVNO_HI(NP_HI)
!C
!C
      DO I = 1, NP_MED
         WAVLN_MED(I) = VSTART + FLOAT(I-1)*DWAVLN  ! Wavelength of medium 
                                                    ! resolution spectrum,
!C         TRAN_MED(I,:)  = 1.0                         ! FWHM=.2 nm, .56-3.1 um.
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
!C         TRAN_STD(I,:)   = 1.0                        ! FWHM=.2 nm, .3-3.1 um.
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


!      DO I = 1, NP_MED
!         FWHM_WAVNO(I) = 10000.*DLT_MED &
!                         /(WAVLN_MED_INDEX(I)*WAVLN_MED_INDEX(I))
!      END DO
!C-----
!C---        print*,'FWHM_WAVNO ,1, 10000,NP_MED = ',FWHM_WAVNO(1),
!C---     &     FWHM_WAVNO(10000), FWHM_WAVNO(NP_MED)
!C-----
!C

!      DO I = 1, NP_MED
!         NCVHF_WAVNO(I) = ( FACDLT * FWHM_WAVNO(I) / DWAVNO + 1.)
!      END DO
!C-----
!C---        print*,'NCVHF_WAVNO ,1, 10000,NP_MED = ',NCVHF_WAVNO(1),
!C---     &     NCVHF_WAVNO(10000), NCVHF_WAVNO(NP_MED)
!C-----

!C Initialize arrays for smoothing medium resolution spectrum (DLT_MED = 0.2 nm,
!C         and point spacing DWAVLN = 0.0001 micron) to coarser spectral
!C         resolution data from imaging spectrometers.

      DO I = 1, NOBS
         NCVHF(I) = ( FACDLT * FWHM(I) / DWAVLN + 1.)
!         write(91,*) I,' RJH: NCVHF=',NCVHF(I),FACDLT,FWHM(I),DWAVLN
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
          FINST2(I)=EXP(-CONST1*(FLOAT(I-NCVHF2)*DWVAVR/DLT2)**2)
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
      END IF

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
      iwndw1c = findMatch(wavobs,nobs,wndow1)
!      write(91,*) 'RJH: IWNDW1=',IWNDW1,iwndw1c,' list[0]=',wavobs(1)
      IWNDW2=FINDMATCH(WAVOBS,NOBS,WNDOW2)
!      write(91,*) 'RJH: IWNDW2=',IWNDW2

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
!      write(91,*) 'RJH: IST1,IED1,IST2,IED2=',IST1,IED1,IST2,IED2,IWP94C,WP94C
      LL=MOD(NBP94,2)
      IF(LL.EQ.0) NBP94=NBP94+1
      NB3HAF=(NBP94-1)/2
      ISTP94=IWP94C-NB3HAF
      IEDP94=IWP94C+NB3HAF
!      write(91,*) 'RJH: ISTP94,IEDP94=',ISTP94,IEDP94

      WT1=(WNDOW2-WP94C)/(WNDOW2-WNDOW1)
      WT2=(WP94C-WNDOW1)/(WNDOW2-WNDOW1)

      IWNDW4=FINDMATCH(WAVOBS,NOBS,WNDOW3)
      IWNDW5=FINDMATCH(WAVOBS,NOBS,WNDOW4)

      WNDOW3=WAVOBS(IWNDW4)
      WNDOW4=WAVOBS(IWNDW5)

!      write(91,*) 'RJH: WNDOW: ',WT1,WT2,IWNDW4,IWNDW5,WNDOW3,WNDOW4
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
!      write(91,*) 'RJH: IST3,IED3,IST4,IED4=',IST3,IED3,IST4,IED4,IW1P14C,W1P14C
      LL=MOD(NB1P14,2)
      IF(LL.EQ.0) NB1P14=NB1P14+1
      NB6HAF=(NB1P14-1)/2
      IST1P14=IW1P14C-NB6HAF
      IED1P14=IW1P14C+NB6HAF
!      write(91,*) 'RJH: IST1P14,IED1P14=',IST1P14,IED1P14

      WT3=(WNDOW4-W1P14C)/(WNDOW4-WNDOW3)
      WT4=(W1P14C-WNDOW3)/(WNDOW4-WNDOW3)
!      write(91,*) 'RJH: WNDOW2: ',WT3,WT4

      ALLOCATE ( ABSCF_CO2(NP_HI,19) )
      ALLOCATE ( ABSCF_N2O(NP_HI,19) )
      ALLOCATE ( ABSCF_CO(NP_HI,19) )
      ALLOCATE ( ABSCF_CH4(NP_HI,19) )
      ALLOCATE ( ABSCF_O2(NP_HI,19) )

       write(filename(1:dln),'(a)') DATPATH(1:dln)
       write(filename(dln+1:),'(a)') 'abscf_gas.nc'
       ncid = ncopn(filename,NCNOWRIT,IRCODE)

       NRHID = NCVID (NCID, 'abscf_co2', IRCODE)
       CALL NCVGT (NCID, NRHID, START, CNT, ABSCF_CO2, IRCODE)
        if (ircode .ne.0) then
          write(*,*) 'Error reading abscf_gas.nc: abscf_co2: rcode=',ircode
          stop 0
       end if

       NRHID = NCVID (NCID, 'abscf_n2o', IRCODE)
       CALL NCVGT (NCID, NRHID, START, CNT, ABSCF_N2O, IRCODE)
        if (ircode .ne.0) then
          write(*,*) 'Error reading abscf_gas.nc: abscf_n2o: rcode=',ircode
          stop 0
       end if

       NRHID = NCVID (NCID, 'abscf_co', IRCODE)
       CALL NCVGT (NCID, NRHID, START, CNT, ABSCF_CO, IRCODE)
        if (ircode .ne.0) then
          write(*,*) 'Error reading abscf_gas.nc: abscf_co: rcode=',ircode
          stop 0
       end if

       NRHID = NCVID (NCID, 'abscf_ch4', IRCODE)
       CALL NCVGT (NCID, NRHID, START, CNT, ABSCF_CH4, IRCODE)
        if (ircode .ne.0) then
          write(*,*) 'Error reading abscf_gas.nc: abscf_ch4: rcode=',ircode
          stop 0
       end if

       NRHID = NCVID (NCID, 'abscf_o2', IRCODE)
       CALL NCVGT (NCID, NRHID, START, CNT, ABSCF_O2, IRCODE)
        if (ircode .ne.0) then
          write(*,*) 'Error reading abscf_gas.nc: abscf_o2: rcode=',ircode
          stop 0
       end if

    CALL NCCLOS(NCID, RCODE)

    END IF ! IFIRST = 1
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
!            write(91,*) 'TRAN_NO2_STD',i,' ',TRAN_NO2_STD(I),totno2,rno2cf(i)
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

       TRAN_HI_OTHERS(:) = 1.0


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
  
       SUMCF(:) = 0.0
       DO J = K_SURF, 19
          SUMCF(:) = SUMCF(:) - ABSCF_CO2(:,J)*DP(J-K_SURF+1)*VMRM(J-K_SURF+1)
       END DO

        TRAN_HI_OTHERS(:)  = TRAN_HI_OTHERS(:)*EXP(SUMCF(:)*Q* 28.966 / &
                              6.0225E+23 / 1.0E-06)
!        do j=1,NP_HI
!           if (TRAN_HI_OTHERS(j).lt.1)  write(91,*) 'CO2: TRAN_HI_OTHERS= ',j,' ',TRAN_HI_OTHERS(j)
!        end do
      END IF

!C--------------------------------------------
!C       For N2O transmittance calculation.

      IF(IN2O.EQ.1) THEN
!          if (first.eq.1) print*,"Doing N2O..."
          DO 324 I=1,NL
            VMRM(I)=0.3*1.0E-06
            VMRM(I)= VMRM(I)*G_OTHER(I)
  324     CONTINUE
  

      SUMCF(:) = 0
       DO J = K_SURF, 19
           SUMCF(:) = SUMCF(:) - ABSCF_N2O(:,J)*DP(J-K_SURF+1)*VMRM(J-K_SURF+1)
       END DO

        TRAN_HI_OTHERS(:)  = TRAN_HI_OTHERS(:)*EXP(SUMCF(:)*Q* 28.966 / &
                              6.0225E+23 / 1.0E-06)
!        if (first.eq.1) then
!        do i=1,NP_HI
!        print*,i," Tran_hi_others=",TRAN_HI_OTHERS(i)
!        end do
!        endif
!        do j=1,NP_HI
!            if (TRAN_HI_OTHERS(j).lt.1) write(91,*) 'N2O: TRAN_HI_OTHERS= ',j,' ',TRAN_HI_OTHERS(j)
!        end do
      END IF

!C--------------------------------------------
!C       For CO transmittance calculation.
      IF(ICO.EQ.1) THEN
!          if (first.eq.1) print*,"Doing CO..."
          DO 325 I=1,NL
            VMRM(I)=0.1*1.0E-06
            VMRM(I)= VMRM(I)*G_OTHER(I)
  325     CONTINUE
  
       SUMCF(:) = 0.0
       DO J = K_SURF, 19
           SUMCF(:) = SUMCF(:) - ABSCF_CO(:,J)*DP(J-K_SURF+1)*VMRM(J-K_SURF+1)
       END DO
        TRAN_HI_OTHERS(:)  = TRAN_HI_OTHERS(:)*EXP(SUMCF(:)*Q* 28.966 / &
                              6.0225E+23 / 1.0E-06)
!        if (first.eq.1) then
!        do i=1,NP_HI
!        print*,i," Tran_hi_others=",TRAN_HI_OTHERS(i)
!        end do
!        endif
!        do j=1,NP_HI
!            if (TRAN_HI_OTHERS(j).lt.1) write(91,*) 'CO: TRAN_HI_OTHERS= ',j,' ',TRAN_HI_OTHERS(j)
!        end do
      END IF

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


       SUMCF(:) = 0.0
       DO J = K_SURF, 19
          SUMCF(:) = SUMCF(:) - ABSCF_CH4(:,J)*DP(J-K_SURF+1)*VMRM(J-K_SURF+1)
       END DO

        TRAN_HI_OTHERS(:)  = TRAN_HI_OTHERS(:)*EXP(SUMCF(:)*Q* 28.966 / &
                              6.0225E+23 / 1.0E-06)
!        do j=1,NP_HI
!           if (TRAN_HI_OTHERS(j).lt.1) write(91,*) 'CH4: TRAN_HI_OTHERS= ',j,' ',TRAN_HI_OTHERS(j)
!        end do

      END IF

!C--------------------------------------------
!C       For O2 transmittance calculation.
      IF(IO2.EQ.1) THEN
            SUMCF(:) = 0
!C***Modified by Bo-Cai Gao on 2/7/2013 to increase O2 absorption
!C---  coefficients by the factor SCL_O2 for wavelengths > 1.2 micron
!C---  in order to model properly the atmospheric O2 band centered
!C---  near 1.265 micron.

          SCL_O2  = 2.60

          DO 327 I=1,NL
            VMRM(I)=0.21
            VMRM(I)= VMRM(I)*G_OTHER(I)
  327     CONTINUE


!C***Modified by Bo-Cai Gao on 2/7/2013 to increase O2 absorption
!C---  coefficients by the factor SCL_O2 for wavelengths > 1.2 micron
!C---  i& < 1.3333 micron in order to model properly the atmospheric
!C---  O2 band centered near 1.265 micron.
       DO J = K_SURF, 19
          SUMCF(:) = SUMCF(:) - ABSCF_O2(:,J)*DP(J-K_SURF+1)*VMRM(J-K_SURF+1)

       END DO
       TRAN_HI_OTHERS(1:9000) = TRAN_HI_OTHERS(1:9000) * EXP(SUMCF(1:9000)*Q* 28.966 / &
                              6.0225E+23 / 1.0E-06)
       TRAN_HI_OTHERS(9001:106600) = TRAN_HI_OTHERS(9001:106600) * EXP(SUMCF(9001:106600)*Q*SCL_O2* 28.966 / &
                              6.0225E+23 / 1.0E-06)
       TRAN_HI_OTHERS(106601:NP_HI) = TRAN_HI_OTHERS(106601:NP_HI) * EXP(SUMCF(106601:NP_HI)*Q* 28.966 / &
                              6.0225E+23 / 1.0E-06)
!        do j=1,NP_HI
!            write(91,*) 'OX: tran_hi_others ',j,' ',TRAN_HI_OTHERS(j)
!            if (j.eq.94924) then
!               write(*,*) 'OX: ',j,' ',EXP(SUMCF(106601:NP_HI)*Q*SCL_O2* 28.966 / &
!                              6.0225E+23 / 1.0E-06),SUMCF(J)
!               do K=K_SURF,19
!               write(*,*) 'OX: ABSCF: ',j,' ', ABSCF_O2(j,K),DP(K-K_SURF+1),VMRM(K-K_SURF+1)
!               end do
!            end if
!        end do

      END IF


!C--------------------------------------------
!C End of calculation of high resolution transmittances for CO2, N2O, CO,
!C     CH4, and O2.
!C--------------------------------------------
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
     
      FIRST = 0
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
!********************************************************************************
!*  TRANCAL combined with TRAN_TABLE by R. Healy 4/28/2015                      *
!*  Name: TRANCAL                                  *
!*  Purpose: This program calculates combined transmittances of H2O, CO2, O3,   *
!*      N2O, CO, CH4, and O2.                                                   *
!*  Parameters: none.                                  *
!*  Algorithm: The calculations were based on the line-by-line absorption       *
!*      parameters supplied by William R. Ridgway of NASA/GSFC.            *
!*  Global output:VAPTT  - geometrically adjusted water vapor total.           *
!*       R094   - channel ratio for 0.94 um region.            *
!*       R114   - channel ratio for 1.14 um region.            *
!*       TRNCAL - total transmittances of all gases that match the     *
!*                         resolutions of imaging spectrometers.               *
!*  Return Codes: none.                                *
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
!*                                         *
!********************************************************************************


      SUBROUTINE TRAN_TABLE

      INCLUDE 'COMMONS_INC_v2.f'

!C  Common variables
      DIMENSION DP(25), PM(25), TM(25), VMRM(25)
      DIMENSION H(25), T(25), P(25), VMR(25)
      DIMENSION WAVOBS(NOBS_MAX),FWHM(NOBS_MAX)
      COMMON /GETINPUT1/ IH2OVP,ICO2,IO3,IN2O,ICO,ICH4,IO2,INO2
      COMMON /GETINPUT3/ H,T,P,VMR,NB,NL,MODEL,IAER,V,TAER55,VRTO3,SNO2
      COMMON /GETINPUT4/ WAVOBS,FWHM
      COMMON /GETINPUT5/ NOBS,IFULLCALC,HSURF,DLT,DLT2
      COMMON /INIT_SPECCAL3/ NH2O
      COMMON /INIT_SPECCAL5/ DP,PM,TM,VMRM
      COMMON /MODEL_ADJ1/ CLMVAP,Q
      COMMON /MODEL_ADJ4/ K_SURF
      COMMON /MODEL_ADJ3/ K_PLANE, DVAP_PLANE, DVAP_LAYER, &
                          DP_PLANE, DP_LAYER, CLMVAPP

      DIMENSION SSH2O(NH2O_MAX)
      DIMENSION VAPTOT(NH2O_MAX), R0P94(NH2O_MAX), R1P14(NH2O_MAX), TRNTBL(NOBS_MAX,NH2O_MAX), &
                 TRAN_KD(NOBS_MAX,NH2O_MAX), DIFF_TRAN(NOBS_MAX,NH2O_MAX),TRH2(NOBS_MAX,NH2O_MAX), &
                 TRNTBLO(NOBS_MAX)
      COMMON /TRAN_TABLE1/ SH2O,VAPTOT,R0P94,R1P14,TRNTBL,TRAN_KD, DIFF_TRAN,TRNTBLO
      ! Split path transmittance
      ! Just need the values for the four table entries calculated on call to get_atrem
      ! rjh - 1/11/2016
      REAL   TRAN_HI_SA(NP_HI,2),TRAN_HI_SAP1(NP_HI,2),TRAN_HI_SB(NP_HI,2),TRAN_HI_SBP1(NP_HI,2), &
              TRNTBL_S(NOBS_MAX,2),TRNTBL_SO(NOBS_MAX,2)
      COMMON /TRAN_TABLES/TRAN_HI_SA,TRAN_HI_SAP1,TRAN_HI_SB,TRAN_HI_SBP1
      COMMON /TRAN_TABLE_L2GEN/TRNTBL_S,TRNTBL_SO
      !COMMON /TRANCAL1/ VAPTT

      REAL   MU,MU0,SSH2O_S(NH2O_MAX,2),VMRM_S(25,2)
      COMMON /GEOMETRY5/MU,MU0,SSH2O_S

      COMMON /GEOMETRY2/ GCO2,GO3,GN2O,GCO,GCH4,GO2,SSH2O,TOTLO3,GGEOM
      REAL G_VAP(25), G_OTHER(25)
      COMMON /GEOMETRY3/ G_VAP, G_OTHER, G_VAP_EQUIV
      COMMON /GEOMETRY4/VAP_SLANT_MDL

      INTEGER start(3)
      INTEGER cnt(3)
      CHARACTER(len=4096) :: filename
      DIMENSION SUMCF(NP_HI), SUMCF_S(NP_HI,2)
      REAL ABSCF_H2O(NP_HI,19)

      REAL, ALLOCATABLE :: TG(:), TKCDF(:,:,:),DIFFT(:,:),SUM_KD(:),F1(:),F2(:),TKCDFC(:,:,:)
      REAL, ALLOCATABLE :: WAVNO_HI(:)
      CHARACTER*31 GNAM, BNAM,LNAM

      REAL, ALLOCATABLE :: TDP(:),TVMRM(:), TRAN_HI(:,:)
      SAVE TDP, TVMRM

      INTEGER FIRST/1/,IDOSMOOTH/0/
      REAL WV/4.1294627/
      LOGICAL DOINLINEKD/.TRUE./
      SAVE FIRST, ABSCF_H2O, TG, TKCDF, DIFFT,N_G, WV, IDOSMOOTH, DOINLINEKD
      data NL/19/,NG/7/

       ALLOCATE( TRAN_HI(NP_HI,NH2O_MAX) )
       ALLOCATE( WAVNO_HI(NP_HI) )

       if (first.eq.1) then

           if (IFULLCALC.eq.1) then
              print*,'ATREM: **WARNING** Full_calc is on.  Processing will be extremely slow...'
           else
              print*,'ATREM: Full_calc is off.  Processing using k-distribution method...'
           endif

           IDOSMOOTH = ICO2 + IN2O + ICO + ICH4 + IO2 + IO3 + INO2

           write(filename(1:dln),'(a)') DATPATH(1:dln)
           write(filename(dln+1:),'(a)') 'abscf_gas.nc'

           start(1) = 1
           start(2) = 1
           start(3) = 1
           cnt(1)   = NP_HI
           cnt(2)   = 19
           cnt(3)   = 0
           ncid = ncopn(filename,NCNOWRIT,IRCODE)
           NRHID = NCVID (NCID, 'abscf_h2o', IRCODE)
           CALL NCVGT (NCID, NRHID, START, CNT, ABSCF_H2O, IRCODE)
            if (ircode .ne.0) then
              write(*,*) 'Error reading abscf_gas.nc: abscf_h2o: rcode=',ircode
              stop 0
           end if
           cnt(2)   = 1
           NRHID = NCVID (NCID, 'waveno', IRCODE)
           CALL NCVGT (NCID, NRHID, START, CNT, WAVNO_HI, IRCODE)
            if (ircode .ne.0) then
              write(*,*) 'Error reading abscf_gas.nc: waveno: rcode=',ircode
              stop 0
           end if

           CALL NCCLOS(NCID, RCODE)

           if (.not.DOINLINEKD) then
               write(filename(1:dln),'(a)') DATPATH(1:dln)
               write(filename(dln+1:),'(a)') 'hico_atrem_h2o_coef7.nc'

               ncid = ncopn(filename,NCNOWRIT,IRCODE)

               CALL NCINQ(NCID,NDIMS,NVARS,NATTS,IRECDIM,IRCODE)

               if (NDIMS.NE.3) then
                  write(*,*) 'Error Wrong number of dims in file'//trim(filename)// &
                  'Expected 4, but got ',NDIMS
                  stop 0
               end if

               CALL NCDINQ(NCID, 1, GNAM, N_G, IRCODE)

               CALL NCDINQ(NCID, 2, BNAM, N_B, IRCODE)

               CALL NCDINQ(NCID, 3, LNAM, N_L, IRCODE)

               !CALL NCDINQ(NCID, 4, LNAM, N_H, IRCODE)

               if (N_B.NE.NOBS) then
                  write(*,*) ' Bands expected (',NOBS,') do not match number of bands in file'// &
                 '(',N_B,') for file='//filename
                  stop 0
               end if

    !           if (N_H.NE.NH2O) then
    !              write(*,*) ' NH2O expected (',NH2O,') do not match NH2O in file'// &
    !             '(',N_H,') for file='//filename
    !              stop 0
    !           end if
            else
                N_G = 7
                N_L = 19
                N_B = NOBS
            endif ! doinlinekd

           ALLOCATE( TG(N_G) )
           ALLOCATE( TKCDF(N_L,N_B,N_G) )
           !ALLOCATE( DIFFT(N_B,N_H) )
           ALLOCATE( F1(N_B) )
           ALLOCATE( F2(N_B) )
           ALLOCATE( SUM_KD(N_B))

           ALLOCATE( TVMRM(N_L) )
           ALLOCATE( TDP(N_L) )

           if (DOINLINEKD) then
               TG(1:7) = (/ 0.0, 0.379, 0.5106, 0.81, 0.9548, 0.9933, 1.0/)
           else

               start(1) = 1
               start(2) = 1
               start(3) = 1
               cnt(1)   = N_G
               cnt(2)   = 0
               cnt(3)   = 0

               NRHID = NCVID (NCID, 'g', IRCODE)
               CALL NCVGT (NCID, NRHID, START, CNT, TG, IRCODE)
                if (ircode .ne.0) then
                  write(*,*) 'Error reading '//filename//': rcode=',ircode
                  stop 0
               end if

               start(1) = 1
               start(2) = 1
               start(3) = 1
               cnt(1)   = N_L
               cnt(2)   = N_B
               cnt(3)   = N_G

               NRHID = NCVID (NCID, 'k_h2o', IRCODE)
               CALL NCVGT (NCID, NRHID, START, CNT, TKCDF, IRCODE)
                if (ircode .ne.0) then
                  write(*,*) 'Error reading '//filename//': rcode=',ircode
                  stop 0
               end if

    !           start(1) = 1
    !           start(2) = 1
    !           start(3) = 1
    !           cnt(1)   = N_B
    !           cnt(2)   = N_H
    !           cnt(3)   = 0

    !           NRHID = NCVID (NCID, 'diffT', IRCODE)
    !           CALL NCVGT (NCID, NRHID, START, CNT, DIFFT, IRCODE)
    !            if (ircode .ne.0) then
    !              write(*,*) 'Error reading '//filename//': rcode=',ircode
    !              stop 0
    !           end if

            endif


      endif

! Moved from INIT_SPECAL because it seems more logical here. - rjh 01/04/2016
!C Initial water vapor VMRs for repeated use in other subroutines and
!C    adjust layered water vapor VMRM with geometrical factors.
      DO I = 1, NL
         VMRM(I) = (VMR(I)+VMR(I+1))/2.0
!C Scale the VMRM by the two-way path geometrical factors. The geometric
!C           factors, G_VAP, varies with atmospheric layer number for
!C           aircraft observational geometries.
         if (ISPLITP.ne.0) then      ! perform the split path calculations
             VMRM_S(I,1) = VMRM(I)*MU0
             if (I.lt.K_PLANE) THEN
                VMRM_S(I,2) = VMRM(I)*MU
             else
                if (I.eq.K_PLANE) THEN
                    VMRM_S(I,2) = VMRM(I)*MU*DVAP_PLANE/DVAP_LAYER
                else
                    VMRM_S(I,2) = 0
                end if
             END IF
             WRITE(91,*) I,' VMRM SOL, OBS=',VMRM_S(I,1),VMRM_S(I,2)
         end if
         VMRM(I) = VMRM(I)*G_VAP(I)
         WRITE(91,*) 'VMRM=',VMRM(I)
      END DO

      if (IFULLCALC.ne.0.or.FIRST.eq.1) then

        !C For each water vapor amount, calculate the geometrically adjusted water
        !C     vapor amount, the channel ratios, and the transmittance spectrum.
              SUMCF(:)    = 0
              IF (ISPLITP.NE.0) THEN
                SUMCF_S(:,1)    = 0
                SUMCF_S(:,2)    = 0
              END IF
              DO J = K_SURF, 19
                 SUMCF(:) = SUMCF(:) - ABSCF_H2O(:,J)*DP(J-K_SURF+1)*VMRM(J-K_SURF+1)
                 IF (ISPLITP.NE.0) THEN
                    SUMCF_S(:,1) = SUMCF_S(:,1) - ABSCF_H2O(:,J)*DP(J-K_SURF+1)*VMRM_S(J-K_SURF+1,1)
                    SUMCF_S(:,2) = SUMCF_S(:,2) - ABSCF_H2O(:,J)*DP(J-K_SURF+1)*VMRM_S(J-K_SURF+1,2)
                 END IF
              END DO
!           do i=1,NP_HI
!               write(91,*) I,' SUM_CF=',SUMCF(I),NH2O,k_SURF,ABSCF_H2O(i,10),VMRM(10-K_surf+1),DP(10-k_surf+1),SSH2O(30),Q
!           end do

              DO  I=1,NH2O
        !C Calculate water vapor transmittances with different scaling factors:
        !C
        !C Multiplying the high resolution water vapor transmittance with combined
        !C    transmittances of CO2, N2O, CO, CH4, and O2:
        !C

               TRAN_HI(:,I) = EXP(SUMCF(:)*SSH2O(I) * Q * 28.966 / &
                            6.0225E+23 / 1.0E-06)            ! *TRAN_HI_OTHERS(:) - Done in TRAN_SMOOTH_OTHERS
        !C Total amount of water vapor (in unit of cm) corresponding to the spectrum.
               VAPTOT(I)= VAP_SLANT_MDL * SSH2O(I)
!         write(91,*) 'VAPTT=',VAPTOT(I),SSH2O(I),VAP_SLANT_MDL

           END DO
!           do i=1,NP_HI
!               write(91,*) I,' TRAN_HI(30)=',tran_hi(i,30),NH2O,k_SURF,SUMCF(i),SSH2O(30),Q
!           end do

           if (ISPLITP.ne.0) then      ! perform the split path calculations
                DO J=1,2
!            write(91,*) 'TRAN_TABLE: SSH2O_S JA = ',   JA_L2,  SSH2O_S(JA_L2,1)
!            write(91,*) 'TRAN_TABLE: SSH2O_S JA+1 = ', JA_L2+1,SSH2O_S(JA_L2+1,2)
!            write(91,*) 'TRAN_TABLE: SSH2O_S JB = ',   JB_L2,  SSH2O_S(JB_L2,1)
!            write(91,*) 'TRAN_TABLE: SSH2O_S JB+1 = ', JB_L2+1,SSH2O_S(JB_L2+1,2)
                    TRAN_HI_SA(:,J) = EXP(SUMCF_S(:,J)  *SSH2O_S(JA_L2,J) * Q * 28.966 / &
                            6.0225E+23 / 1.0E-06)
                    TRAN_HI_SAP1(:,J) = EXP(SUMCF_S(:,J)*SSH2O_S(JA_L2+1,J) * Q * 28.966 / &
                            6.0225E+23 / 1.0E-06)

                    TRAN_HI_SB(:,J) = EXP(SUMCF_S(:,J)  *SSH2O_S(JB_L2,J) * Q * 28.966 / &
                            6.0225E+23 / 1.0E-06)
                    TRAN_HI_SBP1(:,J) = EXP(SUMCF_S(:,J)*SSH2O_S(JB_L2+1,J) * Q * 28.966 / &
                            6.0225E+23 / 1.0E-06)
                END DO
!                    write(91,*) 'SPLIT: JA,JB=',JA_L2,JB_L2,'TRAN_HI=',TRAN_HI(1,JA_L2), &
!                    ' TRAN_HI_SA:',TRAN_HI_SA(1,1),TRAN_HI_SA(1,2), &
!                    TRAN_HI_SAP1(1,1),TRAN_HI_SAP1(1,2),TRAN_HI_SB(1,1),TRAN_HI_SB(1,2), &
!                    TRAN_HI_SBP1(1,1),TRAN_HI_SBP1(1,2)
!                    write(91,*)  'SSH2O JAJB:',SSH2O_S(JA_L2,1),SSH2O_S(JA_L2,2),SSH2O_S(JB_L2,1),SSH2O_S(JB_L2,2), ' SOL=',    &
!                    SUMCF_S(1,1), SSH2O_S(JA_L2,1), SUMCF_S(1,1)  *SSH2O_S(JA_L2,1) * Q * 28.966 / 6.0225E+23 / 1.0E-06,' SEN=', &
!                    SUMCF_S(1,2), SSH2O_S(JA_L2,2), SUMCF_S(1,2)  *SSH2O_S(JA_L2,2) * Q * 28.966 / 6.0225E+23 / 1.0E-06,' SOL=', &
!                    exp(SUMCF_S(1,1)  *SSH2O_S(JA_L2,1) * Q * 28.966 / 6.0225E+23 / 1.0E-06),' SEN=', &
!                    exp(SUMCF_S(1,2)  *SSH2O_S(JA_L2,2) * Q * 28.966 / 6.0225E+23 / 1.0E-06)

           endif

     endif
     if (IFULLCALC.eq.0) then
! The following code calculates the k-distrubition coefficients
! as opposed to reading them in from a netcdf file.
! -rjh 9/23/2015
        if (FIRST.eq.1.AND.DOINLINEKD) then
           ALLOCATE( TKCDFC(N_L,N_B,N_G) )
           write(6,*) 'Calling F kdist'
           call kdist_gas_abs(tkcdfc,abscf_h2o,NP_HI,WAVNO_HI,wavobs,N_B)
           !write(6,*) 'Calling C kdist'
           !call kdistgasabs(tkcdfc,abscf_h2o,WAVNO_HI,wavobs,NP_HI,N_L,N_B)
           TKCDF = TKCDFC
        end if
        !  Perform the k-distribution calculation Trapezoidal integral for transmittance over the NOBS(NBANDS) bands
           DO  I=1,NH2O
              TRAN_KD(:,I) = 1.0

              DO J = K_SURF, 19

                 SUM_KD(:) = 0.0

                 F1(:) = EXP(-TKCDF(J,:,1)  *DP(J-K_SURF+1)*VMRM(J-K_SURF+1)*SSH2O(I))

                 DO K = 1, N_G-1
                    F2(:) = EXP(-TKCDF(J,:,K+1)*DP(J-K_SURF+1)*VMRM(J-K_SURF+1)*SSH2O(I))
                    SUM_KD(:) = SUM_KD(:)+(F1(:) + F2(:))*(TG(K+1) - TG(K))/2.0
                    F1(:) = F2(:)
                 END DO

                 TRAN_KD(:,I) = TRAN_KD(:,I)*SUM_KD(:)

              END DO

        !C Total amount of water vapor (in unit of cm) corresponding to the spectrum.
               VAPTOT(I)= VAP_SLANT_MDL * SSH2O(I)
           END DO

     endif

      if (IFULLCALC.ne.0.or.FIRST.eq.1) then

!C
!C Smooth the high resolution spectra to resolutions of measured spectrum
!C

            CALL TRAN_SMOOTH(TRAN_HI)

            if (first.eq.1) then
                 DO  I=1,NH2O
                     TRH2(:,I) = TRNTBL(:,I)
!                 DO J=1,NOBS
!                    write(6,*) 'RJH: TRH2: ',TRH2(J,I)
!                 END DO
                 END DO
            endif

       endif

       TRNTBLO(:) = 1.0
       TRNTBL_SO(:,:) = 1.0
       ! Only call the TRAN_SMOOTH_OTHERS if the other trace gases were selected
       IF (IDOSMOOTH.GT.0) THEN
         CALL TRAN_SMOOTH_OTHERS
       ENDIF

       if (IFULLCALC.eq.0) then

           if (FIRST.eq.1) then
               do k=1,NOBS
                    print*,K," FAST)TRNTBL=",TRAN_KD(k,60)*DIFF_TRAN(k,1),TRNTBLO(k)
               end do
           end if
           DO  I=1,NH2O
               TRNTBL(:,I) = TRAN_KD(:,I)*DIFF_TRAN(:,I)*TRNTBLO(:)

                if (first.eq.1) then
                 DO J=1,NOBS
                    write(6,*) 'RJH: DIFF_TRAN: ',I,J,TRNTBL(J,I),TRNTBLO(J),TRH2(J,I),TRAN_KD(J,I),DIFF_TRAN(J,I)
                 END DO
                endif
           END DO

       ELSE
           if (FIRST.eq.1) then
               do k=1,NOBS
                    print*,K," SLOW)TRNTBL=",TRNTBL(k,60),TRNTBLO(k)
               end do
           end if

           DO  I=1,NH2O
               TRNTBL(:,I) = TRNTBL(:,I)*TRNTBLO(:)
           END DO

           TRNTBL_S(:,:) = TRNTBL_S(:,:)*TRNTBL_SO(:,:)

       endif

!C Calculate 3-channel ratio values, R094 and R114, from simulated spectra.
       CALL CHNLRATIO
!       CALL channelRatio ! C version

       FIRST=0
       DEALLOCATE(TRAN_HI,WAVNO_HI)
      RETURN
      END

      subroutine linterp(x,y,n,xi,yi)
      real*4 x(n),y(n),xi,yi
      integer n

      klo=1
      khi=n
      k=0
      do while (khi-klo.gt.1)
        k=(khi+klo)/2
        if(x(k).gt.xi)then
          khi=k
        else
          klo=k
        endif
      end do

      if (khi.eq.klo) then  !extrapolate
         klo = klo - 1
         if (klo.le.0) then
            !print*,"Ooops..."
            klo=1
            !stop 0
         endif
         yi = y(khi) + (xi-x(khi))/(x(khi)-x(klo))*(y(khi)-y(klo))
      else
         yi = y(klo) + (xi-x(klo))/(x(khi)-x(klo))*(y(khi)-y(klo))
      endif

      end

      SUBROUTINE kdist_gas_abs(kcdf,abscf,np_hi,waveno,wavobs,nwave)

      ! R. Healy 9/2015
      ! Construct k coefficients from gas transmission coefficients (from Amir Ibrahim's Matlab Code)

      REAL*4 :: abscf(NP_HI,19),waveno(NP_HI),wavobs(nwave),kcdf(19,nwave,7)
      real*4 alayers(NP_HI,19)
      REAL, ALLOCATABLE :: diflam(:),UV_lam(:),IR_lam(:),uv_dlam(:),ir_dlam(:),lam(:), &
                             dwave(:),dwn(:),wn(:),g_i(:),k_i(:),g(:),kg(:),logk(:),kint(:,:), &
                             k7(:,:,:),kk(:)
      INTEGER, ALLOCATABLE :: wavel_window(:,:)
      INTEGER binnum(19)
      real*4 g7(7)/ 0,    0.379,    0.6,    0.81,    0.9548,    0.9933,    1/
      character*1 x /''/
      integer*1 inul
      EQUIVALENCE (inul,x)

      data Q/2.15199993E+25/  ! Q=# of molecules above the surface at one atmosphere (molecules/cm**2)

      do i=1,19

        alayers(:,i) = abscf(:,i)*Q*28.966/ &
                            6.0225E+23 / 1.0E-06
      end do
! wavobs - array of center wavelengths of instrument
! DLAM_cent - delta center wavelengths

      ALLOCATE( diflam(nwave) )

      dlam = 9999
      wmin = 9999
      wmax = -9999
      do i=1,nwave-1
        diflam(i) = wavobs(i+1) - wavobs(i)
        if (dlam.gt.diflam(i)) dlam = diflam(i)
        if (wmin.gt.wavobs(i)) wmin = wavobs(i)
        if (wmax.lt.wavobs(i)) wmax = wavobs(i)
      end do

      diflam(nwave) = diflam(nwave-1)

      if (wmin.gt.wavobs(nwave)) wmin = wavobs(nwave)
      if (wmax.lt.wavobs(nwave)) wmax = wavobs(nwave)

      ndxwuv = (wmin-0.0001 - 0.3)/dlam + 1
      ndxwir = (3.1- (wmax+0.0001))/diflam(nwave) + 1
      ndxtot = nwave+ndxwuv+ndxwir

      ALLOCATE( UV_lam(ndxwuv) )
      ALLOCATE( IR_lam(ndxwir) )
      ALLOCATE( UV_dlam(ndxwuv) )
      ALLOCATE( IR_dlam(ndxwir) )
      ALLOCATE( lam(ndxtot) )
      ALLOCATE( dwave(ndxtot) )
      ALLOCATE( dwn(ndxtot) )
      ALLOCATE( wn(ndxtot) )
      ALLOCATE( wavel_window(ndxtot,2) )

      k = ndxtot
      swav = 0.3

      do i=1,ndxwuv
        UV_lam(i) = swav
        uv_dlam(i) = diflam(1)
        swav = swav + dlam
        lam(k) = UV_lam(i)
        dwn(k) = 10000.*dlam/(lam(k)**2)
         wn(k) = 10000./lam(k)
        k = k - 1
      end do
      swav = wmin
      ks=k
      do i=1,nwave
        lam(k) = swav
        dwn(k) = 10000.*diflam(i)/(lam(k)**2)
        wn(k) = 10000./lam(k)
        swav = swav + diflam(i)
        k = k - 1
      end do

      swav = wmax+0.0001

      do i=1,ndxwir
        IR_lam(i) = swav
        ir_dlam(i) = diflam(nwave)
        swav = swav + diflam(nwave)
        lam(k) = IR_lam(i)
        dwn(k) = 10000.*diflam(nwave)/(lam(k)**2)
         wn(k) = 10000./lam(k)
        k = k - 1
      end do

      ALLOCATE( kint(ndxtot,7) )
      ALLOCATE( k7(19,ndxtot,7) )

      iw = nwave ! index for wavelngths for output array kcdf

      do i=ndxwir+1,ndxtot-ndxwuv

         kwavdn = -999
         kwavup = -999

         do k=1,NP_HI
           if (kwavdn.lt.0.and.waveno(k).gt.(wn(i)-dwn(i)/2.0)) kwavdn = k
           if (kwavup.lt.0.and.waveno(k).gt.(wn(i)+dwn(i)/2.0)) kwavup = k
         end do

         wavel_window(i,1) = kwavdn
         wavel_window(i,2) = kwavup

         if (kwavup.lt.1.or.kwavdn.lt.1) then
            nsamp = 1
         else
            nsamp = kwavup - kwavdn + 1
         end if

         nbins = nsamp
         ALLOCATE( g_i(nbins+1) )
         ALLOCATE( k_i(nbins+1) )
         ALLOCATE( kk(nsamp+1) )

         do k=1,19
            if (kwavdn.lt.1.or.kwavup.lt.1) then
               kk(:) = 0
            else
               kk(1:nsamp) = alayers(kwavdn:kwavup,k)
            end if

            binnum(k) = nbins
            call ecdf(k_i,g_i,binnum(k),kk,nsamp)

            ALLOCATE( logk(nbins+1) )
            logk(:) = log10(k_i(:))

            do l=1,7


              call linterp(g_i,logk,binnum(k),g7(l),kint(i,l))
              k7(k,i,l) = 10**kint(i,l)

              if (i.gt.ndxwir.and.i.le.ndxtot-ndxwuv) then
                  if (isnan(k7(k,i,l))) then
                     k7(k,i,l) = 0
                  endif
                  kcdf(k,iw,l) = k7(k,i,l)
              end if
            end do

            DEALLOCATE(logk)

         end do
         DEALLOCATE (g_i, stat=ialloerr)
         DEALLOCATE (k_i)
         DEALLOCATE (kk)

         if (i.gt.ndxwir.and.i.le.ndxtot-ndxwuv) iw = iw - 1
      end do

      DEALLOCATE (k7)
      DEALLOCATE (kint)
      DEALLOCATE( UV_lam )
      DEALLOCATE( IR_lam )
      DEALLOCATE( UV_dlam )
      DEALLOCATE( IR_dlam )
      DEALLOCATE( lam )
      DEALLOCATE( dwave )
      DEALLOCATE( dwn )
      DEALLOCATE( wn )
      DEALLOCATE( wavel_window )

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

      SUBROUTINE TRAN_SMOOTH(TRAN_HI)

      INCLUDE 'COMMONS_INC_v2.f'

      REAL :: TRAN_HI(NP_HI,NH2O_MAX)
      DIMENSION WAVOBS(NOBS_MAX),FWHM(NOBS_MAX)
      COMMON /GETINPUT4/ WAVOBS,FWHM

      DIMENSION VAPTOT(NH2O_MAX), R0P94(NH2O_MAX), R1P14(NH2O_MAX), TRNTBL(NOBS_MAX,NH2O_MAX), &
                 TRAN_KD(NOBS_MAX,NH2O_MAX), DIFF_TRAN(NOBS_MAX,NH2O_MAX),TRNTBLO(NOBS_MAX)
      COMMON /TRAN_TABLE1/ SH2O,VAPTOT,R0P94,R1P14,TRNTBL,TRAN_KD, DIFF_TRAN,TRNTBLO

      DIMENSION TRAN_IA(NH2O_MAX),TRAN_IAP1(NH2O_MAX)

      COMMON /GETINPUT5/ NOBS,IFULLCALC,HSURF,DLT,DLT2
      REAL, ALLOCATABLE :: TRAN_STD(:,:),TRAN_MED(:,:)

      REAL   TRAN_HI_SA(NP_HI,2),TRAN_HI_SAP1(NP_HI,2),TRAN_HI_SB(NP_HI,2),TRAN_HI_SBP1(NP_HI,2), &
              TRNTBL_S(NOBS_MAX,2),TRNTBL_SO(NOBS_MAX,2)
      COMMON /TRAN_TABLES/TRAN_HI_SA,TRAN_HI_SAP1,TRAN_HI_SB,TRAN_HI_SBP1
      COMMON /TRAN_TABLE_L2GEN/TRNTBL_S,TRNTBL_SO
      REAL   TRAN_MED_INDEX_SA(NP_MED,2),TRAN_MED_INDEX_SB(NP_MED,2), &
              TRAN_MED_INDEX_SAP1(NP_MED,2),TRAN_MED_INDEX_SBP1(NP_MED,2), &
              TRAN_MED_SA(NP_MED,2),TRAN_MED_SAP1(NP_MED,2),TRAN_MED_SB(NP_MED,2),TRAN_MED_SBP1(NP_MED,2), &
              TRAN_STD_SA(NP_STD,2),TRAN_STD_SAP1(NP_STD,2),TRAN_STD_SB(NP_STD,2),TRAN_STD_SBP1(NP_STD,2),             &
              TRAN_IA_SA(2),TRAN_IA_SAP1(2),TRAN_IA_SB(2),TRAN_IA_SBP1(2), &
              TRAN_IAP1_SA(2),TRAN_IAP1_SAP1(2),TRAN_IAP1_SB(2),TRAN_IAP1_SBP1(2), &
              TRNTBLSA(2),TRNTBLSAP1(2),TRNTBLSB(2),TRNTBLSBP1(2)
      COMMON /TRAN_TABLES1/TRAN_MED_INDEX_SA,TRAN_MED_INDEX_SAP1,TRAN_MED_INDEX_SB,TRAN_MED_INDEX_SBP1, &
              TRAN_MED_SA,TRAN_MED_SAP1,TRAN_MED_SB,TRAN_MED_SBP1, &
              TRAN_STD_SA,TRAN_STD_SAP1,TRAN_STD_SB,TRAN_STD_SBP1

      INTEGER FIRST/1/,IA(NOBS_MAX)
      REAL    FINSTR_WAVNO(5000,NP_MED), FWHM_WAVNO(NP_MED)
      INTEGER NCVHF_WAVNO(NP_MED)

      SAVE FIRST, IA, FINSTR_WAVNO, FWHM_WAVNO, NCVHF_WAVNO

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
    IF (FIRST.eq.1) THEN !RJH
      DO 466 J =1, NP_MED
            FWHM_WAVNO(J) = 10000.*DLT_MED &
                         /(WAVLN_MED_INDEX(J)*WAVLN_MED_INDEX(J))

            NCVHF_WAVNO(J) = ( FACDLT * FWHM_WAVNO(J) / DWAVNO + 1.)

             NCVTOT_WAVNO = 2 * NCVHF_WAVNO(J) - 1
    
             SUMINS = 0.0

!        write(91,*) 'NCVHF_WAVNO ',j,' = ',NCVHF_WAVNO(J)
!        write(91,*) 'FWHM_WAVNO ',j,' = ',FWHM_WAVNO(J)
!        write(91,*) 'WAVLN_MED_INDEX ',j,' = ',WAVLN_MED_INDEX(J)

          DO 560 I = NCVHF_WAVNO(J), NCVTOT_WAVNO
             FINSTR_WAVNO(I,J) = &
                EXP( -CONST1*(FLOAT(I-NCVHF_WAVNO(J))*DWAVNO &
                /FWHM_WAVNO(J))**2)
             SUMINS = SUMINS + FINSTR_WAVNO(I,J)
!             if (j.eq.1) write(91,*) 'F:A:finstr_wavno sumins: ',i,sumins,FINSTR_WAVNO(I,J),I-NCVHF_WAVNO(J)
 560      CONTINUE

          DO 565 I = 1, NCVHF_WAVNO(J)-1
             FINSTR_WAVNO(I,J) = FINSTR_WAVNO(NCVTOT_WAVNO-I+1,J)
             SUMINS = SUMINS + FINSTR_WAVNO(I,J)
!             if (j.eq.1) write(91,*) 'F:B:finstr_wavno sumins: ',i,sumins,FINSTR_WAVNO(NCVTOT_WAVNO-I+1,J), &
!                               NCVTOT_WAVNO-I+1
 565      CONTINUE

          SUMINS = SUMINS * DWAVNO

          DO 570 I = 1, NCVTOT_WAVNO
             FINSTR_WAVNO(I,J) = FINSTR_WAVNO(I,J)*DWAVNO/SUMINS
!             write(91,*) 'finstr_wavno ',i,j,FINSTR_WAVNO(I,J),SUMINS
  570      CONTINUE
  

 466  CONTINUE


!C  Index searching...
!C
       DO J=1,NOBS
         !C Calculate instrumental response functions...

          SUMINS = 0.0

          NCVTOT = 2 * NCVHF(J) - 1
!          if (j.eq.1) write(91,*) 'NCVTOT=',ncvtot,ncvhf(j)
          DO 1560 I = NCVHF(J), NCVTOT
             FINSTR(I,J) = &
                EXP( -CONST1*(FLOAT(I-NCVHF(J))*DWAVLN/FWHM(J))**2)
             SUMINS = SUMINS + FINSTR(I,J)
!             if (j.eq.1) write(91,*) I,' F:A:finstr sumins =',i,sumins,FINSTR(I,J),&
!             I-NCVHF(j),fwhm(j),dwavln,FLOAT(I-NCVHF(J))*DWAVLN/FWHM(J),(FLOAT(I-NCVHF(J))*DWAVLN/FWHM(J))**2,&
!             -CONST1*(FLOAT(I-NCVHF(J))*DWAVLN/FWHM(J))**2,EXP( -CONST1*(FLOAT(I-NCVHF(J))*DWAVLN/FWHM(J))**2)
1560      CONTINUE

          DO 1565 I = 1, NCVHF(J)-1
             FINSTR(I,J) = FINSTR(NCVTOT-I+1,J)
             SUMINS = SUMINS + FINSTR(I,J)
!             if (j.eq.1) write(91,*) 'F:B:finstr sumins: ',i,sumins,finstr(i,j),ncvtot-i+1
1565      CONTINUE

          SUMINS = SUMINS * DWAVLN

          DO 1570 I = 1, NCVTOT
             FINSTR(I,J) = FINSTR(I,J)*DWAVLN/SUMINS
1570      CONTINUE


          CALL HUNT(WAVLN_STD, NP_STD, WAVOBS(J), IA(J))
!          write(91,*) J,' RJH: HUNT: ',NP_STD,WAVOBS(J),IA(J)
       END DO

      FIRST = 0
      END IF

       ALLOCATE(TRAN_STD(NP_STD,NH2O_MAX))
       ALLOCATE(TRAN_MED(NP_MED,NH2O_MAX))
!C*** !!!High resolution transmittances of CO2, N2O, CO, CH4, and O2 should
!C         also be calculated somewhere else (wavelength start = 0.56 micron).
!C***  Here assuming TCO2*TN2O*TCO*TCH4*TO2 is already calculated previously,
!C     i.e., TRANS(I) = TRAN_CO2(I)*TRAN_N2O(I)*TRAN_CO(I)*TRAN_CH4(I)*TRAN_O2(I)
!C      and TRAN_HI(I) = TRAN_HI_H2O(I) * TRANS(I), and TRAN_HI_H2O for varying
!C      water vapor amounts is calculated in this subroutine.

     TRAN_MED_INDEX(:,:)  = 0.0
     IF(ISPLITP.NE.0) THEN
         TRAN_MED_INDEX_SA(:,:)    = 0.0
         TRAN_MED_INDEX_SAP1(:,:)  = 0.0
         TRAN_MED_INDEX_SB(:,:)    = 0.0
         TRAN_MED_INDEX_SBP1(:,:)  = 0.0
     END IF


     DO  J =1, NP_MED

          ndx1 = INDEX_MED(J)-(NCVHF_WAVNO(J)-1)
          ndx2 = INDEX_MED(J)+ NCVHF_WAVNO(J)-1
          ! write(91,*) 'TRAN_MED_INDEX: NDX1,NDX2=',ndx1,ndx2
          DO 491 K = ndx1,ndx2
             TRAN_MED_INDEX(J,:) = TRAN_MED_INDEX(J,:) + TRAN_HI(K,:)*  &
                          FINSTR_WAVNO(K-INDEX_MED(J)+NCVHF_WAVNO(J),J)
!              write(91,*) J,K,' TRAN_MED_INDEX:',TRAN_MED_INDEX(J,30),TRAN_HI(K,30),FINSTR_WAVNO(K-INDEX_MED(J)+NCVHF_WAVNO(J),J), &
!                           INDEX_MED(J),NCVHF_WAVNO(J),K-INDEX_MED(J)+NCVHF_WAVNO(J), &
!                           TRAN_HI(K,30)*FINSTR_WAVNO(K-INDEX_MED(J)+NCVHF_WAVNO(J),J)
              IF(ISPLITP.NE.0) THEN
                 TRAN_MED_INDEX_SA(J,:) = TRAN_MED_INDEX_SA(J,:) + TRAN_HI_SA(K,:)*  &
                          FINSTR_WAVNO(K-INDEX_MED(J)+NCVHF_WAVNO(J),J)
                 TRAN_MED_INDEX_SAP1(J,:) = TRAN_MED_INDEX_SAP1(J,:) + TRAN_HI_SAP1(K,:)*  &
                          FINSTR_WAVNO(K-INDEX_MED(J)+NCVHF_WAVNO(J),J)
                 TRAN_MED_INDEX_SB(J,:) = TRAN_MED_INDEX_SB(J,:) + TRAN_HI_SB(K,:)*  &
                          FINSTR_WAVNO(K-INDEX_MED(J)+NCVHF_WAVNO(J),J)
                 TRAN_MED_INDEX_SBP1(J,:) = TRAN_MED_INDEX_SBP1(J,:) + TRAN_HI_SBP1(K,:)*  &
                          FINSTR_WAVNO(K-INDEX_MED(J)+NCVHF_WAVNO(J),J)
!          if (j.eq.8205) then
!            write(91,*) J,K,' SPLIT: TRAN_MED_INDEX:',TRAN_MED_INDEX_SA(J,1),TRAN_HI_SA(K,1), &
!                           FINSTR_WAVNO(K-INDEX_MED(J)+NCVHF_WAVNO(J),J), &
!                           INDEX_MED(J),NCVHF_WAVNO(J),K-INDEX_MED(J)+NCVHF_WAVNO(J), &
!                           TRAN_HI_SA(K,1)*FINSTR_WAVNO(K-INDEX_MED(J)+NCVHF_WAVNO(J),J)
!          endif
              END IF
 491      CONTINUE
     END DO

!C
!C Linear interpolation to get TRAN_MED from TRAN_MED_INDEX:
!C     (Note that WAVLN_MED_INDEX(J) >= WAVLN_MED(J)    )
!C
         TRAN_MED(1,:)      = TRAN_MED_INDEX(1,:)
         TRAN_MED(NP_MED,:) = TRAN_MED_INDEX(NP_MED,:)
         IF(ISPLITP.NE.0) THEN
            TRAN_MED_SA(1,:)        = TRAN_MED_INDEX_SA(1,:)
            TRAN_MED_SA(NP_MED,:)   = TRAN_MED_INDEX_SA(NP_MED,:)
            TRAN_MED_SAP1(1,:)      = TRAN_MED_INDEX_SAP1(1,:)
            TRAN_MED_SAP1(NP_MED,:) = TRAN_MED_INDEX_SAP1(NP_MED,:)
            TRAN_MED_SB(1,:)        = TRAN_MED_INDEX_SB(1,:)
            TRAN_MED_SB(NP_MED,:)   = TRAN_MED_INDEX_SB(NP_MED,:)
            TRAN_MED_SBP1(1,:)      = TRAN_MED_INDEX_SBP1(1,:)
            TRAN_MED_SBP1(NP_MED,:) = TRAN_MED_INDEX_SBP1(NP_MED,:)
         ENDIF

      DO J = 2, NP_MED-1
         IF(WAVLN_MED_INDEX(J).LE.WAVLN_MED(J)) THEN
           TRAN_MED(J,:) = TRAN_MED_INDEX(J,:)
!                write(91,*) ' TRAN_MED_WAVLN: ',j,30,tran_med(j,30),tran_med_index(j,30),&
!                        wavln_med_index(j),wavln_med(j)
           IF(ISPLITP.NE.0) THEN
                TRAN_MED_SA(J,:)   = TRAN_MED_INDEX_SA(J,:)
                TRAN_MED_SAP1(J,:) = TRAN_MED_INDEX_SAP1(J,:)
                TRAN_MED_SB(J,:)   = TRAN_MED_INDEX_SB(J,:)
                TRAN_MED_SBP1(J,:) = TRAN_MED_INDEX_SBP1(J,:)
           ENDIF
         ELSE
           DLT  =  WAVLN_MED_INDEX(J) - WAVLN_MED_INDEX(J-1)
           FJM1 = (WAVLN_MED_INDEX(J) - WAVLN_MED(J))        /DLT
           FJ   = (WAVLN_MED(J)       - WAVLN_MED_INDEX(J-1))/DLT
           TRAN_MED(J,:) = FJM1*TRAN_MED_INDEX(J-1,:) + FJ*TRAN_MED_INDEX(J,:)
!               write(91,*) J,30,' TRAN_MED_INTERP:',TRAN_MED_INDEX(J-1,30),TRAN_MED_INDEX(J,30),FJM1,FJ,DLT

           IF(ISPLITP.NE.0) THEN
                TRAN_MED_SA(J,:)   = FJM1*TRAN_MED_INDEX_SA(J-1,:)   + FJ*TRAN_MED_INDEX_SA(J,:)
                TRAN_MED_SAP1(J,:) = FJM1*TRAN_MED_INDEX_SAP1(J-1,:) + FJ*TRAN_MED_INDEX_SAP1(J,:)
                TRAN_MED_SB(J,:)   = FJM1*TRAN_MED_INDEX_SB(J-1,:)   + FJ*TRAN_MED_INDEX_SB(J,:)
                TRAN_MED_SBP1(J,:) = FJM1*TRAN_MED_INDEX_SBP1(J-1,:) + FJ*TRAN_MED_INDEX_SBP1(J,:)
!                write(91,*) J,'SPLIT TRAN_MED_SA_SOL=',TRAN_MED_SA(J,1),FJM1,FJ,TRAN_MED_INDEX_SA(J-1,1), &
!                TRAN_MED_INDEX_SA(J,1)
           ENDIF
!C---
!C---           print*,j,fjm1,fj
!C---
         END IF
      END DO

!C
!C--- Here multiplying O3 and NO2 spectra and other spectrum at medium resolution:
!C
       DO I = 1, NP_STD
          TRAN_STD(I,:) = 1.
       END DO

       DO I = NPSHIF+1, NP_STD
          TRAN_STD(I,:) = TRAN_STD(I,:)*TRAN_MED(I-NPSHIF,:)
!            write(91,*) 'TRAN_STD: ',30,i,i-NPSHIF,tran_std(i,30),tran_med(i-NPSHIF,30)
       END DO

       IF(ISPLITP.NE.0) THEN
          DO I = 1, NP_STD
           TRAN_STD_SA(I,:)   = 1.
           TRAN_STD_SAP1(I,:) = 1.
           TRAN_STD_SB(I,:)   = 1.
           TRAN_STD_SBP1(I,:) = 1.
          END DO
          DO I = NPSHIF+1, NP_STD
            TRAN_STD_SA(I,:)   = TRAN_STD_SA(I,:)  *TRAN_MED_SA(I-NPSHIF,:)
            TRAN_STD_SAP1(I,:) = TRAN_STD_SAP1(I,:)*TRAN_MED_SAP1(I-NPSHIF,:)
            TRAN_STD_SB(I,:)   = TRAN_STD_SB(I,:)  *TRAN_MED_SB(I-NPSHIF,:)
            TRAN_STD_SBP1(I,:) = TRAN_STD_SBP1(I,:)*TRAN_MED_SBP1(I-NPSHIF,:)
!          write(91,*) I,' TRAN_STD_SA_SOL=',TRAN_STD_SA(I,1),TRAN_MED_SA(I-NPSHIF,1)
          END DO
!          DO I = 1,NP_STD
!           write(91,*) I,' TRAN_STD_SA_SOL=',TRAN_STD_SA(I,1),TRAN_MED_SA(I-NPSHIF,1)
!            write(91,*) I,' TRAN_STD_SAP1_SOL=',TRAN_STD_SAP1(I,1)
!          END DO
       ENDIF

 

!C The 2nd stage of smoothing - smooth the medium resolution spectrum (resolution
!C     of 0.2 nm and point spacing of 0.1 nm) with about 25,000 points to match
!C     the coarser and variable resolution spectrum from imaging spectrometers.
!C
!C Initialize some index parameters:
!C
      DO 1466 J =1, NOBS

             TRNTBL(J,:)       = 0.0
             TRAN_IA(:)        = 0.0
             TRAN_IAP1(:)      = 0.0

       IF(ISPLITP.NE.0) THEN
             TRNTBL_S(J,:)     = 0.0
             TRAN_IA_SA(:)     = 0.0
             TRAN_IAP1_SA(:)   = 0.0
             TRAN_IA_SAP1(:)   = 0.0
             TRAN_IAP1_SAP1(:) = 0.0
             TRAN_IA_SB(:)     = 0.0
             TRAN_IAP1_SB(:)   = 0.0
             TRAN_IA_SBP1(:)   = 0.0
             TRAN_IAP1_SBP1(:) = 0.0
       ENDIF
!C---
!C---        print*,'J= ',j, 'NCVHF =', NCVHF(J), 'NCVTOT=',NCVTOT
    
!C---          IF(j.eq.87) then
!C---           print*, ' J = 87 =', J
!C---           do i = 1, NCVTOT
!C---            print*,i, FINSTR(i)
!C---           end do
!C---          end if

!C
!C---          print*,'J =', J, ' IA =', IA

!
!C  Smoothing...
!C
          DO 1491 K = IA(J)-(NCVHF(J)-1), IA(J)+NCVHF(J)-1
             TRAN_IA(:) = TRAN_IA(:) + TRAN_STD(K,:)* &
                          FINSTR(K-IA(J)+NCVHF(J),J)
!             write(91,*) J,K,ia(j),' TRAN_IA:',TRAN_IA(30),TRAN_STD(K,30), K-IA(J)+NCVHF(J),&
!             FINSTR(K-IA(J)+NCVHF(J),J),NCVHF(J)
             IF(ISPLITP.NE.0) THEN
                 TRAN_IA_SA(:) = TRAN_IA_SA(:)     + TRAN_STD_SA(K,:)*   &
                          FINSTR(K-IA(J)+NCVHF(J),J)
                 TRAN_IA_SAP1(:) = TRAN_IA_SAP1(:) + TRAN_STD_SAP1(K,:)* &
                          FINSTR(K-IA(J)+NCVHF(J),J)
                 TRAN_IA_SB(:) = TRAN_IA_SB(:)     + TRAN_STD_SB(K,:)*   &
                          FINSTR(K-IA(J)+NCVHF(J),J)
                 TRAN_IA_SBP1(:) = TRAN_IA_SBP1(:) + TRAN_STD_SBP1(K,:)* &
                          FINSTR(K-IA(J)+NCVHF(J),J)
!                 write(91,*) K,' TRAN_IA=',TRAN_IA(JA_L2),TRAN_STD(K,JA_L2),FINSTR(K-IA(J)+NCVHF(J),J)
!                 write(91,*) J,K,' TRAN_IA_SA_SOL=',TRAN_IA_SA(1),TRAN_STD_SA(K,1),FINSTR(K-IA(J)+NCVHF(J),J)
!                 write(91,*) K,' TRAN_IA_SAP1_SOL=',TRAN_IA_SAP1(1),TRAN_STD_SAP1(K,1),FINSTR(K-IA(J)+NCVHF(J),J)
             ENDIF

!C---             IF(J.eq.1) then
!C---          print*,'j =', j, 'K = ',K, ' IA= ',IA,
!C---     &       K-IA+NCVHF(J),
!C---     &       FINSTR(K-IA+NCVHF(J))
!C---             End IF

1491      CONTINUE

            IA_P1 = IA(J) + 1
          DO 1492 K = IA_P1-(NCVHF(J)-1), IA_P1+NCVHF(J)-1
             TRAN_IAP1(:) = TRAN_IAP1(:) + TRAN_STD(K,:)* &
                          FINSTR(K-IA_P1+NCVHF(J),J)
             IF(ISPLITP.NE.0) THEN
                TRAN_IAP1_SA(:) = TRAN_IAP1_SA(:)     + TRAN_STD_SA(K,:)*   &
                          FINSTR(K-IA_P1+NCVHF(J),J)
                TRAN_IAP1_SAP1(:) = TRAN_IAP1_SAP1(:) + TRAN_STD_SAP1(K,:)* &
                          FINSTR(K-IA_P1+NCVHF(J),J)
                TRAN_IAP1_SB(:) = TRAN_IAP1_SB(:)     + TRAN_STD_SB(K,:)*   &
                          FINSTR(K-IA_P1+NCVHF(J),J)
                TRAN_IAP1_SBP1(:) = TRAN_IAP1_SBP1(:) + TRAN_STD_SBP1(K,:)* &
                          FINSTR(K-IA_P1+NCVHF(J),J)
 !                write(91,*) J,K,' TRAN_IAP1_SA_SOL=',TRAN_IAP1_SA(1),TRAN_STD_SA(K,1),FINSTR(K-IA(J)+NCVHF(J),J)
             ENDIF
1492      CONTINUE
!C
!C Linear interpolation to get TRNCAL from TRAN_IA and TRAN_IAP1:
!C
           DLT_IA  =  WAVLN_STD(IA_P1) - WAVLN_STD(IA(J))
           FIA     = (WAVLN_STD(IA_P1) - WAVOBS(J)) /DLT_IA
!C          FIA_P1  = (WAVOBS(J)     - WAVLN_STD(IA))/DLT_IA
           FIA_P1  = 1. - FIA
           TRNTBL(J,:) = FIA*TRAN_IA(:) + FIA_P1*TRAN_IAP1(:)
!           do k=1,NH2O_MAX
!           write(91,*) J,K,' TRNTBLA:',TRNTBL(J,k),FIA,TRAN_IA(k),TRAN_IAP1(k), &
!           IA(J),IA_P1,WAVLN_STD(IA(J)),WAVLN_STD(IA_P1),WAVOBS(J),DLT_IA
!           end do
           IF(ISPLITP.NE.0) THEN
            TRNTBLSA(:)   = FIA*TRAN_IA_SA(:)   + FIA_P1*TRAN_IAP1_SA(:)
            TRNTBLSAP1(:) = FIA*TRAN_IA_SAP1(:) + FIA_P1*TRAN_IAP1_SAP1(:)
            TRNTBLSB(:)   = FIA*TRAN_IA_SB(:)   + FIA_P1*TRAN_IAP1_SB(:)
            TRNTBLSBP1(:) = FIA*TRAN_IA_SBP1(:) + FIA_P1*TRAN_IAP1_SBP1(:)
            TRNTBL_S(J,:) = 0.5*(F1A*TRNTBLSA(:) + F2A*TRNTBLSAP1(:) + &
                                 F1B*TRNTBLSB(:) + F2B*TRNTBLSBP1(:) )
!            write(91,*) J,' TRNTBL_TOT=',TRNTBL(J,JA_L2),FIA,FIA_P1,TRAN_IA(JA_L2),TRAN_IAP1(JA_L2)
!            write(91,*) J,IA(j),' TRNTBL_SOL_SA=',TRNTBL_S(J,1),TRNTBLSA(1),TRNTBLSAP1(1),TRAN_IA_SA(1), &
!                             TRAN_IAP1_SA(1),TRAN_IA_SAP1(1),TRAN_IAP1_SAP1(1),F1A,F2A,F1B,F2B
!            write(91,*) J,IA(j),' TRNTBL_SOL_SB=',TRNTBL_S(J,1),TRNTBLSB(1),TRNTBLSBP1(1), &
!                             TRAN_IAP1_SB(1),TRAN_IA_SBP1(1),TRAN_IAP1_SBP1(1),F1A,F2A,F1B,F2B,TRAN_IA_SB(1)
!            write(91,*) J,JA_L2,JB_L2,' TRNTBL_SOL=',TRNTBL_S(J,1),TRNTBL(J,1),F1A,F2A,F1B,F2B,TRNTBLSA(1),TRNTBLSAP1(1), &
!                             TRNTBLSB(1),TRNTBLSBP1(1)
!            write(91,*) J,JA_L2,JB_L2,' TRNTBL_SEN=',TRNTBL_S(J,2),TRNTBL(J,2),F1A,F2A,F1B,F2B,TRNTBLSA(2),TRNTBLSAP1(2), &
!                             TRNTBLSB(2),TRNTBLSBP1(2)
           ENDIF
           IF (IFULLCALC.eq.0) THEN
           do K=1,NH2O_MAX
             !! TRAN_KD(J,K) = TRAN_KD(J,K)*DIFF_TRAN(J,K)
             DIFF_TRAN(J,K) =  (TRNTBL(J,K)-TRAN_KD(J,K))/TRAN_KD(J,K) + 1
             !!DTT =  (TRAN_KD(J,K)-TRNTBL(J,K))/TRNTBL(J,K)
             write(*,*) 'RJH: TRAN2: ',J, K, TRAN_KD(J,K),TRNTBL(J,K),DIFF_TRAN(J,K) !!,100*DTT
           end do

           end if
!C---
!C---       print*,'j=',j,'IA =',IA,'FIA =',FIA,'FIA_P1=',FIA_P1,DLT_IA
!C---
!C
1466  CONTINUE


      DEALLOCATE(TRAN_STD,TRAN_MED)
      RETURN
      END

!********************************************************************************
!*                                                                              *
!*  Name: TRAN_SMOOTH_OTHERS                                                    *
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
!* Added by R. Healy (richard.healy@nasa.gov) 9/29/2015                         *
!*                                                                              *
!* This version smooths the other gases, not H2O                                *
!********************************************************************************

      SUBROUTINE TRAN_SMOOTH_OTHERS

      INCLUDE 'COMMONS_INC_v2.f'

      DIMENSION TRAN_O3_STD(NO3PT)
      COMMON /INIT_SPECCAL16/ TRAN_O3_STD

      DIMENSION TRAN_NO2_STD(NO3PT)
      COMMON /INIT_SPECCAL17/ TRAN_NO2_STD

      DIMENSION WAVOBS(NOBS_MAX),FWHM(NOBS_MAX)
      COMMON /GETINPUT4/ WAVOBS,FWHM

      DIMENSION VAPTOT(NH2O_MAX), R0P94(NH2O_MAX), R1P14(NH2O_MAX), TRNTBL(NOBS_MAX,NH2O_MAX), &
                 TRAN_KD(NOBS_MAX,NH2O_MAX), DIFF_TRAN(NOBS_MAX,NH2O_MAX),TRNTBLO(NOBS_MAX)
      COMMON /TRAN_TABLE1/ SH2O,VAPTOT,R0P94,R1P14,TRNTBL,TRAN_KD, DIFF_TRAN,TRNTBLO
      REAL   TRAN_HI_SA(NP_HI,2),TRAN_HI_SAP1(NP_HI,2),TRAN_HI_SB(NP_HI,2),TRAN_HI_SBP1(NP_HI,2), &
              TRNTBL_S(NOBS_MAX,2),TRNTBL_SO(NOBS_MAX,2)
      COMMON /TRAN_TABLES/TRAN_HI_SA,TRAN_HI_SAP1,TRAN_HI_SB,TRAN_HI_SBP1
      COMMON /TRAN_TABLE_L2GEN/TRNTBL_S,TRNTBL_SO
      REAL   MU,MU0, SSH2O_S(NH2O_MAX,2)
      COMMON /GEOMETRY5/MU,MU0,SSH2O_S

      REAL*4 TRAN_IA,TRAN_IAP1, TRAN_STD_O(NP_STD),TRAN_STD_OS(NP_STD,2),AM(NP_STD),TRAN_IAS(2),TRAN_IAP1S(2)

      COMMON /GETINPUT5/ NOBS,IFULLCALC,HSURF,DLT,DLT2
      REAL TRAN_MED_O(NP_MED),TRAN_MED_INDEX_O(NP_MED)        !Transmittance of medium resolution data
      REAL TRAN_MED_OS(NP_MED,2),TRAN_MED_INDEX_OS(NP_MED,2)
      INTEGER FIRST/1/,IA(NOBS_MAX)
      REAL    FINSTR_WAVNO(5000,NP_MED), FWHM_WAVNO(NP_MED)
      INTEGER NCVHF_WAVNO(NP_MED)

      SAVE FIRST, IA, FINSTR_WAVNO, FWHM_WAVNO, NCVHF_WAVNO

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
    IF (FIRST.eq.1) THEN !RJH
      DO 466 J =1, NP_MED
            FWHM_WAVNO(J) = 10000.*DLT_MED &
                         /(WAVLN_MED_INDEX(J)*WAVLN_MED_INDEX(J))

            NCVHF_WAVNO(J) = ( FACDLT * FWHM_WAVNO(J) / DWAVNO + 1.)

             NCVTOT_WAVNO = 2 * NCVHF_WAVNO(J) - 1
             SUMINS = 0.0

          DO 560 I = NCVHF_WAVNO(J), NCVTOT_WAVNO
             FINSTR_WAVNO(I,J) = &
                EXP( -CONST1*(FLOAT(I-NCVHF_WAVNO(J))*DWAVNO &
                /FWHM_WAVNO(J))**2)
             SUMINS = SUMINS + FINSTR_WAVNO(I,J)
 560      CONTINUE

          DO 565 I = 1, NCVHF_WAVNO(J)-1
             FINSTR_WAVNO(I,J) = FINSTR_WAVNO(NCVTOT_WAVNO-I+1,J)
             SUMINS = SUMINS + FINSTR_WAVNO(I,J)
 565      CONTINUE

          SUMINS = SUMINS * DWAVNO

          DO 570 I = 1, NCVTOT_WAVNO
             FINSTR_WAVNO(I,J) = FINSTR_WAVNO(I,J)*DWAVNO/SUMINS
 570      CONTINUE


 466  CONTINUE


         TRAN_MED_O(:)  = 1.0                         ! FWHM=.2 nm, .56-3.1 um.


!C  Index searching...
!C
       DO J=1,NOBS
         !C Calculate instrumental response functions...

          SUMINS = 0.0

          NCVTOT = 2 * NCVHF(J) - 1

          DO 1560 I = NCVHF(J), NCVTOT
             FINSTR(I,J) = &
                EXP( -CONST1*(FLOAT(I-NCVHF(J))*DWAVLN &
                /FWHM(J))**2)
             SUMINS = SUMINS + FINSTR(I,J)
1560      CONTINUE

          DO 1565 I = 1, NCVHF(J)-1
             FINSTR(I,J) = FINSTR(NCVTOT-I+1,J)
             SUMINS = SUMINS + FINSTR(I,J)
1565      CONTINUE

          SUMINS = SUMINS * DWAVLN

          DO 1570 I = 1, NCVTOT
             FINSTR(I,J) = FINSTR(I,J)*DWAVLN/SUMINS
!            write(91,*) I,J,' FINSTR=',FINSTR(I,J)
1570      CONTINUE


          CALL HUNT(WAVLN_STD, NP_STD, WAVOBS(J), IA(J))
       END DO

      END IF ! IFIRST = 0

     airmass = 1.0/mu0 + 1.0/mu

!C*** !!!High resolution transmittances for CO2, N2O, CO, CH4, and O2 should
!C         also be calculated for wavelength start = 0.56 micron.

     TRAN_MED_INDEX_O(:)  = 0.0
     IF(ISPLITP.NE.0) TRAN_MED_INDEX_OS(:,:)  = 0.0

!     do k=1,NP_HI
!       write(91,*) K,' TRAN_HI_OTHERS (SMOOTH) =',TRAN_HI_OTHERS(K)
!     end do

     DO  J =1, NP_MED

          ndx1 = INDEX_MED(J)-(NCVHF_WAVNO(J)-1)
          ndx2 = INDEX_MED(J)+ NCVHF_WAVNO(J)-1
          !write(91,*) j,' NDX1,NDX2=',ndx1,ndx2,INDEX_MED(J),NCVHF_WAVNO(J)
          DO 491 K = ndx1,ndx2
             TRAN_MED_INDEX_O(J) = TRAN_MED_INDEX_O(J) + TRAN_HI_OTHERS(K)*  &
                          FINSTR_WAVNO(K-INDEX_MED(J)+NCVHF_WAVNO(J),J)
!              write(91,*) J,K,' TRAN_MED_INDEX_O:',TRAN_MED_INDEX_O(J),TRAN_HI_OTHERS(K),&
!                           FINSTR_WAVNO(K-INDEX_MED(J)+NCVHF_WAVNO(J),J), &
!                           INDEX_MED(J),NCVHF_WAVNO(J),K-INDEX_MED(J)+NCVHF_WAVNO(J), &
!                           TRAN_HI_OTHERS(K)*FINSTR_WAVNO(K-INDEX_MED(J)+NCVHF_WAVNO(J),J)
              IF(ISPLITP.NE.0) THEN
                AIRM = -log(TRAN_HI_OTHERS(K))/airmass
                TRAN_HI_OTHERS_SOL = exp(-AIRM/mu0)
                TRAN_HI_OTHERS_SEN = exp(-AIRM/mu)

                TRAN_MED_INDEX_OS(J,1) = TRAN_MED_INDEX_OS(J,1) + TRAN_HI_OTHERS_SOL*  &
                          FINSTR_WAVNO(K-INDEX_MED(J)+NCVHF_WAVNO(J),J)
                TRAN_MED_INDEX_OS(J,2) = TRAN_MED_INDEX_OS(J,2) + TRAN_HI_OTHERS_SEN*  &
                          FINSTR_WAVNO(K-INDEX_MED(J)+NCVHF_WAVNO(J),J)

              ENDIF
!C            if (FIRST.eq.1) print*,j," tran_med_index:", TRAN_MED_INDEX_O(J)
 491      CONTINUE
     END DO

!C
!C Linear interpolation to get TRAN_MED from TRAN_MED_INDEX_O:
!C     (Note that WAVLN_MED_INDEX(J) >= WAVLN_MED(J)    )
!C
         TRAN_MED_O(1)      = TRAN_MED_INDEX_O(1)
         TRAN_MED_O(NP_MED) = TRAN_MED_INDEX_O(NP_MED)
         IF(ISPLITP.NE.0) THEN
            TRAN_MED_OS(1,:)      = TRAN_MED_INDEX_OS(1,:)
            TRAN_MED_OS(NP_MED,:) = TRAN_MED_INDEX_OS(NP_MED,:)
         END IF

      DO J = 2, NP_MED-1
         IF(WAVLN_MED_INDEX(J).LE.WAVLN_MED(J)) THEN
           TRAN_MED_O(J) = TRAN_MED_INDEX_O(J)
           IF(ISPLITP.NE.0) TRAN_MED_OS(J,:) = TRAN_MED_INDEX_OS(J,:)
         ELSE
           DLT  =  WAVLN_MED_INDEX(J) - WAVLN_MED_INDEX(J-1)
           FJM1 = (WAVLN_MED_INDEX(J) - WAVLN_MED(J))        /DLT
           FJ   = (WAVLN_MED(J)       - WAVLN_MED_INDEX(J-1))/DLT
           TRAN_MED_O(J) = FJM1*TRAN_MED_INDEX_O(J-1) + FJ*TRAN_MED_INDEX_O(J)
!           write(91,*) J,'3/31 TRAN_MED_O=',TRAN_MED_O(J),WAVLN_MED_INDEX(J),WAVLN_MED_INDEX(J-1),WAVLN_MED(J),FJM1,FJ,DLT
           IF(ISPLITP.NE.0) &
             TRAN_MED_OS(J,:) = FJM1*TRAN_MED_INDEX_OS(J-1,:) + FJ*TRAN_MED_INDEX_OS(J,:)

         END IF
      END DO

!C
!C--- Here multiplying O3 and NO2 spectra and other spectrum at medium resolution:
!C
       TRAN_STD_O(:) = 1.
       IF(ISPLITP.NE.0) TRAN_STD_OS(:,:) = 1.

       TRAN_STD_O(1:NO3PT) = TRAN_O3_STD(1:NO3PT) * TRAN_NO2_STD(1:NO3PT)
       IF(ISPLITP.NE.0) THEN
           AM(1:NO3PT) = -log(TRAN_STD_O(1:NO3PT))/airmass
           TRAN_STD_OS(1:NO3PT,1) = exp(-AM(1:NO3PT)/mu0)
           TRAN_STD_OS(1:NO3PT,2) = exp(-AM(1:NO3PT)/mu)
       END IF

       TRAN_STD_O(NPSHIF+1:NP_STD) = TRAN_STD_O(NPSHIF+1:NP_STD)*TRAN_MED_O(1:NP_STD-NPSHIF)
       IF(ISPLITP.NE.0) THEN
           TRAN_STD_OS(NPSHIF+1:NP_STD,:) = TRAN_STD_OS(NPSHIF+1:NP_STD,:)*TRAN_MED_OS(1:NP_STD-NPSHIF,:)
!       do j=NPSHIF+1,NP_STD
!        write(91,*) j, ' TRAN_STD_OS_SOL=',TRAN_STD_OS(j,1),TRAN_MED_OS(j-NPSHIF,1)
!        write(91,*) j, ' TRAN_STD_OS_SEN=',TRAN_STD_OS(j,2),TRAN_MED_OS(j-NPSHIF,2)
!       end do
       ENDIF
!C The 2nd stage of smoothing - smooth the medium resolution spectrum (resolution
!C     of 0.2 nm and point spacing of 0.1 nm) with about 25,000 points to match
!C     the coarser and variable resolution spectrum from imaging spectrometers.
!C
!C Initialize some index parameters:
!C
      DO 1466 J =1, NOBS

             TRAN_IA         = 0.0
             TRAN_IAP1       = 0.0
             TRAN_IAS(:)     = 0.0
             TRAN_IAP1S(:)   = 0.0

!
!C  Smoothing...
!C
!          write(91,*) 'TRAN_IA_LOOP:',IA(J)-(NCVHF(J)-1), IA(J)+NCVHF(J)-1,IA(J),NCVHF(J)
          DO 1491 K = IA(J)-(NCVHF(J)-1), IA(J)+NCVHF(J)-1
             TRAN_IA = TRAN_IA + TRAN_STD_O(K)* &
                          FINSTR(K-IA(J)+NCVHF(J),J)
!             write(91,*) J,K,' TRAN_IA_O:',TRAN_IA,TRAN_STD_O(K),FINSTR(K-IA(J)+NCVHF(J),J),IA(J),NCVHF(J)
            IF(ISPLITP.NE.0) THEN
                TRAN_IAS(:) = TRAN_IAS(:) + TRAN_STD_OS(K,:)* &
                          FINSTR(K-IA(J)+NCVHF(J),J)
!             write(91,*) J,K,' TRAN_IAS_SOL:',TRAN_IAS(1),TRAN_STD_OS(K,1),K-IA(J)+NCVHF(J), &
!                          MU0,IA(J),NCVHF(J)
!             write(91,*) J,K,' TRAN_IAS_SEN:',TRAN_IAS(2),TRAN_STD_OS(K,2),K-IA(J)+NCVHF(J), &
!                                       MU,IA(J),NCVHF(J)
            END IF


1491      CONTINUE

            IA_P1 = IA(J) + 1
          DO 1492 K = IA_P1-(NCVHF(J)-1), IA_P1+NCVHF(J)-1
             TRAN_IAP1 = TRAN_IAP1 + TRAN_STD_O(K)* &
                          FINSTR(K-IA_P1+NCVHF(J),J)
           IF(ISPLITP.NE.0) &
             TRAN_IAP1S(:) = TRAN_IAP1S(:) + TRAN_STD_OS(K,:)* &
                          FINSTR(K-IA_P1+NCVHF(J),J)
 1492     CONTINUE
!C
!C Linear interpolation to get TRNCAL from TRAN_IA and TRAN_IAP1:
!C
           DLT_IA  =  WAVLN_STD(IA_P1) - WAVLN_STD(IA(J))
           FIA     = (WAVLN_STD(IA_P1) - WAVOBS(J)) /DLT_IA
!C          FIA_P1  = (WAVOBS(J)     - WAVLN_STD(IA))/DLT_IA
           FIA_P1  = 1. - FIA
           TRNTBLO(J) = FIA*TRAN_IA + FIA_P1*TRAN_IAP1
 !           write(91,*) J,' TRNTBLO:',TRNTBLO(J),FIA,TRAN_IA,TRAN_IAP1, &
 !           IA(J),IA_P1,WAVLN_STD(IA(J)),WAVLN_STD(IA_P1),WAVOBS(J),DLT_IA
           IF(ISPLITP.NE.0) THEN
            TRNTBL_SO(J,:) = FIA*TRAN_IAS(:) + FIA_P1*TRAN_IAP1S(:)
!            write(91,*) 'TRNTBL_SO:',J,TRNTBL_SO(J,1),TRAN_IAS(1), &
!             TRAN_IAP1S(1),FIA,FIA_P1
!            write(91,*) 'TRNTBL_SO:',J,TRNTBL_SO(J,2),TRAN_IAS(2),TRAN_IAP1S(2), &
!             FIA,FIA_P1
           ENDIF

!C---
!C---       print*,'j=',j,'IA =',IA,'FIA =',FIA,'FIA_P1=',FIA_P1,DLT_IA
!C---
!C
1466  CONTINUE


      FIRST = 0

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
      PARAMETER (NOBS_MAX=1024,NH2O_MAX=60)

!C  Common variables
      DIMENSION CONST1(NH2O_MAX),CONST2(NH2O_MAX),CONST3(NH2O_MAX)
      DIMENSION CONST4(NH2O_MAX),CONST5(NH2O_MAX),CONST6(NH2O_MAX)

      COMMON /GETINPUT7/ NB1,NB2,NBP94,NB3,NB4,NB1P14
      COMMON /INIT_SPECCAL6/ IST1,IED1,IST2,IED2,ISTP94,IEDP94
      COMMON /INIT_SPECCAL7/ IST3,IED3,IST4,IED4,IST1P14,IED1P14
      COMMON /INIT_SPECCAL8/ WT1,WT2,WT3,WT4,JA
      DIMENSION VAPTOT(NH2O_MAX), R0P94(NH2O_MAX), R1P14(NH2O_MAX),TRNTBL(NOBS_MAX,NH2O_MAX), &
                 TRAN_KD(NOBS_MAX,NH2O_MAX), DIFF_TRAN(NOBS_MAX,NH2O_MAX),TRNTBLO(NOBS_MAX)
      COMMON /TRAN_TABLE1/ SH2O,VAPTOT,R0P94,R1P14,TRNTBL, TRAN_KD, DIFF_TRAN,TRNTBLO

!C Calculate average of spectra over window and water vapor absorption regions.
      CONST1(:)=0.0
      DO 560 I=IST1,IED1
        CONST1(:)=CONST1(:)+TRNTBL(I,:)
  560 CONTINUE
      CONST1(:)=CONST1(:)/FLOAT(NB1)
      CONST2=0.0
      DO 570 I=IST2,IED2
        CONST2(:)=CONST2(:)+TRNTBL(I,:)
  570 CONTINUE
      CONST2(:)=CONST2(:)/FLOAT(NB2)

      CONST3(:)=0.0
      DO 575 I=ISTP94,IEDP94
        CONST3(:)=CONST3(:)+TRNTBL(I,:)
  575 CONTINUE
      CONST3(:)=CONST3(:)/FLOAT(NBP94)

      R0P94(:)=CONST3(:)/((WT1*CONST1(:)) + (WT2*CONST2(:)))

      CONST4(:)=0.0
      DO 580 I=IST3,IED3 
        CONST4(:)=CONST4(:)+TRNTBL(I,:)
  580 CONTINUE
      CONST4(:)=CONST4(:)/FLOAT(NB3)

      CONST5(:)=0.0
      DO 590 I=IST4,IED4
        CONST5(:)=CONST5(:)+TRNTBL(I,:)
  590 CONTINUE
      CONST5(:)=CONST5(:)/FLOAT(NB4)

      CONST6(:)=0.0
      DO 595 I=IST1P14,IED1P14
        CONST6(:)=CONST6(:)+TRNTBL(I,:)
  595 CONTINUE
      CONST6(:)=CONST6(:)/FLOAT(NB1P14)

      R1P14(:)=CONST6(:)/((WT3*CONST4(:)) + (WT4*CONST5(:)))

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
      PARAMETER (NOBS_MAX=1024)

      DIMENSION XORGN(1050),YORGN(1050),Y2(1050)
      DIMENSION XINT(NOBS_MAX),YINT(NOBS_MAX)
      INTEGER N                               !number of elements in XORGN,YORGN

      COMMON /GETINPUT5/ NOBS,IFULLCALC,HSURF,DLT,DLT2
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

      INTEGER FUNCTION FINDMATCH(LIST,NOBS,ELEM)
      PARAMETER (NOBS_MAX=1024)
      DIMENSION LIST(NOBS_MAX)
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
