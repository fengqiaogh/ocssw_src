********************************************************************************
*									       *
*  Name: SOLAR_IRR							       *
*  Purpose:    obtain solar irradiances above the atmosphere corresponding     *
*              to the measurement time and geographic location. 	       *
*  Parameters: none							       *
*  Global used:   WAVOBS - wavelengths of input imaging spectrometer data      *
*                 FWHM   - FWHM for each of the channels		       *
*                 SOLZNI - solar zenith angle				       *
*                 IDAY   - the day number in the year when the measurement     *
*                          was made    					       *
*  Global output: YIRR   - Solar irradiances just above the atmosphere	       *
*  Return Codes:  none							       *
*  Special Considerations: 						       *
*       1. The solar irradiance curve is contained in the file 'sun_binary'.   *
*          It is now based on the "thuillier_atlas3.dat" in the 0.2 - 2.4      *
*          micron range (Thuillier, G., et al., Solar irradiance reference     *
*          spectra for two solar active levels, Adv. in Space Research, 34,    *
*          256-261, 2004). Above 2.4 micron, the solar curve contained         *
*          in the data file "SUN01kurucz2005.dat" inside the MODTRAN 5.2 code  *
*          & multiplied by a factor of 1.02 is adopted. Previously, ATREMi     *
*          used a solar curve converted from the file "sun_kur" in MODTRAN3.5. *
*       2. The solar irradiance curve must have enough entries starting at     *
*          about .30 so that the interpolation and smoothing will work         *
*          correctly. 							       *
*       3. This subroutine doesn't contain the statement INCLUDE 'COMMONS_INC'.*
*          As a result, the variables, such as NP_HI, NP_MED, WAVNO_HI,        *
*          TRAN_HI, WAVLN_MED, TRAN_MED, and DWAVNO are all local variables    *
*          used for smoothing the 1 cm-1 resolution solar irradiance data from *
*          MODTRAN3.5 to the desired resolutions of imaging spectrometers.     *
*          These variables have different meanings in most of other            *
*          subroutines.          					       *
*									       *
********************************************************************************
      SUBROUTINE  SOLAR_IRR_PC

C  Declarations for common variables
      DIMENSION WAVOBS(1024),FWHM(1024)
      COMMON /GETINPUT4/ WAVOBS,FWHM

      DIMENSION YIRR(1024)
      COMMON /SOLAR_IRR1/YIRR

      COMMON /GETINPUT5/ NOBS,HSURF,DLT,DLT2
      COMMON /GEOMETRY1/ SOLZNI,SOLAZ,OBSZNI,OBSPHI,IDAY


C     Arrays for wavelength positions and FWHM of measured imaging spectrometer
C         data and for smoothing the medium resolution spectrum (FWHM = 0.2 nm,
C         and point spacing = 0.1 nm) to match coarser resolution spectrum
C         of imaging spectrometer data

      PARAMETER (NOBS_MAX = 1024, NINSTR_MAX = 3001)

      REAL    FINSTR(NINSTR_MAX)
      INTEGER NCVHF(NOBS_MAX)


      PARAMETER (NP_HI=49933)      !Number of solar irradiance spectral points,
                                   !  covering 51 cm-1 to 49,983 cm-1 at
                                   !  point spacing of 1.00 cm-1

      PARAMETER (NP_MED = 28001)   !Number of medium resolution spectral
                                   !  points from 0.30 um to 3.1 um at 0.1 nm
                                   !  point spacing and 0.2 nm full width at
                                   !  half maximum (FWHM)

      REAL WAVNO_HI(NP_HI)         !Wavenumber of solar irradiance data at a
                                   !  point spacing of 1.0 cm-1 
      REAL WAVLN_MED(NP_MED)       !Wavelength of medium resolution data
                                   ! from 0.30 to 3.1 um (0.1 nm point spacing).
      REAL TRAN_HI(NP_HI)          !Solar irradiance data at 1.0 cm-1 point
                                   ! spacing
      REAL TRAN_MED(NP_MED)        !Solar irradiance data at medium resolution 
                                   ! from 0.30 to 3.1 um (0.1 nm point spacing).

      REAL VSTART, VEND, DWAVLN    !VSTART: Starting wavelength for calculations
                                   !VEND:   Ending wavelength for calculations
                                   !DWAVLN: Interval between wavelengths
      PARAMETER (VSTART = 0.30, VEND = 3.1, DWAVLN = 0.0001)

      REAL DWAVNO
      PARAMETER (DWAVNO = 1.00)  !Point spacing of input solar irradiance
                                 ! data ( 1.0 cm-1).

      REAL DLT_MED               ! 0.2-nm medium resolution spectrum.
      PARAMETER (DLT_MED = 0.0002)

      REAL FACDLT                !Factor to multiply DLT by to determine the
      PARAMETER (FACDLT = 2.0)   !  range used in the Gaussian function


C Programmer's note: CONST2=DLT*SQRT(3.1415926/CONST1)  CONST2=total
C                    area of gaussian curve = DLT*SQRT(pi/(4.0*ln(2)))
      REAL CONST1
      PARAMETER (CONST1 = 2.7725887)    ! CONST1=4.0*ln(2)=2.7725887

      INTEGER INDEX_MED(NP_MED)
      REAL    WAVLN_MED_INDEX(NP_MED), TRAN_MED_INDEX(NP_MED)
C
      REAL    FINSTR_WAVNO(5000), FWHM_WAVNO(NP_MED)
      INTEGER NCVHF_WAVNO(NP_MED)

      DATA PI,DTORAD /3.1415926,0.0174533/

C--- Temp Code ---
        OPEN(31,file='sun_binary_PC',status='old',
     & form='unformatted',access='direct',recl=4*49933)
        read(31,REC=1) (wavno_hi(I), I = 1, 49933)
        read(31,REC=2) (Tran_hi(I),  I = 1, 49933)
        close(31)
C--- End of Temp Code

      DO I = 1, NP_MED
         WAVLN_MED(I) = VSTART + FLOAT(I-1)*DWAVLN  ! Wavelength of medium
                                                    ! resolution spectrum,
         TRAN_MED(I)  = 1.0                         ! FWHM=.2 nm, .30-3.1 um.
      END DO


C Note: The grids of WAVNO_HI do not match the grids of 10000./WAVLN_MED.
C       INDEX_MED is a specially designed index for finding close matches
C       between  the two kinds of grids.
C
      DO I = 1, NP_MED
         INDEX_MED(I) = ( (10000./WAVLN_MED(I) - 51.)/DWAVNO + 1.)
      END DO
C
C Note:     WAVLN_MED_INDEX(I) is very close to WAVLN_MED(I),
C       and WAVLN_MED_INDEX(I) >= WAVLN_MED(I)
C
      DO I = 1, NP_MED
         WAVLN_MED_INDEX(I) = 10000. /(FLOAT(INDEX_MED(I)-1)*DWAVNO
     &                               + 51.)
      END DO
 
C
      DO I = 1, NP_MED
         FWHM_WAVNO(I) = 10000.*DLT_MED
     &                   /(WAVLN_MED_INDEX(I)*WAVLN_MED_INDEX(I))
         IF(FWHM_WAVNO(I).LT.1.0) FWHM_WAVNO(I) = 1.0
      END DO
C

      DO I = 1, NP_MED
         NCVHF_WAVNO(I) = ( FACDLT * FWHM_WAVNO(I) / DWAVNO + 1.)
      END DO

C Initialize arrays for smoothing medium resolution spectrum (FWHM = 0.2 nm,
C         and point spacing = 0.1 nm) to coarser spectral resolution data
C         from imaging spectrometers.

      DO I = 1, NOBS
         NCVHF(I) = ( FACDLT * FWHM(I) / DWAVLN + 1.)
      END DO


C First stage of smoothing - smooth solar irradiance spectrum with 49,933 
C     points at a point spacing of 1 cm-1 to medium resolution spectrum
C     (resolution of 0.2 nm and point spacing of 0.1 nm) with about
C     25,000 points.
C
C     The smoothing is done in wavenumber domain. For a spectrum with a 
C     constant 0.2 nm resolution in wavelength domain, it has variable
C     resolution in wavenumber domain. This effect is properly taken care of
C     in the design of smoothing functions.
C
C     Because the input solar spectrum is in wavenumber units (cm-1), while 
C     the medium resolution spectrum is in wavelength units, the two kinds of 
C     grids do not automatically match. In order to match the grids, arrays 
C     of INDEX_MED and TRAN_MED_INDEX are specially designed. The desired 
C     medium resolution spectrum, TRAN_MED, at constant 0.1 nm point spacing 
C     and 0.2 nm resolution is obtained through linear interpolation of 
C     TRAN_MED_INDEX array. 
C
      DO 466 J = 1, NP_MED

             TRAN_MED_INDEX(J)  = 0.0
             NCVTOT_WAVNO = 2 * NCVHF_WAVNO(J) - 1
    
             SUMINS = 0.0

          DO 560 I = NCVHF_WAVNO(J), NCVTOT_WAVNO
             FINSTR_WAVNO(I) = 
     &           EXP( -CONST1*(FLOAT(I-NCVHF_WAVNO(J))*DWAVNO
     &           /FWHM_WAVNO(J))**2.0)
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
  

          DO 491 K = INDEX_MED(J)-(NCVHF_WAVNO(J)-1), 
     &                            INDEX_MED(J)+NCVHF_WAVNO(J)-1
             TRAN_MED_INDEX(J) = TRAN_MED_INDEX(J) + TRAN_HI(K)* 
     &                     FINSTR_WAVNO(K-INDEX_MED(J)+NCVHF_WAVNO(J))       
 491      CONTINUE

 466  CONTINUE

C
C Linear interpolation to get TRAN_MED from TRAN_MED_INDEX:
C     (Note that WAVLN_MED_INDEX(J) >= WAVLN_MED(J)    )
C
         TRAN_MED(1)      = TRAN_MED_INDEX(1)  
         TRAN_MED(NP_MED) = TRAN_MED_INDEX(NP_MED)  

      DO J = 2, NP_MED-1
           DLT  =  WAVLN_MED_INDEX(J) - WAVLN_MED_INDEX(J-1)
         IF(DLT.LT.1.0E-06) THEN
           TRAN_MED(J) = TRAN_MED_INDEX(J)
         ELSE
           FJM1 = (WAVLN_MED_INDEX(J) - WAVLN_MED(J))        /DLT
           FJ   = (WAVLN_MED(J)       - WAVLN_MED_INDEX(J-1))/DLT
           TRAN_MED(J) = FJM1*TRAN_MED_INDEX(J-1) + FJ*TRAN_MED_INDEX(J)
         END IF
      END DO

 
C Second stage of smoothing - smooth the medium resolution spectrum (resolution 
C     of 0.2 nm and point spacing of 0.1 nm) with about 28,000 points to match
C     the coarser and variable resolution spectrum from imaging spectrometers.
C
C Initialize some index parameters:
      IA = 1000
C
      DO 1466 J =1, NOBS

             YIRR(J) = 0.0
             TRAN_IA     = 0.0
             TRAN_IAP1   = 0.0

             NCVTOT = 2 * NCVHF(J) - 1
    
C Calculate instrumental response functions...
             SUMINS = 0.0

          DO 1560 I = NCVHF(J), NCVTOT
             FINSTR(I) = 
     &           EXP( -CONST1*(FLOAT(I-NCVHF(J))*DWAVLN
     &           /FWHM(J))**2.0)
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
  
C  Index searching...
C
          CALL HUNT(WAVLN_MED, NP_MED, WAVOBS(J), IA)


C  Smoothing...
C
          DO 1491 K = IA-(NCVHF(J)-1), IA+NCVHF(J)-1
             TRAN_IA = TRAN_IA + TRAN_MED(K)* 
     &                     FINSTR(K-IA+NCVHF(J))       
1491      CONTINUE

            IA_P1 = IA + 1
          DO 1492 K = IA_P1-(NCVHF(J)-1), IA_P1+NCVHF(J)-1
             TRAN_IAP1 = TRAN_IAP1 + TRAN_MED(K)* 
     &                     FINSTR(K-IA_P1+NCVHF(J))       
1492      CONTINUE
C
C Linear interpolation to get YIRR from TRAN_IA and TRAN_IAP1:
C
           DLT_IA  =  WAVLN_MED(IA_P1) - WAVLN_MED(IA)
           FIA     = (WAVLN_MED(IA_P1) - WAVOBS(J)) /DLT_IA
C          FIA_P1  = (WAVOBS(J)     - WAVLN_MED(IA))/DLT_IA
           FIA_P1  = 1. - FIA
           YIRR(J) = FIA*TRAN_IA + FIA_P1*TRAN_IAP1
C
1466  CONTINUE


C--          DO I = 1, NOBS
C--             print*,WAVOBS(I),YIRR(I)
C--          END DO

      COSSOL=COS(SOLZNI)

C The Sun-earth distance correction factor, a function of the day in the year
C (see reference: IQBAL, p.3)

      GAMM=2.0*PI*FLOAT(IDAY-1)/365.
      E0=1.000110+0.034221*COS(GAMM)+0.001280*SIN(GAMM)
     & + 0.000719*COS(2.0*GAMM)+0.000077*SIN(2.0*GAMM)

      DO 20 I=1,NOBS
        YIRR(I)=YIRR(I)*COS(SOLZNI)*E0
   20 CONTINUE

C      WRITE(*,*)'SOLZNI=',SOLZNI,'COSSOL= ',COS(SOLZNI),'GAMM=',GAMM,
C     &'IDAY=',IDAY,' E0= ',E0

C--      WRITE(*,*)'SOLZNI=',SOLZNI,'COSSOL= ',COS(SOLZNI),'GAMM=',GAMM,
C--     &'IDAY=',IDAY,' E0= ',E0
C--      DO 23 I=1,NOBS
C--   23   WRITE(*,*)'I=',I,' WAVOBS(I)=',WAVOBS(I),' YIRR(I)=',YIRR(I),
C--     &YIRR(I)/(COS(SOLZNI)*E0)
 

      RETURN
      END
