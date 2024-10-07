!C
!C---A series of parameters and arrays for initialization of spectral
!C     calculation and smoothing routines.
      INTEGER NP_HI, NP_MED, NP_STD, NH2O_MAX

      PARAMETER (NP_HI=300000)     !Number of high resolution spectral points, 
                                   !  covering 3000 cm-1 to 18,000 cm-1 at
                                   !  point spacing of 0.05 cm-1

      PARAMETER (NP_MED = 25401)   !Number of medium resolution spectral
                                   !  points from 0.56 um to 3.1 um at 0.1 nm
                                   !  point spacing and 0.2 nm full width at
                                   !  half maximum (FWHM)

      PARAMETER (NP_STD = 28001)   !Number of STANDARD medium resolution spec.
                                   !  points from 0.30 um to 3.1 um at 0.1 nm
                                   !  point spacing and 0.2 nm full width at
                                   !  half maximum (FWHM)

      PARAMETER (NPSHIF = 2600)    !NPSHIF = NP_STD - NP_MED
      
      PARAMETER (NH2O_MAX=60)      ! Max Number of column water vapor values in table
      !REAL WAVNO_HI(NP_HI)         !Wavenumber of high resolution data(0.05cm-1)
      REAL WAVLN_MED(NP_MED)       !Wavelength of medium resolution data
                                   ! from 0.56 to 3.1 um (0.1 nm point spacing).
      REAL WAVLN_STD(NP_STD)       !Wavelength of medium resolution data
                                   ! from 0.3 to 3.1 um (0.1 nm point spacing).
      !REAL TRAN_HI(NP_HI,NH2O_MAX)          !Transmittance of high res. data(.05cm-1) for
                                   !  6 gases (H2O, CO2, N2O, CO, CH4, and O2).
      !REAL TRAN_MED(NP_MED,NH2O_MAX)        !Transmittance of medium resolution data
                                   ! from 0.56 to 3.1 um (0.1 nm point spacing).
      !REAL TRAN_STD(NP_STD,NH2O_MAX)        !Transmittance of medium resolution data
                                   ! from 0.30 to 3.1 um (0.1 nm point spacing) 
      REAL TRAN_HI_OTHERS(NP_HI)   !Transmittance of high res. data(.05cm-1) for
                                   !  5 gases (     CO2, N2O, CO, CH4, and O2).
      !REAL A(NP_HI)                !High resolution absorption data (.05cm-1).

!      COMMON /INIT_SPECCAL1/ WAVNO_HI, TRAN_HI, TRAN_HI_OTHERS
      COMMON /INIT_SPECCAL1/ TRAN_HI_OTHERS
!      COMMON /INIT_SPECCAL12/ WAVLN_MED, TRAN_MED, WAVLN_STD, TRAN_STD
      COMMON /INIT_SPECCAL12/ WAVLN_MED, WAVLN_STD


      REAL VSTART, VEND, DWAVLN    !VSTART: Starting wavelength for calculations 
                                   !VEND:   Ending wavelength for calculations
                                   !DWAVLN: Interval between wavelengths
      PARAMETER (VSTART = 0.56, VEND = 3.1, DWAVLN = 0.0001)


      REAL DWVAVR                ! ~ 0.01 um =  wavelength spacing of AVIRIS data

      REAL DWAVNO
      PARAMETER (DWAVNO = 0.05)  !Point spacing of high res. spectra (cm-1).

      REAL DLT_MED               ! 0.2-nm medium resolution spectrum.
      PARAMETER (DLT_MED = 0.0002)
 
      REAL FACDLT                !Factor to multiply DLT by to determine the
      PARAMETER (FACDLT = 2.0)   !  range used in the Gaussian function 
                                 !  calculation for smoothing spectra.

      REAL RADEG
      PARAMETER (RADEG=57.29577951) ! degrees per radian

      INTEGER IH2OLQ               !0 - Indicates that liquid water should not
      PARAMETER (IH2OLQ = 0)       !    be considered in the calculations.
                                   !1 - Indicates that liquid water should be
                                   !    considered in the calculations.
      REAL    RLQAMT               !Amount of "equivalent surface water" in cm,
      PARAMETER (RLQAMT = 0.0)     !    not used yet.

      INTEGER NO3PT                !Number of O3 data points
      PARAMETER (NO3PT = 5001)

      INTEGER NGASTT
      PARAMETER (NGASTT = 7)       !Total number of gases


!C Programmer's note: CONST2=DLT*SQRT(3.1415926/CONST1)  CONST2=total
!C                    area of gaussian curve = DLT*SQRT(pi/(4.0*ln(2)))
      REAL CONST1
      PARAMETER (CONST1 = 2.7725887)    ! CONST1=4.0*ln(2)=2.7725887

      INTEGER INDEX_MED(NP_MED)
      REAL    WAVLN_MED_INDEX(NP_MED), TRAN_MED_INDEX(NP_MED,NH2O_MAX)

      REAL    FINSTR_WAVNO(5000,NP_MED), FWHM_WAVNO(NP_MED)
      INTEGER NCVHF_WAVNO(NP_MED)

      COMMON /INIT_SPECCAL13/ INDEX_MED,WAVLN_MED_INDEX,TRAN_MED_INDEX
      COMMON /INIT_SPECCAL14/ FINSTR_WAVNO, FWHM_WAVNO, NCVHF_WAVNO
 
!C
!C     Arrays for wavelength positions and FWHM of measured imaging spectrometer
!C         data and for smoothing the medium resolution spectrum (FWHM = 0.2 nm,
!C         and point spacing = 0.1 nm) to match coarser resolution spectrum
!C         of imaging spectrometer data

      PARAMETER (NOBS_MAX=1024, NINSTR_MAX = 3001, FILENAME_MAX=4096)

      REAL    FINSTR(NINSTR_MAX,NOBS_MAX)
      INTEGER NCVHF(NOBS_MAX)
      CHARACTER*4096 DATPATH
      INTEGER DLN

      COMMON /INIT_SPECCAL15/ FINSTR, NCVHF
      COMMON /GEOMETRY_L2GEN/ SENZN_L2,SENAZ_L2,SOLZN_L2,WTRVPR, &
             JA_L2,JB_L2,ISPLITP,f1a,f2a,f1b,f2b
      COMMON /INPUT_L2GEN/ DATPATH,DLN

