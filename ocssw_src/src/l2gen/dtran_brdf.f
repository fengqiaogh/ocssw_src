
      subroutine 
     & diff_tran_corr   (Iphase,Sun,View,
     &                  Delphi,Chlor,Ta865,Correct)

C     Iphase  is the aerosol model index (1=O50, 2=O70, ... 16=T99)
C     Iwave   is the wavelength index (1=412nm, 2=443nm, ... 6=670nm)
C     Sun     is the solar zenith angle (DEGREES)
C     View    is the satellite viewing angle (DEGREES)
C     Delphi  is the relative azimuth angle (RADIANS)  <---- NOTE NOTE NOTE
C     Chlor   is the chlorophyll concentration (mg/m^3)
C     Ta      is the aerosol optical depth for the GIVEN WAVELENGTH
C     Tr      is the Rayleigh optical depth for the GIVEN WAVELENGTH
C     Correct is vhe value of (t-t*)/t* for the GIVEN WAVELENGTH


c     The procedure is to apply the standard SeaWiFS/MODIS atmospheric correction 
c     algorithm, which provides tvLw.  During this procedure, the atmospheric 
c     correction algorithm yields two aerosol models (and the associated aerosol 
c     optical depths) and an interpolation parameter (Ratio?) that provides the 
c     fraction of each of the aerosol models to be used in the computation of 
c     the aerosol radiance.  Using the approximation tv = tv* and ts = ts* for 
c     each aerosol model (with the t*'s computed from the aerosol models and 
c     optical thicknesses and already provided in the Gordon and Wang algorithm) 
c     and interpolating between them using the same parameter (Ratio), the 
c     normalized water-leaving radiance and the chlorophyll concentration (Chl) 
c     is computed.  Then, using Chl, the subroutine above diff_tran_corr is 
c     entered for each visible wavelength (412 to 670 nm, SeaWiFS Bands) using 
c     the appropriate optical thicknesses for the each wavelength and the 
c     sun-sensor viewing geometry. This provides the diffuse transmittance 
c     correction factor Correct =(tv-tv*)/tv*.  The corrected diffuse 
c     transmittance is then 
c
c                       (Correct + 1) _ tv*.
c
c
c     This value of tv  is then used with the original tvLw to provide a better 
c     estimate of Lw.  It is expected that using the approximation tv = tv* will 
c     not have much influence on Chl so this probably does not need to be 
c     iterated.  If iteration is required, simply recomputed Chl and enter 
c     diff_tran_corr again.  This procedure will result in more precise 
c     values of Lw. 



      IMPLICIT NONE
      SAVE

      integer NRAD,NPHI,NWAVE,NPHASE,NGAUS,NNG,NG,NUM,NBIG
      integer I, J, M, Ibig, Iwave, Iphase

      PARAMETER (NBIG=10,NRAD=31,NPHI=4)       ! Note NRAD=31 NOT 91
      PARAMETER (NPHASE=16,NWAVE=6)

      PARAMETER (NGAUS=2*NRAD-1, NNG=50, NG=2*NNG)
      PARAMETER (NUM=75)

      REAL Sun, View, Delphi, Ta865, Ta, Tr(NWAVE)
      REAL Correct(NWAVE), Chlor
      REAL  MBRDF(NWAVE,NRAD,NPHI)   ! Same as RAD1 in subroutine, 
                                     ! but has wave index RAD1(NWAVE,NRAD,NPHI)
      REAL  SmMBRDF(NWAVE,NBIG,NPHI) ! Shortened version with Mu averaging
      REAL PHSA(Nphase,Nwave,Nbig,NGAUS,NPHI), PHSR(Nbig,NGAUS,NPHI)

      REAL APHSRADA(NPHASE,NWAVE,NGAUS,NGAUS,NPHI)

      REAL MU(NGAUS),PDIV(NG),PWT(NG), THETA(NGAUS)
      REAL PHSRADA(NGAUS,NGAUS,NPHI), PHSRADR(NGAUS,NGAUS,NPHI)
      REAL TAUA_RAT(NPHASE,NWAVE)

      REAL  RAD1(NRAD,NPHI)
      REAL  tst(2),tdf(2)

      REAL TSTARR(NWAVE,NRAD), TSTARA(NPHASE,NWAVE,NRAD)

      REAL ya1,ya2,yr1,yr2
      REAL za1,za2,zr1,zr2
      REAL x1,z1,v1,w1
      REAL x2,z2,v2,w2
      REAL f1,f2,yl,order,fac,rfres
      REAL aint_p_a, aint_p_r, aint_pl_a, aint_pl_r
      REAL pterm_a, pterm_r, plterm_a, plterm_r
      REAL taur, taua, adelphi
      REAL t_star_r, t_star_a, t_diff_r, t_diff_a
      REAL tstar1, tstar2, tdiff1, tdiff2
      REAL ans1, ans2, ans, slope_ans, slope_tst, slope_tdf
      REAL PI, aindex, fresref, rad, ang, x 

      REAL tstartest1, tstartest2

      INTEGER JP, JDN, JUP, MAXPHI, MGAUS, IVIEW

      CHARACTER INFL1(2)*80,INFL2*80, DUMMY*80

c      common /comfour/ aphsrada, taua_rat
      common /comphase/ phsa, phsr, taua_rat
      COMMON /TSTAR/ TSTARR, TSTARA 

      DATA TR /0.3132, 0.2336, 0.1547, 0.1330, 0.0957, 0.0446/

      RAD(X)=X*PI/180.
      ANG(X)=X*180./PI

      PI=4.0*ATAN(1.0)

      AINDEX = 1.334

      MGAUS  = NGAUS
      MAXPHI = NPHI

      DO I = 1, NGAUS   !NRAD
        THETA(I)    = FLOAT(I-1)*90./Float(NRAD-1)
        MU(I)       = COS(RAD(THETA(I)))
      ENDDO

C     Load the appropriate aerosol quantities 
      call read_partial_phase_integrations ! After the first call to 
                                           ! this it just RETURNs
      call read_tstar              ! After the first call to this it just RETURNs

c       Load the Radiance/BRDF tables

      call Morel_BRDF(Sun, Chlor, MBRDF, SmMBRDF)

C     Find where they are in the tables.

      DO I = 1, NRAD
        IF( VIEW .LT. THETA(I) ) THEN
            IVIEW = I
            GO TO 10
        ENDIF
      ENDDO
10    CONTINUE   ! USE IVIEW and IVIEW -1


399    FORMAT(a)
400    FORMAT( 5(3X, E12.6))


      Ta = Ta865
      DO IWAVE = 1, NWAVE

      DO JUP = IVIEW-1,IVIEW

       adelphi = ang(delphi)
       JDN = MGAUS + 1 - JUP

C      First do the integrals that result from the scattering of Lw toward 
C      the sensor 

   
      YA1 = 0
c      YA2 = 0.
      YR1 = 0.
c      YR2 = 0.
      DO M = 1, MAXPHI
      ORDER = FLOAT(M-1)
      X1 = 0.
c      Z1 = 0.
      V1 = 0.
c      W1 = 0.
       DO I = 1, NBIG
        X2 = PHSA(Iphase,Iwave,I,JUP,M)*SmMBRDF(Iwave,I,M)
c        Z2 = PHSA(Iphase,Iwave,I,JUP,M)
        V2 = PHSR(I,JUP,M)*SmMBRDF(Iwave,I,M)
c        W2 = PHSR(I,JUP,M)
        X1 = X1 + X2
c        Z1 = Z1 + Z2
        V1 = V1 + V2
c       W1 = W1 + W2
       ENDDO

       IF(M.GT.1) THEN 
        X1 = X1*2
        V1 = V1*2
c       Z1 = 0.
c       W1 = 0.
       ENDIF
        YA1 = YA1 + COS(ORDER*DELPHI)*X1  ! Index 1 has Rad1    (for t)
c        YA2 = YA2 + COS(ORDER*DELPHI)*Z1  ! Index 2 is w/o Rad1 (for t*)
        YR1 = YR1 + COS(ORDER*DELPHI)*V1  
c        YR2 = YR2 + COS(ORDER*DELPHI)*W1
      ENDDO

C     Now do the integrals that result from the scattering of Lw back 
C     toward the surface with subsequint Fresnel reflection from the water. 

      ZA1 = 0
c      ZA2 = 0.
      ZR1 = 0.
c      ZR2 = 0.
      DO M = 1, MAXPHI
      ORDER = FLOAT(M-1)
      X1 = 0.
c      Z1 = 0.
      V1 = 0.
c      W1 = 0.
       DO I = 1, NBIG
        X2 = PHSA(Iphase,Iwave,I,JDN,M)*SmMBRDF(Iwave,I,M)
c        Z2 = PHSA(Iphase,Iwave,I,JDN,M)
        V2 = PHSR(I,JDN,M)*SmMBRDF(Iwave,I,M)
c        W2 = PHSR(I,JDN,M)
        X1 = X1 + X2
c        Z1 = Z1 + Z2
        V1 = V1 + V2
c       W1 = W1 + W2
       ENDDO
       IF(M.GT.1) THEN 
        X1 = X1*2
        V1 = V1*2
c       Z1 = 0.
c       W1 = 0.
       ENDIF
        ZA1 = ZA1 + COS(ORDER*DELPHI)*X1  ! Aerosol w    Rad1 for t
c        ZA2 = ZA2 + COS(ORDER*DELPHI)*Z1  ! Aerosol w/o  Rad1 for t*
        ZR1 = ZR1 + COS(ORDER*DELPHI)*V1  ! Rayleigh w   Rad1 for t
c        ZR2 = ZR2 + COS(ORDER*DELPHI)*W1  ! Rayleigh w/o Rad1 for t*
      ENDDO

C     Compute the radiance in the viewing direction from Fourier coefficients

      YL = 0.
      DO M = 1, MAXPHI
        ORDER = FLOAT(M-1)
        FAC = 1.
        IF (M .GT. 1) FAC = 2.
        YL = YL + MBRDF(Iwave,JUP,M)*FAC*COS(ORDER*DELPHI)
      ENDDO

C     Combine the two contributions above the integrals to form the diffuse 
C     reflectance in single scattering

      rfres = fresref(mu(JUP), aindex)

c      aint_p_a  = (YA2 + rfres*ZA2)         /(1.-rfres)
      aint_pl_a = (YA1 + rfres*ZA1)/YL      /(1.-rfres)
c      aint_p_r  = (YR2 + rfres*ZR2)         /(1.-rfres)
      aint_pl_r = (YR1 + rfres*ZR1)/YL      /(1.-rfres)

c      pterm_a  = (1.-aint_p_a/2. )/MU(JUP)     ! aint_p should be *omega0
      plterm_a = (1.-aint_pl_a/2.)/MU(JUP)     !    "
c      pterm_r  = (1.-aint_p_r/2. )/MU(JUP)     
      plterm_r = (1.-aint_pl_r/2.)/MU(JUP)      

c       t_star_a = exp(-pterm_a*ta)
        t_diff_a = exp(-plterm_a*ta*taua_rat(iphase,iwave))
c       t_star_r = exp(-pterm_r*tr)
        t_diff_r = exp(-plterm_r*tr(iwave))
        
            
        if (JUP .EQ. IVIEW-1) THEN 
c            tstar1 = t_star_a * t_star_r
            tstar1 = tstarr(Iwave,Jup)
     & * (tstara(Iphase,Iwave,Jup)**(ta*taua_rat(iphase,iwave)))
            tdiff1 = t_diff_a * t_diff_r
            ans1   = (tdiff1 - tstar1)/(tstar1)
            
        else
c            tstar2 = t_star_a * t_star_r
            tstar2 = tstarr(Iwave,Jup)
     & * (tstara(Iphase,Iwave,Jup)**(ta*taua_rat(iphase,iwave)))
            tdiff2 = t_diff_a * t_diff_r
            ans2   = (tdiff2 - tstar2)/(tstar2)

        endif

        ENDDO ! JUP

        slope_ans   = (ans2  -  ans1)/(Theta(iview)-Theta(iview-1))
        slope_tst   = (tstar2-tstar1)/(Theta(iview)-Theta(iview-1))
        slope_tdf   = (tdiff2-tdiff1)/(Theta(iview)-Theta(iview-1))
        Correct(iwave)   =  ans1   + slope_ans*(view-theta(iview-1))
        
        ENDDO ! IWAVE

C        tst(iphase) = tstar1 + slope_tst*(view-theta(iview-1))
C        tdf(iphase) = tdiff1 + slope_tdf*(view-theta(iview-1))

C        ans = correct

c        print*, tstar1, tstar2
c        print*, adelphi, tst(Iphase), tdf(Iphase), ans 

      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        function fresref(muair,  index)

        real muair, index

        theta1  = acos(muair)
        stheta1 = sin(theta1)
        stheta2 = stheta1/index
        theta2  = asin(stheta2)

        if(theta1 .gt. 0.) then 

                fresref = 0.5*(
     &         ( tan(theta1-theta2)/tan(theta1+theta2) )**2
     &        +( sin(theta1-theta2)/sin(theta1+theta2) )**2
     &             )
        else
                fresref = ( (index-1)/(index+1) )**2
        endif
        return
        end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C     *********************************************************

      subroutine Morel_BRDF(Sun, Chlor, MBRDF, SmMBRDF)

      implicit none

      integer NRAD,NPHI,NSUN,NCHL,NWAVE,NCASE,NBIG
      integer IRAD,IPHI,ISUN,ICHL,IWAVE,ICASE,IBIG
      integer I, M

      PARAMETER (NBIG=10,NRAD=31,NPHI=4)    ! Note NRAD=31 NOT 91
      PARAMETER (NSUN=6,NCHL=6,NWAVE=6,NCASE=NSUN*NCHL*NWAVE)

      REAL  BRDF(NWAVE,NSUN,NCHL,NRAD,NPHI)
      REAL  SmBRDF(NWAVE,NSUN,NCHL,NBIG,NPHI)
      REAL  Theta0(NSUN), Chl(NCHL), Wave(NWAVE), Thetav(NRAD)
      REAL  LChl(NCHL)
      REAL  SUN, Chlor, LChlor, chl_interp, sun_interp
      REAL  INTERP1, INTERP2
      REAL  MBRDF(NWAVE,NRAD,NPHI), SmMBRDF(NWAVE,NBIG,NPHI)

      CHARACTER*80  INFL_FOURIER31
      CHARACTER*80  INFL_FOURIER10
      CHARACTER*20  DUMMY

      CHARACTER filedir*255
      INTEGER   len, lenstr

      INTEGER IONCE

      SAVE

      DATA IONCE /0/

      IF (IONCE .EQ. 0) THEN

        call getenv('OCDATAROOT',filedir)
        if (filedir .eq. '') then
            write(*,*)
     .      '-E- : Environment variable OCDATAROOT undefined'
            call exit(1)
        endif
        len = lenstr(filedir)

      INFL_FOURIER31 = filedir(1:len)//
     .  '/eval/common/dtran_brdf/NEW_Morel_NRAD31-1-EDITED'
      INFL_FOURIER10 = filedir(1:len)//
     .  '/eval/common/dtran_brdf/NEW_Morel_SMALL-1-EDITED'

      OPEN(UNIT=11,FILE=INFL_FOURIER31,STATUS='UNKNOWN')
      OPEN(UNIT=12,FILE=INFL_FOURIER10,STATUS='UNKNOWN')

1     format(a)
2     format(' ', 13(2x,f10.6))
400   FORMAT( 5(3X, E12.6))

      DO iwave = 1, nwave
      DO isun  = 1, nsun
      DO ichl  = 1, nchl
      DO iphi = 1, nphi
         DO irad = 1, nrad
          BRDF(iWAVE,iSUN,iCHL,iRAD,iPHI)=0.
         ENDDO
         DO ibig = 1, nbig
          SmBRDF(iWAVE,iSUN,iCHL,iBIG,iPHI)=0.
         ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C     READ IN THE FILES

      DO iwave = 1, nwave
      DO isun  = 1, nsun
      DO ichl  = 1, nchl
        Read(11, *) wave(iwave), theta0(isun), chl(ichl)
c        Write(6,*)  wave(iwave), theta0(isun), chl(ichl)
        DO iphi = 1, nphi
         Read(11, 1)   DUMMY
c         Write(6, 1)   DUMMY
         Read(11,400) (BRDF(iWAVE,iSUN,iCHL,iRAD,iPHI),irad =1,nrad)
c         Write(6,400) (BRDF(iWAVE,iSUN,iCHL,iRAD,iPHI),irad =1,nrad)
        ENDDO
c        Read(12, *) wave(iwave), theta0(isun), chl(ichl)
        Read(12, *) wave(iwave), theta0(isun), chl(ichl)
c        Write(6,*)  wave(iwave), theta0(isun), chl(ichl)
        DO iphi = 1, nphi
         Read(12, 1)   DUMMY
c         Write(6, 1)   DUMMY
         Read(12,400) (SmBRDF(iWAVE,iSUN,iCHL,iBIG,iPHI),ibig =1,nbig)
c         Write(6,400) (SMBRDF(iWAVE,iSUN,iCHL,iBIG,iPHI),ibig =1,nbig)
        ENDDO
      ENDDO
      ENDDO
      ENDDO

      IONCE = 1
      ENDIF

      DO I = 1, NCHL
       LCHL(I) = alog10(Chl(i))
      ENDDO

C     For sun angle interpolation
      DO I = 1, NSUN
        IF( SUN .LT. THETA0(I) ) THEN
            ISUN = I
            GO TO 10
        ENDIF
      ENDDO
10    CONTINUE   ! USE ISUN and ISUN -1

C     For Chl interpolation
      DO I = 1, NCHL
        IF( Chlor .LT. Chl(I) ) THEN
            ICHL = I
            GO TO 11
        ENDIF
        ICHL = NCHL
      ENDDO
11    CONTINUE   ! USE ICHL and ICHL -1

c       Compute the log10 of the chlorophyll
        LChlor = alog10(Chlor)


C       Interpolate in the tables to get the water radiance

        sun_interp = (SUN   - Theta0(ISUN-1) )/(Theta0(ISUN)-Theta0(ISUN-1))
        chl_interp = (LChlor - LChl(ICHL-1) )/(LChl(ICHL)-LChl(ICHL-1))
c        PRINT*, 'Interpolation terms ', sun_interp, chl_interp
        
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012

        DO Iwave = 1, NWAVE
        DO I     = 1, NRAD
        DO M     = 1, NPHI
            INTERP1 = (1.-sun_interp)*(1.-chl_interp)*
     &                BRDF(iWAVE,iSUN-1,iCHL-1,I,M)
            INTERP1 = INTERP1 + sun_interp*(1.-chl_interp)*
     &                BRDF(iWAVE,iSUN,iCHL-1,I,M)
            INTERP1 = INTERP1 + (1.-sun_interp)*chl_interp*
     &                BRDF(iWAVE,iSUN-1,iCHL,I,M)
            INTERP1 = INTERP1 + sun_interp*chl_interp*
     &                BRDF(iWAVE,iSUN,iCHL,I,M)
            MBRDF(Iwave, I, M) = INTERP1

            IF(I .LE. NBIG) THEN
              INTERP2 = (1.-sun_interp)*(1.-chl_interp)*
     &                SmBRDF(iWAVE,iSUN-1,iCHL-1,I,M)
              INTERP2 = INTERP2 + sun_interp*(1.-chl_interp)*
     &                SmBRDF(iWAVE,iSUN,iCHL-1,I,M)
              INTERP2 = INTERP2 + (1.-sun_interp)*chl_interp*
     &                SmBRDF(iWAVE,iSUN-1,iCHL,I,M)
              INTERP2 = INTERP2 + sun_interp*chl_interp*
     &                SmBRDF(iWAVE,iSUN,iCHL,I,M)
              SmMBRDF(Iwave, I, M) = INTERP2
            ENDIF
c        print*,  'Iwave,I,M,MBRDF = ',Iwave,I,M,MBRDF(Iwave, I, M)
        ENDDO
        ENDDO
        ENDDO

        DO Iwave = 1, NWAVE
        DO I     = 1, NBIG
        DO M     = 1, NPHI
c        print*,  'Iwave,I,M,SmMBRDF = ',Iwave,I,M,SmMBRDF(Iwave, I, M)        
        ENDDO
        ENDDO
        ENDDO


      RETURN
      END 
      


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine read_partial_phase_integrations

      IMPLICIT NONE

      integer NRAD,NPHI,NWAVE,NPHASE,NGAUS,NBIG

      PARAMETER (NRAD=31,NPHI=4,NBIG=10,NGAUS=2*NRAD-1)  !Note NRAD=31 NOT 91
      PARAMETER (NPHASE=16,NWAVE=6)

      integer I, J, JUP, M, Iwave, Iphase, Ibig, IONCE

      REAL PHSA(Nphase,Nwave,Nbig,NGAUS,NPHI), PHSR(Nbig,NGAUS,NPHI)
      REAL TAUA_RAT(NPHASE,NWAVE)

      CHARACTER INFL_AER*512,INFL_RAY*512, DUMMY*80
      CHARACTER INFL_RAT*512

      character filedir*255
      integer   len, lenstr

      common /comphase/ phsa, phsr, taua_rat

      DATA  IONCE /0/

      Save

      IF (IONCE .EQ. 0) THEN

         call getenv('OCDATAROOT',filedir)
         if (filedir .eq. '') then
            write(*,*)
     .      '-E- : Environment variable OCDATAROOT undefined'
            call exit(1)
         endif
         len = lenstr(filedir)
         filedir = filedir(1:len)//'/eval/common/dtran_brdf/'
         len = lenstr(filedir)

      INFL_AER = filedir(1:len)//'Aerosols_Partial_Inegr.dat'
      INFL_RAY = filedir(1:len)//'Rayleigh_Partial_Inegr.dat'
      INFL_RAT = filedir(1:len)//'spec_var_EDITED.dat'

      OPEN(UNIT=11,FILE=INFL_AER,STATUS='OLD')
      OPEN(UNIT=12,FILE=INFL_RAY,STATUS='OLD')
      OPEN(UNIT=15,FILE=INFL_RAT,STATUS='OLD')

399    FORMAT(a)
400    FORMAT( 5(3X, E12.6))


      DO Iphase = 1,NPHASE
       DO Iwave = 1,NWAVE
        DO M = 1, NPHI
         DO JUP = 1, NGAUS  ! Read in for thetav = 0, 3, 6, 9, etc.
           Read(11,399) Dummy
c           Write(6,399) Dummy
           Read(11,400) (PHSA(Iphase,Iwave,Ibig,JUP,M), Ibig = 1, NBIG)
c           Write(6,400) (PHSA(Iphase,Iwave,Ibig,JUP,M), Ibig = 1, NBIG)
         ENDDO      !JUP
        ENDDO       !M
       ENDDO        !Iwave
      ENDDO         !Iphase

      CLOSE(11)

       DO M = 1, NPHI
        DO JUP = 1, NGAUS  ! Read in for thetav = 0, 3, 6, 9, etc.
           Read (12,399) Dummy
c           Write( 6,399) Dummy
           Read (12,400) (PHSR(Ibig,JUP,M), Ibig = 1, NBIG)
c           Write( 6,400) (PHSR(Ibig,JUP,M), Ibig = 1, NBIG)
        ENDDO      !JUP
       ENDDO       !M

       CLOSE(12)

C     READ THE TAUA RATIOS (Lambda to 865)

        do iphase = 1, nphase
         read(15,399) Dummy
C         write(6,399) Dummy
         do iwave = 1, nwave
            read(15,*) j, taua_rat(iphase,iwave)
C            write(6,*) iphase, iwave, j, taua_rat(iphase,iwave)
         enddo
        enddo
        close(15)

       IONCE = 1
       ENDIF
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE read_tstar

      IMPLICIT NONE

      integer NRAD,NPHI,NWAVE,NPHASE
      integer I, J, Iwave, Iphase, IONCE

      PARAMETER (NRAD=31)       ! Note NRAD=31 NOT 91
      PARAMETER (NPHASE=16,NWAVE=6)

      REAL TSTARR(NWAVE,NRAD), TSTARA(NPHASE,NWAVE,NRAD)

      CHARACTER*80 Dummy

      character filedir*255
      integer   len, lenstr

      COMMON /TSTAR/ TSTARR, TSTARA 

      DATA IONCE /0/

      SAVE

 397   Format(' ' , '  Iwave = ', I2)
 398   Format(' ' , ' Iphase = ', I2, '  Iwave = ', I2)
 399   FORMAT(a)
 400   FORMAT( 5(3X, E12.6))
      

      If (IONCE .EQ. 0) THEN

         call getenv('OCDATAROOT',filedir)
         if (filedir .eq. '') then
            write(*,*)
     .      '-E- : Environment variable OCDATAROOT undefined'
            call exit(1)
         endif
         len = lenstr(filedir)
         filedir = filedir(1:len)//'/eval/common/dtran_brdf/'
         len = lenstr(filedir)

       OPEN(UNIT=21,FILE=filedir(1:len)//'tstar_rayleigh.dat',
     .      STATUS='UNKNOWN')
       OPEN(UNIT=22,FILE=filedir(1:len)//'tstar_aerosol.dat',
     .      STATUS='UNKNOWN')

        DO Iphase = 1, Nphase
         DO Iwave = 1, Nwave
           If(Iphase .eq.1) then 
              Read(21,399) Dummy
              Read(21,400) (TSTARR(Iwave, J), J = 1, Nrad-3)
c              Write(6,397) Iwave
c              Write(6,400) (TSTARR(Iwave, J), J = 1, Nrad-3)
           Endif
              Read(22,399) Dummy
              Read(22,400) (TSTARA(Iphase,Iwave, J), J = 1, Nrad-3) 
c              Write(6,398) Iphase, Iwave
c              Write(6,400) (TSTARA(Iphase,Iwave, J), J = 1, Nrad-3)
         ENDDO
        ENDDO   

        Close(21)
        Close(22)
        
        IONCE = 1

        RETURN
      ELSE
        RETURN
      ENDIF

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
