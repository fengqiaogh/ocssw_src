      SUBROUTINE ASAPS(LI,MI,IPLT,ORBIN,IYR,IDAY,SEC,NSTP,CDRG,
     *  TVEC,XVEC)

C $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.1/swfnav/asaps.f,v 1.1 1995/01/17 23:02:01 seawifsd Exp seawifsd $
C $Log: asaps.f,v $
C Revision 1.1  1995/01/17 23:02:01  seawifsd
C Initial revision
C                                                                        

C  This is a subroutine version of the program ASAP, obtained from COSMIC. 
C  The conversion was performed to support the filtering of SeaWiFS GPS 
C  orbit data as part of navigation.  An extensive description of ASAP is 
C  provided in the original program documentation (below).   Several input 
C  parameters were included as calling arguments to override the parameters; 
C  read from the file; added input values to specify time range, and output 
C  arrays to pass back the times and orbit vectors. The format of the input 
C  file was preserved to maintain compatibility with the standalone version.
C
C  Calling Arguments

C  Name         Type    I/O     Description
C
C  LI           I*4      I      Order of gravity field model
C  MI           I*4      I      Degree of gravity field model
C  IPLT         I*4      I      =1, write results to ASCII output files
C  ORBIN(6)     I*4      I      Initial orbital elements (type of elements 
C                                is assumed to be determined by IORB in 
C                                parameter file; 1=mean, 0=osculating
C  IYR          I*4      I      Orbital element epoch year
C  IDAY         I*4      I      Orbital element epoch day-of-year
C  SEC          R*8      I      Orbital element epoch seconds-of-day
C  NSTP         I*4      I      Number of vectors requested (interval is 
C                                specified by STEP in parameter file)
C  CDRG         R*8      I      Atmospheric drag coefficient
C  TVEC(*)      R*8      O      Time tags (offsets from sart of epoch day in seconds)
C  XVEC(6,*)    R*8      O      Orbit vectors
        
C   Conversion performed by Frederick S. Patt, GSC, 20 Dec 93.

C   Modified to add error return code for call to kozsak2.  F.S. Patt, GSC,
C   October 16, 1994.

C   Cleaned up code based on suggestions from Tiger Cheng.  F.S. Patt, GSC,
C   August 15, 1994.

C   Changed TIME common block name to TIMECMN, to eliminate linker 
C   linker warnings.  B. A. Franz, GSC, November 14, 1997.

C      PROGRAM ASAP
C  VERSION 2.03 - 4/18/88
C  PURPOSE   
C    MAIN DRIVER FOR THE ARTIFICIAL SATELLITE ANALYSIS PROGRAM (ASAP).
C    ASAP IS AN EPHEMERIS PROPAGATION PROGRAM FOR ORBITING PLANETARY
C    SPACECRAFTS.  IT USES COWELL'S METHOD OF SPECIAL PERTURBATION.  IT
C    INCLUDES HARMONICS OF UP TO 40X40 FIELD, LUNI-SOLAR GRAVITY, DRAG,
C    AND SOLAR RADIATION PRESSURE.  IT USES A RUNGE-KUTTA 8TH ORDER
C    METHOD WITH STEP SIZE CONTROL FOR NUMERICAL INTEGRATION.  THE
C    PROGRAM IS IN MODULAR FORM SO THAT USERS CAN MODIFY IT TO INCLUDE
C    DIFFERENT I/O MODULES, OR DENSITY MODELS, OR DIFFERENT INTEGRATORS,
C    ETC.  IT ASSUMES A PLANET MEAN  EQUATOR OF EPOCH SYSTEM AND IGNORES
C    POLAR MOTION.  ALL INPUTS ARE EITHER I30 FOR INTEGERS OR D30. FOR
C    DOUBLE PRECISION WITH ONE VALUE TO EACH RECORD.  THE VALUE MUST BE
C    PLACED WITHIN THE FIRST 30 COLUMNS BUT THERE IS NO NEED TO RIGHT
C    JUSTIFY.  COLUMNS 31 TO 80 CAN BE USED FOR COMMENTS.  THE EXCEPTION
C    IS THE INPUT OF THE SPHERICAL HARMONIC COEFFICIENTS.  THIS PROGRAM
C    AND ITS DOCUMENTATION ARE AVAILABLE FROM COSMIC, CAT. #NPO-16731.
C    THERE IS A COMPANION PROGRAM CALLED LONG-TERM PREDICTION (LOP) THAT
C    USES AN AVERAGING METHOD.  LOP IS AVAILABLE FROM COSMIC, CAT.
C    #NPO-17052.
C  INPUT
C    L      = DEGREE OF GRAVITY HARMONICS TO BE INCLUDED
C    M      = ORDER OF GRAVITY HARMONICS TO BE INCLUDED.  MAXIMUM FOR L
C             AND M IS 40.  L CAN BE BIGGER THAN M.
C    IRES   = AN OPTION TO ONLY INCLUDE TESSERALS OF ORDER IRES.  CAN BE
C             USED TO STUDY RESONANCE EFFECTS WHEREBY NON-RESONANT
C             TERRESALS ARE NOT INCLUDED IN THE CALCULATION AND THUS
C             CUTTING PROPAGATION TIME.  DOES NOT AFFECT ZONAL TERMS.
C           = 0, TO TURN OFF THE OPTION.
C    ISUN   = 0, NO SOLAR GRAVITY
C           = 1, WITH SOLAR GRAVITY
C    IMOON  = 0, NO LUNAR GRAVITY
C           = 1, WITH LUNAR GRAVITY
C    IEPHEM = 0, USE TWO-BODY EPHEMERIDES FOR SUN AND MOON
C           = 1, FOR EARTH ORBITING SPACECRAFT ONLY, USE BUILT-IN
C                LUNI-SOLAR EPHEMERIDES.  MUST USE EARTH MEAN EQUATOR
C                AND EQUINOX OF EPOCH SYSTEM
C    IDRAG  = 0, NO DRAG
C           = 1, WITH DRAG
C    IDENS  = 0, USE EXPONENTIAL MODEL
C           = 1, USE STATIC 1977 EARTH MODEL
C    ISRP   = 0, NO SOLAR RADIATION PRESSURE
C           = 1, WITH SOLAR RADIATION PRESSURE
C    IORB     FLAG FOR INPUTING EITHER MEAN OR OSCULATING ORBITAL
C             ELEMENTS.  USES THE VALUE OF C20 TO COMPUTE SHORT PERIOD
C             EFFECTS
C           = 0, INPUT ORBITAL ELEMENTS ORB ARE OSCULATING VALUES
C           = 1, INPUT ORBITAL ELEMENTS ORB ARE MEAN VALUES
C    IPRINT = 0, PRINT AT CONSTANT STEP AS SPECIFIED BY STEP
C           = 1, ALSO PRINT AT PERIAPSIS AND APOAPSIS
C    INODE  = 0, NO NODAL CROSSING PRINT
C           = 1, NODAL CROSSING TIME AND INFO PRINT
C    IPLOT  = 0, NO OUTPUT ASCII FILE FOR PLOTTING
C           = 1, WRITE OUTPUT AT EVERY 'STEP', APSIS AND NODAL CROSSINGS
C                TO FILE 8
C    ORB      OSCULATING OR MEAN ORBITAL ELEMENTS OF THE SPACECRAFT AT
C             TINT.  SEE IORB
C       (1) = A, SEMI-MAJOR AXIS (KM)
C       (2) = E, ECCENTRICITY
C       (3) = I, INCLINATION (DEG)
C       (4) = CAPW, LONGITUDE OF ASCENDING NODE (DEG)
C       (5) = W, ARGUMENT OF PERIAPSIS (DEG)
C       (6) = M, MEAN ANOMALY (DEG)
C    RELERR = RELATIVE ACCURACY OF THE INTEGRATOR.  RECOMMENDED VALUES
C             BETWEEN 1.D-6 TO 1.D-12
C    ABSERR = ABSOLUTE ACCURACY OF THE INTEGRATOR.  RECOMMENDED VALUES
C             BETWEEN 1.D-6 TO 1.D-12
C    STEP   = TIME STEP TO PRINT (SEC).
C    TINT(1)= INITIAL CALENDAR DATE OF RUN (YYYYMMDD.).  FOR EXAMPLE,
C             19880726.D0, FOR 26 JULY 1988.  ALL TIME USED IN THIS
C             PROGRAM ASSUMES EPHEMERIS TIME AND NOT UNIVERSAL TIME.
C        (2)= INITIAL TIME OF DAY OF RUN (HHMMSS.SS...D0).  FOR EXAMPLE,
C             130723.1234D0, FOR 13 HR 7 MIN 23.1234 SEC
C    TFIN   = SAME AS TINT EXCEPT FOR END OF RUN
C    TREF   = SAME AS TINT EXCEPT THIS IS THE TIME CORRESPONDING TO THE
C             POSITION OF THE LUNI-SOLAR EPHEMERIDES (ES AND EM) AND THE
C             PRIME MERIDIAN (PM) INPUT.  IF IEPHEM=1 TO USE BUILT-IN
C             LUNI-SOLAR EPHEMERIDES, THEN THIS IS THE TIME
C             CORRESPONDING TO THE PRIME MERIDIAN INPUT ONLY
C    GE     = PRODUCT OF GRAVITATIONAL CONSTANT AND MASS OF PLANET
C             (KM**3/SEC**2).  RECOMMENDED VALUES (JPL DE118),
C           = 3.9860045D5, FOR EARTH
C           = 3.2485877D5, FOR VENUS
C           = 4.2828287D4, FOR MARS
C    RE     = RADIUS OF PLANET (KM).  RECOMMENDED VALUES (IAU 1982),
C           = 6378.140D0, FOR EARTH
C           = 6051.D0, FOR VENUS
C           = 3393.4D0, FOR MARS
C    RATE   = ROTATION RATE OF THE PLANET (DEG/SEC).  RECOMMENDED VALUES
C             (IAU 1982),
C           = 4.178074216D-3, FOR EARTH
C           = -1.71460706D-5, FOR VENUS
C           = 4.061249803D-3, FOR MARS
C    PM     = LOCATION OF THE PRIME MERIDIAN RELATIVE TO THE INERTIAL
C             X-AXIS AT TREF (DEG).
C    ELLIP  = ELLIPTICITY OF THE REFERENCE ELLIPSOID.  USED BY DRAG
C             TO COMPUTE GEODETIC ALTITUDE FOR ATMOSPHERE DENSITY
C             EVALUATION AND BY OUTPUT ROUTINE TO COMPUTE GEODETIC
C             ALTITUDE.  RECOMMENDED VALUES (IAU 1982),
C           = 0.D0, TO USE A SPHERE
C           = .8182D-1, FOR EARTH
C           = 0.D0, FOR VENUS
C           = .1017D0, FOR MARS
C    RATM     RADIUS OF THE PLANET INCLUDING ATMOSPHERIC BLOCKAGE (KM).
C             IT IS USED TO COMPUTE SHADOW ENTRY AND EXIT IN SOLAR
C             RADIATION PRESSURE COMPUTATION.  RECOMMENDED VALUES,
C           = RE+90 KM FOR VENUS, EARTH, AND MARS
C           = RE FOR MERCURY OR THE MOON
C    RDENS    REFERENCE DENSITY AT REFERENCE HEIGHT RHT TO BE USED BY
C             THE EXPONENTIAL DENSITY MODEL (KG/KM**3)
C    RHT      REFERENCE HEIGHT FOR THE EXPONENTIAL DENSITY MODEL (KM).
C             USE PERIAPSIS ALTITUDE IF POSSIBLE.
C    SHT      SCALE HEIGHT OF THE EXPONENTIAL DENSITY MODEL (KM)
C    ALTMAX   MAXIMUM ALTITUDE TO INCLUDE DRAG PERTURBATION (KM)
C    WT       WEIGHT FACTOR TO BE APPLIED TO THE DENSITY, 1.D0 FOR
C             NOMINAL, 2.D0 FOR TWICE DENSER, ETC.
C    AREAD    EFFECTIVE SPACECRAFT AREA FOR DRAG (KM**2)
C    AREAS    EFFECTIVE SPACECRAFT AREA FOR SOLAR RADIATION PRESSURE
C             (KM**2)
C    SCMASS   EFFECTIVE SPACECRAFT MASS FOR THE LENGTH OF PROPOGATION
C             (KG)
C    CDRAG    DRAG COEFFICIENT, RECOMMENDED VALUES BETWEEN 2.D0 TO 2.2D0
C    CSRP     SOLAR RADIATION PRESSURE CONSTANT (KG/KM-SEC**2).
C             RECOMMENDED VALUES,
C           = G * 4.4D-3, FOR EARTH,
C           = G * 8.4D-3, FOR VENUS,
C           = G * 1.9D-3, FOR MARS,
C             WHERE G < 1 FOR TRANSLUCENT MATERIAL
C                     = 1 FOR BLACK BODY
C                     = 2 FOR PERFECTLY REFLECTIVE MATERIAL
C    GS       PRODUCT OF GRAVITATIONAL CONSTANT AND MASS OF SUN
C             (KM**3/SEC**2).  RECOMMENDED VALUE (JPL DE118),
C           = .13271244D12
C    ES     = ORBITAL ELEMENTS OF THE SUN IN PLANET EQUATOR OF EPOCH,
C             USED IN CALCULATING POINT MASS PERTURBATION DUE TO THE
C             SUN AND SOLAR RADIATION PRESSURE AS WELL.  SEE IEPHEM FOR
C             BUILT-IN LUNI-SOLAR EPHEMERIDES FOR THE EARTH,
C      (1)  = SEMI-MAJOR AXIS (KM)
C      (2)  = ECCENTRICITY
C      (3)  = INCLINATION (DEG)
C      (4)  = LONGITUDE OF THE ASCENDING NODE (DEG), EQUAL TO ZERO IF
C             X-AXIS IS EQUINOX OF EPOCH
C      (5)  = ARGUMENT OF PERIAPSIS (DEG)
C      (6)  = MEAN ANOMALY AT TREF (DEG)
C      (7)  = MEAN MOTION (DEG/SEC)
C    GM       PRODUCT OF GRAVITATIONAL CONSTANT AND MASS OF THE MOON
C             (KM**3/SEC**2).  RECOMMENDED VALUE,
C           = .490279D4, FOR EARTH'S MOON
C    EM       AN ARRAY OF 7 ORBITAL ELEMENTS OF A MOON IN PLANET EQUATOR
C             OF EPOCH SIMILAR TO ES.  SEE IEPHEM FOR BUILT-IN
C             LUNI-SOLAR EPHEMERIDES FOR THE EARTH.
C    N, M, C, S
C             DEGREE, ORDER, OF THE SPHERICAL HARMONIC COEFFICIENTS.
C             CONTINUE TO AS MANY SPHERICAL HARMONICS AS THE FIELD
C             REQUIRES, UP TO 40X40 FIELD.  N SHOULD BE WITHIN COLUMNS
C             1 TO 5, M SHOULD BE WITHIN COLUMNS 6 TO 10, CNM SHOULD BE
C             WITHIN COLUMNS 11 TO 40, AND SNM SHOULD BE WITHIN COLUMNS
C             41 TO 70.  GRAVITY COEFFICIENTS GREATER THAN L AND M ARE
C             NOT COMPUTED IN THE FORCE COMPUTATION.  ONE MAY SET UP A
C             40X40 FIELD HERE ANYWAY EVEN THOUGH A SMALLER FIELD IS
C             RUN.  THE PENALTY IS LONGER DISK READ TIME.
C  CALL SUBROUTINES  
C    JULIAN, KOZSAK, KEPLER, COORD, SETSM, SETTHD, RK78CN, POUT, RK78
C  REFERENCES
C    JPL EM 312/87-153, 20 APRIL 1987
C  ANALYSIS
C    J. H. KWOK - JPL  
C  PROGRAMMER
C    J. H. KWOK - JPL  
C  PROGRAM MODIFICATIONS 
C    NONE
C  COMMENTS
C    NOTE THAT THE SPHERICAL HARMONIC COEFFICIENTS ARE STORED
C    AS C22=C(3,3), ETC.  COEFS. ARE DIMENSIONED FOR 40X40 FIELD.   
C    INPUTS ASSUME KM, SEC, KG, AND DEGREES.
C
C    ALL THE UNIT CONVERSIONS SHOULD BE MADE IN THIS DRIVER
C    PROGRAM SO THAT ALL SUBROUTINES WOULD BE CONSISTENT IN UNITS. 
C   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      SAVE
      DIMENSION TVEC(*),XVEC(6,*),ORBIN(6)
      DIMENSION ORB(6),Y(6),X(6),X1(6)
      DIMENSION TINT(2),TFIN(2),TREF(2)
      CHARACTER*256 FILNM
      LOGICAL FIRST
      COMMON/OPTION/L,M,IRES,ISUN,IMOON,IEPHEM,IDRAG,IDENS,ISRP,IORB
     1  ,IPRINT,INODE,IPLOT
      COMMON/TIMECMN/TI,TF,TR
      COMMON/PLTCON/GE,RE,RATE,PM,AJ2,ELLIP,RATM
      COMMON/ATMCON/RDENS,RHT,SHT,ALTMAX,WT
      COMMON/SPCCON/AREAD,AREAS,SCMASS,CDRAG,CSRP
      COMMON/SUNCON/GS,ES(7),ET(7)
      COMMON/MUNCON/GM,EM(7),EN(7)
      COMMON/HARMON/C(41,41),S(41,41)
      EXTERNAL DER
      DATA NEQ/6/
      DATA DTR,DTS/.1745329251994330D-1,8.64D4/
      DATA TPI/6.283185307179586D0/
      DATA ZERO/0.D0/
      DATA HSTART,HLARGE/60.D0,1.D99/
      DATA FIRST/.TRUE./
      DATA SP/0/
C
C  BEGIN READING INPUT DATA.  THIS BLOCK CAN BE REPLACED
C  IF A NAMELIST ROUTINE IS AVAILABLE
C
      IF (FIRST) THEN
        FIRST = .FALSE.
        FILNM = '$ASAP_PARMS/asap_parms.dat'
        CALL FILENV(FILNM,FILNM)
        OPEN(5,FILE=FILNM,ERR=999)
      END IF

        READ(5,4000)L,M,IRES,ISUN,IMOON,IEPHEM,IDRAG,IDENS,ISRP,IORB
     1    ,IPRINT,INODE,IPLOT
        READ(5,3000)(ORB(I),I=1,6),RELERR,ABSERR,STEP
        READ(5,3000)(TINT(I),I=1,2),(TFIN(I),I=1,2),(TREF(I),I=1,2)
        READ(5,3000)GE,RE,RATE,PM,ELLIP,RATM
        READ(5,3000)RDENS,RHT,SHT,ALTMAX,WT
        READ(5,3000)AREAD,AREAS,SCMASS,CDRAG,CSRP
        READ(5,3000)GS,(ES(I),I=1,7)
        READ(5,3000)GM,(EM(I),I=1,7)
C  
C  BEGIN OUTPUT INPUT DATA
C
C       WRITE(*,4000)L,M,IRES,ISUN,IMOON,IEPHEM,IDRAG,IDENS,ISRP,IORB
C     1   ,IPRINT,INODE,IPLOT
C       WRITE(*,3000)(ORB(I),I=1,6),RELERR,ABSERR,STEP
C       WRITE(*,3000)TINT,TFIN,TREF
C       WRITE(*,3000)GE,RE,RATE,PM,ELLIP,RATM
C       WRITE(*,3000)RDENS,RHT,SHT,ALTMAX,WT
C       WRITE(*,3000)AREAD,AREAS,SCMASS,CDRAG,CSRP
C       WRITE(*,3000)GS,(ES(I),I=1,7)
C       WRITE(*,3000)GM,(EM(I),I=1,7)
C
C  INPUT AND OUTPUT GRAVITY FIELD
C
        DO 10 I=1,41
        DO 10 J=1,41
        C(I,J)=ZERO
   10   S(I,J)=ZERO
        DO 20 N=1,2000
        READ(5,5000,END=30)I,J,C(I+1,J+1),S(I+1,J+1)
   20   CONTINUE
   30   CONTINUE
        REWIND (5)
C       CLOSE (5)
C      END IF

C  Replace values with calling arguments
      L = LI
      M = MI
      IPLOT = IPLT
      CDRAG = CDRG
      DO I=1,6
        ORB(I) = ORBIN(I)
      END DO
      print *,'ASAPS:',l,m,cdrag
      print *,orb
C

C      IF (L.EQ.0.AND.M.EQ.0) GO TO 50
C      DO 40 I=2,L
C      DO 40 J=0,I
C      IF (J.GT.M) GO TO 40
C      WRITE(*,5000)I,J,C(I+1,J+1),S(I+1,J+1)
C   40 CONTINUE
C   50 CONTINUE
      AJ2=-C(3,1)
      IF (IRES.EQ.0) IRES=1
C
C  CONVERT CALENDAR DATE AND TIME TO JULIAN DATE
C
C  TINT and TFIN are not used in subroutine version
C      CALL JULIAN(TINT,TID)
C      CALL JULIAN(TFIN,TFD)
      TID = JD(IYR,1,IDAY) + SEC/DTS - 0.5D0
      TFD = TID + STEP*(NSTP-1)/DTS
      CALL JULIAN(TREF,TRD)
      WRITE(*,6000)TID,TFD,TRD
C
C  CONVERT START CALENDAR DATE TO FLOATING POINT DAY-OF-YEAR
      TDY = IDAY + SEC/DTS
C
C  UNIT CONVERSION FROM DEG TO RADIAN AND DAY TO SECONDS
C
      PM=PM*DTR
      RATE=RATE*DTR
      DO 60 I=3,6
      ORB(I)=ORB(I)*DTR
      ES(I)=ES(I)*DTR
      EM(I)=EM(I)*DTR
   60 CONTINUE
      ES(7)=ES(7)*DTR
      EM(7)=EM(7)*DTR
      TI=TID*DTS
      TF=TI+STEP*(NSTP-1)
C      TF=TFD*DTS
      TR=TRD*DTS
C
C  CONVERT MEAN ELEMENTS TO OSCULATING ELEMENTS IF NECESSARY
C
      IF (IORB.EQ.1) THEN
        CALL KOZSAK2(1,GE,RE,AJ2,ORB,X,X1,IER)
        DO 70 I=1,NEQ
        X1(I) = ORB(I)
   70   ORB(I)=X(I)
      ENDIF
C
C  CHANGE INPUT ORBITAL ELEMENTS TO CARTESIAN COORDINATES
C
      DO 80 I=1,NEQ
   80 Y(I)=ORB(I)
      CALL KEPLER(Y(6),Y(2),EA,SE,CE)
      Y(6)=EA
      CALL COORD(Y,GE,X)
C
C  SET UP ROTATIONAL MATRIX FOR ANALYTICAL EPHEMERIDES FOR THE SUN
C  AND THE MOON.  BESIDES SUN AND SRP, SETSUN IS REQUIRED FOR MOON
C  PERTURBATION BECAUSE THE ECLIPTIC ANGLE IS NEEDED FOR TRANSFORMATION
C  FROM EMO OF DATE TO EME OF DATE
C
      IF (IEPHEM.EQ.1) THEN
        CALL SETSUN(TI,ES)
      ENDIF
      IF (ISUN.EQ.1.OR.ISRP.EQ.1) CALL SETTHD(ES,ET)      
      IF (IMOON.EQ.1.AND.IEPHEM.EQ.0) CALL SETTHD(EM,EN)
C
C  SET UP INTEGRATOR PARAMETERS
C
      H=HSTART
      T=TI
      IVEC = 1
      CALL RK78CN
C
C  SET UP PLOTTING FILE
C
      IF (IPLOT.EQ.1) THEN
        OPEN (8,FILE='8')
        OPEN (9,FILE='9')
        WRITE (9,2000) IYR,TDY
      END IF
C
C  THIS IS SET UP TO PRINT AT FIXED STEP INTERVALS, SEE SP EXPLANATION
C  IF YOU WANT TO DO ADVANCE PROGRAMMING
C
  500 CONTINUE
      IF (IPLOT.EQ.1) CALL POUT(T,X,9)
      TVEC(IVEC) = T - TI + SEC
      DO I=1,6
        XVEC(I,IVEC) = X(I)
      END DO
      TOUT=T+STEP
      IVEC = IVEC + 1   
      CALL RK78(DER,T,TOUT,NEQ,X,H,RELERR,ABSERR,SP)

      IF (TOUT.GE.TF) GO TO 900 
      GO TO 500 
  900 CONTINUE  
      IF (IPLOT.EQ.1) CALL POUT(T,X,9)
      TVEC(IVEC) = T - TI + SEC
      DO I=1,6
        XVEC(I,IVEC) = X(I)
      END DO
C      IF (IORB.EQ.1) THEN
C  CONVERT FINAL VECTOR TO MEAN ELEMENTS FOR OUTPUT
C        CALL EQNOX(X,GE,Y1)
C        X(1) = Y1(1)
C        DO I=2,6
C          X(I) = Y1(I+5)
C        END DO
C       X1(6) = X(5)+X(6)-X1(5)
C        X1(6) = DMOD(X1(6),TPI)
C        CALL KOZSAK2(3,GE,RE,AJ2,X,ORB,X1)
C        DO I=3,6
C          ORB(I) = ORB(I)/DTR
C        END DO
C        WRITE (8,*) ORB
C      END IF
 2000 FORMAT(I10,F14.8)
 3000 FORMAT(BN,D30.16)
 4000 FORMAT(BN,I5)
 5000 FORMAT(BN,2I5,2D30.16)
 6000 FORMAT(1H1,/
     1      ,5X,'RUN STARTS ON JULIAN DATE             = ',D25.16,/
     2      ,5X,'RUN ENDS ON JULIAN DATE               = ',D25.16,/
     3      ,5X,'REFERENCE JULIAN DATE OF PM AND EPHEM = ',D25.16)
      IF (IPLOT.EQ.1) THEN
        OPEN (8,FILE='8')
        OPEN (9,FILE='9')
      END IF
      RETURN
  999 print *,'ERROR OPENING FILE 5'
      RETURN
      END



