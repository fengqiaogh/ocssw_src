      SUBROUTINE POUT(T,X,LUN)  
C  VERSION OF 4/13/87
C  PURPOSE
C    PRINTOUT ROUTINE FOR TRAJECTORY PROGRAM   
C  INPUT 
C    T      = CURRENT TIME (SEC)
C    X      = 6-D CARTESIAN COORD X,Y,Z,XD,YD,ZD (KM,KM/SEC)
C    SEE ROUTINES ASAP AND DER FOR EXPLANATION OF COMMON BLOCK
C    VARIABLES
C  OUTPUT
C    NONE  
C  CALL SUBROUTINES  
C    EQNOX 
C  REFERENCES
C    JPL EM 312/87-153, 20 APRIL 1987
C  ANALYSIS
C    J. H. KWOK - JPL  
C  PROGRAMMER
C    J. H. KWOK - JPL  
C  PROGRAM MODIFICATIONS 
C    NONE  
C  COMMENTS  
C    NONE
C   
C  MODIFICATION HISTORY
C   Changed TIME common block name to TIMECMN, to eliminate linker 
C   linker warnings.  B. A. Franz, GSC, November 14, 1997.
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DIMENSION X(6)
      DIMENSION Y(18)
      COMMON/OPTION/L,M,IRES,ISUN,IMOON,IEPHEM,IDRAG,IDENS,ISRP,IORB
     1  ,IPRINT,INODE,IPLOT
      COMMON/TIMECMN/TI,TF,TR  
      COMMON/PLTCON/GE,RE,RATE,PM,AJ2,ELLIP,RATM
      COMMON/COUNT/NODE,NREV
      COMMON/CASE/ICASE
      DATA HTS,DTH,RTD,PI,TPI,ONE/3.6D3,24.D0,57.29577951308232D0,
     1  3.141592653589793D0,6.283185307179586D0,1.D0/
      DATA SMALL/1.D-3/
      THOUR=(T-TI)/HTS
      DDAY=THOUR/DTH
      IDAY=IDINT(DDAY)
      HOUR=THOUR-IDAY*DTH
      CALL EQNOX(X,GE,Y)
      HA=DMOD(PM+(T-TR)*RATE,TPI)
      ENODE=DMOD(Y(9)-HA+TPI,TPI)
      ALONG=DATAN2(X(2),X(1))-HA
      ALONG=DMOD(ALONG+TPI+TPI,TPI)
      ALAT=DASIN(X(3)/Y(16))
      ARGL=DMOD(Y(14)-Y(9)+TPI,TPI)
      EALT=Y(16)-DSQRT(RE**2*(ONE-ELLIP**2)/(ONE-(ELLIP*DCOS(ALAT))**2))
      PALT=Y(1)*(ONE-Y(7))-RE
      ALT=Y(16)-RE
C
C *** DOUBLE CHECK APSIS CROSSING, THIS IS NEEDED IN CASE SOME
C *** PERTURBATION ON NEARLY CIRCULAR ORBIT IS CHANGING MEAN ANOMALY
C *** BACKWARDS
C
      IF (IPRINT.EQ.1) THEN
      IF (ICASE.EQ.3.OR.ICASE.EQ.4) THEN
        IF (DABS(Y(11)).LT.SMALL.OR.DABS(Y(11)-TPI).LT.SMALL) ICASE=3
        IF (DABS(Y(11)-PI).LT.SMALL) ICASE=4
      ENDIF
      ENDIF
      IF (IPRINT.EQ.1.OR.INODE.EQ.1) THEN
C
C *** CHANGE ICASE IF FIRST PRINT OUT IS A NODAL OR APSIS CROSSING
C
      IF (NODE.EQ.0.AND.(DABS(ARGL).LT.SMALL.OR.DABS(ARGL-TPI).LT.SMALL)
     1) ICASE=1
      IF (NODE.EQ.0.AND.DABS(ARGL-PI).LT.SMALL) ICASE=2
      IF (NREV.EQ.0.AND.(DABS(Y(11)).LT.SMALL.OR.DABS(Y(11)-TPI).LT.SMAL
     1L)) ICASE=3
      IF (NREV.EQ.0.AND.DABS(Y(11)-PI).LT.SMALL) ICASE=4
      ENDIF
      IF (IPRINT.EQ.1.OR.INODE.EQ.1) THEN
      IF (ICASE.EQ.0) GO TO 105
      GO TO (101,102,103,104) ICASE
  101 NODE=NODE+1
      WRITE(*,1001)NODE
      GO TO 105
  102 WRITE(*,1002)NODE
      GO TO 105
  103 NREV=NREV+1
      WRITE(*,1003)NREV
      GO TO 105
  104 WRITE(*,1004)NREV
  105 CONTINUE
      ENDIF
C
C *** START UNIT CONVERSION
C
      Y(6)=Y(6)*RTD
      DO 10 I=8,15
   10 Y(I)=Y(I)*RTD
      Y(18)=Y(18)/HTS
      ARGL=ARGL*RTD
      HA=HA*RTD
      ENODE=ENODE*RTD
      ALONG=ALONG*RTD
      ALAT=ALAT*RTD
c      WRITE(*,1000)IDAY,HOUR
c      WRITE(*,2000)X
      WRITE(LUN,2001)IDAY,HOUR,X
c      WRITE(*,3000)(Y(I),I=1,6)
c      WRITE(*,4000)(Y(I),I=7,15),ARGL
c      WRITE(*,5000)(Y(I),I=16,17),ENODE,ALONG,ALAT,HA,ALT,EALT,PALT
c     1  ,Y(18)
      IF (IPLOT.EQ.1) THEN
        WRITE(8,9000)ICASE,NODE,NREV,DDAY,THOUR,Y(1)
     1  ,(Y(I),I=7,11),ALONG,ALAT,ENODE,EALT,PALT,ALT
        ICASE=0
      ENDIF
 1001 FORMAT(/5X,'PRINTING INFORMATION AT ASCENDING NODE # ',I5)
 1002 FORMAT(/5X,'PRINTING INFORMATION AT DESCENDING NODE #',I5)
 1003 FORMAT(/5X,'PRINTING INFORMATION AT PERIAPSIS OF REV #',I5)
 1004 FORMAT(/5X,'PRINTING INFORMATION AT APOAPSIS OF REV # ',I5)
 1000 FORMAT(/5X,I5,' DAYS',F12.6,' HOURS FROM EPOCH')
 2001 FORMAT(I3,F11.7,3F11.4,3F11.7)
 2000 FORMAT(1P,5X,'CARTESIAN COORD X, Y, Z, XD, YD, ZD',/,6D22.14)
 3000 FORMAT(1P,5X,'EQUINOCTIAL ELEMENTS A, H, K, P, Q, MEAN LONG',/
     1,6D22.14)
 4000 FORMAT(1P,5X,'CLASSICAL ELEMENTS E, I, NODE, W, MA, TA, EA, TRUE L
     1ONG, ECC LONG, ARG OF LAT',/,5D22.14/5D22.14) 
 5000 FORMAT(1P,5X,'OTHER PARAMETERS R, V, ENODE, LONG, LAT, HOUR ANGLE/
     1 ALT, ELLIPSOIDAL ALT, PERIAPSIS ALT, PERIOD',/,6D22.14/6D22.14)
 9000 FORMAT(3I5,2F12.6,F14.5,F12.8,7F12.5,3F14.5)
      RETURN
      END
