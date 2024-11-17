      SUBROUTINE KOZSAK2(IFLAG,GE,RE,AJ2,X,XN,XI,IER)

C $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.1/swfnav/kozsak2.f,v 1.1 1995/01/17 23:02:23 seawifsd Exp seawifsd $
C $Log: kozsak2.f,v $
C Revision 1.1  1995/01/17 23:02:23  seawifsd
C Initial revision
C                                                                        
C  VERSION 4/28/88
C  PURPOSE
C    CONVERTS BETWEEN MEAN AND OSCULATING ORBITAL ELEMENTS
C  INPUT
C    IFLAG   = 1, MEAN TO OSCULATING
C            = 2, OSCULATING TO MEAN
C            = 3, OSCULATING TO MEAN WITH INITIAL GUESS FOR MEAN
C    GE      = PRODUCT OF GRAVITATIONAL CONSTANT * MASS OF
C              PLANET (KM**3/SEC**2)
C    RE      = RADIUS OF PLANET (KM)
C    AJ2     = J2 = -C20
C    X       = AN ARRAY OF 6 ORBITAL ELEMENTS, A, E, I, NODE, W, AND M
C              (KM, RAD)
C    XI      = AN ARRAY OF 6 INITIAL MEAN ELEMENTS
C  OUTPUT ARGUMENTS
C    XN      = AN ARRAY OF 6 ORBITAL ELEMENTS AFTER CONVERSION
C    IER     = ERROR FLAG (SET TO 1 IF OSCULATING-TO-MEAN CONVERSION
C               DOES NOT CONVERGE)
C  CALL SUBROUTINES
C    DELM
C  REFERENCES
C    JPL IOM 312/85.2-927, 23 JANUARY 1985, BY C. UPHOFF
C    JPL EM 312/87-153, 20 APRIL 1987
C  ANALYSIS
C    JOHNNY H. KWOK - JPL
C  PROGRAMMER
C    JOHNNY H. KWOK - JPL
C  MODIFICATIONS
C    ADDED OPTION TO INPUT INITIAL GUESS AT MEAN ELEMENTS FOR 
C    OSCULATING-TO-MEAN CONVERSION AND CONVERGENCE TEST FOR CONVERSION
C    F.S. PATT, GSC, 29 NOVEMBER 1993
C
C    CHANGED MAXIMUM NUMBER OF ITERATIONS TO 200 AND INITIALIZED ERROR
C    RETURN CODE.  F.S. PATT, GSC, OCTOBER 26, 1993
C
C  COMMENTS
C    THIS PROGRAM USES AN ALGORITHM DERIVED BY C. UPHOFF (REF) WHICH
C    USES A COMBINATION OF KOZAI'S AND IZSAK'S THEORY.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(6),XN(6),XI(6),DX(6),XTOL(6)
      LOGICAL DONE
      DATA TPI/6.283185307179586D0/
      DATA ZERO/0.D0/,RTOL/1.D-6/
      DATA XTOL/1.D-5,1.D-8,2*1.D-7,2*1.D-5/
      PI = TPI/2.D0
      P3O2 = 3.D0*TPI/2.D0
C      IMAX=5
C      IF (X(2).LT.1.D-1) IMAX=10
C      IF (X(2).LT.1.D-2) IMAX=20
C      IF (X(2).LT.1.D-3) IMAX=30
      IMAX = 200
      IER = 0
      IF (IFLAG.NE.3) THEN
        DO I=1,6
          XI(I)=X(I)
        END DO
      END IF
      ICOUNT=0
c       print *,' '
c       print *,x
      IF (IFLAG.GT.1) THEN
        DONE = .FALSE.
        DO WHILE ((.NOT.DONE).AND.(ICOUNT.LT.IMAX))
          CALL DELM(RE,GE,AJ2,XI,DX)
          DO I=3,6
            DX(I)=DMOD(DX(I)+P3O2,TPI)-PI
          END DO
c         print *,dx
          ICOUNT = ICOUNT + 1
          FAC = 0.5*(1.D0-ABS(DX(5))/PI)
          DONE = .TRUE.
          DO I=1,6
            XDEL = (X(I)-DX(I)-XI(I))
            IF (I.GE.3) XDEL = DMOD(XDEL+P3O2,TPI)-PI
            RDEL = XDEL/DX(I)
            IF ((ABS(XDEL).GT.XTOL(I)).AND.(ABS(RDEL).GT.RTOL)) 
     *        DONE = .FALSE.
            XN(I)=XI(I) + XDEL*FAC
          END DO
          IF (XN(2).LT.ZERO) XN(2)=ZERO
c      DO 30 I=3,6
c   30 XN(I)=DMOD(XN(I)+TPI,TPI)
          DO I=1,6
            XI(I)=XN(I)
          END DO
        ENDDO
        print *,'KOZSAK2:  Mean elements converged in ',icount,
     *    ' iterations'
        IF (ICOUNT.EQ.IMAX) IER = 1
      ELSE
        CALL DELM(RE,GE,AJ2,XI,DX)
        DO I=1,6
          XN(I)=XI(I)+DX(I)
        ENDDO
      ENDIF
  900 CONTINUE
      DO I=3,6
        XN(I)=DMOD(XN(I)+TPI,TPI)
      END DO
      RETURN
      END
