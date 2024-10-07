      SUBROUTINE DATAINTP(LL,LAT,LON,F1,DT1,F2,DT2,IPT,NBAND,RNG,DEF
     .  ,INTORDER,DUMMY,FOUT, func, int_bad,IR,JC)
C**********************************************************************
C  Data interpolation (spatial and temporal) subroutine.       
C  Spatial interpolation is weighted by distance and temporal interpolation
C  is a simple linear approximation.       
C  [Only valid surrounding data will be used.  If no valid surrounding value 
C  found, the user specified values will be assigned as output values.]
C------------------------------------------------------------
C  Input parameters:
C  -----------------
C  LL        R*4(2)   Actual pixel location, latitude (-90.0 to 90.0) and
C                     longitude (-180. to 180.)
C  LAT       R*4(x)   Surrounding latitudes (x=IPT=number of data points 
C                     surrounding LL), -90.0 to 90.0.
C  LON       R*4(x)   Surrounding longitudes (x=IPT=number of data points
C                     surrounding LL), -180.0 to 180.0.
C  F1        R*4()    Time step 1 surrounding data values.  For the case of 
C                     each data point has only one value (such as P,SST,...),
C                     the callig program declares F1 as a one dimensional array
C                     , F1(IPT).
C                     For multiple values data (e.g. ozone optical thicknesses),
C                     F1 will be declared as a two dimensional array,
C                     F1(IPT,NBAND). 
C  DT1       R*4      Time difference bt. actual time and time step 1. (DT1>=0.)
C                     [Note: actual time is a time bt. time steps 1 & 2].
C  F2        R*4()    Same as F1 except for time step 2.
C  DT2       R*4      Time difference bt. actual time and time step 2. (DT2>=0.)
C  IPT       I*4      Number of data points surrounding position LL at any
C                     single time step.  (1<=IPT<=25)
C  NBAND     I*4      Band number of each data point.  If each data point has
C                     only one corresponding value, set NBAND=1. (1<=NBAND<=15)
C  RNG       R*4(2)   Data thresholds.  Any surrounding data outside the range
C                     (i.e. data less than or data greater than specified
C                     range) will be excluded for data interpolation.
C  DEF       R*4()    User specified default value(s) which will only be used
C                     if no data is valid for interpolation.  If NBAND=1 the
C                     calling program declares DEF as a scalar.  Otherwise,
C                     DEF has the size of NBAND, DEF(NBAND).
C  INTORDER  I*4      Interpolation order.  Valid sequences are:
C                     112: space then time;    
C                     121: time then space;       
C                     110: spatial interpolation only;
C                     101: temporal interpolation only.
C                     [Note: If all data are good, results from 112 and 121 
C                     are identical] 
C  DUMMY     R*4()    A dummy array which has the size as that of F1.
C  IR,JC     I*4      Sizes, Row and column, of output parameter FOUT.
C                     1)If you expect the output as a single value, set
C                       IR=1, JC=1 and FOUT as a scalar;
C                     2)Output data contain NBAND values, set
C                       IR=1, JC=NBAND and FOUT(NBAND);
C                     3)Output data contain IPT values (e.g. INTORDER=101), set
C                       IR=IPT, JC=1 and FOUT(IPT);  
C                     4)Results contain IPT & NBAND values (e.g. INTORDER=101), 
C                       IR=IPT, JC=NBAND and FOUT(IPT,NBAND).
C  Output parameters:
C  -----------------
C  FOUT      R*4()    Interpolated results for actual location and/or time from
C                     surrounding points. (Reference: parameters IR,JC)
C  func      R*4()    Uncertainty estimate - for 112 intorder, the difference 
C                     between the 2 time values
C  int_bad   I*4      flag to indicate that the space interpolation used
C                     the default or one or both of the inputs to the time 
C                     interpolation was a default - cause for ATMWARN
C                     0 if good, 1 if bad
C------------------------------------------------------------
C
C  Created by Eueng-nan Yeh, 11/17/92, GSC
C  W. Robinson, SAIC 6 Dec 2013  add the int_bad argument
C  W. Robinson, SAIC 4 Aug 2016  add an uncertainty estimate and pass back
C
C**********************************************************************
C%%

      INTEGER*4 MBAND
      PARAMETER (MBAND=15)
      INTEGER*4 NBAND,IPT,INTORDER,IR,JC, int_bad, int_bad1, int_bad2
      REAL*4    LL(1),LAT(1),LON(1),F1(IPT,1),F2(IPT,1),RNG(*),DEF(1)
     .         ,FOUT(1),FOUT1(MBAND),FOUT2(MBAND),DUMMY(IPT,1)
     .         ,MIMX(2)
      REAL*4    func(1)
      REAL*8    DT1,DT2

      int_bad = 0
      IF(RNG(2).GT.RNG(1))THEN
        MIMX(2)=RNG(2)       ! maximum      
        MIMX(1)=RNG(1)       ! minimum
      ELSE
        MIMX(1)=RNG(2)
        MIMX(2)=RNG(1)
      END IF 
   
      IF(INTORDER .EQ. 112)THEN            ! 112
        CALL SPACEINT(LL,LAT,LON,F1,IPT,NBAND,MIMX,DEF,FOUT1, int_bad1 )
        CALL SPACEINT(LL,LAT,LON,F2,IPT,NBAND,MIMX,DEF,FOUT2, int_bad2 )
        CALL TIMEINT(FOUT1,FOUT2,DT1,DT2,1,NBAND,MIMX,DEF,FOUT,func,IR,JC)
        IF( ( int_bad1 .EQ. 1 ) .OR. ( int_bad2 .EQ. 1 ) ) int_bad = 1
      ELSE
        ntot = NBAND * IR
        DO 10 i = 1, ntot
          func(i) = 0.
   10     CONTINUE
        IF(INTORDER.EQ.110)THEN            ! 110
          CALL SPACEINT(LL,LAT,LON,F1,IPT,NBAND,MIMX,DEF,FOUT, int_bad )
        ELSE
          if(INTORDER.EQ.101)then          ! 101
            CALL TIMEINT(F1,F2,DT1,DT2,IPT,NBAND,MIMX,DEF,FOUT,func,IR,JC)
          else
            IF(INTORDER.EQ.121)THEN        ! 121
              CALL TIMEINT(F1,F2,DT1,DT2,IPT,NBAND,MIMX,DEF,DUMMY,func
     .          ,IPT,NBAND)
              CALL SPACEINT(LL,LAT,LON,DUMMY,IPT,NBAND,MIMX,DEF,FOUT, int_bad )
            END IF
          end if
        END IF
      END IF
      RETURN
      END

      SUBROUTINE TIMEINT(F1,F2,DT1,DT2,IPT,NBAND,MIMX,DEF,FOUT,func
     .   ,IR,JC)

      INTEGER*4 NBAND,IPT,IR,JC,N,I,K,IJ
      REAL*4    F1(IPT,1),F2(IPT,1),MIMX(*)
     .         ,DEF(1),FOUT(1),func(1)
      REAL*8    DT1,DT2,T1,T2,W1,W2
C:  FOUT(1)=FOUT(# data point, # of band)   
      T1=ABS(DT1)
      T2=ABS(DT2)
      IF(T1.EQ.0. .AND. T2.EQ.0.) T2=1
C: Interpolating data between T1 and T2
      W1=T2/(T1+T2)
      W2=T1/(T1+T2)
      IF (W1.GE.1.) W2=0.
      IF (W2.GE.1.) W1=0.
      DO 200 N=1,NBAND
        IJ=(N-1)*IR
        DO 201 I=1,IPT
          K=IJ+I
           IF(F1(I,N).GE.MIMX(1).AND.F1(I,N).LE.MIMX(2))THEN  ! F1 is OK
            if(F2(I,N).GE.MIMX(1).AND.F2(I,N).LE.MIMX(2))then ! F1 OK & F2 OK
              FOUT(K)=W1*F1(I,N)+W2*F2(I,N)
              func(K)=abs(F1(I,N)-F2(I,N))
            else
              FOUT(K)=F1(I,N)           ! F1 OK & F2 NG
              func(K)=0
            end if
          ELSE                     ! F1 is NG
            if(F2(I,N).GE.MIMX(1).AND.F2(I,N).LE.MIMX(2))then
              FOUT(K)=F2(I,N)           ! F1 NG & F2 OK
              func(K)=0
            else
              FOUT(K)=DEF(N)            ! F1 NG & F2 NG
              func(K)=0
            end if
          END IF
201     CONTINUE
200   CONTINUE    

9000  CONTINUE
      RETURN
      END
c
c
c
      SUBROUTINE SPACEINT(LL,LAT,LON,F,IPT,NBAND,MIMX,DEF,FOUT, int_bad )
C  Rectangular bi-linear interpolation (IPT==4).
C  by Eueng-nan Yeh           GSC/SAIC                 11/17/95
c
      INTEGER*4 NBAND,IPT,I,N,NNG, int_bad
      INTEGER*4 NCC,NC(4),N1(3,4),N2(5,6),N3(2,4)
      REAL*4    LL(*),LAT(*),LON(*),F(IPT,1),MIMX(*),DEF(1)
     $         ,FOUT(1),LLFT(2),URHT(2),DX,DY,DD,XP,YP,A(4),G(4)
      DATA NC/1,10,100,1000/
      DATA N1/1,2,4,10,1,3,100,2,4,1000,1,3/
      DATA N2/11,1,2,4,3,101,1,3,4,2,1001,1,4,2,3,110,2,3,1,4
     $     ,1010,2,4,3,1,1100,3,4,2,1/
      DATA N3/111,4,1110,1,1101,2,1011,3/
      LLFT(1)=LAT(1)
      LLFT(2)=LON(1)
      URHT(1)=LAT(3)
      URHT(2)=LON(3)
C: Rectangular bi-linear interpolation.
      XP=LL(2)-LLFT(2)
      IF(abs(XP) .gt. 180.) XP=sign(360-abs(XP), XP)
      YP=LL(1)-LLFT(1)
      DX=URHT(2)-LLFT(2)
      IF(abs(DX) .gt. 180.) DX=sign(360-abs(DX), DX)
      DY=URHT(1)-LLFT(1)
      DD=DX*DY
      DO 110 N=1,NBAND
        NCC=0
        NNG=0
        int_bad = 0
        DO 115 I=1,IPT   !IPT==4
          G(I)=F(I,N)    !assign F to local variable G
          IF(G(I).LT.MIMX(1) .OR. G(I).GT.MIMX(2))THEN  !F is ng
            NCC=NCC+NC(I)
            NNG=NNG+1
          ENDIF
115     CONTINUE
        IF(NNG .eq. 0)GOTO 9000

        IF(NNG .eq. 2)THEN    !Two bad points
          DO I=1,6   ! WDR replace IPT with 6 to correct
            IF(NCC .eq. N2(1,I))THEN
              G(N2(2,I))=G(N2(4,I))
              G(N2(3,I))=G(N2(5,I))
              GOTO 9000
            ENDIF
          ENDDO
        ELSE
          IF(NNG .eq.1)THEN   !One bad point
            DO I=1,IPT
              IF(NCC .eq. N1(1,I))THEN
                G(I)=(G(N1(2,I))+G(N1(3,I)))/2.
                GOTO 9000
              ENDIF
            ENDDO
          ENDIF
        ENDIF
        IF(NNG .eq. 4)THEN    !Four bad points
          FOUT(N)=DEF(N)
          int_bad = 1
          GOTO 110
        ELSE
          IF(NNG .eq. 3) THEN !Three bad points
            DO I=1,IPT
              IF(NCC .eq. N3(1,I))THEN
                FOUT(N)=G(N3(2,I))
                GOTO 110
              ENDIF
            ENDDO
          ENDIF
        ENDIF
9000    CONTINUE
        A(1)=G(1)
        IF (ABS(DX) .le. 1.E-9) THEN
           A(2) = 0
        ELSE
           A(2)=(G(4)-A(1))/DX
        ENDIF 
        IF (ABS(DY) .le. 1.E-9) THEN
           A(3) = 0
        ELSE
           A(3)=(G(2)-A(1))/DY
        ENDIF
        IF (ABS(DD) .le. 1.E-9) THEN
            A(4) = 0
        ELSE    
            A(4)=(A(1)-G(2)+G(3)-G(4))/DD
        ENDIF
        FOUT(N)=A(1)+A(2)*XP+YP*(A(3)+A(4)*XP)
110   CONTINUE
      RETURN
      END
