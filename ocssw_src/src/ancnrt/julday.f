      SUBROUTINE JULDAY (YYMMDD,JDAY)
c*******************************************************************
c*
c*      julday
c*
c*      PURPOSE:  Converts month and date to julian day
c*
c*      INPUT PARAMETER:
c*         yymmdd     integer     yy*10000 + mm*100 + dd
c*      OUTPUT PARAMETER:
c*         jday       integer     converted julian date
c*
c*      SUBPROGRAMS REFERENCED
c*         TAE:
c*         IIS:
c*         M75SUB:
c*         TPIO:
c*         LOCAL:
c*
c*      I/O ():
c*
c*      COMMON AREAS: none
c**
c*      BY:   Jean Liu                  GSC             3/23/1984
c*                      TO UNIX: BD Schieber, SAIC/GSC 93
c********************************************************************
c
      IMPLICIT NONE
      INTEGER YYMMDD, JDAY, MONTH(12), YY, MM, DD, IT
c
      DATA MONTH /0,31,59,90,120,151,181,212,243,273,304,334/
c
c********************************************************************
C
C
      DD = MOD (YYMMDD, 100)
      IT = YYMMDD / 100
      MM = MOD (IT, 100)
      YY = IT / 100
C
      JDAY = MONTH(MM) + DD
      IF (MOD (YY, 4) .EQ. 0) THEN
          IF (MM .GT. 2) JDAY = JDAY + 1
      ENDIF
C
C
C     EXIT
C
C
      RETURN
      END
