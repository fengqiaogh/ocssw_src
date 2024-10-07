      integer function JULIAN(DAY, MONTH, YEAR)
C****************************************************************
C
C     Converts Gregorian DAY/MONTH/YEAR into Julian day.
C
C     Input parameters:
C     ----------------
C     DAY:    I   Day of MONTH in YEAR; 1 to 31
C     MONTH:  I   Month in YEAR; 1 to 12
C     YEAR:   I   Year in which to convert DAY/MONTH
C
C     Function:
C     --------
C     LEAP:   L   SP GEN_UTILS
C
C     Created by Michael Darzi, GSC, 9/89.
C               To UNIX: BD Schieber, SAIC/GSC, 1993
C****************************************************************
      implicit none
      logical LEAP
      integer DAY, MONTH, YEAR, JDAYS(12,2), II
      data JDAYS/0,31,59,90,120,151,181,212,243,273,304,334,
     .           0,31,60,91,121,152,182,213,244,274,305,335/

      II = 1
      if (LEAP(YEAR)) II = 2
      JULIAN = min(max(1,DAY),31) + JDAYS(min(max(1,MONTH),12),II)

      return
      end
