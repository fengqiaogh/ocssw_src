      subroutine GREGOR(JULDAY, YEAR, MONTH, DAY)
C********************************************************************
C
C     Converts Julian day into Gregorian MONTH/DAY.
C
C     Input parameters:
C     ----------------
C     JULDAY: I   Julian day of YEAR
C     YEAR:   I   Year in which JULDAY occurs
C
C     Output parameters:
C     -----------------
C     DAY:    I   Day of MONTH in YEAR
C     MONTH:  I   Month in YEAR
C
C     Functions:
C     ---------
C     LEAP:   L     SP$GEN_LIB
C
C     Created by Michael Darzi, GSC, 10/89.
C               To UNIX: B. D. Schieber, SAIC/GAC, 93
C********************************************************************

      implicit none
      logical LEAP
      integer DAY, MONTH, YEAR, JULDAY, JDAYS(12,2), II
      data JDAYS/0,31,59,90,120,151,181,212,243,273,304,334,
     .           0,31,60,91,121,152,182,213,244,274,305,335/

      II = 1
      if (LEAP(YEAR)) II = 2

      MONTH = 12
   10 continue
      if (JULDAY .gt. JDAYS(MONTH,II)) go to 20
      MONTH = MONTH - 1
      go to 10

   20 continue
      DAY = JULDAY - JDAYS(MONTH,II)
      return
      end
