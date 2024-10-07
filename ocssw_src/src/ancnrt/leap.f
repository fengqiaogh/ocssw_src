c************************  FUNCTION LEAP  ******************************
c
c     returns as "true" if YEAR is a leap year, else returns as "false".
c
c     Created by Michael Darzi, GSC, 8/89
c               To UNIX: BD Schieber, SAIC/GSC, 93
c***********************************************************************
c
      logical function LEAP(YEAR)
      integer*4 YEAR

      if ((mod(YEAR,400) .eq. 0) .or.
     +    ((mod(YEAR,4).eq.0) .and. (mod(YEAR,100).ne.0))) then
        LEAP = .TRUE.
      else
        LEAP = .FALSE.
      endif

      return
      end
   
