c ====================================================================
c
c leapchk - determines if the input year is a leap year
c
c Synopsis:
c
c       result = leapchk( year )
c
c       Parameter       Scope           Type            Range
c       year            input           integer*4       (+/- 2**31)
c       result          output          logical*1       (.true. or .false.)
c
c Description:
c
c       Function returns true if "year" is a leap year, otherwise it
c       returns false. The year is expressed as a four-digit number
c       (e.g., 1992 or 2006).
c
c Programmer:
c
c       Bryan A. Franz, General Sciences Corp., 1/16/96
c
c Modifications:
c
c
c ====================================================================

      logical*1 function leapchk(year)
c
      implicit none

      integer*4 year
c
      if (mod(year,4) .eq. 0) then
          leapchk = .true.
      else
          leapchk = .false.
      endif
c
      return
      end

