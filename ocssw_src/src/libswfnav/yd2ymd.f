c -------------------------------------------------------------
c Subroutine yd2ymd
c
c  Computes Month and Day of Month from Year and Day of Year
c
c  year                 I       I*4     4-digit year
c  dayOfYear            O       I*4     day of year
c  month                I       I*4     month number (1 - 12)
c  dayOfMonth           I       I*4     day of month
c
c BA Franz, GSC, 12/96
c -------------------------------------------------------------

      subroutine yd2ymd(year,dayOfYear,month,dayOfMonth)
c
      implicit none
c
      save endOfMonth
c
      integer*4 year
      integer*4 month
      integer*4 dayOfMonth
      integer*4 dayOfYear
      integer*4 leap
c
      integer*4 endOfMonth(12,2)
c
      data endOfMonth / 31,59,90,120,151,181,212,243,273,304,334,365,
     .                  31,60,91,121,152,182,213,244,274,305,335,366 /
c
      if (mod(year,4) .eq. 0) then ! Not valid for year 2100
          leap = 2
      else
          leap = 1
      endif
c
      month = 1
      dowhile ( dayOfYear .gt. endOfMonth(month,leap) ) 
          month = month+1
      enddo
c
      if (month .gt. 1) then
          dayOfMonth = dayOfYear - endOfMonth(month-1,leap)
      else
          dayOfMonth = dayOfYear
      endif
c
      return
      end
