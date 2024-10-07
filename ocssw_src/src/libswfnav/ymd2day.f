c ==============================================================
c  Computes Day of Year from Gregorian date
c
c  year                 I       I*4     4-digit year
c  month                I       I*4     month number (1 - 12)
c  dayOfMonth           I       I*4     day of month
c  dayOfYear            O       I*4     day of year
c
c BA Franz, GSC, 12/96
c ==============================================================

      subroutine ymd2day(year,month,dayOfMonth,dayOfYear)
c
      implicit none
c
      save startOfMonth
c
      integer*4 year
      integer*4 month
      integer*4 dayOfMonth
      integer*4 dayOfYear
      integer*4 leap
c
      integer*4 startOfMonth(12,2)
c
      data startOfMonth / 0,31,59,90,120,151,181,212,243,273,304,334,
     .                    0,31,60,91,121,152,182,213,244,274,305,335 /
c
      if (mod(year,4) .eq. 0) then ! Not valid for year 2100
          leap = 2
      else
          leap = 1
      endif
c
      dayOfYear = startOfMonth(month,leap) + dayOfMonth
c
      return
      end
