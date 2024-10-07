c ----------------------------------------------------------------
c Subroutine jul2yds
c
c Converts a julian time to Year, Day, and Seconds of Day.
c
c BA Franz, GSC, 1/97
c ----------------------------------------------------------------
      subroutine jul2yds(jul,year,day,sec)
c
      implicit none
c
      real*8    jul
      integer*4 year
      integer*4 day
      real*8    sec
      integer*4 julday
      integer*4 month
      integer*4 dayOfMonth
c
      julday = int( jul )
      call jddate( julday, year, month, dayOfMonth )
      call ymd2day( year, month, dayOfMonth, day )
      sec = dmod(jul,1.D0) * 86400.D0
c
      return
      end     
