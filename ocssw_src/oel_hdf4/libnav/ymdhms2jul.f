c -------------------------------------------------------------
c Subroutine ymdhms2jul
c
c Convers from Year, Month, Day of Month, Hour, Minute, Second
c to Julian time.
c
c BA Franz, GSC, 1/97
c -------------------------------------------------------------
      subroutine ymdhms2jul(year,month,day,hour,minute,sec,jul)
c
      implicit none
c
      integer*4 year
      integer*4 month
      integer*4 day
      integer*4 hour
      integer*4 minute
      real*8    sec
      real*8    jul
      integer*4 jd
c
      jul = jd(year,month,day) 
     .    + (hour*3600.D0 + minute*60.D0 + sec)/86400.D0
c
      return
      end
