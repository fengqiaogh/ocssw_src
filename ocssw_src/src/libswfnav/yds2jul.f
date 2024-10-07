c ----------------------------------------------------------------
c Subroutine yds2jul
c
c Converts to julian time from Year, Day, and Seconds of Day.
c
c BA Franz, GSC, 1/97
c ----------------------------------------------------------------
      subroutine yds2jul(year,day,sec,jul)
c
      implicit none
c
      real*8    jul
      integer*4 year
      integer*4 day
      real*8    sec
      integer*4 jd
c
      jul = jd(year,1,day) + sec/86400.D0
c
      return
      end     
