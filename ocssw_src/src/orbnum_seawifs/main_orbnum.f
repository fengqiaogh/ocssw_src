      integer*4 year
      integer*4 day
      integer*4 msec

      integer*2 iyear
      integer*2 iday
      real*8    dsec

      real*8    usec

      integer*4 orbnum
      integer*4 orbnum2
      integer*4 inorad
      real*8    pdsec

      real*8    yds2unix

      character*69 cline1
      character*69 cline2

      write(*,*) 'Enter year, day'
      read(*,*) year,day

      orbnum2 = 0

      do msec=0,86400*1000,100
          call getorbit(year,day,msec,orbnum,pdsec,inorad,cline1,cline2)
          if (msec .eq. 0) orbnum2 = orbnum

          if (orbnum .ne. orbnum2) then
              orbnum2 = orbnum
              iyear = year
              iday = day
              dsec = msec/1000.D0
              usec = yds2unix(iyear,iday,dsec)
              write(*,10) iyear,iday,dsec,usec,pdsec,orbnum
 10           format(i4,x,i3,x,f7.2,x,f12.1,x,f7.2,x,i7)
              call exit(0)
          endif
      enddo

      end
