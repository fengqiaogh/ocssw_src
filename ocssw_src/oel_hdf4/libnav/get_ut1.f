      subroutine get_ut1(iyr,iday,ut1c,ier)

c  Routine to retrieve th UTC-UT1 correction from the utcpole.dat file

c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  iyr           I*4     I      Year
c  iday          I*4     I      Day of year
c  ut1c          R*8     O      UT1 - UTC (seconds)
c  ier           I*4     O      Error return code

c  By: Frederick S. Patt, SAIC, July 19, 2011
c
c  Notes: 
c
c  Modification History:
c

      implicit none

      real*8 ut1c, utin
      real*4 tmp1,tmp2,tmp3,tmp4
      integer*4 iyr,iday,ier
      integer*4 jd,jdayt,jdayin/0/
      character*80 utcpolnm,hdr1
      character*108 hdr2
      character*1 tab1,tab2,tab3,tab4,tab5

c  Get path to utcpole.dat file and open file
      utcpolnm = '$OCVARROOT/modis/utcpole.dat'
      call filenv (utcpolnm,utcpolnm)
      print *,utcpolnm
      open (19, file=utcpolnm, status='old', err=999)
      read (19, 1001, err=999) hdr1
 1001 format(a78)
      read (19, 1002, err=999) hdr2
 1002 format(a108)

c  Get truncated Julian day for date

      jdayt = jd(iyr,1,iday) - 2400000
c      print *, jdayt

c  Read file until date is found
      jdayin = 0
      do while (jdayin.ne.jdayt)
        read (19, 1200, end=998) jdayin,tab1,tmp1,tab2,tmp2
     *        ,tab3,tmp3,tab4,tmp4,tab5,utin
 1200   format (i5,2(a1,f9.6,a1,f8.6),a1,f9.6)
      end do
      print *, 'utcpole.dat file record', jdayin, utin
      
      ut1c = utin
c      print *,'UT1C value', ut1c
      close(19)
      return

 998  print *,'End of utcpole file', jdayin, jdayt
      ier = 1
      return

 999  print *,'Error accessing utcpole file'
      ier = 1
      return
      end

