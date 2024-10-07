      subroutine interp_att( msecl, nlin, msecg, attg, attl, iret)
c
c  interp_att( msecl, nlin, msecg, attg, attl, iret)
c
c  Purpose:  interpolate GAC attitude angles to LAC scan line times
c
c  Calling Arguments:
c     
c     Name         Type    I/O     Description
c     --------     ----    ---     -----------
c     msecl        I*4      I      LAC scan line time
c     nlin         I*4      I      number of GAC scan lines
c     msecg(*)     I*4      I      size nlin array of GAC scan line times
c     attg(3,*)    R*4      I      size 3 by nlin array of GAC attitude angles
c     attl(3)      R*4      O      size 3 array of interpolated attitude angles
c     iret         R*4      O      return code (gt 3, number of scans 
c                                   LAC attitude was extrapolated)
c
c  By: F. S. Patt, SAIC GSC, 24 Sep 98
c
c  Notes:  This routine performs linear interpolation of the attitude angles.
c
c  Modification History:
c
c  Modified to allow LAC times within a specified limit (currently 5 msec) 
c  outside the range of GAC time.  F. S. Patt, SAIC GSC, April 22, 1999
c  
c  Changed time limit from previous update to 505 msec to allow for 3 LAC
c  lines after the last GAC line (due to GAC subsampling).
c  F. S. Patt, SAIC GSC, May 11, 1999
c
c  Modified to remove limit on extrapolation and return number of lines 
c  extrapolated by if more than 3.  F. S. Patt, SAIC GSC, April 27, 2000

      implicit none

      real*4 attg(3,*), attl(3)
      real*4 fac, dt
      integer*4 msecg(*), msecl, nlin, iret
      integer*4 ind, i, limdif, lactim
      data ind/0/, limdif/505/, lactim/166/

      iret = 0

c  Search for first GAC time later than LAC time

      dowhile ((msecl.ge.msecg(ind+1)) .and. (ind.lt.(nlin-1)))
         ind = ind + 1
      end do

c  Check for LAC time after end of GAC time array
      if (ind.ge.nlin-1) then
         if ((msecl-msecg(nlin)).gt.limdif) then
            iret = (msecl-msecg(nlin))/lactim
            write(*,*) 'interp_att:  LAC time after available GAC', msecl
         end if
         ind = nlin - 1
      end if

c  Check for LAC time earlier than first GAC time
      if (ind.le.0) then
         if ((msecg(1)-msecl).gt.limdif) then
            iret = (msecg(1)-msecl)/lactim
            write(*,*) 'interp_att:  LAC time before available GAC',msecl
            return
         end if
         ind = 1
      end if

c  Compute linear interpolation factor
      dt = msecg(ind+1) - msecg(ind)
      fac = (msecl - msecg(ind)) / dt

c  Interpolate attitude angles
      do i=1,3
         attl(i) = attg(i,ind)*(1.0 - fac) + attg(i,ind+1)*fac
      end do

      return
      end
 


