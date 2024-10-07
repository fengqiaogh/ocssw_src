      subroutine checkorb( msecl, posl, msecg, posg, nlinl, nling, iret)
c
c  checkorb(  msecl, posl, msecg, posg, nlinl, nling, iret)
c
c  Purpose:  checks for consistency of LAC and GAC orbit
c
c  Calling Arguments:
c     
c     Name         Type    I/O     Description
c     --------     ----    ---     -----------
c     msecl(*)     I*4      I      LAC scan line times
c     posl(3,*)    R*4      I      LAC orbit position vectors (km)
c     msecg(*)     I*4      I      GAC scan line times
c     posg(3,*)    R*4      I      GAC orbit position vectors (km)
c     nlinl        I*4      I      number of LAC scan lines
c     nling        I*4      I      number of GAC scan lines
c     iret         R*4      O      return code (gt 3, number of scans 
c                                   LAC attitude was extrapolated)
c
c  By: F. S. Patt, SAIC, 13 Sep 2002
c
c  Notes:  
c
c  Modification History:
c
c  Increase tolerance to 0.02 km.  F. S. Patt, SAIC, January 10, 2003.

      implicit none

      real*4 posl(3,*), posg(3,*)
      integer*4 msecl(*), msecg(*), nlinl, nling, iret
      integer*4 ig, il, i
      real*4 toldif

      data toldif /0.02/
    
c     Find first GAC and LAC frame with same time

      ig = 1
      il = 1
      iret = 0

      dowhile (msecl(il) .ne. msecg(ig))

c     If LAC time is greater increment GAC index
         if (msecl(il) .gt. msecg(ig)) then
            ig = ig + 1

c       Check for GAC index out of range
            if (ig .gt. nling) then
               iret = -1
               write(*,*) 'CHECKORB:  No matching GAC and LAC times'
               go to 999
            end if

         else
c     Else increment LAC index
             il = il + 1

c       Check for LAC index out of range
            if (il .gt. nlinl) then
               iret = -1
               write(*,*) 'CHECKORB:  No matching GAC and LAC times'
               go to 999
            end if
         end if
      end do

c     Compare LAC and GAC orbit data

      do i=1,3
         if (abs(posl(i,il) - posg(i,ig)) .gt. toldif) then
            write(*,*)'CHECKORB:  Inconsistent GAC and LAC orbit'
            write(*,*) i, il, posl(i,il), ig, posg(i,ig)
            iret = -1
         end if
      end do
      
 999  return
      end
