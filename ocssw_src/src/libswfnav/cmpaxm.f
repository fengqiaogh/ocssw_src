        subroutine cmpaxm( nlines, pos, vel, attxfm )
c
c  cmpaxm( nlines, pos, vel, attxfm )
c
c  Purpose: compute the attitude transform matrix that converts
c           from rotating geocentric to spacecraft frame
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  nlines       I*4      I      Number of lines to process
c  pos          R*4      I      size 3 by nlines GPS orbit position
c                               in km/
c  vel          R*4      I      size 3 by nlines GPS orbit velocity
c                               in km / sec
c  attxfm       R*4      O      size 3 by 3 by nlines attitude transform
c                               matricies 
c
c  By: W. Robinson, GSC, 24 Mar 93
c
c  Notes:  This uses code taken from the original ORIENT program
c          written by Fred Patt
c
c  Modification History:
c
c  W. Robinson, GSC, 17 Apr 93  changed sign of z
c
      implicit none
c
      integer*4 nlines
      real*4 pos(3,nlines), vel(3,nlines), attxfm(3,3,nlines)

      integer*4 i, ilin
      real*4 pm, omf2p, pxy, temp, onm, rd
      real*4 on(3),p(3),v(3),x(3),y(3),z(3)
#include "nav_cnst.fin"

c
c
c       set some computation parameters
c

c       Compute constants for navigation model using Earth radius values
      rd=1.d0/omf2
c
c
      do ilin = 1, nlines

c          Determine local vertical reference axes
c          Uses method of local ellipsoid approximation good to 0.3 arcsecond

         do i=1,3
            v(i) = vel(i,ilin)
            p(i) = pos(i,ilin)
         end do
c
c          modify velocity by earth rotation here WDR 19 Apr 93
c
         v(1) = v(1) - omegae * pos(2,ilin)
         v(2) = v(2) + omegae * pos(1,ilin)

c          Compute X axis as local nadir vector
         pm = sqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3))
         omf2p = (omf2*rem + pm - rem)/pm
         pxy = p(1)*p(1)+p(2)*p(2)
         temp = sqrt(p(3)*p(3) + omf2p*omf2p*pxy)
         x(1) = -omf2p*p(1)/temp
         x(2) = -omf2p*p(2)/temp
         x(3) = -p(3)/temp

c          Compute Z axis along orbit normal
         call crossp(v,x,on)
         onm = sqrt(on(1)*on(1)+on(2)*on(2)+on(3)*on(3))
         z(1) = on(1)/onm
         z(2) = on(2)/onm
         z(3) = on(3)/onm

c          Compute Y axis to complete orthonormal triad
c
         call crossp(z,x,y)  

c          Compute attitude matrix in geocentric frame
c          Store local vertical reference vectors in matrix 
         do i=1,3
            attxfm(1,i,ilin)=x(i)
            attxfm(2,i,ilin)=y(i)
            attxfm(3,i,ilin)=z(i)
         end do

      end do
c
c       and return
c
      return
      end
