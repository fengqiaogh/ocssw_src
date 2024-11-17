       subroutine vec2mean(vec,ge,aj2,xinit,xmean,ier)

c $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.1/swfnav/vec2mean.f,v 1.1 1995/01/17 23:03:07 seawifsd Exp seawifsd $
c $Log: vec2mean.f,v $
c Revision 1.1  1995/01/17 23:03:07  seawifsd
c Initial revision
c                                                                        

c  Purpose:  This subroutine computes mean elements with respect to J2 
c   from an orbit vector.
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  vec(6)       R*8      I      Input orbit vector
c  ge           R*8      I      Earth gravitational constant
c  aj2          R*8      I      J2 gravity term
c  xinit(6)     R*8      I      Initial guess at mean elements
c  xmean(6)     R*8      O      Output mean elements
c  ier          I*4      O      Error flag (mean elements did not converge)
c
c
c  By: Frederick S. Patt, GSC, October 19, 1994
c
c  Notes:  
c
c  Modification History:
c
      implicit none
#include "nav_cnst.fin"

      real*8 vec(6),xinit(6),xmean(6),ge,aj2
      real*8 x(6),tmp(6),y(18)
      integer*4 i,ier

      call eqnox(vec,ge,y)
      x(1) = y(1)
      do i=2,6
          x(i) = y(i+5)
      end do
      tmp(1) = xinit(1)
      tmp(2) = xinit(2)
      do i=3,6
        tmp(i) = xinit(i)/radeg
      end do 
      tmp(6) = x(5) + x(6) - tmp(5)
      call kozsak2(3,ge,re,aj2,x,xmean,tmp,ier)
      do i=3,6
          xmean(i) = xmean(i)*radeg
      end do
      
      return
      end
