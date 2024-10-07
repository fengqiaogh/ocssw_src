        subroutine sunvec(procrng,sunproc,navctl,sun_bod)
c
c  sunvec(procrng,sunproc,navctl,sun_bod)
c
c  Purpose: process sun angles to a sun vector in spacecraft coord system
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  procrng      I*4      I      size 2 by 3 start and end of the active
c                               ranges for the 3 sun sensors
c  sunproc      struct   I      processed sun angle data and flags
c  navctl       struct   I      nav constants, specifically, the transform
c                               matrix from sensor to spacecraft frame
c  sun_bod      R*4      O      size 3 by 3 sensors by # scan lines sun 
c                               vector in spacecraft coord system.  Note
c                               that it is valid only over active range
c
c  By: W. Robinson, GSC, 5 Apr 93
c
c  Notes:  
c
c  Modification History:
c
c  March 30, 1995:  Modified to accept Sun sensor inputs as tangents 
c  instead of angles.  Frederick S. Patt, GSC
c
c  April 23, 1998:  Corrected bug which caused the same bias to be used for
c  both sun sensor angles.  F. S. Patt, SAIC GSC.

      implicit none
c
#include "nav_cnst.fin"
#include "sunpr_s.fin"
#include "navctl_s.fin"
c
      type(sunproc_struct) :: sunproc(3)
      type(navctl_struct) :: navctl
c
      integer*4 procrng(2,3)
      real*4 sun_bod(3,3,*)
c
      integer*4 isens, ilin
      real*4 vec(3), tan1, tan2
c
c
c       start, perform the calculation for each sensor
c
      do isens = 1, 3
c
c          only treat the active range for each sensor
c
         if( procrng(1,isens) .gt. 0 )then
            do ilin = procrng(1,isens), procrng(2,isens)
c
c                compute the sun vector in the sensor frame
c                (Wertz eq 7-23, our boresight is x instead of z)
c
               if( sunproc(isens)%flag(ilin) .eq. 0 ) then
                  tan1 = sunproc(isens)%ang(1,ilin) * 
     *              navctl%sun_scal(1,isens) - navctl%sun_bias(1,isens)
                  tan2 = sunproc(isens)%ang(2,ilin) *
     *              navctl%sun_scal(2,isens) - navctl%sun_bias(2,isens)
c
c                  note that vec(1) is x (2) is y, (3) is z
                  vec(1) = 1. / sqrt( tan1 * tan1 + tan2 * tan2 + 1. )
                  vec(3) = vec(1) * tan2
                  vec(2) = vec(1) * tan1
c
c                   rotate the vector to the spacecraft ref frame
c
                  call matvec( navctl%sun_mat(1,1,isens), vec, 
     1                      sun_bod(1,isens,ilin) )
               end if
            end do
         end if
      end do
c
      return
      end
