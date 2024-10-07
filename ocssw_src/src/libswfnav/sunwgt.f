        subroutine sunwgt(sun_bod,procrng,sunproc,nlines,wsun,wsunfl)
c
c  sunwgt(sun_bod,procrng,sunproc,nlines,wsun,wsunfl)
c
c  Purpose: blend or weight the sun vector from the 3 sensors
c           into one consistent set
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  sun_bod      R*4      I      size 3 by 3 sensor by nlines sun vector 
c                               data 
c  procrng      I*4      I      size 2 by 3 sensors low and high active 
c                               ranges for each sensor
c  sunproc      struct   I      contains flags for the sensor and line,
c                               0 = good, 1 = bad
c  nlines       I*4      I      number of lines of data
c  wsun         R*4      O      size 3 by nlines weighted sun vector data
c  wsunfl       I*4      O      size nlines flags for the weighted sun vectors
c
c  By: W. Robinson, GSC, 6 Apr 93
c
c  Notes:  
c
c  Modification History:
c
c  March 30, 1995:  Modified weighting scheme to account for Sun sensor
c  data as tangents of angles and to allow for possible swapping of axes.  
c  Frederick S. Patt, GSC
c
c  October 16, 1997.  Modified the maximum sensor value to correspond to 
c  a restricted range of useful sensor data (tangent value of +/- 1.9).
c  Frederick S. Patt, GSC

      implicit none
c
#include "nav_cnst.fin"
#include "sunpr_s.fin"
c
      type(sunproc_struct) :: sunproc(3)
c
      integer*4 procrng(2,3), nlines, wsunfl(nlines)
      real*4 sun_bod(3,3,nlines), wsun(3,nlines), sumwgt
c
      integer*4 ivec, ilin, isens
      real*4 weight,sunmax
      data sunmax/1.9/
c
c
c       start, loop through the lines
c
      do ilin = 1, nlines
c
c         initialize weight, final sun vectors and set final flag to bad
c
        sumwgt = 0
        wsunfl(ilin) = 0
        do ivec = 1, 3
          wsun(ivec,ilin) = 0
        end do
c
c         loop over each sensor
c
        do isens = 1, 3
c
c             compute the weight for this sensor and line as 
c             the position from the end of the sensor range

            if( sunproc(isens)%flag(ilin) .eq. 0 ) then

            weight = sunmax - sqrt(sunproc(isens)%ang(1,ilin)**2 +
     *           sunproc(isens)%ang(2,ilin)**2)
c
c             get the sum of the weights, and weight * vector
c
            sumwgt = sumwgt + weight
            do ivec = 1, 3
              wsun(ivec,ilin) = wsun(ivec,ilin) + 
     1               sun_bod(ivec,isens,ilin) * weight
            end do
            wsunfl(ilin) = 0
c
          end if
        end do
c
c         divide the sum by the weight
c
        if( sumwgt .gt. 0 ) then
          do ivec = 1, 3
            wsun(ivec,ilin) = wsun(ivec,ilin) / sumwgt
          end do
        else
          wsunfl(ilin) = 1
        end if
c
      end do
c
      return
      end
