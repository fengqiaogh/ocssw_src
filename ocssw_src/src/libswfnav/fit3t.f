      subroutine fit3t(meas, nmeas, nquant, flag, nper, time, timdif,
     1     measout, flgout )
c
c  fit3t(meas, nmeas, nquant, flag, nper, time, timdif, measout, flgout )
c
c  Purpose: fit data taken at certain times to a regular time interval
c           using a least squares fit to a cubic polynomial
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  meas         R*4      I      size nquant by nmeas array of measured
c                               quantitys
c  nmeas        I*4      I      number of measurements in array
c  nquant       I*4      I      number of quantities in the array
c  flag         I*4      I      flag array for meas: 0- good, 1- bad
c  nper         I*4      I      expansion factor from meas to measout
c  time         R*8      I      times to fit to, also the first, 
c                               nper+1th... are the tag times for
c                               the input data
c  timdif       R*4      I      offset from the tag time that the 
c                               measurement was made of the quantities
c  measout      R*4      O      size nquant by nmeas * nper array of
c                               fitted measurements
c  flgout       I*4      O      size nmeas array of output flags
c
c  By: W. Robinson, GSC, 1 Apr 93
c
c  Notes:  
c
c  Modification History:
c
c  Eliminate unused variable flgout.  F.S. Patt, GSC, August 16, 1996.
c
c  Added code to set flgout output variable.  F.S. Patt, GSC, Oct. 3, 1996.
c
c  Modified to compute outputs only within range of valid input points.
c  F. S. Patt, SAIC GSC, Feb. 4, 1999.
c
c  Modified to include a check for poorly conditioned polynomial by means of
c  a check on the inverted matrix.   F. S. Patt, SAIC GSC, Mar. 10, 1999.

      implicit none
#include "nav_cnst.fin"
c
      integer*4 nmeas, nquant, nper
      real*4 meas(nquant,nmeas), measout(nquant,nmeas*nper),
     *     timdif(nmeas)
      integer*4 iret, flag(nmeas), 
     *     flgout(nmeas*nper)
      real*8 time(nmeas*nper)
c
c       note that summ is sum of measurements, sumt is sum of times...
      real*8 sumt, sumt2, sumt3, sumt4, sumt5, sumt6, summ(20), 
     1     summt(20), summt2(20), summt3(20), x(4), m(4,4), a(4)
      real*8 newtim(maxlin), t, tp, tdiff, t1
      real*4 tolmat, tm3
      integer*4 imeas, iquant, nsum
      data tolmat/100./

c
c     initialize output flags
c
      do imeas = 1, nmeas*nper
         flgout(imeas) = 1
      end do
c
c
c       start, create the actual measurement time from the
c       time and the time difference
c
      do imeas = 1, nmeas
        newtim(imeas) = time( (imeas - 1) * nper + 1 ) +
     1       timdif(imeas) - time(1)
      end do
c
c       create components to go into the 4 simultaneous equations
c       which make up the solution to the least squares fit
c
      sumt = 0
      sumt2 = 0
      sumt3 = 0
      sumt4 = 0
      sumt5 = 0
      sumt6 = 0
      nsum = 0
      tdiff = 999.
      do iquant = 1,nquant
        summ(iquant) = 0
        summt(iquant) = 0
        summt2(iquant) = 0
        summt3(iquant) = 0
      end do  
c
      do imeas = 1, nmeas
        if (nsum.gt.0) tdiff = newtim(imeas) - tp
        if( (flag(imeas) .eq. 0) .and. (tdiff .gt. 0.1) ) then
          sumt = sumt + newtim(imeas)
          sumt2 = sumt2 + newtim(imeas) * newtim(imeas)
          sumt3 = sumt3 + newtim(imeas) **3 
          sumt4 = sumt4 + newtim(imeas) **4
          sumt5 = sumt5 + newtim(imeas) **5
          sumt6 = sumt6 + newtim(imeas) **6
c
c           there are nquant solutions that will be done here
c
          do iquant = 1, nquant
            summ(iquant) = summ(iquant) + meas(iquant,imeas)
            summt(iquant) = summt(iquant) + meas(iquant,imeas) * 
     1                      newtim(imeas)
            summt2(iquant) = summt2(iquant) + meas(iquant,imeas) *
     1                      newtim(imeas) * newtim(imeas)
            summt3(iquant) = summt3(iquant) + meas(iquant,imeas) *
     1                      newtim(imeas) **3
          end do
          nsum = nsum + 1
          if (nsum .eq. 1) t1 = newtim(imeas)
          tp = newtim(imeas)
        end if
      end do
c
c     check for minimum number of samples
c
      if (nsum.gt.3) then
c
c       set up the matrix for the computations
c
      m(1,1) = nsum 
      m(2,1) = sumt
      m(3,1) = sumt2 
      m(4,1) = sumt3
      m(1,2) = m(2,1)
      m(2,2) = sumt2
      m(3,2) = sumt3
      m(4,2) = sumt4
      m(1,3) = m(3,1)
      m(2,3) = m(3,2)
      m(3,3) = sumt4 
      m(4,3) = sumt5
      m(1,4) = m(4,1)
      m(2,4) = m(4,2)
      m(3,4) = m(4,3)
      m(4,4) = sumt6
c
c       loop and solve for each quantity
c
      do iquant = 1,nquant
c
c         create the right side of the simultaneous equation 
c
        x(1) = summ(iquant) 
        x(2) = summt(iquant)
        x(3) = summt2(iquant)
        x(4) = summt3(iquant)
c
c         invert and solve for coefficients of cubic polynomial
c         x = m a    find a - the coefficients
c         As m is the inverse, only create the inverse once, 
c         after that, multiply the inverted matrix by the right
c         hand sides
c
        if( iquant .eq. 1 ) then
          call invert(m, x, 4, 4, a, iret )

c     check for poorly conditioned solution
          tm3 = sqrt(m(4,4))*(tp - t1)**3
          if (tm3 .gt. tolmat) iret = -1

        else
          call matvec2( m, 4, x, a )
        end if
c
        if( iret .eq. 0 ) then
c
c           fit the output measurements
c
           do imeas = 1, nmeas * nper
              t = time(imeas) - time(1)
              if ((t .ge. t1) .and. (t .le. tp)) then
                 measout(iquant,imeas) = a(1) + a(2) * t + 
     *                a(3) * t * t + a(4) * t * t * t
                 flgout(imeas) = 0
              end if
           end do
        else
           print *,' An error occured in the inversion process'
        end if
      end do

      else
         print *,' Insufficient samples for smoothing'

      end if
c
c       and end
c
  990 continue
      return
      end
