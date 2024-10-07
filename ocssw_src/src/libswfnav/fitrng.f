      subroutine fitrng(meas, nmeas, nquant, flag, nper, 
     1     measout, flgout )
c
c  fitrng(meas, nmeas, nquant, flag, nper, measout, flgout )
c
c  Purpose: fit data to a finer sampling over a range
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
c  measout      R*4      O      size nquant by nmeas * nper array of
c                               fitted measurements
c  flgout       I*4      O      flag array for measout: 0- good, 1- bad
c
c  By: W. Robinson, GSC, 25 Mar 93
c
c  Notes:  
c
c  Modification History:
c
c  Eliminated redundant call to fndflg.  F.S. Patt, GSC, August 16, 1996.
c

      implicit none
c
      integer*4 nmeas, nquant, nper
      real*4 meas(nquant,nmeas), measout(nquant,nmeas*nper)
      integer*4 flag(nmeas), flgout(nmeas*nper)
c
      real*4 del(20)
      integer*4 imeas, iquant, i1, i2
      logical end
c
c
c       Use interpolation for now 
c
c       fill the output flag array with bad values
c
      do imeas = 1,nmeas * nper
         flgout(imeas) = 1
      end do
c
c       move the array values from the meas array to 
c       the output array
c
      do imeas = 1,nmeas
         flgout( ( imeas - 1 ) * nper + 1 ) = flag( imeas )
         do iquant = 1, nquant
            measout( iquant, ( imeas - 1 ) * nper + 1 ) = 
     1              meas( iquant, imeas )
         end do
      end do
c
c       extrapolate any lines required at the start of the segment
c
c          use first 2 good measurements to extrapolate or interpolate below
c
      call fndflg(flag, nmeas, 1, i1 )
      call fndflg(flag, nmeas, ( i1 + 1 ), i2 )
c
      if( flag(1) .eq. 1) then
c
         do iquant = 1, nquant
            del(iquant) = ( meas(iquant,i2) - meas(iquant,i1) ) / 
     1            ((i2 - i1) * nper)
         end do
c
         do imeas = 1, (i1 - 1) * nper
            flgout( imeas ) = 0
            do iquant = 1, nquant
               measout(iquant, imeas ) = 
     1                   meas(iquant,i1) - del(iquant) *
     1                  ( (i1 - 1) * nper + 1 - imeas )
            end do
         end do
      end if
c
c       interpolate through the available measurements, start 
c        with first 2 found above
c
      end = .FALSE.
      do while( .not. end )
         call fndflg(flag, nmeas, (i1 + 1), i2 )
         if( i2 .le. 0 ) then
            end = .TRUE.
         else
            do iquant = 1, nquant
               del(iquant) = ( meas(iquant,i2) - meas(iquant,i1) ) / 
     1               ((i2 - i1) * nper)
            end do
c
            if( ( i2 * nper - i1 * nper ) .gt. 1 ) then
c
c                there are spaces to fill in output array
c
               do imeas = (i1 - 1) * nper + 2, (i2 - 1) * nper
                  flgout(imeas) = 0
                  do iquant = 1, nquant
                     measout(iquant,imeas) = 
     1                          meas(iquant,i1) + del(iquant) * 
     1                         ( imeas - ( i1 - 1 ) * nper - 1 )
                  end do
               end do
            end if
c
c             find next pair to interpolate
            i1 = i2
c            call fndflg(flag, nmeas, ( i1 + 1 ), i2 )
c            if( i2 .le. 0 ) end = .TRUE.
         end if
      end do
c
c       extrapolate times to the end of the segment
c
      if( ( flag(nmeas) .eq. 1 )  .or.  ( nper .ne. 1 ) ) then
         do imeas = (i1 - 1 ) * nper + 2, nmeas * nper
            flgout(imeas) = 0
            do iquant = 1, nquant
               measout(iquant,imeas) = meas(iquant,i1) + del(iquant) * 
     1                   ( imeas - ( i1 - 1 ) * nper - 1 )
            end do
         end  do
      end if
c
c       and end
c
  990 continue
      return
      end
