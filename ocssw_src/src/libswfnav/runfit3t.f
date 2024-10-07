      subroutine runfit3t(nfpts, nskip, measlcl, meas, nmeas, nquant, 
     1     flag, nper, time, timdif, measout, flgout, iret )
c $Header$
c $Log$
c 
c
c
c  runfit3t(nfpts, nskip, measlcl, meas, nmeas, nquant, 
c        flag, nper, time, timdif, measout, flgout, iret )
c
c  Purpose: fit a set of data points taken at a set of times by fitting
c       a cubic polynomial to overlaping sections of the data points
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  nfpts        I*4      I      number of points to use in the polynomial
c                               fit, must be >= 4
c  nskip        I*4      I      skip or number of points of input data
c                               to move to do the next fit. nskip must
c                               be <= nfpts
c  measlcl      R*4      I      work array of size nquant, nfpts
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
c  flgout       I*4      O      flag array for measout: 0- good, 1- bad
c  iret         I*4      O      return code: 0 if all good, -1 if bad 
c                               inputs
c
c  By: W. Robinson, GSC, 3 Apr 93
c
c  Notes:  The routine uses routine fit3t to do the polynomial fit.
c          this routine will just feed it the sub-sections of the
c          input data and smoothly merge the fits in the overlap
c          area.  the overlap is nfpts - nskip.  Also note that if
c          any sub-range contains fewer than 4 good points on which
c          to do the fit, the output flags for that range will be
c          set to bad = 1
c
c  Modification History:
c
c  Disabled check for nskip gt nfpts for processing short scenes
c   F. S. Patt, GSC, March 3, 1994
c
c  Minor "clean-up" of code.  F. S. Patt, January 12, 1995.
c
c  Eliminated unused variable flglcl in call to fit3t and corrected logic
c  in calculation of ngood1 and ngood2.  F.S. Patt, GSC, August 15, 1996.
c
c  Added flglcl back to fit3t arguments and included logic to check flags
c  in returned fitted data.  F.S. Patt, GSC, October 6, 1997.
c
c  Added a check in averaging logic for a current value in the output array.
c  F.S. Patt, December 24, 1997.

      implicit none
#include "nav_cnst.fin"
c
      integer*4 nquant, nmeas, nper
      real*4 meas(nquant,nmeas), measout(nquant,nmeas*nper),
     1     timdif(nmeas), measlcl(nquant,nmeas*nper)
      integer*4 nfpts, nskip, flag(nmeas),
     1     flgout(nmeas*nper), flglcl(maxlin), iret
      real*8 time(nmeas*nper)
c
      integer*4 ilin, iquant, line, ovrlap, ifits, nfits, 
     1     istfit, isttim, istfill, ngood,
     1     ngood1, ngood2, nfp, nskp
      real*4 fracinc
      logical lastdone, do_end, badrange
c
c
c       start, make sure that the input controls are valid
c
      nfp = nfpts
      nskp = nskip

      if (nfp.gt.nmeas) nfp = nmeas
      iret = 0
      if( nfp .lt. 4 ) then
        print *,' Error in setup of runfit3t: # fit points must be > 3'
        iret = -1
        go to 990
      end if
c
c       initialize the output flags as bad 
c
      do ilin = 1, nmeas * nper
        flgout(ilin) = 1
      end do
c
c       find the number of fits that can be fully
c       done over the range nmeas
c
      nfits = ( nmeas - nfp ) / nskp + 1
c
c       most likely, there are some points left over at the end
c       that will have to be fit.  The fit range can not be
c       moved by nskip in this case, so it must be flagged
c       and handled differently
c
      if( ( nfits * nskp + nfp - nskp ) .lt. nmeas ) then
        do_end = .true.
        nfits = nfits + 1
      else
        do_end = .false.
      end if
c
c       compute the overlap and indicate that the previous fit
c       could not be done (as the initialization)
c
      ovrlap = nfp - nskp
      lastdone = .false.
c
c       loop over the fit ranges
c
      do ifits = 1, nfits
c
c         find the start point in the input data to fit for and
c         the corrosponding start point in the output ( and in the
c         time array). if the end condition is set, compute start of
c         fit range so the range ends at last measurement.
c         also, re-set the overlap if the end is reached
c
        if( ifits .eq. nfits  .and.  do_end ) then
          istfit = nmeas - nfp + 1
          ovrlap = ( (nfits-1) * nskp + ( nfp - nskp ) ) - istfit
        else
          istfit = ( ifits - 1 ) * nskp + 1
        end if
        isttim = ( istfit - 1 ) * nper + 1
c
c         check for at least 4 unflagged points in the fit range
c         and 2 good points in each half of the range.  If this
c         is not satisfied, do not do that range
c         also, add check for unflagged endpoints
c
        ngood = 0
        badrange = .false.
        ngood1 = 0
        ngood2 = 0
        do ilin = istfit, istfit + nfp - 1
c
          if( flag(ilin) .eq. 0 ) ngood = ngood + 1
c
          if( ilin .le. ( istfit + (nfp - 1) / 2 )  ) then
            if( flag(ilin) .eq. 0 ) ngood1 = ngood1 + 1
          else
            if( flag(ilin) .eq. 0 ) ngood2 = ngood2 + 1
          end if
c
        end do
c
        if( ( ngood .lt. 4 ) .or. ( ngood1 .lt. 2 ) .or.
c     1      ( flag(istfit) .ne. 0 ) .or. 
c     1      ( flag( istfit + nfp - 1 ) .ne. 0 ) .or.
     1      ( ngood2 .lt. 2 )  ) badrange = .true.
c
        if( badrange ) then
          lastdone = .false.
        else
c
c           code for good fit if 4 or more points exist in the range.
c           fit over the current range
c
          call fit3t(meas(1,istfit), nfp, nquant, flag(istfit), nper,
     1     time(isttim), timdif(istfit), measlcl, flglcl)
c
c           merge the overlap region if any
c
          if( lastdone .and. ovrlap .gt. 0 ) then
c
c             set the weighting as the fraction of the overlap
c
            fracinc = 1 / float( ovrlap * nper )
c
c             re-set the output measurements in the overlap by 
c             weighting with the current measurements
c
            do ilin = 1, ovrlap * nper
              line = ilin - 1 + isttim
              if (flglcl(ilin).eq.0) then
c
c              check to be sure there is already a value to average
c
                 if (flgout(line).eq.0) then
                    do iquant = 1, nquant
                       measout(iquant,line) = measout(iquant,line) *
     1                      ( 1 - fracinc * ilin ) + 
     1                      measlcl(iquant,ilin) * fracinc * ilin
                    end do
                 else
                    do iquant = 1, nquant
                       measout(iquant,line) = measlcl(iquant,ilin)
                    end do
                 end if
                 flgout(line) = 0
              end if
            end do
c
c             set the start index for transfer of the rest of the
c             local measurements (measlcl) to the output (measout)
c
            istfill = ovrlap * nper + 1
          else
c
c             set the start transfer index at 1
c
            istfill = 1
          end if
c
c           note that the fit was done and transfer local measurements 
c           to the output array
c
          lastdone = .true.
c
          do ilin = istfill, nper * nfp
             if (flglcl(ilin).eq.0) then
                line = ilin - 1 + isttim
                do iquant = 1, nquant
                   measout(iquant,line) = measlcl(iquant,ilin)
                end do
                flgout(line) = 0
             end if 
          end do
        end if
      end do
c
c       and end
c
  990 continue
      return
      end
