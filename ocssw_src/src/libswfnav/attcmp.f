      subroutine attcmp(nlines, gaclac, tlm, pos, vel, timref, 
     1     time, navctl, navqc, attxfm, navblk, tiltpr, tiltfl)
c
c  attcmp(nlines, gaclac, tlm, pos, vel, timref, time, navctl, 
c         navqc, attxfm, navblk)
c
c  Purpose: Compute attitude from the spacecraft sensor data
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  nlines       I*4      I      number of scan lines to process
c  gaclac       I*4      I      flag for GAC or LAC data type,
c  tlm          struct   I      telemetry block data
c  pos          R*4      I      size 3 by nlines orbit position in km
c  vel          R*4      I      size 3 by nlines orbit velocity in km/sec
c  timref       R*8      I      size 3 reference time at start line
c                               of data: year, day, sec
c  time         R*8      I      array of time in seconds relative to
c                               timref for every scan line
c  navctl       struct   I      controls for processing and sensor, 
c                               offset data
c  navqc        struct   I      navigation quality control parameters
c  attxfm       R*4      I      size 3 by 3 by nlines attitude transform
c                               matricies
c  navblk       struct   O      navigation block containing navigation
c                               information
c  tiltpr       R*4      O      size nlines processed tilt data 
c                               for each line and for two tilt measurements
c  tiltfl       I*4      O      flags for goodness of tilts 0 - good,
c                               1 - bad
c
c  By: W. Robinson, GSC,  1 Apr 93
c
c  Notes:  
c
c  Modification History:
c
c  Added rs to call to l_sun.  F.S. Patt, GSC, August 15, 1996.
c
c  Modified to compute sen_mat and scan_ell using zero tilt if tilt flags 
c  are set.  F. S. Patt, GSC, October 9, 1996.
c
c  Minor modifications for compatibility with Sun OS.  B. A. Franz, GSC,
c  November 14, 1997.
c
c  Added initial index to some array references for Sun OS compatibility.
c  B. A. Franz, GSC, November 14, 1997.
c
c  Added check for large yaw angle (indicating large uncertainty) and set
c  nav flag.  F. S. Patt, SAIC GSC, December 11, 1997.
c  
c  Modified to work with new (Kalman filter) version of attdet.f.
c  F. S. Patt, SAIC GSC, June 2, 1998.
c
c  Modified to add two-pass (forward and backward) filtering.
c  F. S. Patt, SAIC GSC, June 4, 1998.
c
c  Fixed bug which caused incorrect initialization for flagged frames at
c  start of scene.  F. S. Patt, SAIC GSC, June 16, 1998.
c
c  Added test for final yaw covariance greater than a limit, indicating 
c  poor geometry.  F. S. Patt, SAIC GSC, September 24, 1998.
c
c  Fixed bug which caused the last scan line not to be checked for covariance
c  greater that the limit.   F. S. Patt, SAIC GSC, December 10, 1998.
c
c  Fixed to compute nadbodfl from widphfl before writing to file.
c  F. S. Patt, SAIC GSC, February 5, 1999.
c
c  Modified to set nflag(8) (navwarn) for nflag(5) = 1 (high yaw uncertainty)
c  or nflag(7) = 2 (tilt change).  F. S. Patt, SAIC, October 10, 2001.
c
c  Modified to recompute yaw covariance from forward and backward passes 
c  before checking against limit.  F. S. Patt, SAIC, March 21, 2003.
c
c  Modified check on input flags to check for valid orbit vector only, in order
c  to allow attitude determination even with flagged time (e.g., time shift).
c  F. S. Patt, SAIC, October 30, 2007.


      implicit none
#include "tlm_str.fin"
#include "navctl_s.fin"
#include "navqc_s.fin"
#include "navblk_s.fin"
c
      type(tlm_struct) tlm
      type(navctl_struct) navctl
      type(navqc_struct) navqc
      type(navblk_struct) navblk(maxlin)
      type(navblk_struct) navblk1
c
      integer*4 nlines, gaclac, isens
      real*4 pos(3,nlines), vel(3,nlines), attxfm(3,3,nlines)
      real*8 timref(3), time(nlines)
c
      integer*4 ilin, i, l, sunbodfl(maxlin), nadbodfl(maxlin), year,
     1     day, delday, tiltfl(maxlin), widphfl(maxlin,2)
      real*4 sun_bod(3,maxlin), nad_bod(3,maxlin), tiltpr(maxlin),
     1     widphse(2,maxlin,2), yaw_lim, covar(3), delta_t, 
     1     att_prev(3), covpass1(3,maxlin), maxcov
      real*8 sec, rs
      logical first
      data yaw_lim/6.0/, att_prev/3*0.0/, maxcov/1.0e-6/
c
c       start, compute the array of spacecraft sun vectors for each line
c
      call suncomp( tlm, time, nlines, gaclac, navqc, navctl, 
     1     sun_bod, sunbodfl )
c
c       debug printout of sun body vectors
      if (navctl%lvdbug.gt.1) then
        open(unit=21,file='sun_bod.out',status = 'unknown')
        do ilin = 1, nlines
          write(21,1200)(sun_bod(isens,ilin),isens=1,3),
     1     sunbodfl(ilin)
 1200   format(2x,3f12.7,i7)
        end do
        close(unit = 21)
      end if
c
c       compute the array of spacecraft nadir vectors for each line
c
      call earcomp(tlm, time, pos, vel, nlines, gaclac, navqc,
     1     navctl, nad_bod, nadbodfl, widphse, widphfl )
c
c       prepare the tilt values reported
c
      call tiltcomp( nlines, tlm, timref, time, gaclac, navqc, 
     1  tiltpr, tiltfl )
c
c       (future) compute the array of spacecraft magnetic vectors
c
c      call magcomp(
c
c       Loop through lines and compute the attitude
c       First pass of Kalman filter
c
      first = .true.
      do ilin = 1, nlines
c
c          compute the model sun vector
c          (get actual year, day, sec from ref time and time array)
c
         year = timref(1)
         day = timref(2)
         sec = timref(3)
         delday = 0
         call ydsadd( year, day, sec, delday, time(ilin) )
c
         call l_sun( year, day, sec, navblk(ilin)%sun_ref(1), rs )
c
c          access the local vertical from the attitude transform matrix
c
         navblk(ilin)%l_vert(1) = attxfm(1,1,ilin)
         navblk(ilin)%l_vert(2) = attxfm(1,2,ilin)
         navblk(ilin)%l_vert(3) = attxfm(1,3,ilin)
c
c          determine attitude from the model and measured vectors
c
         if( (widphfl(ilin,1).ne.0) .and. (widphfl(ilin,2).ne.0)) then
           navblk(ilin)%nflag(4) = 1
         end if
c
         if( sunbodfl(ilin) .ne. 0 ) then
           navblk(ilin)%nflag(3) = 1
         end if
c
         if( navblk(ilin)%nflag(2) .eq. 0 ) then
c
           if (first) then
              first = .false.
              delta_t = 0.0
           else
              delta_t = time(ilin) - time(ilin-1)
           end if

           call attdet(navblk(ilin), navctl,  nad_bod(1,ilin),
     1          attxfm(1,1,ilin), pos(1,ilin), vel(1,ilin), 
     1          widphse(1,ilin,1), widphfl(ilin,1), widphse(1,ilin,2), 
     1          widphfl(ilin,2), sun_bod(1,ilin), sunbodfl(ilin), 
     1          att_prev, delta_t, covar)

c  Store attitude for next line and save covariance for second pass
           do i=1,3
              att_prev(i) = navblk(ilin)%att_ang(i)
              covpass1(i,ilin) = covar(i)
           end do

        end if
      end do

c  Check last scan for yaw angle or covariance outside of limit
      if ((covar(1).ge.maxcov).or.
     1     (abs(navblk(nlines)%att_ang(1)).gt.yaw_lim)) then
         navblk(nlines)%nflag(8) = 1
         navblk(nlines)%nflag(5) = 1
      end if

c  Now make a second pass backwards
      do l=1,nlines-1
         ilin = nlines - l

c  Already have sun_ref and flags from first pass

         if( navblk(ilin)%nflag(2) .eq. 0 ) then

            do i=1,3
               navblk1%sun_ref(i) = navblk(ilin)%sun_ref(i)
            end do
            delta_t = time(ilin) - time(ilin+1)

            call attdet(navblk1, navctl,  nad_bod(1,ilin),
     1           attxfm(1,1,ilin), pos(1,ilin), vel(1,ilin), 
     1           widphse(1,ilin,1), widphfl(ilin,1), widphse(1,ilin,2),
     1           widphfl(ilin,2), sun_bod(1,ilin), sunbodfl(ilin), 
     1           att_prev, delta_t, covar)
            
            do i=1,3
               att_prev(i) = navblk1%att_ang(i)
            end do

c  Compute weighted average of attitude from two passes
            do i=1,3
               navblk(ilin)%att_ang(i) = 
     1              ( navblk(ilin)%att_ang(i)*covar(i) + 
     1              navblk1%att_ang(i)*covpass1(i,ilin)) /
     1              (covar(i) + covpass1(i,ilin))
               covpass1(i,ilin) = 1./(1./covpass1(i,ilin) +
     1              1./covar(i))
            end do

c            write (65,*) (covpass1(i,ilin),i=1,3)
                        
c  Check for yaw angle or covariance outside of limit
            if ((covpass1(1,ilin).ge.maxcov).or.
     1           (abs(navblk(ilin)%att_ang(1)).gt.yaw_lim)) then
               navblk(ilin)%nflag(8) = 1
               navblk(ilin)%nflag(5) = 1
            end if

         end if
      end do
c
c       debug printout of nadir  body vectors
      if (navctl%lvdbug.gt.1) then
        open(unit=21,file='nad_bod.out',status = 'unknown')
        do ilin = 1, nlines
           if ((widphfl(ilin,1).ne.0) .and. (widphfl(ilin,2).ne.0)) then
              nadbodfl(ilin) = 1
           else
              nadbodfl(ilin) = 0
           end if
           write(21,1200) (nad_bod(isens,ilin),isens=1,3), 
     1          nadbodfl(ilin)
        end do
        close(unit = 21)
      end if
c
c       add scan ellipse coefficients and tilt flags to the nav block
c
      do ilin = 1, nlines

         navblk(ilin)%nflag(7) = tiltfl(ilin)
         if ( tiltfl(ilin) .eq. 1) then
            navblk(ilin)%nflag(1) = 1
         else if ( tiltfl(ilin) .eq. 2) then
            navblk(ilin)%nflag(8) = 1
         end if

         call ellxfm( attxfm(1,1,ilin), navblk(ilin)%att_ang(1), 
     1     tiltpr(ilin), pos(1,ilin), navctl,  
     1     navblk(ilin)%sen_mat(1,1), navblk(ilin)%scan_ell(1) )

      end do
c
c       and return
c
      return
      end
