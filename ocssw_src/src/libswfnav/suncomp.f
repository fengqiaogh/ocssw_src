        subroutine suncomp(tlm, time, nlines, gaclac, navqc,
     1     navctl, sun_bod, sunbodfl )
c
c  suncomp(tlm, time, nlines, gaclac, navqc,
c           navctl, sun_bod, sunbodfl )
c
c  Purpose: Create one set of sun vectors from the raw sun sensor data
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  tlm          struct   I      telemetry containing the sun sensor
c                               data
c  time         R*8      I      size nlines time of each line in seconds
c                               past the reference time
c  nlines       I*4      I      number of lines in this segment
c  gaclac       I*4      I      flag for GAC or LAC data type,
c                               1 = GAC, 0 = LAC
c  navqc        struct   I      quality controls including sun angle
c                               and sun angle rate change tolerences
c  navctl       struct   I      navigation controls including sun-sensor
c                               to spacrcraft transform matricies
c  sun_bod      R*4      O      size 3 by nlines weighted sun vectors 
c                               in S/C frame
c  sunbodfl     I*4      O      quality flags for sun_bod vectors, 0 = good
c
c  By: W. Robinson, GSC, 6 Apr 93
c
c  Notes:  this routine takes the 2 sun angle measurements from each
c       of the 3 sun sensors for each telemetry record, performs
c       quality control on the angles, converts the angles into vectors
c       in the spacecraft coordinate system and merges the 3 sensor
c       sets into one consistent set of sun vectors.
c
c  Modification History:
c
c  Modified logic to set larger smoothing interval for LAC and HRPT data.
c  F. S. Patt, SAIC GSC, December 24, 1997.
c
c  Removed check for end of valid data within 3 frames of end of interval, 
c  since this occasionally caused problems and is no longer needed with the
c  Kalman filter attitude algorithm.  F. S. Patt, SAIC GSC, August 6, 1998

      implicit none
c
#include "tlm_str.fin"
#include "navqc_s.fin"
#include "navctl_s.fin"
c
      type(tlm_struct) :: tlm
      type(navqc_struct) :: navqc
      type(navctl_struct) :: navctl
      integer*4 nlines, gaclac, sunbodfl(maxlin)
      real*8 time(nlines)
      real*4 sun_bod(3,maxlin)
c
#include "sunpr_s.fin"
      type(sunproc_struct) :: sunproc(3)
      integer*4 sunrng(2,3), nper, numact, iret, procrng(2,3)
      integer*4 isens, i, k, nfpts, nskip
      real*4 measlcl(2,maxlin), sun_bod3(3,3,maxlin)
c
c       initially set all output flags
      do isens = 1,3
        do i = 1, nlines
          sunproc(isens)%flag(i) = 1
        end do
      end do
c
c
c       start, determine the active ranges for each sensor
c       (assume 1 quasi-contiguous range)
c
      do isens = 1, 3
         call actrng( tlm%ntlm, tlm%sun(isens)%active, sunrng(1,isens) )
      end do
c
c       for the active ranges, check the 2 angles to be within tolerence
c
      call suntol( tlm%sun, navqc, sunrng )
c
c       check the remaining unflagged angles for consistency
c
      call suncnst( gaclac, navqc, sunrng, tlm%sun )
        
c
c       perform a running 3rd order polynomial fit over the active ranges
c
      nfpts = navctl%nsfpts
      nskip = navctl%nsskip 
      if( gaclac .eq. 1 ) then
         nper = 5
      else
         nper = 1
         nfpts = nfpts*20
         nskip = nskip*20
      end if
c
      do isens = 1, 3
         if( sunrng(1,isens) .gt. 0 ) then
            if (sunrng(1,isens) .le. 3) sunrng(1,isens) = 1
            numact = sunrng(2,isens) - sunrng(1,isens) + 1
c
            procrng(1,isens) = ( sunrng(1,isens) - 1 ) * nper + 1
            procrng(2,isens) = sunrng(2,isens) * nper
c
            call runfit3t( nfpts, nskip, measlcl, 
     1          tlm%sun(isens)%ang(1,sunrng(1,isens)), numact, 
     1          2, tlm%sun(isens)%flag(sunrng(1,isens)), 
     1          nper, time(procrng(1,isens)), 
     1          tlm%sun(isens)%deltim(sunrng(1,isens)), 
     1          sunproc(isens)%ang(1,procrng(1,isens)),
     1          sunproc(isens)%flag(procrng(1,isens)), iret )
c
         else
            procrng(1,isens) = -1
         end if
      end do
c
c
c       the sensor data is now checked and fitted to every data line.
c       Now, convert the data to the spacecraft ref frame
c
c       convert the angles to vectors in the spacecraft frame
c
      call sunvec( procrng, sunproc, navctl, sun_bod3 )

      if (navctl%lvdbug .gt. 2) then 
        do i=1,nlines
          do isens=1,3
            write (67,*) i,(sun_bod3(k,isens,i),k=1,3)
          end do
        end do
      end if
c
c       blend the vectors from the 3 sensors into 1 list of vectors
c
      call sunwgt( sun_bod3, procrng, sunproc, nlines, 
     1     sun_bod, sunbodfl )
c
c       and end
c
      return
      end
