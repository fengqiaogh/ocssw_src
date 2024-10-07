        subroutine earcomp(tlm, time, pos, vel, nlines, gaclac, navqc,
     1     navctl, nad_bod, nadbodfl, widphse, widphfl )
c
c  earcomp(tlm, time, pos, vel, nlines, gaclac, navqc,
c           navctl, nad_bod, nadbodfl, widphse, widphfl )
c
c  Purpose: process the earth sensor data into line-by-line satellite 
c           nadir vectors
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  tlm          struct   I      telemetry containing the sun sensor
c                               data
c  time         R*8      I      size nlines time of each line in seconds
c                               past the reference time
c  pos          R*4      I      size 3 by nlines GPS smoothed position 
c                               data for the orbit (km)
c  vel          R*4      I      size 3 by nlines GPS smoothed velocity
c                               data for the orbit (km/sec)
c  nlines       I*4      I      number of lines in this segment
c  gaclac       I*4      I      flag for GAC or LAC data type,
c                               1 = GAC, 0 = LAC
c  navqc        struct   I      quality controls including sun angle
c                               and sun angle rate change tolerences
c  navctl       struct   I      navigation controls including sun-sensor
c                               to spacrcraft transform matricies
c  nad_bod      R*4      O      size 3 by nlines weighted sun vectors 
c                               in S/C frame
c  nadbodfl     I*4      O      quality flags for sun_bod vectors, 0 = good
c  widphse      R*4      O      size 2 by nlines by 2 sensor processed 
c                               width and phase
c  widphfl      I*4      O      size nlines by 2 sensor width and phase flags
c
c  By: W. Robinson, GSC, 13 Apr 93
c
c  Notes:  
c
c  Modification History:
c
c  Modified logic to set larger smoothing interval for LAC and HRPT data.
c  F. S. Patt, GSC, December 24, 1997.
c
      implicit none
c
#include "tlm_str.fin"
#include "navqc_s.fin"
#include "navctl_s.fin"
c
      type(tlm_struct) :: tlm
      type(navqc_struct) :: navqc
      type(navctl_struct) :: navctl
      integer*4 nlines, gaclac, nadbodfl(maxlin), widphfl(maxlin,2)
      real*8 time(nlines)
      real*4 nad_bod(3,maxlin), pos(3,nlines), 
     1     vel(3,nlines), widphse(2,maxlin,2)
c
      integer*4 earrng(2,2), isens, nper, numact, iret, procrng(2,2),
     1     ilin, nfpts, nskip
      real*4 measlcl(2,maxlin)
c
c       initially set all output flags
      do isens = 1, 2
        do ilin = 1, nlines
          widphfl(ilin,isens) = 1
        end do
      end do
c
c       start, determine the active ranges for each sensor
c       (assume 1 quasi-contiguous range)
c
      do isens = 1, 2
         call actrng( tlm%ntlm, tlm%earth(isens)%active, 
     1         earrng(1,isens) )
      end do
c
c       for the active ranges, check the 2 angles to be within tolerence
c
      call eartol( tlm%earth, navqc, earrng )
c
c       check the remaining unflagged angles for consistency
c
      call earcnst( gaclac, navqc, earrng, tlm%earth )
c
c       perform a running 3rd order polynomial fit over the active ranges
c
      nfpts = navctl%nefpts
      nskip = navctl%neskip 
      if( gaclac .eq. 1 ) then
         nper = 5
      else
         nper = 1
         nfpts = nfpts*20
         nskip = nskip*20
      end if
c
      do isens = 1, 2

c       If valid data within 3 frames of start or end of data,
c        extend range to end
         if (earrng(1,isens) .le. 3) earrng(1,isens) = 1
         if (earrng(2,isens) .ge. (tlm%ntlm-3)) 
     *        earrng(2,isens) = tlm%ntlm
         if( earrng(1,isens) .gt. 0 ) then
            numact = earrng(2,isens) - earrng(1,isens) + 1
c
            procrng(1,isens) = ( earrng(1,isens) - 1 ) * nper + 1
            procrng(2,isens) = earrng(2,isens) * nper
c
            call runfit3t( nfpts, nskip, measlcl, 
     1          tlm%earth(isens)%widphse(1,earrng(1,isens)), numact, 2,  
     1          tlm%earth(isens)%flag(earrng(1,isens)), nper, 
     1          time(procrng(1,isens)), 
     1          tlm%earth(isens)%deltim(earrng(1,isens)), 
     1          widphse(1,procrng(1,isens),isens),
     1          widphfl(procrng(1,isens),isens), iret )
c
         else
            procrng(1,isens) = -1
         end if
      end do
c
c
c       the sensor data is now checked and fitted to every data line.
c       Now, create the set of nadir vectors assuming a reference
c       yaw of 0
c
c      yawest = 0.
c      do ilin = 1, nlines
c        call earth( pos(1,ilin), vel(1,ilin), widphse(1,ilin,1), 
c     1     widphfl(ilin,1), widphse(1,ilin,2), widphfl(ilin,2), 
c     1     yawest, navctl, nad_bod(1,ilin), nadbodfl(ilin) )
c      end do
c
c       and end
c
      return
      end
