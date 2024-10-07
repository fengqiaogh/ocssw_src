      subroutine scpar(nlines, attxfm, att_ang, attangfl,
     1     timref, time, pos, tiltpr, tiltfl, navctl, navblk)

c $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.1/swfnav/scpar.f,v 1.1 1995/01/17 23:02:39 seawifsd Exp seawifsd $
c $Log: scpar.f,v $
c Revision 1.1  1995/01/17 23:02:39  seawifsd
c Initial revision
c
c
c  scpar(nlines, attxfm, att_ang, attangfl, timref, time, navblk)
c
c  Purpose: compute additional spacecraft parameters for path
c           where attitude is provided by the spacecraft
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  nlines       I*4      I      number of lines to process
c  attxfm       R*4      I      size 3 by 3 by nlines attitude transform
c                               matricies
c  att_ang      R*4      I      3 by nlines array of spacecraft yaw,
c                               roll and pitch
c  attangfl     I*4      I      array of flags for the attitude values
c  timref       R*8      I      size 3 reference time of the start of 
c                               the current segment
c  time         R*8      I      size nlines time of each line in seconds 
c                               past the reference time
c  pos          R*4      I      size 3 spacecraft position: x, y, z
c  tiltpr       R*4      I      size nlines processed tilt data
c                               for each line and for two tilt measurements
c  tiltfl       I*4      I      flags for goodness of tilts 0 - good,
c                               1 - bad but tilt assumed = 0.
c  navctl       struct   I      controls for navigation processing
c  navblk       struct   O      navigation block containing finished nav
c                               information
c
c  By: W. Robinson, GSC, 29 Mar 93
c
c  Notes:  
c
c  Modification History:
c
c  Changed declaration of sec to R*8.  F. S. Patt, January 12, 1995.
c
c  Added rs to call to l_sun.  F.S. Patt, GSC, August 15, 1996.
c
c  Modified to compute sen_mat and scan_ell using zero tilt if tilt flags
c  are set.  F. S. Patt, GSC, October 9, 1996.
c
c  Added initial index to some array references for Sun OS compatibility.
c  B. A. Franz, GSC, November 14, 1997.
c

      implicit none
#include "navblk_s.fin"
#include "navctl_s.fin"
#ifdef LINUX
#include "nav_cnst.fin"
#endif
c
      integer*4 nlines
c
#ifdef LINUX
      type(navblk_struct) :: navblk(maxlin)
#else
      type(navblk_struct) :: navblk(nlines)
#endif
      type(navctl_struct) :: navctl
c
      integer*4  attangfl(nlines), tiltfl(nlines), j
      real*4 tiltpr(nlines)
      real*4 attxfm(3,3,nlines), att_ang(3,nlines), pos(3,nlines)
      real*8 timref(3), time(nlines), sec, rs
c
      integer*4 ilin, iyr, iday
c
c
c       start, loop through all the lines
c
      do ilin = 1,nlines
c
c          pull the local verical from the attitude transform matrix
c
         navblk(ilin)%l_vert(1) = attxfm(1,1,ilin)
         navblk(ilin)%l_vert(2) = attxfm(1,2,ilin)
         navblk(ilin)%l_vert(3) = attxfm(1,3,ilin)
c
c          get the sun ref angle relative to the earth
c
         iyr = timref(1)
         iday = timref(2)
         sec = timref(3)
         call ydsadd(iyr, iday, sec, 0, time(ilin) )
         call l_sun(iyr, iday, sec, navblk(ilin)%sun_ref(1), rs )
c
c          set navigation flags for bad angles or tilt
c
         navblk(ilin)%nflag(7) = tiltfl(ilin)
         navblk(ilin)%nflag(5) = attangfl(ilin)
         if ((navblk(ilin)%nflag(5) .eq. 1) .or. 
     *          (navblk(ilin)%nflag(7) .eq. 1)) then
            navblk(ilin)%nflag(1) = 1
         end if
c
c          get the scan ellipse coefficients and sensor xfm matrix
c
         call ellxfm( attxfm(1,1,ilin), att_ang(1,ilin), 
     1      tiltpr(ilin), pos(1,ilin), navctl, 
     1      navblk(ilin)%sen_mat, navblk(ilin)%scan_ell )
        
c
         do j = 1, 3
            navblk(ilin)%att_ang(j) = att_ang(j,ilin)
         end do
c
      end do
c
c       and exit
c
      return
      end
