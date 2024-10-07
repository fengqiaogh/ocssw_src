      program interp_nav
c
c
c  Purpose:  This is the main routine for a program which extracts navigation
c  for a LAC scene from the GAC file covering the same time period.  The 
c  attitude angles from the GAC data are interpolated to the LAC scan times 
c  and used to replace the attitude-dependent fields in the LAC file
c
c  Command-line Arguments:
c
c     Name         Type    Description
c     --------     ----    -----------
c     lacfile      char    LAC file name
c     gacfile      char    GAC file name
c
c  By: Frederick S. Patt, SAIC GSC, 24 Sep. 98
c
c  Notes:  Linear interpolation is used for the attitude angles.  
c  The ECEF-to-orbital transformation matrices for the LAC scan lines are 
c  extracted from the original sensor transformation matrices using the 
c  tilt angle, sensor alignment matrix and original attitude angles.
c  The tilt angle and orbit position are assumed to be correct and are 
c  re-used from the LAC data.
c
c  Modification History:
c
c  Modified to stop processing and return error code if LAC data is out of 
c  range of the GAC file provided.  F. S. Patt, SAIC GSC, April 19, 1999.
c
c  Modified to check for GAC file start time on previous day. 
c  F. S. Patt, SAIC GSC, April 22, 1999.
c
c  Modified to exit with non-zero return codes for all error conditions.
c  F. S. Patt, SAIC GSC, May 11, 1999.
c
c  Modified to return number of extrapolated frames > 3 as exit code.
c  F. S. Patt, SAIC GSC, April 27, 2000.
c
c  Fixed bug which caused the same tilt angle to be used for all scan lines.
c  F. S. Patt, SAIC GSC, June 6, 2000.
c
c  Added check for orbit data consistency between LAC and GAC.
c  F. S. Patt, SAIC, September 13, 2002.
c
c  Changed data type for data start day from 2 to 4 bytes, to allow use of 
c  Julian day.  
c  F. S. Patt, SAIC, November 5, 2002.


        implicit none
c
#include "nav_cnst.fin"
#include "navctl_s.fin"
      type(navctl_struct) :: navctl
c
      integer*4 nling, nlinl, iret, amode, prod_ID, sfend, i
      integer*4 msecg(maxlin), msecl(maxlin), msecday, maxret
      integer*4 igday, ilday
      real*4 posg(3,maxlin), posl(3,maxlin), tilt(maxlin)
      real*4 attg(3,maxlin), attl(3,maxlin), coef(6,maxlin)
      real*4 attxfm(3,3), smat(3,3,maxlin)
      character*80 lacfile, gacfile
      data msecday/86400000/

c  Initialize navigation constants
      call cdata
      
c  Read control parameter file
      call readctl(navctl, iret)

      if (iret.ne.0) then
         write (*,*) 'Error reading navctl.dat'
         go to 999
      end if
      
c  Get LAC file name
      call getarg(1,lacfile)
      write( 6, 300 ) lacfile
 300  format( ' interp_nav.f: LAC file =',/,a,/ )

c  Get GAC file name
      call getarg(2,gacfile)
      write( 6, 301 ) gacfile
 301  format( ' interp_nav.f: GAC file =',/,a,/ )

c  Open GAC file 
      amode = 1
      call get_l1a_open( gacfile, amode, prod_ID, igday, nling, iret )

      if (iret.ne.0) then
         write(*,*) 'Error opening GAC file'
         go to 999
      end if

c  Read GAC data
      call get_l1a_data(  prod_ID, nling, msecg, posg, smat, attg, 
     *     tilt, iret )

      if (iret.ne.0) then
         write(*,*) 'Error reading GAC'
         go to 999
      end if

c  Close GAC file
      iret = sfend(prod_ID)
      
c  Open LAC file 
      amode = 3
      call get_l1a_open( lacfile, amode, prod_ID, ilday, nlinl, iret )

      if (iret.ne.0) then
         write(*,*) 'Error opening LAC file'
         go to 999
      end if

c  Read LAC data
      call get_l1a_data( prod_ID, nlinl, msecl, posl, smat, attl, 
     *     tilt, iret )

      if (iret.ne.0) then
         write(*,*) 'Error reading LAC file'
         go to 999
      end if

c  Check for GAC file starting on previous day
      if ( (ilday-igday) .eq. 1 ) then 
         do i=1,nlinl
            msecl(i) = msecl(i) + msecday
         end do
      end if

c  Check for inconsistency between LAC and GAC orbit data
      call checkorb( msecl, posl, msecg, posg, nlinl, nling, iret)
      if (iret.ne.0) then
         write(*,*) 'Orbit data differences in LAC and GAC data'
         go to 999
      end if
      

      maxret = 0

c  Do for each LAC scan line
      do i=1,nlinl

c    Extract the ECEF-to-orbital transformation matrix
         call get_xfm( smat(1,1,i), navctl, tilt(i), attl(1,i),
     *        attxfm)

c    Interpolate the GAC attitude angles to the LAC scan line time 
         call interp_att(msecl(i), nling, msecg, attg, attl(1,i), iret)

         if (iret.ne.0) then
            if (maxret.lt.iret) maxret = iret
c            write(*,*) 'Error interpolating attitude at time',msecl(i)
c            go to 999
         end if

c    Recompute the sensor transformation matrix and the ellipse coefficients
         call ellxfm( attxfm, attl(1,i), tilt(i), posl(1,i), navctl, 
     *       smat(1,1,i), coef(1,i))

      end do
      write(*,*) maxret, ' LAC lines extrapolated'

c  Write recomputed data to LAC file

      call put_l1a_data( prod_ID, nlinl, smat, attl, coef, iret)
      if (iret.ne.0) then
         write(*,*) 'Error writing data to LAC file'
      end if

c  Write recomputed metadata to LAC file

      call put_l1a_metadata( prod_ID, nlinl, posl, smat, coef, iret)
      if (iret.ne.0) then
         write(*,*) 'Error writing metadata to LAC file'
      end if

c  Close LAC file

      iret = sfend(prod_ID)
      stop 

 999  iret = sfend(prod_ID)
      write(*,*) 'interp_nav exiting with error'
      call exit(1)
      
      end






