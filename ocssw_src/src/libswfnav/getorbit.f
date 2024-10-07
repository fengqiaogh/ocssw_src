      subroutine getorbit(iyr,iday,msec,iorbno,inorad,cline1,cline2)
c
c $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.5/swfnav/getorbit.f,v 1.3 1996/04/24 20:30:48 seawifsd Exp seawifsd $
c $Log: getorbit.f,v $
c Revision 1.3  1996/04/24 20:30:48  seawifsd
c comments out 4 lines using 'd'(instead of 'c') based on the recommendation
c from Bob that the orbit number returned by this function will be consistent
c with the orbit number in the schedule file.
c
c Revision 1.2  1996/01/05 20:12:34  seawifsd
c Corrected element-file search to handle year boundaries
c
c Revision 1.1  1995/01/17 23:02:17  seawifsd
c Initial revision
c
c Revision 1.1  1994/08/30 20:53:29  seawifsd
c Initial revision
c
c Revision 1.1  1994/07/13 19:07:24  seawifst
c Initial revision
c
c Revision 1.3  1994/05/23  20:08:16  seawifst
c New version from Watson. <Mark>
c
c Revision 1.3  1994/04/18  17:10:13  seawifst
c Added $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.5/swfnav/getorbit.f,v 1.3 1996/04/24 20:30:48 seawifsd Exp seawifsd $ $Log: getorbit.f,v $
c Added $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.4/swfnav/getorbit.f,v 1.2 1996/01/05 20:12:34 seawifsd Exp seawifsd $ Revision 1.3  1996/04/24 20:30:48  seawifsd
c Added $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.4/swfnav/getorbit.f,v 1.2 1996/01/05 20:12:34 seawifsd Exp seawifsd $ comments out 4 lines using 'd'(instead of 'c') based on the recommendation
c Added $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.4/swfnav/getorbit.f,v 1.2 1996/01/05 20:12:34 seawifsd Exp seawifsd $ from Bob that the orbit number returned by this function will be consistent
c Added $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.4/swfnav/getorbit.f,v 1.2 1996/01/05 20:12:34 seawifsd Exp seawifsd $ with the orbit number in the schedule file.
c Added $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.4/swfnav/getorbit.f,v 1.2 1996/01/05 20:12:34 seawifsd Exp seawifsd $
c Added $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.5/swfnav/getorbit.f,v 1.3 1996/04/24 20:30:48 seawifsd Exp seawifsd $ Revision 1.2  1996/01/05 20:12:34  seawifsd
c Added $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.5/swfnav/getorbit.f,v 1.3 1996/04/24 20:30:48 seawifsd Exp seawifsd $ Corrected element-file search to handle year boundaries
c Added $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.5/swfnav/getorbit.f,v 1.3 1996/04/24 20:30:48 seawifsd Exp seawifsd $
c Added $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.5/swfnav/getorbit.f,v 1.3 1996/04/24 20:30:48 seawifsd Exp seawifsd $ Revision 1.1  1995/01/17 23:02:17  seawifsd
c Added $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.5/swfnav/getorbit.f,v 1.3 1996/04/24 20:30:48 seawifsd Exp seawifsd $ Initial revision
c Added $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.5/swfnav/getorbit.f,v 1.3 1996/04/24 20:30:48 seawifsd Exp seawifsd $
c Added $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.5/swfnav/getorbit.f,v 1.3 1996/04/24 20:30:48 seawifsd Exp seawifsd $ Revision 1.1  1994/08/30 20:53:29  seawifsd
c Added $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.5/swfnav/getorbit.f,v 1.3 1996/04/24 20:30:48 seawifsd Exp seawifsd $ Initial revision
c Added $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.5/swfnav/getorbit.f,v 1.3 1996/04/24 20:30:48 seawifsd Exp seawifsd $
c Added $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.5/swfnav/getorbit.f,v 1.3 1996/04/24 20:30:48 seawifsd Exp seawifsd $ Revision 1.1  1994/07/13 19:07:24  seawifst
c Added $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.5/swfnav/getorbit.f,v 1.3 1996/04/24 20:30:48 seawifsd Exp seawifsd $ Initial revision
c Added $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.5/swfnav/getorbit.f,v 1.3 1996/04/24 20:30:48 seawifsd Exp seawifsd $
c Revision 1.3  1994/05/23  20:08:16  seawifst
c New version from Watson. <Mark>
c
c  Orbit number fraction is computed by first computing true
c  anomaly instead of mean anomaly, and adding to the argument
c  of perigee. 
c
c Revision 1.4  1994/07/28  20:08:16  seawifst
c New version from Watson. <Mark>
c
c  NORAD Two-line elements are returned from the subroutine.
c  They are as character*69.
c
c
c.......................................................................
c Originally: rdnor.f - 20 August 92
c
c Altered and revised to: rnorad.f
c                         Kenneth S. Lambert
c                         Sea Grant - University of Maryland
c                         15 July 1993
c
c Altered and revised to: getorbit.f
c                         Watson Gregg
c                         Code 902.3
c                         NASA/GSFC
c                         December 1993
c  Modified for better orbit number calculation, using argument
c  of perigee and mean anomaly, 4/18/94, WWG
c
c  Modified to compute and use true anomaly in orbit number
c  calculation, 5/93, Robert Woodward, GSC
c
c  NORAD Two-line elements are returned from the subroutine.
c  They are as character*69, 7/28/94, WWG
c
c  Corrected element-file search to handle year boundaries
c  1/3/96, BAF
c
c  Minor modifications for Sun OS compatibility. removed lenstr and 
c  filenv subroutines, as they are available as separate modules, 
c  add parens around parameter statement, and reversed save and 
c  implicit statements.  B. A. Franz, 11/14/97.
c
c  Given year, day of year, and time (GMT, in milliseconds), computes
c  orbit number for SeaStar/SeaWiFS, for use by the SeaWiFS Data
c  Processing System.  Requires a recent NORAD Two-Line Element
c  File for SeaStar.  I recommend that the element file be no more 
c  than 1 week old, but if the file is updated monthly only a few (1-3)
c  minutes of error should result.  A few minutes is of no 
c  consequence to HRPT stations, since the error occurs at the
c  crossing of the equator on ascending node (always night), but
c  may result in error of lunar calibration files created by the
c  SDPS.
c
c
c  Calling argument list:
c   Input
c    iyr    (I4,1) - year of data
c    iday   (I4,1) - day of year of data
c    msec   (I4,1) - time of day in milliseconds (GMT) of data
c   Output
c    iorbno  (I4,1) - orbit number for the data file
c    inorad  (I4,1) - error return flag:
c                         .eq. 0 - proper return
c                         .ne. 0 - error reading input format
c    cline1  (A69) - first line of NORAD Two-Line Elements
c    cline2  (A69) - second line of NORAD Two-Line Elements
c
c   Local Variables
c    ifound (I4,1) - flag describing whether satellite was found
c    iepyr  (I4,1) - four digit epoch year of element set
c    epday  (R8,1) - epoch day (day of year & fraction: XXX.XXXXXXXX)
c    elem   (R8,6) - element array:
c                        (1) semi-major axis (km)
c                        (2) eccentricity 
c                        (3) inclination (deg)
c                        (4) right ascension of ascending node (deg)
c                        (5) argument of perigee (deg)
c                        (6) mean anomaly (deg)
c    iorbit (i4,1) - orbit number (revolution number) for epoch
c    pi     (R8,1) - Pi
c    radeg  (R8,1) - Radians to degrees conversion factor
c    re     (R8,1) - Earth equatorial radius (km)
c    gm     (R8,1) - Earth gravitational constant (km3/s2)
c    bj2    (R8,1) - 2nd Spherical Harmonic Term
c
c.......................................................................

      implicit real*8 (a-h,o-z)
      save    !required for IRIX compilers 
      parameter (MaxFileSrch = 30)
      logical*1 found
      character*69 cline1,cline2
      character*128 noradf
      real*8 elem(6)
      real*8 rmanom8,reanom8,dkeplr
      data isatnum/24883/
      data imon/1/

c  Constants
      pi = dacos(-1.0D0)
      re = 6378.137D0 
      radeg = 180.0D0/pi
      gm = 398600.5D0
      bj2 = 1082.63D-6
c
c  Variable declarations
      pi2 = 2.0D0*pi                            ! 2*pi
      inorad = 0                                ! error flag

c  Open NORAD Two-Line Element file; search back in time to find 
c   most recent file

      jyr     = iyr
      jday    = iday
      icnt    = 0       ! Count of number of days searched
      found   = .FALSE.

c     !
c     ! Begin search
c     !
      do while ( (icnt .lt. MaxFileSrch) .and. (.not. found) )

          icnt = icnt+1

c         !
c         ! If the file day goes to zero, convert to previous year
c         !
          if (jday .lt. 1) then

              jyr = iyr-1

c             !
c             ! Check for leap year before setting new day number
c             !
              rlp  = jyr/4.0
              ilp  = rlp
              rchk = rlp - ilp

              if (rchk .eq. 0.0) then   
                  jday = 366    ! leap year
              else
                  jday = 365    ! non-leap year
              endif

          endif

c         !
c         ! Try to open the file
c         ! 
          call selfilnm(jyr,jday,noradf)
          open(4,file=noradf,status='old',form='formatted',err=300)
              found = .TRUE.
300       continue

          jday = jday-1

      enddo

c     !
c     ! Verify that a file was found.  Exit on error.
c     !
      if (.not. found) then

          write(6,*)
          write(6,*)'                  ERROR  '
          write(6,*)
          write(6,*)' RECENT NORAD TWO LINE ELEMENT FILE NOT FOUND'
          write(6,*)'               FOR SEASTAR: '
          write(6,*)' File does not exist or most recent file'
          write(6,*)'          exceeds 1 month old'
          write(6,*)'       Obtain updated element file'

          inorad = 1
          return

      endif

c     !
c     ! Search through file and read two line elements. look for satellite no.
c     !
      ieof   = 0
      ifound = 0
      do while ((ieof .eq. 0) .and. (ifound .eq. 0))
          read(4,'(a)',end=100,err=101) cline1
          if (cline1(1:1) .eq. '1') then
              read(4,'(a)',end=100,err=101) cline2
              if (cline2(1:1) .eq. '2') then
                  read(cline1(3:7),'(i5)') isat
                  if (isat .eq. isatnum) ifound = 1
                  read(cline1(10:11),'(i2)') idyr
                  read(cline1(19:32),'(i2,f12.8)') iepyr2, epday
                  read(cline1(34:43),'(f9.8)')rndot2
                  read(cline2(8:63),'(2f9.4,i8,2f9.4)') elem(3),elem(4),
     .                                                  iecc,
     .                                                  elem(5),elem(6)
                  read(cline2(52:63),'(f12.8)') xmm
                  read(cline2(64:68),'(i5)') iorbit
              else
                  backspace(4)
              endif
          endif

      enddo

c  Satellite found
      if (ifound .eq. 1) then

c  Convert integer eccentricity (from assumed decimal) to real.
        elem(2) = dble(iecc)*1.0D-7

c  Compute semi-major axis using method of SGP4
        if (xmm .gt. 1.0D-10) then          !trap for divide by zero
          xno = xmm*pi2/1440.0D0            !mean motion in radians/min
          tothrd = 2.0D0/3.0D0
          xincl = elem(3)/radeg             !inclination
          eo = elem(2)
          ck2 = abs(0.5D0*bj2)               !1/2*j2*ae**2
          xke = sqrt(gm/re**3)*60.0D0     
          a1 = (xke/xno)**tothrd
          cosio = cos(xincl)
          theta2 = cosio*cosio
          x3thm1 = 3.0D0*theta2-1.0D0
          eosq = eo*eo
          betao2 = 1.0D0-eosq
          betao = sqrt(betao2)
          del1 = 1.5D0*ck2*x3thm1/(a1*a1*betao*betao2)
          ao = a1*(1.0D0-del1*(0.5D0*tothrd+del1*
     *      (1.0D0+134.0D0/81.0D0*del1)))
          delo = 1.5D0*ck2*x3thm1/(ao*ao*betao*betao2)
          xnodp = xno/(1.0D0+delo)
          aodp = ao/(1.0D0-delo) !semi-major axis in Earth radii
          elem(1) = aodp*re        !semi-major axis in km
        else
          inorad = 2
        endif
c        write(6,*)'semi-major axis = ',elem(1)
c  Compute period using method of Casey and Way, 1991 (IEEE Trans.
c   Geoscience and Remote Sensing)
        dsma2 = elem(1)*elem(1)
        dsma3 = dsma2*elem(1)
        clsspd = 2.0*pi*sqrt(dsma3/gm)   !classical period
        cosarg = 4.0D0*cosio*cosio - 1.0
        pcorr = 1.0 - 3.0/2.0*bj2*re*re/dsma2*cosarg
        pdsec = clsspd*pcorr
c        write(6,*)'pdsec = ',pdsec
c
c  Output iepyr as four digits, keep within range of 1950 - 2049
        if (iepyr2 .lt. 50) then
          iepyr = 2000 + iepyr2
        else
          iepyr = 1900 + iepyr2
        endif
c
c  Convert to Julian day
c   Compute floating point days since Jan 1.5, 2000
c    Note that the Julian day starts at noon on the specified date
        sec = msec*0.001D0
        tdata = jd(iyr,imon,iday) - 2451545.0D0 
     *          + (sec-43200.D0)/86400.D0
c        write(6,*)'iyr,iday,msec,sec,tdata = ',iyr,iday,msec,sec,tdata
        iepday = epday
        fracday = epday - iepday
        sec = fracday*86400.0D0
        tepoch = jd(iepyr,imon,iepday) - 2451545.0D0 
     *           + (sec-43200.D0)/86400.D0
c  Compute orbit number as function of time difference
c  New as of May, 1994 -- use true anomaly instead of mean anomaly
c        fracorb = (elem(5) + elem(6))/360.0D0
        rmanom8 = elem(6)/radeg   !mean anomaly in radians
        reanom8 = dkeplr(rmanom8,elem(2))  !eccentric anomaly
        a = ((1-elem(2))/(1+elem(2)))**0.5
        rtanom = 2.0 * datan(dtan((reanom8/2.0)/a)) * radeg  !true anom.
        if (rtanom .lt. .0)rtanom = rtanom + 360.0
        fracorb = (elem(5) + rtanom)/360.0D0
c        if (fracorb .ge. 1.0D0)then
c         rorbit = iorbit + (fracorb-1.0D0)
c        else
         rorbit = iorbit + fracorb
c        endif
        diff = (tdata-tepoch)*86400.0D0
        rorb = diff/pdsec + rorbit
        iorbno = rorb
      else !SeaStar satellite  number not found in NORAD TLE file
        print *, 'SeaStar not found in NORAD TLE file...'
        inorad = 1
      endif
      go to 400


100   ieof = 1
      print *, 'SeaStar not found in NORAD TLE file...'
      inorad = 1
      go to 400

101   print *, 'Error reading ',radfilnm
      goto 400

400   close(4)
      return
      end
c
c ******************************************************************
      subroutine selfilnm(myr,idoy,fnme)
c
c  Given year, day of year, returns appropriate file name for
c  daily NORAD Two-Line Element file for SeaStar.
c
      character*128 fnme
      character cday1*1,cday2*2,cday3*3,cday*3,c1*1,cyr*2
      character*10 chead
      character*4 cdat
      character*1 zero
c
c      chead = '/ftp/pub/mission-ops/sel'
      chead = '$NORAD/sel'
      cdat = '.dat'
      zero = '0'
c
c  Convert day
      jday = idoy
      if (jday .lt. 10)then
       write(cday1,'(i1)')jday
       cday = zero//zero//cday1
      else if (jday .lt. 100)then
       write(cday2,'(i2)')jday
       cday = zero//cday2
      else
       write(cday3,'(i3)')jday
       cday = cday3
      endif
c
c  Convert year
      if (myr .le. 1999)then
       itmp = myr - 1900
       write(cyr,'(i2)')itmp
      else if (myr .eq. 2000)then
       cyr = zero//zero
      else if (myr .lt. 2010)then
       itmp = myr - 2000
       write(c1,'(i1)')itmp
       cyr = zero//c1
      else
       itmp = myr - 2000
       write(cyr,'(i2)')itmp
      endif
c
c  File name
      fnme = chead//cday//cyr//cdat
      call filenv(fnme,fnme)
c
      return
      end
c
C ******************************************************************
      FUNCTION DKEPLR(M,E)
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE ! REQUIRED FOR IRIX
      REAL*8 M,PI2/6.283185307179586D0/,TOL/0.5D-15/
C
C SUBROUTINE TO SOLVE KEPLER'S EQUATION
C KEPLER'S EQUATION RELATES GEOMETRY OR POSITION IN ORBIT PLANE
C  TO TIME.
C
C  M - MEAN ANOMALY (O<M<2PI)
C  E - ECCENTRICITY
C  EA - ECCENTRIC ANOMALY
C
      EA=0
      IF(M)1,2,1
    1 EA=M + E*DSIN(M)
      DO 22 I=1,12
      OLDEA=EA
      FE=EA-E*DSIN(EA)-M
      EA=EA-FE/(1-E*DCOS(EA-0.5D0*FE))
C TEST FOR CONVERGENCE
      DELEA=DABS(EA-OLDEA)
      IF(DELEA.LE.TOL)GO TO 2
   22 CONTINUE
    2 EA=DMOD(EA,PI2)
      DKEPLR=EA
      RETURN
      END
C
