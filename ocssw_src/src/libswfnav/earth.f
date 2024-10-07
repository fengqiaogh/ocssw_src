       subroutine earth( pos, vel, widphse1, widphfl1, widphse2,
     1     widphfl2, yawest, navctl, nad_bod, nadbodfl, ssdec)
c
c  earth( pos, vel, widphse1, widphfl1, widphse2,
c       widphfl2, yawest, navctl, nad_bod, nadbodfl)
c
c  Purpose: determine the nadir body vectors from the earth
c           sensor data
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  pos          R*4      I      size 3 orbit Position Vector (km)
c  vel          R*4      I      size 3 orbit Velocity Vector (km/sec)
c  widphse1     R*4      I      size 2 processed earth width and 
c                               phase for sensor 1, in degrees
c  widphfl1     I*4      I      flag for earth sensor 1 - 0 good
c  widphse2     R*4      I      size 2 processed earth width and 
c                               phase for sensor 2, in degrees
c  widphfl2     I*4      I      flag for earth sensor 2 - 0 good
c  yawest       R*4      I      Estimated Spacecraft Yaw in degrees
c  navctl       struct   I      controls including specification of
c                               cone angle, crossing biases, matrix
c                               for transform of sensor frame to 
c                               spacecraft frame
c  nad_bod      R*4      O      size 3 Unit Nadir Vector in Spacecraft Frame
c  nadbodfl     I*4      O      flag for the the nadir vector, 1 - bad
c  ssdec        R*4      O      Sine of the solar declination
c
c  By: W. Robinson, GSC, 12 Apr 93 from version by J Whiting
c
c  Notes: This was taken from the original earth without much mods.
c         due to the flags for earth sensors and deterministic nadir
c         solution added, it may be good to re-write this from scratch
c               WR  
c
c  Modification History:
c       W. Robinson, GSC, 16 Apr 93  This is a change in body of code
c               from original: flip indicies on array hv so that vector
c               routines may be called with this array.
c       F. Patt, GSC, 23 Mar 95:  Fixed scanner azimuth calculation to
c               allow for arbitrary scanner orientation
c       F. Patt, GSC, 23 Sep 97:  Added horizon sensor in- and out-crossing
c               angle calibrations as 4th-order polynomials
c       F. Patt, GSC, 16 Oct 97:  Corrected a bug in the oblateness modeling.
c       F. Patt, SAIC GSC, 24 Apr 98:  Added horizon crossing biases to 
c               initial calculation of rho at scanner axis azimuth
c       F. Patt, SAIC GSC, 11 May 98:  Modified calculation of horizon angle
c               by adjusting flattening factor to fit apparent horizon 
c               scanner triggering height.
c       F. Patt, SAIC GSC, 19 June 98:  Further modified calculation of the
c               horizon angle to include a correction based on the solar
c               declination angle and the subsatellite latitude.
c       F. Patt, SAIC GSC, 5 Feb 99:  Set nad_bod to zero if flags are set.
c       F. Patt, SAIC GSC, 23 Aug 1999:  Modified horizon sensor calibrations
c               to use vendor-provided calibration tables; modified horizon 
c               angle calculation seasonal correction based on analysis of 
c               sensor data with revised calibration.
c       F. Patt, SAIC GSC, 9 Dec 1999:  Modified horizon angle calculation 
c               seasonal correction based on further analysis of sensor data.
c       F. Patt, SAIC GSC, 20 Dec 1999:  Modified flattening factor used for
c               horizon angle calculation based on analysis of island targets.
c       F. Patt, SAIC GSC, 7 Jan 2000:  Incorporated new form of horizon angle 
c                calculation to allow seasonal shifting of ellipsoid in place 
c                of adjustment of subsatellite radius.
c       F. Patt, SAIC GAC, 4 Apr 2001:  Modified flattening factor to include
c                seasonal correction based on island target results.

      implicit none
c
#include "nav_cnst.fin"
#include "navctl_s.fin"
c
      type(navctl_struct) :: navctl
c
      integer*4 nadbodfl, widphfl1, widphfl2
      real*4 pos(3),vel(3),yawest,nad_bod(3),
     1     widphse1(2), widphse2(2)
c
      real*8 xnad(3)

      real*4 xang(4), v2(3), xang1(4)
      real*4 sca1, sca2, cca1, cca2, hv0(3), hv(3,4)
      real*4 sc0x(3), sc0y(3), sc0z(3), sc0zm, pxy2, posm
      real*4 n(3), nm, ndoty, ndotz, yawref
      real*4 aa1, aa2, at1, da1, da2, at2, rh1, rh2, xaz(4)
      real*4 lat, clat, ffrom1, ffrom2, clat2, re2
      real*4 slat, cxaz, cxaz2, rho(4), crho
      real*4 rho1, rho2, hv1(3), hv2(3), or1(3), or2(3), or3(3),
     1     norm1, norm2, norm3, vmag, nad1(3), ff, ssdec, ssdfac(3)
      real*4 aa, ab, ac, af, ah, ai, aj, dx, ftmp
      real*8 mv(3), state(3,3), xnadm
      real*8 a, b, c
      integer*2 i,j,k, fl1, fl2, iax
      integer*4 ierr,i3
c
      data rho/4*0.0/
      data ssdfac/0.0,11.5,-4.5/
      data i3/3/
c
c
c
c       start, set local flags from the input flags and see 
c       if this index should be done
c
      fl1 = widphfl1
      fl2 = widphfl2
c
      if( ( fl1 .eq. 1 ) .and. ( fl2 .eq. 1 ) ) then
        nadbodfl = 1
        do i=1,3
           nad_bod(i) = 0.0
        end do
      else
        nadbodfl = 0

c           Initialize:
c           Convert scanner alignment angles to transformation matrices
c
c           Convert Earth width and split-to-index to in- & 
c           out-crossing angles (xang)
c
        if( fl1 .eq. 0 ) then
           xang1(1) = widphse1(2) - widphse1(1)/2.d0 
           xang1(2) = widphse1(2) + widphse1(1)/2.d0 
        else 
           xang1(1) = -135.
           xang1(2) = -135.
        end if
c
        if( fl2 .eq. 0 ) then
           xang1(3) = widphse2(2)  - widphse2(1)/2.d0
           xang1(4) = widphse2(2)  + widphse2(1)/2.d0
        else 
           xang1(3) = 135.
           xang1(4) = 135.
        end if

c           Apply horizon angle calibrations

        call earcal(xang1, xang)

        do i=1,4
           xang(i) = xang(i)/radeg
           xang1(i) = xang1(i)/radeg
        end do 

c           Compute horizon vectors in scanner coordinates 
c           and immediately transform to spacecraft coordinates
c           hv = [cos(xang)sin(ca), sin(xang)sin(ca), cos(ca)]
c
        if( fl1 .eq. 0 ) then
          sca1 = sin( navctl%ear1sca / radeg  )
          cca1 = cos( navctl%ear1sca / radeg  )
          do i=1,2
            hv0(1) = cos(xang(i))*sca1
            hv0(2) = sin(xang(i))*sca1
            hv0(3) = cca1
            call matvec( navctl%ear_mat(1,1,1), hv0, hv(1,i) )
          enddo
        end if
c
c           note that we define both scanner systems approx paallel
c           within the limits of the orientation
c
        if( fl2 .eq. 0 ) then
          sca2 = sin( navctl%ear2sca / radeg )
          cca2 = cos( navctl%ear2sca / radeg )
          do i=3,4
            hv0(1) = cos(xang(i))*sca2
            hv0(2) = sin(xang(i))*sca2
            hv0(3) = cca2
            call matvec( navctl%ear_mat(1,1,2), hv0, hv(1,i) )
          enddo
        end if

c           Compute model nadir-to-horizon angles (rho) 
c           corresponding to vectors:
c           Compute (zero attitude) s/c axes from pos. & vel. vectors
c
        v2(1) = vel(1) - omegae*pos(2)
        v2(2) = vel(2) + omegae*pos(1)
        v2(3) = vel(3)

        call crossp(pos,v2,sc0z)
        sc0zm = sqrt(sc0z(1)*sc0z(1)+sc0z(2)*sc0z(2)+sc0z(3)*sc0z(3))
        pxy2 = pos(1)*pos(1)+pos(2)*pos(2)
        posm = sqrt(pxy2+pos(3)*pos(3))
        do i=1,3 
          sc0x(i) = -pos(i)/posm
           sc0z(i) = sc0z(i)/sc0zm
        enddo
        call crossp(sc0z,sc0x,sc0y)

c           Compute local North vector (=<pos>X<Earth Centered Z-axis>X<pos>)
c
        n(1) = -pos(1)*pos(3)
        n(2) = -pos(2)*pos(3)
        n(3) =  pxy2
        nm = sqrt(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))
        do i=1,3 
          n(i) = n(i)/nm
        enddo

c           Compute s/c (geocentric) latitude
c
        lat = asin(pos(3)/posm)

c           Compute sub-satellite Earth radius (Wertz p102, eq 4-23)
c
        clat = cos(lat)

c           Flattening factor with seasonal correction
        ff = 1.d0/1.85d2 - ssdec/6.4d2
c        ff = 1.d0/2.0d2
        ffrom1 = 1.d0 - ff
        ffrom2 = 2.d0 - ff
        ftmp = ffrom2*ff/(ffrom1*ffrom1)
        clat2 = clat*clat

c           Set-up to compute rho values (Wertz p102, eq 4-24)
c
        re2 = (re + ssdfac(3)*ssdec)**2
        slat = pos(3)/posm
        aa = 1.0
        ab = 1.0 + ftmp*clat*clat
        ac = 1.0 + ftmp*slat*slat
        af = 2.*ftmp*slat*clat
        ah = posm*af
        ai = 2.*posm*ac
        aj = posm*posm*ac - re2

c         Add seasonal correction to ellipsoid

        dx = ssdfac(1) + ssdfac(2)*ssdec
        ah = ah - 2.*dx*clat/(ffrom1*ffrom1)
        ai = ai - 2.*dx*slat/(ffrom1*ffrom1)
        aj = aj - (2.*dx*pos(3) - dx*dx)/(ffrom1*ffrom1)
c           Compute yaw reference azimuth (angle betw. north & orbit y-axis)
c
        ndoty = n(1)*sc0y(1)+n(2)*sc0y(2)+n(3)*sc0y(3)
        ndotz = n(1)*sc0z(1)+n(2)*sc0z(2)+n(3)*sc0z(3)
        yawref = -atan2(ndotz,ndoty)

c           Compute scanner axis azimuth (scanner axis align.: 
c           yaw~= +90 & +270)
c
        da1 = atan2(navctl%ear_mat(3,3,1),navctl%ear_mat(2,3,1))
        aa1 = yawref + yawest/radeg + da1
        da2 = atan2(navctl%ear_mat(3,3,2),navctl%ear_mat(2,3,2))
        aa2 = yawref + yawest/radeg + da2

c           Compute in- & out-xing azimuths (xaz) relative to North
c           xaz = aa +/- atan(tan(ca)sin(ew/2)) (approx)
c
        if( fl1 .eq. 0 ) then

c          First compute rho at scanner axis azimuth 
           cxaz = cos(aa1)
           cxaz2 = cxaz*cxaz
           a = ai*ai - 4.d0*aj*ac
           b = (4.d0*aj*af - 2.d0*ah*ai)*cxaz
           c = (ah*ah+4.d0*aj*(aa-ab))*cxaz2 - 4.d0*aj*aa
           rh1 = atan(2.d0*a/(sqrt(b*b - 4.d0*a*c) - b))
c            Need to include average crossing bias
           rh1 = rh1 + (navctl%e1biasic + navctl%e1biasoc)/(2.d0*radeg)

c          Compute azimuth offsets to horizon crossings
          at1 = acos((sca1**2*cos(xang(2)-xang(1)) 
     *          + cca1**2-cos(rh1)**2) / sin(rh1)**2)/2.
          xaz(1) = aa1 + at1
          xaz(2) = aa1 - at1
        end if
c
        if( fl2 .eq. 0 ) then

c          First compute rho at scanner axis azimuth 
           cxaz = cos(aa2)
           cxaz2 = cxaz*cxaz
           a = ai*ai - 4.d0*aj*ac
           b = (4.d0*aj*af - 2.d0*ah*ai)*cxaz
           c = (ah*ah+4.d0*aj*(aa-ab))*cxaz2 - 4.d0*aj*aa
           rh2 = atan(2.d0*a/(sqrt(b*b - 4.d0*a*c) - b))
c            Need to include average crossing bias
           rh2 = rh2 + (navctl%e2biasic + navctl%e2biasoc)/(2.d0*radeg)

c          Compute azimuth offsets to horizon crossings
          at2 = acos((sca2**2*cos(xang(4)-xang(3))
     *          + cca2**2 - cos(rh2)**2) / sin(rh2)**2)/2.
          xaz(3) = aa2 + at2
          xaz(4) = aa2 - at2
        end if

c           Compute nadir-to-horizon angles (rho):
c
        if( fl1 .eq. 0 ) then
          do i=1,2
            cxaz = cos(xaz(i))
            cxaz2 = cxaz*cxaz
            a = ai*ai - 4.d0*aj*ac
            b = (4.d0*aj*af - 2.d0*ah*ai)*cxaz
            c = (ah*ah+4.d0*aj*(aa-ab))*cxaz2 - 4.d0*aj*aa
            rho(i) = atan(2.d0*a/(sqrt(b*b - 4.d0*a*c) - b))
          enddo
        end if
c
        if( fl2 .eq. 0 ) then
          do i=3,4
            cxaz = cos(xaz(i))
            cxaz2 = cxaz*cxaz
            a = ai*ai - 4.d0*aj*ac
            b = (4.d0*aj*af - 2.d0*ah*ai)*cxaz
            c = (ah*ah+4.d0*aj*(aa-ab))*cxaz2 - 4.d0*aj*aa
            rho(i) = atan(2.d0*a/(sqrt(b*b - 4.d0*a*c) - b))
          enddo
        end if

c               Add horizon scanner triggering biases to angles

        rho(1) = rho(1) + navctl%e1biasic/radeg
        rho(2) = rho(2) + navctl%e1biasoc/radeg
        rho(3) = rho(3) + navctl%e2biasic/radeg
        rho(4) = rho(4) + navctl%e2biasoc/radeg

        if (navctl%lvdbug.gt.2) then
           do i=1,4
              write (68,*) (hv(j,i),j=1,3),rho(i), xang(i)
           end do
        end if

c           depending on whether one or 2 earth sensors are 
c           available, process differently
c
        if( fl1 .eq. 1  .or.  fl2 .eq. 1 ) then
c
c             use deterministic approach. First, move the good 
c             horizon and rho vectors to new locations
c
          if( fl1 .eq. 0 ) then
c
            do iax = 1, 3
              hv1(iax) = hv(iax,1)
              hv2(iax) = hv(iax,2)
            end do
            rho1 = rho(1)
            rho2 = rho(2)
c
          else
c
            do iax = 1, 3
              hv1(iax) = hv(iax,3)
              hv2(iax) = hv(iax,4)
            end do
            rho1 = rho(3)
            rho2 = rho(4)
c
          end if
c
c             Construct an orthonormal system out of the 
c             horizon vectors: sum, difference and cross product
c
          do iax = 1, 3
            or1(iax) = hv1(iax) + hv2(iax)
            or3(iax) = hv1(iax) - hv2(iax)
          end do
c
          call crossp( hv1, hv2, or2 )
c
          norm1 = vmag( or1 )
          norm2 = vmag( or2 )
          norm3 = vmag( or3 )
c
          do iax = 1, 3
            or1(iax) = or1(iax) / norm1
            or2(iax) = or2(iax) / norm2
            or3(iax) = or3(iax) / norm3
          end do
c
c             dot the nadir vector into the or1 and or3 axis
c             to get the components of the nadir in this system
c             ( remember, cos(rho) = nadir * horiz vec )
c
          nad1(1) = ( cos( rho1 ) + cos( rho2 ) ) / norm1
          nad1(2) = ( cos( rho1 ) - cos( rho2 ) ) / norm3
          nad1(3) = -sqrt( 1 - nad1(1) * nad1(1) - nad1(2) * nad1(2) )
c
c             the square root has 2 solutions. pick the one that has 
c             the greatest x component: thus, if or2(1) is < 0,
c             use the negative solution
c
c          if( or2(1) .lt. 0 ) nad1(3) = -nad1(3)
c
c             transform the nadir vector back to spacecraft system
c
          do iax = 1, 3
            nad_bod(iax) = nad1(1) * or1(iax) + 
     1                     nad1(2) * or3(iax) + 
     1                     nad1(3) * or2(iax)
          end do
c
c
        else
c
c             in this case, three or four valid vectors, use 
c             least-squares algorithm:
c             Initialize (upper right of) state matrix & 
c             measurement vector as zeroes
c
          do i=1,3 
            mv(i)=0.d0
            do j=i,3
              state(i,j)=0.d0
            enddo
          enddo

c             For each (valid) horizon vector:
c
          do i=1,4 

c             Update state matrix (upper right half) using 
c             horizon vector (hv(i)):
c               [SM] = [SM] + <HV>#<HV>     (# denotes the outer product)
c
            do j=1,3
              do k=j,3
                state(j,k) = state(j,k) + hv(j,i)*hv(k,i)
              enddo
            enddo

c               Update measurement vector using horizon 
c               vector (hv) & angle (rho):
c               <MV> = <MV> + cos(rho)*<HV>

            crho = cos(rho(i))
            do j=1,3
              mv(j) = mv(j) + crho*hv(j,i)
            enddo
          enddo

c             Update lower left hand of (symmetric) state matrix
c
          state(2,1) = state(1,2)
          state(3,1) = state(1,3)
          state(3,2) = state(2,3)

c             Invert state matrix and solve state equations for nadir vector:
c           <nad> = inv[SM]*<MV>
c
          call invert(state,mv,i3,i3,xnad,ierr)

c             Normalize the nadir vector to be safe.
c
          xnadm = dsqrt( xnad(1) * xnad(1) + 
     1                   xnad(2) * xnad(2) + 
     1                   xnad(3) * xnad(3))
          do i=1,3 
            nad_bod(i) = xnad(i)/xnadm
          enddo
c
        end if
      end if
c
      return
      end

      subroutine earcal(xang,yang)
      
c
c  earcal( xang, yang)
c
c  Purpose:  apply calibrations to the horizon scanner crossing angles
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c     xang(4)   R*4      I      Input horizon scanner crossing angles
c                               (order is HSA in, out, HSB in, out
c     yang(4)   R*4      O      Output (calibrated) angles
c
c  By: F. S. Patt, SAIC GSC, May 4, 1998
c
c  Notes: This subroutine uses the look-up tables provided by OSC (and,
c     presumably, the manufacturer).
c
c  Modification History:
c
c  Modified to use cubic Lagrange instead of linear interpolation.
c  F. S. Patt, SAIC GSC, August 16, 1999
c
c  Modified to read all calibration tables from a single file.
c  F. S. Patt, SAIC GSC, August 18, 1999
      
      implicit none

      real*4 xang(4), yang(4)
      real*4 hscal(541,4), hsref(541), xfac, xf(4)
      integer*4 i, j, lenstr
      character*256 hsfile
      logical first
      
      common /calcom/ hscal, hsref, first

      data first/.true./

c  If first call, open files and read tables


      if (first) then
         first = .false.
         hsfile = '$NAVCTL/'
         call filenv(hsfile, hsfile)
         hsfile = hsfile(1:lenstr(hsfile)) // 'hs_cal.dat'
         open(1,file = hsfile)
         do j=1,541
            read(1,*) hsref(j),(hscal(j,i), i=1,4)
         end do

         close(1)
      end if

c  Find location in lookup table for each angle
      do i=1,4
         j = 3
         if (xang(i).lt.0) xang(i) = xang(i) + 360.
         dowhile ((xang(i).gt.hsref(j)).and.(j.lt.540))
           j = j + 1
         end do
         
c     Compute interpolation factors
         xfac = (xang(i)-hsref(j-1))/(hsref(j)-hsref(j-1))
         xf(1) = -xfac*(xfac-1.0)*(xfac-2.0)/6.
         xf(2) = (xfac+1.0)*(xfac-1.0)*(xfac-2.0)/2.
         xf(3) = -(xfac+1.0)*xfac*(xfac-2.0)/2.
         xf(4) = (xfac+1.0)*xfac*(xfac-1.0)/6.
         
         yang(i) = xf(1)*hscal(j-2,i) + xf(2)*hscal(j-1,i)
     *        + xf(3)*hscal(j,i) + xf(4)*hscal(j+1,i)
        
      end do

      return
      end
