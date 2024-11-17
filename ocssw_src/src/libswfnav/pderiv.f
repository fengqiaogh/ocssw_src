      subroutine pderiv(ngps,igyr,igday,gpsec,vecs,orbinit,
     *    iyinit,idinit,secinit,driv)
c
c  Purpose:  This routine computes partial derivatives of the orbit positions
c               with respect to the initial mean elements.  Methods and 
c               algorithms are described in TBD
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  ngps         I*4      I      Number of vectors
c  igyr         I*4      I      Time tag year
c  igday        I*4      I      Time tag day-of-year
c  gpsec(ngps)  R*8      I      Time tag seconds-of-day
c  vecs(6,ngps) R*8      I      Input orbit vectors
c  orbinit(6)   R*8      I      Initial orbital elements
c  iyinit       I*4      I      Orbital element epoch year
c  idinit       I*4      I      Orbital element epoch day-of-year
c  secinit      R*8      I      Orbital element epoch seconds-of-day
c  driv(6,3,ngps) R*8    O      Partial derivatives of orbit position vectors
c
c
c  By: Frederick S. Patt, GSC, December 21, 1993
c
c  Notes:  
c
c  Modification History:
c
      implicit none
#include "nav_cnst.fin"

      real*8 gpsec(maxlin),vecs(6,maxlin),driv(6,3,maxlin),orbinit(6)
      real*8 on(3), v(3), p(3), gm, aj2, ddif, secinit
      real*8 pm, vm, om, oxy, tmp, xmm, t, oa, am, dldi, dmdi
      real*8 sini, cosi, sinl, cosl, sino, coso, sinm, cosm
      integer*4 ngps, igyr, igday, iyinit, idinit
      integer*4 i, j, k, jd
      data gm/398600.5d0/,aj2/1.08263d-3/


c  Compute day offset of data from element epoch
      ddif = (jd(igyr,1,igday) - jd(iyinit,1,idinit))*864.d2

c  Compute constant for precession rates
      tmp = 1.5d0*aj2*sqrt(gm)*re**2/orbinit(1)**3.5

c  Compute mean motion with J2 correction
      xmm = sqrt(gm/orbinit(1)**3)*(1.d0 + 1.5*aj2*(re/orbinit(1))**2
     *   *(4.d0*cos(orbinit(3)/radeg)**2 - 1.d0))

c  Loop through orbit vectors
      do i=1,ngps

c  Compute time offset from element epoch
        t = gpsec(i) - secinit + ddif

c  Load vector in local arrays
c  Orbit velocities must be corrected for Earth rotation rate
        p(1) = vecs(1,i)
        p(2) = vecs(2,i)
        p(3) = vecs(3,i)
        v(1) = vecs(4,i) - omegae*p(2)
        v(2) = vecs(5,i) + omegae*p(1)
        v(3) = vecs(6,i)

c  Compute position and velocity magnitudes
        pm = sqrt(p(1)*p(1) + p(2)*p(2) + p(3)*p(3))
        vm = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))

c  Compute instantaneous orbit normal
        on(1) = p(2)*v(3) - p(3)*v(2)
        on(2) = p(3)*v(1) - p(1)*v(3)
        on(3) = p(1)*v(2) - p(2)*v(1)
        om = sqrt(on(1)*on(1) + on(2)*on(2) + on(3)*on(3))
        oxy = sqrt(on(1)*on(1) + on(2)*on(2))

c  Compute sines and cosines of inclination, RA of ascending node and 
c   orbit angle
        sini = oxy/om
        cosi = on(3)/om
        sinl = on(1)/oxy
        cosl = -on(2)/oxy
        coso = p(1)*cosl + p(2)*sinl
        sino = p(3)/sini
        oa = atan2(sino,coso)

c  Compute anomaly (includes correction for perigee precession)
        am = oa - orbinit(5)/radeg - tmp*t*(2.-2.5*sini*sini)
        sinm = sin(am)
        cosm = cos(am)

c  Compute partial derivatives

c   With respect to semimajor axis
        do j=1,3 
          driv(1,j,i) = (p(j) - 1.5d0*v(j)*t)/pm
        end do 

c   With respect to eccentricity
        do j=1,3
          driv(2,j,i) = 2.d0*sinm*pm*v(j)/vm - cosm*p(j)
        end do

c   With respect to RA of ascending node 
        driv(4,1,i) = -p(2)
        driv(4,2,i) = p(1)
        driv(4,3,i) = 0.d0

c  With respect to argument of perigee
        do j=1,3
          driv(5,j,i) = -orbinit(2)*(2.d0*cosm*pm*v(j)/vm + sinm*p(j))
        end do

c  With respect to mean anomaly
        do j=1,3
          driv(6,j,i) = v(j)/xmm
        end do

c  With respect to inclination
c   Includes corrections to mean motion and ascending node rates
        dmdi = -8.*sini*cosi*tmp
        dldi = sini*tmp
        do j=1,3
          driv(3,j,i) = sino*on(j)/om + driv(4,j,i)*t*dldi 
     *      + driv(6,j,i)*t*dmdi
        end do

c  Convert all angle derivatives to degrees
        do  j=1,3
          do k=3,6
            driv(k,j,i) = driv(k,j,i)/radeg
          end do
        end do

c  End of main loop
      end do

      return
      end
