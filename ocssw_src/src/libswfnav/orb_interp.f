        subroutine orb_interp(orbit,nlines,timref,time,posi,veli,orbfl)
c
c  prolog(orbit,timref,time,posi,veli)
c
c  Purpose: Interpolate orbit position and velocity vectors to scan line times
c  
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  orbit        struct   I      Input structure of filtered orbit data
c  nlines       I*4      I      Number of scan line times
c  timref       R*8      I      Size 3 reference time at start line 
c                                of data:  year, day, sec
c  time         R*8      I      Array of time in seconds relative to 
c                                timref for every scan line
c  posi         R*4      O      Size 3 by nlines interpolated position
c  veli         R*4      O      Size 3 by nlines interpolated velocity
c  orbfl        I*4      O      Size nlines orbit vector flags (0=good)
c
c
c  By: Frederick S. Patt, GSC, August 10, 1993
c
c  Notes:  Method uses cubic polynomial to match positions and velocities
c   at input data points.
c
c  Modification History:
c
c  Added check for array index exceeded number of scan lines in check for
c  scan lines preceding first scan line.  F.S. Patt, SAIC GSC, May 23, 2000.
c
      implicit none
#include "nav_cnst.fin"
#include "orbit_s.fin"
      type(orbit_struct) :: orbit
c
      integer*4 ind, i, j, i1, nlines, iy, id, jd
      real*8 timref(3), time(nlines)
      real*8 a0(3),a1(3),a2(3),a3(3), t, dt, x, x2, x3, ddifs
      real*4 posi(3,nlines), veli(3,nlines)
      integer*4 orbfl(nlines)
      logical first       

      ind = 1
      i1 = 1
      first = .true.
c  Make sure that orbit vector precedes first scan line
      iy = timref(1)
      id = timref(2)
      ddifs = (jd(iy,1,id) - jd(orbit%iyr,1,orbit%iday))*864.d2
      t = ddifs + timref(3)+time(i1)
      dowhile ((t.lt.orbit%torb(ind)).and.(i1.le.nlines))
        do j=1,3
          posi(j,i1) = 0.0
          veli(j,i1) = 0.0
        end do
        orbfl(i1) = 1
        print *, 'Scan line times before available orbit data'
        i1 = i1 + 1
        t = ddifs + timref(3)+time(i1)
      end do
      
      do i=i1,nlines
        t = ddifs + timref(3)+time(i)
        if ((t.gt.orbit%torb(ind+1)).or.first) then
          dowhile (t.gt.orbit%torb(ind+1))
            ind = ind + 1
            if (ind.gt.orbit%nvec) then
              i1 = 1
              go to 990
            end if
          end do
          first = .false.
c  Set up cubic interpolation
          dt = orbit%torb(ind+1) - orbit%torb(ind)
          do j=1,3
            a0(j) = orbit%pos(j,ind)
            a1(j) = orbit%vel(j,ind)*dt
            a2(j) = 3.d0*orbit%pos(j,ind+1) - 3.d0*orbit%pos(j,ind) 
     *              - 2.d0*orbit%vel(j,ind)*dt - orbit%vel(j,ind+1)*dt
            a3(j) = 2.d0*orbit%pos(j,ind) - 2.d0*orbit%pos(j,ind+1) 
     *              + orbit%vel(j,ind)*dt + orbit%vel(j,ind+1)*dt
          end do
        end if

c  Interpolate orbit position and velocity components to scan line time
        x = (t - orbit%torb(ind))/dt
        x2 = x*x
        x3 = x2*x
        do j=1,3
          posi(j,i) = a0(j) + a1(j)*x + a2(j)*x2 + a3(j)*x3
          veli(j,i) = (a1(j) + 2.*a2(j)*x + 3.*a3(j)*x2)/dt
        end do
        orbfl(i) = 0

      end do
      return

990   print *, 'Scan line times after available orbit data'
      do i=i1,nlines
        do j=1,3
          posi(j,i) = 0.0
          veli(j,i) = 0.0
        end do
      orbfl(i) = 1
      end do
      return

      end       
