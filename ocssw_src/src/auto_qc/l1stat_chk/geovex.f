        subroutine geovex(pos,rm,coef,sun,geovec )
c*******************************************************************
c
c   geovex
c
c   purpose: modification of geovec to get vector from earth
c                center to the satellite subpoint
c
c
c   Parameters: (in calling order)
c
c  Name         Type    I/O     Description
c  ----------   ------  ---     --------------------------
c  pos(3)       R*4      I      Orbit Position Vector (km)
c  rm(3,3)      R*4      I      Sensor Orientation Matrix
c  coef(6)      R*4      I      Scan path coefficients
c  sun(3)       R*4      I      Sun unit vector in geocentric rotating
c                                reference frame
c  geovec(3)    R*4      O      vector to the point
c      
c
c   Modification history:
c      Programmer        Date            Description of change
c      ----------        ----            ---------------------
c      W. Robinson       5 Mar 96        Original development from geovec
c
c*******************************************************************
c
c       Subprograms Called:
c
c       CROSSP          Compute cross product of two vectors
c
c*******************************************************************
 
        real pos(3),coef(6),rm(3,3),sun(3)
        real xlat(1),xlon(1),solz(1),sola(1),senz(1),sena(1)
        real geovec(3),no(3),up(3),ea(3),rmtq(3)
        real*8 pi,radeg,re,rem,f,omf2,omegae
        real*8 sinc
        logical first/.true./
        common /gconst/pi,radeg,re,rem,f,omf2,omegae
        data sinc/0.0015835d0/
        data ea/0.0,0.0,0.0/

c  Compute sensor-to-surface vectors for the center angle
          a = coef(1)
          b = coef(4)
          c = coef(6)
          r = b*b-4.d0*c*a  !begin solve quadratic equation

c  Check for scan past edge of Earth
          if (r.lt.0.) then
            do j = 1, 3
              geovec(j) = 999.
            end do
          else
c  Solve for magnitude of sensor-to-pixel vector and compute components
            q = (-b-sqrt(r))/(2.d0*a)
            Qx = q

c  Transform vector from sensor to geocentric frame
            do j=1,3
              rmtq(j) = Qx*rm(1,j)
              geovec(j) = rmtq(j) + pos(j)
            end do
          end if

        return
        end
