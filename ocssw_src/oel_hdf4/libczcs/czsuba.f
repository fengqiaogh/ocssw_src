      subroutine czsuba(rad,altcor,tilt,xnu,pitch,roll,yaw, 
     * rpsi,rths,rth) 
c 
c  Internal navigation routine to compute trackline effective 
c  bearing to pixel in radians (rpsi), effective scan angle 
c  modified by tilt (rths), and earth surface distance in radians
c  from ground point to pixel (rth).
c
      implicit real*8 (a-h,o-z)
      real*4 tilt,pitch,roll,yaw
c 
      pid2 = atan(1.0)*2.0 
      pi = acos(-1.0) 
      pi3d2 = pid2*3.0 
      rnu = xnu/rad 
      t = tilt*0.5 
      rt = t/rad 
c  If tilt = 0, it's easy 
      if (tilt .eq. 0.0)then 
       if (xnu .ge. 0.0)then 
        rpsi = pid2 
       else 
        rpsi = pi3d2 
       endif 
       rths = abs(rnu)
       go to 2 
      endif 
c  If not, ... 
c   Calculate theta-n and phi-n from Scripps Geolocation Algorithm 
c   Report, Wilson, et al 
      root = 1.0/sqrt(2.0) 
      cosc = root*(sin(rt)*cos(rnu)+cos(rt)) 
      rthr = 2.0*acos(cosc) 
      rphim = atan(sin(rnu)/(sin(rt)-cos(rt)*cos(rnu)))+pi 
c   Calculate theta and phi 
      rthp = acos(-sin(rthr)*cos(rphim)) 
      rtan = sin(rthr)/cos(rthr)
      rphip = atan(rtan*sin(rphim))
      if (tilt .lt. 0.0)then
       rphip = pi+rphip 
      endif
c   Roll, pitch, yaw corrections 
      if (roll .ne. 0.0 .or. pitch .ne. 0.0 .or. yaw .ne. 0.0)then 
       u = sin(rthp)*cos(rphip) 
       v = sin(rthp)*sin(rphip) 
       w = cos(rthp) 
       y = yaw/rad 
       r = roll/rad 
       p = pitch/rad
       xp = (1.0+y*r*p)*u + (-y+r*p)*v + p*w
       yp = y*u + v - r*w 
       zp = (-p+y*r)*u + (-p*y+r)*v + w 
       rths = acos(zp) 
       rpsi = atan2(yp,xp) 
      else 
       rpsi = rphip 
       rths = rthp 
      endif 
c 
c  Earth distance 
2     sina = altcor*sin(rths)
      a = asin(sina)
      rth = a - rths 
c 
      return 
      end 
