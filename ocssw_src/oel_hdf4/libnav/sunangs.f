      subroutine sunangs(iyr,iday,gmt,xlon,ylat,sunz,suna)
c
c  Given year, day of year, time in hours (GMT) and latitude and 
c  longitude, returns an accurate solar zenith and azimuth angle.  
c  Based on IAU 1976 Earth ellipsoid.  Method for computing solar 
c  vector and local vertical from Patt and Gregg, 1993, Int. J. 
c  Remote Sensing.
c
c  Subroutines required: sun2000
c                        gha2000
c                        jd
c
      save  !required for IRIX compilers
      real suni(3),sung(3),up(3),no(3),ea(3)
      real*8 pi,radeg,re,rem,f,omf2,omegae
      real*8 sec,gha,ghar,day
      common /gconst/pi,radeg,re,rem,f,omf2,omegae
      data ea /0.0,0.0,0.0/
c
c  Compute sun vector
c   Compute unit sun vector in geocentric inertial coordinates
      sec = gmt*3600.0D0
      call sun2000(iyr,iday,sec,suni,rs)

c   Get Greenwich mean sidereal angle
      day = iday + sec/86400.0D0
      call gha2000(iyr,day,gha)
      ghar = gha/radeg

c   Transform Sun vector into geocentric rotating frame
      sung(1) = suni(1)*cos(ghar) + suni(2)*sin(ghar)
      sung(2) = suni(2)*cos(ghar) - suni(1)*sin(ghar)
      sung(3) = suni(3)
c
c  Convert geodetic lat/lon to Earth-centered, earth-fixed (ECEF) 
c  vector (geodetic unit vector)
      rlon = xlon/radeg
      rlat = ylat/radeg
      cosy = cos(rlat)
      siny = sin(rlat)
      cosx = cos(rlon)
      sinx = sin(rlon)
      up(1) = cosy*cosx
      up(2) = cosy*sinx
      up(3) = siny
c
c  Compute the local East and North unit vectors  
      upxy = sqrt(up(1)*up(1)+up(2)*up(2))
      ea(1) = -up(2)/upxy
      ea(2) = up(1)/upxy
      no(1) = up(2)*ea(3) - up(3)*ea(2)  !cross product
      no(2) = up(3)*ea(1) - up(1)*ea(3)
      no(3) = up(1)*ea(2) - up(2)*ea(1)

c  Compute components of spacecraft and sun vector in the
c  vertical (up), North (no), and East (ea) vectors frame
      sunv = 0.0
      sunn = 0.0
      sune = 0.0
      do j = 1,3
       sunv = sunv + sung(j)*up(j)
       sunn = sunn + sung(j)*no(j)
       sune = sune + sung(j)*ea(j)
      enddo
c
c  Compute the solar zenith and azimuth
      sunz = radeg*atan2(sqrt(sunn*sunn+sune*sune),sunv)
c  Check for zenith close to zero
      if (sunz .gt. 0.05D0)then 
       suna = radeg*atan2(sune,sunn)
      else 
       suna = 0.0D0
      endif
      if (suna .lt. 0.0D0)suna = suna + 360.0D0
c
      return
      end
