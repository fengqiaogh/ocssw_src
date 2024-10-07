      subroutine subsat(pi,rad,ylatp,xlonp,xinclp,altcor, 
     * tilt,xnu,pitch,roll,yaw,ylats,xlons) 
c 
c  Calculates the satellite ground point latitude and longitude 
c  Input 
c   rad    = radians/degrees conversion 
c   ylatp  = latitude of pixel 
c   xlonp  = longitude of pixel 
c   xinclp = equatorial trackline inclination = 9.28 deg. 
c   altcor = altitude correction factor = (R+H)/R = 1.149 
c   xnu    = scan angle to pixel for zero tilt 
c   tilt   = scan mirror tilt angle 
c   pitch  = pitch angle of satellite 
c   roll   = roll angle of satellite 
c   yaw    = yaw of spacecraft 
c  Output 
c   ylats  = latitude of spacecraft 
c   xlons  = longitude of spacecraft 
c
      implicit real*8 (a-h,o-z)
      real*4 ylatp,xlonp,tilt,pitch,roll,yaw
c
      call czsuba(rad,altcor,tilt,xnu,pitch,roll,yaw,rpsi,rths,rth)
c
      rlatp = ylatp/rad
      rlonp = xlonp/rad 
      rinclp = xinclp/rad 
c 
c  Set up quadratic variables from Scripps Algorithm Report; CZCS
c  Geolocation Algorithms
      a = sin(rlatp)-sin(rth)*sin(rpsi)*sin(rinclp) 
      b = cos(rth) 
      c = sin(rth)*cos(rpsi) 
      d = cos(rinclp)*cos(rinclp) 
      ab = a*b 
      xd = b*b + c*c 
      x = ab*ab - xd*(a*a - d*c*c) 
      xp = 0.0 
      if (x .gt. 1.0E-10)xp = sqrt(x) 
c 
c  Get the +- values from the quadratic 
      x1 = (ab+xp)/xd 
      x2 = (ab-xp)/xd
      rlats1 = asin(x1)
      rlats2 = asin(x2)
c 
c  If tilt is positive use the negative quadratic and vice versa
      if (tilt .lt. 0.0)then
       rlats = rlats1 
      else 
       rlats = rlats2 
      endif 
      ylats1 = rlats1*rad 
      ylats2 = rlats2*rad 
      ylats = rlats*rad 
c 
c  Comment out the following because not used
cc  Calculate the true bearing from north to pixel
c      rinclo = asin(sin(rinclp)/cos(rlats))
c      if (rinclo .gt. rpsi)then
c       rpsip = 2.0*pi + (rpsi - rinclo)
c      else
c       rpsip = rpsi - rinclo
c      endif
cc
cc  Calculate the longitude difference of G.P. to pixel 
c      x = sin(rpsip)*sin(rth)/cos(rlatp) 
c      dlam = asin(x) 
c      xlam = dlam*rad 
c      xlons = xlonp + xlam 
c 
      return 
      end 
