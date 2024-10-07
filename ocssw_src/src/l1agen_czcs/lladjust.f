      subroutine lladjust(tilt,roll,pitch,yaw,senz,xlon,ylat)
 
c  Adjusts navigation for CZCS pixel longitudes and latitudes
c   This version applies an along-track (orbit) and yaw adjustment
c   Adjustments are assumes to be small angles
 
      real senz(1968),xlon(1968),ylat(1968)
      real vll(3),vecn(3),veco(3),dvo(3),dvy(3)
      real*8 th,tmp1,tmp2,sini
      data xincl/99.28/
      data orbadj/0.046/,yawadj/0.18/
 
      real*8 pi,radeg,re,rem,f,omf2,omegae
      common /gconst/pi,radeg,re,rem,f,omf2,omegae
 
 
c  First find nadir pixel
      szm = 90.
      i = 900
      dowhile (senz(i).lt.szm)
         szm = senz(i)
         i = i + 1
      end do
      i = i - 1
 
c  Compute orbit normal unit vector
      
      sini = sin(xincl/radeg)
      cosi = cos(xincl/radeg)
 
c    Orbit angle from ascending node
c      th = radeg*asin(sin(ylat(i)/radeg)/sini)
      tmp1 = ylat(i)/radeg
      tmp2 = dsin(tmp1)/sini
      if (tmp2.gt.1.d0) then
         th = 90.d0
      else if (tmp2.lt.(-1.d0)) then
         th = -90.d0
      else
         th = radeg*dasin(tmp2)
      end if
  
c    Longitude of ascending node
      xlnn = xlon(i) - radeg*atan(tan(th/radeg)*cosi)
 
c    Orbit normal vector
      veco(1) = sini*sin(xlnn/radeg)
      veco(2) = -sini*cos(xlnn/radeg)
      veco(3) = cosi
 
c  Compute unit vector along spacecraft nadir
c    Determine along-track angle from zenith and tilt
      if ((tilt+pitch).gt.0) then
         tal = senz(i) - tilt - pitch
      else
         tal = -senz(i) - tilt - pitch
      end if
      th = th - tal
 
c    Compute nadir vector
      vecn(1) = cos(th/radeg)*cos(xlnn/radeg) 
     *     -  sin(th/radeg)*cosi*sin(xlnn/radeg) 
      vecn(2) = sin(th/radeg)*cosi*cos(xlnn/radeg) 
     *     + cos(th/radeg)*sin(xlnn/radeg)
      vecn(3) = sin(th/radeg)*sini
 
c    Scale both vectors by adjustment angle
      do i=1,3
         vecn(i) = vecn(i)*yawadj/radeg
         veco(i) = veco(i)*orbadj/radeg
      end do
      
c  Adjust pixel lons/lats
 
c    For every pixel
      do i=1,1968
 
c      Convert to a vector
         vll(1) = cos(xlon(i)/radeg)*cos(ylat(i)/radeg)
         vll(2) = sin(xlon(i)/radeg)*cos(ylat(i)/radeg)
 
        vll(3) = sin(ylat(i)/radeg)
 
c      Compute adjustment vectors
         call crossp(vecn,vll,dvy)
         call crossp(veco,vll,dvo)
 
c      Adjust vector
         do j=1,3
            vll(j) = vll(j) + dvy(j) + dvo(j)
         end do
 
c      Compute adjusted lon and lat
         xlon(i) = radeg*atan2(vll(2),vll(1))
         ylat(i) = radeg*asin(vll(3))
 
      end do
 
      return
      end

