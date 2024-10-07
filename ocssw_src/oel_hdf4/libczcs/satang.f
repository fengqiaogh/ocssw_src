      subroutine satang(pi,rad,tilt,roll,pitch,yaw,xlon,ylat,senz,
     *sena)
c
c  Given lons and lats and nominal satellite info, computes
c  satellite (sensor) zenith and azimuth angles
c
      implicit real*8 (a-h,o-z)
      parameter(npix=1968)
      parameter(Re=6378.0,Hs=955.0)
      parameter(altcor=(Re+Hs)/Re)
      real*4 xlon(npix),ylat(npix)
      real*4 senz(npix),sena(npix)
      real*4 tilt,roll,pitch,yaw
      data xinclp /9.28/
      data difov /0.04/   !IFOV in degrees
c
c  Satellite ground point (sub-satellite point)
      xnu = 0.0 - 10.0*difov  !stay a litle off center to keep precision
      ip = npix/2 - 10
      call subsat(pi,rad,ylat(ip),xlon(ip),xinclp,altcor,tilt,xnu,
     * pitch,roll,yaw,ylats,xlons)
      rinclp = xinclp/rad
      rlats = ylats/rad
c
c  Compute satellite zenith and azimuth
      do ip = 1,npix
       xnu = -39.36 + 0.04*float(ip)
       call czsuba(rad,altcor,tilt,xnu,pitch,roll,yaw,rpsi,rths,rth)
c 
c  Calculate the true bearing from north to pixel 
       tmp = sin(rinclp)/cos(rlats)
       tmp1 = min(tmp,1.0D0)
c       rinclo = asin(sin(rinclp)/cos(rlats))
       rinclo = asin(tmp1)
       if (rinclo .gt. rpsi)then
        rpsip = 2.0*pi + (rpsi - rinclo)
       else
        rpsip = rpsi - rinclo 
       endif
       ssca = rpsip*rad
c
c  Re-orient spacecraft azimuth with respect pixel
       rlatp = ylat(ip)/rad
       if (abs(rth) .gt. 1.0E-4)then
        cosca = (sin(rlats)-cos(rth)*sin(rlatp))/(sin(rth)*cos(rlatp))
        eps = abs(cosca)
        if (eps .ge. 1.1)then
c         write(6,*)'Error in acos argument in azimuth'
c         write(6,*)'argument = ',cosca
         if (cosca .lt. 0.0)cosca=-1.0
         if (cosca .gt. 0.0)cosca=1.0
        else if (eps .gt. 1.0)then
         if (cosca .lt. 0.0)cosca=-1.0
         if (cosca .gt. 0.0)cosca=1.0
        endif
        sca = acos(cosca)*rad
        if (ssca .le. 180.0)sca = 360.0 - sca
       else
c        if (rpsi .lt. pi)sca = rpsip*rad+180.0
c        if (rpsi .ge. pi)sca = rpsip*rad-180.0
        if (rpsi .lt. pi)sca = 180.0
        if (rpsi .ge. pi)sca = 0.0
       endif
       sena(ip) = sca
c 
c  Calculate spacecraft zenith angle 
c  re and Hs are Earth radius and satellite altitude 
       r = altcor*sin(abs(rths)) 
       temp = asin(r) 
       senz(ip) = temp*rad 
      enddo
c
      return
      end
