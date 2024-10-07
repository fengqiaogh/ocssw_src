        subroutine geonav_lonlat(pos,rm,coef,sun,nsta,ninc,npix,
     *          xlat,xlon)

c  This subroutine performs navigation of a scanning sensor on the 
c  surface of an ellipsoid based on an input orbit position vector and 
c  sensor orientation matrix.  It uses a closed-form algorithm for 
c  determining points on the ellipsoidal surface which involves 
c  determining the intersection of the scan plan with the ellipsoid.  
c  The scan angle values are based on the SeaWiFS LAC scan parameters,  
c  i.e., 1285 pixels centered at nadir and 1.5835 milliradians per 
c  pixel.  The scan angles navigated are determined by the scan 
c  parameters passed to the routine, i.e., the start pixel, the number  
c  of pixels and the increment.  (For GAC data these values are 147, 
c  248 and 4, respectively.)  The reference ellipsoid is set according 
c  to the scan intersection coefficients in the calling sequence; a 
c  useful model is an equatorial radius of 6378.137 km. and a 
c  flattening factor of 1/298.257, used by both the Geodetic Reference 
c  System (GRS) 1980 and the World Geodetic System (WGS) 1984.
c
c  It then computes geometric parameters using the pixel locations on
c  the Earth, the spaecraft position vector and the unit Sun vector in
c  the geocentric rotating reference frame.  The outputs are arrays of
c  geodetic latitude and longitude, solar zenith and azimuth and sensor
c  zenith and azimuth.  The azimuth angles are measured from local
c  North toward East.  Flag values of 999. are returned for any pixels
c  whose scan angle is past the Earth's horizon.

c  Calling Arguments

c  Name         Type    I/O     Description
c
c  pos(3)       R*4      I      Orbit Position Vector (km)
c  rm(3,3)      R*4      I      Sensor Orientation Matrix
c  coef(6)      R*4      I      Scan path coefficients
c  sun(3)       R*4      I      Sun unit vector in geocentric rotating
c                                reference frame
c  nsta         I*4      I      Scan start pixel number
c                                SeaWifs LAC:  1
c                                SeaWifs GAC:  147
c  ninc         R*4      I      Pixel number increment
c                                SeaWifs LAC:  1
c                                SeaWifs GAC:  4
c  npix         I*4      I      Pixels per scan (maximum value of 1285)
c                                SeaWiFS LAC:  1285
c                                SeaWiFS GAC:  248
c  xlat(*)      R*4      O      Pixel geodetic latitudes
c  xlon(*)      R*4      O      Pixel geodetic longitudes

c       Subprograms Called:

c       CROSSP          Compute cross product of two vectors


c       Program written by:     Frederick S. Patt
c                               General Sciences Corporation
c                               October 20, 1992
c
c       Modification History:   
c  Modified to correspond to paper: "Exact closed-form geolocation
c  algorithm for Earth survey sensors", International Journal of
c  Remote Sensing, Patt and Gregg 1993, by W. Gregg, 4/5/93.
c
c  Corrected sinc (radians per pixel) to 0.0015897 correspond to measured
c  value from SBRS data book.  F. S. Patt, 9/22/97
c  
c  Added out-of-plane correction to look vector (first-order approximation) 
c  and updated radians per pixel value based on navigation assessment 
c  results.  F. S. Patt, GSC, 12/8/97

        real pos(3),coef(6),rm(3,3),sun(3)
        real xlat(1),xlon(1)
        real geovec(3),no(3),up(3),ea(3),rmtq(3)
        real*8 pi,radeg,re,rem,f,omf2,omegae
        real*8 sinc,elev,sinl,cosl,h
        logical first/.true./
        common /navicom/sina(1285),cosa(1285),sinl,cosl
        common /gconst/pi,radeg,re,rem,f,omf2,omegae
        data sinc/0.0015911d0/
        data ea/0.0,0.0,0.0/

c  If first call, store array of sines and cosines of scan angles for 
c  future use
        if (first) then
          first = .false.
          call cdata

c  Compute elevation (out-of-plane) angle
          elev = sinc*1.2
          sinl = sin(elev)
          cosl = cos(elev)
          do i=1,1285
            sina(i) = sin((i-643)*sinc)*cosl
            cosa(i) = cos((i-643)*sinc)*cosl
          end do
        end if

c  Compute correction factor for out-of-plane angle
        h = (rm(2,1)*pos(1)+rm(2,2)*pos(2)+rm(2,3)*pos(3)/omf2)*2.d0

c  Compute sensor-to-surface vectors for all scan angles
        do i=1,npix
          in = ninc*(i-1) + nsta
          a = coef(1)*cosa(in)*cosa(in)+coef(2)*cosa(in)*sina(in)
     *          +coef(3)*sina(in)*sina(in)
          b = coef(4)*cosa(in)+coef(5)*sina(in)
          c = coef(6)
          r = b*b-4.d0*c*a  !begin solve quadratic equation

c  Check for scan past edge of Earth
          if (r.lt.0.) then
            xlat(i) = 999.
            xlon(i) = 999.
          else
c  Solve for magnitude of sensor-to-pixel vector and compute components
            q = (-b-sqrt(r))/(2.d0*a)

c  Add out-of-plane correction 
            q = q*(1.d0 + sinl*h/sqrt(r))
            
            Qx = q*cosa(in)
            Qy = q*sinl
            Qz = q*sina(in)

c  Transform vector from sensor to geocentric frame
            do j=1,3
              rmtq(j) = Qx*rm(1,j) + Qy*rm(2,j) + Qz*rm(3,j)
              geovec(j) = rmtq(j) + pos(j)
            end do

c    Compute geodetic latitude and longitude
            tmp = sqrt(geovec(1)*geovec(1)+geovec(2)*geovec(2))*omf2
            xlat(i) = radeg*atan2(geovec(3),tmp)
            xlon(i) = radeg*atan2(geovec(2),geovec(1))

          endif


        end do

        return
        end
