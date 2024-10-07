        subroutine gha2000(iyr,day,gha)

c  This subroutine computes the Greenwich hour angle in degrees for the
c  input time.  It uses the model referenced in The Astronomical Almanac
c  for 1984, Section S (Supplement) and documented in "Exact 
c  closed-form geolocation algorithm for Earth survey sensors", by 
c  F.S. Patt and W.W. Gregg, Int. Journal of Remote Sensing, 1993.
c  It includes the correction to mean sideral time for nutation
c  as well as precession.

c  Calling Arguments

c  Name         Type    I/O     Description
c
c  iyr          I*4      I      Year (four digits)
c  day          R*8      I      Day (time of day as fraction)
c  gha          R*8      O      Greenwich hour angle (degrees)


c       Subprograms referenced:
c
c       JD              Computes Julian day from calendar date
c       EPHPARMS        Computes mean solar longitude and anomaly and
c                        mean lunar lontitude and ascending node
c       NUTATE          Compute nutation corrections to lontitude and 
c                        obliquity
c       
c
c       Program written by:     Frederick S. Patt
c                               General Sciences Corporation
c                               November 2, 1992
c
c       Modification History:
c
        implicit real*8 (a-h,o-z)
        common /nutcm/dpsi,eps,nutime
        common /gconst/pi,radeg,re,rem,f,omf2,omegae
        data imon/1/,nutime/-9999/

c  Compute days since J2000
        iday = day
        fday = day - iday
        jday = jd(iyr,imon,iday)
        t = jday - 2451545.5d0 + fday
        
c  Compute Greenwich Mean Sidereal Time (degrees)
        gmst = 100.4606184d0 + 0.9856473663d0*t + 2.908d-13*t*t

c  Check if need to compute nutation correction for this day
        nt = t
        if (nt.ne.nutime) then
          nutime = nt
          call ephparms(t,xls,gs,xlm,omega)
          call nutate(t,xls,gs,xlm,omega,dpsi,eps)
        end if

c  Include apparent time correction and time-of-day
        gha = gmst + dpsi*cos(eps/radeg) + fday*360.d0
        gha = dmod(gha,360.d0)
        if (gha.lt.0.d0) gha = gha + 360.d0

        return
        end

        subroutine ephparms(t,xls,gs,xlm,omega)

c  This subroutine computes ephemeris parameters used by other Mission
c  Operations routines:  the solar mean longitude and mean anomaly, and
c  the lunar mean longitude and mean ascending node.  It uses the model
c  referenced in The Astronomical Almanac for 1984, Section S 
c  (Supplement) and documented and documented in "Exact closed-form 
c  geolocation algorithm for Earth survey sensors", by F.S. Patt and 
c  W.W. Gregg, Int. Journal of Remote Sensing, 1993.  These parameters 
c  are used to compute the solar longitude and the nutation in 
c  longitude and obliquity.

c  Calling Arguments

c  Name         Type    I/O     Description
c
c  t            R*8      I      Time in days since January 1, 2000 at 
c                                12 hours UT
c  xls          R*8      O      Mean solar longitude (degrees)
c  gs           R*8      O      Mean solar anomaly (degrees)
c  xlm          R*8      O      Mean lunar longitude (degrees)
c  omega        R*8      O      Ascending node of mean lunar orbit 
c                                (degrees)
c
c
c       Program written by:     Frederick S. Patt
c                               General Sciences Corporation
c                               November 2, 1992
c
c       Modification History:
c
        implicit real*8 (a-h,o-z)

c  Sun Mean Longitude           
        xls = 280.46592d0 + 0.9856473516d0*t
        xls = dmod(xls,360.d0)
 
c  Sun Mean Anomaly             
        gs = 357.52772d0 + 0.9856002831d0*t 
        gs = dmod(gs,360.d0)

c  Moon Mean Longitude          
        xlm = 218.31643d0 + 13.17639648d0*t 
        xlm = dmod(xlm,360.d0)

c  Ascending Node of Moon's Mean Orbit  
        omega = 125.04452d0 - 0.0529537648d0*t 
        omega = dmod(omega,360.d0)
        
        return
        end

        subroutine nutate(t,xls,gs,xlm,omega,dpsi,eps)

c  This subroutine computes the nutation in longitude and the obliquity
c  of the ecliptic corrected for nutation.  It uses the model referenced
c  in The Astronomical Almanac for 1984, Section S (Supplement) and 
c  documented in "Exact closed-form geolocation algorithm for Earth 
c  survey sensors", by F.S. Patt and W.W. Gregg, Int. Journal of 
c  Remote Sensing, 1993.  These parameters are used to compute the 
c  apparent time correction to the Greenwich Hour Angle and for the 
c  calculation of the geocentric Sun vector.  The input ephemeris 
c  parameters are computed using subroutine ephparms.  Terms are 
c  included to 0.1 arcsecond.

c  Calling Arguments

c  Name         Type    I/O     Description
c
c  t            R*8      I      Time in days since January 1, 2000 at 
c                                12 hours UT
c  xls          R*8      I      Mean solar longitude (degrees)
c  gs           R*8      I      Mean solar anomaly   (degrees)
c  xlm          R*8      I      Mean lunar longitude (degrees)
c  Omega        R*8      I      Ascending node of mean lunar orbit 
c                                (degrees)
c  dPsi         R*8      O      Nutation in longitude (degrees)
c  Eps          R*8      O      Obliquity of the Ecliptic (degrees)
c                                (includes nutation in obliquity)
c
c
c       Program written by:     Frederick S. Patt
c                               General Sciences Corporation
c                               October 21, 1992
c
c       Modification History:
c
        implicit real*8 (a-h,o-z)
        common /gconst/pi,radeg,re,rem,f,omf2,omegae
        
c  Nutation in Longitude
        dpsi = - 17.1996D0*sin(omega/radeg) 
     *          + 0.2062D0*sin(2.0D0*omega/radeg)
     *          - 1.3187D0*sin(2.0D0*xls/radeg) 
     *          + 0.1426D0*sin(gs/radeg) 
     *          - 0.2274D0*sin(2.0D0*xlm/radeg) 
        
c  Mean Obliquity of the Ecliptic       
        epsm = 23.439291d0 - 3.560d-7*t 

c  Nutation in Obliquity 
        deps = 9.2025D0*cos(omega/radeg) + 0.5736D0*cos(2.0D0*xls/radeg)

c  True Obliquity of the Ecliptic 
        eps = epsm + deps/3600.d0

        dpsi = dpsi/3600.d0

        return
        end
