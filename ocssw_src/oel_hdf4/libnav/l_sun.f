        subroutine l_sun(iyr,iday,sec,sunr,rs)
c
c  Computes unit Sun vector in geocentric rotating coodinates, using 
c  subprograms to compute inertial Sun vector and Greenwich hour angle
c
c
c       Arguments:
c
c       Name    Type    I/O     Description
c       --------------------------------------------------------
c       IYR     I*4      I      Year, four digits (i.e, 1993)
c       IDAY    I*4      I      Day of year (1-366)
c       SEC     R*8      I      Seconds of day 
c       SUNR(3) R*4      O      Unit Sun vector in geocentric rotating 
c                                coordinates
c       RS      R*4      O      Earth-to-Sun distance (AU)
c
c       Subprograms referenced:
c
c       SUN2000         Computes inertial Sun vector
c       GHA2000         Computes Greenwich sidereal angle
c
c       Coded by:  Frederick S. Patt, GSC, September 29, 1992
c
c       Modification History:
c
c       Modifified to use new Sun and hour angle routines
c       Frederick S. Patt, November 3, 1992
c
c       Removed internal jd() function, since it is available as an
c       independent module.  B. A. Franz, GSC, November 14, 1997.

        implicit real*8 (a-h,o-z)
        real*4 sunr(3),su(3),rs
        common /gconst/pi,radeg,re,rem,f,omf2,omegae
        
c  Get unit Sun vector in geocentric inertial coordinates
        call sun2000(iyr,iday,sec,su,rs)

c  Get Greenwich mean sideral angle
        day = iday + sec/864.d2 
        call gha2000(iyr,day,gha)
        ghar = gha/radeg

c  Transform Sun vector into geocentric rotating frame
        sunr(1) = su(1)*cos(ghar) + su(2)*sin(ghar)
        sunr(2) = su(2)*cos(ghar) - su(1)*sin(ghar)
        sunr(3) = su(3)

        return
        end

        subroutine sun2000(iyr,iday,sec,sun,rs)
c
c  This subroutine computes the Sun vector in geocentric inertial 
c  (equatorial) coodinates.  It uses the model referenced in The 
c  Astronomical Almanac for 1984, Section S (Supplement) and documented
c  in "Exact closed-form geolocation algorithm for Earth survey
c  sensors", by F.S. Patt and W.W. Gregg, Int. Journal of Remote
c  Sensing, 1993.  The accuracy of the Sun vector is approximately 0.1 
c  arcminute.
c
c       Arguments:
c
c       Name    Type    I/O     Description
c       --------------------------------------------------------
c       IYR     I*4      I      Year, four digits (i.e, 1993)
c       IDAY    I*4      I      Day of year (1-366)
c       SEC     R*8      I      Seconds of day 
c       SUN(3)  R*4      O      Unit Sun vector in geocentric inertial 
c                                coordinates of date
c       RS      R*4      O      Magnitude of the Sun vector (AU)
c
c       Subprograms referenced:
c
c       JD              Computes Julian day from calendar date
c       EPHPARMS        Computes mean solar longitude and anomaly and
c                        mean lunar lontitude and ascending node
c       NUTATE          Compute nutation corrections to lontitude and 
c                        obliquity
c
c       Coded by:  Frederick S. Patt, GSC, November 2, 1992
c       Modified to include Earth constants subroutine by W. Gregg,
c               May 11, 1993.


        implicit real*8 (a-h,o-z)
        real*4 sun(3),rs
        common /nutcm/dpsi,eps,nutime
        common /gconst/pi,radeg,re,rem,f,omf2,omegae
        data xk/0.0056932/              !Constant of aberration 
        data imon/1/

c   Compute floating point days since Jan 1.5, 2000 
c    Note that the Julian day starts at noon on the specified date
        t = jd(iyr,imon,iday) - 2451545.0d0 + (sec-43200.d0)/86400.d0

c  Compute solar ephemeris parameters
        call ephparms(t,xls,gs,xlm,omega)

c  Check if need to compute nutation corrections for this day
        nt = t
        if (nt.ne.nutime) then
          nutime = nt
          call nutate(t,xls,gs,xlm,omega,dpsi,eps)
        end if

c  Compute planet mean anomalies
c   Venus Mean Anomaly  
        g2 = 50.40828D0 + 1.60213022D0*t
        g2 = dmod(g2,360.d0)

c   Mars Mean Anomaly           
        g4 = 19.38816D0 + 0.52402078D0*t
        g4 = dmod(g4,360.d0)

c  Jupiter Mean Anomaly 
        g5 = 20.35116D0 + 0.08309121D0*t
        g5 = dmod(g5,360.d0)

c  Compute solar distance (AU)
        rs = 1.00014D0 - 0.01671D0*cos(gs/radeg) 
     *       - 0.00014D0*cos(2.0D0*gs/radeg)

c  Compute Geometric Solar Longitude 
        dls =   (6893.0D0 - 4.6543463D-4*t)*sin(gs/radeg) 
     *          + 72.0D0*sin(2.0D0*gs/radeg) 
     *          - 7.0D0*cos((gs - g5)/radeg)
     *          + 6.0D0*sin((xlm - xls)/radeg) 
     *          + 5.0D0*sin((4.0D0*gs - 8.0D0*g4 + 3.0D0*g5)/radeg) 
     *          - 5.0D0*cos((2.0D0*gs - 2.0D0*g2)/radeg)
     *          - 4.0D0*sin((gs - g2)/radeg) 
     *          + 4.0D0*cos((4.0D0*gs - 8.0D0*g4 + 3.0D0*g5)/radeg) 
     *          + 3.0D0*sin((2.0D0*gs - 2.0D0*g2)/radeg)
     *          - 3.0D0*sin(g5/radeg) 
     *          - 3.0D0*sin((2.0D0*gs - 2.0D0*g5)/radeg)  !arcseconds

        xlsg = xls + dls/3600.d0

c  Compute Apparent Solar Longitude; includes corrections for nutation 
c   in longitude and velocity aberration
        xlsa = xlsg + dpsi - xk/rs

c   Compute unit Sun vector 
        sun(1) = cos(xlsa/radeg)
        sun(2) = sin(xlsa/radeg)*cos(eps/radeg)
        sun(3) = sin(xlsa/radeg)*sin(eps/radeg)
c       print *,' Sunlon = ',xlsg,xlsa,eps

        return
        end
