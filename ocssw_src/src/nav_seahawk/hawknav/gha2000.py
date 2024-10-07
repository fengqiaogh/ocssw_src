def gha2000(iyr,day):
    #      converted from gha2000.f
    # This subroutine computes the Greenwich hour angle in degrees for the
    # input time.  It uses the model referenced in The Astronomical Almanac
    # for 1984, Section S (Supplement) and documented in "Exact 
    # closed-form geolocation algorithm for Earth survey sensors", by 
    # F.S. Patt and W.W. Gregg, Int. Journal of Remote Sensing, 1993.
    # It includes the correction to mean sideral time for nutation
    # as well as precession.
    
    # Calling Arguments
    
    # Name         Type    I/O     Description
    # iyr          I*4      I      Year (four digits)
    # day          R*8      I      Day (time of day as fraction)
    # gha          R*8      O      Greenwich hour angle (degrees)
    
    
    #      Subprograms referenced:
    #      JD              Computes Julian day from calendar date
    #      EPHPARMS        Computes mean solar longitude and anomaly and
    #                       mean lunar lontitude and ascending node
    #      NUTATE          Compute nutation corrections to lontitude and 
    #                       obliquity
    #      
    #      Program written by:     Frederick S. Patt
    #                              General Sciences Corporation
    #                              November 2, 1992
    # Liang Hong, 3/24/2020, array calculation
    
    #common /nutcm/dpsi,eps,nutime
    #common /gconst/pi,radeg,re,rem,f,omf2,omegae
    #data imon/1/,nutime/-9999/
    imon = 1
    
    import numpy as np
    from hawknav.jd import jd
    
    # Compute days since J2000
    iday = np.fix(day)
    fday = day - iday
    jday = jd(iyr,imon,iday)
    t = jday - 2451545.5 + fday
    
    # Compute Greenwich Mean Sidereal Time (degrees)
    gmst = 100.4606184 + 0.9856473663*t + 2.908e-13*np.multiply(t,t)
    
    # Check if need to compute nutation correction for this day
    [xls,gs,xlm,omega] = ephparms(t)
    [dpsi,eps,epsm] = nutate(t,xls,gs,xlm,omega)
    
    # Include apparent time correction and time-of-day
    gha = gmst + np.multiply(dpsi, np.cos(np.deg2rad(eps))) + fday*360.0
    gha = gha + 360.0
    gha = np.mod(gha,360)
    
    return gha

def ephparms(t):

    # This subroutine computes ephemeris parameters used by other Mission
    # Operations routines:  the solar mean longitude and mean anomaly, and
    # the lunar mean longitude and mean ascending node.  It uses the model
    # referenced in The Astronomical Almanac for 1984, Section S 
    # (Supplement) and documented and documented in "Exact closed-form 
    # geolocation algorithm for Earth survey sensors", by F.S. Patt and 
    # W.W. Gregg, Int. Journal of Remote Sensing, 1993.  These parameters 
    # are used to compute the solar longitude and the nutation in 
    # longitude and obliquity.

    # Calling Arguments

    # Name         Type    I/O     Description

    # t            R*8      I      Time in days since January 1, 2000 at 
    #                               12 hours UT
    # xls          R*8      O      Mean solar longitude (degrees)
    # gs           R*8      O      Mean solar anomaly (degrees)
    # xlm          R*8      O      Mean lunar longitude (degrees)
    # omega        R*8      O      Ascending node of mean lunar orbit 
    #                               (degrees)


    #      Program written by:     Frederick S. Patt
    #                              General Sciences Corporation
    #                              November 2, 1992
    # Liang Hong, 3/24/2020, array calculation
    
    import numpy as np
    
    # Sun Mean Longitude           
    xls = 280.46592 + 0.9856473516*t
    xls = np.mod(xls,360)

    # Sun Mean Anomaly             
    gs = 357.52772 + 0.9856002831*t 
    gs = np.mod(gs,360)

    # Moon Mean Longitude          
    xlm = 218.31643 + 13.17639648*t 
    xlm = np.mod(xlm,360)

    # Ascending Node of Moon's Mean Orbit  
    omega = 125.04452 - 0.0529537648*t 
    omega = np.mod(omega,360)

    return [xls,gs,xlm,omega]
    
def nutate(t,xls,gs,xlm,omega):
    
    # This subroutine computes the nutation in longitude and the obliquity
    # of the ecliptic corrected for nutation.  It uses the model referenced
    # in The Astronomical Almanac for 1984, Section S (Supplement) and 
    # documented in "Exact closed-form geolocation algorithm for Earth 
    # survey sensors", by F.S. Patt and W.W. Gregg, Int. Journal of 
    # Remote Sensing, 1993.  These parameters are used to compute the 
    # apparent time correction to the Greenwich Hour Angle and for the 
    # calculation of the geocentric Sun vector.  The input ephemeris 
    # parameters are computed using subroutine ephparms.  Terms are 
    # included to 0.1 arcsecond.

    # Calling Arguments

    # Name         Type    I/O     Description
    
    # t            R*8      I      Time in days since January 1, 2000 at 
    #                               12 hours UT
    # xls          R*8      I      Mean solar longitude (degrees)
    # gs           R*8      I      Mean solar anomaly   (degrees)
    # xlm          R*8      I      Mean lunar longitude (degrees)
    # Omega        R*8      I      Ascending node of mean lunar orbit 
    #                               (degrees)
    # dPsi         R*8      O      Nutation in longitude (degrees)
    # Eps          R*8      O      Obliquity of the Ecliptic (degrees)
    #                               (includes nutation in obliquity)
    
    
    #      Program written by:     Frederick S. Patt
    #                              General Sciences Corporation
    #                              October 21, 1992
    
    #      Modification History:
    
    import numpy as np
    
    xls_rad = np.deg2rad(xls)
    gs_rad = np.deg2rad(gs)
    xlm_rad = np.deg2rad(xlm)
    omega_rad = np.deg2rad(omega)
    
    # Nutation in Longitude
    dpsi = - 17.1996 * np.sin(omega_rad) + 0.2062 * np.sin(2.0*omega_rad) - 1.3187 * np.sin(2.0*xls_rad) + 0.1426 * np.sin(gs_rad) - 0.2274 * np.sin(2.0*xlm_rad)
    
    # Mean Obliquity of the Eclipti#      
    epsm = 23.439291 - 3.560e-7*t 
    
    # Nutation in Obliquity 
    deps = 9.2025 * np.cos(omega_rad) + 0.5736 * np.cos(2.0*xls_rad)
    
    # True Obliquity of the Ecliptic 
    eps = epsm + deps/3600.0
    
    dpsi = dpsi/3600.0
    
    return [dpsi,eps,epsm]
