def l_sun(iyr,iday,sec):
    # Computes unit Sun vector in geocentric rotating coodinates, using 
    # subprograms to compute inertial Sun vector and Greenwich hour angle
    # converted from l_sum.pro
    # Liang Hong, 2017/6/26
    
    import math
    import numpy as np    
    from hawknav.gha2000 import gha2000
    
    # Get unit Sun vector in geocentric inertial coordinates
    [su,rs] = sun2000(iyr,iday,sec)

    # Get Greenwich mean sideral angle
    day = iday + sec/86400.0 
    gha = gha2000(iyr,day)
    ghar = np.deg2rad(gha)
        
    # Transform Sun vector into geocentric rotating frame
    n = np.size(day)    
    sunr = np.zeros((3,n))
    sunr[0][:] = su[0][:]*np.cos(ghar) + su[1][:]*np.sin(ghar)
    sunr[1][:] = su[1][:]*np.cos(ghar) - su[0][:]*np.sin(ghar)
    sunr[2][:] = su[2][:]
    
    return [sunr,rs]
    
def sun2000(iyr,iday,sec):
    # This subroutine computes the Sun vector in geocentric inertial 
    # (equatorial) coodinates.  It uses the model referenced in The 
    # Astronomical Almanac for 1984, Section S (Supplement) and documented
    # for the SeaWiFS Project in "Constants and Parameters for SeaWiFS
    # Mission Operations", in TBD.  The accuracy of the Sun vector is
    # approximately 0.1 arcminute.
    
    import math
    import numpy as np    
    from hawknav.jd import jd
    from hawknav.gha2000 import gha2000,ephparms,nutate
    
    xk = 0.0056932        #Constant of aberration 
    imon = np.ones(np.shape(iyr))
    
    # common nutcm,dpsi,eps,nutime
        
    #  Compute floating point days since Jan 1.5, 2000 
    #   Note that the Julian day starts at noon on the specified date
    jday = jd(iyr,imon,iday)
    t = jday - 2451545.0 + (sec-43200.0)/86400.0
    if np.isscalar(t):
        n = 1
    else:
        n = len(t)
    sun = np.zeros((3,n))
    
    # Compute solar ephemeris parameters
    [xls,gs,xlm,omega] = ephparms(t)
    
    # Check if need to compute nutation corrections for this day
    [dpsi,eps,epsm] = nutate(t,xls,gs,xlm,omega)
    
    # Compute planet mean anomalies
    #  Venus Mean Anomaly     
    g2 = 50.40828 + 1.60213022*t     
    g2 = g2 % 360
    
    #  Mars Mean Anomaly         
    g4 = 19.38816 + 0.52402078*t     
    g4 = g4 % 360
    
    # Jupiter Mean Anomaly 
    g5 = 20.35116 + 0.08309121*t     
    g5 = g5 % 360
    
    # Compute solar distance (AU)
    rs = 1.00014-0.01671*np.cos(np.deg2rad(gs))-0.00014*np.cos(2.0*np.deg2rad(gs))
    
    # Compute Geometric Solar Longitude 
    dls =     (6893.0 - 4.6543463e-4*t)*np.sin(np.deg2rad(gs)) \
             + 72.0*np.sin(2.0*np.deg2rad(gs))         \
             - 7.0*np.cos(np.deg2rad(gs - g5))         \
             + 6.0*np.sin(np.deg2rad(xlm - xls))         \
             + 5.0*np.sin(np.deg2rad(4.0*gs - 8.0*g4 + 3.0*g5)) \
             - 5.0*np.cos(np.deg2rad(2.0*gs - 2.0*g2))     \
             - 4.0*np.sin(np.deg2rad(gs - g2))         \
             + 4.0*np.cos(np.deg2rad(4.0*gs - 8.0*g4 + 3.0*g5)) \
             + 3.0*np.sin(np.deg2rad(2.0*gs - 2.0*g2))     \
             - 3.0*np.sin(np.deg2rad(g5))             \
             - 3.0*np.sin(np.deg2rad(2.0*gs - 2.0*g5))

    xlsg = xls + dls/3600.0
    
    # Compute Apparent Solar Longitude; includes corrections for nutation 
    #  in longitude and velocity aberration
    xlsa = xlsg + dpsi - xk/rs
    
    #  Compute unit Sun vector 
    sun[0][:] = np.cos(np.deg2rad(xlsa))
    sun[1][:] = np.sin(np.deg2rad(xlsa))*np.cos(np.deg2rad(eps))
    sun[2][:] = np.sin(np.deg2rad(xlsa))*np.sin(np.deg2rad(eps))

    return [sun,rs]
