def j2000_to_mod(iy,idy,sec):
    # Get J2000 to MOD (precession) transformation

    # Arguments:
    # Name 		Type 	I/O 	Description
    # --------------------------------------------------------
    # iy      	 I   	 I  	Year, four digits
    # idy     	 I   	 I  	Day of year
    # sec     	R*8  	 I  	Seconds of day
    #    j2mod(3,3)	 R 	 O  	J2000 to MOD matrix
    # ported from IDL (Fred Patt), Liang Hong, 9/20/2019 
    # Liang Hong, 3/24/2020, array calculation

    import numpy as np
    from hawknav.jd import jd

    t = jd(iy,np.ones(np.shape(iy)),idy) - 2451545.50 + sec/86400.0
    t = t/36525.0

    j2mod = np.zeros((np.size(iy),3,3))
    zeta0 = np.deg2rad(( 2306.2181*t + 0.302*t*t + 0.018*t*t*t ))/3600.0
    thetap = np.deg2rad(( 2004.3109*t - 0.4266*t*t - 0.04160*t*t*t ))/3600.0
    xip = np.deg2rad(( 2306.2181*t + 1.095*t*t + 0.018*t*t*t ))/3600.0

    j2mod[:,0,0] = -np.sin(zeta0)*np.sin(xip) + np.cos(zeta0)*np.cos(xip)*np.cos(thetap)
    j2mod[:,1,0] = -np.cos(zeta0)*np.sin(xip) - np.sin(zeta0)*np.cos(xip)*np.cos(thetap)
    j2mod[:,2,0] = -np.cos(xip)*np.sin(thetap)
    j2mod[:,0,1] = np.sin(zeta0)*np.cos(xip) + np.cos(zeta0)*np.sin(xip)*np.cos(thetap)
    j2mod[:,1,1] = np.cos(zeta0)*np.cos(xip) - np.sin(zeta0) * np.sin(xip) * np.cos(thetap)
    j2mod[:,2,1] = -np.sin(xip)*np.sin(thetap)
    j2mod[:,0,2] = np.cos(zeta0)*np.sin(thetap)
    j2mod[:,1,2] = -np.sin(zeta0)*np.sin(thetap)
    j2mod[:,2,2] = np.cos(thetap)

    return j2mod
