def tle2orb(TLEline1,TLEline2,t_start,t_end,t_interval):    
    # derived from sgp4_rot which rotateS the SGP4 orbit vectors to the ECR frame
    # converted from: sgp4_rot.f
    # inherited from TLE2ECR.py
    # Liang Hong, 2/26/2020
    # input: TLEline1,TLEline2, two lines from the TLE file
    # e.g. 1 43820U 18099BQ  20019.13664524  .00000344  00000-0  36192-4 0  9993
    #      2 43820  97.7183  93.1246 0011684 330.0533  30.0019 14.95507336 60873
    # input: t_start, in UTC
    # input: t_end, in UTC
    # input: t_interval, in seconds
    # output: orb, ECR data [doy.partialDay pos[0-2], vel[0-2]]
    # Liang Hong, 4/8/2020
    
    import sgp4
    import time
    import numpy as np
    import datetime
    from sgp4.api import Satrec
    from sgp4.api import jday    
    from hawknav.jd import jd
    
    OMEGAE = 7.29211585494e-5
    
    satellite = Satrec.twoline2rv(TLEline1,TLEline2)
    time_tuple = time.gmtime(t_start)
    jd, fr = jday(time_tuple.tm_year, time_tuple.tm_mon, time_tuple.tm_mday, time_tuple.tm_hour, time_tuple.tm_min, time_tuple.tm_sec)
    ts = np.arange(t_start,t_end,t_interval)    # time stamps in UTC seconds
    dt_array = (ts-t_start)/86400.0
    fr_array = fr+dt_array
    jd_array = jd*np.ones(np.size(fr_array))
    jd_array = jd_array + np.modf(fr_array)[1]
    fr_array = np.modf(fr_array)[0]    
    e, r, v = satellite.sgp4_array(jd_array, fr_array)
    nrec = np.size(jd_array)
    
    # Compute days since J2000
    t = jd_array - 2451545 + fr_array
    
    # Compute Greenwich Mean Sidereal Time    (degrees)
    gmst = 100.4606184 + 0.9856473663*t + 2.908e-13*t*t
    
    # Add time-of-day to get GHA
    gha = gmst + fr_array*360
    gha = np.mod(gha,360)
    
    # Rotate position and velocity vectors
    cogha = np.cos(np.deg2rad(gha))
    sigha = np.sin(np.deg2rad(gha))
    posr = np.zeros((nrec,3))
    posr[:,0] = r[:,0]*cogha + r[:,1]*sigha
    posr[:,1] = r[:,1]*cogha - r[:,0]*sigha
    posr[:,2] = r[:,2]
    
    # Include Earth rotation rate in velocity
    velr = np.zeros((nrec,3))
    velr[:,0] = v[:,0]*cogha + v[:,1]*sigha + OMEGAE*posr[:,1]
    velr[:,1] = v[:,1]*cogha - v[:,0]*sigha - OMEGAE*posr[:,0]
    velr[:,2] = v[:,2]
    
    yd = [datetime.datetime.utcfromtimestamp(t).timetuple().tm_yday for t in ts]
    tt = yd + fr_array
    
    orb = np.concatenate((tt[:,None],posr,velr),axis=1)
    return orb