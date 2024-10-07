def orb2lla(orb):
    # Compute geodetic longitude, latitude and altitude 
    # from Cartesian orbit vector in ECR coordinates
    # Uses geodetic approximation from Patt and Gregg, IJRS, 1993.
    # input: orb: position vector[0,1,2]
    # reference: orb2lla.pro by Fred Patt
    # ported to Python by Liang Hong, 2/13/2020
    
    import numpy as np
    
    RE = 6378.137           # Earth equatorial radius (km)
    REM = 6371.0            # Earth mean radius (km)
    F = 1.0/298.257       # Earth flattening factor
    OMF2 = (1.0-F)**2
        
    # Compute longitude
    lon = np.rad2deg(np.arctan2(orb[:,1],orb[:,0]))
    
    # Compute geodetic latitude 
    rad = np.sqrt(np.square(orb[:,0])+np.square(orb[:,1])+np.square(orb[:,2]))
    omf2p = (OMF2*REM + rad - REM)/rad
    pxy = np.square(orb[:,0])+np.square(orb[:,1])
    temp = np.sqrt(np.square(orb[:,2]) + omf2p*omf2p*pxy)
    lat = np.rad2deg(np.arcsin(np.array(orb[:,2])/temp))

    # Compute altitude
    clatg = np.cos(np.arctan(OMF2*np.tan(np.deg2rad(lat))))
    rl = RE*(1.0-F)/np.sqrt(1.0-(2.0-F)*F*np.square(clatg))
    alt = rad - rl
    
    return lon,lat,alt


