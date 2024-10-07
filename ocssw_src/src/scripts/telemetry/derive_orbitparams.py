import numpy as np

def derive_orbitparams(posr, velr):

    # adapted from Fred Patt's IDL routine get_ephem_from_hkt.pro
    omegae = 7.29211585494e-5
    for i in np.arange(len(posr)):
        velr[i][0] += posr[i][1]*omegae
        velr[i][1] -= posr[i][0]*omegae

    # adapted from Fred Patt's IDL routine orb2lla.pro
    re = 6378.137      # Earth equatorial radius (km)
    rem = 6371.0       # Earth mean radius (km)
    f = 1./298.257     # Earth flattening factor
    omf2 = (1.0-f) * (1.0-f)

    xyz = np.array([p/1000.0 for p in posr], dtype=np.float64)  # convert meters to km
    x, y, z = [xyz[:,i] for i in np.arange(3)]  # separate coords

    # Compute longitude
    lon = np.arctan2(y, x)

    # Compute geodetic latitude
    rad = np.linalg.norm(xyz,axis=1)  # Euclidean distance
    omf2p = (omf2*rem + rad - rem)/rad
    pxy = x*x + y*y
    temp = np.sqrt(z*z + omf2p*omf2p*pxy)
    lat = np.arcsin(z/temp)

    # Compute altitude
    clatg = np.cos(np.arctan(omf2*np.tan(lat)))
    rl = re*(1.0-f)/np.sqrt(1.0-(2.0-f)*f*clatg*clatg)
    alt = rad - rl

    # return as dictionary
    orbitparams = {}
    orbitparams['posr'] = posr
    orbitparams['velr'] = velr
    orbitparams['lon']  = np.rad2deg(lon)
    orbitparams['lat']  = np.rad2deg(lat)
    orbitparams['alt']  = alt * 1000.0  # convert from km to meters

    return orbitparams
