def ned2ecr(lon,lat):
    # Program to compute the transformation from North, East, Down to ECR
    #  at a given geodetic latitude and longitude
    # input: lon, lat
    # output: n2e
    # Ported from ned2ecr.pro by Fred Patt
    # Liang Hong, 2/19/2020
    
    import numpy as np
    
    sla = np.sin(np.deg2rad(lat))
    cla = np.cos(np.deg2rad(lat))
    slo = np.sin(np.deg2rad(lon))
    clo = np.cos(np.deg2rad(lon))
    n = [-sla*clo,-sla*slo,cla]
    e = [-slo,clo,0.0]
    d = [-cla*clo,-cla*slo,-sla]
    
    n2e = [n, e, d]
    
    return n2e
 
