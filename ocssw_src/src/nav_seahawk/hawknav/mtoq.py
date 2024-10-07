def mtoq(rm):
    # Convert direction cosine matrix to equivalent quaternion
    # input: rm
    # output: q
    # Ported from mtoq.pro by Fred Patt
    # Liang Hong, 2/18/2020
    
    import numpy as np
    
    q = np.zeros(4)
    e = np.zeros(3)

    # Compute Euler angle
    cphi = (rm[0,0]+rm[1,1]+rm[2,2]-1.0)/2.0
    if (np.abs(cphi) < 0.98):
        phi = np.arccos(cphi)
    else:
        ssphi = ((rm[1,0]-rm[0,1])**2 + (rm[0,2]-rm[2,0])**2 + (rm[2,1]-rm[1,2])**2)/4.0
        phi = np.arcsin(np.sqrt(ssphi))
        if (cphi < 0):
            phi = np.pi - phi
    
    # Compute Euler axis
    e[0] = (rm[2,1]-rm[1,2])/(np.sin(phi)*2.0)
    e[1] = (rm[0,2]-rm[2,0])/(np.sin(phi)*2.0)
    e[2] = (rm[1,0]-rm[0,1])/(np.sin(phi)*2.0)
    e = e/np.sqrt(np.sum(e*e))

    # Compute quaternion
    q[0] = e[0]*np.sin(phi/2.0)
    q[1] = e[1]*np.sin(phi/2.0)
    q[2] = e[2]*np.sin(phi/2.0)
    q[3] = np.cos(phi/2.0)
    
    return q