def triad(b1,b2,r1,r2):
    # Program to perform Triad algorithm
    # input: b1,b2,r1,r2
    # output: rbmat
    # Ported from triad.pro by Fred Patt
    # Liang Hong, 2/18/2020
    
    import numpy as np
    
    xb = np.zeros((3,3))
    xr = np.zeros((3,3))
    
    xb[0,:] = b1
    tb = np.cross(b1,b2)
    xb[2,:] = tb/np.sqrt(np.sum(tb*tb))
    xb[1,:] = np.cross(xb[2,:],xb[0,:])
    
    xr[0,:] = r1
    tr = np.cross(r1,r2)
    xr[2,:] = tr/np.sqrt(np.sum(tr*tr))
    xr[1,:] = np.cross(xr[2,:],xr[0,:])
    
    ## rbmat = xb#invert(xr)
    rbmat = np.dot(np.linalg.inv(xr),xb)
    
    return rbmat

