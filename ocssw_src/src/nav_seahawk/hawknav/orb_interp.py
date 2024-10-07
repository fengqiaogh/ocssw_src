def orb_interp(torb,p,v,t_out):
    # 
    # orb_interp(torb,p,v,t_out)
    #
    # Purpose: Interpolate orbit position and velocity vectors to scan line times
    # 
    #
    # Calling Arguments:
    #
    # Name         Type    I/O     Description
    # --------     ----    ---     -----------
    # torb(*)     double      I      Array of orbit vector times in seconds of day
    # p(*,3)       float     I       Array of orbit position vectors for
    #                               each time torb
    # v(*,3)       float     I    Array of orbit velocity vectors for
    #                               each time torb
    # t_out(*)    double     I    Array of time in seconds of day 
    #                 for every scan line
    # posii(*,3)    float     O    Array of interpolated positions
    # veli(*,3)    float     O    Array of interpolated velocities
    # orbfl(*)    int     O    Interpolated orbit flags (0=good)
    #
    #
    # By: Frederick S. Patt, GSC, August 10, 1993
    #
    # Notes:  Method uses cubic polynomial to match positions and velocities
    #  at input data points.
    #
    # Modification History:
    #
    # Created IDL version from FORTRAN code.  F.S. Patt, SAIC, November 29, 2006
    # python code by Liang Hong, 2018/3/28
    # updated by Liang Hong, 2020/4/6, converted loop to matrix calculation

    import numpy as np
    
    norb = len(torb)
    if np.isscalar(t_out):
        t_out = [t_out]
    nlines = len(t_out)
    posi = np.zeros((nlines,3))
    veli = np.zeros((nlines,3))
    orbfl = np.zeros(nlines)
        
    # find 'Scan line times out of available orbit data'
    ind = np.searchsorted(torb,t_out)-1
    ind_invalid = np.where((ind<0 )| (ind>norb-2))
    ind[ind_invalid] = 0
       
    # find direct match of output position & velocity from input without interpolation
    ind_direct = np.where(t_out[:, None] == torb[None, :])
    if np.size(ind_direct)>0:
        posi[ind_direct[0],:] = p[ind_direct[1],:]
        veli[ind_direct[0],:] = v[ind_direct[1],:]
        orbfl[ind_direct[0]] = 0
        
    # cubic interpolation
    dt = torb[ind+1] - torb[ind]
    a0 = p[ind,:]
    a1 = v[ind,:]*dt[:,None]
    a2 = 3.0*p[ind+1,:] - 3.0*p[ind,:] - 2.0*v[ind,:]*dt[:,None] - v[ind+1,:]*dt[:,None]
    a3 = 2.0*p[ind,:] - 2.0*p[ind+1,:] + v[ind,:]*dt[:,None] + v[ind+1,:]*dt[:,None]   
    
    x = (t_out - torb[ind])/dt
    x2 = x*x
    x3 = x2*x
    posi = a0 + a1*x[:,None] + a2*x2[:,None] + a3*x3[:,None]
    veli = (a1 + 2.0*a2*x[:,None] + 3.0*a3*x2[:,None])/dt[:,None]
        
    posi[ind_invalid,:] = 0.0
    veli[ind_invalid,:] = 0.0
    orbfl[ind_invalid] = 1
    
    return (posi,veli,orbfl)
           
