def propagate(rates):
    # Procedure to propagate angular rates to generate quaternions
    # input: rates
    # output: quats
    # reference: propagate.pro by Fred Patt
    # ported to Python by Liang Hong, 2/19/2020
    # updated by Liang Hong, 12/3/2020
    
    import numpy as np
    from hawknav.qprod import qprod
    
    nr = np.shape(rates)[0]
    quats = np.zeros((nr+1,4))
    
    quats[0,:] = [0.0,0.0,0.0,1.0]
    qr = np.zeros((nr,4))
    
    qr[:,0:3] = np.sin(rates/2.0)
    qr[:,3] = np.sqrt(1.0 - np.sum(qr[:,0:3]**2,axis=1))
    
    for i in range(0,nr):
        qt = qprod(quats[i,:],qr[i,:])
        qt = qt/np.sqrt(np.sum(qt*qt))   # LH , 12/3/2020
        quats[i+1,:] = qt
    
    return quats
