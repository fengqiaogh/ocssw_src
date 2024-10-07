def get_qr(q_in):
    # input: q_in
    # output: qr
    # Ported from get_qr.pro by Fred Patt
    # Liang Hong, 2/19/2020
    
    import numpy as np
    from hawknav.qinv import qinv
    from hawknav.qprod import qprod
    
    nq = np.shape(q_in)[0]
    qi = qinv(q_in[0:nq-1,:])
    qr = qprod(qi,q_in[1:nq,:])
    
    j=np.where(qr[:,3] < 0)
    if (np.size(j)>0):
        qr[j,:] = -qr[j,:]
        
    return qr

