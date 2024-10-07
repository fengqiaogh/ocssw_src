def qinv(q):
    # input: q
    # output: qin
    # Liang Hong, 2/19/2020
    
    import numpy as np
    qin = np.copy(q)
    try:
        qin[:,0:3] = -qin[:,0:3]
    except:
        qin[0:3] = -qin[0:3]
    
    return qin
