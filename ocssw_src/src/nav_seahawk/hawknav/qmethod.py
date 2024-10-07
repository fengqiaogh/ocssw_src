def qmethod(vecb, veci, wgt):
    # Program to determine attitude using Davenport's Q-method
    # input: vecb, veci, wgt
    # output: qib
    # reference: qmethod.pro by Fred Patt
    # Liang Hong, 2/20/2020
    
    import numpy as np
    from numpy import linalg as LA
    
    w = np.copy(vecb)
    v = np.copy(veci)
    w = w*np.transpose([np.sqrt(wgt),np.sqrt(wgt),np.sqrt(wgt)])
    v = v*np.transpose([np.sqrt(wgt),np.sqrt(wgt),np.sqrt(wgt)])
    
    b = np.dot(np.transpose(v),w)
    s = np.transpose(b) + b
    z = [b[2,1]-b[1,2],b[0,2]-b[2,0],b[1,0]-b[0,1]]
    sig = np.sum(np.diagonal(b))
    
    k = np.zeros((4,4))
    k[0:3,0:3] = s
    for i in range(0,3):
        k[i,i] = k[i,i] - sig
    k[0:3,3] = z
    k[3,0:3] = z
    k[3,3] = sig
    
    eigval,evec = LA.eig(k)  # eigen vector differs from "eval = eigenql(k,eigenvectors = evec)" IDL output in signs
    tmpind = np.where(eigval == np.max(eigval))
    qib = np.squeeze(evec[:,tmpind])
    
    return qib
