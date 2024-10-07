'''
Created on Nov 19, 2015

@author: rhealy
'''
import numpy as np

def bore_sight(dth):
    import math as m
    
    cx=m.cos(dth[0])
    sx=m.sin(dth[0])
    cy=m.cos(dth[1])
    sy=m.sin(dth[1])
    cz=m.cos(dth[2])
    sz=m.sin(dth[2])
    
#    print('bore_sight:',cx,sx,cy,sy,cz,sz)
    
    rot = np.zeros((3,3))

    rot[0,:]=[cy*cz,cy*sz,-sy]
    rot[1,:]=[sx*sy*cz-cx*sz,sx*sy*sz+cx*cz,cy*sx]
    rot[2,:]=[cx*sy*cz+sx*sz,cx*sy*sz-sx*cz,cy*cx]
    
    return rot

if __name__ == '__main__':
    rad2Deg = 180.0/np.pi
    deg2Rad = np.pi/180.0
    as2Deg = 1.0/3600.0
    as2Rad = as2Deg*deg2Rad
    
    r_hico_bs = np.zeros((3))
    r_hico_bs[0] = -0.9957
    r_hico_bs[1] = 0.0268
    r_hico_bs[2] = -0.0128
    rot = bore_sight(r_hico_bs*deg2Rad)
#    print('rot=',str(rot))
    r_hico_to_iss = np.zeros((3,3))
    r_hico_to_iss[0][0] = -1.0
    r_hico_to_iss[1][1] = -1.0
    r_hico_to_iss[2][2] = 1.0
#    print('r_iss:',r_hico_to_iss)
#    print(np.dot(r_hico_to_iss,rot))