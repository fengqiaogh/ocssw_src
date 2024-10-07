def j2000_to_ecr(iy,idy,sec,strUTCpole):
    # Get J2000 to ECEF transformation matrix
    # Arguments:
    # Name 		Type 	I/O 	Description
    # --------------------------------------------------------
    # iy      	 I   	 I  	Year, four digits
    # idy     	 I   	 I  	Day of year
    # sec     	R*8  	 I  	Seconds of day
    # ecmat(3,3)	 R 	 O  	J2000 to ECEF matrix
    # strUTCpole  string, I    Earth motion file, 'utcpole.dat'
    # Ported from j2000_to_ecr.pro by Fred Patt
    # Liang Hong, 2/14/2020
    # Liang Hong, 3/24/2020, array calculation
    
    from hawknav.j2000_to_mod import j2000_to_mod
    from hawknav.get_nut import get_nut
    from hawknav.get_ut1 import get_ut1
    from hawknav.gha2000 import gha2000
    import numpy as np
    
    daysec = 86400.0     # Second per day
    
    # Get transformation from J2000 to mean-of-date inertial
    
    j2mod = np.squeeze(j2000_to_mod(iy,idy,sec))
    
    # Get nutation and UT1-UTC (once per run)
    xnut = np.squeeze(get_nut(iy,idy))
    ut1utc = get_ut1(iy,idy,strUTCpole)
    
    # Compute Greenwich hour angle for time of day
    
    day = idy + (sec+ut1utc)/daysec
    gha = gha2000(iy,day)
    gham = np.zeros((np.size(iy),3,3))
    gha = np.deg2rad(gha)
    gham[:,0,0] = np.cos(gha)
    gham[:,1,1] = np.cos(gha)
    gham[:,2,2] = 1.0
    gham[:,1,0] = np.sin(gha)
    gham[:,0,1] = -np.sin(gha)
    
    # Combine all transformations
    ## IDL code: ecmat = gham#np.transpose(xnut)#j2mod
    if np.size(iy)==1:
        xnut = np.squeeze(xnut)
        gham = np.squeeze(gham)
        #ecmat = np.dot(j2mod,np.dot(np.transpose(xnut),gham))
        ecmat = np.linalg.multi_dot([j2mod,np.transpose(xnut),gham])
    else:
        ecmat = np.matmul(j2mod,np.matmul(xnut.swapaxes(1,2),gham))
    
    return ecmat
