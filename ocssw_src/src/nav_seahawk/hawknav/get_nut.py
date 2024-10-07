def get_nut(iyr,day):
    #  Program to get nutation matrix for TOD-to-MOD conversion
    
    #  Calling Arguments:
    
    #  Name         Type    I/O     Description
    #  --------     ----    ---     -----------
    #  iyr          I*4      I      Year
    #  day         I*4      I      Day of year
    #  xnut(3,3)    R*8      O      Nutation matrix
    # ported from IDL (Fred Patt), Liang Hong, 9/20/2019
    # Liang Hong, 3/24/2020, array calculation
    
    import numpy as np
    from hawknav.jd import jd
    from hawknav.gha2000 import ephparms,nutate
    
    xnut = np.zeros((np.size(iyr),3,3))
    imon = 1
    iday = np.round(day)
    t = jd(iyr,imon,iday) - 2451545.50 + day - iday
    [xls,gs,xlm,omega] = ephparms(t)
    [dpsi,eps,epsm] = nutate(t,xls,gs,xlm,omega)
    
    dpsi_rad = np.deg2rad(dpsi)
    epsm_rad = np.deg2rad(epsm)
    eps_rad = np.deg2rad(eps)
        
    xnut[:,0,0] = np.cos(dpsi_rad)
    xnut[:,0,1] = -np.sin(dpsi_rad)*np.cos(epsm_rad)
    xnut[:,0,2] = -np.sin(dpsi_rad)*np.sin(epsm_rad)
    xnut[:,1,0] = np.sin(dpsi_rad)*np.cos(eps_rad)
    xnut[:,1,1] = np.cos(dpsi_rad)*np.cos(eps_rad)*np.cos(epsm_rad) \
            + np.sin(eps_rad)*np.sin(epsm_rad)
    xnut[:,1,2] = np.cos(dpsi_rad)*np.cos(eps_rad)*np.sin(epsm_rad) \
            - np.sin(eps_rad)*np.cos(epsm_rad)
    xnut[:,2,0] = np.sin(dpsi_rad)*np.sin(eps_rad)
    xnut[:,2,1] = np.cos(dpsi_rad)*np.sin(eps_rad)*np.cos(epsm_rad) \
            - np.cos(eps_rad)*np.sin(epsm_rad)
    xnut[:,2,2] = np.cos(dpsi_rad)*np.sin(eps_rad)*np.sin(epsm_rad) \
            + np.cos(eps_rad)*np.cos(epsm_rad)
    
    return xnut
