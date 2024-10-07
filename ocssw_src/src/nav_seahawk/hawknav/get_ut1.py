def get_ut1(iyr,idy,strUTCpole):
    # Program to get UT1-UTC from utcpole.dat file
    # input: iyr, idy, strUTCpole (Earth motion file, 'utcpole.dat')
    # output: ut1utc = UT1-UTC(s)
    # Ported from get_ut1.pro by Fred Patt
    # Liang Hong, 2/25/2020
    # Liang Hong, 3/24/2020, array calculation
    
    from hawknav.jd import jd
    import numpy as np
    
    mjd = jd(iyr,1,idy) - 2400000
    
    # line:  MJD	x(arc sec)	x error		y(arc sec)	y error		UT1-UTC(s)	UT error	qual
    utcpole = np.loadtxt(strUTCpole,skiprows=2, usecols=[0,5])
    
    # readf,ulun,ijd,xn,xe,yn,ye,ut1utc,ute
    #ut1utc = utcpole[np.where(utcpole[:,0]==mjd),1]
    tmpind = np.where(np.in1d(utcpole[:,0], mjd))[0]
    ut1utc = utcpole[tmpind,1]
    
    return ut1utc