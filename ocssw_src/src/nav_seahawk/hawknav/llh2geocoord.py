def llh2geocoord(iyr,iday,lon,lat,hgt):
    # output: geomat
    # format data as input to geomag70.c 
    # Liang Hong, 2/24/2020
    
    from hawknav.jddate import jddate
    from hawknav.jd import jd
    import numpy as np
    
    nl = np.size(lon)
    njday = jd(iyr,1,iday)
    iy,mon,idy = jddate(njday)
    
    # row:  iyr,mon,idy,fix(hgt(i)),lat(i),lon(i),format='(i4,2(",",i02)," D K",i3,2f9.3)'
    geomat = np.column_stack((iyr*np.ones(nl),mon*np.ones(nl),idy*np.ones(nl),np.ones(nl),np.fix(hgt),lat,lon))
    
    return geomat