def remake_orbit_objects(norb,nc_fid,ngid):
    # Program to redefine orbit objects in Hawkeye L1A file
    # input: 
    #  - norb: number of orbits
    #  - nc_fid: netcdf file id
    #  - ngid: group id
    # Liang Hong, 2/26/2020, note: there is a library bug in netCDF4 module giving runtime error when renamevariable
    # Liang Hong, 2/24/2021, removed rpy_offset from navigation_data group attribute
        
    ## Rename existing data objects
    #ngid.renameVariable('orb_time','old_time')
    #ngid.renameVariable('orb_pos','old_pos')
    #ngid.renameVariable('orb_vel','old_vel')
    
    # set rol, pitch, yaw initial offsets to 0
    # nc_fid.groups['navigation_data']
    # ngid.rpy_offset = '0,0,0'
        
    # Create dimension
    if 'x' not in nc_fid.dimensions.keys():
        nc_fid.createDimension('x',norb)
        nc_fid.createDimension('y',3)

    # Create data objects
    orb_time_new = nc_fid.createVariable('/navigation_data/orb_time','f8',('x'),fill_value=-999)
    orb_pos_new = nc_fid.createVariable('/navigation_data/orb_pos',float,('x','y'),fill_value=-999999)
    orb_vel_new = nc_fid.createVariable('/navigation_data/orb_vel',float,('x','y'),fill_value=-999)

    # Set attributes 
    orb_time_new.long_name = "Orbit vector time (seconds of day)"
    orb_pos_new.long_name = "Orbit position vectors"
    orb_vel_new.long_name = "Orbit velocity vectors"
    orb_time_new.valid_min = 0.0
    orb_pos_new.valid_min = -7000.0
    orb_vel_new.valid_min = -8.0
    orb_time_new.valid_max = 86400.0
    orb_pos_new.valid_max = 7000.0
    orb_vel_new.valid_max = 8.0
    orb_time_new.units = "seconds"
    orb_pos_new.units = "kilometers"
    orb_vel_new.units = "kilometers/second"
        
    return
