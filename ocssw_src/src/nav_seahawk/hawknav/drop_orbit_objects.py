def drop_orbit_objects(l1afile,l1arenav):
    # drop orb_time,orb_pos and orb_vel from navigation_data group in the seahawk l1a file
    # in preparation of adding updated navigation data to same variables
    # Liang Hong, 3/11/2020
    # Liang Hong, 2/19/2021, added missing group attributes to destination file
    import netCDF4
    
    with netCDF4.Dataset(l1afile) as src, netCDF4.Dataset(l1arenav, "w") as dst:
        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
        # copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        # copy all file data except for the excluded
        for name, group in src.groups.items():
            dst.createGroup(name)
            dst[name].setncatts(src[name].__dict__)
            for vname,variable in group.variables.items():
                if vname.find('orb_')>-1:
                    continue
                gvname = name+'/'+vname
                x = dst.createVariable(gvname, variable.datatype, variable.dimensions)
                # copy variable attributes all at once via dictionary
                dst[gvname].setncatts(src[gvname].__dict__)
                dst[gvname][:] = src[gvname][:]
            