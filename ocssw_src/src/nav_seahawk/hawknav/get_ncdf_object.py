def get_ncdf_object(file,group,object):
    # Program to read a data object from a NetCDF file
    # Liang Hong, 2/11/2020
    
    from netCDF4 import Dataset
    
    nc_fid = Dataset(file,'r')
    ngid = nc_fid.groups[group]
    data = ngid.variables[object][:].data
    nc_fid.close()
    
    return data 

