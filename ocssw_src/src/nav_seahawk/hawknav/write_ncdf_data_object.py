def write_ncdf_data_object(ngid,dname,data,start=0):
    #  procedure to write data to a NetCDF data object

    #	Arguments
    #     
    #     	Name    Type 	  I/O 	 Description
    #     	----	---- 	  --- 	 -----------
    #       ngid     int        I     Group ID for data object
    #       dname   string     I     Data object name
    #       data               I     Data array to write
    # 	start   int(*)     I     Optional array of dimension offsets
    # Ported from write_ncdf_data_object.pro by Fred Patt
    # Liang Hong, 2/25/2020
    
    import numpy as np
    
    # Open data object
    did = ngid.variables[dname]
    orig_data = np.copy(did[:].data)
       
    # Check for start array
    if (start==0):
        if orig_data.ndim==1:
            did[0:np.size(data)]=data
        else:
            did[0:np.shape(data)[0],:]=data
    else:
        if orig_data.ndim==1:
            did[:] = np.stack(orig_data[0:start],data)
        else:
            did[0:np.shape(data)[0],:] = np.stack(orig_data[0:start,:],data)
    return