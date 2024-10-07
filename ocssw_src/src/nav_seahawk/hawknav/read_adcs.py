def fss_read_bytes(datain):
    # read fine sun sensor bytes
    import numpy as np
    import struct
    
    datain = np.flip(datain,1)
    dataout = struct.unpack('<%df'%(np.size(datain)/4),np.ascontiguousarray(datain))
    dataout = np.flip(np.reshape(dataout,(-1,2)),1)
    
    return dataout

def read_adcs_from_l1a(l1afile):
    # Program to read ADCS data from a Hawkeye L1A file
    # Assumes file path/name in varible l1afile
    # Inputs:
    #  - l1afile
    # Outputs:
    #  - ptime,sunv,mag2,gyro,rw,stime,sunxp,sunzp,sunzn,atime,quat
    # Ported from renav_hawkeye.pro by Fred Patt.
    # Liang Hong, 2/12/2020
    # Liang Hong, 3/25/2020, set sun sensor zp/zn readings to 0 when readings stay unchanged
    # Liang Hong, 4/27/2020, return empty arrays when no valid values in L1A
    
    import numpy as np
    from functools import reduce
    from hawknav.get_ncdf_object import get_ncdf_object
    
    group = 'navigation_data'
    
    # find data indices with valid time tag
    atime = get_ncdf_object(l1afile, group, 'att_time')
    ptime = get_ncdf_object(l1afile, group, 'processed_sensor_time')
    stime = get_ncdf_object(l1afile, group, 'sensor_time')
    kq = np.where(atime >= 0)
    kp = np.where(ptime >= 0)
    ks = np.where(stime >= 0)
    ind_valid_ts = reduce(np.intersect1d, (kq,kp,ks))
    if np.size(ind_valid_ts)==0:
        # no valid readings from L1A file
        print("No valid navigation time tags.")
        return [],[],[],[],[],[],[],[],[],[],[]
        
    # Attitude quaternions     
    quat = get_ncdf_object(l1afile, group, 'att_quat')    
    atime = atime[ind_valid_ts]
    quat = np.squeeze(quat[ind_valid_ts,:])
    
    # Processed sensor telemetry    
    sunv = get_ncdf_object(l1afile, group, 'sun_vector') 
    mag1 = get_ncdf_object(l1afile, group, 'mag1')
    mag2 = get_ncdf_object(l1afile, group, 'mag2')
    gyro = get_ncdf_object(l1afile, group, 'gyro_rates')
    rw = get_ncdf_object(l1afile, group, 'rwheels')
    ptime = np.squeeze(ptime[ind_valid_ts])
    sunv = np.squeeze(sunv[ind_valid_ts,:])
    mag2 = np.squeeze(mag2[ind_valid_ts,:])
    gyro = np.squeeze(gyro[ind_valid_ts,:])
    rw = np.squeeze(rw[ind_valid_ts,:])
    
    # FSS angles and error code    
    sensor = get_ncdf_object(l1afile, group, 'sensor_bus_telemetry')
    stime = stime[ind_valid_ts]
    sensor = np.squeeze(sensor[ind_valid_ts,:])
    
    # Extract FSS angles
    nks = np.size(ind_valid_ts)
    fssxp = fss_read_bytes(np.copy(sensor[:,0:8]))
    fsszp = fss_read_bytes(np.copy(sensor[:,9:17]))
    fsszn = fss_read_bytes(np.copy(sensor[:,18:26]))
    fsserr = sensor[:,[8,17,26]]
    
    # Convert FSS angles to vectors
    sunxp = np.zeros((nks,3))
    tana = np.tan(np.deg2rad(fssxp[:,0]))
    tanb = np.tan(np.deg2rad(fssxp[:,1]))
    sunxp[:,0] = 1/np.sqrt(tana**2 + tanb**2 + 1.0)
    sunxp[:,1] = -tanb*sunxp[:,0]
    sunxp[:,2] = tana*sunxp[:,0]
    
    sunzp = np.zeros((nks,3))
    if len(np.argwhere(np.diff(fsszp[:,0])==0))/nks<0.5:
        # use sun sensor zp when <50% readings changing during the time, 3/25/2020
        tana = np.tan(np.deg2rad(fsszp[:,0]))
        tanb = np.tan(np.deg2rad(fsszp[:,1]))
        sunzp[:,2] = 1/np.sqrt(tana**2 + tanb**2 + 1.0)
        sunzp[:,1] = tanb*sunzp[:,2]
        sunzp[:,0] = tana*sunzp[:,2]
    
    sunzn = np.zeros((nks,3))
    if len(np.argwhere(np.diff(fsszn[:,0])==0))/nks<0.5:
        # use sun sensor zp when <50% readings changing during the time, 3/25/2020
        tana = np.tan(np.deg2rad(fsszn[:,0]))
        tanb = np.tan(np.deg2rad(fsszn[:,1]))
        sunzn[:,2] = 1/np.sqrt(tana**2 + tanb**2 + 1.0)
        sunzn[:,1] = (tana-tanb)*sunzn[:,2]/np.sqrt(2)
        sunzn[:,0] = (tana+tanb)*sunzn[:,2]/np.sqrt(2)
        sunzn[:,2] = -sunzn[:,2]
    
    return ptime,sunv,mag2,gyro,rw,stime,sunxp,sunzp,sunzn,atime,quat

