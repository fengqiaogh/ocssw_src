#!/usr/bin/env python3
# Amir Ibrahim - November 2017

# IMPORTS
from datetime import datetime
import time
import os
import re
import sys
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import requests
from netCDF4 import Dataset
import argparse
import json
import collections
from pathlib import Path
import tempfile

__version__ = '1.3.0_2020-12-11'
global verbose

# Define the sensor to generate the table

########################################
# TODO: use $OCDATAROOT/common/SensorInfo.json as starter dict
# read file

#sensorInfoFile = os.path.join(os.environ['OCDATAROOT'] , 'common','SensorInfo.json')
#with open(sensorInfoFile, 'r') as myfile:
#    sensorDefs=myfile.read()
#
## parse file
#sensors = json.loads(sensorDefs)
#print(sensors)

sensors = collections.defaultdict(dict)

sensors['modisa']['name'] = "MODIS-Aqua"
sensors['modisa']['rsr'] = "aqua_modis_RSR.nc"
sensors['modist']['name'] = "MODIS-Terra"
sensors['modist']['rsr'] = "terra_modis_RSR.nc"
sensors['seawifs']['name'] = "SeaWiFS"
sensors['seawifs']['rsr'] = "orbview-2_seawifs_RSR.nc"
sensors['meris']['name'] = "MERIS"
sensors['meris']['rsr'] = "envisat_meris_RSR.nc"
sensors['octs']['name'] = "OCTS"
sensors['octs']['rsr'] = "adeos_octs_RSR.nc"
sensors['czcs']['name'] = "CZCS"
sensors['czcs']['rsr'] = "nimbus-7_czcs_RSR.nc"
sensors['viirsn']['name'] = "VIIRS-SNPP"
sensors['viirsn']['rsr'] = "suomi-npp_viirs_RSR.nc"
sensors['viirsj1']['name'] = "VIIRS-JPSS1"
sensors['viirsj1']['rsr'] = "jpss-1_viirs_RSR.nc"
sensors['aviris']['name'] = "AVIRIS"
sensors['aviris']['rsr'] = "aviris_RSR.nc"
sensors['oci']['name'] = "OCI"
sensors['oci']['rsr'] = "pace_oci_RSR.nc"
sensors['s3aolci']['name'] = "OLCI-S3A"
sensors['s3aolci']['rsr'] = "sentinel-3a_olci_RSR.nc"
sensors['s3bolci']['name'] = "OLCI-S3B"
sensors['s3bolci']['rsr'] = "sentinel-3b_olci_RSR.nc"
sensors['hico']['name'] = "HICO"
sensors['hico']['rsr'] = "iss_hico_RSR.nc"
sensors['oli']['name'] = "Landsat-8 OLI"
sensors['oli']['rsr'] = "landsat-8_oli_RSR.nc"
sensors['goci']['name'] = "GOCI"
sensors['goci']['rsr'] = "coms_goci_RSR.nc"
sensors['ocm1']['name'] = "OCM-1"
sensors['ocm1']['rsr'] = "irs-p4_ocm_RSR.nc"
sensors['ocm2']['name'] = "OCM-2"
sensors['ocm2']['rsr'] = "oceansat-2_ocm-2_RSR.nc"
sensors['osmi']['name'] = "OSMI"
sensors['osmi']['rsr'] = "kompsat_osmi_RSR.nc"
sensors['misr']['name'] = "MISR"
sensors['misr']['rsr'] = "terra_misr_RSR.nc"

########################################


def gaussian(x, mu, sig):
    ''' Calculate gaussian filter
        Input: wavelength, band cetner, FWHM
        Output: Guassian filter
    '''
    sig = sig/2.355
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def get_hico_Rsr (all_wl):
    ''' Calculate HICO RSR - for testing and comparing to operational code
        Input: Wavelength (nm) at which the RSR will be calculates,
               for this exercise, use the water vapor tx wavelength
        Output: hico band center, RSR
    '''
    hico_wl=np.array([352.528,444.176,535.824,627.472,719.12,810.768,902.416,994.064,
    358.256,449.904,541.552,633.2,724.848,816.496,908.144,999.79,
    363.984,455.632,547.28,638.928,730.576,822.224,913.872,1005.52,
    369.712,461.36,553.008,644.656,736.304,827.952,919.6,1011.25,
    375.44,467.088,558.736,650.384,742.032,833.68,925.328,1016.98,
    381.168,472.816,564.464,656.112,747.76,839.408,931.056,1022.7,
    386.896,478.544,570.192,661.84,753.488,845.136,936.784,1028.43,
    392.624,484.272,575.92,667.568,759.216,850.864,942.512,1034.16,
    398.352,490,581.648,673.296,764.944,856.592,948.24,1039.89,
    404.08,495.728,587.376,679.024,770.672,862.32,953.968,1045.62,
    409.808,501.456,593.104,684.752,776.4,868.048,959.696,1051.34,
    415.536,507.184,598.832,690.48,782.128,873.776,965.424,1057.07,
    421.264,512.912,604.56,696.208,787.856,879.504,971.152,1062.8,
    426.992,518.64,610.288,701.936,793.584,885.232,976.88,1068.53,
    432.72,524.368,616.016,707.664,799.312,890.96,982.608,1074.26,
    438.448,530.096,621.744,713.392,805.04,896.688,988.336,1079.98])
    hico_wl = np.sort(hico_wl[:,np.newaxis],axis=0)
    fwhm = []
    for wl in hico_wl:
        if(wl<745):
            fwhm = np.append(fwhm,10)
        elif(wl>745):
            fwhm = np.append(fwhm,20)
    fwhm = fwhm[:,np.newaxis]
    Rsr = gaussian(all_wl, hico_wl, fwhm)
    #Rsr = Rsr.transpose()
    hico_wl = hico_wl.ravel()
    return hico_wl, Rsr

def read_sensor_RSR (sensor):
    ''' Reads a sensor spectral response (RSR) from the OCW through https
        Input: Sensor name (string)
        Output: normalized Rsr (numpy 2-d array - wavelength x sensor band),
                wavelength,
                sensor band labels,
    '''
    global verbose

    defpath = 'https://oceancolor.gsfc.nasa.gov/docs/rsr/'
    webpath = 'undefined'

    try:
        webpath = defpath + sensors[sensor]['rsr']
        if verbose:
            print("Reading spectral response function from %s" % webpath)
    except:
        print("Unkown sensor")
        sys.exit()

    tmpdir = tempfile.TemporaryDirectory()
    tmpfile = Path(tmpdir.name + '/' + sensors[sensor]['rsr'])
    r = requests.get(webpath, stream=True)

    if (r.status_code == requests.codes.ok):
        with open(tmpfile, 'wb') as fd:
            for chunk in r.iter_content(chunk_size=128):
                fd.write(chunk)

        ncid = Dataset(tmpfile,'r')
        wavelength = ncid['wavelength'][:]
        Rsr = ncid['RSR'][:][:]
        labels = ncid['bands'][:]
        return wavelength, Rsr, labels
    else:
        print("Received error: %s" % r.status_code)
        print("while attempting to retrieve %s" % webpath)
        sys.exit()


def read_TPVMR ():
    ''' Reads temprature, pressure, and vmr profiles from nc file in /common for
        models: 0:Tropical,1:Mid Latitude Summer,2:Mid Latitude Winter,
        3:Subarctic Summer,4:Subarctic Winter,5:US Standard 1962,6:User Defined Model
        Input: None
        Output: dp, vmrm
    '''
    abscf_path = os.environ.get('OCSSWROOT') + '/share/common/'
    tpvmr = Dataset(abscf_path + 'atrem_tpvmr.nc')
    p = np.array(tpvmr.variables['p'])
    #t = np.array(tpvmr.variables['t'])
    vmr = np.array(tpvmr.variables['vmr'])
    vmr = vmr.transpose()*1e-6
    p = p.transpose()/1013
    vmrm= np.empty([len(vmr)-1,7])
    DP= np.empty([len(vmr)-1,7])
    for i in range(0,len(vmr)-1):
        vmrm[i] = 0.5*(vmr[i] + vmr[i+1])
        DP[i] = p[i] - p[i+1]
    return DP, vmrm, p

def read_gas_abscf ():
    ''' Reads gas (water vapor) absorption coeff. from a nc file in /common
        Input: None
        Output: waveno, A
    '''
    abscf_path = os.environ.get('OCSSWROOT') + '/share/common/'
    a = Dataset(abscf_path + 'abscf_gas.nc')
    waveno = np.array(a.variables['waveno'])
    A = np.array(a.variables['abscf_h2o'])
    return waveno, A

def read_gas_abscf_co2 ():
    ''' Reads gas (CO2) absorption coeff. from a nc file in /common
        Input: None
        Output: waveno, A
    '''
    abscf_path = os.environ.get('OCSSWROOT') + '/share/common/'
    a = Dataset(abscf_path + 'abscf_gas.nc')
    waveno = np.array(a.variables['waveno'])
    A = np.array(a.variables['abscf_co2'])
    return waveno, A

def read_gas_abscf_o2 ():
    ''' Reads gas (O2) absorption coeff. from a nc file in /common
        Input: None
        Output: waveno, A
    '''
    abscf_path = os.environ.get('OCSSWROOT') + '/share/common/'
    a = Dataset(abscf_path + 'abscf_gas.nc')
    waveno = np.array(a.variables['waveno'])
    A = np.array(a.variables['abscf_o2'])
    return waveno, A

def read_gas_abscf_n2o ():
    ''' Reads gas (N2O) absorption coeff. from a nc file in /common
        Input: None
        Output: waveno, A
    '''
    abscf_path = os.environ.get('OCSSWROOT') + '/share/common/'
    a = Dataset(abscf_path + 'abscf_gas.nc')
    waveno = np.array(a.variables['waveno'])
    A = np.array(a.variables['abscf_n2o'])
    return waveno, A

def read_gas_abscf_co ():
    ''' Reads gas (CO) absorption coeff. from a nc file in /common
        Input: None
        Output: waveno, A
    '''
    abscf_path = os.environ.get('OCSSWROOT') + '/share/common/'
    a = Dataset(abscf_path + 'abscf_gas.nc')
    waveno = np.array(a.variables['waveno'])
    A = np.array(a.variables['abscf_co'])
    return waveno, A

def read_gas_abscf_ch4 ():
    ''' Reads gas (CH4) absorption coeff. from a nc file in /common
        Input: None
        Output: waveno, A
    '''
    abscf_path = os.environ.get('OCSSWROOT') + '/share/common/'
    a = Dataset(abscf_path + 'abscf_gas.nc')
    waveno = np.array(a.variables['waveno'])
    A = np.array(a.variables['abscf_ch4'])
    return waveno, A

def LBL_trans_H2O (waveno, A, dp, vmrm, wv):
    ''' Calculate LBL water vapor transmittance
        Input: wave number, absorption coefficients, delta-pressure,
        volume mixing ratio (ppm), water vapor amount (cm)
        Output: LBL water vapor transmittance
    '''
    T = np.empty([int(3e5),len(wv),6])
    #wv = wv[:,np.newaxis].transpose()
    Q = 2.15199993E+25
    for i in range(6):
        DP = dp[:,i]
        DP = DP[:,np.newaxis]
        VMRM = vmrm[:,i]
        VMRM = VMRM[:,np.newaxis]
        C = (Q*28.966 / 6.0225e23 / 1e-6)*DP*VMRM
        C = C/ (np.sum(Q*VMRM*DP)/3.34e22)   #----------> Normalization
        for w in range(len(wv)):
            T[:,w,i] = np.prod(np.exp(-A*np.reshape(C, [len(C), 1]) * wv[w]), axis=0)
    return T

def LBL_trans_O2 (waveno, A):
    ''' Calculate LBL O2 transmittance
        Input: wave number, absorption coefficients, delta-pressure,
        volume mixing ratio (ppm), water vapor amount (cm)
        Output: LBL water vapor transmittance
    '''
    abscf_path = os.environ.get('OCSSWROOT') + '/share/common/'
    tpvmr = Dataset(abscf_path + 'atrem_tpvmr.nc')
    p = np.array(tpvmr.variables['p'])
    pp = p[5,:]/1013
    DP = np.empty(19)
    for i in range(19):
        DP[i] = pp[i] - pp[i+1]
    T = np.empty([int(3e5)])
    #wv = wv[:,np.newaxis].transpose()
    Q = 2.15199993E+25
    #CO2 = np.linspace(340,400,10)
    VMRM = (0.21)
    C = (Q*28.966) / 6.0225e23 / 1e-6
    A[:,9001:106600] = A[:,9001:106600]*2.6
    T = np.prod(np.exp(-A*VMRM*C*np.reshape(DP,[-1,1])), axis=0)
    return T

def LBL_trans_N2O (waveno, A):
    ''' Calculate LBL N2O transmittance
        Input: wave number, absorption coefficients, delta-pressure,
        volume mixing ratio (ppm), water vapor amount (cm)
        Output: LBL water vapor transmittance
    '''
    abscf_path = os.environ.get('OCSSWROOT') + '/share/common/'
    tpvmr = Dataset(abscf_path + 'atrem_tpvmr.nc')
    p = np.array(tpvmr.variables['p'])
    pp = p[5,:]/1013
    DP = np.empty(19)
    for i in range(19):
        DP[i] = pp[i] - pp[i+1]
    T = np.empty([int(3e5)])
    #wv = wv[:,np.newaxis].transpose()
    Q = 2.15199993E+25
    #CO2 = np.linspace(340,400,10)
    VMRM = 0.3*1.0E-06
    C = (Q*28.966) / 6.0225e23 / 1e-6
    T = np.prod(np.exp(-A*VMRM*C*np.reshape(DP,[-1,1])), axis=0)
    return T

def LBL_trans_CO (waveno, A):
    ''' Calculate LBL CO transmittance
        Input: wave number, absorption coefficients, delta-pressure,
        volume mixing ratio (ppm), water vapor amount (cm)
        Output: LBL water vapor transmittance
    '''
    abscf_path = os.environ.get('OCSSWROOT') + '/share/common/'
    tpvmr = Dataset(abscf_path + 'atrem_tpvmr.nc')
    p = np.array(tpvmr.variables['p'])
    pp = p[5,:]/1013
    DP = np.empty(19)
    for i in range(19):
        DP[i] = pp[i] - pp[i+1]
    T = np.empty([int(3e5)])
    #wv = wv[:,np.newaxis].transpose()
    Q = 2.15199993E+25
    #CO2 = np.linspace(340,400,10)
    VMRM = 0.1*1.0E-06
    C = (Q*28.966) / 6.0225e23 / 1e-6
    T = np.prod(np.exp(-A*VMRM*C*np.reshape(DP,[-1,1])), axis=0)
    return T

def LBL_trans_CO2 (waveno, A, CO2):
    ''' Calculate LBL CO2 transmittance
        Input: wave number, absorption coefficients, delta-pressure,
        volume mixing ratio (ppm), water vapor amount (cm)
        Output: LBL water vapor transmittance
    '''
    abscf_path = os.environ.get('OCSSWROOT') + '/share/common/'
    tpvmr = Dataset(abscf_path + 'atrem_tpvmr.nc')
    p = np.array(tpvmr.variables['p'])
    pp = p[5,:]/1013
    DP = np.empty(19)
    for i in range(19):
        DP[i] = pp[i] - pp[i+1]
    T = np.empty([int(3e5)])
    #wv = wv[:,np.newaxis].transpose()
    Q = 2.15199993E+25
    #CO2 = np.linspace(340,400,10)
    VMRM = CO2*1.0E-06
    C = (Q*28.966) / 6.0225e23 / 1e-6
    T = np.prod(np.exp(-A*VMRM*C*np.reshape(DP,[-1,1])), axis=0)
    return T

def LBL_trans_CH4 (waveno, A):
    ''' Calculate LBL CH4 transmittance
        Input: wave number, absorption coefficients, delta-pressure,
        volume mixing ratio (ppm), water vapor amount (cm)
        Output: LBL water vapor transmittance
    '''
    abscf_path = os.environ.get('OCSSWROOT') + '/share/common/'
    tpvmr = Dataset(abscf_path + 'atrem_tpvmr.nc')
    p = np.array(tpvmr.variables['p'])
    pp = p[5,:]/1013
    DP = np.empty(19)
    for i in range(19):
        DP[i] = pp[i] - pp[i+1]
    T = np.empty([int(3e5)])
    #wv = wv[:,np.newaxis].transpose()
    Q = 2.15199993E+25
    #CO2 = np.linspace(340,400,10)
    VMRM = 1.8*1.0E-06
    C = (Q*28.966) / 6.0225e23 / 1e-6
    T = np.prod(np.exp(-A*VMRM*C*np.reshape(DP,[-1,1])), axis=0)
    return T

def calc_wv_trans_sensor (waveno, T, wavelength, Rsr, bands, wv):
    ''' Calculate the water vapor transmittance within a sensor band after applying sensor RSR
        Input: wave number, LBL transmittance, RSR wavelength, RSR, band centers, water vapor
        Output: water vapor transmittance per sensor band, Interpolated RSR to LBL wavelength,
                LBL wavelength, testing variable (to be removed)
        '''
    F0_path = os.environ.get('OCSSWROOT') + '/share/common/'
    F0_file  = F0_path + 'Thuillier_F0.dat'
    F0 = np.array(pd.read_csv(F0_file, skiprows= 15, delimiter='\s', names=['λ','F0'], engine='python'))
    trans_sensor = np.empty([6,len(bands),len(wv)])
    Rsr_f = interp1d(wavelength,Rsr,kind='nearest',fill_value="extrapolate")
    F0_f = interp1d(F0[:,0],F0[:,1],kind='nearest',fill_value="extrapolate")
    all_wl = np.append(1e7/waveno,np.arange(np.min(1e7/waveno)-0.1,299.91,-0.1))
    Rsr_int = Rsr_f(all_wl)
    F0_int = F0_f(all_wl)
    for i in range(6):
        for j,_ in enumerate(bands):
            for k, _ in enumerate(wv):
                trans_sensor[i,j,k] = np.trapz((all_wl),
                                               np.append(T[:,k,i],np.repeat(1,len(all_wl)-len(waveno)))
                                               *Rsr_int[j,:]*F0_int)/np.trapz(all_wl,F0_int*Rsr_int[j,:])
    z=np.append(T[:,k,i],np.repeat(1,len(all_wl)-len(waveno)))
    return trans_sensor, Rsr_int, all_wl, z

def calc_trans_sensor (waveno, T, wavelength, Rsr, bands):
    ''' Calculate the  transmittance within a sensor band after applying sensor RSR
        Input: wave number, LBL transmittance, RSR wavelength, RSR, band centers, water vapor
        Output: water vapor transmittance per sensor band, Interpolated RSR to LBL wavelength,
                LBL wavelength, testing variable (to be removed)
        '''
    F0_path = os.environ.get('OCSSWROOT') + '/share/common/'
    F0_file  = F0_path + 'Thuillier_F0.dat'
    F0 = np.array(pd.read_csv(F0_file, skiprows= 15, delimiter='\s', names=['λ','F0'], engine='python'))
    trans_sensor = np.empty([len(bands)])
    Rsr_f = interp1d(wavelength,Rsr,kind='nearest',fill_value="extrapolate")
    F0_f = interp1d(F0[:,0],F0[:,1],kind='nearest',fill_value="extrapolate")
    all_wl = np.append(1e7/waveno,np.arange(np.min(1e7/waveno)-0.1,299.91,-0.1))
    Rsr_int = Rsr_f(all_wl)
    F0_int = F0_f(all_wl)
    for i,_ in enumerate(bands):
        trans_sensor[i] = np.trapz((all_wl),
            np.append(T[:],np.repeat(1,len(all_wl)-len(waveno)))
            *Rsr_int[i,:]*F0_int)/np.trapz(all_wl,F0_int*Rsr_int[i,:])

    return trans_sensor

def create_nc(sensor, ofile, waterVapor, bands, models, history=None):
    """
    This function opens and defines products for output as a netcdf file
    """

    if verbose:
        print("Writing netCDF file: %s" % ofile)
    # Create the netCDF file and add dimensions
    dataset = Dataset(ofile, 'w',format='NETCDF4')
    dataset.createDimension('n_water_vapor', len(waterVapor))
    dataset.createDimension('nmodels', len(models))
    dataset.createDimension('nwavelengths', len(bands))

    # Set global attributes
    dataset.title = "Atmospheric gas transmittance table"
    dataset.description = "Atmospheric gas transmittance table for " + sensors[sensor]['name'] + " generated for a US Standard Atmosphere"
    dataset.date_created = (datetime.fromtimestamp(time.time()).strftime('%Y-%m-%dT%H:%M:%SZ')).encode('ascii')
    dataset.creator_name = "Amir Ibrahim, NASA/GSFC/OBPG"
    dataset.creator_email = "data@oceancolor.gsfc.nasa.gov"
    dataset.creator_url = "https://oceancolor.gsfc.nasa.gov"
    dataset.product_name=ofile.encode("ascii")
    dataset.history=history.encode("ascii")
#    dataset.source = 'generate_gas_tables.py'

    # Create the table data sets
    waterVaporSDS = dataset.createVariable('water_vapor', np.float32, ('n_water_vapor',))
    modelSDS = dataset.createVariable('model', str, 'nmodels')
    bandSDS = dataset.createVariable('wavelength', np.float32, ('nwavelengths',))
    wvTransSDS = dataset.createVariable('water_vapor_transmittance', np.float32,('nmodels','nwavelengths','n_water_vapor'))
    co2TransSDS = dataset.createVariable('carbon_dioxide_transmittance', np.float32,('nwavelengths'))
    coTransSDS = dataset.createVariable('carbon_monoxide_transmittance', np.float32,('nwavelengths'))
    ch4TransSDS = dataset.createVariable('methane_transmittance', np.float32,('nwavelengths'))
    n2oTransSDS = dataset.createVariable('nitrous_oxide_transmittance', np.float32,('nwavelengths'))
    o2TransSDS = dataset.createVariable('oxygen_transmittance', np.float32,('nwavelengths'))
    
    # Add some variable attributes
    waterVaporSDS.long_name = 'water vapor concentration'
    waterVaporSDS.units = 'cm'
    modelSDS.long_name = 'US Standard atmosphere model'
    modelSDS.units = ''
    bandSDS.long_name = 'wavelength'
    bandSDS.units = 'nm'

    wvTransSDS.long_name = 'water vapor transmittance'
    wvTransSDS.valid_min = 0
    wvTransSDS.valid_max = 1

    co2TransSDS.long_name = 'carbon dioxide transmittance'
    co2TransSDS.valid_min = 0
    co2TransSDS.valid_max = 1

    coTransSDS.long_name = 'carbon monoxide transmittance'
    coTransSDS.valid_min = 0
    coTransSDS.valid_max = 1

    ch4TransSDS.long_name = 'methane transmittance'
    ch4TransSDS.valid_min = 0
    ch4TransSDS.valid_max = 1

    n2oTransSDS.long_name = 'nitrous oxide transmittance'
    n2oTransSDS.valid_min = 0
    n2oTransSDS.valid_max = 1

    o2TransSDS.long_name = 'oxygen transmittance'
    o2TransSDS.valid_min = 0
    o2TransSDS.valid_max = 1

    # Write the static data
    waterVaporSDS[:] = waterVapor
    modelSDS[:] = models
    bandSDS[:] = bands

    return dataset

def write_transmittance(ncid,product,transmittance):
    ncid[product][:] = transmittance

def close_nc(ncid):
    # Close the file
    ncid.close()

def main():
    """
    Primary driver of the program; get command line arguments, check the files
    specified and kick off the processing
    """
    global verbose
    parser = argparse.ArgumentParser(description=\
        'Generates gas transmittance tables and writes them into a netCDF file')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('sensor', type=str,
                        choices=['modisa',
                                 'modist',
                                 'viirsn',
                                 'viirsj1',
                                 'seawifs',
                                 'meris',
                                 's3aolci',
                                 's3bolci',
                                 'hico',
                                 'oli',
                                 'goci',
                                 'ocm1',
                                 'ocm2',
                                 'osmi',
                                 'aviris',
                                 'octs',
                                 'czcs',
                                 'oci',
                                 'misr'],
                        default='modis-aqua', help='Sensor for which to generate the table')
    parser.add_argument('--output_file', type=str, default='gas_transmittance.nc', help='output netCDF LUT filename; default is <SENSOR>_gas_transmittance.nc')
    parser.add_argument('--verbose', '-v', action='store_true')

    args = parser.parse_args()

    verbose = args.verbose

    ofile = args.output_file
    if ofile == 'gas_transmittance.nc':
        ofile = '_'.join([args.sensor,ofile])

    if verbose:
        print("Generating %s ..." % ofile)

    # these are the water vapor values at which the transmittance are calculated (From Bo-Cai LOWTRAN tables)
    wv = np.array([0,0.002,0.004,0.008,0.012,0.016,0.02,0.024,0.028,0.032,0.036,0.04,0.08,0.12,0.16,0.2,0.24,0.28,
    0.32,0.36,0.4,0.44,0.48,0.52,0.56,0.6,0.64,0.68,0.72,0.76,0.8,0.84,0.88,0.92,0.96,1,1.04,1.08,
    1.12,1.16,1.2,1.24,1.28,1.32,1.36,1.4,1.44,1.48,1.52,1.56,1.6,1.64,1.68,1.72,1.76,1.8,1.84,1.88,
    1.92,1.96,2,2.04,2.08,2.12,2.16,2.2,2.24,2.28,2.32,2.36,2.4,2.44,2.48,2.52,2.56,2.6,2.64,2.68,
    2.72,2.76,2.8,2.84,2.88,2.92,2.96,3,3.04,3.08,3.12,3.16,3.2,3.24,3.28,3.32,3.36,3.4,3.44,3.48,
    3.52,3.56,3.6,3.64,3.68,3.72,3.76,3.8,3.84,3.88,3.92,3.96,4,4.2,4.4,4.6,4.8,5,5.2,5.4,5.6,5.8,
    6,6.2,6.4,6.6,6.8,7,7.2,7.4,7.6,7.8,8,8.4,8.8,9.2,9.6,10,10.4,10.8,11.2,11.6,12,12.4,12.8,13.2,
    13.6,14,14.4,14.8,15.2,15.6,16,16.4,16.8,17.2,17.6,18,18.4,18.8,19.2,19.6,20,20.4,20.8,21.2,
    21.6,22,22.4,22.8,23.2,23.6,24,24.8,25.6,26.4,27.2,28,28.8,29.6,30.4,31.2,32,32.8,33.6,34.4,
    35.2,36,36.8,37.6,38.4,39.2,40,42,44,46,48,50,52,54,56,58,60,64,68,72,76,80,88,96,104,112,120,
    128,136,144,152,160,168,176,188,200])
    #wv = np.array([1,2])

    # These are the standard models
    models = np.array(['Tropical','Mid Latitude Summer','Mid Latitude Winter','Subarctic Summer',
              'Subarctic Winter','US Standard 1962'],dtype='str')

    # Load sensor RSR
    wavelength, Rsr, bands = read_sensor_RSR(args.sensor)

    # Open a netCDF file that will contain the per sensor water vapor transmittance
    history = ''
    history = "{} {}".format(parser.prog,args.sensor)
    if args.output_file != 'gas_transmittance.nc':
        history += " --output_file {}".format(ofile)

    ncid = create_nc(args.sensor, ofile, wv, bands, models, history=history)

    # water vapor transmittance
    # Load delta-pressure and volume mixing ratio of water vapor for all 6 standard models
    if verbose:
        print("Loading delta-pressure and volume mixing ratio of water vapor for all 6 standard models")
    dp,vmrm,p = read_TPVMR()

    # Load water vapor absorption coefficients and wave numbers
    if verbose:
        print("Loading water vapor absorption coefficients and wave numbers")
    waveno,A = read_gas_abscf()

    # Calculate the LBL water vapor transmittance
    if verbose:
        print("Calculating water vapor transmittance")
    T_wv = LBL_trans_H2O(waveno,A,dp,vmrm,wv)

    # Apply the sensor RSR on the LBL to get the water vapor transmittance within the sensor's band
    if verbose:
        print("Convolving water vapor transmittance with sensor specrtal response")
    sensor_trans_w, R,all_wl, z = calc_wv_trans_sensor(waveno, T_wv, wavelength, Rsr, bands, wv)

    if verbose:
        print("Writing water vapor transmittance")
    write_transmittance(ncid,'water_vapor_transmittance',sensor_trans_w)

    # CO2 transmittance
    if verbose:
        print("Loading carbon dixoide absorption coefficients and wave numbers")
    waveno, A_co2 = read_gas_abscf_co2()

    if verbose:
        print("Calculating carbon dioxide transmittance")
    T_co2 = LBL_trans_CO2(waveno,A_co2, 400)

    if verbose:
        print("Convolving carbon dioxide transmittance with sensor specrtal response")
    sensor_trans_co2 = calc_trans_sensor(waveno, T_co2, wavelength, Rsr, bands)

    if verbose:
        print("Writing carbon dioxide transmittance")
    write_transmittance(ncid,'carbon_dioxide_transmittance',sensor_trans_co2)

    # O2 transmittance
    if verbose:
        print("Loading oxygen absorption coefficients and wave numbers")
    waveno, A_o2 = read_gas_abscf_o2()

    if verbose:
        print("Calculating oxygen transmittance")
    T_o2 = LBL_trans_O2(waveno,A_o2)

    if verbose:
        print("Convolving oxygen transmittance with sensor specrtal response")
    sensor_trans_o2 = calc_trans_sensor(waveno, T_o2, wavelength, Rsr, bands)

    if verbose:
        print("Writing oxygen transmittance")
    write_transmittance(ncid,'oxygen_transmittance',sensor_trans_o2)

    # NO2 transmittance
    if verbose:
        print("Loading nitrous oxide absorption coefficients and wave numbers")
    waveno, A_n2o = read_gas_abscf_n2o()

    if verbose:
        print("Calculating nitrous oxide transmittance")
    T_n2o = LBL_trans_N2O(waveno,A_n2o)

    if verbose:
        print("Convolving nitrous oxide transmittance with sensor specrtal response")
    sensor_trans_no2 = calc_trans_sensor(waveno, T_n2o, wavelength, Rsr, bands)

    if verbose:
        print("Writing nitrous oxide transmittance")
    write_transmittance(ncid,'nitrous_oxide_transmittance',sensor_trans_no2)

    # CO transmittance
    if verbose:
        print("Loading carbon monoxide absorption coefficients and wave numbers")
    waveno, A_co = read_gas_abscf_co()

    if verbose:
        print("Calculating carbon monoxide transmittance")
    T_co = LBL_trans_CO(waveno,A_co)
 
    if verbose:
        print("Convolving carbon monoxide transmittance with sensor specrtal response")
    sensor_trans_co = calc_trans_sensor(waveno, T_co, wavelength, Rsr, bands)

    if verbose:
        print("Writing carbon monoxide transmittance")
    write_transmittance(ncid,'carbon_monoxide_transmittance',sensor_trans_co)

    # CH4 transmittance
    if verbose:
        print("Loading methane absorption coefficients and wave numbers")
    waveno, A_ch4 = read_gas_abscf_ch4()

    if verbose:
        print("Calculating methane transmittance")
    T_ch4 = LBL_trans_CH4(waveno,A_ch4)

    if verbose:
        print("Convolving methane transmittance with sensor specrtal response")
    sensor_trans_ch4 = calc_trans_sensor(waveno, T_ch4, wavelength, Rsr, bands)

    if verbose:
        print("Writing methane transmittance")
    write_transmittance(ncid,'methane_transmittance',sensor_trans_ch4)

    close_nc(ncid)
    if verbose:
        print("Processing complete!")

# The following allows the file to be imported without immediately executing.
if __name__ == '__main__':
    sys.exit(main())
