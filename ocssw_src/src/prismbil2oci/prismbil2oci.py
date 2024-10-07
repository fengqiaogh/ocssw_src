#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
written by J.Scott on 2017/10/03 (joel.scott@nasa.gov)
'''
import sys
import argparse
import os
import logging
from collections import OrderedDict
from PRISM_module import rd_prism_hdr, rd_L1B_bil
import numpy as np
import operator
import re
import subprocess
from netCDF4 import Dataset
from datetime import datetime, timedelta

__version__ = '1.1'

#==========================================================================================================================================

def ParseCommandLine(args):
    '''
    Defines and parses command line arguments
    '''
    ocdataroot = os.getenv('OCDATAROOT')
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, 
        description=' Converts PRISM data from the CORAL project to netCDF proxy data in the format of PACE-OCI ', 
        add_help=True)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('-i', '--ifile', nargs=1, required=True, type=str, 
        help=' PRISM L1 prm*img.hdr file (required) ')
    parser.add_argument('-o', '--ofile', nargs=1, required=True, type=str, 
        help=' Output netCDF filename (required) ')
    parser.add_argument('--cdlfile', nargs=1, type=str, 
        default=os.path.join(ocdataroot, 'prism', 'OCIP_Level-1C_Data_Structure.cdl'), 
        help=' L1C format spec CDL file name ')
    parser.add_argument('-d', '--debug', action='store_true', default=False)
    parsedArgs=parser.parse_args(args)
    return parsedArgs

#==========================================================================================================================================

def SetLogger(pargs):
    '''
    Initiates a logger instance
    '''
    lgrName = 'prismbil2oci_%s_T_%s' % (datetime.date(datetime.now()), datetime.time(datetime.now()))
    lgr = logging.getLogger(lgrName)
    fmt = '%(message)s'
    if pargs.debug:
        level = logging.DEBUG
        fmt = '%(asctime)s - %(name)s - %(levelname)s -\
               [%(module)s..%(funcName)s..%(lineno)d] -\
               %(message)s'
        formatter = logging.Formatter(fmt)
        fh = logging.FileHandler('%s.log' % lgrName)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        lgr.addHandler(fh)

    level = logging.INFO
    formatter = logging.Formatter(fmt)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    lgr.addHandler(ch)
    lgr.setLevel(level)
    lgr.debug('Logger initialized')
    return lgr

#==========================================================================================================================================

def read_prism_L1(logger, fname_in):
    '''
    A subroutine to read PRISM L1 data
    '''
    L1_hdr   = OrderedDict()
    L1_data  = OrderedDict()

    #create list of hdr fnames
    ls_hdr = []
    ls_hdr.append(fname_in) #img.hdr

    if os.path.exists(fname_in.rsplit('_',1)[0] + '_loc_ort.hdr'):
        ls_hdr.append(fname_in.rsplit('_',1)[0] + '_loc_ort.hdr') #loc_ort.hdr
    else:
        ls_hdr.append(fname_in.rsplit('_',1)[0] + '_loc.hdr') #loc.hdr

    if os.path.exists(fname_in.rsplit('_',1)[0] + '_obs_ort.hdr'):
        ls_hdr.append(fname_in.rsplit('_',1)[0] + '_obs_ort.hdr') #obs_ort.hdr
    else:
        ls_hdr.append(fname_in.rsplit('_',1)[0] + '_obs.hdr') #obs.hdr

    #create list of bil fnames
    ls_bil = [fname.rsplit('.')[0] for fname in ls_hdr]

    for fname_hdr,fname_bil in zip(ls_hdr,ls_bil):
        if not os.path.exists(fname_hdr):
            logger.warning('Missing PRISM header file: {:}'.format(fname_hdr))

        if not os.path.exists(fname_bil):
            logger.warning('Missing PRISM BIL file: {:}'.format(fname_bil))

        logger.info('Reading HDR file ' + fname_hdr)
        [hdr_type, dict_hdr] = rd_prism_hdr(fname_hdr)

        #save all hdr content
        L1_hdr[hdr_type] = OrderedDict()
        for key in dict_hdr:
            L1_hdr[hdr_type][key] = dict_hdr[key]

        logger.info('Reading BIL file ' + fname_bil)
        dict_data = rd_L1B_bil(fname_bil, dict_hdr)
        L1_data.update(dict_data)

    return [L1_hdr, L1_data]

#==========================================================================================================================================

def subsample_spectral(band_data, band_weights):
    '''
    a subroutine to subsample via weighted averaging
    spectral data about a desired band center

    Inputs:    band_data    -- 3-dim numpy array of shape-order (number_of_bands(_to_average), number_of_scans, pixels_per_scan)
               band_weights -- list or 1-dim numpy array of length (number_of_bands(_to_average))

    Outputs:   data         -- 2-dim numpy array of shape-order (number_of_scans, pixels_per_scan)
               width        -- floating point value of sub-sampled FWHM value
    '''
    width = np.sum(band_weights)
    
    data_numerator = np.empty_like(band_data)
    for i,weight in zip(range(data_numerator.shape[0]), band_weights):
        data_numerator[i,:,:] = np.multiply(band_data[i,:,:], weight)
    
    data = np.divide(np.sum(data_numerator, axis=0), width)
    
    return [data, width]

#==========================================================================================================================================

def resample_prism2oci(bands_prism, widths_prism, lt_prism):
    '''
    A subroutine to resample PRISM data to OCI targets to create OCIP data
    '''

    #map bands and FWHMs from PRISM to OCI for 350:890nm
    index_ocip  = list(range(0,191))
    bands_ocip  = list(bands_prism[0:191])
    widths_ocip = list(widths_prism[0:191])

    #map bands and FWHMs from PRISM to OCI for 940nm
    lis_940bands  = list(range(202,215))
    band_array = np.ndarray((len(lis_940bands),1,1), dtype='f8')
    band_array[:,0,0] = [bands_prism[i] for i in lis_940bands]
    [data,width] = subsample_spectral(band_array, [widths_prism[i] for i in lis_940bands])
    bands_ocip.append(data[0,0])
    widths_ocip.append(width)

    #map bands and FWHMs from PRISM to OCI for 1038nm
    lis_1038bands = list(range(239,246))
    band_array = np.ndarray((len(lis_1038bands),1,1), dtype='f8')
    band_array[:,0,0] = [bands_prism[i] for i in lis_1038bands]
    [data,width] = subsample_spectral(band_array, [widths_prism[i] for i in lis_1038bands])
    bands_ocip.append(data[0,0])
    widths_ocip.append(width)

    #resample Lt
    lt_ocip = np.ndarray((len(bands_ocip),lt_prism.shape[1],lt_prism.shape[2]), dtype="f8")
    lt_ocip[0:len(index_ocip),:,:] = lt_prism[index_ocip,:,:]

    [data,width] = subsample_spectral(lt_prism[lis_940bands,:,:], [widths_prism[i] for i in lis_940bands])
    lt_ocip[-2,:,:] = data[:,:]

    [data,width] = subsample_spectral(lt_prism[lis_1038bands,:,:], [widths_prism[i] for i in lis_1038bands])
    lt_ocip[-1,:,:] = data[:,:]

    return [bands_ocip, widths_ocip, lt_ocip]

#==========================================================================================================================================

def read_OCIP_L1C_CDL(fname_in, num_scans, num_pix, num_bands):
    '''
    A subroutine to read a OCIP L1C CDL file and fill the proper dimensions
    '''
    with open(fname_in,'r') as fid:
        lines = fid.readlines()

    #define dims and chunksizes in CDL file
    lines = [re.sub('NUMSSCANS',  str(num_scans), line) for line in lines]
    lines = [re.sub('NUMSPIXELS', str(num_pix),   line) for line in lines]
    lines = [re.sub('NUMSBANDS',  str(num_bands), line) for line in lines]

    return lines

#==========================================================================================================================================

def PRISMmain(args):
    '''
    Coordinates calls to PRISM L1 BIL to L1C netCDF conversion process
    '''
    #parse input arguments and instantiate logger
    pArgs = ParseCommandLine(args)
    mainLogger = SetLogger(pArgs)
    mainLogger.info("prismbil2oci.py %s" % __version__)
    ifile = ''.join(pArgs.ifile)
    ofile = ''.join(pArgs.ofile)
    cfile = ''.join(pArgs.cdlfile)

    #load L1 data
    [dsL1_hdr, dsL1_data] = read_prism_L1(mainLogger, ifile)

    #resample PRISM band centers and widths to more closely match OCI
    [bands_ocip, widths_ocip, lt_ocip] = resample_prism2oci([float(d) for d in dsL1_hdr['rdn_img']['wavelength']], 
                                                            [float(d) for d in dsL1_hdr['rdn_img']['fwhm'][:len(dsL1_hdr['rdn_img']['wavelength'])]], 
                                                            dsL1_data['lt'])

    #read and parse CDL file
    cdl_lines = read_OCIP_L1C_CDL(cfile, dsL1_data['utc time'].shape[0], dsL1_data['utc time'].shape[1], len(bands_ocip))

    #mk netCDF4 file via NCGEN subprocess call
    if '/' in ofile:
        nc_fname = ofile.rsplit('/',1)[-1]
    else:
        nc_fname = ofile

    mainLogger.info('Creating netCDF4 file {:}'.format(ofile))
    pid = subprocess.run('ncgen -b -v3 -o ' + ofile, input=''.join(cdl_lines), encoding='ascii', 
                         shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if pid.stderr:
        logger.warning('Error in ncgen call: {:}'.format(pid.stderr))

    nc_fid = Dataset(ofile, 'a')


    #fill global attributes
    mainLogger.info('Filling global attributes')
    nc_fid.product_name = nc_fname
    nc_fid.prism_version = dsL1_hdr['rdn_img']['prism_version']
    dt_now = datetime.now()
    nc_fid.date_created = dt_now.strftime("%Y-%m-%dT%H:%M:%SZ")
    dt_start = dsL1_hdr['rdn_img']['start datetime']
    nc_fid.flight_line = dt_start.strftime("prm%Y%m%dt%H%M%S")
    nc_fid.time_coverage_start = dt_start.strftime("%Y-%m-%dT%H:%M:%SZ")
    dt_end = dt_start + timedelta(microseconds=np.nanmax(dsL1_data['utc time']))
    nc_fid.time_coverage_end = dt_end.strftime("%Y-%m-%dT%H:%M:%SZ")
    nc_fid.northernmost_latitude = float(np.nanmax(dsL1_data['latitude (wgs-84)']))
    nc_fid.southernmost_latitude = float(np.nanmin(dsL1_data['latitude (wgs-84)']))
    nc_fid.easternmost_longitude = float(np.nanmax(dsL1_data['longitude (wgs-84)']))
    nc_fid.westernmost_longitude = float(np.nanmin(dsL1_data['longitude (wgs-84)']))
    nc_fid.geospatial_lat_max = float(np.nanmax(dsL1_data['latitude (wgs-84)']))
    nc_fid.geospatial_lat_min = float(np.nanmin(dsL1_data['latitude (wgs-84)']))
    nc_fid.geospatial_lon_max = float(np.nanmax(dsL1_data['longitude (wgs-84)']))
    nc_fid.geospatial_lon_min = float(np.nanmin(dsL1_data['longitude (wgs-84)']))
    nc_fid.day_night_flag = "Day"
    nc_fid.earth_sun_distance_correction = float(np.nanmean(dsL1_data['earth-sun distance (au)']))

    #fill sbp-group vars
    mainLogger.info('Filling and writing sensor_band_parameters')
    nc_fid.groups['sensor_band_parameters'].variables['wavelength'][:] = bands_ocip
    nc_fid.groups['sensor_band_parameters'].variables['fwhm'][:] = widths_ocip

    #fill sla-group vars
    mainLogger.info('Filling and writing scan_line_attributes')
    dsL1_data['scan_start_time'][np.isnan(dsL1_data['scan_start_time'])] = nc_fid.groups['scan_line_attributes'].variables['scan_start_time']._FillValue
    nc_fid.groups['scan_line_attributes'].variables['scan_start_time'][:] = dsL1_data['scan_start_time']

    dsL1_data['scan_end_time'][np.isnan(dsL1_data['scan_end_time'])] = nc_fid.groups['scan_line_attributes'].variables['scan_end_time']._FillValue
    nc_fid.groups['scan_line_attributes'].variables['scan_end_time'][:] = dsL1_data['scan_end_time']

    #fill nd-group vars
    mainLogger.info('Filling and writing navigation_data')
    dsL1_data['longitude (wgs-84)'][np.isnan(dsL1_data['longitude (wgs-84)'])] = nc_fid.groups['navigation_data'].variables['longitude']._FillValue
    nc_fid.groups['navigation_data'].variables['longitude'][:,:] = dsL1_data['longitude (wgs-84)']

    dsL1_data['latitude (wgs-84)'][np.isnan(dsL1_data['latitude (wgs-84)'])] = nc_fid.groups['navigation_data'].variables['latitude']._FillValue
    nc_fid.groups['navigation_data'].variables['latitude'][:,:] = dsL1_data['latitude (wgs-84)']

    dsL1_data['elevation (m)'][np.isnan(dsL1_data['elevation (m)'])] = nc_fid.groups['navigation_data'].variables['height']._FillValue
    nc_fid.groups['navigation_data'].variables['height'][:,:] = dsL1_data['elevation (m)']

    dsL1_data['path length (m)'][np.isnan(dsL1_data['path length (m)'])] = nc_fid.groups['navigation_data'].variables['range']._FillValue
    nc_fid.groups['navigation_data'].variables['range'][:,:] = dsL1_data['path length (m)']

    dsL1_data['to-sensor zenith (0 to 90 degrees from zenith)'][np.isnan(dsL1_data['to-sensor zenith (0 to 90 degrees from zenith)'])] = nc_fid.groups['navigation_data'].variables['sensor_zenith']._FillValue
    nc_fid.groups['navigation_data'].variables['sensor_zenith'][:,:] = dsL1_data['to-sensor zenith (0 to 90 degrees from zenith)']

    dsL1_data['to-sensor azimuth (0 to 360 degrees cw from n)'][np.isnan(dsL1_data['to-sensor azimuth (0 to 360 degrees cw from n)'])] = nc_fid.groups['navigation_data'].variables['sensor_azimuth']._FillValue
    nc_fid.groups['navigation_data'].variables['sensor_azimuth'][:,:] = dsL1_data['to-sensor azimuth (0 to 360 degrees cw from n)']

    dsL1_data['to-sun zenith (0 to 90 degrees from zenith)'][np.isnan(dsL1_data['to-sun zenith (0 to 90 degrees from zenith)'])] = nc_fid.groups['navigation_data'].variables['solar_zenith']._FillValue
    nc_fid.groups['navigation_data'].variables['solar_zenith'][:,:] = dsL1_data['to-sun zenith (0 to 90 degrees from zenith)']

    dsL1_data['to-sun azimuth (0 to 360 degrees cw from n)'][np.isnan(dsL1_data['to-sun azimuth (0 to 360 degrees cw from n)'])] = nc_fid.groups['navigation_data'].variables['solar_azimuth']._FillValue
    nc_fid.groups['navigation_data'].variables['solar_azimuth'][:,:] = dsL1_data['to-sun azimuth (0 to 360 degrees cw from n)']

    #fill od-group vars
    mainLogger.info('Filling and writing observation_data')
    lt_ocip[np.isnan(lt_ocip)] = nc_fid.groups['observation_data'].variables['Lt']._FillValue
    nc_fid.groups['observation_data'].variables['Lt'][:,:,:] = lt_ocip

    #close nc file
    nc_fid.close()
    mainLogger.info('Successfully created and filled: {:}'.format(ofile))

    return

if __name__ == "__main__":
    PRISMmain(sys.argv[1:])
