#!/usr/bin/env python3

"""
@author: aakmal

Based on IDL code by wdrobbins:
cams_monthly_average.pro 
monthly_mean_trace_gas_concentrations_on_merra_grid.pro
ncdf_read_all.pro  
run.pro

Information from original IDL version:

file info.txt

information on the CAMS gas profile data

Original source:
  https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/
    cams-global-greenhouse-gas-inversion?tab=form

for CH4, select Methane,concentration, surface air sample and satellite,
    daily mean, latest version, 2019
for CO2, select Carbon dioxide, concentration, surface air sample,
    instantaneous, latest version, 2020
for N2O, select Nitrous oxide, concentration, surface air sample,
    instantaneous, latest version, 2020

You will get 12 monthly files of the gas data on a hybrid sigma vertical 
coordinate system with 39 levels, n2o 96 in lat, 97 in lon, 8 time / day
co2 same, ch4 34 levels,90 in lat, 120 in lon, 1 time / day

When distributing the data, include the following
'Generated using Copernicus Atmosphere Monitoring Service information [Year]’ 
with [Year] 2019 for CH$ and N2O and 2020 for CO2

Local storage for data:
The ECMWF files and final average file are under area:
~/data/ancillary/ecmwf/cams/base_ecmwf_data
in the following sub-directories:  CH4_2019  CO2_2020  N2O_2019

Program:

The program to make the average files is in
~/src/idl_lib/cams_avg/lcl_code
To run:
- edit the file monthly_mean_trace_gas_concentrations_on_merra_grid.pro
and insert the gas to process: 'co2', 'n2o', or 'ch4'
- invoke IDL
- compile the routines:
   .compile monthly_mean_trace_gas_concentrations_on_merra_grid
   .compile ncdf_read_all
   .compile cams_monthly_average
   .compile run
- Make the average:
  run

final file will be in the CH4_2019  CO2_2020  N2O_2019 directories 
respectively, names:
  cams73_latest_ch4_conc_surface_satellite_dm_2019.nc
  cams73_latest_co2_conc_surface_inst_2020.nc
  cams73_latest_n2o_conc_surface_inst_2019.nc

They are averaged for each month and on the 42 level pressure coordinates 
used by MERRA2

"""

import netCDF4 as nc4
import numpy as np
import numpy.ma as ma
import datetime as dt
from datetime import timedelta, timezone
from dateutil.relativedelta import relativedelta
import sys

__version__ = '1.0 (2025-04-24)'

#locate a value (val) in a monotonically decreasing array (vec)
def value_locate(vec, val):
    indices = np.where(val >= vec)
   #print(indices)
    if indices[0].size == 0:
            return vec.size-1
    else:
            return indices[0][0]-1

#read a monthly CAMS netCDF file and return a dict holding all info
def ncdf_read_all(cdffile):
    print("reading ...")
    nc_dataset = nc4.Dataset(cdffile, mode='r', format='NETCDF4')
    sCDFinfo = {'ndims': len(nc_dataset.dimensions)}
    sCDFinfo['nvars'] = len(nc_dataset.variables)
    sCDFinfo['ngatts'] = len(nc_dataset.ncattrs())
    sCDFinfo['recdim'] = -1

    loc = -1
    for i in nc_dataset.dimensions.keys():
        loc += 1
        if(nc_dataset.dimensions[i].isunlimited):
                sCDFinfo['recdim'] = loc

    global_atts = nc_dataset.ncattrs()
    sGatts = {}
    cnt = 0
    for i in global_atts:
        sGatts[i] = getattr(nc_dataset, i)
        cnt +=1

    sVars = {}
    variables = nc_dataset.variables.keys()
    for i in variables:
        sVars[i] = {'data': nc_dataset[i][:], 'atts': {}}
        var_atts = nc_dataset[i].ncattrs()
        for j in var_atts:
            sVars[i]['atts'][j] = getattr(nc_dataset[i], j)

    s = {'info': sCDFinfo, 'gatts': sGatts, 'vars': sVars}

    return s

#create an month average profile on MERRA grid
def cams_monthly_average(s, merra_pressure, gas):

    print("averaging ...")

    nlat = (s['vars']['latitude']['data']).size
    nlon = (s['vars']['longitude']['data']).size
    ntime = (s['vars']['time']['data']).size

#Define the missing value used for grid points that are below the surface 
    MISSING_VALUE=-999.0

    npress_merra = merra_pressure.size
    merra_pressure = 100.0*merra_pressure #Convert to Pa, since gas concentration files use Pa as their pressure co-ordinate

#Get ap and bp coefficients for the hybrid sigma co-ordinate of the concentraion array
#pressure = ap + bp*surface_pressure

    if(gas == 'co2' or gas == 'n2o'):
        nedge=s['vars']['ap']['data'].size
        nmid=s['vars']['level']['data'].size
        edge=np.arange(float(nedge)) 
        mid=s['vars']['level']['data']+0.5
        ap=np.interp(mid,edge,s['vars']['ap']['data'])
        bp=np.interp(mid,edge,s['vars']['bp']['data'])
        srfc_prs=s['vars']['Psurf']['data']
    elif(gas == 'ch4'):
        nmid=s['vars']['hyam']['data'].size
        ap   = s['vars']['hyam']['data']
        bp   = s['vars']['hybm']['data']
        srfc_prs = s['vars']['ps']['data']

#Calculate average in time of surface pressure so that concentrations can be averaged on a fixed pressure grid
    avg_srfc_prs = np.average(srfc_prs, axis=0)
    
   #print("srfc_prs: ", srfc_prs.shape)
   #print("avg_srf_prs: ", avg_srfc_prs.shape)
    
   #print("ap: ", ap.shape)
   #print("bp: ", bp.shape)

   #print("nlat, nlon, ntime, npress_merra: ",nlat, nlon, ntime, npress_merra)
   #print(s['vars']['ps']['atts']['add_offset'], s['vars']['ps']['atts']['scale_factor'])

#Get gas concentration arrays -  pressure, time, pressure, lat, lon; scale/offset already applied by netCDF4 library
    xgas = s['vars'][gas.upper()]['data']
#Create an array to store average gas concentration on the original horizontal and vertical grid
    xgas_average           = np.empty((nmid, nlat, nlon), dtype=np.float32)
    xgas_average_merra_prs = np.empty((npress_merra, nlat, nlon), dtype=np.float32) 

   #print("xgas, xgas_average, xgas_average_merra_prs: ", xgas.shape, xgas_average.shape, xgas_average_merra_prs.shape)
   #print("xgas: ", type(xgas))
   #print("xgas_average: ",type(xgas_average))
   #print("xgas_average_merra_prs: ",type(xgas_average_merra_prs))

    for ilat in range(nlat):
        for ilon in range(nlon):
            #Average sigma pressure profile for this location
            average_pressure=ap + bp*avg_srfc_prs[ilat,ilon]
            #Array to average instantaneous profiles into
            average_profile=np.zeros((nmid))
            #Time average
            for itime in range(ntime):
                pressure=ap + bp*srfc_prs[itime, ilat, ilon]
                inst_profile = np.interp(average_pressure, pressure[np.argsort(pressure)], xgas[itime,np.argsort(pressure), ilat, ilon])
                average_profile=average_profile+inst_profile

            average_profile=average_profile/float(ntime)
            xgas_average[:,ilat,ilon]=average_profile

            #Check where surface pressure lies on Merra pressure grid
            iz=value_locate(merra_pressure,average_pressure[0])	
            #Create an array for the merra profile filled with missing values
            merra_profile=np.full(npress_merra,MISSING_VALUE)
            if iz < 0:
                iz=0 #Can interpolate over full Merra pressure range
            #interpolate gas concentrations to merra grid with missing values for pres. levels below surface pressure
            merra_profile[iz:npress_merra-1] = np.interp(merra_pressure[iz:npress_merra-1], average_pressure[np.argsort(average_pressure)], average_profile[np.argsort(average_pressure)])
            #Save profile into gridded array
            xgas_average_merra_prs[:,ilat, ilon] = merra_profile

#   merra_pressure = merra_pressure/100.0 #Convert back to hPa, which is standar MERRA usage
    return xgas_average_merra_prs


def main(args):
    print("tracegas_avg", __version__)
    
    if len(args) < 3:
        print("Usage: tracegas_avg gas year [ base data directory (default: ./base_ecmwf_data) ]")
        exit()
    elif (args[1] != 'ch4' and args[1] != 'co2' and args[1] != 'n2o'):
        print("gas must be: ch4, co2 or ch4")
        exit()

    gas     = str(args[1])
   #gas     = 'ch4'
   #gas     = 'co2'
   #gas     = 'n2o'

    if len(args) == 4:
        src_dir = str(args[-1])
    else:
        src_dir = './base_ecmwf_data'

#URL that the data are downloaded from
   #creator_url="https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-global-greenhouse-gas-inversion?tab=overview"
    creator_url = "https://oceandata.sci.gsfc.nasa.gov" ;

    MISSING_VALUE=-999.0

    merra_pressure=np.array([1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 725, 700, 650, 600, 550, 500, 450, 
    400, 350, 300, 250, 200, 150, 100, 70, 50, 40, 30, 20, 10, 7, 5, 4, 3, 2, 1, 0.7, 0.5, 0.4, 0.3, 0.1])

    npress = merra_pressure.size
    
    year = args[2]
   #year = 2019
   #year = 2020
   #dat_dir = src_dir + '/CH4_2019/'
   #dat_dir = src_dir + '/CO2_2020/'
   #dat_dir = src_dir + '/N2O_2020/'
    dat_dir = src_dir + '/' + gas.upper() + '_' + str(year) + '/'

   #stub    = 'cams73_latest_'+gas+'_conc_surface_satellite_dm_2019'
   #stub    = 'cams73_latest_'+gas+'_conc_surface_inst_2020'
    if gas == 'ch4' :
        stub = 'cams73_latest_'+gas+'_conc_surface_satellite_dm_' + str(year)
    else:
        stub = 'cams73_latest_'+gas+'_conc_surface_inst_' + str(year)

    creation_time=dt.datetime.now(dt.timezone.utc).isoformat()
    print("creation_time: ", creation_time)

    outfile =  dat_dir + stub + ".nc"
#   outfile = 'example.nc'
    print("Gas: ", gas.upper())
    print("Writing to: ", outfile)

    nmon = 12 

    s = {}
    start_time = {}
    end_time = {}

    for imon in range(1, nmon+1):
        fname = dat_dir + stub + "{:02d}".format(imon) + ".nc"
        s[imon] = ncdf_read_all(fname)    #dict holding full data structures for each month

        epoch_start_ts = s[imon]['vars']['time']['atts']['units']
        time_data = s[imon]['vars']['time']['data']

#Get start and end times of data being averaged
        if( gas == 'ch4'):
            a=dt.datetime.strptime(epoch_start_ts,'hours since %Y-%m-%d %H:%M')
            end_time[imon] = (a + relativedelta(months=imon)).replace(tzinfo=timezone.utc).isoformat()
        else:
            a=dt.datetime.strptime(epoch_start_ts,'hours since %Y-%m-%d %H:%M:%S')
            end_time[imon] = (a + relativedelta(months=1)).replace(tzinfo=timezone.utc).isoformat()

        start_time[imon] = (a + timedelta(hours=time_data[0])).replace(tzinfo=timezone.utc).isoformat()

    nlat = (s[1]['vars']['latitude']['data']).size
    nlon = (s[1]['vars']['longitude']['data']).size

    xgas_atts=s[1]['vars'][gas.upper()]['atts']



    avg_conc = {} #dict holding variable names for average concentrations for each month

# Create a new NetCDF file (use 'w' for write mode)
    dataset = nc4.Dataset(outfile, 'w', format='NETCDF4')

# Define dimensions

    month_id     = dataset.createDimension('time', nmon)
    latitude_id  = dataset.createDimension('latitude', nlat)
    longitude_id = dataset.createDimension('longitude', nlon)
    pressure_id  = dataset.createDimension('pressure', npress)

# Add global attributes
    dataset.title = "Monthly mean averages of CAMS " + gas.upper() + " trace gas concentrations for PACE mission"
   #dataset.creation_time = creation_time
    dataset.date_created = creation_time
    dataset._FillValue = MISSING_VALUE
    dataset.creator_url = creator_url
    dataset.standard_name_vocabulary = "CF Standard Name Table v36" 
    dataset.creator_name = "NASA/GSFC/OBPG" 
    dataset.creator_email = "data@oceancolor.gsfc.nasa.gov" 

    dataset.time_coverage_start = start_time[1]
    dataset.time_coverage_end   = end_time[nmon]
    dataset.geospatial_lat_min=s[1]['vars']['latitude']['data'][0] 
    dataset.geospatial_lat_max=s[1]['vars']['latitude']['data'][-1] 
    dataset.geospatial_lon_min=s[1]['vars']['longitude']['data'][0] 
    dataset.geospatial_lon_max=s[1]['vars']['longitude']['data'][-1]
    dataset.geospatial_vertical_min=str(merra_pressure[-1])+" hPa" 
    dataset.geospatial_vertical_max=str(merra_pressure[0])+" hPa" 
    dataset.geospatial_vertical_positive="down"
    dataset.processing_level="L4"
    dataset.publisher_name = "NASA/GSFC/OBPG" 
    dataset.publisher_url = "https://oceandata.sci.gsfc.nasa.gov" 
    dataset.publisher_email = "data@oceancolor.gsfc.nasa.gov" 
    dataset.time_coverage_duration="M"
    dataset.time_coverage_resolution="M"

    for global_attr in s[1]['gatts'].keys():
        if global_attr != 'title' and global_attr != 'source' and global_attr != 'release' and global_attr != 'author' :
            setattr(dataset, global_attr, s[1]['gatts'][global_attr])

    dataset.history = dataset.history + "; " + dataset.date_created + ": " + "tracegas_avg " + gas + " " +  year

    if gas == 'ch4' :
        dataset.comment = s[1]['gatts']['title'] + "; source: " + s[1]['gatts']['source'] + \
        "; release: " + s[1]['gatts']['release'] + "; author: " +s[1]['gatts']['author']
    else:
        dataset.comment = s[1]['gatts']['title'] + "; source: " + s[1]['gatts']['source'] 
        
    dataset.Conventions = "CF-1.8, ACDD-1.3" 
    dataset.institution = "NASA Goddard Space Flight Center, Ocean Ecology Laboratory, Ocean Biology Processing Group" 

# Create variables
    conc = gas.upper()
    latitude  = dataset.createVariable('latitude', np.float32, ('latitude',))
    longitude = dataset.createVariable('longitude', np.float32, ('longitude',))
    pressure  = dataset.createVariable('pressure', np.float32, ('pressure',))
    time      = dataset.createVariable('time', np.uint8, ('time',))
    conc_avg  = dataset.createVariable(conc, np.float32, ('time', 'pressure', 'latitude', 'longitude'),
            fill_value= MISSING_VALUE)

#   for imon in range(1,nmon+1):
#       gas_conc=gas.upper()+"_"+ "{:02d}".format(imon)
#       avg_conc[imon]  = dataset.createVariable(gas_conc, np.float32,('pressure', 'latitude', 'longitude'), 
#       fill_value= MISSING_VALUE)

# Add attributes to variables 
    setattr(latitude, 'long_name', 'latitude')
    for attr in s[1]['vars']['latitude']['atts'].keys():
        setattr(latitude, attr, s[1]['vars']['latitude']['atts'][attr])
    setattr(latitude, 'coverage_content_type', 'coordinate')

    setattr(longitude, 'long_name', 'longitude')
    for attr in s[1]['vars']['longitude']['atts'].keys():
        setattr(longitude, attr, s[1]['vars']['longitude']['atts'][attr])
    setattr(longitude, 'coverage_content_type', 'coordinate')

    setattr(pressure, 'long_name', 'MERRA2 Pressure Levels')
    setattr(pressure, 'units', 'hPa')
    setattr(pressure, 'valid_min', np.float32(0.1))
    setattr(pressure, 'valid_max', 1000.0)
    setattr(pressure, 'coverage_content_type', 'coordinate')

    setattr(time, 'long_name', 'Month of year')
    setattr(time, 'units', 'month')
    setattr(time, 'valid_min', np.uint8(1))
    setattr(time, 'valid_max', np.uint8(12))
    setattr(time, 'coverage_content_type', 'coordinate')

    for attr in xgas_atts.keys():
        if( (attr != "add_offset") and (attr != "scale_factor") and (attr != 'cell_measures') and (attr != 'cell_methods')):
            setattr(conc_avg, attr, xgas_atts[attr])
    setattr(conc_avg, 'coverage_content_type', 'modelResult')


# Write data to variables
    latitude[:]  = s[1]['vars']['latitude']['data']
    longitude[:] = s[1]['vars']['longitude']['data']
    pressure[:] = merra_pressure
    time[:] = np.arange(1, 13, dtype=np.uint8)

    for imon in range(1,nmon+1):
        print("Month: ", imon)
        print("start_time: ", start_time[imon])
        print("end_time:   ", end_time[imon])
        conc_avg[imon-1,:,:] = cams_monthly_average(s[imon], merra_pressure, gas)

# Close the dataset to save changes
    dataset.close()

if __name__ == "__main__":
    main(sys.argv)
