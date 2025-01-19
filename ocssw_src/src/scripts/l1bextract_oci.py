#! /usr/bin/env python3

# Extractor for L1B_OCI files

import argparse
import netCDF4
import pathlib
from os.path import basename
import numpy as np
import sys
import time
import datetime
from shutil import copyfile as cp
from seadasutils.pixlin_utils import pixlin
from seadasutils.setupenv import env
#from seadasutils.MetaUtils import readMetadata
from seadasutils.netcdf_utils import ncsubset_vars
#import seadasutils.ProcUtils as ProcUtils

versionStr = "1.2 (2025-01-13)"
global pixDim
global scanDim

class extract:

    def __init__(self, ifile, ofile=None,
                 north=None, south=None, west=None, east=None,
                 spixl=None, epixl=None, sline=None, eline=None,
                 verbose=False):
        # inputs
        self.ifile = pathlib.Path(ifile)
        self.ofile = pathlib.Path(ofile)
        self.north = north
        self.south = south
        self.west = west
        self.east = east
        self.spixl = spixl
        self.epixl = epixl
        self.sline = sline
        self.eline = eline
        self.verbose = verbose

        # defaults
        self.runtime = None
        self.attrs = None

        # unused, but needed by setupenv.py
        self.dirs = {}
        self.ancdir = None
        self.curdir = False
        self.sensor = None
        env(self)  # run setupenv

    def runextract(self, subset):
        srcfile = self.ifile
        if srcfile.exists():
            dstfile = self.ofile
            if self.verbose:
                print('Extracting', srcfile)

            ncsubset_vars(srcfile, dstfile, subset, timestamp=self.runtime)

            # Read infile as netCDF dataset
            infile = netCDF4.Dataset(srcfile, 'r')
            # Read and write outfile as netCDF4 dataset
            outfile = netCDF4.Dataset(dstfile, 'r+')
        
            #_____________________________________________________________________________________________________________________________________
            #                                                                                 |
            # Modify time_coverage_start/end of output file:                                  |
            #_________________________________________________________________________________|
            
            # Read infile as netCDF dataset
            infile = netCDF4.Dataset(srcfile, 'r')
            # Read and write outfile as netCDF4 dataset
            outfile = netCDF4.Dataset(dstfile, 'r+')
         
            # Number of seconds at infile time_coverage_start/end
            infile_start_sec: float = infile.groups['scan_line_attributes'].variables['time'][0]
            infile_end_sec = infile.groups['scan_line_attributes'].variables['time'][infile.dimensions[scanDim].size - 1]
            # Number of seconds at outfile time_coverage_start/end
            outfile_start_sec = outfile.groups['scan_line_attributes'].variables['time'][0]
            outfile_end_sec = outfile.groups['scan_line_attributes'].variables['time'][outfile.dimensions[scanDim].size - 1]

            # Take infile time_coverage_start/end
            infile_start_time = infile.time_coverage_start
            infile_end_time = infile.time_coverage_end

            # Extract year, month, day, hours, minutes, seconds from infile time coverage start/end
            start_form = datetime.datetime.strptime(infile_start_time[0:20] + '000', '%Y-%m-%dT%H:%M:%S.%f')
            end_form = datetime.datetime.strptime(infile_end_time[0:20] + '000', '%Y-%m-%dT%H:%M:%S.%f')
            # Extract year,...,seconds from epoch
            epoch = datetime.datetime.strptime('1970 01 01 00 00 00.000','%Y %m %d %H %M %S.%f')
            # Determine difference in time from infile time_coverage_start to epoch
            diff_start = start_form - epoch
            # Determine difference in time from infile time_coverage_end to epoch
            diff_end = end_form - epoch

            # Calculate the number of seconds contained in the time difference in previous step
            diff_sec_start = diff_start.total_seconds()
            # Calculate the number of seconds contained in the time difference in previous step
            diff_sec_end = diff_end.total_seconds()

            # Seconds between infile/outfile time_coverage_start/end
            diff_infile_outfile_start = outfile_start_sec - infile_start_sec
            diff_infile_outfile_end = outfile_end_sec - infile_end_sec

            # Add the input/output file time_coverage_start/end difference to the infile time_coverage_start/end # of seconds
            outfile_tot_start_sec = diff_sec_start + diff_infile_outfile_start
            outfile_tot_end_sec = diff_sec_end + diff_infile_outfile_end

            # Create time structure for the outfile time_coverage_start/end
            outfile_start_time_since = time.gmtime(outfile_tot_start_sec)
            outfile_end_time_since = time.gmtime(outfile_tot_end_sec)

            # Extract year, month, day, hours, minutes, seconds from outfile start/end time structs
            ostart_y = outfile_start_time_since.tm_year
            ostart_mon = "{0:0=2d}".format(outfile_start_time_since.tm_mon)
            ostart_d = "{0:0=2d}".format(outfile_start_time_since.tm_mday)
            ostart_h = "{0:0=2d}".format(outfile_start_time_since.tm_hour)
            ostart_min = "{0:0=2d}".format(outfile_start_time_since.tm_min)
            ostart_s = "{0:0=2d}".format(outfile_start_time_since.tm_sec)
            
            oend_y = outfile_end_time_since.tm_year
            oend_mon = "{0:0=2d}".format(outfile_end_time_since.tm_mon)
            oend_d = "{0:0=2d}".format(outfile_end_time_since.tm_mday)
            oend_h = "{0:0=2d}".format(outfile_end_time_since.tm_hour)
            oend_min = "{0:0=2d}".format(outfile_end_time_since.tm_min)
            oend_s = "{0:0=2d}".format(outfile_end_time_since.tm_sec)

            # Change outfile time_coverage_start/end
            outfile.time_coverage_start = str(ostart_y) + '-' + str(ostart_mon) + '-' + str(ostart_d) + 'T' + str(ostart_h) + ':' + str(ostart_min) + ':' + str(ostart_s)
            outfile.time_coverage_end = str(oend_y) + '-' + str(oend_mon) + '-' + str(oend_d) + 'T' + str(oend_h) + ':' + str(oend_min) + ':' + str(oend_s) 
            
            #_____________________________________________________________________________________________________________________________________
            #                                                                                 |
            # Gring Calculation:                                                              |
            #_________________________________________________________________________________|
            
            outfile.set_auto_mask(False)

            scans = outfile.dimensions[scanDim].size - 1
            pixels = outfile.dimensions[pixDim].size - 1
            
            latitude = outfile.groups['geolocation_data'].variables['latitude']
            longitude = outfile.groups['geolocation_data'].variables['longitude']
            
            lon_min = longitude[0, 0]
            lon_max = longitude[scans, pixels]
            lat_min = latitude[0, 0]
            lat_max = latitude[scans, pixels]
            
            # Degrees to separate latitude points, here it is 20
            lat_add = float((lat_max - lat_min) / 20)
        
            # Place one point in the middle of longitude
            lat_l = []
            lat_r = []
            lon_l = []
            lon_r = []

            # add latitude right values
            if lat_add > 0 and int(scans/lat_add) != 0:
                for i in range(0, scans - 1, int(scans/lat_add)):
                    lat_r.append(i)
                    # add longitude left/right values
                    lon_l.append(0)
                    lon_r.append(pixels)
            else:
                 for i in range(0, scans - 1, 1):
                    lat_r.append(i)
                    # add longitude left/right values
                    lon_l.append(0)
                    lon_r.append(pixels)

            # add latitude left values
            lat_l = list(reversed(lat_r))
            
            # add longitude up values
            lon_u = [pixels, (pixels/2), 0]
            # add latitude up values
            lat_u = [scans, scans, scans]
            # add longitude down values
            lon_d = list(reversed(lon_u)) 
            # add longitude up values
            lat_d = [0, 0, 0]
            
            lat_values_u = []
            lon_values_u = []
            lat_values_d = []
            lon_values_d = []
            lat_values_l = []
            lon_values_l = []
            lat_values_r = []
            lon_values_r = []
        
            num = 0
            for i in range(len(lat_u)):
                lat_values_u.append(float(lat_u[i]))
                lon_values_u.append(float(lon_u[i]))
                num += 1
            for i in range(len(lat_l)):
                lat_values_l.append(float(lat_l[i]))
                lon_values_l.append(float(lon_l[i]))
                num += 1
            for i in range(len(lat_d)):
                lat_values_d.append(float(lat_d[i]))
                lon_values_d.append(float(lon_d[i]))
                num += 1
            for i in range(len(lat_r)):
                lat_values_r.append(float(lat_r[i]))
                lon_values_r.append(float(lon_r[i]))
                num += 1
        
            p_seq = []
            
            for i in range(num):
                p_seq.append(np.dtype('int32').type(i + 1))
            
            args_lat = (lat_values_u, lat_values_l, lat_values_d, lat_values_r)
            args_lon = (lon_values_u, lon_values_l, lon_values_d, lon_values_r)
            lat_values = np.concatenate(args_lat)
            lon_values = np.concatenate(args_lon)
            
            g_lat = []
            g_lon = []
            
            for i in range(0,len(lat_values)):
                g_lat.append(latitude[int(lat_values[i])][int(lon_values[i])])
                g_lon.append(longitude[int(lat_values[i])][int(lon_values[i])])
                
            outfile.setncattr('GRingPointLatitude', g_lat)
            outfile.setncattr('GRingPointLongitude', g_lon)
            outfile.setncattr('GRingPointSequenceNo', p_seq)   
            
            #_____________________________________________________________________________________________________________________________________
            #                                                                                 |
            # Geolocation lat/lon min/max update:                                              |
            #_________________________________________________________________________________|
            def minmax(arr):
                return arr.min(), arr.max()

            lat_min, lat_max = latitude[0, 0], latitude[scans, pixels]
            lon_min, lon_max = longitude[0, 0], longitude[scans, pixels]

            outfile.setncattr('geospatial_lat_min', lat_min)
            outfile.setncattr('geospatial_lat_max', lat_max)
            outfile.setncattr('geospatial_lon_min', lon_min)
            outfile.setncattr('geospatial_lon_max', lon_max)
            
        return 0

    def getpixlin(self):
        if self.verbose:
            print("north={} south={} west={} east={}".
                  format(self.north, self.south, self.west, self.east))

        # run lonlat2pixline
        pl = pixlin(geofile=self.ifile,
                    north=self.north, south=self.south,
                    west=self.west, east=self.east,
                    verbose=self.verbose)
        pl.lonlat2pixline(zero=False)  # using 1-based indices
        self.spixl, self.epixl, self.sline, self.eline = \
        (pl.spixl, pl.epixl, pl.sline, pl.eline)
        return pl.status

    def run(self):
        global scanDim
        # convert to zero-based index
        self.spixl, self.epixl, self.sline, self.eline = \
        (v-1 for v in (self.spixl, self.epixl, self.sline, self.eline))

        dims = list(infile.dimensions.keys())

        # extract file
        subset = {'pixels':[self.spixl, self.epixl],
            'scans': [self.sline, self.eline]}
        if scanDim == 'number_of_scans':
            # prior to V3 dimension change
            subset = {'SWIR_pixels':[self.spixl, self.epixl],
                    'ccd_pixels': [self.spixl, self.epixl],
                    'number_of_scans': [self.sline, self.eline]}
        elif 'SWIR_pixels' in dims and 'scans' in dims:
            # V3 non-baseline configuration
            subset = {'SWIR_pixels':[self.spixl, self.epixl],
                    'pixels': [self.spixl, self.epixl],
                    'scans': [self.sline, self.eline]}
        self.runextract(subset)

        #return status
        return 0

def chk_pixl(args, infile):
    global pixDim
    global scanDim
    if args.epixl == -1:
        args.epixl = infile.dimensions[pixDim].size
    if args.eline == -1:
        args.eline = infile.dimensions[scanDim].size
    return args.spixl, args.epixl, args.sline, args.eline
   
if __name__ == "__main__":
    print("l1bextract_oci", versionStr)
    status = 0
    # parse command line
    parser = argparse.ArgumentParser(
        description='Extract specified area from OCI Level 1B files.',
        epilog='Specify either geographic limits or pixel/line ranges, not both.')
    parser.add_argument('-v', '--verbose', help='print status messages',
                        action='store_true')
    parser.add_argument('ifile',
                        help='Level 1B input file')
    parser.add_argument('ofile', nargs='?',
                        help='output file')

    group1 = parser.add_argument_group('geographic limits')
    group1.add_argument('-n', '--north', type=float, help='northernmost latitude')
    group1.add_argument('-s', '--south', type=float, help='southernmost latitude')
    group1.add_argument('-w', '--west', type=float, help='westernmost longitude')
    group1.add_argument('-e', '--east', type=float, help='easternmost longitude')

    group2 = parser.add_argument_group('pixel/line ranges (1-based)')
    group2.add_argument('--spixl', type=int, help='start pixel', default = 1)
    group2.add_argument('--epixl', type=int, help='end pixel', default = -1)

    group2.add_argument('--sline', type=int, help='start line', default = 1)
    group2.add_argument('--eline', type=int, help='end line', default = -1)
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1) 
    args = parser.parse_args()

    pixel_bounds_specified = not (args.eline == -1 and args.sline == 1 and args.epixl == -1 and args.spixl == 1)
    
    # Read infile as netCDF dataset
    infile = netCDF4.Dataset(args.ifile, 'r')
    dims = list(infile.dimensions.keys())
    pixDim = 'pixels'
    scanDim = 'scans'
    if 'number_of_scans' in dims:
        pixDim = 'SWIR_pixels'
        scanDim = 'number_of_scans'

    args.spixl, args.epixl, args.sline, args.in_eline = chk_pixl(args, infile)
             
    # initialize
    this = extract(ifile=args.ifile,
                   ofile=args.ofile,
                   north=args.north,
                   south=args.south,
                   west=args.west,
                   east=args.east,
                   spixl=args.spixl,
                   epixl=args.epixl,
                   sline=args.sline,
                   eline=args.eline,
                   verbose=args.verbose)

    # input value checks
    goodlatlons = None not in (this.north, this.south, this.west, this.east)
    if (goodlatlons and pixel_bounds_specified):
        print("ERROR: Specify either geographic limits or pixel/line ranges, not both.")
        sys.exit(1)
    elif (goodlatlons and not pixel_bounds_specified):
        status = this.getpixlin()
        if status not in (0, 110):
            print("No extract; lonlat2pixline status =", status)
            exit(status)
        status = this.run()
    elif (pixel_bounds_specified and not goodlatlons):
        status = this.run()
    else:
        print("No extract; subset not specified")
        status = 1

    exit(status)
