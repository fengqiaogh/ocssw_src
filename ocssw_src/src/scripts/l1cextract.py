#! /usr/bin/env python3

# Extractor for L1C files

import argparse
import netCDF4
import pathlib
from os.path import basename
import numpy as np
import sys
import time
import calendar
import datetime
import re
import math
from shutil import copyfile as cp
import xml.etree.ElementTree as ET
from seadasutils.pixlin_utils import pixlin
from seadasutils.setupenv import env
#from seadasutils.MetaUtils import readMetadata
from seadasutils.netcdf_utils import nccopy, nccopy_grp, ncsubset_vars
#import seadasutils.ProcUtils as ProcUtils

versionStr = "1.3 (2024-09-23)"

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

            ncsubset_vars(srcfile, dstfile, subset, timestamp=self.runtime, verbose=self.verbose)

            # Read infile as netCDF dataset
            infile = netCDF4.Dataset(srcfile, 'r')
            # Read and write outfile as netCDF4 dataset
            outfile = netCDF4.Dataset(dstfile, 'r+')
            
            # Update nadir_bin 
            try:
                if infile.nadir_bin:
                    nadir_bin = np.dtype('int32').type(infile.nadir_bin)
                    if (nadir_bin > self.epixl):
                        outfile.nadir_bin = np.dtype('int32').type(-1)
                    else:
                        outfile.nadir_bin = np.dtype('int32').type(nadir_bin - (self.spixl + 1))
            except:
                pass

            #_____________________________________________________________________________________________________________________________________
            #                                                                                 |
            # Add extract_pixel/line_start/stop:                                              |
            #_________________________________________________________________________________|
            
            if 'extract_pixel_start' in infile.ncattrs():
                outfile.extract_pixel_start = np.dtype('int32').type(infile.extract_pixel_start + self.spixl + 1)
                outfile.extract_pixel_stop  =  np.dtype('int32').type(infile.extract_pixel_stop + self.epixl + 1)
                outfile.extract_line_start  =  np.dtype('int32').type(infile.extract_line_start + self.sline + 1)
                outfile.extract_line_stop   =  np.dtype('int32').type(infile.extract_line_stop + self.eline + 1)
            else:
                outfile.extract_pixel_start = np.dtype('int32').type(self.spixl + 1)
                outfile.extract_pixel_stop  =  np.dtype('int32').type(self.epixl + 1)
                
                outfile.extract_line_start  =  np.dtype('int32').type(self.sline + 1)
                outfile.extract_line_stop   =  np.dtype('int32').type(self.eline + 1)
            
            #_____________________________________________________________________________________________________________________________________
            #                                                                                 |
            # Modify time_coverage_start/end of output file:                                  |
            #_________________________________________________________________________________|

            # Read infile as netCDF dataset
            infile = netCDF4.Dataset(srcfile, 'r')
            # Read and write outfile as netCDF4 dataset
            outfile = netCDF4.Dataset(dstfile, 'r+')
        
            # Account for different file formats
            token = 0
            try:
                if 'nadir_view_time' in infile.groups['bin_attributes'].variables:
                    # Number of seconds at infile time_coverage_start/end
                    infile_start_sec = infile.groups['bin_attributes'].variables['nadir_view_time'][0]
                    # catch negative start time error
                    if (np.ma.is_masked(infile_start_sec)): raise ValueError("iFile contains a negative start time for the variable nadir_view_time.")
                    infile_end_sec = infile.groups['bin_attributes'].variables['nadir_view_time'][infile.dimensions['bins_along_track'].size - 1]
                    # Number of seconds at outfile time_coverage_start/end
                    outfile_start_sec = outfile.groups['bin_attributes'].variables['nadir_view_time'][0]
                    outfile_end_sec = outfile.groups['bin_attributes'].variables['nadir_view_time'][outfile.dimensions['bins_along_track'].size - 1] 

                    # Take infile time_coverage_start/end
                    infile_start_time = infile.time_coverage_start
                    infile_end_time = infile.time_coverage_end

                    # Extract year, month, day, hours, minutes, seconds from infile time coverage start/end
                    start_form = datetime.datetime.strptime(infile_start_time[0:19], '%Y-%m-%dT%H:%M:%S')
                    end_form = datetime.datetime.strptime(infile_end_time[0:19], '%Y-%m-%dT%H:%M:%S')
                    # Extract year,...,seconds from epoch
                    token = 1
            except ValueError:
                print("iFile contains a negative start time for the variable nadir_view_time.")
                sys.exit(1)
            except AttributeError:
                if 'time_nadir' in infile.groups['geolocation_data'].variables:
                    #infile.groups['geolocation_data']['time_nadir'].valid_max = 86400.00
                    #infile.groups['geolocation_data']['time_nadir'].valid_min = 0.00
           
                    # Number of seconds at infile time_coverage_start/end
                    infile_start_sec: float = infile.groups['geolocation_data'].variables['time_nadir'][0]
                    infile_end_sec = infile.groups['geolocation_data'].variables['time_nadir'][infile.dimensions['bins_along_track'].size - 1]
                    # Number of seconds at outfile time_coverage_start/end
                    outfile_start_sec = outfile.groups['geolocation_data'].variables['time_nadir'][0]
                    outfile_end_sec = outfile.groups['geolocation_data'].variables['time_nadir'][outfile.dimensions['bins_along_track'].size - 1] 

                    # Take infile time_coverage_start/end
                    infile_start_time = infile.time_coverage_start
                    infile_end_time = infile.time_coverage_end

                    # Extract year, month, day, hours, minutes, seconds from infile time coverage start/end
                    start_form = datetime.datetime.strptime(infile_start_time[0:19], '%Y-%m-%dT%H:%M:%S')
                    end_form = datetime.datetime.strptime(infile_end_time[0:19], '%Y-%m-%dT%H:%M:%S')
                    # Extract year,...,seconds from epoch
                    token = 1
            except:
                pass
            if token == 1:

                # Calculate the number of seconds contained in the time difference in previous step
                diff_sec_start = start_form.timestamp()
                # Calculate the number of seconds contained in the time difference in previous step
                diff_sec_end = end_form.timestamp()

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
                outfile.time_coverage_start = str(ostart_y) + '-' + str(ostart_mon) + '-' + str(ostart_d) + 'T' + str(ostart_h) + ':' + str(ostart_min) + ':' + str(ostart_s) + 'Z'
                outfile.time_coverage_end = str(oend_y) + '-' + str(oend_mon) + '-' + str(oend_d) + 'T' + str(oend_h) + ':' + str(oend_min) + ':' + str(oend_s) + 'Z'
            
            #_____________________________________________________________________________________________________________________________________
            #                                                                                 |
            # Gring Calculation:                                                              |
            #_________________________________________________________________________________|
            outfile.set_auto_mask(False)
            try:
                if 'latitude' in outfile.groups['geolocation_data'].variables:
                    bins_along_track = outfile.dimensions['bins_along_track'].size - 1
                    bins_across_track = outfile.dimensions['bins_across_track'].size - 1
                    
                    latitude = outfile.groups['geolocation_data'].variables['latitude']
                    longitude = outfile.groups['geolocation_data'].variables['longitude']
                    
                    lon_min = longitude[0, 0]
                    lon_max = longitude[bins_along_track, bins_across_track]
                    lat_min = latitude[0, 0]
                    lat_max = latitude[bins_along_track, bins_across_track]
                    
                    # Degrees to separate latitude points, here it is 20
                    lat_add = float((lat_max - lat_min) / 20)
                
                    # Place one point in the middle of longitude
                    lat_l = []
                    lat_r = []
                    lon_l = []
                    lon_r = []
        
                    # add latitude right values
                    for i in range(0, bins_along_track - 1, int(bins_along_track/lat_add)):
                        lat_r.append(i)
                        # add longitude left/right values
                        lon_l.append(0)
                        lon_r.append(bins_across_track) 
    
                    # add latitude left values
                    lat_l = list(reversed(lat_r))
                    
                    # add longitude up values
                    lon_u = [bins_across_track, (bins_across_track/2), 0] 
                    # add latitude up values
                    lat_u = [bins_along_track, bins_along_track, bins_along_track]
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
                        
                    # Format geospatial_bounds to WKT format
                    geospatial_bounds_str = "POLYGON((" + ', '.join([f"{lon} {lat}" for lon, lat in zip(g_lon, g_lat)]) + '))'
                    outfile.setncattr('geospatial_bounds', geospatial_bounds_str)

                    # Update geospatial lat/lon min/max
                    outfile.setncattr('geospatial_lat_max', str(max(g_lat)))
                    outfile.setncattr('geospatial_lat_min', str(min(g_lat)))
                    outfile.setncattr('geospatial_lon_max', str(max(g_lon)))
                    outfile.setncattr('geospatial_lon_min', str(min(g_lon)))

                    # Set gringpoints as global attributes and exclude last repeated pnts
                    outfile.setncattr('gringpointlatitude', ', '.join(map(str, g_lat[:-1])))
                    outfile.setncattr('gringpointlongitude', ', '.join(map(str, g_lon[:-1])))
                    outfile.setncattr('gringpointsequence', ', '.join(map(str, p_seq[:-1])))

            except:
                pass
  
        return 0

    def getpixlin(self):
        if self.verbose:
            print("north={} south={} west={} east={}".
                  format(self.north, self.south, self.west, self.east))

        # run lonlat2pixline
        pl = pixlin(file=self.ifile,
                    north=self.north, south=self.south,
                    west=self.west, east=self.east,
                    verbose=self.verbose)
        pl.lonlat2pixline(zero=False)  # using 1-based indices
        self.spixl, self.epixl, self.sline, self.eline = \
        (pl.spixl, pl.epixl, pl.sline, pl.eline)
        return pl.status

    def run(self):
        # convert to zero-based index
        self.spixl, self.epixl, self.sline, self.eline = \
        (v-1 for v in (self.spixl, self.epixl, self.sline, self.eline))

        # extract file
        subset = {'bins_across_track':[self.spixl, self.epixl],
                  'bins_along_track':   [self.sline, self.eline]}
        self.runextract(subset)

        #return status
        return 0

def chk_pixl(args, infile):
    if args.epixl == -1:
        args.epixl = infile.dimensions['bins_across_track'].size
    if args.eline == -1:
        args.eline = infile.dimensions['bins_along_track'].size
    return args.spixl, args.epixl, args.sline, args.eline
    
def chk_spex_width(args, infile):
        args.spixl = 245
        args.epixl = 273
        if args.eline == -1:
            args.eline = infile.dimensions['bins_along_track'].size
        return args.spixl, args.epixl, args.sline, args.eline
        
if __name__ == "__main__":
    print("l1cextract", versionStr)

    # parse command line
    parser = argparse.ArgumentParser(
        description='Extract specified area from OLCI Level 1C files.',
        epilog='Specify either geographic limits or pixel/line ranges, not both.')
    parser.add_argument('-v', '--verbose', help='print status messages',
                        action='store_true')
    parser.add_argument('ifile',
                        help='Level 1C input file')
    parser.add_argument('ofile',
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
    
    group3 = parser.add_argument_group('spex width (overwrites any pixel ranges or geographic limits)')
    group3.add_argument('--spex_width', help='spex width', action='store_true', default=None)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1) 
    args = parser.parse_args()
    
    # Read infile as netCDF dataset
    infile = netCDF4.Dataset(args.ifile, 'r')
    
    if args.spex_width == None:
        token = 0
    else:
        if args.spex_width != None:
            token = 1
    
    if token == 1:
        args.spixl, args.epixl, args.sline, args.eline = chk_spex_width(args, infile)
    else:
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
    goodindices = None not in (this.spixl, this.epixl, this.sline, this.eline)
    if (goodlatlons and goodindices):
        print("ERROR: Specify either geographic limits or pixel/line ranges, not both.")
        sys.exit(1)
    elif goodlatlons:
        status = this.getpixlin()
        if status not in (0, 110):
            print("No extract; lonlat2pixline status =", status)
            exit(status)
    elif goodindices:
        pass

    status = this.run()

    exit(status)
