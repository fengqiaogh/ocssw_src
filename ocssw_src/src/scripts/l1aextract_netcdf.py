#! /usr/bin/env python3

# Extractor for L1C files

import argparse
import datetime
from dateutil import tz
import netCDF4
import numpy as np
from os.path import basename
import pathlib
from shutil import copyfile as cp
import sys

from seadasutils.setupenv import env
from seadasutils.netcdf_utils import ncsubset_vars

versionStr = "1.0_2024-11-06"

class extract:

    def __init__(self, ifile, ofile=None,
                 spixl=None, epixl=None, sline=None, eline=None,
                 verbose=False):
        # inputs
        self.ifile = pathlib.Path(ifile)
        self.ofile = pathlib.Path(ofile)
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
            utc = tz.tzutc()
            scantime = outfile.variables['time'][:]
            stime = datetime.datetime.fromtimestamp(scantime[0],tz=utc)
            etime = datetime.datetime.fromtimestamp(scantime[-1],tz=utc)
            # Change outfile time_coverage_start/end
            outfile.time_coverage_start = str(stime.strftime('%Y-%m-%dT%H:%M:%S') + 'Z')
            outfile.time_coverage_end = str(etime.strftime('%Y-%m-%dT%H:%M:%S') + 'Z')

    def run(self):
        # convert to zero-based index
        self.spixl, self.epixl, self.sline, self.eline = \
        (v-1 for v in (self.spixl, self.epixl, self.sline, self.eline))

        # extract file
        subset = {'pixels':[self.spixl, self.epixl],
                  'scans':   [self.sline, self.eline]}
        infile = netCDF4.Dataset(self.ifile, 'r')
        dims = list(infile.dimensions.keys())
        infile.close()
        if 'nsamp' in dims:
            subset = {'nsamp':[self.spixl, self.epixl],
                      'lines':[self.sline, self.eline]}

        self.runextract(subset)

def chk_pixl(args, infile):
    if args.epixl == -1:
        args.epixl = infile.dimensions['pixels'].size
    if args.eline == -1:
        args.eline = infile.dimensions['scans'].size
    return args.spixl, args.epixl, args.sline, args.eline
        
if __name__ == "__main__":
    print("l1aextract_netcdf", versionStr)

    # parse command line
    parser = argparse.ArgumentParser(
        description='Extract specified area from an L1A netCDF file')
    parser.add_argument('-v', '--verbose', help='print status messages',
                        action='store_true')
    parser.add_argument('ifile',
                        help='SeaWiFS Level 1A input file')
    parser.add_argument('ofile',
                        help='output file')

    group2 = parser.add_argument_group('pixel/line ranges (1-based)')
    group2.add_argument('--spixl', type=int, help='start pixel', default = 1)
    group2.add_argument('--epixl', type=int, help='end pixel', default = -1)

    group2.add_argument('--sline', type=int, help='start line', default = 1)
    group2.add_argument('--eline', type=int, help='end line', default = -1)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1) 
    args = parser.parse_args()
    
    # Read infile as netCDF dataset
    infile = netCDF4.Dataset(args.ifile, 'r')
    
    args.spixl, args.epixl, args.sline, args.in_eline = chk_pixl(args, infile)
             
    # initialize
    this = extract(ifile=args.ifile,
                   ofile=args.ofile,
                   spixl=args.spixl,
                   epixl=args.epixl,
                   sline=args.sline,
                   eline=args.eline,
                   verbose=args.verbose)

    this.run()