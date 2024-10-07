#! /usr/bin/env python3

# Extractor for OLCI Sentinel 3 and MERIS SAFE formatted L1B files

import argparse
from datetime import datetime, timedelta
import netCDF4
import pathlib
import sys
import time
from shutil import copyfile as cp
import xml.etree.ElementTree as ET

from seadasutils.netcdf_utils import ncsubset_vars
from seadasutils.pixlin_utils import pixlin
from seadasutils.setupenv import env

versionStr = "1.0 (2023-05-04)"

def parseManifest(manifestfile):
    tree = ET.parse(manifestfile)
    root = tree.getroot()
    radfiles = []
    tiefiles = []
    engfiles = []
    for child in root.find('dataObjectSection'):
        filename = pathlib.PurePath(child.find('byteStream').find('fileLocation').attrib['href']).name
        if 'radiance' in filename:
            radfiles.append(filename)
        elif 'tie_' in filename:
            tiefiles.append(filename)
        elif 'time_coordinates' in filename:
            continue
        else:
            engfiles.append(filename)

    return (radfiles, tiefiles, engfiles)

def minmax(arr):
    return arr.min(), arr.max()

def epoch2000(usec):
    # format Epoch 2000 time (microseconds since 2000-01-01)
    base = datetime(2000, 1, 1, 0, 0, 0)
    t = base + timedelta(microseconds=int(usec))
    return t.strftime('%Y-%m-%dT%H:%M:%S.%fZ')

class extract:

    def __init__(self, idir, odir=None,
                 north=None, south=None, west=None, east=None,
                 spixl=None, epixl=None, sline=None, eline=None,
                 verbose=False):
        # inputs
        self.idir = pathlib.Path(idir)
        self.odir = pathlib.Path(odir)
        self.north = north
        self.south = south
        self.west = west
        self.east = east
        self.spixl = spixl
        self.epixl = epixl
        self.sline = sline
        self.eline = eline
        self.verbose = verbose
        self.geofile = self.idir / 'geo_coordinates.nc'
        self.timefile = self.idir / 'time_coordinates.nc'
        self.tiefile = self.idir / 'tie_geo_coordinates.nc'
        self.manifest = self.idir / 'xfdumanifest.xml'

        # defaults
        self.runtime = None
        self.attrs = None

        # unused, but needed by setupenv.py
        self.dirs = {}
        self.ancdir = None
        self.curdir = False
        self.sensor = None
        env(self)  # run setupenv

    def runextract(self, files, subset):
        # subset each file
        for filename in files:
            srcfile = self.idir / filename
            if srcfile.exists():
                dstfile = self.odir / filename
                if self.verbose:
                    print('Extracting', srcfile)

                ncsubset_vars(srcfile, dstfile, subset, timestamp=self.runtime)

                # update global attributes
                with netCDF4.Dataset(dstfile, mode='a') as dst:
                    dst.setncatts(self.attrs)
        return 0

    def getpixlin(self):
        if self.verbose:
            print("north={} south={} west={} east={}".
                  format(self.north, self.south, self.west, self.east))

        # run lonlat2pixline
        pl = pixlin(geofile=self.geofile,
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

        # check/create output directory
        if not self.odir:
            self.odir = '.'.join([self.idir, 'subset'])
        pathlib.Path(self.odir).mkdir(parents=True, exist_ok=True)

        # find tie file endpoints
        with netCDF4.Dataset(self.tiefile, 'r') as src:
            npixl = src.dimensions['tie_columns'].size
            nline = src.dimensions['tie_rows'].size
            dpixl = getattr(src, 'ac_subsampling_factor', 1)  # tie_col_pts
            dline = getattr(src, 'al_subsampling_factor', 1)  # tie_row_pts
        spixl, epixl = [self.spixl, self.epixl + dpixl - 1] // dpixl
        sline, eline = [self.sline, self.eline + dline - 1] // dline

        # Make sure tie files have enough points for interpolation.
        # l1_olci.c uses type gsl_interp_cspline, which requires
        # at least 3 pixels in each dimension.
        points_required = 3

        pixls_needed = points_required - (epixl - spixl + 1)
        if pixls_needed > 0:
            spixl = max(0, spixl - pixls_needed // 2)
            epixl = min(spixl + points_required, npixl) - 1
            if epixl == npixl - 1:
                spixl = npixl - points_required

        lines_needed = points_required - (eline - sline + 1)
        if lines_needed > 0:
            sline = max(0, sline - lines_needed // 2)
            eline = min(sline + points_required, nline) - 1
            if eline == nline - 1:
                sline = nline - points_required

        # scale to full-resolution coordinates
        # will be aligned with tie files
        spixl, epixl = [dpixl*v for v in (spixl, epixl)]
        sline, eline = [dline*v for v in (sline, eline)]
        if self.verbose:
            print("spixl={} epixl={} sline={} eline={}".
                  format(spixl+1, epixl+1, sline+1, eline+1))

        # find new start, stop times
        with netCDF4.Dataset(self.timefile, 'r') as src:
            ts = src['time_stamp'][[sline, eline]]
            start_time = epoch2000(ts[0])
            stop_time = epoch2000(ts[1])

        # find new lat/lon ranges
        with netCDF4.Dataset(self.geofile, 'r') as src:
            lat_min, lat_max = minmax(src['latitude']
                                    [sline:eline, spixl:epixl])
            lon_min, lon_max = minmax(src['longitude']
                                    [sline:eline, spixl:epixl])

        # define global attributes
        self.attrs = {'start_time': start_time,
                      'stop_time':  stop_time,
                      'geospatial_lat_min': lat_min,
                      'geospatial_lat_max': lat_max,
                      'geospatial_lon_min': lon_min,
                      'geospatial_lon_max': lon_max }
        self.runtime = time.gmtime()  # same for all files

        # extract time_coordinates
        subset = {'rows':   [sline, eline]}
        status = self.runextract([self.timefile.name], subset)

        # extract full-resolution files
        subset = {'columns':[spixl, epixl],
                  'rows':   [sline, eline]}
        (radfiles, tiefiles, engfiles) = parseManifest(self.manifest)
        status += self.runextract(radfiles + engfiles, subset)

        # extract lower-resolution (tie) files
        subset = {'tie_columns':[spixl, epixl] // dpixl,
                  'tie_rows':   [sline, eline] // dline}
        status += self.runextract(tiefiles, subset)

        return status


if __name__ == "__main__":
    print("l1bextract_safe_nc", versionStr)

    # parse command line
    parser = argparse.ArgumentParser(
        description='Extract specified area from OLCI Level 1B files.',
        epilog='Specify either geographic limits or pixel/line ranges, not both.')
    parser.add_argument('-v', '--verbose', help='print status messages',
                        action='store_true')
    parser.add_argument('idir',
                        help='directory containing OLCI Level 1B files')
    parser.add_argument('odir', nargs='?',
                        help='output directory (defaults to "idir.subset")')

    group1 = parser.add_argument_group('geographic limits')
    group1.add_argument('-n', '--north', type=float, help='northernmost latitude')
    group1.add_argument('-s', '--south', type=float, help='southernmost latitude')
    group1.add_argument('-w', '--west', type=float, help='westernmost longitude')
    group1.add_argument('-e', '--east', type=float, help='easternmost longitude')

    group2 = parser.add_argument_group('pixel/line ranges (1-based)')
    group2.add_argument('--spixl', type=int, help='start pixel')
    group2.add_argument('--epixl', type=int, help='end pixel')
    group2.add_argument('--sline', type=int, help='start line')
    group2.add_argument('--eline', type=int, help='end line')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    # initialize
    this = extract(idir=args.idir,
                   odir=args.odir,
                   north=args.north,
                   south=args.south,
                   west=args.west,
                   east=args.east,
                   spixl=args.spixl,
                   epixl=args.epixl,
                   sline=args.sline,
                   eline=args.eline,
                   verbose=args.verbose)

    # file checks
    if not this.idir.exists():
        print("ERROR: Directory '{}' does not exist. Exiting.".format(this.idir))
        sys.exit(1)
    if not this.timefile.exists():
        print("ERROR: Timestamp file ({}) not found!".format(this.timefile))
        sys.exit(1)
    if not this.geofile.exists():
        print("ERROR: Geolocation file ({}) not found!".format(this.geofile))
        sys.exit(1)
    if not this.tiefile.exists():
        print("ERROR: Tie file ({}) not found!".format(this.tiefile))
        sys.exit(1)

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
    else:
        print("ERROR: Specify all values for either geographic limits or pixel/line ranges.")
        sys.exit(1)

    # run
    status = this.run()

    # copy the manifest in case we ever need it
    cp(this.manifest, this.odir / 'xfdumanifest.xml')

    exit(status)
