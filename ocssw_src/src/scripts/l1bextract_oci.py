#! /usr/bin/env python3

# Extractor for L1B_OCI files

import argparse
import pathlib
import sys
from datetime import datetime, timedelta

import netCDF4
import numpy as np

from seadasutils.netcdf_utils import ncsubset_vars
from seadasutils.setupenv import env

versionStr = "1.3_2025-10-17"


class extract:
    def __init__(
        self,
        ifile,
        ofile=None,
        spixl=None,
        epixl=None,
        sline=None,
        eline=None,
        verbose=False,
    ):
        # inputs
        self.ifile = pathlib.Path(ifile)
        self.ofile = pathlib.Path(ofile)
        self.spixl = spixl
        self.epixl = epixl
        self.sline = sline
        self.eline = eline
        self.verbose = verbose

        # unused, but needed by setupenv.py
        self.dirs = {}
        self.ancdir = None
        self.curdir = False
        env(self)  # run setupenv

    def runextract(self):
        #
        # open input file
        infile = netCDF4.Dataset(self.ifile, "r")

        # adjust end line/pixel as needed
        dims = infile.dimensions
        npixl = dims["pixels"].size
        nline = dims["scans"].size
        if self.epixl == -1:
            self.epixl = npixl
        else:
            self.epixl = min(self.epixl, npixl)
        if self.eline == -1:
            self.eline = nline
        else:
            self.eline = min(self.eline, nline)

        # convert to zero-based index
        self.spixl, self.epixl, self.sline, self.eline = (
            v - 1 for v in (self.spixl, self.epixl, self.sline, self.eline)
        )

        # extract region
        subset = {
            "pixels": [self.spixl, self.epixl],
            "scans": [self.sline, self.eline],
        }
        ncsubset_vars(self.ifile, self.ofile, subset, verbose=self.verbose)

        # reopen output file to update metadata
        outfile = netCDF4.Dataset(self.ofile, "r+")

        # update temporal extents
        timevar = outfile.groups["scan_line_attributes"].variables["time"]
        daystart = datetime.strptime(timevar.units, "seconds since %Y-%m-%d")
        stime = daystart + timedelta(seconds=np.min(timevar))
        etime = daystart + timedelta(seconds=np.max(timevar))
        outfile.time_coverage_start = stime.isoformat(timespec="milliseconds") + "Z"
        outfile.time_coverage_end = etime.isoformat(timespec="milliseconds") + "Z"

        # find new spatial extents
        geovars = outfile.groups["geolocation_data"].variables
        lat = geovars["latitude"]
        lon = geovars["longitude"]

        south_lat = np.min(lat)
        north_lat = np.max(lat)
        west_lon = np.min(lon)
        east_lon = np.max(lon)
        if east_lon - west_lon > 180:  # granule crosses 180
            west_lon = np.min(lon[lon >= 0])  # western part in eastern hemisphere
            east_lon = np.max(lon[lon < 0])  # eastern part in western hemisphere

        def corners(a):
            return [a[-1, -1], a[-1, 0], a[0, 0], a[0, -1], a[-1, -1]]

        p_lat = corners(lat)
        p_lon = corners(lon)
        polygon = ", ".join([f"{lon:.5f} {lat:.5f}" for lon, lat in zip(p_lon, p_lat)])

        # ...update spatial global attributes
        outfile.geospatial_lat_min = south_lat
        outfile.geospatial_lat_max = north_lat
        outfile.geospatial_lon_min = west_lon
        outfile.geospatial_lon_max = east_lon
        outfile.geospatial_bounds = f"POLYGON(({polygon}))"  #  WKT format

        # add extract_pixel/line_start/stop (1-based)
        outfile.extract_pixel_start = np.int32(self.spixl + 1)
        outfile.extract_pixel_stop = np.int32(self.epixl + 1)
        outfile.extract_line_start = np.int32(self.sline + 1)
        outfile.extract_line_stop = np.int32(self.eline + 1)
        if "extract_pixel_start" in infile.ncattrs():
            outfile.extract_pixel_start += np.int32(infile.extract_pixel_start - 1)
            outfile.extract_pixel_stop += np.int32(infile.extract_pixel_start - 1)
            outfile.extract_line_start += np.int32(infile.extract_line_start - 1)
            outfile.extract_line_stop += np.int32(infile.extract_line_start - 1)


# // end class extract


if __name__ == "__main__":
    print("l1bextract_oci", versionStr)

    # parse command line
    parser = argparse.ArgumentParser(
        description="Extract specified area from OCI Level 1B files"
    )
    parser.add_argument(
        "-v", "--verbose", help="print status messages", action="store_true"
    )
    parser.add_argument("ifile", help="Level 1B input file")
    parser.add_argument("ofile", help="extracted L1B file")

    group2 = parser.add_argument_group("pixel/line ranges (1-based)")
    group2.add_argument("--spixl", type=int, help="start pixel", default=1)
    group2.add_argument("--epixl", type=int, help="end pixel", default=-1)

    group2.add_argument("--sline", type=int, help="start line", default=1)
    group2.add_argument("--eline", type=int, help="end line", default=-1)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    # initialize
    this = extract(
        ifile=args.ifile,
        ofile=args.ofile,
        spixl=args.spixl,
        epixl=args.epixl,
        sline=args.sline,
        eline=args.eline,
        verbose=args.verbose,
    )

    # run extracts and update metadata
    this.runextract()
