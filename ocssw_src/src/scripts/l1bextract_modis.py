#! /usr/bin/env python3

# Extractor for MODIS Collection 7 L1B and GEO files

import argparse
import pathlib
import sys
from datetime import datetime, timezone

import netCDF4
import numpy as np

from seadasutils.netcdf_utils import ncsubset_vars
from seadasutils.setupenv import env

versionStr = "1.0.1_2025-10-17"


class extract:
    def __init__(
        self,
        ifile,
        ofile=None,
        ifilegeo=None,
        ofilegeo=None,
        spixl=None,
        epixl=None,
        sline=None,
        eline=None,
        verbose=False,
    ):
        # inputs
        self.ifile = pathlib.Path(ifile)
        self.ofile = pathlib.Path(ofile)
        self.ifilegeo = pathlib.Path(ifilegeo)
        self.ofilegeo = pathlib.Path(ofilegeo)
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
        # open input files
        l1b_in = netCDF4.Dataset(self.ifile, "r")
        geo_in = netCDF4.Dataset(self.ifilegeo, "r")

        # TODO: verify that L1B and GEO match

        # adjust end line/pixel as needed
        dims = l1b_in["HDFEOS/SWATHS/MODIS_SWATH_Type_L1B"].dimensions
        npixl = dims["Max_EV_frames"].size
        nline = dims["10*nscans"].size
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

        # preserve integer number of scans
        pixlrange = [self.spixl, self.epixl]
        scanrange = [int(self.sline / 10), int(self.eline / 10)]
        linerange = [10 * scanrange[0], 10 * (scanrange[1] + 1) - 1]
        self.sline, self.eline = linerange

        # extract L1B
        subset = {
            "Max_EV_frames": pixlrange,
            "10*nscans": linerange,
            "2*nscans": [2 * scanrange[0], 2 * (scanrange[1] + 1) - 1],
            "1KM_geo_dim": [int(self.spixl / 5), int(self.epixl / 5)],
            "phony_dim_9": scanrange,
            "phony_dim_12": scanrange,
        }
        if self.verbose:
            print(f"L1B {subset=}")
        ncsubset_vars(self.ifile, self.ofile, subset, verbose=self.verbose)

        # extract GEO
        subset = {
            "mframes": pixlrange,
            "mframes*2": [2 * pixlrange[0], 2 * (pixlrange[1] + 1) - 1],
            "nscans*10": linerange,
            "nscans*20": [2 * linerange[0], 2 * (linerange[1] + 1) - 1],
            "phony_dim_5": scanrange,
            "phony_dim_11": pixlrange,
        }
        if self.verbose:
            print(f"GEO {subset=}")
        ncsubset_vars(self.ifilegeo, self.ofilegeo, subset, verbose=self.verbose)

        # reopen output files to update metadata
        l1b_out = netCDF4.Dataset(self.ofile, "r+")
        geo_out = netCDF4.Dataset(self.ofilegeo, "r+")

        # ...update dimensional global attributes
        nscans = np.int32(scanrange[1] - scanrange[0] + 1)
        nframes = np.int32(pixlrange[1] - pixlrange[0] + 1)
        l1b_out.setncattr("Number of Scans", nscans)
        geo_out.setncattr("Number of Scans", nscans)
        l1b_out.setncattr("Max Earth View Frames", nframes)
        geo_out.setncattr("Max Earth Frames", nframes)

        # find new temporal extents
        tai93 = np.take(geo_out.variables["EV start time"], [0, -1])
        tai93_offset = datetime(1993, 1, 1) - datetime(1970, 1, 1)
        times = [  # (convert TAI93 to UTC but ignore leap seconds)
            datetime.fromtimestamp(t, tz=timezone.utc) + tai93_offset for t in tai93
        ]

        # ...update temporal global attributes
        for outfile in [l1b_out, geo_out]:
            outfile.RangeBeginningDate = times[0].strftime("%Y-%m-%d")
            outfile.RangeBeginningTime = times[0].strftime("%H:%M:%S.%f")
            outfile.RangeEndingDate = times[1].strftime("%Y-%m-%d")
            outfile.RangeEndingTime = times[1].strftime("%H:%M:%S.%f")

        # find new spatial extents
        group = geo_out["HDFEOS/SWATHS/MODIS_Swath_Type_GEO/Geolocation Fields"]
        lat = group["Latitude"][:]
        lon = group["Longitude"][:]

        south_lat = np.double(np.min(lat))
        north_lat = np.double(np.max(lat))
        west_lon = np.min(lon)
        east_lon = np.max(lon)
        if east_lon - west_lon > 180:  # granule crosses 180
            west_lon = np.min(lon[lon >= 0])  # western part in eastern hemisphere
            east_lon = np.max(lon[lon < 0])  # eastern part in western hemisphere
        west_lon = np.double(west_lon)
        east_lon = np.double(east_lon)

        def corners(a):
            vals = [a[0, 0], a[0, -1], a[-1, -1], a[-1, 0]]
            return [np.double(v) for v in vals]

        # ...update spatial global attributes
        geo_out.NorthBoundingCoordinate = north_lat
        geo_out.SouthBoundingCoordinate = south_lat
        geo_out.EastBoundingCoordinate = east_lon
        geo_out.WestBoundingCoordinate = west_lon
        for outfile in [l1b_out, geo_out]:
            outfile.geospatial_lat_min = south_lat
            outfile.geospatial_lat_max = north_lat
            outfile.geospatial_lon_min = west_lon
            outfile.geospatial_lon_max = east_lon
            outfile.GRingPointLatitude = corners(lat)
            outfile.GRingPointLongitude = corners(lon)

        # add extract_pixel/line_start/stop (1-based)
        for infile, outfile in zip([l1b_in, geo_in], [l1b_out, geo_out]):
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
    print("l1bextract_modis", versionStr)

    # parse command line
    parser = argparse.ArgumentParser(
        description="Extract specified area from a MODIS L1B/GEO NetCDF4 or HDF5 file pair"
    )
    parser.add_argument(
        "-v", "--verbose", help="print status messages", action="store_true"
    )
    parser.add_argument("ifile", help="Level 1B input file")
    parser.add_argument("ofile", help="extracted L1B file")
    parser.add_argument("ifilegeo", help="Geolocation input file")
    parser.add_argument("ofilegeo", help="extracted GEO file")

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
        ifilegeo=args.ifilegeo,
        ofilegeo=args.ofilegeo,
        spixl=args.spixl,
        epixl=args.epixl,
        sline=args.sline,
        eline=args.eline,
        verbose=args.verbose,
    )

    # run extracts and update metadata
    this.runextract()
