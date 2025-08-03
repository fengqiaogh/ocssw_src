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

version_str = "1.1 (2025-06-09)"


def parse_manifest(manifestfile):
    tree = ET.parse(manifestfile)
    root = tree.getroot()
    radfiles = []
    tiefiles = []
    engfiles = []
    for child in root.find("dataObjectSection"):
        filename = pathlib.PurePath(
            child.find("byteStream").find("fileLocation").attrib["href"]
        ).name
        if "radiance" in filename:
            radfiles.append(filename)
        elif "tie_" in filename:
            tiefiles.append(filename)
        elif "time_coordinates" in filename:
            continue
        else:
            engfiles.append(filename)

    return (radfiles, tiefiles, engfiles)


def get_min_and_max(arr):
    return arr.min(), arr.max()


def epoch_2000(usec):
    # format Epoch 2000 time (microseconds since 2000-01-01)
    base = datetime(2000, 1, 1, 0, 0, 0)
    t = base + timedelta(microseconds=int(usec))
    return t.strftime("%Y-%m-%dT%H:%M:%S.%fZ")


class Extractor:

    def __init__(
        self,
        idir,
        odir=None,
        north=None,
        south=None,
        west=None,
        east=None,
        spixl=None,
        epixl=None,
        sline=None,
        eline=None,
        verbose=False,
    ):
        # inputs
        self.idir = pathlib.Path(idir)
        if odir is None:
            self.odir = self.idir / "subset"
        else:
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
        self.geofile = self.idir / "geo_coordinates.nc"
        self.timefile = self.idir / "time_coordinates.nc"
        self.tiefile = self.idir / "tie_geo_coordinates.nc"
        self.manifest = self.idir / "xfdumanifest.xml"

        # defaults
        self.runtime = None
        self.attrs = None

        # unused, but needed by setupenv.py
        self.dirs = {}
        self.ancdir = None
        self.curdir = False
        self.sensor = None
        env(self)  # run setupenv

    def extract(self, files, subset):
        # subset each file
        for filename in files:
            srcfile = self.idir / filename
            if srcfile.exists():
                dstfile = self.odir / filename
                if self.verbose:
                    print("Extracting", srcfile)

                ncsubset_vars(srcfile, dstfile, subset, timestamp=self.runtime)

                # update global attributes
                with netCDF4.Dataset(dstfile, mode="a") as dst:
                    dst.setncatts(self.attrs)
        return 0

    def get_indices_from_lat_lon(self):
        if self.verbose:
            print(
                "north={} south={} west={} east={}".format(
                    self.north, self.south, self.west, self.east
                )
            )

        # run lonlat2pixline
        pl = pixlin(
            geofile=self.geofile,
            north=self.north,
            south=self.south,
            west=self.west,
            east=self.east,
            verbose=self.verbose,
        )
        pl.lonlat2pixline(zero=False)  # using 1-based indices
        self.spixl, self.epixl, self.sline, self.eline = (
            pl.spixl,
            pl.epixl,
            pl.sline,
            pl.eline,
        )
        return pl.status

    def run(self):
        # convert to zero-based index
        self.spixl, self.epixl, self.sline, self.eline = (
            v - 1 for v in (self.spixl, self.epixl, self.sline, self.eline)
        )

        # check/create output directory
        if not self.odir:
            self.odir = ".".join([self.idir, "subset"])
        pathlib.Path(self.odir).mkdir(parents=True, exist_ok=True)

        # find tie file endpoints
        with netCDF4.Dataset(self.tiefile, "r") as src:
            npixl = src.dimensions["tie_columns"].size
            nline = src.dimensions["tie_rows"].size
            dpixl = getattr(src, "ac_subsampling_factor", 1)  # tie_col_pts
            dline = getattr(src, "al_subsampling_factor", 1)  # tie_row_pts
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
        spixl, epixl = [dpixl * v for v in (spixl, epixl)]
        sline, eline = [dline * v for v in (sline, eline)]
        if self.verbose:
            print(
                "spixl={} epixl={} sline={} eline={}".format(
                    spixl + 1, epixl + 1, sline + 1, eline + 1
                )
            )

        # find new start, stop times
        with netCDF4.Dataset(self.timefile, "r") as src:
            ts = src["time_stamp"][[sline, eline]]
            start_time = epoch_2000(ts[0])
            stop_time = epoch_2000(ts[1])

        # find new lat/lon ranges
        with netCDF4.Dataset(self.geofile, "r") as src:
            lat_min, lat_max = get_min_and_max(
                src["latitude"][sline:eline, spixl:epixl]
            )
            lon_min, lon_max = get_min_and_max(
                src["longitude"][sline:eline, spixl:epixl]
            )

        # define global attributes
        self.attrs = {
            "start_time": start_time,
            "stop_time": stop_time,
            "geospatial_lat_min": lat_min,
            "geospatial_lat_max": lat_max,
            "geospatial_lon_min": lon_min,
            "geospatial_lon_max": lon_max,
        }
        self.runtime = time.gmtime()  # same for all files

        # extract time_coordinates
        subset = {"rows": [sline, eline]}
        status = self.extract([self.timefile.name], subset)

        # extract full-resolution files
        subset = {"columns": [spixl, epixl], "rows": [sline, eline]}
        (radfiles, tiefiles, engfiles) = parse_manifest(self.manifest)
        status += self.extract(radfiles + engfiles, subset)

        # extract lower-resolution (tie) files
        subset = {
            "tie_columns": [spixl, epixl] // dpixl,
            "tie_rows": [sline, eline] // dline,
        }
        status += self.extract(tiefiles, subset)

        return status

    def verify_input(self) -> str:
        # file checks
        if not self.idir.exists():
            return "ERROR: Directory '{}' does not exist. Exiting.".format(self.idir)
        if not self.timefile.exists():
            return "ERROR: Timestamp file ({}) not found!".format(self.timefile)
        if not self.geofile.exists():
            return "ERROR: Geolocation file ({}) not found!".format(self.geofile)
        if not self.tiefile.exists():
            return "ERROR: Tie file ({}) not found!".format(self.tiefile)

        # input value checks
        using_geospatial_coords = None not in (
            self.north,
            self.south,
            self.west,
            self.east,
        )
        using_cartesian_coords = None not in (
            self.spixl,
            self.epixl,
            self.sline,
            self.eline,
        )
        if using_geospatial_coords and using_cartesian_coords:
            return "ERROR: Specify either geographic limits or pixel/line ranges, not both."
        elif using_geospatial_coords and not using_cartesian_coords:
            status = self.get_indices_from_lat_lon()
            if status not in (0, 110):
                return "No extract; lonlat2pixline status ={}".format(status)
        elif using_cartesian_coords and not using_geospatial_coords:
            for coord in (self.spixl, self.epixl, self.sline, self.eline):
                if coord < 1:
                    return "ERROR: Pixel/line ranges must be >= 1"
            pass
        else:
            return "ERROR: Specify all values for either geographic limits or pixel/line ranges."

        return None  # All good


if __name__ == "__main__":
    print("l1bextract_safe_nc", version_str)

    # parse command line
    parser = argparse.ArgumentParser(
        description="Extract specified area from OLCI Level 1B files.",
        epilog="Specify either geographic limits or pixel/line ranges, not both.",
    )
    parser.add_argument(
        "-v", "--verbose", help="print status messages", action="store_true"
    )
    parser.add_argument("idir", help="directory containing OLCI Level 1B files")
    parser.add_argument(
        "odir", help='output directory (defaults to "idir.subset")', nargs="?"
    )

    group1 = parser.add_argument_group("geographic limits")
    group1.add_argument("-n", "--north", type=float, help="northernmost latitude")
    group1.add_argument("-s", "--south", type=float, help="southernmost latitude")
    group1.add_argument("-w", "--west", type=float, help="westernmost longitude")
    group1.add_argument("-e", "--east", type=float, help="easternmost longitude")

    group2 = parser.add_argument_group("pixel/line ranges (1-based)")
    group2.add_argument("--spixl", type=int, help="start pixel")
    group2.add_argument("--epixl", type=int, help="end pixel")
    group2.add_argument("--sline", type=int, help="start line")
    group2.add_argument("--eline", type=int, help="end line")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    # initialize
    extractor = Extractor(
        idir=args.idir,
        odir=args.odir,
        north=args.north,
        south=args.south,
        west=args.west,
        east=args.east,
        spixl=args.spixl,
        epixl=args.epixl,
        sline=args.sline,
        eline=args.eline,
        verbose=args.verbose,
    )

    error_message = extractor.verify_input()
    if error_message is not None:
        print(error_message)
        sys.exit(1)

    # run
    status = extractor.run()

    # copy the manifest in case we ever need it
    cp(extractor.manifest, extractor.odir / "xfdumanifest.xml")

    exit(status)
