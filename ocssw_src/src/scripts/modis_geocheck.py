#! /usr/bin/env python3

import argparse
import sys
import os
from modis.modis_utils import modis_env
import modis.modis_GEO_utils as modisGEO
from seadasutils.setupenv import env


if __name__ == "__main__":

    version = "1.1"

    # Read commandline options...
    parser = argparse.ArgumentParser(prog="modis_atteph")
    parser.add_argument('--version', action='version', version='%(prog)s ' + version)
    parser.add_argument("geofile", nargs='?',
                      help="Input GEO file", metavar="GEOFILE")  
    parser.add_argument("--threshold", dest='geothresh', default=95, type=float,
                      help="percentage of geo-populated pixels required to pass geocheck validation test", metavar="THRESHOLD")
    parser.add_argument("-v", "--verbose", action="store_true",
                      default=False, help="print status messages")

    args = parser.parse_args()
    if args.geofile:
        if not os.path.exists(args.geofile):
            print ("*** ERROR: Provided geolocation file does not exist.")
            print ("*** Validation test failed for geolocation file:", geofile)
            sys.exit(1)
    else:
        parser.print_help()
        sys.exit(1)

    # kluge: use geofile as l1afile for setup
    m = modisGEO.modis_geo(filename=args.geofile, geofile=args.geofile,
                           geothresh=args.geothresh,
                           verbose=args.verbose
    )
    env(m)
    modis_env(m)
    m.geochk()
    sys.exit(0)
