#! /usr/bin/env python3

import seadasutils.anc_utils as anc_utils
import argparse
from shlex import quote
from modis.modis_utils import buildpcf, modis_env
import modis.modis_l1aextract_utils as ex
import modis.modis_GEO_utils as ga
from seadasutils.setupenv import env
import sys

if __name__ == "__main__":

    version = "1.1"

    # Read commandline options...

    parser = argparse.ArgumentParser(prog="modis_L1A_extract")
    parser.add_argument('--version', action='version', version='%(prog)s ' + version)

    parser.add_argument("filename", nargs='?',
                      help="Input L1A file to extract", metavar="L1AFILE")    
    parser.add_argument("-p", "--parfile", 
                      help="Parameter file containing program inputs", metavar="PARFILE")
    parser.add_argument("-o", "--output",
                      help="Output L1A extract filename - defaults to L1AFILE.sub", metavar="EXTRACTFILE")
    parser.add_argument("-g", "--geofile",
                      help="INPUT GEOFILE filename - defaults to basename of L1AFILE +'.GEO'", metavar="GEOFILE")
    parser.add_argument("-n", "--north", type=float,
                      help="Northernmost desired latitude", metavar="NORTH")
    parser.add_argument("-s", "--south", type=float,
                      help="Southernmost desired latitude", metavar="SOUTH")
    parser.add_argument("-w", "--west", type=float,
                      help="Westernmost desired longitude", metavar="WEST")
    parser.add_argument("-e", "--east", type=float,
                      help="Easternmost desired longitude", metavar="EAST")
    parser.add_argument("--extract_geo",
                      help="extract geolocation filename", metavar="EXGEO")
    parser.add_argument("--att1",
        help="Input attitude  filename 1 (chronological)", metavar="ATT1")
    parser.add_argument("--att2",
        help="Input attitude  filename 2 (chronological)", metavar="ATT2")
    parser.add_argument("--att3",
        help="Input attitude  filename 3 (chronological)", metavar="ATT3")
    parser.add_argument("--eph1",
        help="Input ephemeris filename 1 (chronological)", metavar="EPH1")
    parser.add_argument("--eph2",
        help="Input ephemeris filename 2 (chronological)", metavar="EPH2")
    parser.add_argument("--eph3",
        help="Input ephemeris filename 3 (chronological)", metavar="EPH3")
    ancdb_help_text = "Use a custom filename for ancillary database. If " \
                      "full path not given, ANCDB is assumed to exist "\
                      "(or will be created) under " + \
                      anc_utils.DEFAULT_ANC_DIR_TEXT + "/log/. If " + \
                      anc_utils.DEFAULT_ANC_DIR_TEXT + "/log/ does not " \
                      "exist, ANCDB is assumed (or will be created) " \
                      "under the current working directory"
    parser.add_argument("--ancdb", default='ancillary_data.db',help=ancdb_help_text, metavar="ANCDB")
    parser.add_argument("--ancdir",
        help="Use a custom directory tree for ancillary files", metavar="ANCDIR")
    parser.add_argument("-v", "--verbose", action="store_true",
                      default=False, help="print status messages")
    parser.add_argument("--log", action="store_true",
                      default=False, help="Save processing log file(s)")

    args = parser.parse_args()

    if args.parfile is None and args.filename is None:
        parser.print_help()
        sys.exit(1)

    parfile = None
    if args.parfile:
        parfile = quote(args.parfile)

    geofile = None
    if args.geofile:
        geofile = quote(args.geofile)

    ofile = "%s.sub" % args.filename
    if args.output:
        ofile = quote(args.output)

    m = ex.extract(filename=quote(args.filename),
                   parfile=parfile,
                   geofile=geofile,
                   outfile=ofile,
                   log=args.log,
                   north=args.north,
                   south=args.south,
                   west=args.west,
                   east=args.east,
                   verbose=args.verbose)

    env(m)
    modis_env(m)
    m.chk()
    status = m.run()
    if not status:
    # Create geolocation file for extract
        g = ga.modis_geo(filename=ofile,
                         geofile=args.extract_geo,
                         ancdb=args.ancdb,
                         ancdir=args.ancdir,
                         a1=args.att1,
                         a2=args.att2,
                         a3=args.att3,
                         e1=args.eph1,
                         e2=args.eph2,
                         e3=args.eph3,
                         log=args.log,
                         verbose=args.verbose)

        env(g)
        modis_env(g)
        g.atteph()
        buildpcf(g)

        g.run()
        sys.exit(0)
    else:
        sys.exit(status)

