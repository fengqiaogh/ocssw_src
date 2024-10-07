#! /usr/bin/env python3
"""
Wrapper program to produce MODIS GEO files.
"""
import argparse
from shlex import quote
import sys
import seadasutils.anc_utils as anc_utils
from modis.modis_utils import buildpcf, modis_env
import modis.modis_GEO_utils as modisGEO
from seadasutils.setupenv import env


def main():
    """
    Driver function for the program.
    """
    version = "1.1"

    # Read commandline options...
    parser = argparse.ArgumentParser(prog="modis_GEO")
    parser.add_argument('--version', action='version', version='%(prog)s ' + version)
    parser.add_argument("filename", nargs='?',
                      help="Input L1A file", metavar="L1AFILE")  
    parser.add_argument("-p", "--parfile",
                      help="Parameter file containing program inputs", metavar="PARFILE")
    parser.add_argument("-o", "--output",
                      help="Output GEO filename - defaults to '(A|T)YYYYDDDHHMMSS.GEO'", metavar="GEOFILE")
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
    parser.add_argument("--ancdb",
        help=ancdb_help_text, metavar="ANCDB",default='ancillary_data.db')
    parser.add_argument("--ancdir",
        help="Use a custom directory tree for ancillary files", metavar="ANCDIR")
    parser.add_argument("-c", "--curdir", action="store_true",
        default=False, help="Download ancillary files directly into current working directory")
    parser.add_argument("-r", "--refreshDB", action="store_true", default=False,
                      help="Remove existing database records and re-query for ancillary files")
    parser.add_argument("-f", "--force-download", action="store_true", dest='force', default=False,
                      help="Force download of ancillary files, even if found on hard disk")
    parser.add_argument("--disable-download", action="store_false", default=True,
                      help="Disable download of ancillary files not found on hard disk")
    parser.add_argument("--timeout", type=float, default=10.0, metavar="TIMEOUT",
                      help="set the network timeout in seconds")
    parser.add_argument("--threshold", type=float, default=95.0,
                      help="percentage of geo-populated pixels required to pass geocheck validation test",
                      metavar="THRESHOLD")
    parser.add_argument("-d", "--enable-dem", action="store_true", default=False,
                      help="Enable MODIS terrain elevation correction")
    parser.add_argument("-v", "--verbose", action="count",
                      default=0, help="print status messages")
    parser.add_argument("--log", action="store_true",
                      default=False, help="Save processing log file(s)")


    args = parser.parse_args()

    if args.parfile is None and args.filename is None:
        parser.print_help()
        sys.exit(1)

    parfile = None
    if args.parfile:
        parfile = quote(args.parfile)

    outputfile = None
    if args.output:
        outputfile = quote(args.output)
    

    m = modisGEO.modis_geo(filename=quote(args.filename),
                           parfile=parfile,
                           geofile=outputfile,
                           a1=args.att1,
                           a2=args.att2,
                           a3=args.att3,
                           e1=args.eph1,
                           e2=args.eph2,
                           e3=args.eph3,
                           terrain=args.enable_dem,
                           geothresh=args.threshold,
                           ancdir=args.ancdir,
                           curdir=args.curdir,
                           ancdb=args.ancdb,
                           refreshDB=args.refreshDB,
                           forcedl=args.force,
                           download=args.disable_download,
                           log=args.log,
                           verbose=args.verbose,
                           timeout=args.timeout)

    env(m)
    modis_env(m)
    m.chk()
    m.utcleap()
    try:
        m.atteph()
    except SystemExit:
        print ("Failed to identify/retrieve attitude/ephemeris required to geolocate %s; exiting." % args.filename)
        raise
    buildpcf(m)
    m.run()
    return 0

if __name__ == "__main__":
    sys.exit(main())
