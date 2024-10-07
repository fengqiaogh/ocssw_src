#! /usr/bin/env python3

import argparse
import os
import sys
import seadasutils.anc_utils as ga
from seadasutils.setupenv import env
from seadasutils.ProcUtils import check_sensor

if __name__ == "__main__":
    version = "2.2"
    
    # Read commandline options...
    usage = '''
    %(prog)s L1A_file
             or
    %(prog)s -m modisa|aqua|modist|terra -s YYYYDDDHHMMSS -e YYYYDDDHHMMSS

    NOTE: Currently NO2 climatological data is used for OBPG operational
          processing, so to match OBPG distributed data products, the default
          behaviour disables NO2 searching.

    This program queries an OBPG server and optionally downloads the optimal
    ancillary data files for modis_GEO processing. If an input file
    is specified the start and end times are determined automatically, otherwise
    a start time must be provided by the user.

    A text file (with the extension '.atteph') is created containing parameters
    that can be directly used for optimal optimal modis_GEO processing

    EXIT STATUS:
        0 - all is well in the world
        1 - predicted attitude selected
        2 - predicted ephemeris selected
        4 - no attitude found
        8 - no ephemeris found
        16 - invalid mission
        12 - no att/eph files currently exist corresponding to the start
             time and therefore no .atteph parameter text file was created
        99 - an error was encountered; no .atteph parameter text file was created
    '''
    
    parser = argparse.ArgumentParser(prog="modis_atteph",usage=usage)
    parser.add_argument('--version', action='version', version='%(prog)s ' + version)
    parser.add_argument("filename", nargs='?',
                      help="Input L1 file", metavar="L1FILE")  
    parser.add_argument("-s", "--start",
                      help="Granule start time (YYYYDDDHHMMSS)",
                      metavar="START")
    parser.add_argument("-e", "--stop",
                      help="Granule stop time (YYYYDDDHHMMSS)", metavar="STOP")

    ancdb_help_text = "Use a custom file for ancillary database. If full " \
                      "path not given, ANCDB is assumed to exist (or " \
                      "will be created) under " + ga.DEFAULT_ANC_DIR_TEXT + \
                      "/log/. If " + ga.DEFAULT_ANC_DIR_TEXT + "/log/ does " \
                                                               "not exist, ANCDB is assumed (or will be created) " \
                                                               " under the current working directory"

    ancdb_help_text = "Use a custom filename for ancillary database. If " \
                      "full path not given, ANCDB is assumed to exist "\
                      "(or will be created) under " + \
                      ga.DEFAULT_ANC_DIR_TEXT + "/log/. If " + \
                      ga.DEFAULT_ANC_DIR_TEXT + "/log/ does not " \
                      "exist, ANCDB is assumed (or will be created) " \
                      "under the current working directory"
    parser.add_argument("--ancdb", default='ancillary_data.db',help=ancdb_help_text, metavar="ANCDB")
    parser.add_argument("--ancdir",
        help="Use a custom directory tree for ancillary files", metavar="ANCDIR")
    parser.add_argument("-c", "--curdir", action="store_true",
        default=False, help="Download ancillary files directly into current working directory")
    parser.add_argument("-r", "--refreshDB", action="store_true", default=False,
                      help="Remove existing database records and re-query for ancillary files")
    parser.add_argument("--disable-download", action="store_false", dest="download",default=True,
                      help="Disable download of ancillary files not found on hard disk")
    parser.add_argument("--timeout", type=float, default=10.0, metavar="TIMEOUT",
                      help="set the network timeout in seconds")
    parser.add_argument("-m", "--mission", help="MODIS mission - A(qua) or T(erra)", metavar="MISSION")

    parser.add_argument("-v", "--verbose", action="count",
                      default=0, help="print status messages")
    parser.add_argument("--noprint", action="store_false", default=True,
                      help="Suppress printing the resulting list of files to the screen")
    parser.add_argument("-f", "--force-download", action="store_true", dest='force', default=False,
                      help="Force download of ancillary files, even if found on hard disk")

    args = parser.parse_args()
    if args.filename is None and (args.mission and args.start is None):
        parser.print_help()
        sys.exit(32)


    m = ga.getanc(filename=args.filename,
                  start=args.start,
                  stop=args.stop,
                  ancdir=args.ancdir,
                  curdir=args.curdir,
                  ancdb=args.ancdb,
                  sensor=args.mission,
                  download=args.download,
                  refreshDB=args.refreshDB,
                  atteph=True,
                  verbose=args.verbose,
                  timeout=args.timeout)

    env(m)
    m.chk()
    if m.sensor is None and os.path.exists(args.filename):
        m.sensor = check_sensor(args.filename)
    if args.filename and m.finddb():
        m.setup()
    elif args.start and m.finddb():
        m.setup()
    else:
        m.setup()
        m.findweb()

    m.locate(forcedl=args.force)
    m.write_anc_par()
    m.cleanup()

    sys.exit(m.db_status)
