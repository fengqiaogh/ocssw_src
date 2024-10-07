#! /usr/bin/env python3

"""
Program to check for updated ancillary data files and download them
as appropriate.
"""

import argparse
import sys
import seadasutils.anc_utils as ga
from seadasutils.setupenv import env

def main():
    """
    The main function for the program. Gets and checks command line options, instantiates 
    a getanc object (from the class defined in anc_utils) and then calls the methods to 
    get the ancillary data.
    """

    version = "2.2"

    # Read commandline options...
    usage = """
    %(prog)s [OPTIONS] FILE
          or
    -s,--start YYYYDDDHHMMSS [-e,--end YYYDDDHHMMSS]  [OPTIONS]

      FILE  Input L1A or L1B file

    NOTE: Currently NO2 climatological data is used for OBPG operational
          processing, so to match OBPG distributed data products, the default
          behaviour disables NO2 searching.

    This program queries an OBPG server and optionally downloads the optimal
    ancillary data files for Level-1 to Level-2 processing. If an input file
    is specified the start and end times are determined automatically, otherwise
    a start time must be provided by the user.

    A text file (with the extension '.anc') is created containing parameters
    that can be directly used as input to the l2gen program for optimal Level-1
    to Level-2 processing, e.g.:

         l2gen ifile=<infile> ofile=<outfile> par=<the *.anc text file>

    EXIT STATUS:
        0  : all optimal ancillary files exist and are present on the locally
        99 : an error was encountered; no .anc parameter text file was created
        31 : no ancillary files currently exist corresponding to the start
             time and therefore no .anc parameter text file was created
      1-30 : bitwise value indicating one or more files are not optimal:

             bit 0 set = missing one or more MET files
             bit 1 set = missing one or more OZONE files
             bit 2 set = no SST file found
             bit 3 set = no NO2 file found
             bit 4 set = no ICE file found

    e.g. STATUS=11 indicates there are missing optimal MET, OZONE, and NO2 files

    """

    parser = argparse.ArgumentParser(prog="getanc",usage=usage)
    parser.add_argument('--version', action='version', version='%(prog)s ' + version)
    parser.add_argument("filename", nargs='?',
                      help="Input L1 file", metavar="L1FILE")  
    parser.add_argument("-s", "--start",
                      help="Time of the first scanline (if used, no input file is required)",
                      metavar="START")
    parser.add_argument("-e", "--stop",
                      help="Time of last scanline", metavar="STOP")

    ancdb_help_text = "Use a custom filename for ancillary database. If " \
                      "full path not given, ANCDB is assumed to exist "\
                      "(or will be created) under " + \
                      ga.DEFAULT_ANC_DIR_TEXT + "/log/. If " + \
                      ga.DEFAULT_ANC_DIR_TEXT + "/log/ does not " \
                      "exist, ANCDB is assumed (or will be created) " \
                      "under the current working directory"
    parser.add_argument("--ancdb", default='ancillary_data.db',help=ancdb_help_text, metavar="ANCDB")
    parser.add_argument("-o", "--ofile", help="output ancillary par file", metavar="ANC_FILE")
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
    parser.add_argument("-m", "--mission", help="Mission name", metavar="MISSION")

    parser.add_argument("-i", "--ice", action="store_false", default=True,
                      help="Do not search for sea-ice ancillary data")
    parser.add_argument("-n", "--no2", action="store_true", default=False, 
                      help="Search for NO2 ancillary data")
    parser.add_argument("-t", "--sst", action="store_false", default=True,
                      help="Do not search for SST ancillary data")
    parser.add_argument("-v", "--verbose", action="count",
                      default=0, help="print status messages")
    parser.add_argument("--noprint", action="store_false", default=True,
                      help="Suppress printing the resulting list of files to the screen")
    parser.add_argument("-f", "--force-download", action="store_true", dest='force', default=False,
                      help="Force download of ancillary files, even if found on hard disk")
    parser.add_argument("-u", "--use_filename", action="store_true", default=False, 
                      help="Use filename to call API instead of deriving start time")

    args = parser.parse_args()
    if args.filename is None and args.start is None:
        parser.print_help()
        sys.exit(32)

    g = ga.getanc(filename=args.filename,
                  start=args.start,
                  stop=args.stop,
                  ancdir=args.ancdir,
                  ancdb=args.ancdb,
                  anc_file=args.ofile,
                  curdir=args.curdir,
                  sensor=args.mission,
                  opt_flag=5,  # defaults to retrieving met, ozone, sst, and ice data
                  verbose=args.verbose,
                  printlist=args.noprint,
                  download=args.download,
                  timeout=args.timeout,
                  refreshDB=args.refreshDB,
                  use_filename=args.use_filename)

    if args.sst is False:
        g.set_opt_flag('sst', off=True)
    if args.no2:
        g.set_opt_flag('no2')
    if args.ice is False:
        g.set_opt_flag('ice', off=True)

    env(g)
    g.chk()
    if (args.filename and g.finddb()) or (args.start and g.finddb()):
        g.setup()
    else:
        g.setup()
        g.findweb()
    g.locate(forcedl=args.force)
    g.write_anc_par()
    g.cleanup()
    return(g.db_status)

if __name__ == "__main__":
    exit(main())
