#! /usr/bin/env python3

"""
Wrapper program for running the l1bgen program on MODIS L1A files.
"""
import argparse
from shlex import quote
import sys
from modis.modis_utils import buildpcf, modis_env
import modis.modis_L1B_utils as modisL1B
from seadasutils.setupenv import env


def main():
    """
    This is the primary driver function for the modis_L1B.py program.
    """
    version = "1.1"

    # Read commandline options...
    parser = argparse.ArgumentParser(prog="modis_L1B")
    parser.add_argument('--version', action='version', version='%(prog)s ' + version)
    parser.add_argument("filename", nargs='?',
                      help="Input L1A file", metavar="L1AFILE")  
    parser.add_argument("geofile", nargs='?',
                      help="INPUT GEOFILE filename - defaults to basename of L1AFILE +'.GEO'", metavar="GEOFILE")
    parser.add_argument("-p", "--parfile",
                      help="Parameter file containing program inputs", metavar="PARFILE")
    parser.add_argument("-o", "--okm",
                      help="Output L1B 1KM filename - defaults to '(A|T)YYYYDDDHHMMSS.L1B_LAC'", metavar="1KMFILE")

    parser.add_argument("-k", "--hkm",
        help="Output MODIS L1B HKM HDF filename", metavar="HKMFILE")
    parser.add_argument("-q", "--qkm",
        help="Output MODIS L1B QKM HDF filename", metavar="QKMFILE")
    parser.add_argument("-c", "--obc",
        help="Output MODIS L1B OBC HDF filename", metavar="OBCFILE")

    parser.add_argument("-l", "--lutver",
        help="L1B LUT version number", metavar="LUTVER")
    parser.add_argument("-d", "--lutdir",
        help="Path of directory containing LUT files", metavar="LUTDIR")

    parser.add_argument("-x", "--del-okm", action="store_const", const=1,
        default=0, help="Delete 1km  resolution L1B file")
    parser.add_argument("-y", "--del-hkm",  action="store_const", const=2,
        default=0, help="Delete 500m resolution L1B file")
    parser.add_argument("-z", "--del-qkm",  action="store_const", const=4,
        default=0, help="Delete 250m resolution L1B file")
    parser.add_argument("--keep-obc",  action="store_const", const=0,
        default=8, help="Save onboard calibration file")

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
    okm = None
    if args.okm:
        okm = quote(args.okm)
    hkm = None
    if args.hkm:
        hkm = quote(args.hkm)
    qkm = None
    if args.qkm:
        qkm = quote(args.qkm)
    obc = None
    if args.obc:
        obc = quote(args.obc)

    delfiles = args.del_okm + args.del_hkm + args.del_qkm + args.keep_obc

    l1b_instance = modisL1B.ModisL1B(inp_file=quote(args.filename),
                                   parfile=parfile,
                                   geofile=geofile,
                                   okm=okm,
                                   hkm=hkm,
                                   qkm=qkm,
                                   obc=obc,
                                   lutver=args.lutver,
                                   lutdir=args.lutdir,
                                   delfiles=delfiles,
                                   log=args.log,
                                   verbose=args.verbose
    )
    env(l1b_instance)
    modis_env(l1b_instance)
    l1b_instance.chk()
    buildpcf(l1b_instance)
    l1b_instance.run()
    return 0

if __name__ == "__main__":
    sys.exit(main())
