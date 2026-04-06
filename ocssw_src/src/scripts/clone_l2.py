#! /usr/bin/env python3


import argparse
import sys

import netCDF4

from seadasutils.netcdf_utils import nccopy_grp, nccopy_var

"""
Make a copy of an l2gen-style file, empty of geophysical data.

Copy the following from a L2 file:
- global attributes
- recursive copy of groups:
    sensor_band_parameters, scan_line_attributes, navigation_data, processing_control
- group geophysical_data attributes
- values geophysical_data/l2_flags

That is, everything except geophysical products.
"""

__version__ = "1.0.0_2026-02-10"
BAD_FIL = 101
BAD_GRP = 102
BAD_VAR = 103


def clone_l2(
    ifile,
    ofile,
    verbose=False,
):

    status = 0

    if verbose:
        print(f"Cloning {ifile} to {ofile}...")

    # open files
    try:
        src = netCDF4.Dataset(ifile, "r")
        src.set_auto_mask(False)
    except FileNotFoundError:
        print(f"File {ifile} not found")
        return BAD_FIL
    try:
        dst = netCDF4.Dataset(ofile, "w")
        dst.set_auto_mask(False)
    except:
        print(f"Error writing file {ofile}")
        return BAD_FIL

    # copy global attributes and dimensions
    nccopy_grp(src, dst, verbose=verbose, copygrps=False)

    # recursively copy specified groups
    to_copy = [
        "sensor_band_parameters",
        "scan_line_attributes",
        "navigation_data",
        "processing_control",
    ]
    for grpname in to_copy:
        if grpname not in src.groups:
            print(f"Group {grpname} not found")
            return BAD_GRP
        dstgrp = dst.createGroup(grpname)
        nccopy_grp(src.groups[grpname], dstgrp, verbose=verbose, copygrps=True)

    # copy empty geophysical_data group
    grpname = "geophysical_data"
    if grpname not in src.groups:
        print(f"Group {grpname} not found")
        return BAD_GRP

    # copy empty geophysical_data group
    dstgrp = dst.createGroup(grpname)
    nccopy_grp(src.groups[grpname], dstgrp, verbose=verbose, copyvars=False)

    # copy l2_flags
    varname = "l2_flags"
    if varname in src.groups[grpname].variables:
        srcvar = src.groups[grpname].variables[varname]
        nccopy_var(srcvar, dstgrp, verbose=verbose)
    else:
        print(f"Variable {grpname}/{varname} not found")
        return BAD_VAR

    return status


def main():
    print("clone_l2", __version__)

    # parse command line
    parser = argparse.ArgumentParser(
        description="Make a copy of an l2gen product file, empty of all geophysical data except l2_flags",
        epilog="""
EXIT Status:
0   : Success
1   : Fatal error
101 : File open error
102 : Group not found
103 : Variable not found
""",
    )
    parser.add_argument(
        "-v", "--verbose", help="print status messages", action="store_true"
    )
    parser.add_argument("ifile", help="Level 2 input file, from l2gen")
    parser.add_argument("ofile", help='"Empty" Level 2 output file')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    # run
    status = clone_l2(
        ifile=args.ifile,
        ofile=args.ofile,
        verbose=args.verbose,
    )
    return status


if __name__ == "__main__":
    sys.exit(main())
