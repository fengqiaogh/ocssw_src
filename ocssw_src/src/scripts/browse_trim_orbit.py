#! /usr/bin/env python3

import argparse
import numpy as np
import netCDF4 as nc
import sys


def find_lines(ifile, verbose=False):
    status = 0

    try:
        l2file = nc.Dataset(ifile,'r')
    except:
        status = 1
        if verbose:
            print(f'Failed to open {ifile}')
        return status
    
    lat = (l2file['navigation_data']['latitude'][:]).data
    lon = (l2file['navigation_data']['longitude'][:]).data

    # Get the gradient for the latitude array, compute the median of the y axis (lines)
    grad_y, grad_x = np.gradient(lat)
    grad = np.median(grad_y, axis=1)

    # Get the sign of each element
    # Returns 1 for positive, -1 for negative, and 0 for zero
    signs = np.sign(grad)
    # Compare adjacent signs using slicing
    # This creates a boolean array where True indicates a change in sign
    sign_changes = signs[:-1] != signs[1:]

    # Use np.where() to get the indices where the condition is True
    # We add 1 to the indices because the comparison checks the *transition* 
    # between index i and i+1, so the change occurs *at* index i+1.
    indices_of_change = np.where(sign_changes)[0] + 1

    sline=1
    eline=len(grad)
    if len(indices_of_change):
        if len(indices_of_change) == 1:
            if indices_of_change[0] < eline / 2:
                sline = indices_of_change[0]
            else:
                eline = indices_of_change[0]
        elif len(indices_of_change) == 2:
                sline = indices_of_change[0]
                eline = indices_of_change[1]
        else:
            status = 1
    
    if not status:
        cpixl = int((lon.shape)[1]/2)
        cline = int((eline - sline) / 2) + sline
        lon0 = lon[cline,cpixl]
        lat0 = lat[cline,cpixl]
        if (sline != 1) or (eline != len(grad)):
            print(f'sline={sline}\neline={eline}\n#\n')
        else:
            status = 110
        print(f'central_meridian={lon0:.4f}\nlat_0={lat0:.4f}')

    return status

def main():

    # parse command line
    parser = argparse.ArgumentParser(
        description="Find start and end lines for an L2 file that keep the file from crossing the poles",
        epilog="""
EXIT Status:
0   : Success
1   : Fatal error - more than two transitions found
110 : No extract needed
""",
    )
    parser.add_argument("ifile", help="Level 2 input file")
    parser.add_argument("-v", "--verbose", help="print status messages", action="store_true")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    # run
    status = find_lines(args.ifile, verbose=args.verbose)
    return status


if __name__ == "__main__":
    sys.exit(main())
