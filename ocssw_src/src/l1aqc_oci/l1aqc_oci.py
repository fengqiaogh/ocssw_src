#!/usr/bin/env python3

import sys
import argparse
from netCDF4 import Dataset

__version__ = '0.1.0_2020-12-04'

def get_ncdf_object(file,group,object):
    # read data object from NetCDF file   
     
    from netCDF4 import Dataset
    
    nc_fid = Dataset(file,'r')
    ngid = nc_fid.groups[group]
    data = ngid.variables[object][:].data
    nc_fid.close()
    
    return data 


def l1aqc_oci(args):
# procedure to get start, end times and more from Level-1a files for OCI

    try:
        nc_fid = Dataset(args.input_file,'r')
    except:
        print("Problem reading L1A file.")
        return 1
        
    status = 0
    
#  Generate output file if user desired it
    output = None
    if args.output: 
        if args.verbose:
            print("Writing output file: %s" % args.output)
        output = open(args.output,'w')
        output.write("# Info for %s\n" % args.input_file)
    
    try:
        stime = getattr(nc_fid, 'time_coverage_start')
        print("start_time=%sZ" % stime)
        if output:
            output.write("start_time=%sZ\n" % stime)
    except:
        return 120

    try:
        etime = getattr(nc_fid, 'time_coverage_end')
        print("stop_time=%sZ" % etime)

        if output:
            output.write("\nstop_time=%sZ\n" % etime)
            output.close()
    except:
        return 120
        
    try:
        nRec = getattr(nc_fid, 'number_of_filled_scans')
        print("number_of_records=%s" % nRec)

        if output:
            output.write("\number_of_records=%s\n" % nRec)
            output.close()
    except:
        return 120
    
    # quality flag will be an evolving development
    quality_flag = 0
    print("quality_flag=%d" % quality_flag)
    if output:
        output.write("\quality_flag=%d\n" % quality_flag)
        output.close()    
    
    nc_fid.close()
    return status

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,description=\
        'Reads OCI L1A science data and reports related info',epilog="""
EXIT Status:
0   : All is well in the world
1   : File not accessible
120 : Problem reading info from metadata
""")
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('input_file', type=str, help='path to the input L1a file')
    parser.add_argument('--output', type=str, help='path to the optional output text file')
    parser.add_argument('--verbose', '-v', action='store_true')

    args = parser.parse_args()

    status = l1aqc_oci(args)
    if status:
        sys.exit(status)

if __name__ == '__main__':
    sys.exit(main())
