#!/usr/bin/env python3

import sys
import os
import argparse
from l0info_utils import read_cfe_header, print_utils_version
from l0info_hkt import l0info_hkt
from l0info_oci import l0info_oci
from l0info_harp import l0info_harp
from l0info_spex import l0info_spex


__version__ = '1.7.3 (2023-08-10)'

__instids__ = ['hkt','oci','harp','spex']

def l0info_pace(args):

    try:
        fh = open(args.input_file, mode='rb')
    except:
        print("%s not found." % args.input_file)
        return 101
    
    rcode = 0
    
    print_utils_version()
    
    if args.verbose:
        print("Opened PACE data file: %s" % args.input_file)
        
    bDSB = False if args.noCFEheader else True
        
    # Get data type from CFE file header
    if bDSB:
        ctime, dtype = read_cfe_header(fh,args.verbose)

        if dtype==-1:
            fh.close()
            return 110
        
        # Check for EOF
        if fh.tell() == os.fstat(fh.fileno()).st_size:
            print('No packets in file.')
            fh.close()
            return 110
        
        fh.seek(0)
        try:
            strInst = __instids__[dtype]
        except:
            strInst = ''
    else:
        try:
            strInst = args.instrument.lower()
        except:
            strInst = ''
                
    output = None
    if args.output: 
        if args.verbose:
            print("Writing output file: %s" % args.output)
        output = open(args.output,'w')
        # output.write("# Info for %s\n" % args.input_file)
        
    print("datatype=%s\n"%strInst.upper())
    # Check if the CRC check requested
    crc_check = args.crc_check
    if strInst!='oci' and crc_check:
        print("Warning: CRC check is only availible for OCI Science data")
    # Call appropriate function depending on data type
    if strInst == 'hkt':    # PACE HKT data
        if output: output.write("datatype=HKT\n")
        rcode = l0info_hkt(args,fh,output)
    
    elif strInst == 'oci': # OCI science data
        # if output: output.write("datatype=OCI\n")
        rcode = l0info_oci(args,fh,output,bDSB,crc_check)

    elif strInst == 'harp': # HARP science data
        if output: output.write("datatype=HARP\n")
        rcode = l0info_harp(args,fh,output)

    elif strInst == 'spex': # SPEXone science data
        # if output: output.write("datatype=SPEX\n")
        rcode = l0info_spex(args,fh,output)
    
    else: 
        print('Invalid file or data type from CFE header. Or specify from [hkt,oci,harp,spex].')
        rcode = 102    

    fh.close()    
    if output:
        output.close()
        
    return rcode


def main():
    print("l0info_pace", __version__)
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,description=\
        'Reads OCI L0 HKT/science data and reports start/stop times',epilog="""
EXIT Status:
0   : All is well in the world
1   : Dunno, something horrible occurred
101 : File open error 
102 : Invalid file or instrument from CFE header
103 : Problem reading packet information
104 : Invalid packet [header/datatype/length]
110 : No [Science] packets found
111 : Bad image data
120 : Problem reading time from ancillary packet
131 : OCI Checksum failure
13X : OCI science data error
""")
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('input_file', type=str, help='path to the input L0 file')
    parser.add_argument('--output', type=str, help='path to the optional output text file')
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('--noCFEheader', action='store_true', help="data without CFE header")
    parser.add_argument('--instrument', type=str, help='hkt,oci,harp or spex')
    parser.add_argument('--crc_check', nargs='?', const=1, type=bool,default=False, help='validate crc checksum')
    args = parser.parse_args()        

    status = l0info_pace(args)
    if status:
        sys.exit(status)

if __name__ == '__main__':
    sys.exit(main())
