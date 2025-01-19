#!/usr/bin/env python3
#
import argparse
import os
import re
import requests
import sys
import textwrap
from seadasutils.ProcUtils import httpdl, compare_checksum
from urllib.parse import urlparse
from pathlib import Path, PurePath
from datetime import datetime
import time

def retrieveURL(request,localpath='.', uncompress=False, verbose=0,force_download=False, appkey=False, checksum=False,timestamp=False):
    if args.verbose:
        print("Retrieving %s" % request.rstrip())

    server = "oceandata.sci.gsfc.nasa.gov"
    parsedRequest = urlparse(request)
    netpath = parsedRequest.path

    if parsedRequest.netloc:
        server = parsedRequest.netloc
    else:
        if not re.match(".*getfile",netpath):
            netpath = '/ob/getfile/' + netpath

    joiner = '?'
    if (re.match(".*getfile",netpath)) and appkey:
        netpath = netpath + joiner +'appkey=' + appkey
        joiner = '&'

    if parsedRequest.query:
        netpath = netpath + joiner + parsedRequest.query

    status = httpdl(server, netpath, localpath=localpath, uncompress=uncompress, verbose=verbose,force_download=force_download,timestamp=timestamp)
    
    if checksum and not uncompress:
        cksumURL = 'https://'+server + '/checkdata/' + parsedRequest.path
        dnldfile = localpath / parsedRequest.path
        if compare_checksum(dnldfile,requests.get(cksumURL).text):
            print("The file %s failed checksum test" % parsedRequest.path)
            status = 1

    return status


if __name__ == "__main__":
    # parse command line
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description='Download files archived at the OB.DAAC',
        epilog=textwrap.dedent('''
Provide one of either filename, --filelist or --http_manifest.

NOTE: For authentication, a valid .netrc file in the user home ($HOME) directory\nor a valid appkey is required.

    Example .netrc:
    machine urs.earthdata.nasa.gov login USERNAME password PASSWD\n

    An appkey can be obtained from:
    https://oceandata.sci.gsfc.nasa.gov/appkey/
'''
    ))
    parser.add_argument('-v', '--verbose', help='print status messages',
                        action='count',default=0)
    parser.add_argument('filename', nargs='?', help='name of the file (or the URL of the file) to retreive')
    parser.add_argument('--filelist',
                        help='file containing list of filenames to retreive, one per line')    
    parser.add_argument('--http_manifest',
                        help='URL to http_manifest file for OB.DAAC data order')
    parser.add_argument('--odir',
                        help='full path to desired output directory; \ndefaults to current working directory: %s' % Path.cwd(),
                        default=Path.cwd())
    parser.add_argument('--uncompress',action="store_true",
                        help='uncompress the retrieved files (if compressed)',
                        default=False)
    parser.add_argument('--checksum',action="store_true",
                        help='compare retrieved file checksum; cannot be used with --uncompress',
                        default=False)
    parser.add_argument('--timestamp',action="store_true",
                        help='change timestamp of retrieved file to time of creation',
                        default=False)
    parser.add_argument('--failed',help='filename to contain list of files that failed to be retrieved')
    parser.add_argument('--appkey',help='value of the users application key')
    parser.add_argument('--force',action='store_true',
                        help='force download even if file already exists locally',
                        default=False)
    args = parser.parse_args()

    filelist = []

    if args.http_manifest:
        status = retrieveURL(args.http_manifest,verbose=args.verbose,force_download=True,appkey=args.appkey)
        if status:
            print("There was a problem retrieving %s (received status %d)" % (args.http_manifest,status))
            sys.exit("Bailing out...")
        else:
            with open('http_manifest.txt') as flist:
                for filename in flist:
                    filelist.append(filename.rstrip())
    elif args.filename:
        filelist.append(args.filename)
    elif args.filelist:
        with open(os.path.expandvars(args.filelist)) as flist:
            for filename in flist:
                filelist.append(os.path.expandvars(filename.rstrip()))

    if not len(filelist):
        parser.print_usage()
        sys.exit("Please provide a filename (or list file) to retrieve")

    if args.uncompress and args.checksum:
        parser.print_usage()
        sys.exit("--uncompress is incompatible with --checksum")

    outpath = Path.resolve(Path.expanduser(Path(os.path.expandvars(args.odir))))

    if args.verbose:
        print("Output directory: %s" % outpath)

    failed = None
    if args.failed:
        failed = open(args.failed, 'w')

    for request in filelist:
        status = retrieveURL(request,localpath=outpath, uncompress=args.uncompress,
                             verbose=args.verbose,force_download=args.force,
                             appkey=args.appkey,checksum=args.checksum,timestamp=args.timestamp)
        if status:
            if status == 304: 
                if args.verbose:
                    print("%s is not newer than local copy, skipping download" % request)
            else:
                print("There was a problem retrieving %s (received status %d)" % (request,status))
                if failed:
                    failed.write(request)
                    failed.write("\n")

    if failed:
        failed.close()
