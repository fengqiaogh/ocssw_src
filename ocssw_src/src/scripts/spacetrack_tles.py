#! /usr/bin/env python3
# Retrieve TLE files from the  https://www.space-track.org website

import argparse
import os
import shutil
import logging
import requests
import sys
import json

__version__ = "2.0-20200212"

norad = {
    'TERRA':       25994,
    'AQUA':        27424,
    'CALIPSO':     29108,
    'CLOUDSAT':    29107,
    'PROBA':       26958,
    'NPP':         37849,
    'JPSS1':       43013,
    'LANDSAT-8':   39084,
    'SEAHAWK-1':   43820,
    'SENTINEL-2A': 40697,
    'SENTINEL-2B': 42063,
    'SENTINEL-3A': 41335,
    'SENTINEL-3B': 43437,
    'OCEANSAT-2':  35931,
    'GCOM-C':      43065
}


parser = argparse.ArgumentParser(prog='spacetrack_tles',
        description='This script retrieves Two-Line Elements from www.space-track.org for a given satellite name or NORAD catalog ID')
parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
parser.add_argument('satellite',nargs='?', help='Satellite TLE to retrieve.  "ALL" will return TLEs for all in the list above.', type=str,
                        choices=['TERRA',
                                 'AQUA',
                                 'CALIPSO',
                                 'CLOUDSAT',
                                 'PROBA',
                                 'NPP',
                                 'JPSS1',
                                 'LANDSAT-8',
                                 'SEAHAWK-1',
                                 'SENTINEL-2A',
                                 'SENTINEL-2B',
                                 'SENTINEL-2B',
                                 'SENTINEL-3A',
                                 'SENTINEL-3B',
                                 'OCEANSAT-2',
                                 'GCOM-C',
                                 'ALL'])
parser.add_argument('--ofile','-o', type=str, default=None, help="Name of file to store retreived TLE")
parser.add_argument('--catID', '-c', type=int, default=None, help='NORAD Catalog ID of satellite to retrieve; ingored if a satellite name is provided.')
parser.add_argument('--username', '-u', type=str, default=None, help="Space-track.org username")
parser.add_argument('--password', '-p', type=str, default=None, help="Space-track.org password")
parser.add_argument('--copyTCE', default=False, action='store_true', help="Copy TLE to latest '<satellite>.tce' file")
parser.add_argument('--directory','-d', type=str, default=os.getcwd(), help="Direcotry to write TLEs; default is the current working directory")
parser.add_argument('--norad_catIDs', '-n', type=str, default=os.path.join(os.environ['OCDATAROOT'],'common','norad_catalog.json'), help="JSON file in the form {'satellite':<catID>}")
parser.add_argument('--verbose', '-v', action='count', default=0)

args = parser.parse_args()

if args.verbose > 1:
    logging.basicConfig(level=logging.DEBUG)

identity = args.username
password = args.password

if not identity:
    netrc = os.path.join(os.environ['HOME'],'.netrc')
    with open(netrc,'r') as cred:
        for line in cred:
            if line.startswith('machine www.space-track.org'):
                parts = line.split()
                identity = parts[3]
                password = parts[5]

credentials={'identity':identity,'password':password}

if not identity or not password:
    sys.exit("Authentication credentials required! Check your $HOME/.netrc or pass in via --username/--password")

baseTLEdir = args.directory

if args.norad_catIDs:
    if os.path.isfile(args.norad_catIDs):
        try:
            with open(args.norad_catIDs) as fj:
                norad = json.load(fj)
        except json.JSONDecodeError as e:
            print(e)
            sys.exit("Failed to read %s" % args.norad_catIDs)
    elif args.norad_catIDs == os.path.join(os.environ['OCDATAROOT'],'common','norad_catalog.json') and not os.path.isfile(args.norad_catIDs):
        if args.verbose:
            print("%s does not exist, using default set" % args.norad_catIDs)
    else:
        sys.exit("%s does not exist...exiting" % args.norad_catIDs)

tleIDs = None
if args.satellite:
    if args.satellite == 'ALL':
        tleIDs = ','.join(str(tle) for tle in norad.values())
    else:
        tleIDs = str(norad[args.satellite])

if args.catID:
        tleIDs = str(args.catID)

if not tleIDs:
    parser.print_help()
    sys.exit("\nEither provide a satellite name (or ALL) or as --catID")

baseurl = "https://www.space-track.org"
loginurl = '/'.join([baseurl, "ajaxauth","login"])
tleurl = '/'.join([baseurl,"basicspacedata","query","class","tle_latest","ORDINAL","1","NORAD_CAT_ID",tleIDs,"format","tle"])

#Get login cookie
stConn = requests.Session()

#Get TLEs
tlehash = {}

try:
    req = stConn.post(loginurl, data=credentials)
    if req.status_code != 200:
        sys.exit("Received response: %d from space-track.org...unable to log in...exiting..." % req.status_code)

    if args.verbose:
        print("Contacting Space-Track with the following request\n\t%s" % tleurl)

except requests.exceptions.RequestException as e:
    print(e)
    sys.exit(1)

try:
    tlereq = stConn.get(tleurl)
    if tlereq.status_code != 200:
        sys.exit("Received response: %d from space-track.org...exiting..." % tlereq.status_code)

    for line in tlereq.text.splitlines():
        if len(line) < 10:
            continue
        parts = line.split()
        tleID = int(parts[1].strip('U'))
        if tleID in tlehash:
            tlehash[tleID].append(line)
        else:
            tlehash[tleID] = [line]

    tlereq.close()
except requests.exceptions.RequestException as e:
    print(e)
    sys.exit(1)

#Process retrieved TLEs
for tleID in tlehash.keys():

    bird = 'NORAD_' + str(tleID)
    for sat, cat in norad.items():
        if cat == tleID:
            bird = sat

    tle = tlehash[tleID]
    if len(tle) != 2:
        print ("Whoops! TLE retrieval error for %s" % bird)
    else:
        tleDate = tle[0].split()[3]
        year = str(int(tleDate[0:2]) + 2000)
        doy = tleDate[2:5]
        tleFile = os.path.join(baseTLEdir,bird+doy+year+".dat")
        if args.ofile:
            tleFile = args.ofile

        if args.verbose:
            print("Writing TLE for NORAD catalog ID %d to %s" % (tleID,tleFile))

        with open(tleFile,'a') as f:
            f.write(tle[0])
            f.write("\n")
            f.write(tle[1])
            f.write("\n")

        if args.copyTCE:
            tceFile = os.path.join(baseTLEdir,bird+".tce")
            if args.ofile:
                tceFile, ext = os.path.splitext(tleFile)
                tceFile = tceFile + ".tce"

            if args.verbose:
                print("Copying TLE: %s to %s" % (tleFile,tceFile))

            shutil.copyfile(tleFile,tceFile)
