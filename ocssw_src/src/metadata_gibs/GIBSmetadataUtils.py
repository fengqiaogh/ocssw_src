#!/usr/bin/env python3
#import hashlib
import datetime
from os import path
import sys
import argparse
from netCDF4 import Dataset

__version__ = '1.1.0_2019-10-28'

def gring_to_geopolygon(lats, lons):
    merged = [None]*(len(lons)+len(lats))
    merged[::2] = lons
    merged[1::2] = lats

# yeah, not very pythonic, but it works mate ;)
    gPolygon = ''
    i=-1
    for element in merged:            
        delimiter = ' '
        if i%2:
            delimiter = ','
        if i==-1:
            delimiter = ''
        gPolygon = gPolygon + delimiter + str(element)
        i = i+1
    gPolygon = gPolygon + ',' + str(merged[0]) + ' ' + str(merged[1])
    return gPolygon

#############################################################################

def create_granule_metadata(browseFileName, fileBaseName, starttime, stoptime, datadays):
    '''This method writes the xml metadata for the GIBS imagery'''

    productionTime = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')
    xmlbody = '''<ImageryMetadata
    xmlns="http://www.w3.org/2005/Atom"
    xmlns:georss="http://www.georss.org/georss/10">
    <ProviderProductId>%s</ProviderProductId>
    <ProductionDateTime>%s</ProductionDateTime>
    <DataStartDateTime>%s</DataStartDateTime>
    <DataEndDateTime>%s</DataEndDateTime>
    <DataDay>%s</DataDay>
</ImageryMetadata>
    ''' % (path.basename(browseFileName),
           productionTime,
           starttime,
           stoptime,datadays[0].strftime("%Y%j"))

    return xmlbody

#############################################################################

#def generateSHA1(filepath):
#    BLOCKSIZE = 65536
#    sha1 = hashlib.sha1()
#    with open(filepath, 'rb') as afile:
#        buf = afile.read(BLOCKSIZE)
#        while len(buf) > 0:
#            sha1.update(buf)
#            buf = afile.read(BLOCKSIZE)
#    return sha1.hexdigest()
#    sha1sum "$1"|cut -d' ' -f1
#    cmd = 'sha1sum %s' % filepath
#    sha1 = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
#    return sha1.split(' ')[0]

#############################################################################

def main():
    """
    Write XML metadata file for GIBS browse imagery
    """
    parser = argparse.ArgumentParser(description=\
        'Generate XML metadata file for GIBS browse imagery')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('-i','--input_file', type=str,  required=True)
    parser.add_argument('-t','--tiff_file', type=str,  required=True)
    parser.add_argument('-o','--output_file', type=str)
    parser.add_argument('-d','--datadays', nargs='+', type=str, help="dataday(s) in YYYDDD or YYYY-MM-DD format")
    parser.add_argument('-v','--verbose', action='store_true')

    args = parser.parse_args()
         
    rootgrp = Dataset(args.input_file, "r", format="NETCDF4")
    starttime = rootgrp.getncattr("time_coverage_start")
    stoptime = rootgrp.getncattr("time_coverage_end")
#    navgrp = rootgrp.createGroup("navigation_data")
#    gRingLats = navgrp.getncattr("gringpointlatitude"); 
#    gRingLons = navgrp.getncattr("gringpointlongitude");
#    gPolygon = gring_to_geopolygon(gRingLats,gRingLons)

    datadays = []
    if args.datadays:
        for dataday in args.datadays:
            if len(dataday) == 7:
                datadays.append(datetime.datetime.strptime(dataday, "%Y%j"))
            elif len(dataday) == 10:
                datadays.append(datetime.datetime.strptime(dataday, "%Y-%m-%d"))
            else:
                print("Could not parse data day: %s" % dataday)
                sys.exit(1)
    else:
        datadays.append(datetime.datetime.strptime(starttime[:10], "%Y-%m-%d"))

    outputFile = args.tiff_file + '.xml'
    if args.output_file:
        outputFile = args.output_file
        
    if args.verbose:
        print("StartTime: %s" % starttime)
        print("StopTime: %s" % stoptime)
        print("Input file: %s" % args.input_file)
        print("TIFF file: %s" % args.tiff_file)
        print("Output file: %s" % outputFile)
        for cnt,dataday in enumerate(datadays):
            print("DataDay[%d]: %s" % (cnt,dataday.strftime("%Y%j")))
#        print("gRingLats: %s" % gRingLats)
#        print("gRingLons: %s" % gRingLons)
#        print("gPolygon: %s" % gPolygon)

    xmlbody = create_granule_metadata(args.tiff_file, args.input_file, starttime, stoptime, datadays)

    xmlfile = open(outputFile,'w')
    xmlfile.write(xmlbody)
    xmlfile.close()

    if args.verbose:
        print("XML for %s:" % outputFile)
        print(xmlbody)

# The following allows the file to be imported without immediately executing.
if __name__ == '__main__':
    sys.exit(main())

