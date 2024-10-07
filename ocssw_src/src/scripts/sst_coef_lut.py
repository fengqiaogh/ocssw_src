#!/usr/bin/env python3

# This script converts the ASCII SST coefficient tables provided by RSMAS
# into a netCDF4 file.  Since this is a very infrequent occurance, you may
# need to modify this as appropritate when new tables are provided.

import pandas as pd
from netCDF4 import Dataset
import datetime
import argparse
import sys

__version__ = '1.0.0_2019-04-05'

def main():
    """
    Primary driver of the program; get command line arguments, check the files
    specified and kick off the processing
    """

    parser = argparse.ArgumentParser(description=\
        'Converts ASCII "Latband" SST coefficient LUTs to netCDF.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('--nlatbands', default=7, type=int, help='number of latbands in LUT')
    parser.add_argument('--ncoefficients', default=7, type=int, help='number of coefficients in LUT')
    parser.add_argument('input_file', type=str, help='path to the input ASCII LUT')
    parser.add_argument('output_file', type=str, help='path to the output netCDF LUT')
    parser.add_argument('--lut_version', type=str, help='version of the LUT being converted')
    parser.add_argument('--lut_description', default="11 um", type=str,help='string to identify sst type to apply these coefficients (11 um, 4 um, Triple Window)')
    parser.add_argument('--fields', nargs='+', help="field names for ASCII table; default:['sensor','month','minlat','maxlat','a0','a1','a2','a3','a4','a5','a6','count']")
    parser.add_argument('--mission', default="modis-aqua", help="mission for which the LUTs apply (modis-aqua, modis-terra, viirs-npp")
    parser.add_argument('--verbose', '-v', action='store_true')

    args = parser.parse_args()

    numberOfLatbands = args.nlatbands
    numberOfCoefficients = args.ncoefficients
    inputFilename = args.input_file
    outputFilename = args.output_file
    lutdescription = args.lut_description
    fields = ['sensor','month','minlat','maxlat','a0','a1','a2','a3','a4','a5','a6','count']
    if args.fields:
        fields = args.fields

    create_date = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')

    skiprows = 0
    with open(inputFilename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                skiprows = skiprows + 1

    if args.verbose:
        print("Processing %s ..." % inputFilename)
        print("...found %d header lines to skip" % skiprows)
        print("Creating to %s ..." % outputFilename)
        print("with coefficients for %s" % lutdescription)

    data = pd.read_csv(inputFilename, skiprows=skiprows, delim_whitespace=True, names=fields)

    #cols = ['month', 'lat1', 'lat2', 'a0', 'a1', 'a2', 'a3', 'a6', 'a5', 'a4', 'count']
    data['latband'] = data.index % numberOfLatbands

    bounds = data[['minlat','maxlat'] ]
    coefficientList = ['a0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6' ]
    if numberOfCoefficients != 7:
        coefficientList = []
        for coef in range (0,numberOfCoefficients):
            coefficientList.append('a'+str(coef))
    coefficients = data[coefficientList]

    #print(bounds[bounds.index.isin([0, 1,2,3, 4,5,6])])
    #print(coefficients)

    rootgrp = Dataset(outputFilename, "w", format="NETCDF4")
    rootgrp.createDimension("month", 12)
    rootgrp.createDimension("latband", numberOfLatbands)
    rootgrp.createDimension("coefficient", numberOfCoefficients)
    rootgrp.createDimension("latitude",2)

    instrument = 'MODIS'
    platform = 'Aqua'

    if args.mission == 'modis-terra':
        platform = 'Terra'

    if args.mission == 'viirs-npp':
        instrument = 'VIIRS'
        platform = 'SNPP'

    rootgrp.title = "%s-%s SST coefficients" % (instrument,platform)
    rootgrp.instrument = instrument
    rootgrp.platform = platform
    rootgrp.description = "Latitude band derived %s SST coefficients for %s-%s" % (lutdescription,instrument,platform)
    if args.lut_version:
        rootgrp.version = args.lut_version
        if args.verbose:
            print("LUT version: %s" %  args.lut_version)

    rootgrp.date_created = create_date
    rootgrp.product_name = outputFilename

    coef = rootgrp.createVariable("coefficients","f4",("month","latband","coefficient",),fill_value=-32767.0)
    coef.long_name = "SST (%s) coefficients" % lutdescription
    coef.equation = "SST = C0 + C1*BT11 + C2*(BT11-BT12)*bsst + C3*sec(abs(satz)) + C4*satz + C5*satz^2 + C6*mirrorside"
    coef.description = "SST (%s) coefficients based on latitude band and month of the year" % lutdescription
    coef.comment = 'Month Index: 0=Jan, 1=Feb, 2=Mar, 3=Apr, 4=May, 5=Jun, 6=Jul, 7=Aug, 8=Sep, 9=Oct, 10=Nov, 11=Dec'

    bound = rootgrp.createVariable("latbands","f4",("latband","latitude",),fill_value=-999.)
    bound.valid_min = -90.
    bound.valid_max = 90.
    bound.long_name = 'Latitude Bands for SST coefficients'
    bounds.units = 'degrees_north'

    latbandList = []
    for latband in range (0,numberOfLatbands):
        latbandList.append(latband)
    b = bounds[bounds.index.isin(latbandList)].values
    bound[:] = b

    coef[:] = coefficients.values

    rootgrp.close()
    if args.verbose:
        print("Done!")

# The following allows the file to be imported without immediately executing.
if __name__ == '__main__':
    sys.exit(main())
