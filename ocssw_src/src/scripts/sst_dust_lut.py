#!/usr/bin/env python3

# This script converts the ASCII SST coefficient tables provided by RSMAS
# into a netCDF4 file.  Since this is a very infrequent occurance, you may
# need to modify this as appropritate when new tables are provided.

import pandas as pd
from netCDF4 import Dataset
import datetime
import argparse
import sys
import os

__version__ = '1.0.0_2019-05-18'

def main():
    """
    Primary driver of the program; get command line arguments, check the files
    specified and kick off the processing
    """

    parser = argparse.ArgumentParser(description=\
        'Converts ASCII SST dust correction coefficient file to netCDF.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('input_file', type=str, help='path to the input ASCII LUT')
    parser.add_argument('output_file', type=str, help='path to the output netCDF LUT')
    parser.add_argument('--lut_version', type=str, help='LUT version',default='1.0')
    parser.add_argument('--fields', nargs='+', help="field names for ASCII table; default:['c0','c1','c2','c3','c4','c5','c6','c7','c8','c9','c10','c11','c12','dsdi_thresh']")
    parser.add_argument('--mission', default="modis-aqua", help="mission for which the LUTs apply (modis-aqua, modis-terra, viirs-npp")
    parser.add_argument('--verbose', '-v', action='store_true')

    args = parser.parse_args()

    numberOfCoefficients = 13
    inputFilename = args.input_file
    outputFilename = args.output_file
    fields = ['c0','c1','c2','c3','c4','c5','c6','c7','c8','c9','c10','c11','c12','dsdi_thresh','dust_thresh']
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

#    data = pd.read_csv(inputFilename, skiprows=skiprows, delim_whitespace=True, names=fields)
    data = pd.read_csv(inputFilename, skiprows=skiprows, names=fields)

    #cols = ['month', 'lat1', 'lat2', 'a0', 'a1', 'a2', 'a3', 'a6', 'a5', 'a4', 'count']

    coefficientList = ['c0','c1','c2','c3','c4','c5','c6','c7','c8','c9','c10','c11','c12' ]
    if numberOfCoefficients != 13:
        coefficientList = []
        for coef in range (0,numberOfCoefficients):
            coefficientList.append('c'+str(coef))
    coefficients = data[coefficientList]

    print(data)
    #print(bounds[bounds.index.isin([0, 1,2,3, 4,5,6])])
    #print(coefficients)

    rootgrp = Dataset(outputFilename, "w", format="NETCDF4")
    rootgrp.createDimension("coefficient", numberOfCoefficients)

    instrument = 'MODIS'
    platform = 'Aqua'

    if args.mission == 'modis-terra':
        platform = 'Terra'

    if args.mission == 'viirs-npp':
        instrument = 'VIIRS'
        platform = 'SNPP'

    rootgrp.title = "%s-%s SST dust correction coefficients" % (instrument,platform)
    rootgrp.instrument = instrument
    rootgrp.platform = platform
    rootgrp.description = "SST dust correction coefficients for %s-%s" % (instrument,platform)
    if args.lut_version:
        rootgrp.version = args.lut_version
        if args.verbose:
            print("LUT version: %s" %  args.lut_version)

    rootgrp.date_created = create_date
    rootgrp.product_name = outputFilename
    rootgrp.dsdi_threshold = data['dsdi_thresh']
    rootgrp.dust_threshold = data['dust_thresh']
    history = ' '.join(sys.argv)
    history = history.replace(os.environ["OCSSWROOT"],"$OCSSWROOT")
    rootgrp.history = history.encode('ascii')


    coef = rootgrp.createVariable("coefficients","f4",("coefficient",),fill_value=-32767.0)
    coef.long_name = "SST dust correction coefficients"
#    coef.equation = "SST = C0 + C1*BT11 + C2*(BT11-BT12)*bsst + C3*sec(abs(satz)) + C4*satz + C5*satz^2 + C6*mirrorside"
    coef.description = "SST dust correction coefficients"

    coef[:] = coefficients.values

    rootgrp.close()
    if args.verbose:
        print("Done!")

# The following allows the file to be imported without immediately executing.
if __name__ == '__main__':
    sys.exit(main())
