#!/usr/bin/env python3

import argparse
from osgeo import gdal
from affine import Affine
import sys
__version__ = '1.1'

# CONUS tile extents - Albers conic projection
# tile-index        (in meters)
#   h   v       ulx/llx     uly/lly
#   0   0       -2565585    3314805
#   32  21       2384415      14805
#
# for the 4x4 grid, the llx/lly are 2834415 and -285195
# ...unless I have x and y backward...
# the numbers passed to xy2rc work, so *meh*
parser = argparse.ArgumentParser(prog='cyan_conus_tiles',
        description='This script generates 54 GeoTIFF tiles in a 9x6 grid from an input full\nCONUS GeoTIFF\n' +
        'The output tile names are based on the input file name with\n'+
        '_row_col.tif appended')
parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
parser.add_argument('ifile',help='input CONUS GeoTIFF file to slice and dice')
parser.add_argument('--verbose', '-v', action='store_true')

args = parser.parse_args()

nXtiles = 9
nYtiles = 6

sourcefile = args.ifile
sourcefilebase = sourcefile.rstrip('.tif')
#Open source dataset
srcDS = gdal.Open(sourcefile)
#Find our starting point
T0 = Affine.from_gdal(*srcDS.GetGeoTransform())
xy2rc = lambda xm, ym: (ym, xm) * ~T0

xoff,yoff = xy2rc(3314805,-2565585)
xoff = int(xoff)
yoff = int(yoff)
endx,endy = xy2rc(-285195,2834415)
endx = int(endx)
endy = int(endy)

xstride = int((endx-xoff)/nXtiles)
ystride = int((endy-yoff)/nYtiles)

for y in range(0,nYtiles):
    ytile = y + 1
    uly = y * ystride + yoff

    for x in range(0,nXtiles):
        xtile = x +1
        ulx = x * xstride + xoff
        srcwin = [ulx,uly,xstride,ystride]
        tileFile = sourcefilebase + '_' + str(xtile) + '_' + str(ytile) + '.tif'
        if args.verbose:
            print('generating %s' % tileFile)
        transop = gdal.TranslateOptions(creationOptions=['COMPRESS=LZW'],srcWin=srcwin)
        gdal.Translate(tileFile,srcDS, options=transop)

sys.exit(0)
