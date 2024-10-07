#!/usr/bin/env python3

import argparse
from osgeo import gdal
from affine import Affine
import sys

__version__ = '1.1'
# CONUS tile extents - Albers conic projection
# tile-index        (in meters)
#   h   v       ulx/lrx     uly/lry
#   0   0       -851715	     2474325
#   16  13      1698285      374325
#
# for the 4x4 grid, the lrx/lry are 1683285 and 404325
# ...unless I have x and y backward...
# the numbers passed to xy2rc work, so *meh*
parser = argparse.ArgumentParser(prog='cyan_conus_tiles',
        description='This script generates 12 GeoTIFF tiles in a 4x3 grid from an input full\nAlaskan GeoTIFF\n' +
        'The output tile names are based on the input file name with\n'+
        '_row_col.tif appended')
parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
parser.add_argument('ifile',help='input Alaskan GeoTIFF file to slice and dice')
parser.add_argument('--verbose', '-v', action='store_true')

args = parser.parse_args()

nXtiles = 4
nYtiles = 3

sourcefile = args.ifile
sourcefilebase = sourcefile.rstrip('.tif')
#Open source dataset
srcDS = gdal.Open(sourcefile)
#Find our starting point
T0 = Affine.from_gdal(*srcDS.GetGeoTransform())
xy2rc = lambda xm, ym: (ym, xm) * ~T0

xoff,yoff = xy2rc(2474325,-851715)
xoff = int(xoff)
yoff = int(yoff)
endx,endy = xy2rc(404325,1683285)
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
            print("generating %s" % tileFile)
        transop = gdal.TranslateOptions(creationOptions=['COMPRESS=LZW'],srcWin=srcwin)
        gdal.Translate(tileFile,srcDS, options=transop)

sys.exit(0)
