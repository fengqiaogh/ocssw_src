#!/usr/bin/env python3
import argparse
import sys
import re
import subprocess
import os
import seadasutils.MetaUtils as MetaUtils
from seadasutils.ParamUtils import ParamProcessing

def getBinRes(resolution):
    resvalue = 2000

    if "km" in resolution:
        resvalue = float(resolution.strip('km')) * 1000.
    elif "m" in resolution:
        resvalue = float(resolution.strip('m'))
    elif "deg" in resolution:
        resvalue = float(resolution.strip('deg')) * 111312.

    if resvalue <= 99.0 :
        return 'HH'
    elif resvalue < 250.0:
        return 'HQ'
    elif resvalue < 500.0:
        return 'Q'
    elif resvalue < 1150.0:
        return 'H'
    elif resvalue < 2300.0:
        return '1'
    elif resvalue < 4600.0:
        return '2'
    elif resvalue < 9200.  :
        return '4'
    else:
        return '9'

if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,description='''\
     This program takes a product from a L2 file, maps it using a Plate Carree cylindrical projection,
      and produces a gray scale PGM or color PPM file.

      The arguments can be specified on the commandline, or put into a parameter file,
      or the two methods can be used together, with commandline over-riding.''',add_help=True)
    parser.add_argument('--parfile', type=str, help='input parameter file')
    parser.add_argument('ifile', nargs='?', type=str, help='input file (L1A/B file or file with list)')
    parser.add_argument('--geofile', nargs='?', type=str,help='input file geolocation file name')
    parser.add_argument('--ofile', type=str, help='output file name; default = <ifile>.L2_MAP.<oformat extension>')
    parser.add_argument('--product', type=str, help='product(s) to map; comma separated')
    parser.add_argument('--resolution', type=str, default='2.0km', help='''\

    resolution
    -----------------
         #.#:  resolution in meters
       #.#km:  resolution in kilometers
      #.#deg:  resolution in degrees
    ''')
    parser.add_argument('--oformat', choices=['netcdf4','hdf4','png','ppm','tiff'], type=str, default="png", help=('''\
     output file format
     --------------------------------------------------------
     netcdf4: netCDF4 file, can contain more than one product
     hdf4:    HDF4 file (old SMI format)
     png:     PNG image file
     ppm:     PPM image file
     tiff:    TIFF file with georeference tags
     '''))
    parser.add_argument('--north', type=float, help=('Northern most Latitude (-999=file north)'))
    parser.add_argument('--south', type=float, help=('Southern most Latitude (-999=file south)'))
    parser.add_argument('--east', type=float, help=('Eastern most Latitude (-999=file east)'))
    parser.add_argument('--west', type=float, help=('Western most Latitude (-999=file west)'))
    parser.add_argument('--projection', type=str,default='platecarree',  help='''\
        enter a proj.4 projection string or choose one of the following
        predefined projections:
        --------------------------------------------------------------
        smi:       Standard Mapped image, cylindrical projection, uses
                   central_meridian.  n,s,e,w default to whole globe.
                   projection="+proj=eqc +lat_0=<central_meridian>"
        platecarree: Plate Carree image, cylindrical projection, uses
                   central_meridian
                   projection="+proj=eqc +lat_0=<central_meridian>"
        mollweide: Mollweide projection
                   projection="+proj=moll +lat_0=<central_meridian>"
        lambert:   Lambert conformal conic projection
                   projection="+proj=lcc +lat_0=<central_meridian>"
        albersconic: Albers equal-area conic projection
                   projection="+proj=aea +lat_0=<central_meridian>"
        mercator:  Mercator cylindrical map projection
                   projection="+proj=merc +lat_0=<central_meridian>"
        ease2:     Ease Grid 2 projection
                   projection="+proj=cea +lon_0=0 +lat_ts=30 +ellps=WGS84
                         +datum=WGS84 +units=m +lat_0=<central_meridian>"
        conus:     USA Contiguous Albers Equal Area-Conic projection, USGS"
                   projection="+proj=aea +lat_1=29.5 +lat_2=45.5
                         +lat_0=23.0 +lon_0=-96 +x_0=0 +y_0=0
                         +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
    ''')
    parser.add_argument('--central_meridian', type=float, help='central meridian for projection in deg east.  Only used for smi, mollweide and raw projection')
    parser.add_argument('--palfile', type=str, help='palette filename. Default uses file for product in product.xml')
    parser.add_argument('--fudge', type=float,default=1.0, help='fudge factor used to modify size of L3 pixels')
    parser.add_argument('--threshold', type=float,default=0, help='minimum percentage of filled pixels before an image is generated')
    parser.add_argument('--datamin', type=float, help=('minimum value for scaling (default from product.xml)'))
    parser.add_argument('--datamax', type=float, help=('maximum value for scaling (default from product.xml)'))
    parser.add_argument('--scaletype', choices=['linear','log','arctan'],type=float, help='''\
        data scaling type (default from product.xml)
        ---------------------------------------------
        linear:  linear scaling
        log:     logarithmic scaling
        arctan:  arc tangent scaling

    ''')
    parser.add_argument('--keep-intermediates', action='store_true',default=False, help='Do not delete the intermediate L2/L3B files produced')
    parser.add_argument('-v','--verbose',action='count',default=0, help="Let's get chatty.  Two -v's a better than one.")

    args=parser.parse_args()

    if args.parfile:
        param_proc = ParamProcessing(parfile=args.parfile)
        param_proc.parseParFile()
        for param in (param_proc.params['main'].keys()):
            value = param_proc.params['main'][param]
            if hasattr(args,param):
                    setattr(args,param,value)
            elif param == 'file':
                args.ifile = value

    if not args.ifile:
        parser.error("You must specify an input file either on the command line or in a parameter file")

    if not args.product:
        parser.error("You must specify at least one product to map, either on the command line or in a parameter file")

    if not args.ofile:
        extension = args.oformat
        if extension == 'netcdf4':
            extension = 'nc'
        setattr(args,'ofile',args.ifile + '.L2_MAP.' + extension)

    default_opts=["product", "datamin", "datamax", "scaletype", "palfile", "north", "south", "east", "west", "central_meridian"]
    l2bin_opts = ["product", "resolution"]
    l3map_opts = ["product", "resolution", "oformat", "north", "south", "west", "east", "projection", "central_meridian", "palfile", "fudge", "threshold"]

    tmpfile = args.ifile + '.L3b'

    if args.verbose > 1:
        print(args)

    # Build the l2bin command line
    clo = ["l2bin"]
    clo.append('flaguse=ATMFAIL,CLDICE,BOWTIEDEL')
    clo.append("ifile=" + args.ifile)
    clo.append("ofile=" + tmpfile)
    clo.append("area_weighting=1")
    for l2bin_opt in l2bin_opts:
        if getattr(args,l2bin_opt):
            if l2bin_opt == 'resolution':
                clo.append(l2bin_opt + "=" + getBinRes(str(getattr(args,l2bin_opt))))
            elif l2bin_opt == 'product':
                clo.append("l3bprod" + "=" + str(getattr(args,l2bin_opt)))
            else:
                clo.append(l2bin_opt + "=" + str(getattr(args,l2bin_opt)))

    if args.verbose:
        print("l2bin parameters:")
        print(clo)

    try:
        if args.verbose > 1:
            subprocess.run(clo,check=True)
        else:
            subprocess.run(clo,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=True)
        if args.verbose:
            print("Generated bin file %s" % tmpfile)

    except subprocess.CalledProcessError as e:
        print("Process error ({0}): message:{1}".format(str(e.returncode), e.output))
        sys.exit()

    # Build the l3mapgen command line
    clo = ["l3mapgen"]
    clo.append('ifile=' + tmpfile)
    clo.append('interp=nearest')
    for opt in l3map_opts:
        if getattr(args,opt):
            if type(getattr(args,opt)) is bool: # handle boolean options
                clo.append(opt + "=1" )
            elif opt == 'product':
                if args.oformat == 'netcdf4':
                    clo.append("product=" + args.product)
            else:
                clo.append(opt + "=" + str(getattr(args,opt)))

    if args.verbose:
        print("l3mapgen parameters:")
        print(clo)

    if args.oformat == 'netcdf4':
        try:
            clo.append('ofile=' + args.ofile)
            if args.verbose > 1:
                subprocess.run(clo,check=True)
            else:
                subprocess.run(clo,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=True)
            if args.verbose:
                print("Generated map file %s" % args.ofile)
        except subprocess.CalledProcessError as e:
            print("Process error ({0}): message:{1}".format(str(e.returncode), e.output))
            sys.exit()
    else:
        for prod in args.product.split(','):
            prod_clo = clo
            prod_clo.append("product=" + prod)
            ofile = args.ofile

            if len(args.product.split(',')) > 1:
                ofile = args.ofile.replace(args.oformat,'') + prod + '.' + args.oformat

            prod_clo.append('ofile=' + ofile)
            try:
                if args.verbose > 1:
                    subprocess.run(prod_clo,check=True)
                else:
                    subprocess.run(prod_clo,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=True)
                if args.verbose:
                    print("Generated map file %s" % ofile)
            except subprocess.CalledProcessError as e:
                print("Process error ({0}): message:{1}".format(str(e.returncode), e.output))
                sys.exit()

    if not args.keep_intermediates:
        os.remove(tmpfile)
        if args.verbose:
            print("Deleted %s" % tmpfile)