#!/usr/bin/env python3
import argparse
from datetime import datetime
import logging
from math import ceil, floor
import os
import re
from shlex import quote
import subprocess
import sys

import seadasutils.MetaUtils as MetaUtils
from seadasutils.ParamUtils import ParamProcessing
from seadasutils.setupenv import build_executable_path
from seadasutils.SensorUtils import by_sensor
import mlp.get_obpg_file_type

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

def get_ftype(filepath):
    level = 1
    cmd = ['ncdump', '-h', filepath]
    result = subprocess.run(cmd, shell=False, capture_output=True, text=True)
    levelstr = None
    if not result.returncode:
        for line in result.stdout.splitlines():
            if 'processing_level' in line:
                levelstr = line
                break
    if levelstr:
        if re.search(r'"L2',levelstr):
            level = 2
        elif re.search(r'"L3',levelstr):
            level = 3

    return level

def whoops(errormsg, status=1):
    logging.error(errormsg)
    sys.exit(status)

def execute_command(command):
    """
    Execute a process on the system
    """
    logging.info("Running: " + ' '.join(command))

    result = subprocess.run(command, shell=False, capture_output=True, text=True)
    std_out = result.stdout
    err_out = result.stderr

    if result.returncode:
        logging.debug(std_out)
        if err_out:
            logging.error(err_out)
        whoops("Exiting: {command} returned a status {status}".format(command=' '.join(command),status=result.returncode), result.returncode)
    else:
        if std_out:
            logging.debug(std_out)
        if err_out:
            logging.debug(err_out)

def run_clo_program(program, inputs):
    """
    Set up and run a program that accepts keyword=value command line arguemnts ala CLO
    """
    logging.info('\n###### {program} ######\n'.format(program=program))
    cmd = [str(build_executable_path(program))]

    params = ["%s=%s" % (k,v) for k, v in inputs.items()]

    cmd.extend(params)
    execute_command(cmd)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,description='''\
generate mapped output from a SeaDAS supported satellite data files
arguments can be specified on the commandline or in a parameter file
the two methods can be used together, with commandline over-riding the parfile''',add_help=True)
    parser.add_argument('--parfile', '-p', type=str, help='input parameter file')
    parser.add_argument('--ifile', '-i', type=str, help='input file or text file list of files')
    parser.add_argument('--geofile', '-g', type=str,help='geolocation file or text file list of files')
    parser.add_argument('--ofile', '-o', type=str, help='output file name; default: <ifile>.MAP.<oformat ext>')
    parser.add_argument('--logfile','-l', type=str, default='mapgen_<timestamp>.log', help='''\
log file
default: mapgen_<timestamp>.log
<timestamp> is in seconds since Jan 1, 1970 00:00:00
this file is deleted if verbose is not set and no errors
occur during processing''')

    parser.add_argument('--use_rgb', action='store_true', default=False, help='''\
generate an RGB image output
default: a pseudo-true color image with bands to use
         controlled by --product_rgb option''')
    prodGroup = parser.add_mutually_exclusive_group()
    prodGroup.add_argument('--product', type=str, help='product(s) to map; comma separated')
    prodGroup.add_argument('--product_rgb', type=str,help='''\
comma separated string of RGB products
e.g., product_rgb=rhos_645,rhos_555,rhos_469
default:  sensor specific, see
$OCDATAROOT/<sensor>/l1mapgen_defaults.par''')

    parser.add_argument('--resolution', '-r', type=str, default='2.0km', help='''\
     #.#:  width of a pixel in meters
   #.#km:  width of a pixel in kilometers
  #.#deg:  width of a pixel in degrees''')
    parser.add_argument('--oformat', choices=['netcdf4','png','ppm','tiff'], type=str, default="png", help='''\
netcdf4: Network Common Data Form v4 file
           can contain more than one product
png:     Portable Network Graphics format image
ppm:     Portable PixMap format image
tiff:    Tagged Image File Format with georeference tags''')
    parser.add_argument('--use_transparency','-t',action='store_true', default=False, help='make missing data transparent\nonly valid for color PNG and TIFF output')
    parser.add_argument('--north', '-n', type=float, help=('northern-most latitude; default: input file max lßatitude'))
    parser.add_argument('--south', '-s', type=float, help=('southern-most latitude; default: input file min latitude'))
    parser.add_argument('--east', '-e', type=float, help=('eastern-most latitude; default: input file max longitude'))
    parser.add_argument('--west', '-w', type=float, help=('western-most latitude; default: input file min longitude'))
    parser.add_argument('--projection', type=str,default='platecarree',  help='''\
 "proj" projection string or one of the following:
    platecarree: Plate Carree (cylindrical) projection
      projection="+proj=eqc +lat_0=<central_meridian>"
    mollweide:   Mollweide projection
      projection="+proj=moll +lat_0=<central_meridian>"
    lambert:     Lambert conformal conic projection
      projection="+proj=lcc +lat_0=<central_meridian>"
    albersconic: Albers equal-area conic projection
      projection="+proj=aea +lat_0=<central_meridian>"
    mercator:    Mercator cylindrical map projection
      projection="+proj=merc +lat_0=<central_meridian>"
    ease2:       Ease Grid 2 projection
      projection="+proj=cea +lon_0=0 +lat_ts=30 +ellps=WGS84
      +datum=WGS84 +units=m +lat_0=<central_meridian>"''')
    parser.add_argument('--central_meridian', type=float, help='central meridian to use for projection in degrees east')
    parser.add_argument('--palfile', type=str, help='palette filename\ndefault: see $OCDATAROOT/common/product.xml')
    parser.add_argument('--fudge', type=float,default=1.0, help='factor used to modify pixel search radius for mapping')
    parser.add_argument('--datamin', type=float, help='minimum value for scaling (default from product.xml)')
    parser.add_argument('--datamax', type=float, help='maximum value for scaling (default from product.xml)')
    parser.add_argument('--scale_type', choices=['linear','log','arctan'],type=str, help='data scaling method (default from product.xml)')
    parser.add_argument('--threshold', type=float, help='minimum percentage of filled pixels for image generation\ndefault: 0')
    parser.add_argument('--trimNSEW', action='store_false',default=True, help='do not trim output to match input NSEW range')
    parser.add_argument('--write_projtext', action='store_true',default=False, help='write projection information to a text file (for mapgen_overlay script)')
    parser.add_argument('--keep-intermediates', action='store_true',default=False, help='do not delete the intermediate L2/L3B files produced')
    parser.add_argument('--verbose', '-v',action='count',default=0, help="let's get chatty; each occurrence increases verbosity\ndefault: error\n-v info -vv debug'")

    args=parser.parse_args()

    sanitize = ['ifile','ofile','logfile','geofile','parfile','product','resolution']
    for arg in vars(args):
        if arg in sanitize:
            value = getattr(args,arg)
            if value:
                setattr(args,arg,quote(value))

    if 'timestamp' in args.logfile:
        setattr(args,'logfile',''.join(['mapgen_',datetime.now().strftime('%s'),'.log']))
    # Set up the logging levels
    levels = [logging.ERROR, logging.INFO, logging.DEBUG]
    # limit verbosity to the number of levels, in case someone thinks if 2 v's are good 20 is better
    verbosity = max(0, min(args.verbose, len(levels)-1))

    logging.basicConfig(filename=args.logfile,
        format='-%(levelname)s- %(asctime)s - %(message)s',
        level=levels[verbosity])
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    if not os.getenv("OCSSWROOT"):
        whoops("You must define the OCSSW enviroment for me to help you...")

    if args.parfile:
        param_proc = ParamProcessing(parfile=args.parfile)
        param_proc.parseParFile()
        for param in (param_proc.params['main'].keys()):
            value = param_proc.params['main'][param]
            if param in sanitize:
                value = quote(value)
            if hasattr(args,param):
                setattr(args,param,value)
            elif param == 'file':
                args.ifile = value

    if not args.ifile:
        parser.print_help()
        whoops("\nYou must specify an input file either on the command line or in a parameter file")

    if not os.path.isfile(args.ifile):
        whoops("\nYou must specify a valid input file either on the command line or in a parameter file")

    if not args.ofile:
        extension = args.oformat
        if extension == 'netcdf4':
            extension = 'nc'
        setattr(args,'ofile',args.ifile + '.MAP.' + extension)        

    geo_opts = ["north", "south", "east", "west"]
    l2bin_opts = ["resolution"]
    l3map_opts = ["product", "use_rgb","product_rgb", "oformat", "resolution",
                  "projection", "central_meridian", "palfile", "fudge", "threshold",
                  "datamin", "datamax", "scale_type", "palfile", "use_transparency",
                  "trimNSEW","write_projtext"]
    tmpfile_l2gen = "mapgen_l2gen_list.txt"
    tmpfile_l2bin = "mapgen_tmp.l2bin.nc"

    if args.verbose > 1:
        msg = "Parameters passed to or assumed by mapgen:\n"
        for key,value in (vars(args)).items():
            msg += "\t{k}={v}\n".format(k=key,v=value)
        logging.debug(msg)

    ftype = None
    ifiles = []
    geofiles = []
    if 'text' in MetaUtils.get_mime_data(args.ifile) and (not re.search('MTL.txt', args.ifile)) and (not re.search('xfdumanifest.xml', args.ifile)):
        with open(args.ifile,'r') as lstfile:           
            ifile_list = lstfile.readlines()
            for ifile in ifile_list:
                if not os.path.isfile(ifile.strip()):
                    whoops("Hmmm...seems I don't recognize this input: {ifile}".format(ifile=ifile))

                ifiles.append(ifile.strip())
                if ftype: 
                    if get_ftype(ifile.strip()) != ftype:
                        whoops("Seems you are mixing input types.  While I might be able to make sense of your request, it'd be nicer to pass me just on type of file at a time...")
                else:
                    ftype = get_ftype(ifile.strip())
    else:
        ifiles.append(args.ifile)
        ftype = get_ftype(args.ifile)

    if not args.use_rgb and not args.product:
        whoops("My cousin has the artificial intelligence in the family...I'll need you to tell me what product(s) you wish to map...")
    if ftype == 3 and len(ifiles) > 1:
        whoops("I can only map one Level 3 file at a time...")
    
    # initialize the command-line option dictionary
    clo = {}
    # For L1 inputs we need to process to L2 for the rhos/t
    if ftype == 1:
        if not args.use_rgb:
            whoops("You gave me Level 1 data, but didn't tell me you wanted RGB output, and I don't want to assume... If you want another product, I'll need Level 2 inputs")

        if args.oformat == 'netcdf4':
            whoops("netCDF output is not supported for RGB images")

        # we only care about geolocation files if we need to run l2gen...
        if args.geofile:
            if 'text' in MetaUtils.get_mime_data(args.geofile):
                with open(args.geofile,'r') as geolistfile:
                    geofiles = geolistfile.readlines()
            else:
                geofiles.append(args.geofile)
            if len(geofiles) != len(ifiles):
                whoops("Number of provided geofiles does not match number of provided ifiles...")
        
        l2filelst = open(tmpfile_l2gen,'w')
        for i, l1file in enumerate(ifiles):
            if args.verbose > 1:
                logging.info("{cnt}: {ifile}".format(cnt=i, ifile=l1file))

            # Build the l2gen command
            l1file = l1file.strip()
            clo.clear()
            clo['ifile']=l1file
            if args.geofile:
                clo['geofile'] = geofiles[i].strip()

            ofile = "{filebase}_{cnt:02}.nc".format(filebase=tmpfile_l2gen.replace('_list.txt',''),cnt=i+1)
            clo['ofile'] = ofile
            if args.product_rgb:
                clo['l2prod'] = args.product_rgb
            else:
                file_typer = mlp.get_obpg_file_type.ObpgFileTyper(l1file)
                # ftype_string, sensor = file_typer.get_file_type()
                sensor = file_typer.get_sensor()
                
                dirs_by_sensor = by_sensor(sensor)
                sensor_dir = dirs_by_sensor['dir']
                sensor_sub_dir = dirs_by_sensor.get('subdir')
                if sensor_sub_dir:
                    filename_l3mapgen_defaults = os.path.join(os.getenv('OCDATAROOT'),
                                    sensor_dir, sensor_sub_dir, 'l3mapgen_defaults.par')
                    if not os.path.exists(filename_l3mapgen_defaults):
                        filename_l3mapgen_defaults = os.path.join(os.getenv('OCDATAROOT'),
                                    sensor_dir, 'l3mapgen_defaults.par')
                else: 
                    filename_l3mapgen_defaults = os.path.join(os.getenv('OCDATAROOT'),
                                    sensor_dir, 'l3mapgen_defaults.par')
                param_proc = ParamProcessing(parfile=filename_l3mapgen_defaults)
                param_proc.parseParFile()
                for param in (param_proc.params['main'].keys()):
                    if param == 'product_rgb':
                        clo['l2prod'] = args.product_rgb = param_proc.params['main']['product_rgb']
            clo['atmocor'] = '0'
            clo['proc_sst'] = '0'

            for geo_opt in geo_opts:
                if getattr(args,geo_opt) is not None:
                    clo[geo_opt] = str(getattr(args,geo_opt))

            if args.verbose > 1:
                msg = "l2gen parameters for {ofile}\n".format(ofile=ofile)
                for key,value in (vars(args)).items():
                    msg += "\t{k}={v}\n".format(k=key,v=value)
                logging.debug(msg)

            run_clo_program('l2gen', clo)
            logging.info("Generated L2 file for input: {ifile}".format(ifile=l1file))

            l2filelst.write(ofile+'\n')
        l2filelst.close()

    # Bin the inputs
    if ftype != 3:
        # Build the l2bin command line
        clo.clear()
        if args.use_rgb:
            clo['flaguse'] = 'ATMFAIL,BOWTIEDEL'
        else:
            if args.product and re.search(r'sst.?',args.product):
                clo['suite'] = 'SST'
                clo['flaguse'] = 'LAND,BOWTIEDEL'
            else:
                clo['flaguse'] = 'ATMFAIL,CLDICE,BOWTIEDEL'
                clo['l3bprod'] = args.product
        if ftype == 1:
            clo['ifile'] = tmpfile_l2gen
        else:
            clo['ifile'] = args.ifile

        clo['ofile'] = tmpfile_l2bin
        clo['area_weighting'] = '1'
        clo['prodtype'] = "regional"

        if args.north is not None:
            clo['latnorth'] = ceil(float(args.north))
        if args.south is not None:
            clo['latsouth'] = floor(float(args.south))
        if args.west is not None:
            clo['lonwest'] = floor(float(args.west))
        if args.east is not None:
            clo['loneast'] = ceil(float(args.east))
        for l2bin_opt in l2bin_opts:
            if getattr(args,l2bin_opt) is not None:
                if l2bin_opt == 'resolution':
                    clo[l2bin_opt] =  getBinRes(str(getattr(args,l2bin_opt)))
                else:
                    clo[l2bin_opt] =  str(getattr(args,l2bin_opt))

        if args.verbose > 1:
            msg = "l2bin parameters:\n"
            for key,value in (vars(args)).items():
                msg += "\t{k}={v}\n".format(k=key,v=value)
            logging.debug(msg)

        run_clo_program('l2bin', clo)
        logging.info("Generated bin file {binfile}".format(binfile=tmpfile_l2bin))

    # Build the l3mapgen command line
    clo.clear()
    if not args.use_rgb:
        clo['mask_land'] = 1
    if ftype != 3:
        clo['ifile'] = tmpfile_l2bin
    else:
        clo['ifile'] = args.ifile

    clo['interp'] = 'nearest'

    for opt in l3map_opts + geo_opts:
        if getattr(args,opt) is not None:
            # set boolean options if defined
            if type(getattr(args,opt)) is bool:
                if getattr(args,opt) is True:
                    clo[opt] = '1'
            else:
                clo[opt] = str(getattr(args,opt))

    if args.product:
        clo['product'] = args.product

    if args.oformat == 'netcdf4':
        clo['ofile'] = args.ofile
    else:
        if args.product:
            if len(args.product.split(',')) > 1:
                ofile = args.ofile.replace(args.oformat,'') + 'PRODUCT' + '.' + args.oformat
                clo['ofile'] = ofile
            else:
                 clo['product'] = args.product
                 clo['ofile'] = args.ofile
        else:
            clo['ofile'] = args.ofile

    # Generate the maps!
    if args.verbose > 1:
        msg = "l3mapgen parameters:\n"
        for key,value in (vars(args)).items():
            msg += "\t{k}={v}\n".format(k=key,v=value)
        logging.debug(msg)

    run_clo_program('l3mapgen', clo)
    logging.info("Generated map file(s)")

    # Clean up our crumbs...
    if not args.keep_intermediates:
        if os.path.exists(tmpfile_l2gen):
            l2filelst = open(tmpfile_l2gen,'r')
            for l2f in l2filelst.readlines():
                l2f = l2f.strip()
                if os.path.exists(l2f):
                    if args.verbose:
                        logging.info("removing: {rmfile}".format(rmfile=l2f))
                    os.remove(l2f)
            l2filelst.close()
            os.remove(tmpfile_l2gen)
        if os.path.exists(tmpfile_l2bin):
            if args.verbose:
                logging.info("removing: {rmfile}".format(rmfile=tmpfile_l2bin))
            os.remove(tmpfile_l2bin)

    # Remove the logfile if empty (because verbosity was zero and no errors were thrown...)
    if os.stat(args.logfile).st_size == 0:
        os.remove(args.logfile)