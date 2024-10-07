#!/usr/bin/env python3

import argparse
import subprocess
from pathlib import Path
import logging
import sys
import os
import json
from shlex import quote

logger = logging.getLogger('viirs_l1_benchmark')
workingDir = {
    'SNPP':'viirs-snpp-sample',
    'NOAA20':'viirs-noaa20-sample'
}

global artifacts
global packageroot

def build_executable_path(prog_name):
    """
    Returns the path to the program named in prog_name.
    None is returned if the program is not found.
    """
    global packageroot
    prog_path = None

    if Path.exists(packageroot / 'bin' / prog_name):
        prog_path = packageroot / 'bin' / prog_name
    else:
        err_msg = "Cannot find program %s" % prog_name
        logging.error(err_msg)
        sys.exit(err_msg)

    return prog_path

def execute_command(command,dieOnError=False):
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
            if dieOnError:
                print("Exiting: {command} returned a status {status}".format(command=' '.join(command),status=result.returncode), result.returncode)
                sys.exit(1)
    else:
        if std_out:
            logging.debug(std_out)
        if err_out:
            logging.debug(err_out)

def run_clo_program(platform,program):
    """
    Set up and run a program that accepts keyword=value command line arguemnts ala CLO
    """
    global artifacts
    logging.info('###### %s ######\n' % program)
    cmd = [str(build_executable_path(program))]

    params = ["%s=%s" % (k,v) for k, v in artifacts[platform][program].items()]

    cmd.extend(params)
    execute_command(cmd,dieOnError=True)

def run_positional_program(platform,program):
    """
    Set up and run a program that accepts positional command line arguemnts
    """
    global artifacts
    logging.info('###### %s ######\n' % program)
    cmd = [str(build_executable_path(program))]

    cmd.extend(artifacts[platform][program])
    execute_command(cmd,dieOnError=True)

def run_verify(platform,tolerance=None):
    """
    Verify the generated products against a standard set
    """
    global artifacts
    logging.info('###### Verifying outputs ######\n')
    prog = build_executable_path('nccmp')

    nccmp_tolerance = ''
    if tolerance:
        nccmp_tolerance = "-T %f" % tolerance
    seen = {}
    skipthese = ['cmn_lut_file','geo_lut_file','polar_wander_file','leapsec_file',
                 'static_lut_file','rsb_dynamic_lut_file','dnb_dynamic_lut_file','straylight_lut_file']

    for program in ('calibrate_viirs','geolocate_viirs'):
        pdict = artifacts[platform][program]
        
        for k,v in pdict.items():
            if k in skipthese:
                continue
            if not v in seen:
                seen[v] = True
                logging.info("Comparing %s output %s ..." % (k,v))
                cmd = [str(prog), '-m','-g','-d','-f','-C','10','-G','date_created,ProductionTime,software_version', nccmp_tolerance, ''.join(['../baseline/',v]), v]
                execute_command(cmd)

if __name__ == "__main__":
    global packageroot
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='''\
    This program runs a benchmark the VIIRS L1 code suite including: 
        l1agen_viirs, calibrate_viirs and geolocate_viirs
    SNPP and NOAA-20 are supported

    This script requires that the benchmark data files have previously been downloaded
    The working directory for this script should be the directory in which the benchmark data have been untar'd
    (yeah, lazy programmers)
    ''', add_help=True)
    parser.add_argument('--mission', type=str, help='Mission to test.  If not set, both SNPP and NOAA20 are run',
                    choices=['SNPP','NOAA20'], default=None)
    parser.add_argument('--packageroot',type=str, help='''Base directory for package install 
    This can also be set with the OCSSWROOT or VIIRSL1ROOT environment variables
    If both are set, OCSSWROOT is used''')
    parser.add_argument('--logfile','-l', type=str, default='viirs_l1_benchmark.log', help="output log filename")
    parser.add_argument('--check', '-c', action="store_true", default=False,
                    help='verify the generated products against the standard artifacts')
    parser.add_argument('--tolerance','-t',type=float, default=None, help="Tolerance (in percentage) to allow when running artifact verification")
    parser.add_argument('--artifacts','-a', type=str, default="viirs_l1_benchmark.json", help="JSON file containing the test artifacts to generate")
    parser.add_argument('--verbose', '-v', action="count", default=0,
                    help='each occurrence increases verbosity (default is ERROR): -v=WARNING -vv=INFO -vvv=DEBUG')
    

    args = parser.parse_args()

    sanitize = ['packageroot','logfile','artifacts']
    for arg in vars(args):
        if arg in sanitize:
            value = getattr(args,arg)
            if value:
                setattr(args,arg,quote(value))

    with open(args.artifacts) as f:
        artifacts = json.load(f)

    # Set up the logging levels
    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    logging.basicConfig(filename=args.logfile,
        format='-%(levelname)s- %(asctime)s - %(message)s',
        level=levels[args.verbose])
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    rootvar = None
    if args.packageroot:
        rootvar = args.packageroot
    else:
        rootvar = os.environ.get('OCSSWROOT')
        if not rootvar:
            rootvar = os.environ.get('VIIRSL1ROOT')
    try:
        packageroot = Path(rootvar)
        if (not packageroot.exists()) or (not packageroot.is_dir()):
            errormsg = "package root variable does not exist or is not a directory path!"
            logging.error(errormsg)
            sys.exit(errormsg)
    except:
        errormsg = "package root variable is not defined!"
        logging.error(errormsg)
        sys.exit(errormsg)

    missions = ['SNPP','NOAA20']
    if args.mission:
        missions = [args.mission]

    cwd = os.getcwd()
    for mission in missions:
        logging.info("Processing platform: %s",mission)
        os.chdir(os.path.join(cwd,workingDir[mission]))
        
        run_positional_program(mission,'l1agen_viirs')
        run_clo_program(mission,'calibrate_viirs')
        run_clo_program(mission,'geolocate_viirs')

        if args.check:
            run_verify(mission,args.tolerance)
