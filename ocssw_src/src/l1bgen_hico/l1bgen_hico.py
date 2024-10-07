#! /usr/bin/env python3
from hico.HicoL0toL1B import HicoL0toL1b
from hico.makePosVelQuatCsv import MakePosVelQuatCSV
from hico.cproc_hico import hico_geo
from hico.exceptions import PVQException, USGNCCoarseTimeException
from hico.exceptions import PVQInterpolation, EmptyCSVException
from hico.exceptions import BadFormatCSVException
from netCDF4 import Dataset
from datetime import datetime as DT
import argparse
import sys
import logging
import os
import pickle
from collections import namedtuple as NT
from numpy import ones as npones
from numpy import float32

__version__ = "2.1"


def ParseCommandLine(args):
    '''Specifies command line arguments and parses command line accordingly.'''
    ocvarroot = os.getenv('OCVARROOT')
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description="Generates a HICO L1B file from an L0 file.")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('-i', '--l0file', type=str, required='True',
                        help='iss*.hico.bil (required)')
    parser.add_argument('-c', '--csvfile', type=str, help='iss*.hico.csv')
    parser.add_argument('-r', '--hdrfile', type=str, help='iss*.hico.hdr')
    parser.add_argument('-l', '--lpsfile', type=str,
                        default=os.path.join(ocvarroot, 'modis', 'leapsec.dat'),
                        help=' leapsec file ')
    parser.add_argument('-e', '--earthfile', type=str,
                        default=os.path.join(ocvarroot, 'hico', 'finals.data'),
                        help=' earth-orient file ')
    parser.add_argument('-b', '--boresight', nargs=3, type=float,
                        default=([-0.9957, 0.0268, -0.0128]),
                        help=('Process Bore Sight Parameters'))
    parser.add_argument('-p', '--pvqcsv', type=str,
                        help='iss*hico_pos_vel_quat.csv (Must be a CSV file) ')
    parser.add_argument('-n', '--navoffset', action='store_true', default=False)
    parser.add_argument('-o', '--ofile', type=str,
                        help=' Output netCDF filename ')
    parser.add_argument('-d', '--debug', action='store_true', default=False)
    parsedArgs = parser.parse_args(args)
    l0basename = parsedArgs.l0file.split('.bil')[0]
    if not parsedArgs.csvfile:
        parsedArgs.csvfile = '%s.csv' % l0basename
    if not parsedArgs.hdrfile:
        parsedArgs.hdrfile = '%s.hdr' % l0basename
    if not parsedArgs.pvqcsv:
        parsedArgs.pvqcsv = '%s_pos_vel_quat.csv' % l0basename
    if not parsedArgs.ofile:
        if "L0M" in parsedArgs.l0file:
            parsedArgs.ofile = '%s.L1BM.nc' % l0basename
        else:
            parsedArgs.ofile = '%s.L1b.nc' % l0basename
    return parsedArgs


def SetLogger(pargs):
    lgrName = 'l1bgen_hico_%s_T_%s' % (DT.date(DT.now()), DT.time(DT.now()))
    lgr = logging.getLogger(lgrName)
    fmt = '%(message)s'
    if pargs.debug:
        level = logging.DEBUG
        fmt = '%(asctime)s - %(name)s - %(levelname)s -\
               [%(module)s..%(funcName)s..%(lineno)d] -\
               %(message)s'
        formatter = logging.Formatter(fmt)
        fh = logging.FileHandler('%s.log' % lgrName)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        lgr.addHandler(fh)

    level = logging.INFO
    formatter = logging.Formatter(fmt)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    lgr.addHandler(ch)
    lgr.setLevel(level)
    lgr.debug('Logger initialized')
    return lgr


def ComputeScaleOffset(dataMin, dataMax, bits=16):
    scale_factor = (dataMax - dataMin) / (2 ** bits - 1)
    add_offset = dataMin + 2 ** (bits - 1) * scale_factor
    return(scale_factor, add_offset)


def GetLocation(scene_id, logger):
    root = os.getenv('OCSSWROOT')
    locfile = os.path.join(root, 'var/hico/HICO_ID_SCENE_NAME_dict.pkl')
    try:
        with open(locfile, 'rb') as pklf:
            loc_dict = pickle.load(pklf)
            location = loc_dict.get(str(scene_id), 'unknown')

    except FileNotFoundError as e:
        logger.warning(e)

    finally:
        if not location:
            location = 'unknown'
        return location


def FillNC(root_grp_ptr, scene_location):
    retGps = NT("returnGroups",  "calGrp, productsGrp, navGrp, slaGrp, periodGrp")
    root_grp_ptr.createDimension('samples', 512)
    root_grp_ptr.createDimension('scan_lines', 2000)
    root_grp_ptr.createDimension('bands', 128)
    root_grp_ptr.instrument = 'HICO'
    root_grp_ptr.institution = 'NASA Goddard Space Flight Center'
    root_grp_ptr.resolution = '100m'
    root_grp_ptr.location_description = scene_location
    root_grp_ptr.license = 'https://science.nasa.gov/earth-science/earth-science-data/data-information-policy/'
    root_grp_ptr.naming_authority = 'gov.nasa.gsfc.sci.oceandata'
    root_grp_ptr.date_created = DT.strftime(DT.utcnow(), '%Y-%m-%dT%H:%M:%SZ')
    root_grp_ptr.creator_name = 'NASA/GSFC'
    root_grp_ptr.creator_email = 'data@oceancolor.gsfc.nasa.gov'
    root_grp_ptr.publisher_name = 'NASA/GSFC'
    root_grp_ptr.publisher_url = 'https://oceancolor.gsfc.nasa.gov'
    root_grp_ptr.publisher_email = 'data@oceacolor.gsfc.nasa.gov'
    root_grp_ptr.processing_level = 'L1B'
    nav_grp = root_grp_ptr.createGroup('navigation')
    nav_vars = list()
    nav_vars.append(nav_grp.createVariable('sensor_zenith', 'f4', ('scan_lines', 'samples',)))
    nav_vars.append(nav_grp.createVariable('solar_zenith', 'f4', ('scan_lines', 'samples',)))
    nav_vars.append(nav_grp.createVariable('sensor_azimuth', 'f4', ('scan_lines', 'samples',)))
    nav_vars.append(nav_grp.createVariable('solar_azimuth', 'f4', ('scan_lines', 'samples',)))
    nav_vars.append(nav_grp.createVariable('longitude', 'f4', ('scan_lines', 'samples',)))
    nav_vars.append(nav_grp.createVariable('latitude', 'f4', ('scan_lines', 'samples',)))
    for var in nav_vars:
        var.units = 'degrees'
        if var.name.endswith('zenith'):
            var.valid_min = 0
        else:
            var.valid_min = -180
        var.valid_max = 180
        var.long_name = var.name.replace('_', ' ').rstrip('s')
    retGps.navGrp = nav_grp
    retGps.productsGrp = root_grp_ptr.createGroup('products')
    lt = retGps.productsGrp.createVariable('Lt', 'u2', ('scan_lines',
                                                        'samples', 'bands'))
    lt.scale_factor = float32([0.02])
    lt.add_offset = float32(0)
    lt.units = "W/m^2/micrometer/sr"
    # lt.valid_range = nparray([0, 16384], dtype='u2')
    lt.long_name = "HICO Top of Atmosphere"
    lt.wavelength_units = "nanometers"
    # lt.createVariable('fwhm', 'f4', ('bands',))
    lt.fwhm = npones((128,), dtype='f4') * -1
    # wv = lt.createVariable('wavelengths', 'f4', ('bands',))
    lt.wavelengths = npones((128,), dtype='f4')
    lt.wavelength_units = "nanometers"
    retGps.slaGrp = root_grp_ptr.createGroup('scan_line_attributes')
    retGps.slaGrp.createVariable('scan_quality_flags', 'u1', ('scan_lines',
                                                              'samples'))
    # Create metadata group and sub-groups
    meta_grp = root_grp_ptr.createGroup('metadata')
    pl_info_grp = meta_grp.createGroup("FGDC/Identification_Information/Platform_and_Instrument_Identification")
    pl_info_grp.Instrument_Short_Name = "hico"
    prc_lvl_grp = meta_grp.createGroup("FGDC/Identification_Information/Processing_Level")
    prc_lvl_grp.Processing_Level_Identifier = "Level-1B"
    retGps.periodGrp = meta_grp.createGroup("FGDC/Identification_Information/Time_Period_of_Content")
    # fill HICO group
    retGps.calGrp = meta_grp.createGroup("HICO/Calibration")
    return retGps


def HicoMain(args):
    '''
    Coordinates calls to hico l0 to l1b processes and manages data recording
    in a netcdf file.
    '''
    # Get command line
    pArgs = ParseCommandLine(args)
    mainLogger = SetLogger(pArgs)
    mainLogger.info("l1bgen_hico %s" % __version__)
    pvqcsv = ''.join(pArgs.pvqcsv)
    navoffset = pArgs.navoffset
    earth_orient_file = ''.join(pArgs.earthfile)
    leap_sec_file = ''.join(pArgs.lpsfile)
    # instantiate, initialize and populate the hico level converter (hlc)
    hc = HicoL0toL1b(pArgs, parentLoggerName=mainLogger.name)
    # pass the hlc object to the pvq maker
    try:
        if not os.path.exists(pvqcsv):

            MakePosVelQuatCSV(hc, outputFileName=pvqcsv,
                              doNavTimeCorrection=navoffset)
        assert os.path.exists(pvqcsv)
    # convert raw count to data
        hc.ConvertRaw2Rad()
        # geolocation
        mainLogger.info("Running geolocation...")

        hicogeo = hico_geo(pvqcsv, earth_orient_file, leap_sec_file,
                           pArgs.boresight)
        mainLogger.info("Writing to netcdf...")
        # Record processed data
        sceneID = str(hc.L0.header['ID'])
        scene_location = GetLocation(sceneID, mainLogger)
        with Dataset(pArgs.ofile, 'w', format='NETCDF4') as root_grp:
            ncGroups = FillNC(root_grp, scene_location=scene_location)
            hc.WriteRadFile(ncGroups.productsGrp, ncGroups.periodGrp)
            hicogeo.write_geo_nc(ncGroups.navGrp, ncGroups.calGrp)
            mainLogger.info("NC file created.")
    except BadFormatCSVException:
        mainLogger.error("Failed to complete:")
        mainLogger.error("CSV file is badly formatted")
        mainLogger.error("unable to resolve issue.")
    except AssertionError:
        mainLogger.error("Failed to create PVQ file.")
    except EmptyCSVException:
        mainLogger.error("Failed to complete:")
        mainLogger.error("CSV file is empty.")
    except PVQException:
        mainLogger.error("Failed to complete:")
        mainLogger.error("Failed geolocation")
    except UserWarning as uw:
        mainLogger.warning(uw)
    except USGNCCoarseTimeException:
        mainLogger.error("Failed to complete:")
        mainLogger.error("No overlap between coarse PV & epoch timestamps.")
        mainLogger.error("Likely due to erroneous time data in CSV file.")
        mainLogger.error("Cannot create *.hico_pos_vel_quat.csv file")
    except PVQInterpolation as pvqinterr:
        mainLogger.error("interpolation error due to %s field in CSV file",
                         pvqinterr)
    except Exception as e:
        mainLogger.error("Failed to complete:")
        mainLogger.error("Undiagnosed failure")
        mainLogger.error("%s" %e, exc_info=True)
    else:
        mainLogger.info("Finished processing %s" % pArgs.l0file)


if __name__ == '__main__':
    HicoMain(sys.argv[1:])
