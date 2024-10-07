#! /usr/bin/env python3

import sys
import numpy as np
import datetime
import os
import argparse
import netCDF4
import logging
from io import BytesIO

from telemetry.PacketUtils import readFileHeader, pad_packet
from telemetry.timestamp import *
from telemetry.derive_orbitparams import derive_orbitparams
from telemetry import ccsdspy

__version__ = '1.5.4 (2024-05-28)'

ignored_apids = [
    636,  # OCI Ancillary packet
    700,  # OCI science packet with spin number but no time field
    720,  # OCI SWIR science packet
    751,  # HARP science packet
    848,  # SPEX science packet
]  # probably don't need this.

# define required packets
att = 108
orb = 128
tilt = 198
good_apids = [att, orb, tilt]

# define packet categories
# keep these in APID order
pktypes = {}
pktypes['SC']   = {'maxapid': 549, 'maxlen': 2048, 'name': 'Spacecraft'}
pktypes['OCI']  = {'maxapid': 749, 'maxlen': 1618, 'name': 'OCI'}
pktypes['HARP'] = {'maxapid': 799, 'maxlen':   40, 'name': 'HARP2'}
pktypes['SPEX'] = {'maxapid': 849, 'maxlen':  298, 'name': 'SPEXone'}

# constant valid max for seconds of day
MAX_SECONDS_OF_DAY = 172800.0   # covers 2 days 


def main():
    print("hktgen_pace", __version__)

    # Read command line options
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description='Convert S-band housekeeping telemetry to NetCDF4',
        epilog='''
EXIT Status:
0   : Success
1   : Fatal error
101 : Non-fatal file open error
102 : Invalid file or instrument from CFE header
103 : Invalid packet header
104 : Invalid packet [header/datatype/length]
110 : No valid packets found
120 : Multiple warnings; see log
'''
    )
    parser.add_argument('ifile', type=str,
                        help='path to S-band telemetry (HSK) file OR list of input files, '
                        'one per line, in chronological order')
    parser.add_argument("-o", "--output", metavar="ofile", dest="ofile",
                        help="output NetCDF4 file; defaults to PACE.yyyymmddThhmmss.HKT.nc")
    parser.add_argument('-s', '--skip_raw', action='store_true',
                        default=False, help="don't save raw housekeeping data")
    parser.add_argument('-v', '--verbose', action='store_true',
                        default=False, help="print status messages")
    parser.add_argument("--pversion", metavar="pversion", dest="pversion",
                        help="Processing version of the program")
    parser.add_argument("--doi", metavar="doi", dest="doi",
                        help="Identifier_product_doi")

    args = parser.parse_args()
    status = 0

    loglevel = logging.INFO if args.verbose else logging.WARNING
    logging.basicConfig(format='%(levelname)s: %(message)s', level=loglevel)

    # read packet definitions
    pktDir = os.path.join(os.getenv('OCDATAROOT'), 'telemetry', 'pace', 'hkt')
    if not os.path.exists(pktDir):
        logging.error(f'ERROR: The directory {pktDir} does not exist.')
        return 1
    pktDict = {}
    for apid in good_apids:
        pktDict[apid] = {}
        pktDict[apid]['recordlist'] = []
        csvfile = os.path.join(f'{pktDir}',f'APID{apid}.csv')
        if not os.path.exists(csvfile):
            logging.error(f'ERROR: The file {csvfile} does not exist.')
            return 1
        pktDict[apid]['packetdef'] = ccsdspy.FixedLength.from_file(csvfile)

    # add conversions  (TODO: specify conversions in csv for ccsdspy)
    pace_tsade_enc_tilt = ccsdspy.converters.LinearConverter(0.0075, -245.76)
    pktDict[tilt]['packetdef'].add_converted_field('TILT_ENCPOS', 'TILT_ENCPOS',
                                                   pace_tsade_enc_tilt)

    # Is input file tlm or list of files?
    filelist = []
    infile = os.path.expandvars(args.ifile)

    try:
        with open(infile, mode='rt') as flist:  # try to read as list of files
            try:
                input_list = True
                for ifile in flist:
                    filelist.append(os.path.expandvars(ifile.rstrip()))
            except UnicodeDecodeError:
                input_list = False  # contains invalid byte - infile is binary

    except IOError as e:
        logging.error(f'{e}; exiting')
        return 1

    if not len(filelist):  # if that didn't work, it's a single HSK file
        filelist.append(infile)

    # Initialize record and packet lists
    for val in pktypes.values():
        val['data'] = []
    oversize_packets = bytearray()
    alltimes = []

    # Step through all input files
    for filename in filelist:
        logging.info(f'Reading {filename}')

        try:
            ifile = open(filename, mode='rb')
        except IOError as e:
            status = 120 if status > 0 else 101
            logging.warning(f'{e}; continuing')
            continue      # input_list errors already caught above

        # Read any file header(s)
        filehdr = readFileHeader(ifile)
        if filehdr:
            logging.info(filehdr)
            logging.info('')

            # Is it the right kind of file?
            desired = (filehdr['subtype'] == 101 and
                       filehdr['length'] == 64 and
                       filehdr['SCID'] == b'PACE' and
                       filehdr['processorid'] in (1, 2, 30))

        if not filehdr or not desired:
            status = 120 if status > 0 else 102
            if input_list:
                logging.warning(f'File {filename} has invalid header; continuing')
                continue    # go to next file
            else:
                logging.error(f'File {filename} has invalid header; returning')
                return status

        # read CCSDS packets
        for packet in ccsdspy.utils.iter_packet_bytes(ifile, include_primary_header=True):

            data   = packet[6:]
            header = ccsdspy.utils.read_primary_headers(BytesIO(packet))
            for k, v in header.items():
                header[k] = v[0]  # remove outer array
            logging.info(header)

            # check for invalid header
            if header['CCSDS_PACKET_LENGTH'] > 16378:
                status = 120 if status > 0 else 103
                logging.warning(
                    f'File {filename} contains invalid CCSDS packet header: {header}')
                continue  # go to next packet

            # check for truncated data
            if len(data) < header['CCSDS_PACKET_LENGTH'] + 1:
                status = 120 if status > 0 else 104
                logging.warning(f'File {filename} has unexpected EOF: '
                                f'expected {header["CCSDS_PACKET_LENGTH"]+1} more bytes, '
                                f'got {len(data)}')
                break    # done with this file

            # check for invalid timestamp
            if (header['CCSDS_SECONDARY_FLAG'] == 1
                and data[0] > 112   # after  2017-07-18T05:49:15Z
                and data[0] < 192): # before 2060-01-28T16:50:35Z
                ts = readTimestamp(data)
                alltimes.append(ts)
                header['timestamp'] = ts
                # keep going to save packet, regardless of timestamp

            # parse useful fields
            apid = header['CCSDS_APID']
            if apid in good_apids:
                myDict = (pktDict[apid]['packetdef']).load(BytesIO(packet))
                for key, val in myDict.items():
                    myDict[key] = val[0]  # remove outer array
                myDict['timestamp'] = decode_timestamp(myDict['seconds'],
                                                       myDict['subsecs'])
                # TODO: implement custom converter in CCSDSPy
                pktDict[apid]['recordlist'].append(myDict)
                logging.info(f"taiTime={myDict['taiTime']} seconds={myDict['seconds']} "
                             f"subsecs={myDict['subsecs']} timestamp={myDict['timestamp']}")

            # accumulate raw packets by category
            if not args.skip_raw and apid not in ignored_apids:
                for val in pktypes.values():
                    if apid <= val['maxapid']:
                        if len(packet) > val['maxlen']:
                            oversize_packets += packet
                        else:
                            packet = pad_packet(packet, val['maxlen'])
                            val['data'].append(packet)
                        break  # add packet to only one category

        # end (for packet in packets)

        # close input file
        ifile.close()

    # end (for filename in filelist)

    # get start and end times
    if len(alltimes) == 0:
        logging.warning('No input packets with valid times')
        return 110

    stime = tai58_as_datetime(alltimes[0])
    etime = tai58_as_datetime(alltimes[-1])
    daystart = stime.replace(hour=0, minute=0, second=0, microsecond=0)
    timeunits = f"seconds since {daystart.strftime('%Y-%m-%d')}"

    # construct product name
    prod_name = stime.strftime('PACE.%Y%m%dT%H%M%S.HKT.nc')
    if args.ofile is None:
        args.ofile = prod_name

    # create new netcdf file
    try:
        logging.info(f'Writing {args.ofile}')
        ofile = netCDF4.Dataset(args.ofile, 'w')
    except BaseException:
        logging.error("Cannot write file \"%s\": exiting." % args.ofile)
        return 1

    # define dimensions
    ofile.createDimension('vector_elements', 3)
    ofile.createDimension('quaternion_elements', 4)

    # write raw housekeeping data to file
    if not args.skip_raw:
        group = ofile.createGroup('housekeeping_data')
        for key, val in pktypes.items():
            npkts = len(val['data'])
            if npkts > 0:
                maxlen = val['maxlen']
                logging.info(
                    f"{npkts}\t{key}_HKT_packets\tx {maxlen}\t= {npkts*maxlen}\tbytes")
                ofile.createDimension(f'{key}_hkt_pkts', None)
                ofile.createDimension(f'max_{key}_packet', maxlen)
                var = group.createVariable(f'{key}_HKT_packets', 'u1',
                                           (f'{key}_hkt_pkts', f'max_{key}_packet'))
                var.long_name = f"{val['name']} housekeeping telemetry packets"
                var[:] = val['data']
        if len(oversize_packets) > 0:
            logging.info(f'{len(oversize_packets)}\toversize_packets')
            ofile.createDimension('os_pkts', None)
            var = group.createVariable('oversize_packets', 'u1', ('os_pkts'))
            var.long_name = 'Buffer for packets exceeding maximum size'
            var[:] = oversize_packets

    # write navigation data to file
    group = ofile.createGroup('navigation_data')

    if len(pktDict[att]['recordlist']) > 0:     # ATTITUDE
        ofile.createDimension('att_records', None)

        var = group.createVariable(  # att_time
            'att_time', 'f8', ('att_records'), fill_value=-32767.)
        var.long_name = "Attitude sample time (seconds of day)"
        var.valid_min = 0
        var.valid_max = MAX_SECONDS_OF_DAY
        var.units = timeunits
        var[:] = [seconds_since(rec['Est_Time'], basetime=daystart)
                  for rec in pktDict[att]['recordlist']]

        var = group.createVariable(  # att_quat
            'att_quat', 'f4', ('att_records', 'quaternion_elements'), fill_value=-32767.)
        var.long_name = "Attitude quaternions (J2000 to spacecraft)"
        var.valid_min = -1
        var.valid_max = 1
        var.units = "seconds"
        var[:] = [rec['q_EciToBrf_Est'] for rec in pktDict[att]['recordlist']]

        var = group.createVariable(  # att_rate
            'att_rate', 'f4', ('att_records', 'vector_elements'), fill_value=-32767.)
        var.long_name = "Attitude angular rates in spacecraft frame"
        var.valid_min = np.array((-0.004), 'f4')
        var.valid_max = np.array(( 0.004), 'f4')
        var.units = "radians/second"
        var[:] = [rec['w_EciToBrf_Brf_Est']
                  for rec in pktDict[att]['recordlist']]

    if len(pktDict[orb]['recordlist']) > 0:     # EPHEMERIS (orbit)
        ofile.createDimension('orb_records', None)

        # calculate orbit parameters
        posr = [np.matmul(rec['DCM_ecef2eci'], rec['scPosJ2000'])
                for rec in pktDict[orb]['recordlist']]
        velr = [np.matmul(rec['DCM_ecef2eci'], rec['scVelJ2000'])
                for rec in pktDict[orb]['recordlist']]
        orbitparams = derive_orbitparams(posr, velr)

        var = group.createVariable(  # orb_time
            'orb_time', 'f8', ('orb_records'), fill_value=-32767.)
        var.long_name = "Orbit vector time (seconds of day)"
        var.valid_min = 0
        var.valid_max = MAX_SECONDS_OF_DAY
        var.units = timeunits
        var[:] = [seconds_since(rec['FSWTime'], basetime=daystart)
                  for rec in pktDict[orb]['recordlist']]

        var = group.createVariable(  # orb_pos
            'orb_pos', 'f4', ('orb_records', 'vector_elements'), fill_value=-9999999)
        var.long_name = "Orbit position vectors (ECR)"
        var.valid_min = -7200000
        var.valid_max = 7200000
        var.units = "meters"
        var[:] = orbitparams['posr']

        var = group.createVariable(  # orb_vel
            'orb_vel', 'f4', ('orb_records', 'vector_elements'), fill_value=-32767.)
        var.long_name = "Orbit velocity vectors (ECR)"
        var.valid_min = -7600
        var.valid_max = 7600
        var.units = "meters/second"
        var[:] = orbitparams['velr']

        var = group.createVariable(  # orb_lon
            'orb_lon', 'f8', ('orb_records'), fill_value=-32767.)
        var.long_name = "Orbit longitude (degrees East)"
        var.valid_min = -180
        var.valid_max = 180
        var.units = "degrees_east"
        var[:] = orbitparams['lon']

        var = group.createVariable(  # orb_lat
            'orb_lat', 'f8', ('orb_records'), fill_value=-32767.)
        var.long_name = "Orbit latitude (degrees North)"
        var.valid_min = -90
        var.valid_max = 90
        var.units = "degrees_north"
        var[:] = orbitparams['lat']

        var = group.createVariable(  # orb_alt
            'orb_alt', 'f8', ('orb_records'), fill_value=-32767.)
        var.long_name = "Orbit altitude"
        var.valid_min = 670000
        var.valid_max = 710000
        var.units = "meters"
        var[:] = orbitparams['alt']

    if len(pktDict[tilt]['recordlist']) > 0:     # TILT
        ofile.createDimension('tilt_records', None)

        var = group.createVariable(  # tilt_time
            'tilt_time', 'f8', ('tilt_records'), fill_value=-32767.)
        var.long_name = "Tilt sample time (seconds of day)"
        var.valid_min = 0
        var.valid_max = MAX_SECONDS_OF_DAY
        var.units = timeunits
        var[:] = [seconds_since(rec['timestamp'], basetime=daystart)
                  for rec in pktDict[tilt]['recordlist']]

        var = group.createVariable(  # tilt
            'tilt', 'f4', ('tilt_records'), fill_value=-32767.)
        var.long_name = "Tilt angle"
        var.valid_min = np.array((-20.1), 'f4')
        var.valid_max = np.array(( 20.1), 'f4')
        var.units = "degrees"
        var[:] = [rec['TILT_ENCPOS'] for rec in pktDict[tilt]['recordlist']]

        var = group.createVariable(  # tilt
            'tilt_flag', 'u1', ('tilt_records'), fill_value=255)
        var.long_name = "Tilt quality flag"
        var.flag_values = np.array([0, 1], 'u1')
        var.flag_meanings = "Valid Not_initialized"
        var[:] = [1 - rec['TILT_HOME_INIT']
                  for rec in pktDict[tilt]['recordlist']]

    # write global metadata

    # static
    ofile.title = "PACE HKT Data"
    ofile.instrument = "Observatory"
    ofile.processing_version = "V1.0"
    ofile.institution = "NASA Goddard Space Flight Center, Ocean Biology Processing Group"
    ofile.license = "https://science.nasa.gov/earth-science/earth-science-data/data-information-policy/"
    ofile.naming_authority = "gov.nasa.gsfc.oceancolor"
    ofile.creator_name = "NASA/GSFC/OBPG"
    ofile.creator_email = "data@oceancolor.gsfc.nasa.gov"
    ofile.creator_url = "https://oceancolor.gsfc.nasa.gov"
    ofile.project = "Ocean Biology Processing Group"
    ofile.publisher_name = "NASA/GSFC/OB.DAAC"
    ofile.publisher_email = "data@oceancolor.gsfc.nasa.gov"
    ofile.publisher_url = "https://oceancolor.gsfc.nasa.gov"
    ofile.processing_level = "L0"
    ofile.Conventions = "CF-1.8 ACDD-1.3"
    ofile.standard_name_vocabulary = 'CF Standard Name Table v79'
    ofile.keywords_vocabulary = "NASA Global Change Master Directory (GCMD) Science Keywords"
    # dynamic
    ofile.product_name = args.ofile
    ofile.time_coverage_start = datetime_repr(stime)
    ofile.time_coverage_end = datetime_repr(etime)
    ofile.history = ' '.join([v for v in sys.argv])
    if input_list:
        ofile.source = ','.join([os.path.basename(f) for f in filelist])
    ofile.date_created = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')

    # if pversion and doi are specified, write the, to the global attributes 
    if args.pversion != None:
        ofile.processing_version = args.pversion

    if args.doi != None:
        ofile.identifier_product_doi_authority = "https://dx.doi.org"
        ofile.identifier_product_doi = args.doi

    # write processing control data into file
    group = ofile.createGroup('processing_control')
    group.setncattr("software_name", "hktgen_pace")
    group.setncattr("software_version", __version__)
    group.setncattr("hsk_files", ",".join(filelist))

    # write input parameters into file, including default files that the
    # program used and not specified by the user
    group1 = group.createGroup("input_parameters")
    group1.setncattr("ifile", args.ifile)
    group1.setncattr("ofile", args.ofile)
    group1.setncattr("skip_raw", "True" if args.skip_raw  else "False")

    # close file
    ofile.close()

    if status:
        logging.warning('Exiting with status code %s' % status)
    return status


if __name__ == '__main__':
    sys.exit(main())
