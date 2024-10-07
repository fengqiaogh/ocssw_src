#!/usr/bin/env python3
from telemetry import ccsdspy
from io import BytesIO
from telemetry.timestamp import *
from telemetry.PacketUtils import readFileHeader

__version__ = '3.0.0_2024-02-26'



def get_timestamp(all_packets, is_stime=True):
    '''
    Given a CCSDSPY packet list, get the start or end time
    from the packets

    PARAMETERS:
        all_packets : List of bytes
        is_stime : True for start time, False for end time 

    RETURNS:
        datetime

    '''
    # reading backwards or forward range object
    range_order = range(len(all_packets)) if is_stime else reversed(range(len(all_packets)))

    time = None
    for i in range_order:
        packet = all_packets[i]
        data = packet[6:]
        header = ccsdspy.utils.read_primary_headers(BytesIO(packet))

        for k, v in header.items():
            header[k] = v[0]  # remove outer array

        #check for invalid header
        if header['CCSDS_PACKET_LENGTH'] > 16378:
            continue  # go to next packet

        # check for truncated data
        if len(data) < header['CCSDS_PACKET_LENGTH'] + 1:
            break    # done with this file

        if (header['CCSDS_SECONDARY_FLAG'] == 1
            and data[0] > 112   # after  2017-07-18T05:49:15Z
            and data[0] < 192): # before 2060-01-28T16:50:35Z
            timestamp = readTimestamp(data)
            time = tai58_as_datetime(timestamp)
            break # only need the first valid time for forward and backwards reading
        
    return time.strftime('%Y-%m-%dT%H:%M:%S.%f')



def l0info_hkt(args, fh, output):
    # procedure to get start and end times from PACE S-band HKT data file name
    print("Running l0info_hkt (version: %s) \n" % __version__)
    
    if args.verbose:
        print("Reading PACE S-band data file.")

    status = 104
    try:
        # read the CFE header and advance the reading pointer if it exists
        headerInfo = readFileHeader(fh)

        # if not a PACE file, raise exception
        if (headerInfo["SCID"] != b"PACE"):
            print('Not a valid PACE L0 file')
            status = 102 # invalid instrument 
            raise Exception

        # get all the packets
        all_packets = ccsdspy.utils.split_packet_bytes(fh, include_primary_header=True)

        # get the start and stop times 
        stime = get_timestamp(all_packets, True)
        etime = get_timestamp(all_packets, False)

        print("start_time=%s" % stime)
        print("stop_time=%s" % etime)
        
        if output:
            output.write("start_time=%s\n" % stime)
            output.write("stop_time=%s\n" % etime)
     
    except Exception as e:
        return status


    # status 0 if no issues
    return 0
