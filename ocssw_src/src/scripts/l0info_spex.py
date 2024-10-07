#!/usr/bin/env python3

# procedure to get start and end times from Level-0 fpacket files for SPEXone
#  The times are extracted from the first fpacket of each image fpacket group.
#  Reference: SPEXone TMTC Database, SPX1-TN-004

__version__ = '1.3.0_2023-01-23'

from l0info_utils import read_packet, read_apid, get_anc_packet_time, is_bad_packet

def is_bad_apid(apid,fh):
    if apid == -1:
        return False
    if (apid < 800) or (apid > 849):
        print("Invalid packet header at byte %d."%fh.tell())
        return True
    else:
        return False

def l0info_spex(args, fh, output):
    print("Running l0info_spex (version: %s) \n" % __version__)
       
    status = 0
    if args.verbose:
        print("Reading SPEXone science data.")

    # Skip CFE header
    fh.seek(64)
    bSPW = True

    # Find first science fpacket
    apid = 0
    packetLength = 8
    firstp = 0
    while ((apid != 848) or (not firstp)) and (packetLength > 7):
        fpacket, packetLength = read_packet(fh, bSPW)
        apid = read_apid(fpacket)
        if is_bad_apid(apid,fh):  return 104
        if fpacket: firstp = (fpacket[2] & 64)>>6
    
    if is_bad_packet(packetLength,fh):  return 104
        
    # If no science packets, rewind and get times from HKT packets
    if (packetLength==0):
        print("No science packets found in file.")
        if output:
            output.write("datatype=SPEX\n")
        status = 110
        fh.seek(64)
        apid = 0
        packetLength = 8
        csec = 0
        while ((apid != 800) or (csec<1000000000)) and (packetLength > 7):
            fpacket, packetLength = read_packet(fh, bSPW)
            apid = read_apid(fpacket)
            csec = int.from_bytes(fpacket[6:10],'big') 
            
        if is_bad_packet(packetLength,fh):  return 104
        
        if fpacket:
            mpacket = fpacket
            try:
                stime = get_anc_packet_time(mpacket[6:12],0,'msec')
                print("start_time=%s" % stime.strftime('%Y-%m-%dT%H:%M:%S.%f'))
                if output:
                    output.write("start_time=%s\n" % stime.strftime('%Y-%m-%dT%H:%M:%S.%f'))
            except:
                return 120
                
            while packetLength > 7:
                fpacket, packetLength = read_packet(fh, bSPW)
                apid = read_apid(fpacket)
                if (apid == 800) and fpacket: mpacket = fpacket
            
            if is_bad_packet(packetLength,fh):  return 104
            
            etime = get_anc_packet_time(mpacket[6:12],0,'msec')
    else:
        # Determine data type from packet fields
        fmode = fpacket[125] # frame mode
        omode = fpacket[130] # output mode
        imrlen = int.from_bytes(fpacket[296:300],'big') # swap_endian(long(packet[296:299],0)); image record length
        dtype = -1
        if (fmode == 2) and (omode == 1) and (imrlen > 0): dtype = 1  # Science
        if (fmode == 1) and (omode == 3) and (imrlen > 0): dtype = 2  # Diagnostic

        if dtype == 1:
            # Science data
            print("datatype=SPEX")
            if output:
                output.write("datatype=SPEX_SCI\n")
        elif dtype == 2:
            # Diagnostic data
            print("datatype=SPEX_DIAG")
            if output:
                output.write("datatype=SPEX_DIAG\n")
        else:
            print("Invalid SPEXone data type")
            status = 101
            
        # Get start time
        try:
            stime = get_anc_packet_time(fpacket[300:306],0,'msec')
            print("start_time=%s" % stime.strftime('%Y-%m-%dT%H:%M:%S.%f'))
            if output:
                output.write("start_time=%s\n" % stime.strftime('%Y-%m-%dT%H:%M:%S.%f'))
        except:
            return 120

        # Read to end of file
        epacket = fpacket
        while packetLength > 0:
            apid = 0
            firstp = 0
            while ((apid != 848) or (not firstp)) and (packetLength > 7):
                fpacket, packetLength = read_packet(fh, bSPW)
                apid = read_apid(fpacket)
                if is_bad_apid(apid,fh):  return 104
                if fpacket: firstp = (fpacket[2] & 64)>>6
            
            if is_bad_packet(packetLength,fh):  return 104    
                                            
            if fpacket: epacket = fpacket
        
        etime = get_anc_packet_time(epacket[300:306],0,'msec')

    print("stop_time=%s" % etime.strftime('%Y-%m-%dT%H:%M:%S.%f'))
    if output:
        output.write("stop_time=%s\n" % etime.strftime('%Y-%m-%dT%H:%M:%S.%f'))

    return status
