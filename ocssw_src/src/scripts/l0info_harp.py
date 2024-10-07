#!/usr/bin/env python3

# procedure to get start and end times from Level-0 fpacket files for HARP2
#  The times are extracted from the detector image headers.  The detector number 
#  (1, 2, 3 or All) is also determined
#  Reference: HARP2 Telemetry Manual

__version__ = '1.5.6_2024-03-06'

import numpy as np
import os
import copy
from datetime import datetime, timedelta
from l0info_utils import read_packet, read_apid, get_anc_packet_time, is_bad_packet

# str_imhead = "HARP" # [72, 65, 82, 80]    # Header start mark ASCII "HARP", \x48\x41\x52\x50
b_imhead = b'HARP'

def is_fill_image(fpacket,packetLength,fh,firstfill):
    if not np.any(np.frombuffer(fpacket[22:packetLength-1],dtype=np.uint8)-255):
        # Routine to output location of first HARP2 image fill data
        if firstfill:
            firstfill = False
            print("Image fill data at byte %d."%fh.tell())                    
        return True, firstfill
    else:
        return False, firstfill

def is_bad_apid(apid,fh):
    if apid == -1:
        return False
    if (apid < 750) or (apid > 799):
        print("Invalid packet or header at byte %d."%fh.tell())
        return True
    else:
        return False

def l0info_harp(args, fh, output):
    print("Running l0info_harp (version: %s) \n" % __version__)
    
    status = 0
    if args.verbose:
        print("Reading HARP2 science data file.")
        
    file_size = os.fstat(fh.fileno()).st_size
        
    # Skip CFE header
    fh.seek(64)
    bSPW = True
    
    # Find first science fpacket (APID 751)
    apid = 0
    packetLength = 8
    firstfill = True
    
    
    while (apid != 751) and (packetLength > 7): 
        fpacket, packetLength = read_packet(fh, bSPW)
        apid = read_apid(fpacket)
        if is_bad_apid(apid,fh):  return 104
    
    if packetLength > 7:
        if is_bad_packet(packetLength,fh):  return 104
        isFillImage, firstfill = is_fill_image(fpacket,packetLength,fh,firstfill) 
        
        # Get detector number from AP header
        detnum = (fpacket[12] & 48)>>4
    
        # Find an image header
        imhp = -1
        while (imhp == -1) and (fh.tell() < file_size):
            imhp = fpacket.find(b_imhead)
                    
            if imhp == -1: # If no header in current fpacket, find next science fpacket
                apid = 0
                packetLength = 8
                while (apid != 751) and (packetLength > 7):
                    fpacket, packetLength = read_packet(fh, bSPW)
                    apid = read_apid(fpacket)
                    if is_bad_apid(apid,fh):  return 101
                
                if packetLength > 7:
                    if is_bad_packet(packetLength,fh):  return 104
                    isFillImage, firstfill = is_fill_image(fpacket,packetLength,fh,firstfill) 
                    
                    # Check for more than one detector in file
                    detn = (fpacket[12] & 48)>>4
                    if detn != detnum: detnum = 999
                    
        if fh.tell() >= file_size:
            print('No image headers found in file.')                
      
        if args.verbose:
            pctr = (fpacket[2] % 64)*256 + fpacket[3]
            print(f"Image pointer {imhp}, detector {detnum} in packet {pctr}")
    else:
        print("No stored science packets found in file.")
        
    # If no stored science packets, check for real-time
    if fh.tell() >= file_size:
        status = 110
        etime = datetime(2000,1,1)
        stime = datetime(3000,1,1)
        
        fh.seek(64)
        apid = 0
        packetLength = 8
        while (fh.tell() < file_size) and (apid != 757) and (packetLength > 7):
            fpacket, packetLength = read_packet(fh, bSPW)
            apid = read_apid(fpacket)
            if is_bad_apid(apid,fh):  return 101
        
        if is_bad_packet(packetLength,fh):  return 104
        
        if fpacket:
            mpacket = fpacket
            isRealTime = False
            try:
                ctime = get_anc_packet_time(fpacket[19:25],0,'msec')
                if ctime<stime:
                    stime = ctime
            except Exception as e:
                print(e)
                return 120
            
            while fpacket and (packetLength > 7):
                fpacket, packetLength = read_packet(fh, bSPW)
                apid = read_apid(fpacket)
                if is_bad_apid(apid,fh):  return 101
                
                if apid == 757: 
                    mpacket = fpacket
                    isRealTime = True
                ctime = get_anc_packet_time(mpacket[19:25],0,'msec')
                # print(f"ctime={ctime.strftime('%Y-%m-%dT%H:%M:%S.%f')}")
                if ctime<stime:
                    stime = ctime
                if ctime>etime:
                    etime = ctime
            
            if is_bad_packet(packetLength,fh):  return 104

            if isRealTime:
                try:
                    print("detector=REAL\n")
                    print("start_time=%s" % stime.strftime('%Y-%m-%dT%H:%M:%S.%f'))                
                    print("stop_time=%s" % etime.strftime('%Y-%m-%dT%H:%M:%S.%f'))
                    if output:
                        output.write("start_time=%s\n" % stime.strftime('%Y-%m-%dT%H:%M:%S.%f'))
                        output.write("stop_time=%s\n" % etime.strftime('%Y-%m-%dT%H:%M:%S.%f'))
                        output.write("detector=REAL\n")
                except Exception as e:
                    print(e)
                    return 120
                #print(f"status={status}")    
                return status

    # If no science packets, rewind and get times from HKT packets (APID 750)
    if fh.tell() >= file_size:
        status = 110
        etime = datetime(2000,1,1)
        stime = datetime(3000,1,1)
        
        fh.seek(64)
        apid = 0
        packetLength = 8
        ccsec = 0
        while ((fh.tell() < file_size) and (apid != 750) and (packetLength > 7)) or (ccsec < 2000000000):
            fpacket, packetLength = read_packet(fh, bSPW)
            apid = read_apid(fpacket)
            if is_bad_apid(apid,fh):  return 101
            
            # Check for invalid time stamp
            ccsec = int.from_bytes(fpacket[6:10],'big') # swap_endian(long(packet[6:9],0))
            
        if is_bad_packet(packetLength,fh):  return 104
        
        if fpacket:
            mpacket = fpacket
            try:
                ctime = get_anc_packet_time(mpacket[6:12],0,'msec')
                if ctime<stime:
                    stime = ctime
            except Exception as e:
                print(e)
                return 120
            
            while fpacket and (packetLength > 7):
                fpacket, packetLength = read_packet(fh, bSPW)
                apid = read_apid(fpacket)
                if apid == 750:  mpacket = fpacket
                ctime = get_anc_packet_time(mpacket[6:12],0,'msec')
                # print(f"ctime={ctime.strftime('%Y-%m-%dT%H:%M:%S.%f')}")
                if ctime<stime:
                    stime = ctime
                if ctime>etime:
                    etime = ctime
            
            if is_bad_packet(packetLength,fh):  return 104
            
        try:
            print("start_time=%s" % stime.strftime('%Y-%m-%dT%H:%M:%S.%f'))
            if output:
                output.write("start_time=%s\n" % stime.strftime('%Y-%m-%dT%H:%M:%S.%f'))
            print("stop_time=%s" % etime.strftime('%Y-%m-%dT%H:%M:%S.%f'))
            if output:
                output.write("stop_time=%s\n" % etime.strftime('%Y-%m-%dT%H:%M:%S.%f'))
        except Exception as e:
            print(e)
            return 120
        #print(f"status={status}")    
        return status
    
    # Get start time
    try:
        stime = get_anc_packet_time(fpacket[imhp+7:imhp+13],0,'msec')
        detnum = (fpacket[12] & 48)>>4
    except:
        return 120

    # Read to end of file
    # Look for next image header
    packetLength = 8
    last_pct = 0
    etime = stime
    
    # when updating end time based on new image header's time,
    # new - curr should not exceed 1 day. Average time
    # diff between images is ~ 0:00:00.221710. You can lower the 
    # threshold more if needed
    MAX_TIME_DIFF = timedelta(days=1.0) 

    while fh.tell() < file_size:
        apid = 0

        # Find next science fpacket
        while (fh.tell() < file_size) and (apid != 751) and (packetLength > 7):
            fpacket, packetLength = read_packet(fh, bSPW)
            apid = read_apid(fpacket)
            if is_bad_apid(apid,fh):  return 104
        
        if fpacket:
            if is_bad_packet(packetLength,fh):  return 104
            isFillImage, firstfill = is_fill_image(fpacket,packetLength,fh,firstfill) 
            
            # Look for image header
            imh = -1
            while (imh == -1) and fh.tell() < file_size:
                if args.verbose:
                    imh_pct = int(fh.tell()/file_size*100)
                    if imh_pct - last_pct > 9:
                        print("%d%% of packets read." % imh_pct)
                        last_pct = imh_pct
                imh = fpacket.find(b_imhead)
                                
                if (imh == -1): # If no header in current fpacket, find next science fpacket
                    apid = 0
                    while (fh.tell() < file_size) and (apid != 751) and (packetLength > 7):
                        fpacket, packetLength = read_packet(fh, bSPW)
                        apid = read_apid(fpacket)
                        if is_bad_apid(apid,fh):  return 104
                    if (fh.tell() < file_size):
                        if is_bad_packet(packetLength,fh):  return 104
                        isFillImage, firstfill = is_fill_image(fpacket,packetLength,fh,firstfill)
                    
                    # if is_bad_image(fpacket,packetLength,fh):  return 111                    
            if fpacket and (imh > -1):  # fh.tell() < file_size:
                # Check for more than one detector in file
                detn = (fpacket[12] & 48)>>4
                if detn != detnum: detnum = 999
                
                imhp = imh
                if args.verbose:
                    pctr = (fpacket[2] % 64)*256 + fpacket[3]
                    print(f"Image pointer {imhp}, detector {detn} in packet {pctr}")
                    
                if imhp < 8202:
                    ctime = get_anc_packet_time(fpacket[imhp+7:imhp+13],0,'msec')
                    # if ctime<stime:
                    #     stime = ctime
                    if ctime>etime and (ctime - etime) < MAX_TIME_DIFF:
                        etime = ctime
                   
            

    if args.verbose:
        print("100% of packets read.")
    
    try:
        if imhp <0:
            status = 110
        print("start_time=%s" % stime.strftime('%Y-%m-%dT%H:%M:%S.%f'))
        print("stop_time=%s" % etime.strftime('%Y-%m-%dT%H:%M:%S.%f'))
        if output:
            output.write("start_time=%s\n" % stime.strftime('%Y-%m-%dT%H:%M:%S.%f')) 
            output.write("stop_time=%s\n" % etime.strftime('%Y-%m-%dT%H:%M:%S.%f'))
    except Exception as e:
        print(e)
        return 120
    
    if detnum == 999:
        print("detector=all")
        if output:
            output.write("detector=all\n")
    else:
        print("detector=%d"%detnum)
        if output:
            output.write("detector=%d\n"%detnum)
    
    return status

