#!/usr/bin/env python3


__version__ = '2.8.0_2023-08-09'

__dtypes__ = ['','','_DARK','_SOL','_SPCA','_LIN','_LUN','_DIAG','_STAT',
    '_SPEC','','_SNAP-X','_SNAP-I','','_LUN-ST','_RAW']

import numpy as np  
from l0info_utils import read_packet, read_apid, get_anc_packet_time, is_bad_packet  



CRC_TABLE = [0x0000, 0x1021, 0x2042, 0x3063, 0x4084, 0x50a5,
             0x60c6, 0x70e7, 0x8108, 0x9129, 0xa14a, 0xb16b,
             0xc18c, 0xd1ad, 0xe1ce, 0xf1ef, 0x1231, 0x0210,
             0x3273, 0x2252, 0x52b5, 0x4294, 0x72f7, 0x62d6,
             0x9339, 0x8318, 0xb37b, 0xa35a, 0xd3bd, 0xc39c,
             0xf3ff, 0xe3de, 0x2462, 0x3443, 0x0420, 0x1401,
             0x64e6, 0x74c7, 0x44a4, 0x5485, 0xa56a, 0xb54b,
             0x8528, 0x9509, 0xe5ee, 0xf5cf, 0xc5ac, 0xd58d,
             0x3653, 0x2672, 0x1611, 0x0630, 0x76d7, 0x66f6,
             0x5695, 0x46b4, 0xb75b, 0xa77a, 0x9719, 0x8738,
             0xf7df, 0xe7fe, 0xd79d, 0xc7bc, 0x48c4, 0x58e5,
             0x6886, 0x78a7, 0x0840, 0x1861, 0x2802, 0x3823,
             0xc9cc, 0xd9ed, 0xe98e, 0xf9af, 0x8948, 0x9969,
             0xa90a, 0xb92b, 0x5af5, 0x4ad4, 0x7ab7, 0x6a96,
             0x1a71, 0x0a50, 0x3a33, 0x2a12, 0xdbfd, 0xcbdc,
             0xfbbf, 0xeb9e, 0x9b79, 0x8b58, 0xbb3b, 0xab1a,
             0x6ca6, 0x7c87, 0x4ce4, 0x5cc5, 0x2c22, 0x3c03,
             0x0c60, 0x1c41, 0xedae, 0xfd8f, 0xcdec, 0xddcd,
             0xad2a, 0xbd0b, 0x8d68, 0x9d49, 0x7e97, 0x6eb6,
             0x5ed5, 0x4ef4, 0x3e13, 0x2e32, 0x1e51, 0x0e70,
             0xff9f, 0xefbe, 0xdfdd, 0xcffc, 0xbf1b, 0xaf3a,
             0x9f59, 0x8f78, 0x9188, 0x81a9, 0xb1ca, 0xa1eb,
             0xd10c, 0xc12d, 0xf14e, 0xe16f, 0x1080, 0x00a1,
             0x30c2, 0x20e3, 0x5004, 0x4025, 0x7046, 0x6067,
             0x83b9, 0x9398, 0xa3fb, 0xb3da, 0xc33d, 0xd31c,
             0xe37f, 0xf35e, 0x02b1, 0x1290, 0x22f3, 0x32d2,
             0x4235, 0x5214, 0x6277, 0x7256, 0xb5ea, 0xa5cb,
             0x95a8, 0x8589, 0xf56e, 0xe54f, 0xd52c, 0xc50d,
             0x34e2, 0x24c3, 0x14a0, 0x0481, 0x7466, 0x6447,
             0x5424, 0x4405, 0xa7db, 0xb7fa, 0x8799, 0x97b8,
             0xe75f, 0xf77e, 0xc71d, 0xd73c, 0x26d3, 0x36f2,
             0x0691, 0x16b0, 0x6657, 0x7676, 0x4615, 0x5634,
             0xd94c, 0xc96d, 0xf90e, 0xe92f, 0x99c8, 0x89e9,
             0xb98a, 0xa9ab, 0x5844, 0x4865, 0x7806, 0x6827,
             0x18c0, 0x08e1, 0x3882, 0x28a3, 0xcb7d, 0xdb5c,
             0xeb3f, 0xfb1e, 0x8bf9, 0x9bd8, 0xabbb, 0xbb9a,
             0x4a75, 0x5a54, 0x6a37, 0x7a16, 0x0af1, 0x1ad0,
             0x2ab3, 0x3a92, 0xfd2e, 0xed0f, 0xdd6c, 0xcd4d,
             0xbdaa, 0xad8b, 0x9de8, 0x8dc9, 0x7c26, 0x6c07,
             0x5c64, 0x4c45, 0x3ca2, 0x2c83, 0x1ce0, 0x0cc1,
             0xef1f, 0xff3e, 0xcf5d, 0xdf7c, 0xaf9b, 0xbfba,
             0x8fd9, 0x9ff8, 0x6e17, 0x7e36, 0x4e55, 0x5e74,
             0x2e93, 0x3eb2, 0x0ed1, 0x1ef0
             ]
CRC_LENGTH: int = 4

science_oci_format_header = 0x2bc
shift_bit = 0xffff


def compute_CRC_checksum(data: list, length: int):
    # byte_0 = 0x34  - Logical Address
    # byte_1 = 0x02  - protocol ID
    CRC_BUFFER: list[int] = [ 0x34,0x02] # placing these two bytes in the beggining (SpaceWire Header)
    crcA = 0xFFFF
    crcB = 0xFFFF
    CRC_BUFFER += data[:length]
    size = (length + 2) // 4 
    for i in range(size):
        crcA = CRC_TABLE[((crcA >> 8) ^ CRC_BUFFER[4 * i + 0])
                         & 0xFF] ^ (crcA << 8)
        crcA = crcA & shift_bit
        crcA = CRC_TABLE[((crcA >> 8) ^ CRC_BUFFER[4 * i + 1])
                         & 0xFF] ^ (crcA << 8)
        crcA = crcA & shift_bit
        crcB = CRC_TABLE[((crcB >> 8) ^ CRC_BUFFER[4 * i + 2])
                         & 0xFF] ^ (crcB << 8)
        crcB = crcB & shift_bit
        crcB = CRC_TABLE[((crcB >> 8) ^ CRC_BUFFER[4 * i + 3])
                         & 0xFF] ^ (crcB << 8)
        crcB = crcB & shift_bit
    return (crcA << 16) | crcB



def get_band_dims(apacket):

# extract aggregation information and compute data dimensions from ancillary packet

#       Arguments
#
#       Name    Type            I/O      Description
#       ----    ----            ---      -----------
#    apacket    byte(*)         I    Ancillary packet array
#    ncp    int         O    Number of CCD band pixels
#    nbb    int         O    Number of blue CCD bands
#    nrb    int         O    Number of red CCD bands
#    nsp    int         O    Number of SWIR band pixels
#    ndc    int         O    Number of CCD dark collect pixels
#    nds    int         O    Number of SWIR dark collect pixels
#    btaps    int         O    Blue CCD tap enable flags
#    rtaps    int         O    Red CCD tap enable flags
#    itable  struct         O    Spatial aggregation table        

#    orig: Frederick S. Patt, SAIC, 1 August 2018

    # Create spatial aggregation table structure
    itable = np.zeros(10, dtype={'names':('dtype', 'iagg', 'lines'),'formats':('i4','i4','i4')})

    # Extract spatial aggregation table and compute numbers of pixels
    ioff = 36
    nagg = [1,2,4,8]
    ncp = 0
    nsp = 0
    ndc = 0
    nds = 0
    for i in range(0,10):
        itable['dtype'][i] = apacket[ioff+3]%16
        itable['iagg'][i] = apacket[ioff+2]%4
        itable['lines'][i] = apacket[ioff]*256 + apacket[ioff+1]
        ioff = ioff + 4
        if ((itable['dtype'][i] > 0) and (itable['dtype'][i] <= 12) and (itable['dtype'][i] != 10)):
            if (itable['dtype'][i] == 2):
                ndc = ndc + itable['lines'][i]/nagg[itable['iagg'][i]] 
                nds = nds + itable['lines'][i]/8
            ncp = ncp + itable['lines'][i]/nagg[itable['iagg'][i]]
            nsp = nsp + itable['lines'][i]/8
    if (ncp == 0):
        ncp = 1
        nsp = 1

    if (ndc == 0):
        ndc = 1
        nds = 1
  
    ioff = ioff + 4

    # Extract spectral aggregation and compute numbers of bands
    #  Tap enable flags
    btap = apacket[ioff+2]*256 + apacket[ioff+3] 
    rtap = apacket[ioff]*256 + apacket[ioff+1]
    btaps = np.zeros(16)
    rtaps = np.zeros(16)
    
    #  Tap aggregation factors
    bagg = int.from_bytes(apacket[ioff+8:ioff+12],'big')
    ragg = int.from_bytes(apacket[ioff+4:ioff+8],'big')
    
    baggs = np.zeros(16)
    raggs = np.zeros(16)
    
    #  Compute number of bands for enabled taps
    nbb = 0
    nrb = 0
    ken = 1
    kag = 3
    lag = 1
    
    for i in range(15,-1,-1):
        btaps[i] = np.bitwise_and(btap, ken)/ken
        if (btaps[i]):
            baggs[i] = nagg[int(np.bitwise_and(bagg, kag)/lag)]
            nbb = nbb + 32/baggs[i]
        
        rtaps[i] = np.bitwise_and(rtap, ken)/ken
        if (rtaps[i]):
            raggs[i] = nagg[int(np.bitwise_and(ragg, kag)/lag)]
            nrb = nrb + 32/raggs[i]
        
        ken = ken*2
        kag = kag*4
        lag = lag*4
        
    return ncp,nbb,nrb,nsp,ndc,nds,btaps,rtaps,itable,baggs,raggs

def anc_compare(apacket0,apacket):

# function to compare spatial and spectral data collection parameters
#  from two OCI ancillary data packets

#        Arguments
# 
#        Name        Type    I/O      Description
#        ----        ----    ---      -----------
#     apacket0(*)    byte     I    Ancillary packet from first spin
#     apacket(*)    byte     I    Ancillary packet from next spin

#  Returns 1 (TRUE) if packets agree, 0 if not
    
    anc_compare = 1
    strmsg = ""
    # Compare spatial data collection fields
    ioff = 36
    ilen = 40
    if apacket[ioff:ioff+ilen] != apacket0[ioff:ioff+ilen]:
        c_time = get_anc_packet_time(apacket)
        strmsg += "\nSpatial table change at %s" % c_time.strftime('%Y-%m-%dT%H:%M:%S.%f')
        anc_compare = 0

    # Compare spectral data collection fields
    joff = 80
    jlen = 12
    if apacket[joff:joff+jlen] != apacket0[joff:joff+jlen]:
        c_time = get_anc_packet_time(apacket)
        strmsg += "\nSpectral table change at %s" % c_time.strftime('%Y-%m-%dT%H:%M:%S.%f')
        anc_compare = 0

    return anc_compare, strmsg

def is_bad_apid(apid,fh):
    if apid == -1:
        return False
    if (apid < 550) or (apid > 749):
        print("Invalid packet header at byte %d."%fh.tell())
        return True
    else:
        return False

def get_oci_data_type(fpacket):
    # To determine the OCI data type from and ancillary packet
    # Returns: dtype (1 - 12) if valid ancillary packet, -1 otherwise

    # Check APID
    apid = read_apid(fpacket)
    if (apid != 636): return -1

    dtypes = np.zeros(10)-1
    ioff = 36
    for i in range(0,10):
        dtypes[i] = fpacket[ioff+3] % 16
        ioff = ioff + 4
    kd = np.argwhere((dtypes != 2) & (dtypes != 10)) # Exclude dark and no processing types
    dtype = np.max(dtypes[kd])
    
    return int(dtype)


def l0info_oci(args, fh, output, bDSB, crcSumCheck):
    # procedure to get start and end times from Level-0 packet files for OCI
    print("Running l0info_oci (version: %s) \n" % __version__)
    
    dtype = -1
    str_stime = ''
    str_etime = ''
    
    status = 0
    if args.verbose:
        print("Reading OCI science data file.")
        
    bSPW = False
    lpoint = 0
    if bDSB:
        # chead = fh.read(64)
        lpoint = 64
        bSPW = True
    fh.seek(lpoint)
    
    # Get first ancillary packet
    apid = 0
    packetLength = 8
    if crcSumCheck:
        print("CRC checksum will be perfomed for each science packet")
    while (apid != 636) and (packetLength > 7):
        try:
            fpacket, packetLength = read_packet(fh, bSPW)
            apid = read_apid(fpacket)
            if is_bad_apid(apid,fh):  return 104

            # check CRC sum
            if fpacket and crcSumCheck:
                phead_format = int.from_bytes(fpacket[:2], "big")
                if  phead_format == science_oci_format_header:
                    dataCompute = list(fpacket[0:-CRC_LENGTH])
                    crcCalcs = compute_CRC_checksum(dataCompute, len(dataCompute))
                    crcExpected = int.from_bytes(fpacket[-CRC_LENGTH:], "big")
                    if crcExpected != crcCalcs:
                        print(f"Error : CRC checksum failed! Expected {hex(crcExpected)}, calculated {hex(crcCalcs)} for packet with apid {apid}")
                        return 131
        except:            
            return 103
    if is_bad_packet(packetLength,fh):  return 104
    if fpacket:
        print("CRC check finished.")
        
    # If no ancillary packets, rewind and get times from HKT packets
    if fpacket is None: 
        print('No stored science packets found in file.')
        status = 110
        
        fh.seek(lpoint)
        apid = 0
        packetLength = 8
        ccsec = 0
        firstraw = 0
        # while ((apid < 550) or (apid > 750) or (apid == 558)) and (packetLength > 8):
        while (ccsec < 1900000000) and (packetLength>7):
            fpacket, packetLength = read_packet(fh, bSPW)
            if fpacket:
                apid = read_apid(fpacket)
            else:
                continue
            ccsec = int.from_bytes(fpacket[6:10],'big')
        # if is_bad_packet(packetLength,fh):  return 104       
        if fpacket:
            mpacket = fpacket
            try:
                stime = get_anc_packet_time(mpacket[6:12],0,'msec')
                str_stime = stime.strftime('%Y-%m-%dT%H:%M:%S.%f')
                # print("start_time=%s" % str_stime)
            except:
                return 120

            while fpacket and (packetLength>7):
                fpacket, packetLength = read_packet(fh, bSPW)
                if fpacket:
                    apid = read_apid(fpacket)
                else:
                    continue
                    
                # Check for raw mode
                if (apid == 705) and (not firstraw):
                    dtype = 15  # '_RAW'
                    firstraw = 1
                
                if (apid >= 550) and (apid < 750) and (apid!=705):  mpacket = fpacket
            
            if is_bad_packet(packetLength,fh):  return 104
            print("start_time=%s" % str_stime)
            etime = get_anc_packet_time(mpacket[6:12],0,'msec')
        else:
            if output:
                output.write("datatype=OCI\n")
            print("No packets with valid time tags in file.")
            return status
    
    else:
        # Get data type
        dtype = get_oci_data_type(fpacket)
        if (dtype<1) or (dtype>12):
            print("Invalid data type %d."%dtype)
            return 104
                    
        # Get start time
        try:
            stime = get_anc_packet_time(fpacket,28,'usec')
            str_stime = stime.strftime('%Y-%m-%dT%H:%M:%S.%f')
            print("start_time=%s" % str_stime)
            epacket = fpacket
            etime = get_anc_packet_time(epacket,28,'usec')
        except:
            return 120        
        
        # Read to end of file
        apacket = fpacket
        acomp = 0
        nagg = np.array([1,2,4,8])
        spnp = int.from_bytes(apacket[24:28],'big')
        
        while fpacket:
            etime = get_anc_packet_time(apacket,28,'usec')
            if args.verbose:
                # spnum = swap_endian(long(apacket(24:27),0))
                spnum = int.from_bytes(apacket[24:28],'big')
                tmpstr = ""
                if ((spnum-spnp) > 10):
                    tmpstr += "\nSpin number gap: %d, %d" % (spnp, spnum)
                spnp = spnum
                if (not acomp):
                    tmpstr += "\n\nInstrument configuration change at spin %d, %s" % (spnum,etime)
                    ncp,nbb,nrb,nsp,ndc,nds,btaps,rtaps,itable,baggs,raggs = get_band_dims(apacket)
                    iz = np.where((itable['dtype'] != 0) & (itable['dtype'] != 10))[0]
                    tmpstr += "\nData Type       %s" % str(itable['dtype'][iz])
                    iagg = nagg[itable['iagg'][iz]]
                    tmpstr += "\nNumber of pixels %s" % str(itable['lines'][iz]/iagg)
                    tmpstr += "\nAggregation     %s" % str(iagg)
                    tmpstr += "\nBlue bands: %d   Red bands: %d" %(nbb, nrb)
                    tmpstr += "\nBlue aggregation per tap: %s" % str(baggs)
                    tmpstr += "\nRed aggregation per tap: %s\n" % str(raggs)
                    print(tmpstr)
                    
                    apacket = fpacket
                    
            apid = 0
            while (apid != 636) and (packetLength>0):
                fpacket, packetLength = read_packet(fh, bSPW)     
                apid = read_apid(fpacket)
                if is_bad_apid(apid,fh):  return 104
            
            if is_bad_packet(packetLength,fh):  return 104
            if fpacket:  epacket = fpacket
        
        try:
            etime = get_anc_packet_time(epacket,28,'usec')
        except:
            return 120
    
    str_etime = etime.strftime('%Y-%m-%dT%H:%M:%S.%f')
    print("stop_time=%s" % str_etime)
    
    datatype_name = ''
    if dtype > 0:
        try:
            datatype_name = 'OCI' + __dtypes__[dtype]
            if datatype_name=="OCI_LUN-ST":
                print("Unexpected datatype [%s] is detected." % datatype_name)
                status = 130
        except:
            print("OCI data type name reading error.")
            status = 130
    
    if datatype_name=='':
        datatype_name = 'OCI'
    
    print("datatype=%s" % datatype_name)
    if output:
        output.write("datatype=%s\n" % datatype_name)
        if str_stime!='':
            output.write("start_time=%s\n" % str_stime)
        if str_etime!='':
            output.write("stop_time=%s\n" % str_etime)
    
    return status
