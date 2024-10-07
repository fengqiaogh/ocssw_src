from datetime import datetime, timezone, timedelta
from seadasutils.time_utils import get_leap_seconds
import julian

__version__ = '1.3.0_2022-12-23'

# seconds between 1970-01-01 and 1958-01-01
epoch_offset = (datetime(1958,1,1,0,0,0,tzinfo=timezone.utc) - datetime(1970,1,1,0,0,0,tzinfo=timezone.utc)).total_seconds()

def print_utils_version():
    print("Running l0info_utils (version: %s) \n" % __version__)

def read_cfe_header(filehandle, bVerbose=False):
    # read header from a DSB file and print contents
    # Reference PACE-SYS-ICD_0010, :Observatory to Ground Data Interface Requirements/Control Document: 
    # Volume II  S-band CFDP Specification, section 3.5.1
    
    chead = filehandle.read(64)
    
    ftime = int.from_bytes(chead[24:28],'big') + int.from_bytes(chead[28:32],'big')/2**32
    ctime = (datetime.utcfromtimestamp(ftime+epoch_offset) - timedelta(seconds=get_leap_seconds(ftime))).isoformat(sep='T', timespec='microseconds')
    if bVerbose: print('File creation date/time: %s' % ctime)
    
    # Check for valid cFE header
    dtype = -1  
    scid =  "".join([chr(x) for x in chead[12:16]])
    if scid != 'PACE':
        print('Not a valid PACE L0 file')
        return ctime, dtype

    # Check for science data file and determine instrument
    subtype = "".join([chr(x) for x in chead[4:8]])
    if subtype == ' SCI':
        instid = "".join([chr(x) for x in chead[20:24]])
        instids = ['',' OCI','HARP','SPEX']
        dtype = instids.index(instid)
    else:
    # Check for S-band file subtype
        subtype = int(chead[7])
        if subtype == 101:   dtype = 0
    
    return ctime, dtype
    
def read_dsb_header(filehandle, bVerbose=False):
    # Routine to read header from a DSB file and print contents
    # Reference PACE-SYS-ICD_0010, :Observatory to Ground Data Interface Requirements/Control Document: 
    # Volume II  S-band CFDP Specification, section 3.5.2

    dhead = filehandle.read(76)

    ftime = int.from_bytes(dhead[0:4],'big') + int.from_bytes(dhead[4:8],'big')/2**32
    etime = (datetime.utcfromtimestamp(ftime+epoch_offset) - timedelta(seconds=get_leap_seconds(ftime))).isoformat(sep='T', timespec='microseconds')
    if bVerbose: 
        print('File close date/time:    %s' % etime)
        print('File path/name:         %s'% "".join([chr(x) for x in dhead[12:44]]))

    return etime

def read_packet(filehandle, bSPW = True):
    # procedure to read a single CCSDS packet from a file.
    # bSPW to read spacewire header

    # Read spacewire header if needed
    if bSPW:
      spwh = filehandle.read(2)
      if not spwh:    # EOF
          return None, 0

    # Read packet header
    phead = filehandle.read(6)
    if not phead:    # EOF
        return None, 0
    elif len(phead)<6:
        return None, 0

    # Get length of packet body
    packetLength = phead[4]*256 + phead[5] + 1
    try:
        pbod = filehandle.read(packetLength)
        if not pbod:
            return None, 0
    except:
        print('Packet reading error.')
        return None, 0
    
    # Concatenate header and body

    packet = phead+pbod
    packetLength = packetLength + 6

    return packet, packetLength
    
def read_apid(fpacket):
    if fpacket:
        return fpacket[1] + 256*(fpacket[0] % 8)
    else:
        return -1
    
def get_anc_packet_time(apacket,toff,str_subsec):
    # procedure to unpack and convert the CCSDS segmented time code from the OCI ancillary packet

    #  Get day count since Jan. 1, 1958 (Julian day 2436205)
    # toff = 28
    sec58 = int.from_bytes(apacket[toff:toff+4],'big') 
    sec58 -= get_leap_seconds(sec58)
    
    day58 = int(sec58/86400)
    sec58 -= day58*86400
    
    usec = 0
    msec = 0
    
    if str_subsec == 'usec':
        # Get microseconds
        usec = int.from_bytes(apacket[toff+4:toff+8],'big')
        usec = usec/4096.0/1000000. #20 MSBs used, last 12 are spares
    elif str_subsec == 'msec':
        # Get milliseconds
        msec = int.from_bytes(apacket[toff+4:toff+6],'big')
        msec = msec/65536.0

    dt = julian.from_jd(day58 + 2436205.0) - timedelta(hours=12)
    dt += timedelta(seconds= sec58 + msec + usec)
    return dt
    
def is_bad_packet(packetLength,fh):
    if (packetLength>0) and (packetLength<8):
        print(f"packetLength={packetLength}")  # debug LH
        print("Invalid packet header at byte %d."%fh.tell())
        return True
    else:
        return False
