import sys
import numpy as np
import datetime

TAI58_OFFSET = datetime.datetime(1970, 1, 1) - datetime.datetime(1958, 1, 1)
LEAPSEC = 37  # TODO: look up by year

def decode_timestamp(seconds,subseconds):
    return float(seconds) + float(subseconds) / pow(2,8*subseconds.itemsize)

CCSDS_timestamp = np.dtype([
    ('seconds',   '>u4'),   # Time since start of epoch in seconds (TAI58)
    ('subsecs',   '>u2'),   # fractional seconds
])                          # total length = 6 bytes

def parse_CCSDS_timestamp(timestr):
    return decode_timestamp(timestr['seconds'], timestr['subsecs'])

def readTimestamp(data):
    tmp = np.frombuffer(data, dtype=CCSDS_timestamp, count=1)
    ts = parse_CCSDS_timestamp(tmp[0])
    return ts

def tai58_as_datetime(tai58):
    # convert tai58 (seconds since 1958, including leapsecs) to datetime object
    dt = datetime.datetime.utcfromtimestamp(tai58 - LEAPSEC) - TAI58_OFFSET
    return dt  # seconds since Jan 1, 1970

def seconds_since(tai58, basetime=None):
    # translate TAI58 to seconds since a baseline time (defaults to start of day)
    dt = tai58_as_datetime(tai58)
    if not basetime:
        basetime = dt.replace(hour=0, minute=0, second=0, microsecond=0)
    return (dt - basetime).total_seconds()

def datetime_repr(dt):
    return dt.strftime('%Y-%m-%dT%H:%M:%S.%fZ')

#-------------------------------------------------------------------------------

def getDict(structured_array):
    '''
    convert numpy structured array to dict of scalars and numpy arrays
    '''
    myDict = {}

    for key in structured_array.dtype.names:
        if ( key.upper().startswith('PADDING') or
             key.upper().startswith('SPARE') ):
            continue  # skip meaningless fields
        myDict[key] = structured_array[key][0]

    return myDict

def pad_packet(data, length):
    # pad or truncate bytearray to expected length
    packet = bytearray(length)
    if len(data) > length:
        packet[:] = data[0:length]
    else:
        packet[0:len(data)] = data
    return packet

#-------------------------------------------------------------------------------

cfeHeaderFields = np.dtype([
    ('filetype', 'S4'),              # "type of file = "cFE1"
    ('subtype', '>u4'),              # 101d = downloaded PACE DS files
    ('length', '>u4'),               # length of CFE file header = 64 bytes
    ('SCID', 'S4'),                  # Spacecraft ID
    ('APID', '>u4'),                 # Application ID
    ('processorid', '>u4'),          # 1=Spacecraft, 2=OCI  (30 in ETE files)
    ('fileopen_seconds', '>u4'),     # Time when file created
    ('fileopen_subseconds', '>u4'),  # Time when file created (fractional seconds)
    ('description', 'S32'),          # = "DS data storage file"
])                                   # total length = 64 bytes

dsfHeaderFields = np.dtype([         # CFE/S Data Storage header
    ('fileclose_seconds', '>u4'),    # Time when file closed
    ('fileclose_subseconds', '>u4'), # Time when file closed (fractional seconds)
    ('filetable_index', '>u2'),      # 0=events, 1=housekeeping
    ('filename_type', '>u2'),        # 1=count, 2=time (PACE uses count)
    ('filename', 'S64'),             # fully qualified file path
])                                   # total length = 76 bytes

def readFileHeader(filehandle):

    # is a CFE header present?
    pos = filehandle.tell()
    fourchars = filehandle.read(4)
    filehandle.seek(pos) # rewind to original position
    if fourchars != b'cFE1':
        return None

    # read CFE header
    data = filehandle.read(cfeHeaderFields.itemsize)
    tmp = np.frombuffer(data, dtype=cfeHeaderFields, count=1)
    myDict = getDict(tmp)
    myDict['fileopen_ts'] = decode_timestamp(myDict['fileopen_seconds'],
                                             myDict['fileopen_subseconds'])
    myDict['fileopen'] = datetime_repr(tai58_as_datetime(myDict['fileopen_ts']))

    # append DS header, if it's there.
    if myDict['description'].startswith(b'DS'):
        data = filehandle.read(dsfHeaderFields.itemsize)
        tmp = np.frombuffer(data, dtype=dsfHeaderFields, count=1)
        myDict.update(getDict(tmp))
        myDict['fileclose_ts'] = decode_timestamp(myDict['fileclose_seconds'],
                                                  myDict['fileclose_subseconds'])
        myDict['fileclose'] = datetime_repr(tai58_as_datetime(myDict['fileclose_ts']))

    return myDict

#-------------------------------------------------------------------------------
