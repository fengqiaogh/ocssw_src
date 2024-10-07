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
