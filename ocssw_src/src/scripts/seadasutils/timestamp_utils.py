from seadasutils.MetaUtils import readMetadata
import seadasutils.ProcUtils as ProcUtils
import datetime
import re

def goci_timestamp(arg):
    """
        Determine the start time, stop time, and platform of a GOCI L1B file.
    """

    meta = readMetadata(arg)
    if 'Sensor name' in meta:
        sat_name = meta['Sensor name'].lower()
        stime = meta['Scene Start time']
        etime = meta['Scene end time']
        sdt_obj = datetime.datetime.strptime(stime[0:20], "%d-%b-%Y %H:%M:%S")
        edt_obj = datetime.datetime.strptime(etime[0:20], "%d-%b-%Y %H:%M:%S")
        start_time = "{}-{}-{}T{}:{}:{}".format(stime[7:11], sdt_obj.strftime("%m"), stime[0:2], stime[12:14], stime[15:17], stime[18:20])
        end_time = "{}-{}-{}T{}:{}:{}".format(etime[7:11], edt_obj.strftime("%m"), etime[0:2], etime[12:14], etime[15:17], etime[18:20])
    return start_time, end_time, sat_name

def hawkeye_timestamp(arg):
    """
        Determine the start time, stop time, and platform of a HAWKEYE L1A file.
    """

    meta = readMetadata(arg)
    if 'instrument' in meta:
        sat_name = meta['instrument'].lower()
        stime = meta['time_coverage_start']
        etime = meta['time_coverage_end']
        start_time = stime[0:19]
        end_time = etime[0:19]
    return start_time, end_time, sat_name

def goci_timestamp(arg):
    """
        Determine the start time, stop time, and platform of a GOCI L1B file.
    """

    meta = readMetadata(arg)
    if 'Sensor name' in meta:
        sat_name = meta['Sensor name'].lower()
        stime = meta['Scene Start time']
        etime = meta['Scene end time']
        sdt_obj = datetime.datetime.strptime(stime[0:20], "%d-%b-%Y %H:%M:%S")
        edt_obj = datetime.datetime.strptime(etime[0:20], "%d-%b-%Y %H:%M:%S")
        start_time = "{}-{}-{}T{}:{}:{}".format(stime[7:11], sdt_obj.strftime("%m"), stime[0:2], stime[12:14], stime[15:17], stime[18:20])
        end_time = "{}-{}-{}T{}:{}:{}".format(etime[7:11], edt_obj.strftime("%m"), etime[0:2], etime[12:14], etime[15:17], etime[18:20])
    return start_time, end_time, sat_name

def hawkeye_timestamp(arg):
    """
        Determine the start time, stop time, and platform of a HAWKEYE L1A file.
    """

    meta = readMetadata(arg)
    if 'instrument' in meta:
        sat_name = meta['instrument'].lower()
        stime = meta['time_coverage_start']
        etime = meta['time_coverage_end']
        start_time = stime[0:19]
        end_time = etime[0:19]
    return start_time, end_time, sat_name

def hico_timestamp(arg):
    """
        Determine the start time, stop time, and platform of a HICO L1B file.
    """

    meta = readMetadata(arg)
    if 'instrument' in meta:
        sat_name = meta['instrument'].lower()
        sdate = meta['Beginning_Date']
        edate = meta['Ending_Date']
        stime = meta['Beginning_Time']
        etime = meta['Ending_Time']
        start_time = '-'.join([sdate[0:4],sdate[4:6],sdate[6:8]]) + 'T' + ':'.join([stime[0:2],stime[2:4],stime[4:len(stime)]])
        end_time = '-'.join([edate[0:4],edate[4:6],edate[6:8]]) + 'T' + ':'.join([etime[0:2],etime[2:4],etime[4:len(etime)]])
    return start_time, end_time, sat_name

def meris_timestamp(arg):
    """
        Determine the start time, stop time, and platform of a MERIS L1B file.
    """

    meta = readMetadata(arg)
    if 'instrument' in meta:
        sat_name = meta['instrument'].lower()
        stime = meta['startTime']
        etime = meta['stopTime']
        start_time = stime[0:19]
        end_time = etime[0:19]
    return start_time, end_time, sat_name   

def msi_timestamp(arg):
    """
        Determine the start time, stop time, and platform of a OLCI L1B file.
    """

    meta = readMetadata(arg)
    if 'platform' in meta:
        sat_name = meta['platform'].lower()
        stime = meta['startTime']
        start_time = stime[0:19]
        end_time = None
    return start_time, end_time, sat_name 

def ocm2_timestamp(arg):
    """
        Determine the start time, stop time, and platform of a OCM2 L1B file.
    """

    meta = readMetadata(arg)
    sat_name = 'ocm2'
    stime = meta['Start Time']
    etime = meta['End Time']
    sdt_obj = datetime.datetime.strptime(stime[0:13], "%Y%j%H%M%S")
    edt_obj = datetime.datetime.strptime(etime[0:13], "%Y%j%H%M%S")
    start_time = "{}-{}-{}T{}:{}:{}".format(stime[0:4], sdt_obj.strftime("%m"), sdt_obj.strftime("%d"), stime[7:9], stime[9:11], stime[11:13])
    end_time = "{}-{}-{}T{}:{}:{}".format(etime[0:4], edt_obj.strftime("%m"), edt_obj.strftime("%d"), etime[7:9], etime[9:11], etime[11:13])
    return start_time, end_time, sat_name 

def olci_timestamp(arg):
    """
        Determine the start time, stop time, and platform of a OLCI L1B file.
    """

    meta = readMetadata(arg)
    if 'platform' in meta:
        sat_name = meta['platform'].lower()
        stime = meta['startTime']
        etime = meta['stopTime']
        start_time = stime[0:19]
        end_time = etime[0:19]
    return start_time, end_time, sat_name

def etm_timestamp(arg):
    """
        Determine the start time, stop time, and platform of a L7ETM L1B file.
    """

    meta = readMetadata(arg)
    if re.search('ETM', meta['SENSOR_ID']):
        sat_name = 'etm'
    stime = meta['SCENE_CENTER_TIME']
    sdate = meta['DATE_ACQUIRED']
    start_time = sdate[0:10] + 'T' + stime[1:9]
    end_time = None
    return start_time, end_time, sat_name

def sgli_timestamp(arg):
    """
        Determine the start time, stop time, and platform of a SGLI L1B file.
    """

    meta = readMetadata(arg)
    sat_name = 'sgli'
    stime = meta['Scene_start_time']
    etime = meta['Scene_end_time']
    start_time = '-'.join([stime[0:4], stime[4:6], stime[6:8]]) + 'T' + stime[9:17]
    end_time = '-'.join([etime[0:4], etime[4:6], etime[6:8]]) + 'T' + etime[9:17]
    return start_time, end_time, sat_name

def tm_timestamp(arg):
    """
        Determine the start time, stop time, and platform of a L5TM L1B file.
    """

    meta = readMetadata(arg)
    sat_name = 'tm'
    stime = meta['startTime']
    sdate = meta['startDate']
    start_time= sdate + "T" + stime[0:8]
    end_time = None
    return start_time, end_time, sat_name

def l9_timestamp(arg):
    """
        Determine the start time, stop time, and platform of a L5TM L1B file.
    """

    meta = readMetadata(arg)
    if re.search('OLI', meta['SENSOR_ID']):
        sat_name = meta['SENSOR_ID'].strip()
    stime = meta['SCENE_CENTER_TIME']
    sdate = meta['DATE_ACQUIRED']
    start_time = sdate[0:10] + 'T' + stime[1:9]
    end_time = None
    return start_time, end_time, sat_name