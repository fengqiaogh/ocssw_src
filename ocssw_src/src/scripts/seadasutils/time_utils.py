"""
Module to hold time utility classes and functions.
"""

from datetime import datetime, timezone, timedelta
import os

def convert_month_day_to_doy(mon, dom, yr):
    """
    Returns a day of year computed from the provided month (mon parameter),
    day of month(dom parameter), and year (yr parameter).
    """
    date_obj = datetime(int(yr), int(mon), int(dom))
    doy = date_obj.timetuple().tm_yday
    return doy

def get_leap_seconds(taitime,epochyear=1958):
    '''
    Return the number of elapsed leap seconds given a TAI time in seconds
    Requires tai-utc.dat
    '''
    
    import julian

    epochsecs = (datetime(epochyear,1,1,0,0,0,tzinfo=timezone.utc) - datetime(1970,1,1,0,0,0,tzinfo=timezone.utc)).total_seconds()
    taidt = datetime.utcfromtimestamp(taitime + epochsecs)
    taiutc = os.path.join(os.environ['OCVARROOT'], 'common', 'tai-utc.dat')
    leapsec = 0
    with open(taiutc, "r") as tdat:
        for line in tdat:
            rec = line.rstrip().split(None, 7)
            dt = julian.from_jd(float(rec[4]))
            if dt < taidt:
                leapsec = float(rec[6])

    return leapsec
