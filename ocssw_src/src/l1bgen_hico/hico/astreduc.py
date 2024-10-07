'''
Created on Nov 24, 2015

@author: rhealy (richard.healy@nasa.gov)

'''


import sys
from datetime import datetime as DT
from datetime import timedelta as TD
import numpy as np

MMM={'JAN':1,'FEB':2,'MAR':3,'APR':4,'MAY':5,'JUN':6,'JUL':7,'AUG':8,'SEP':9,'OCT':10,'NOV':11,'DEC':12}

edf_delimiter = [2,2,2,9,3,9,9,1,9,9,3,10,10,1,7,7]
edf_usecols = [0,1,2,5,8,11,14]
edf_dtype = [np.int32,np.int32,np.int32,np.float64,np.float64,np.float64,np.float64]
edf_names = ['year','month','day','pmx','pmy','ut1mutc','loda']

ls_delimiter = [6,3,3,5,10,10,12]
ls_usecols = [0,1,2,6]
ls_dtype = [np.int32,(np.str_,3),np.int32,np.float64]
ls_names = ['yfo','mco','dfo','taimutco']

rad2Deg = 180.0/np.pi
deg2Rad = np.pi/180.0
as2Deg = 1.0/3600.0
as2Rad = as2Deg*deg2Rad
"""
if __name__ == '__main__':
    import os
    try:
        ocvarroot = os.environ['OCVARROOT']
    except KeyError:
        print("OCSSW environement variable OCVARROOT not set.")
        print("Typically it is set to $OCSSW/run/var.")
        sys.exit()
    filein = "{}/hico/{}".format(ocvarroot,'nutation.dat')

    nut = initReduc(filein)

    print(nut.head(10))

class HicoSec2Date():
    def __init__(self, delta):
        # From MJM 2012/10/24 Adjusting the time from the epoch (2009-01-01-00:00:00 UTC).
        # Input is GPS seconds since the epoch. Output is UTC date/time.
        # adapted to python by RJH @ NASA 11/24/2015
        mlen = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
        dinc=delta
        idinc=np.floor(delta)
        #Seconds since epoch of the (only) leap second so far
        ls_2012_06_30 = 86400*((365*3)+(31+29+31+30+31+30))+1
        # use epoch definition as start
        self.yyyy=2009
        self.mm=1
        self.dd=1
        self.hour= np.dtype('int32')
        self.min = np.dtype('int32')
        self.hour=0
        self.min=0
        self.sec=0.e0
        if idinc > ls_2012_06_30:
            #date after leap second, so adjust this to new epoch
            self.yyyy=2012
            self.mm=7
            self.dd=1
            idinc=idinc-ls_2012_06_30
            dinc=dinc-ls_2012_06_30
            self.hour=np.int(0)
            self.min=np.int(0)
            self.sec=0.e0

        mm = self.mm
        mmlen=mlen[mm - 1]
        jdinc=np.dtype('int32')
        while idinc >= 86400*mmlen:
            if mm == 2:
                if isLeapYear(self.yyyy) == 1:
                    mmlen = mlen[mm-1]+1
                else:
                    mmlen = mlen[mm-1]
            else:
                mmlen = mlen[mm-1]
            jdinc = 86400*mmlen
            idinc=idinc-jdinc
            dinc=dinc - 86400*mmlen
            if mm != 12:
                mm=mm+1
            else:
                mm=1
                self.yyyy=self.yyyy+1

        self.mm = mm
# at this point, we have less a full month. Let's figure out how many whole days
        nd=np.int(idinc/86400)
        self.dd=self.dd+nd
        idinc=idinc - 86400*nd
        dinc=dinc - 86400*nd
# Less than 24 hours left
        self.hour=np.int(idinc/3600)
        idinc=idinc - self.hour*3600
        dinc=dinc - self.hour*3600
# Less than 1 hour left
        self.min=np.int(idinc/60)
        idinc=idinc - self.min*60
        dinc=dinc - self.min*60
# Less than 1 minute left
        self.sec=dinc"""

def get_astrodate(delta_seconds):
    d_epoch = DT(year=2009, month=1, day=1)
    delta_t = TD(seconds=delta_seconds)
    return d_epoch + delta_t

def isLeapYear(yyyy):
    if (yyyy % 4) == 0:
        if (yyyy % 100) == 0:
            leap=1
        else:
            if (yyyy % 400) == 0:
                leap=1
            else:
                leap=0
    else:
        leap=0

    return leap

def increment_seconds_to_date(delta):
    # From MJM 2012/10/24 Adjusting the time from the epoch (2009-01-01-00:00:00 UTC).
    # Input is GPS seconds since the epoch. Output is UTC date/time.
    # adapted to python by RJH @ NASA 11/24/2015
    mlen = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    dinc=delta
    idinc=int(np.floor(delta))
    #Seconds since epoch of the (only) leap second so far
    ls_2012_06_30 = 86400*((365*3)+(31+29+31+30+31+30))+1
    #else use epoch definition as start
    yyyy=2009
    mm=1
    dd=1
    if idinc > ls_2012_06_30:
        #date after leap second, so adjust this to new epoch
        yyyy=2012
        mm=7
        dd=1
        idinc=idinc-ls_2012_06_30
        dinc=dinc-ls_2012_06_30

    mmlen=mlen[mm - 1]
    jdinc=np.dtype('int32')
    while idinc >= 86400*mmlen:
        if mm == 2:
            if isLeapYear(yyyy) == 1:
                mmlen = mlen[mm-1]+1
            else:
                mmlen = mlen[mm-1]
        else:
            mmlen = mlen[mm-1]
        jdinc = 86400*mmlen
        dinc=dinc - 86400*mmlen
        idinc=idinc-jdinc
        if mm != 12:
            mm=mm+1
        else:
            mm=1
            yyyy=yyyy+1
# at this point, we have less a full month. Let's figure out how many whole days
    nd=int(idinc/86400)
    dd=dd+nd
    idinc=idinc - 86400*nd
    dinc=dinc - 86400*nd
# Less than 24 hours left
    hour=int(idinc/3600)
    idinc=idinc - hour*3600
    dinc=dinc - hour*3600
# Less than 1 hour left
    minute=int(idinc/60)
    idinc=idinc - min*60
    dinc=dinc - min*60
# Less than 1 minute left
    sec=dinc

    return yyyy,mm,dd,hour,minute,sec

def initReduc(filename):
    import pandas as pd
    Convrt= 0.0001e0/3600.0e0 # 0.0001 # to deg
    nutation_data = pd.read_csv(filename,skiprows=110,skipfooter=552-216,skipinitialspace=True,sep=' ',
                               names=['a1','a2','a3','a4','a5','A','B','C','D','NUM'], engine='python')
    nutation_data['A'] = nutation_data['A']*Convrt
    nutation_data['B'] = nutation_data['B']*Convrt
    nutation_data['C'] = nutation_data['C']*Convrt
    nutation_data['D'] = nutation_data['D']*Convrt

    return nutation_data

def getEDFData(filename,year,mon,day):
    return np.concatenate([ [row[3]*as2Rad,row[4]*as2Rad,row[5],row[6]/1000.] for row in parseEDFfile(filename) if (row[0] == year and row[1] == mon and row[2] == day )])

def parseEDFfile(filename,skiphdr=0):
    return np.genfromtxt(filename,skip_header=skiphdr,
                         delimiter = edf_delimiter,
                         usecols = edf_usecols,
                         dtype = edf_dtype,
                         names = edf_names)

def getLSData(filename,year,mon,day):
    mdate=year*10000 + mon*100 + day
    lsdat = np.concatenate([ [row[3]] for row in parseLSfile(filename,1) if (mdate > (row[0]*10000 + MMM[row[1]]*100 + row[2] ) )])

    return lsdat[-1]

def parseLSfile(filename,skiphdr=0):
    return np.genfromtxt(filename,skip_header=skiphdr,
                         delimiter = ls_delimiter,
                         usecols = ls_usecols,
                         dtype = ls_dtype,
                         names = ls_names)

''''
converted to python by R. Healy 11/27/2015
* ----------------------------------------------------------------------------
*
*                           SUBROUTINE POLARM
*
*  this function calulates the transformation matrix for polar motion.
*    the units for polar motion are input in rad because it's more
*    efficient to do this in the main routine.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    xp          - Polar motion coefficient       rad
*    yp          - Polar motion coefficient       rad
*
*  Outputs       :
*    PM          - Polar motion transformation (pef-ecef)
*
*  Locals        :
*    None.
*
*  Coupling      :
*    ROT2        - Rotation about the second axis
*    ROT3        - Rotation about the third axis
*
*  References    :
*    Vallado       2007, 228
*
* ----------------------------------------------------------------------------
'''
def polarm(xp,yp):


    cosxp = np.cos(xp)
    sinxp = np.sin(xp)
    cosyp = np.cos(yp)
    sinyp = np.sin(yp)

    pm = np.zeros((3,3),dtype=np.float64)

    pm[0][0] =  cosxp
    pm[0][1] =  sinxp * sinyp
    pm[0][2] =  sinxp * cosyp

    pm[1][0] =  0.0
    pm[1][1] =  cosyp
    pm[1][2] = -sinyp

    pm[2][0] = -sinxp
    pm[2][1] =  cosxp * sinyp
    pm[2][2] =  cosxp * cosyp

    return pm

def quat_to_rot(q):
    '''
    #*****************************************************************************80
    # Converted to Python by R. Healy 11/30/2015
    #  from q_to_r.f90
    #
    # Very slightly modified my Marcos Montes.
    # MJM Note: scalar is q(1).
    ## ROTATION_QUAT2MAT_3D converts rotation from quaternion to matrix form in 3D.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    27 July 1999
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Reference:
    #
    #    James Foley, Andries van Dam, Steven Feiner, John Hughes,
    #    Computer Graphics, Principles and Practice,
    #    Second Edition,
    #    Addison Wesley, 1990.
    #
    #  Parameters:
    #
    #    Input, real ( kind = 8 ) Q(4), the quaternion representing the rotation.
    #
    #    Output, real ( kind = 8 ) a[3,3), the rotation matrix.
    #
    '''
    sin_phi = np.sqrt ( np.sum ( np.square(q[1:4]) ) )

    cos_phi = q[0]

    angle = 2.0 * np.arctan2( sin_phi, cos_phi )

    if sin_phi == 0.0:
        v1 = 1.0
        v2 = 0.0
        v3 = 0.0
    else:
        v1 = q[1] / sin_phi
        v2 = q[2] / sin_phi
        v3 = q[3] / sin_phi

    ca = np.cos ( angle )
    sa = np.sin ( angle )

    a = np.zeros((3,3))
    a[0,0] = v1 * v1 + ca * ( 1.0 - v1 * v1 )
    a[0,1] = ( 1.0 - ca ) * v1 * v2 - sa * v3
    a[0,2] = ( 1.0 - ca ) * v1 * v3 + sa * v2

    a[1,0] = ( 1.0 - ca ) * v2 * v1 + sa * v3
    a[1,1] =  v2 * v2 + ca * ( 1.0 - v2 * v2 )
    a[1,2] = ( 1.0 - ca ) * v2 * v3 - sa * v1

    a[2,0] = ( 1.0 - ca ) * v3 * v1 - sa * v2
    a[2,1] = ( 1.0 - ca ) * v3 * v2 + sa * v1
    a[2,2] = v3 * v3 + ca * ( 1.0 - v3 * v3 )

    return a


def hms2ut(Hr, Minute, Sec):
    '''
    * Converted to python by R. Healy 11/30/2015
    * broken into two functions from fortran sub HMS_SEC:
            hms2ut
            ut2hms
    * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    *
    *                              SUBROUTINE HMS_SEC
    *
    *  this subroutine converts Hours, Minutes and Seconds into seconds from the
    * beginning of the day.
    *
    *  Author        : David Vallado                  719 - 573 - 2600    1 Mar 2001
    *
    *  Inputs          Description                    Range / Units
    * Hr - Hours                          0 .. 24
    * minute - Minutes                        0 .. 59
    * Sec - Seconds                        0.0D0 .. 59.99D0
    * Direction - Which set of vars to output    FROM  TOO
    *
    *  OutPuts       :
    *    Sec - Seconds                        0.0D0 .. 86400.0D0
    *
    *  Locals        :
    *    Temp - Temporary variable
    *
    *  Coupling      :
    *    None.
    *
    * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
'''
    UTSec= Hr*3600.0 + Minute*60.0 + Sec

    return UTSec

def ut2hms(UTSec):
    '''
    * Converted to python by R. Healy 11/30/2015
    * broken into two functions from fortran sub HMS_SEC:
            hms2ut
            ut2hms
    * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    *
    *                              SUBROUTINE HMS_SEC
    *
    *  this subroutine converts Hours, Minutes and Seconds into seconds from the
    * beginning of the day.
    *
    *  Author        : David Vallado                  719 - 573 - 2600    1 Mar 2001
    *
    *  Inputs          Description                    Range / Units
    * Hr - Hours                          0 .. 24
    * minute - Minutes                        0 .. 59
    * Sec - Seconds                        0.0D0 .. 59.99D0
    * Direction - Which set of vars to output    FROM  TOO
    *
    *  OutPuts       :
    *    Sec - Seconds                        0.0D0 .. 86400.0D0
    *
    *  Locals        :
    *    Temp - Temporary variable
    *
    *  Coupling      :
    *    None.
    *
    * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
'''

    Temp= UTSec / 3600.0
    Hr  = np.floor( Temp )
    Minute = np.floor( (Temp - Hr)*60.0 )
    Sec = (Temp - Hr - Minute/60.0 ) * 3600.0

    return Hr, Minute, Sec


def JDay        ( Year,Mon,Day,Hr,minute, Sec ):
    '''
    * Converted to python by R. Healy
    * -----------------------------------------------------------------------------
    *
    *                           SUBROUTINE JDay
    *
    *  this subroutine finds the Julian date given the Year, Month, Day, and Time.
    *
    *  Author        : David Vallado                  719-573-2600    1 Mar 2001
    *
    *  Inputs          Description                    Range / Units
    *    Year        - Year                           1900 .. 2100
    *    Mon         - Month                          1 .. 12
    *    Day         - Day                            1 .. 28,29,30,31
    *    Hr          - Universal Time Hour            0 .. 23
    *    minute         - Universal Time minute             0 .. 59
    *    Sec         - Universal Time Sec             0.0D0 .. 59.999D0
    *    WhichType   - Julian .or. Gregorian calender   'J' .or. 'G'
    *
    *  Outputs       :
    *    JD          - Julian Date                    days from 4713 BC
    *
    *  Locals        :
    *    B           - Var to aid Gregorian dates
    *
    *  Coupling      :
    *    None.
    *
    *  References    :
    *    Vallado       2007, 189, Alg 14, Ex 3-14
    * -----------------------------------------------------------------------------
    '''

    jday = 367.0 * Year \
             - np.int( (7* (Year+np.int ( (Mon+9)/12) ) ) * 0.25 ) \
             + np.int( 275*Mon / 9 ) \
             + Day + 1721013.5 \
             + ( (Sec/60.0 + minute ) / 60.0 + Hr ) / 24.0
    return jday

class UT2time():
    '''
    * Converted to python by R. Healy 11/30/2015
    * ------------------------------------------------------------------------------
    *
    *                           SUBROUTINE CONVTIME
    *
    *  this subroutine finds the time parameters and Julian century values for inputs
    *    of UTC or UT1. Numerous outputs are found as shown in the local variables.
    *    Because calucations are in UTC, you must include TimeZone IF ( you enter a
    *    local time, otherwise it should be zero.
    *
    *  Algorithm     : A file of record contains the timing data
    *                  Seeks are performed to obtain the data
    *                    Data starts Jan 1, 1980, thus JD = 2444238.5D0 in the code
    *                  Calculate the answer depending on initial time type
    *
    *  Author        : David Vallado                  719-573-2600    1 Mar 2001
    *
    *  Inputs          Description                    Range / Units
    *    Year        - Year                           1900 .. 2100
    *    Mon         - Month                          1 .. 12
    *    Day         - Day                            1 .. 28,29,30,31
    *    Hr          - Universal Time Hour            0 .. 23
    *    minute      - Universal Time minute          0 .. 59
    *    SEC         - Universal Time SEC             0.0D0 .. 59.999D0
    *    TimeZone    - Offset to UTC from local SITE  0 .. 23 hr
    *    TypeUTIn    - Type of input UT               1 (UT1), else UTC
    *    DUT1        - Delta of UTC - UT1             SEC
    *
    *  Outputs       :
    *    DAT         - Delta of TAI-UTC              SEC [MJM: THIS IS AN INPUT]
    *    xp          - Polar motion coefficient       arcsec [MJM: INPUT]
    *    yp          - Polar motion coefficient       arcsec [MJM: INPUT]
    *    UT1         - Universal time                 SEC
    *    TUT1        - Julian centuries of UT1
    *    JDUT1       - Julian Date of UT1             days from 4713 BC
    *    UTC         - Coordinated Universal Time     SEC
    *    TAI         - Atomic time                    SEC
    *    TDT         - Terrestrial Dynamical time     SEC
    *    TTDT        - Julian centuries of TDT
    *    JDTDT       - Julian Date of TDT             days from 4713 BC
    *    TDB         - Terrestrial Barycentric time   SEC
    *    TTDB        - Julian centuries of TDB
    *    JDTDB       - Julian Date of TDB             days from 4713 BC
    *
    *  Locals        :
    *    HrTemp      - Temporary hours                hr
    *    MinTemp     - Temporary miNutes              minute
    *    SecTemp     - Temporary seconds              SEC
    *    LocalHr     - Difference to local time       hr
    *    JD          - Julian Date of request         days from 4713 BC
    *    ME          - Mean Anomaly of the Earth      rad
    *    TimeFile    - File of record with time data
    *    CurrTimeRec - Current Time record
    *
    *  Coupling      :
    *    HMS_SEC     - Conversion between hr-minute-SEC .and. seconds
    *    jday   - Find the Julian date
    *
    *  References    :
    *    vallado       2007, 201, alg 16, ex 3-7
    *
    * ------------------------------------------------------------------------------
    '''
    def __init__(self, Year, Mon,Day, Hr, Minute, Sec, dUT1,lsdat,TimeZone=0,TypeUTIn='UTC'):

        self.year = Year
        self.mon  = Mon
        self.day  = Day
        self.hr   = Hr
        self.minute = Minute
        self.sec    = Sec
        self.typeUTIn = TypeUTIn
        self.timezone = TimeZone
        self.localhr  = TimeZone + Hr
        if TypeUTIn == 'UT1' :
            self.ut1 = hms2ut(self.localhr,Minute,Sec)
            self.jdut1 = JDay(Year,Mon,Day, self.localhr, Minute, Sec)
            self.tut1  = (self.jdut1 - 2451545.0 )/ 36525.0
            self.utc   = self.ut1 - dUT1
        else:
            self.utc = hms2ut(self.localhr,Minute,Sec)
            self.ut1 = self.utc + dUT1
            hrTemp,minTemp,secTemp = ut2hms(self.ut1)
            self.jdut1 = JDay(Year,Mon,Day, np.int(hrTemp), np.int(minTemp),secTemp)
            self.tut1  = (self.jdut1 - 2451545.0 )/ 36525.0

        self.tai = self.utc + lsdat
        self.tt = self.tai + 32.184
        hrTemp,minTemp,secTemp = ut2hms(self.tt)
        self.jdtt = JDay(Year,Mon,Day, hrTemp, minTemp, secTemp)
        self.ttt = (self.jdtt - 2451545.0 )/ 36525.0

        #self.me = 357.5277233 + 35999.05034*self.ttt   # approx - should do with TTDB
        self.me = np.fmod((357.5277233 + 35999.05034*self.ttt),360.0)*deg2Rad
        self.tdb = self.tt + 0.001658 * np.sin(self.me) + 0.00001385*np.sin(2.0*self.me)
        hrTemp,minTemp,secTemp = ut2hms(self.tdb)
        self.jdtdb = JDay(Year,Mon,Day, hrTemp, minTemp, secTemp)
        self.ttdb  =  (self.jdtdb - 2451545.0)/ 36525.0

def precession(TTT):
    '''
#     * ----------------------------------------------------------------------------
#     *
#     *                           SUBROUTINE PRECESSION
#     *
#     *  this function calulates the transformation matrix for precession.
#     *
#     *  Author        : David Vallado                  719-573-2600    1 Mar 2001
#     *
#     *  Inputs          Description                    Range / Units
#     *    TTT         - Julian Centuries of TT         centuries
#     *
#     *  Outputs       :
#     *    Prec        - Precession transformation (eci-mod)
#     *
#     *  Locals        :
#     *    TTT2        - TTT squared
#     *    TTT3        - TTT cubed
#     *    Zeta        - PRECESSION ANGLE               rad
#     *    z           - PRECESSION ANGLE               rad
#     *    Theta       - PRECESSION ANGLE               rad
#     *
#     *  Coupling      :
#     *    none
#     *
#     *  References    :
#     *    Vallado       2007, 228
#     *
#     * ----------------------------------------------------------------------------
    '''
# MJM note: Reid Reynolds uses TTDB not TTT.

    # --------------------- PRECESSION angles ---------------------
    # MJM These are in arcsec; look like values in Reid Reynolds memo
    zeta= (( 0.017998*TTT + 0.30188)*TTT + 2306.2181)*TTT
    theta = (( - 0.041833*TTT - 0.42665)*TTT + 2004.3109)*TTT
    z     = (( 0.018203*TTT + 1.09468)*TTT + 2306.2181)*TTT

    zeta = zeta  * deg2Rad / 3600.0 # convert to radians
    theta= theta * deg2Rad / 3600.0
    z    = z     * deg2Rad / 3600.0

    coszeta   = np.cos(zeta)
    sinzeta   = np.sin(zeta)
    costheta = np.cos(theta)
    sintheta = np.sin(theta)
    cosz    = np.cos(z)
    sinz    = np.sin(z)

    prec=np.zeros((3,3))
    # ----------------- form matrix  J2000 to MOD -----------------
    prec[0,0] =  coszeta * costheta * cosz - sinzeta * sinz
    prec[0,1] = -sinzeta * costheta * cosz - coszeta * sinz
    prec[0,2] = -sintheta * cosz

    prec[1,0] =  coszeta * costheta * sinz + sinzeta * cosz
    prec[1,1] = -sinzeta * costheta * sinz + coszeta * cosz
    prec[1,2] = -sintheta * sinz

    prec[2,0] =  coszeta * sintheta
    prec[2,1] = -sinzeta * sintheta
    prec[2,2] =  costheta

    return prec

def nutation(TTT,nut_data):
    '''
    # * ----------------------------------------------------------------------------
    # *
    # *                           SUBROUTINE NUTATION
    # *
    # *  this function calulates the transformation matrix for nutation.
    # *
    # *  Author        : David Vallado                  719-573-2600    1 Mar 2001
    # *
    # *  Inputs          Description                    Range / Units
    # *    TTT         - Julian Centuries of TT
    # *
    # *  Outputs       :
    # *    Nut         - Nutation transformation (mod-tod)
    # *    DeltaPsi    - NUTATION ANGLE                 rad
    # *    TrueEps     - True obliquity of the ecliptic rad
    # *    Omega       -                                rad
    # *    MeanEps     - Mean obliquity of the ecliptic rad
    # *
    # *  Locals        :
    # *    TTT2        - TTT squared
    # *    TTT3        - TTT cubed
    # *    MeanEps     - Mean obliquity of the ecliptic rad
    # *    l           -                                rad
    # *    ll          -                                rad
    # *    F           -                                rad
    # *    D           -                                rad
    # *    DeltaEps    - Change in obliquity            rad
    # *
    # *  Coupling      :
    # *    none
    # *
    # *  References    :
    # *    Vallado       2007, 228
    # *
    # * ----------------------------------------------------------------------------
    '''
    # ---- Determine coefficients for IAU 1980 NUTATION Theory ----
    TTT2= TTT*TTT
    TTT3= TTT2*TTT
    TTT4= TTT2*TTT2

# Meaneps is in arcseconds first

    MeanEps = ((0.001813*TTT - 0.00059)*TTT - 46.8150)*TTT + \
                84381.448

    MeanEps = np.fmod(MeanEps/3600.0 , 360.0)  # degrees, [0-360)
    MeanEps = MeanEps * deg2Rad                # radians

    l    =  134.96340251 + ( 1717915923.2178*TTT +          \
           31.8792*TTT2 + 0.051635*TTT3 - 0.00024470*TTT4 ) \
           / 3600.0
    l1   =  357.52910918 + (  129596581.0481*TTT -           \
            0.5532*TTT2 - 0.000136*TTT3 - 0.00001149*TTT4 ) \
           / 3600.0
    F    =   93.27209062 + ( 1739527262.8478*TTT -           \
           12.7512*TTT2 + 0.001037*TTT3 + 0.00000417*TTT4 ) \
           / 3600.0
    D    =  297.85019547 + ( 1602961601.2090*TTT -           \
            6.3706*TTT2 + 0.006593*TTT3 - 0.00003169*TTT4 ) \
           / 3600.0
    Omega=  125.04455501 + (   -6962890.2665*TTT +           \
            7.4722*TTT2 + 0.007702*TTT3 - 0.00005939*TTT4 ) \
           / 3600.0
# ! Above are in degrees and they correspond to:
# ! l=mean anomaly of moon
# ! l1 = mean anomaly of sun
# ! F = mean latitude of moon
# ! D = mean elongation of sun
# ! omega = ascending node of the moon
# !     Because the coefficients (of above in degrees) are all integers,
# !     we can safely MOD to get into range of 0-360; but not necessary.
    l    = np.fmod( l , 360.0 )     * deg2Rad
    l1   = np.fmod( l1 , 360.0 )    * deg2Rad
    F    = np.fmod( F , 360.0 )     * deg2Rad
    D    = np.fmod( D , 360.0 )     * deg2Rad
    Omega= np.fmod( Omega , 360.0 ) * deg2Rad

    DeltaPsi= 0.0
    DeltaEps= 0.0

    for i in range(len(nut_data)-1,-1,-1):
        TempVal= nut_data['a1'][i]*l + nut_data['a2'][i]*l1 + nut_data['a3'][i]*F + \
                nut_data['a4'][i]*D + nut_data['a5'][i]*Omega
        DeltaPsi= DeltaPsi + (nut_data['A'][i]+ nut_data['B'][i]*TTT) * \
                  np.sin( TempVal )
        DeltaEps= DeltaEps + (nut_data['C'][i]+ nut_data['D'][i]*TTT) * \
                  np.cos( TempVal )

# --------------- Find NUTATION Parameters --------------------
# ! MJM begin
# !     Now, Reid Reynolds document with the SAME values says they are in
# !     microarcseconds and microarcseconds per century (RAR80 terms).
# !     other references show these units (for same values!!!) as 0.1mas
# !          deltaPsi = deltaPsi* 1.d-4 / 3600. ! convert to as, then deg
# !          deltaEps = deltaEps* 1.d-4 / 3600. ! convert to as, then deg
# ! ACTUALLY THESE ARE CONVERTED WHEN THEY ARE READ!!!! See INITREDUC
# ! MJM end
    DeltaPsi = np.fmod(DeltaPsi , 360.0 ) * deg2Rad
    DeltaEps = np.fmod(DeltaEps , 360.0 ) * deg2Rad
    TrueEps  = MeanEps + DeltaEps

# Checked order against Reid Reynolds. Multiplication is correct.
    cospsi  = np.cos(DeltaPsi)
    sinpsi  = np.sin(DeltaPsi)
    coseps  = np.cos(MeanEps)
    sineps  = np.sin(MeanEps)
    costrueeps = np.cos(TrueEps)
    sintrueeps = np.sin(TrueEps)

    nut = np.zeros((3,3))
    nut[0,0] =  cospsi
    nut[0,1] = -coseps * sinpsi
    nut[0,2] = -sineps * sinpsi

    nut[1,0] =  costrueeps * sinpsi
    nut[1,1] =  costrueeps * coseps * cospsi + sintrueeps * sineps
    nut[1,2] =  costrueeps * sineps * cospsi - sintrueeps * coseps

    nut[2,0] =  sintrueeps * sinpsi
    nut[2,1] =  sintrueeps * coseps * cospsi - sineps * costrueeps
    nut[2,2] =  sintrueeps * sineps * cospsi + costrueeps * coseps

    return DeltaPsi, TrueEps, MeanEps, Omega, nut

def sidereal(jdut1,DeltaPsi,MeanEps,Omega,LOD,terms):
    # * ----------------------------------------------------------------------------
    # *
    # *                           SUBROUTINE SIDEREAL
    # *
    # *  this function calulates the transformation matrix that accounts for the
    # *    effects of nutation. Notice that deltaspi should not be moded to a
    # *    positive number because it is multiplied rather than used in a
    # *    trigonometric argument.
    # *
    # *  Author        : David Vallado                  719-573-2600    1 Mar 2001
    # *
    # *  Inputs          Description                    Range / Units
    # *    JDUT1       - Julian Date of UT1             days from 4713 BC
    # *    DeltaPsi    - NUTATION ANGLE                 rad
    # *    MeanEps     - Mean obliquity of the ecliptic rad
    # *    Omega       -                                rad
    # *    LOD         - Excess length of day           sec
    # *    terms       - number of terms to include with ast 0, 2
    # *
    # *  Outputs       :
    # *    St          - Sidereal Time transformation (tod-pef)
    # *    StDot       - Sidereal Time rate transformation (tod-pef)
    # *    Omegaearth  - rotation of the earth          rad
    # *
    # *  Locals        :
    # *    GST         - Mean Greenwich SIDEREAL Time   0 to 2Pi rad
    # *    AST         - Apparent GST                   0 to 2Pi rad
    # *    Hr          - hour                           hr
    # *    minute         - minutes                        minute
    # *    SEC         - seconds                        SEC
    # *
    # *  Coupling      :
    # *    none
    # *
    # *  References    :
    # *    Vallado       2007, 228
    # *
    # * ----------------------------------------------------------------------------

    OmegaEarth = np.zeros(3)

    Conv1 = np.pi / (180.0*3600.0)

    # ------------------------ Find Mean GST ----------------------
    gst= gstime( jdut1 )

    # ------------------------ Find Mean AST ----------------------
    if terms > 0 and jdut1 > 2450449.5:
        ast = gst + DeltaPsi* np.cos(MeanEps) \
            + 0.00264*Conv1*np.sin(Omega)     \
            + 0.000063*Conv1*np.sin(2.0*Omega)
    else:
        ast = gst + DeltaPsi* np.cos(MeanEps)

    st = np.zeros((3,3))
    st[0,0] =  np.cos(ast)
    st[0,1] =  np.sin(ast)
    st[0,2] =  0.

    st[1,0] = -np.sin(ast)
    st[1,1] =  np.cos(ast)
    st[1,2] =  0.

    st[2,0] =  0.
    st[2,1] =  0.
    st[2,2] =  1.

    # ------------ compute sidereal time rate matrix --------------
    ThetaSa   =  7.29211514670698e-05 * (1.0 - LOD/86400.0)
    OmegaEarth[2] = ThetaSa

    stdot = np.zeros((3,3))
    stdot[0,0] = -OmegaEarth[2] * np.sin(ast)
    stdot[0,1] =  OmegaEarth[2] * np.cos(ast)
    stdot[0,2] =  0.0

    stdot[1,0] = -OmegaEarth[2] * np.cos(ast)
    stdot[1,1] = -OmegaEarth[2] * np.sin(ast)
    stdot[1,2] =  0.0

    stdot[2,0] =  0.0
    stdot[2,1] =  0.0
    stdot[2,2] =  0.0

    return st,stdot,OmegaEarth

def gstime(jd):
    # * -----------------------------------------------------------------------------
    # *
    # *                           FUNCTION GSTIME
    # *
    # *  this function finds the Greenwich sidereal time (iau-82).
    # *
    # *  Author        : David Vallado                  719-573-2600    1 Mar 2001
    # *
    # *  Inputs          Description                    Range / Units
    # *    JD          - Julian Date                    days from 4713 BC
    # *
    # *  OutPuts       :
    # *    GSTIME      - Greenwich SIDEREAL Time        0 to 2Pi rad
    # *
    # *  Locals        :
    # *    Temp        - Temporary variable for reals   rad
    # *    TUT1        - Julian Centuries from the
    # *                  Jan 1, 2000 12 h epoch (UT1)
    # *
    # *  Coupling      :
    # *
    # *  References    :
    # *    vallado       2007, 193, Eq 3-43
    # * -----------------------------------------------------------------------------


    TUT1= ( jd - 2451545.0 ) / 36525.0
# !!! MJM Commented out this implmenetation !!!!!
#         T2= - 6.2D-6*TUT1*TUT1*TUT1
#      &        + 0.093104D0*TUT1*TUT1
#      &        + (876600.0D0*3600.0D0 + 8640184.812866D0)*TUT1
#      &        + 67310.54841D0
#         T2= DMOD( T2*Deg2Rad/240.0D0,TwoPi ) ! 360/86400 = 1/240, to deg, to rad
# !!! END MJM comment out
    TUT1=0.13090052920192846
# Begin MJM add (based on Reid Reynolds), better numerically
    temp= 280.46061837 + 360.98564736629 * (jd - 2451545.0) + \
        ( 0.000387933 - (TUT1/38710000.0))*TUT1**2
    temp=np.fmod(temp*deg2Rad,2*np.pi )
# END MJM added

    # ------------------------ Check quadrants --------------------
    if temp < 0.0:
        temp += 2*np.pi

    return temp

def ecef2latlon(xyz):
    #     ! Need to check this set of equations; it seems pretty good,
    #     ! but now there are different equations on Wikipedia. See the papers
    #     ! I downloaded to the ReferenceFrames folder.
    #
    #     ! MJM copied algorithm from web page
    #     ! http://en.wikipedia.org/wiki/Geodetic_system
    #     ! COnverted to IDL 2009 Sept. 29.
    #     ! Converted to F90 2012 Oct 25
    #     ! 2012/10/31 Converted to Vermeille's algorithm since at least I have a referece. But performance looks
    #     ! exactly the same.
    #     ! 2012/10/31 Changed output order to LON (deg), LAT (DEG), height (m)
    #     !
    #     ! ASSUMES WGS 84 ellipsoid
    #     ! Results are GEODETIC latitude and longitude
    #     ! These will be different from Geocentric latitude
    a=6378137.0
    a2=a**2
    f=1./298.257223563
    e2= 2*f-f**2
    e=np.sqrt(e2)
    e4=e2**2
# ! Vermeille's Algorithm (Journal of Geodesy (2002) 76:451-454
# !    e=sqrt(e2)
# !    a2=a**2

    p=(xyz[0,:]**2 + xyz[1,:]**2)/a2
    q=(1-e2)/a2*xyz[2,:]**2
    r=(p+q-e4)/6
    s=e4*p*q/4/r**3
    t=(1. + s + np.sqrt(s*(2. + s)))**(1./3.)
    u=r*(1. + t + 1./t)
    v=np.sqrt(u**2+q*e4)
    w=e2*(u+v-q)/2/v
    k=np.sqrt(u+v+w**2)-w
    d=k*np.sqrt(xyz[0,:]**2 + xyz[1,:]**2)/(k+e2)

    tmp=np.sqrt(d**2 + xyz[2,:]**2)

    llh = np.zeros((3,len(xyz[0,:])))
    llh[0,:] = np.arctan2(xyz[1,:],xyz[0,:])     # longitude in rad
    llh[1,:] = 2*np.arctan2(xyz[2,:],(d + tmp))  # latitude in rad
    llh[2,:] = tmp*(k+e2-1)/k                    # height in same units as "a", i.e., meters

    return llh

def wgs84_intercept(rsc,u):
    # rsc ! spacecraft radius
    # u   ! unit vector pointing toward earth, a lot of them
    ff=1.0/298.257223563   # WGS-84
    re=6378137.            # WGS-84
    #     ! F is really the diagonal elements of a 3x3 matrix, all rest are 0
    #     ! but we only need the diagonal for the WGS-84 ellipsoid fo rvery simple
    #     ! mathematics, below.
    F=[1.,1.,1./(1.-ff)**2]

    # ! Using F as above (but see the comments - really F as the 3x3 matrix)
    # ! If Rg is a vector from the center to the surface of the earth, then
    # ! the equation of the WGS-84 ellipsoid is transpose(Rg)F*Rg=Re*Re
    #
    # ! The equation of a ray with unit direction (vector) u from the spacecraft
    # ! (at vector position Rsc) to the earth is the vector equation:
    # !  Rg = Rsc + s*u ; s is length of ray from the spacecraft to earth;
    # ! substitue this vector equation into the one above, and use the form
    # ! of the various vectors and matrices to get the equations below,
    # ! which yields a simple quadratic equation for s. For us, Rsc is one
    # ! location, and we have 512 unit direction vectors (elements of u).
    #
    # ! precalculate a few things, using the magic rules of sum () and  matmul()

    c=sum(F*rsc**2)-re**2  # scalar, POSITIVE since spacecraft is above earth surface

    b=2.*np.dot(F*rsc,u) # ns elements; negative since we're looking at earth and this is

    #     ! essentially a dot prod of spacecraft vector (away from center of earth) and view vector
    #     ! (towards earth).

    a=np.dot(F,np.square(u))  # ns elements, positive, since it is like u^2.

    det= np.square(b) - 4*np.dot(a,c)

    if (det < 0).any():
        print('ERROR IN WGS84_intercept: invalid answer. no intercept')
        print('CHECK INPUT!')
        print((det < 0))
        sys.exit


    # ! Note that -b is positive. So the closest root, the one with the smallest s,
    # ! that is the smallest distance from the spacecraft (the one on the spacecraft
    # ! side of the earth) is the one with the negative sign before the sqrt(). The OTHER one (positive)
    # ! is when the ray emerges from the earth.
    # ! Now the distance along the ray, from the spacecraft to the near intercept is:
    s= (-b -np.sqrt(det))/(2. * a ) # ns elements
    # ! Once you know s, rout is the the location in ECEF of the intercept with the Earth,
    # ! traveling a distance s from rsc in the direction(s) u.
    nr = len(u[:,0])
    ns = len(u[0,:])
    rout = np.zeros((nr,ns))
    #rout= spread(rsc,dim=2,ncopies=ns) + spread(s,dim=1,ncopies=3)*u # 3xns
    sout=np.tile(s[:,np.newaxis],(1,nr)).T*u
    rscout = np.tile(rsc[np.newaxis,:],(ns,1)).T
    rout=np.tile(rsc[np.newaxis,:],(ns,1)).T + sout #*u
    # Now we need an ECEF->LATLONH conversion
    out =  ecef2latlon(rout)

    # ! The below essentially use terms from the ECEF -> ENU conversion on the surface of an oblate spheroid,
    # ! our WGS-84 ellipsoid
    # ! The normal to the ellipsoid has the normal snorm
    snorm = np.zeros((nr,ns))
    snorm[0,:] =np.cos(out[0,:])*np.cos(out[1,:])
    snorm[1,:] =np.sin(out[0,:])*np.cos(out[1,:])
    snorm[2,:] =np.sin(out[1,:])
    # The cos(view zenith) is formed by the pointing vector from the ground to the spacecraft
    # dotted into the surface normal; this is (for one location) sum(-u(:,i)*snorm(:,i))
    view_zen=np.arccos(np.sum(-u*snorm,axis=0))*rad2Deg # deg from zenith

    uDotE=np.sin(out[0,:])*u[0,:] - np.cos(out[0,:])*u[1,:]
    uDotN=np.cos(out[0,:])*np.sin(out[1,:])*u[0,:] + np.sin(out[0,:])*np.sin(out[1,:])*u[1,:] - \
        np.cos(out[1,:])*u[2,:]
    view_az=np.arctan2(uDotE,uDotN)*rad2Deg # deg from N, clockwise
    out=out*rad2Deg

    return out,view_zen,view_az

def solar_geometry(iyr,imn,idy,xh,xm,xs,xlat,xlong):
    #     ! 2012-10-31 MJM Made able to handle an array of xlat and xlong, where xlat and xlong
    #     !                are each vectors
    #     ! 2002-12-11 MJM Turned into array valued function. Now should only
    #     !                necessary solar geometry calculations. This is
    #     !                called many nlines times when given appropriate info for
    #     !                each line.
    #     ! xlat: positive NORTH
    #     ! xlon: positive WEST
    #     ! 2001-05-16 MJM Solar factors only.

    sol_zen_az = np.zeros((2,len(xlat)))


    tt=2*np.pi*((xh)/24.+xm/1440.+xs/86400.)

    xlatr = xlat * deg2Rad
    xlongr = xlong * deg2Rad

    # at this TIME there is only ONE astronomical solar position
    dec,haz = suncor(idy,imn,iyr,tt)
    solaz,el = hazel(haz+tt-xlongr,dec,xlatr)

    #!!!C---Note: DEC, SOLAZ,and EL DEC are in radians

    # !    IF (EL .LE. 0.) THEN
    # !       print*,IYR,Imn,idy,IH,IM,IS,xlat,xlong
    # !       print*,'el=',el
    # !       call fatal_error(file_name,module_name,subroutine_name,&
    # !            &'ERROR: Sun is below the horizon!!!'//&
    # !            &' Check input date, time, latitude and longitude.')
    # !    ENDIF

    solzni = np.pi/2.0 - el
    sol_zen_az[0,:]=solzni #! radians
    sol_zen_az[1,:]=solaz  #!radians

    return sol_zen_az

def suncor(iday,month,iyr,tz):

    jd=julian(iday,month,iyr)
    fjd=0.5 + tz/(2.0*np.pi)
    ras,dec,gsdt,bzero,pzero,solong = solcor(jd,fjd)
    haz=gsdt-ras-tz

    return dec,haz

def julian(iday,month,iyr):

    md=[0,31,59,90,120,151,181,212,243,273,304,334]

    jyr=iyr-1600
    I1=np.int(jyr/400)               #Number of complete 400 yr spans; each have 97 leap years
    I2=np.int((jyr-400*I1)/100)      #Number of 100 year spans, taking away the 400 year spans above
    I3=np.int((jyr-400*I1-100*I2)/4) #Number of four year spans, taking away complete 400 and 100 above
                                    #    JD(AD1600) normal yr  400 yr  100yr  leap yr this century
    jday=2305447 + 365*jyr + 97*I1 + 24*I2 +I3 #Counts IYR if IYR a leap year

    jday=jday+md[month-1]+iday

    if month == 2:
        leap = 0
        if isLeapYear(jyr):
            leap = 1

        #Already counted leap day above, so take out if jan/feb
        jday = jday - leap
    return jday

def hazel(h,d,xlat):
    sne = np.sin(d)*np.sin(xlat)+np.cos(d)*np.cos(xlat)*np.cos(h)
    e=np.arcsin(sne)
    sna = np.cos(d)*np.sin(h)
    csa=(np.sin(xlat)*np.cos(h)*np.cos(d)-np.sin(d)*np.cos(xlat))
    a=np.arctan2(sna,csa)+np.pi

    return a,e

def solcor(jd,fjd):

    jyr=365.25
    d=(jd-2415020)+fjd
    iyr=np.int(d/jyr)

    g=-.026601523+.01720196977*d-1.95e-15*d*d -2*np.pi*iyr
    xlms=4.881627938+.017202791266*d+3.95e-15*d*d-2*np.pi*iyr
    obl=.409319747-6.2179e-9*d
    ecc=.01675104-1.1444e-9*d
    f=d-jyr*iyr
    dgsdt=1.739935476+2*np.pi*f/jyr+1.342027e-4*d/jyr
    gsdt=dgsdt+2*np.pi*(fjd-0.5)
    xlts=xlms+2.*ecc*np.sin(g)+1.25*ecc*ecc*np.sin(2.*g)
    sndc=np.sin(xlts)*np.sin(obl)
    ddecs=np.arcsin(sndc)
    decs=ddecs
    csra=np.cos(xlts)/np.cos(ddecs)
    ras=np.arccos(csra)
    if np.sin(xlts) < 0.:
        ras=2*np.pi-ras
    omega=1.297906+6.66992e-7*d
    thetac=xlts-omega
    bzro=np.arcsin(.126199*np.sin(thetac))
    p=-np.arctan(np.cos(xlts)*np.tan(obl))-np.arctan(.127216*np.cos(thetac))
    xlmm=np.arctan2(.992005*np.sin(thetac),np.cos(thetac))
    jdr=jd-2398220
    frot=(jdr+fjd)/25.38-np.int((jdr+fjd)/25.38)
    solong=xlmm-2*np.pi*frot+np.pi-3.e-4
    if solong <  0.:
        solong += 2*np.pi
    elif solong > 2*np.pi:
        solong -= 2*np.pi

    return ras,decs,gsdt,bzro,p,solong

def frac2DegMinSec(fracd,fracm,fracs):
#     ! A subroutine to get degrees, min, sec from
#     ! fractional pieces. The output are all real, although a
#     ! strong case can be made for degrees and minutes to be
#     ! integers.
#     ! MJM Unknown date
    d_all = fracd+fracm/60.+fracs/3600.
    d_hold = np.floor(fracd+fracm/60.+fracs/3600.)
    m_all = (d_all-d_hold)*60.
    m_hold = np.floor(m_all)

    return [d_hold,m_hold,(m_all-m_hold)*60.]
