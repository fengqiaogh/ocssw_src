# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 14:01:27 2015
Reads Hico L0 bil file.
Header (first 256 bytes) parsed into header dictionary.
Data
@author: EKarako:ylu:
"""
import sys
import struct
import numpy as np
import datetime as dtime
from math import sqrt, atan, atan2, pi
import logging
import time


class HicoL0Reader(object):

    def __init__(self, l0FileName, parentLoggerName):

        self.data = {}
        self.header = {}
        self.logger = logging.getLogger('%s.HicoL0Reader' % parentLoggerName)
        with open(l0FileName, 'rb') as fh:
            self.__ParseFillHeaderDict(fh)
            self.__ParseFillDataDict(fh)

    @staticmethod
    def __BinaryCodedDec2Int(x):
        return((x & 0xf0) >> 4) * 10 + (x & 0x0f)

    def __HicoDateTime(self):
        def __timeWrapAround(*args):
            hr, mn, secs = (arg for arg in args)
            mn += secs // 60
            secs = secs % 60
            hr += mn // 60
            mn = mn%60
            return [hr, mn, secs]
        t_start = time.time()
        # Might have to add a safety to ensure date/time data is in allowable range.
        centerDate = []
        centerTime = []
        # dateList.append(self.__BinaryCodedDec2Int(self.header['FFyearMSB']) * 100 + \
        # self.__BinaryCodedDec2Int(self.header['FFyearLSB']))
        year = self.header['FFyearMSB'] * 100 + self.header['FFyearLSB']
        centerDate.append(year)
        # dateList.append(self.__BinaryCodedDec2Int(self.header['FFmonth']))
        centerDate.append(self.header['FFmonth'])
        # dateList.append(self.__BinaryCodedDec2Int(self.header['FFday']))
        centerDate.append(self.header['FFday'])
        # timeList.append(self.__BinaryCodedDec2Int(self.header['FFhours']))
        centerTime.append(self.header['FFhours'])
        # timeList.append(self.__BinaryCodedDec2Int(self.header['FFmin']))
        centerTime.append(self.header['FFmin'])

        # pps wraps at 2**16
        if self.header['LFpps'] >= self.header['FFpps']:
            LFpps_ = self.header['LFpps']
        else:
            LFpps_ = (self.header['LFpps'] + 2**16)
        time_imageInt = LFpps_ + self.header['LFppsSub'] * 16.762e-6 - \
            (self.header['FFpps'] + self.header['FFppsSub'] * 16.762e-6)

        # secs = self.__BinaryCodedDec2Int(self.header['FFsec']) + \
        secs = self.header['FFsec'] + \
            ((int(self.header['Word08'] & 0x0f) << 16) +
             int(self.header['FFsubLSB'])) * 1.0e-6 + 0.5 * time_imageInt
        centerTime.append(secs)
        startDate = centerDate[:]
        endDate = centerDate[:]
        startTime = centerTime[:]
        endTime = centerTime[:]
        startTime[2] -= 0.5 * time_imageInt
        endTime[2] += 0.5 * time_imageInt
        startTime = __timeWrapAround(*startTime)
        endTime = __timeWrapAround(*endTime)

        tempDate = dtime.date(centerDate[0], centerDate[1], centerDate[2])
        iday = tempDate.timetuple().tm_yday
        fhr = centerTime[0] + (centerTime[1] + centerTime[2] / 60) / 60
        resDict = {'startTime': startTime, 'centerTime': centerTime,
                   'endTime': endTime, 'startDate': startDate, 'centerDate': centerDate,
                   'endDate': endDate, 'iday': iday, 'fhr': fhr}
        tot_time = time.time() - t_start
        self.logger.debug('time taken: %f' % tot_time)
        return resDict

    @staticmethod
    def __LatLonH(*arg):

        if len(arg) == 1:
            xyz, = arg
            x, y, z = xyz
        elif len(arg) == 3:
            x, y, z = arg
        else:
            print("!!-> Incorrect arg. # @ LatLonH")
            sys.exit(1)
        esa_m = 6378137.0  # earth semimajor axis (radius?) in meters
        recFlat = 1 / 298.257223563  # reciprocal flattening
        seminax = esa_m * (1 - recFlat)  # semi-minor axis
        fEcc_2 = 2 * recFlat * (1 - recFlat)  # first eccentricity, squared
        sEcc_2 = recFlat * (2 - recFlat) / ((1 - recFlat) ** 2)  # second eccentricity, squared
        r2 = x**2 + y**2
        r = sqrt(r2)
        ee2 = esa_m**2 - seminax**2
        ff = 54 * seminax ** 2 * z**2
        g = r2 + (1 - fEcc_2) * z ** 2 - fEcc_2*ee2
        c = (fEcc_2**2 * ff * r2) / (g**3)
        s = 1.0 + c + sqrt(c * (c + 2)) ** (1.0/3.0)
        p = ff / (2.0 * (s + 1.0 / (s) + 1.0)**2 * g**2)
        q = sqrt(1.0 + 2 * fEcc_2**2 * p)
        ro = -(fEcc_2*p*r) / (1+q) + sqrt((esa_m**2 / 2.0) * (1+1.0/q) -
                                          ((1-fEcc_2) * p * z**2) / (q * (1+q)) - p*r2/2.0)
        tmp = (r - fEcc_2 * ro)**2
        u = sqrt(tmp+z**2)
        v = sqrt(tmp+(1-fEcc_2)*z**2)
        z0 = (seminax**2 * z) / (esa_m * v)
        h = u*(1 - seminax ** 2 / (esa_m * v))
        phi = atan((z + sEcc_2*z0)/r)*180.0/pi
        lambd = atan2(y, x) * 180.0 / pi
        returnList = [phi, lambd, h]
        return returnList

    def __ParseFillHeaderDict(self, fHandle):
        """
        Goes through 256-byte embedded HICO header data.
        Note that header is big-endian while data is little-endian
        """
        t_start = time.time()
        keys = ['Sync', 'Len', 'FFpre', 'FFyearMSB', 'FFyearLSB', 'FFmonth', 'FFday',
                'FFhours', 'FFmin', 'FFsec', 'Word08', 'FFsubLSB', 'FFpps', 'FFppsSub',
                'LFpps', 'LFppsSub', 'roiw', 'roih', 'roix', 'roiy', 'hbin', 'vbin',
                'ROport', 'ROspeed', 'IcuSwVersion', 'CamClearMode', 'TotalFrames',
                'Dark1', 'Scene', 'Dark2', 'ID', 'ExpTime']
        values = struct.unpack('>4sL8b10H6b5HL', fHandle.read(56))  # 33 items
        d1 = dict(zip(keys, values))
        if d1['Sync'] != b'HICO':
            sys.exit("File is not L0. Exiting...")
        keys2 = ['FFpos', 'FFvel', 'FFquat']
        FFpos = np.fromfile(fHandle, dtype='>f4', count=3)  # big endian float array
        FFvel = np.fromfile(fHandle, dtype='>f4', count=3)
        FFquat = np.fromfile(fHandle, dtype='>f4', count=4)

        keys3 = ['FFgncSec', 'FFgncSubSec', 'FFgncQual', 'TotalErr']
        values3 = struct.unpack('>L2bH', fHandle.read(8))

        keys4 = ['LFpos', 'LFvel', 'LFquat']
        LFpos = np.fromfile(fHandle, dtype='>f4', count=3)  # big endian float array
        LFvel = np.fromfile(fHandle, dtype='>f4', count=3)
        LFquat = np.fromfile(fHandle, dtype='>f4', count=4)

        keys5 = ['LFgncSec', 'LFgncSubSec', 'LFgncQual', 'EventFlags', 'Dark1Offset',
                 'SceneOffset', 'Dark2Offset', 'Again', 't_on', 't_dark', 't_obs', 't_ang',
                 'TriggerCount', 'EMgain', 'Dark1Tel', 'Dark1Pos', 'Dark1TriggerCount',
                 'SceneTel', 'ScenePos', 'SceneTriggerCount', 'Dark2Tel', 'Dark2Pos',
                 'Dark2TriggerCount', 't_motor', 'Spare100', 'Spare101', 'Spare102',
                 'ErrorLogCount']
        values5 = struct.unpack('>L2bH4L7Hh2Hh2Hh6H', fHandle.read(64))

        keys6 = ['Errors']
        Errors = np.fromfile(fHandle, dtype='>i2', count=24)
        d2 = dict(zip(keys2, (FFpos, FFvel, FFquat)))
        d3 = dict(zip(keys3, values3))
        d4 = dict(zip(keys4, (LFpos, LFvel, LFquat)))
        d5 = dict(zip(keys5, values5))
        d6 = dict(zip(keys6, Errors))
        self.header.update(d1)
        self.header.update(d2)
        self.header.update(d3)
        self.header.update(d4)
        self.header.update(d5)
        self.header.update(d6)

        # BCD correction:
        self.header['FFpre'] = self.__BinaryCodedDec2Int(self.header['FFpre'])
        self.header['FFyearMSB'] = self.__BinaryCodedDec2Int(self.header['FFyearMSB'])
        self.header['FFyearLSB'] = self.__BinaryCodedDec2Int(self.header['FFyearLSB'])
        self.header['FFmonth'] = self.__BinaryCodedDec2Int(self.header['FFmonth'])
        self.header['FFday'] = self.__BinaryCodedDec2Int(self.header['FFday'])
        self.header['FFhours'] = self.__BinaryCodedDec2Int(self.header['FFhours'])
        self.header['FFmin'] = self.__BinaryCodedDec2Int(self.header['FFmin'])
        self.header['FFsec'] = self.__BinaryCodedDec2Int(self.header['FFsec'])
        tot_time = time.time() - t_start
        self.logger.debug('time taken: %f' % tot_time)
        return None

    def __ParseFillDataDict(self, fHandle):
        # this is to fill the L0 dictionary (hash)
        t_start = time.time()
        m2ft = 0.3048
        dtDict = self.__HicoDateTime()
        ffLatLonH = self.__LatLonH(self.header['FFpos'] * m2ft)
        lfLatLonH = self.__LatLonH(self.header['LFpos'] * m2ft)
        ffLatLonH[2] /= 1e3  # converts m to KM.
        lfLatLonH[2] /= 1e3
        nl = int(self.header['TotalFrames'])
        nb = int(self.header['roih'])
        ns = int(self.header['roiw'])
        count = nl * nb * ns
        pic0 = np.reshape(np.fromfile(fHandle, dtype='<i2', count=count),
                          (nl, nb, ns))
        self.data = {"ns": ns, "nb": nb, "nl": nl, "offset": 2561,
                     "file_type": 'ENVI Standard', "data_type": 12, "interleave": 'bil',
                     "sensor_type": self.header['Sync'], "byte order": 0,
                     "exposure_time": self.header['TriggerCount'] * 1.117460459,
                     "n_predark": int(self.header['Dark1']),
                     "n_image": int(self.header['Scene']),
                     "n_postdark": int(self.header['Dark2']),
                     "start_date": dtDict['startDate'], "start_time": dtDict['startTime'],
                     "center_date": dtDict['centerDate'], "center_time": dtDict['centerTime'],
                     "end_date": dtDict['endDate'], "end_time": dtDict['endTime'],
                     "yearDay": dtDict['iday'], "floatHour": dtDict['fhr'],
                     "ScenePointingAngle": (self.header['ScenePos'] + 395) * 0.004,
                     "FFVelMag": sqrt(sum((self.header['FFvel'] * m2ft - 3) ** 2)),
                     "LFVelMag": sqrt(sum((self.header['LFvel'] * m2ft - 3) ** 2)),
                     "FFLatLonH": ffLatLonH,
                     "LFLatLonH": lfLatLonH,
                     "FFquat": self.header["FFquat"],
                     "LFquat": self.header["LFquat"],
                     "Dark1Pos": self.header["Dark1Pos"],
                     "Dark2Pos": self.header["Dark2Pos"],
                     "pic0": pic0
                     }
        tot_time = time.time() - t_start
        self.logger.debug('time taken: %f' % tot_time)
        return None
