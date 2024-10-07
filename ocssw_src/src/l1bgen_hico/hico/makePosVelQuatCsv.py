from datetime import datetime as DT
from datetime import timedelta as TDel
import pickle
from astropy.time import Time
from scipy.interpolate import UnivariateSpline as UVS
import pandas as pd
import socket
import sys
import numpy as np
from hico.auxiliary import ConvLst2DT
from hico.exceptions import USGNCCoarseTimeException
from hico.exceptions import PVQInterpolation
from hico.exceptions import EmptyCSVException
from hico.exceptions import BadFormatCSVException
import os
import getpass
import platform
import re

__version__ = "0.1"
__author__ = "R. Healy & E. Karakoylu (erdem.m.karakoylu@nasa.gov)"
__date__ = "2/3/2017"


class MakePosVelQuatCSV:

    def __init__(self, hicoPtr, **kwargs):
        gps_time = DT(1980, 1, 6, 0, 0, 0)
        sec_per_count = 1.117460459e-6
        secOffset = TDel(seconds=15)
        self.__subSeconds2Seconds = 16.72e-6
        self.__timeWrapAround = 65536
        self.__time0Factor = 0.05
        self.__time1Factor = 62.5e-9
        do_nav_time_correction = kwargs.pop('doNavTimeCorrection', False)
        self.beginTime = DT.now()
        self.paramsDict = {}
        self.paramsDict['d_tisspvq'] = kwargs.pop('delta_tisspvq', 0)
        self.paramsDict['d_ticugps'] = kwargs.pop('delta_ticugps', 0)
        self.paramsDict['orientation'] = hicoPtr.params['issOrientation']
        self.paramsDict['gps_secs_090101'] = (DT(2009, 1, 1, 0, 0, 0) + secOffset
                                              - gps_time).total_seconds()
        self.paramsDict['trigger_pulse_width'] = 0.860e-3
        self.paramsDict['start_date_time'] = ConvLst2DT(hicoPtr.L0.data['start_date'],
                                                        hicoPtr.L0.data['start_time'])
        self.paramsDict['end_date_time'] = ConvLst2DT(hicoPtr.L0.data['end_date'],
                                                      hicoPtr.L0.data['end_time'])
        self.paramsDict['start_date_time'] -= secOffset
        self.paramsDict['end_date_time'] += secOffset
        self.paramsDict['thetas'] = hicoPtr.L0.data['ScenePointingAngle']
        self.paramsDict['nls'] = hicoPtr.L0.data['n_image']
        ffpps_all = hicoPtr.L0.header['FFpps'] +\
            hicoPtr.L0.header['FFppsSub'] * self.__subSeconds2Seconds
        lfpps_all = hicoPtr.L0.header['LFpps'] +\
            hicoPtr.L0.header['LFppsSub'] * self.__subSeconds2Seconds
        if ffpps_all > lfpps_all:
            ffpps_all += self.__timeWrapAround
        self.paramsDict['ffpps_all'] = ffpps_all
        self.paramsDict['lfpps_all'] = lfpps_all
        self.paramsDict['exptimes'] = hicoPtr.L0.header['TriggerCount'] * sec_per_count
        self.paramsDict['odrc_time_offset'] = self.__GetODRCTimeOffset(hicoPtr.inpDict['hdrFile'])
        self.paramsDict['nav_time_offset'] = 0
        self.paramsDict['jday_end'] = Time(self.paramsDict['end_date_time']).jd
        self.paramsDict['jday_start'] = Time(self.paramsDict['start_date_time']).jd
        self.paramsDict['exptimes'] = hicoPtr.L0.header['TriggerCount'] * sec_per_count
        self.paramsDict['cexp'] = '_expDEF'
        rootName = os.path.basename(hicoPtr.inpDict['l0File']).split('.bil')[0]
        self.paramsDict['pvqFileName'] = kwargs.pop('outputFileName',
                                                    '%s_pos_vel_quat.csv' % rootName)
        self.paramsDict['anglFileName'] = '%s_LonLatViewAngles.bil' % rootName
        self.paramsDict['csvName'] = hicoPtr.inpDict['csvFile']
        self.paramsDict['n_pixels'] = hicoPtr.L0.data['ns']
        self.__ReadCSVFile()
        if do_nav_time_correction:
            self.__GetNavOffsetTimeCorrection(rootName)
        self.__FillPoVeQuDF()
        self.finishTime = DT.now()
        self.__WriteHeader2CSV()
        self.__WriteData2CSV()

    def __ReadCSVFile(self):
        usecols = ['USGNC_PS_Pointing_Coarse_Time_Tag',
                   'USGNC_PS_Pointing_Inert_Posn_VectorX',
                   'USGNC_PS_Pointing_Inert_Posn_VectorY',
                   'USGNC_PS_Pointing_Inert_Posn_VectorZ',
                   'USGNC_PS_Pointing_Inert_Vel_VectorX',
                   'USGNC_PS_Pointing_Inert_Vel_VectorY',
                   'USGNC_PS_Pointing_Inert_Vel_VectorZ',
                   'USGNC_PS_Pointing_Current_Inert_Att_Quatrn_0',
                   'USGNC_PS_Pointing_Current_Inert_Att_Quatrn_1',
                   'USGNC_PS_Pointing_Current_Inert_Att_Quatrn_2',
                   'USGNC_PS_Pointing_Current_Inert_Att_Quatrn_3',
                   'USGNC_PS_PD_Fine_Pointing_Fine_Time_Tag',
                   'HSTCLOCKTIME0', 'HSTATTITUDEQUATS',
                   'HSTATTITUDESTATUSMODE',
                   'HSTATTITUDEQUATX', 'HSTATTITUDEQUATY', 'HSTATTITUDEQUATZ',
                   'HSTATTITUDETIME0', 'HSTATTITUDETIME1', '_ISSGPSTIME',
                   'ISSTIMESUBSECOND', 'ICUTIMEGPSSECONDS',
                   'ICUTIMEISSSUBSECOND', 'ICUTIMEHWREGSECOND',
                   'ICUTIMEHWREGSUBSECOND', 'ISSTIMECENTURY',
                   'ISSTIMEYEAR', 'ISSTIMEMONTH', 'ISSTIMEDAY', 'ISSTIMEHOUR',
                   'ISSTIMEMINUTE', 'ISSTIMESECOND']

        try:
            self.dfCSV = pd.read_csv(self.paramsDict['csvName'],
                                     usecols=usecols)
        except ValueError:
            print("badly formatted csv file - attempting to resolve")
            ocdataroot = os.getenv('OCDATAROOT')
            hicocaldir = os.path.join(ocdataroot, 'hico/cal')
            col_name_file = os.path.join(hicocaldir, 'csv_columns.pkl')
            try:
                with open(col_name_file, 'rb') as f:
                    col_names = pickle.load(f)
                try:
                    self.dfCSV = pd.read_csv(self.paramsDict['csvName'],
                                             skiprows=1, names=col_names,
                                             usecols=usecols)
                    print('csv format issue resolved')
                except ValueError:
                    raise BadFormatCSVException
            except FileNotFoundError:
                raise BadFormatCSVException
        else:
            if self.dfCSV.empty:
                raise EmptyCSVException


    def __WriteHeader2CSV(self):
        with open(self.paramsDict['pvqFileName'], 'w') as fpvq:
            print('\nThis file name, %s' % os.path.basename(self.paramsDict['pvqFileName']),
                  file=fpvq)
            print('Source CSV file, %s' % self.paramsDict['csvName'], file=fpvq)
            print('Subset CSV file, None', file=fpvq)
            print('Expected Lon/Lat/View angle filename, %s' % self.paramsDict['anglFileName'],
                  file=fpvq)
            print('Epoch, 2009 Jan 01 00:00:00 UTC', file=fpvq)
            print('Requested number of pixels, %d' % self.paramsDict['n_pixels'], file=fpvq)
            print('Distance Unit, Feet', file=fpvq)
            print('Central Body, Earth', file=fpvq)
            print('CoordinateSystem, J2000', file=fpvq)
            print('Theta (degrees from stowed position), %.16f' % self.paramsDict['thetas'],
                  file=fpvq)
            print('Code name, %s' % __name__, file=fpvq)
            print('Code version, %s' % __version__, file=fpvq)
            print('Code date, %s' % __date__, file=fpvq)
            print('Code author, %s' % __author__, file=fpvq)
            print('Code executed on computer, %s' % socket.gethostname(), file=fpvq)
            print('Code executed by username, %s' % getpass.getuser(), file=fpvq)
            print('Code run under Python osfamily, %s' % os.name, file=fpvq)
            print('Code run under Python os, %s' % platform.release(), file=fpvq)
            print('Code start time, %s' % self.beginTime.strftime("%a %b %d %H:%M:%S"), file=fpvq)
            print('Code end time, %s' % self.finishTime.strftime("%a %b %d %H:%M:%S"), file=fpvq)
            print('Exposure interval (frame time), %f' % self.paramsDict['exptimes'], file=fpvq)
            print('ISS orientation, %s' % self.paramsDict['orientation'], file=fpvq)
            print('Trigger pulse width (s), %.3e' % self.paramsDict['trigger_pulse_width'],
                  file=fpvq)
            print('ODRC broadcast time - gps time (s), %s' % self.paramsDict['odrc_time_offset'],
                  file=fpvq)
            print('delta_ticugps (s), %.6f' % 0, file=fpvq)
            print('delta_tisspvq (s), %.6f' % 0, file=fpvq)
            print(file=fpvq)

    def __WriteData2CSV(self):
        with open(self.paramsDict['pvqFileName'], 'a') as fpvq:
            self.dfPVQ.to_csv(fpvq, index=False)

    def __GetPVQIdx(self, relTimes):
        crstmtagfld = 'USGNC_PS_Pointing_Coarse_Time_Tag'
        fntmtagfld = 'USGNC_PS_PD_Fine_Pointing_Fine_Time_Tag'
        locs_goodt = ((self.dfCSV[crstmtagfld] >= relTimes[0] - 10) &
                      (self.dfCSV[crstmtagfld] <= relTimes[-1] + 10))
        idx = locs_goodt[locs_goodt].index
        if idx.size == 0:
            raise(USGNCCoarseTimeException)
        u_usgnc = self.dfCSV.loc[idx, crstmtagfld].duplicated()
        udx_usgnc = u_usgnc.index[np.logical_not(u_usgnc)]
        u_usgnc_coarse = self.dfCSV.loc[udx_usgnc, crstmtagfld].values
        u_usgnc_fine = self.dfCSV.loc[udx_usgnc, fntmtagfld].values
        t_issposvelquat = u_usgnc_coarse + u_usgnc_fine / 256 + self.paramsDict['d_tisspvq']
        return udx_usgnc, t_issposvelquat

    def __GetStarTracking(self, hicoTimes, hstField):
        hstattitudetime = self.dfCSV.HSTATTITUDETIME0 * self.__time0Factor +\
            self.dfCSV.HSTATTITUDETIME1 * self.__time1Factor
        hstclocktime = self.dfCSV.HSTCLOCKTIME0 * self.__time0Factor +\
            self.dfCSV.HSTATTITUDETIME1 * self.__time1Factor
        hstqt_pps = self.__GetIssPosVelQuatXYZ(hstclocktime, hstattitudetime, hstField)
        t_issgps = self.dfCSV._ISSGPSTIME + 1e-6 * self.dfCSV.ISSTIMESUBSECOND -\
            self.paramsDict['odrc_time_offset'] + self.paramsDict['nav_time_offset']
        interpolator = self.__PVQInterpolator(t_issgps, hstqt_pps)
        return interpolator(hicoTimes)

    def __SplineIfNeeded(self, hicoRelTimes, fieldName):
        if self.dfCSV[fieldName].unique().size > 1:
            return self.__GetStarTracking(hicoRelTimes, fieldName)
        else:
            return np.ones((self.dfPVQ.shape[0],),
                           dtype=self.dfCSV[fieldName].dtype) * self.dfCSV[fieldName].unique()

    def __GetNavOffsetTimeCorrection(self, rootname):
        '''distance units are in km, time in s, speed in km/s'''
        hico_ground_speed = 6.9 #km/s
        rootnameparts = rootname.split('.')
        scene_ref = 'H%s%s' % (rootnameparts[1], rootnameparts[3])
        nav_offset_file = os.path.join(os.getenv('OCVARROOT'), 'hico', 'navigation_offsets.txt')
        regex_pattern = re.compile(r'H\d+\s(-*\d*\.\d+|-*\d+)\s*(-*\d*\.\d+|-*\d+)*')
        with open(nav_offset_file) as nof:
            for line in nof:
                if scene_ref in line:
                    linematch = regex_pattern.findall(line)
                    if linematch:
                        along_track = float(linematch[0][0])
                        across_track = float(linematch[0][1])
                        self.paramsDict['nav_offset_km'] = {'along_track': along_track,
                                                            'across_track': across_track,
                                                            'scene': scene_ref}
                        self.paramsDict['nav_time_offset'] = along_track / hico_ground_speed
                        # TODO add logging
                    break



    def __FillPoVeQuDF(self):
        fields = ['SecondsSinceEpoch', 'ISSPOSX', 'ISSPOSY', 'ISSPOSZ', 'ISSVELX', 'ISSVELY',
                  'ISSVELZ', 'ISSQX', 'ISSQY', 'ISSQZ', 'ISSQS', 'STQX', 'STQY', 'STQZ', 'STQS',
                  'HST_ATTITUDE_STATUS']
        vecPosXfld = 'USGNC_PS_Pointing_Inert_Posn_VectorX'
        vecPosYfld = 'USGNC_PS_Pointing_Inert_Posn_VectorY'
        vecPosZfld = 'USGNC_PS_Pointing_Inert_Posn_VectorZ'
        vecVelXfld = 'USGNC_PS_Pointing_Inert_Vel_VectorX'
        vecVelYfld = 'USGNC_PS_Pointing_Inert_Vel_VectorY'
        vecVelZfld = 'USGNC_PS_Pointing_Inert_Vel_VectorZ'
        vecQuat0fld = 'USGNC_PS_Pointing_Current_Inert_Att_Quatrn_0'
        vecQuat1fld = 'USGNC_PS_Pointing_Current_Inert_Att_Quatrn_1'
        vecQuat2fld = 'USGNC_PS_Pointing_Current_Inert_Att_Quatrn_2'
        vecQuat3fld = 'USGNC_PS_Pointing_Current_Inert_Att_Quatrn_3'
        hstAttXfld = 'HSTATTITUDEQUATX'
        hstAttYfld = 'HSTATTITUDEQUATY'
        hstAttZfld = 'HSTATTITUDEQUATZ'
        hstAttSfld = 'HSTATTITUDEQUATS'
        hstAttStatfld = 'HSTATTITUDESTATUSMODE'
        self.dfPVQ = pd.DataFrame(columns=fields)
        hicoRelTimes = self.__GetSecsSinceEpoch()
        self.dfPVQ['SecondsSinceEpoch'] = hicoRelTimes - self.paramsDict['gps_secs_090101']
        pvIdx, pvTimes = self.__GetPVQIdx(hicoRelTimes)
        self.dfPVQ['ISSPOSX'] = self.__GetIssPosVelQuatXYZ(hicoRelTimes, pvTimes, vecPosXfld, pvIdx)
        self.dfPVQ['ISSPOSY'] = self.__GetIssPosVelQuatXYZ(hicoRelTimes, pvTimes, vecPosYfld, pvIdx)
        self.dfPVQ['ISSPOSZ'] = self.__GetIssPosVelQuatXYZ(hicoRelTimes, pvTimes, vecPosZfld, pvIdx)
        self.dfPVQ['ISSVELX'] = self.__GetIssPosVelQuatXYZ(hicoRelTimes, pvTimes, vecVelXfld, pvIdx)
        self.dfPVQ['ISSVELY'] = self.__GetIssPosVelQuatXYZ(hicoRelTimes, pvTimes, vecVelYfld, pvIdx)
        self.dfPVQ['ISSVELZ'] = self.__GetIssPosVelQuatXYZ(hicoRelTimes, pvTimes, vecVelZfld, pvIdx)
        self.dfPVQ['ISSQX'] = self.__GetIssPosVelQuatXYZ(hicoRelTimes, pvTimes, vecQuat1fld, pvIdx)
        self.dfPVQ['ISSQY'] = self.__GetIssPosVelQuatXYZ(hicoRelTimes, pvTimes, vecQuat2fld, pvIdx)
        self.dfPVQ['ISSQZ'] = self.__GetIssPosVelQuatXYZ(hicoRelTimes, pvTimes, vecQuat3fld, pvIdx)
        self.dfPVQ['ISSQS'] = self.__GetIssPosVelQuatXYZ(hicoRelTimes, pvTimes, vecQuat0fld, pvIdx)
        self.dfPVQ['STQX'] = self.__SplineIfNeeded(hicoRelTimes, hstAttXfld)
        self.dfPVQ['STQY'] = self.__SplineIfNeeded(hicoRelTimes, hstAttYfld)
        self.dfPVQ['STQZ'] = self.__SplineIfNeeded(hicoRelTimes, hstAttZfld)
        self.dfPVQ['STQS'] = self.__SplineIfNeeded(hicoRelTimes, hstAttSfld)
        self.dfPVQ['HST_ATTITUDE_STATUS'] = self.__SplineIfNeeded(hicoRelTimes, hstAttStatfld)

    @staticmethod
    def __PVQInterpolator(timeInGrid, pvqOutGrid):
        ditp = pd.DataFrame(dict(time=timeInGrid, pvq=pvqOutGrid))
        ditp.drop_duplicates(subset='time', inplace=True)
        ditp.sort_values('time', inplace=True)
        return UVS(ditp.time.values, ditp.pvq.values)

    @staticmethod
    def __GetODRCTimeOffset(hdrFilePath):
        try:
            with open(hdrFilePath, 'rt') as fhdr:
                for line in fhdr.readlines():
                    if 'odrc_time_offset' in line:
                        return float(line.split('=')[1])
            return 0
        except FileNotFoundError as e:
            sys.exit(e)


    def __ExtractTimeDataFromCSV(self):
        """
        Extracts columns for the rest of processing
        replaces half-baked date and time columns with proper datetime
        column.
        """
        def GetDateTime(row):
            row = row.astype(int)
            year = row['ISSTIMECENTURY'] * 100 + row['ISSTIMEYEAR']
            return DT(year, row['ISSTIMEMONTH'], row['ISSTIMEDAY'],
                      row['ISSTIMEHOUR'], row['ISSTIMEMINUTE'], row['ISSTIMESECOND'])

        hsdat = pd.DataFrame(columns=['datetime', ])
        cols2copy = ['ICUTIMEGPSSECONDS', 'ICUTIMEISSSUBSECOND',
                     'ICUTIMEHWREGSECOND', 'ICUTIMEHWREGSUBSECOND']
        hsdat = self.dfCSV[cols2copy].copy()
        #hsdat['datetime'] = self.dfCSV.apply(GetDateTime, axis=1)
        hsdat['datetime'] = self.dfCSV.filter(regex='ISSTIME',
                                              axis=1).apply(GetDateTime,
                                                            axis=1)
        try:
            np.testing.assert_array_equal(hsdat.datetime.values, np.sort(hsdat.datetime.values))
        except AssertionError:
            print("Caution: datetime column unordered")
        finally:
            return hsdat

    def __GetSecsSinceEpoch(self):
        def CheckStuff(self, hsdat, hwreg):
            badhwregs = (hsdat.ICUTIMEHWREGSECOND[0] > hsdat.ICUTIMEHWREGSECOND).values
            if badhwregs.any():
                hsdat.loc[badhwregs, 'ICUTIMEHWREGSECOND'] += self.__timeWrapAround
            if self.paramsDict['lfpps_all'] < self.paramsDict['ffpps_all']:
                self.paramsDict['ffpps_all'] += self.__timeWrapAround
            if self.paramsDict['lfpps_all'] < hwreg[0]:
                self.paramsDict['lfpps_all'] += self.__timeWrapAround
            if self.paramsDict['ffpps_all'] < hwreg[0]:
                self.paramsDict['ffpps_all'] += self.__timeWrapAround

        def GetHwReg(hsdat, ffpps_all):
            hwreg = hsdat.ICUTIMEHWREGSECOND +\
                    hsdat.ICUTIMEHWREGSUBSECOND * 16.72e-6
            locs_lt = hwreg <= ffpps_all
            idx_lt = locs_lt.index[locs_lt]
            # the first of these is the tick just after the camera stops
            # hwreg_locs_gt = hwreg[hwreg >= lfpps_all]
            return hwreg, idx_lt

        hsdat = self.__ExtractTimeDataFromCSV()
        hwreg, idx_lt = GetHwReg(hsdat, self.paramsDict['ffpps_all'])

        d_ticugps = 0
        t_icugps = hsdat.loc[:, 'ICUTIMEGPSSECONDS'].values +\
            1.e-6 * hsdat.loc[:, 'ICUTIMEISSSUBSECOND'].values

        t_icugps -= self.paramsDict['odrc_time_offset']
        t_icugps += d_ticugps  # 2012/01/17 for testing
        hwstart = (self.paramsDict['ffpps_all'] - hwreg[idx_lt[-2]]
                   ) / (hwreg[idx_lt[-1]] - hwreg[idx_lt[-2]])
        time_start = hwstart*(t_icugps[idx_lt[-1]] -
                              t_icugps[idx_lt[-2]]) + t_icugps[idx_lt[-2]]
        scanrange = np.arange(self.paramsDict['nls'], dtype=float)
        testrange = np.multiply(scanrange, self.paramsDict['exptimes'])
        hico_times = np.add(np.add(time_start, testrange),
                            (self.paramsDict['trigger_pulse_width'] + 0.5 *
                             self.paramsDict['exptimes']))
        return hico_times

    def __GetIssPosVelQuatXYZ(self, hicoTimes, pvtimes, inputField, pvidx=np.zeros(1)):
        if pvidx.any():
            pvGrid = self.dfCSV.loc[pvidx, inputField].values
        else:
            pvGrid = self.dfCSV[inputField].values
        try:
            f = self.__PVQInterpolator(pvtimes, pvGrid)
        except ValueError:
            raise(PVQInterpolation(inputField))
        issPvOut = f(hicoTimes)
        return issPvOut

    def __GetSTQXYZS(self):
        pass

    def __GetHSTAttitudeStat(self):
        pass
