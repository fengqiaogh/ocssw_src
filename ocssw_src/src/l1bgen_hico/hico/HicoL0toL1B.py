#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 14:32:11 2015
Transaltion of NRL's IDL code hico_raw_to_radiance.pro
@author: Erdem Karako:ylu:

NOTES ON THE IDL=>PYTHON TRANSLATION
    There a number of substantial dangling assignments (variables assigned to that are not
    subsequently used -- this is apparently the cleaned-up version.)

    Lines-1333->1359: Nothing useful, a couple of assignments to file_ and geo_history.
        ->deal with that later.

"""
import os
import sys
import time
import functools
import math as m
import calendar as cal
import logging
import pickle

import pandas as pd
from datetime import datetime as DT
from datetime import timedelta as TDel
import numpy as np
from scipy import interpolate

from hico.HicoL0Reader import HicoL0Reader



def TimeIt(func):
    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        start = time.time()
        result = func(self, *args, **kwargs)
        tot_time = time.time() - start
        self.logger.debug('time taken by %s: %.2f secs.' % (func.__name__, tot_time))
        return result
    return wrapper

def LoopDrkLine(data, startLine, stopLine):
    '''
    This is here because I was trying to jit it with numba. tbc...
    '''
    win_size = 5
    win_n_sd = 3.0
    for drk_line in range(startLine, stopLine+1):
        wmin = max(startLine, drk_line - win_size)
        wmax = min(stopLine, drk_line + win_size)
        nwin = wmax - wmin + 1
        win_ind = np.arange(nwin) + wmin
        win_ind = win_ind[np.where(win_ind != drk_line)]
        win_mat = data[win_ind]
        win_sd = win_mat.std(axis=0, dtype='f4', ddof=1)
        offset = win_n_sd * win_sd
        win_median = np.median(win_mat, axis=0).astype('f4')
        minCrit = win_median - offset.astype('f4')
        maxCrit = win_median + offset.astype('f4')
        data[drk_line] = np.where((data[drk_line] < minCrit) |
                                  (data[drk_line] > maxCrit),
                                  win_median, data[drk_line])


class HicoL0toL1b:
    """
    Object to assemble the data for converting a HICO L0 file to
    L1b format.
    Optional keyword arguments:
    hdrFile: header file name
      => default: <L0 basename>.hdr
    csvFile: csv file name
      => default: <L0 basename>.csv
    outFile: output file name
      => default: <L0 basename>.nc
    bandScaleFactorFile: as the name suggests
      => default: Gao_scaling_Final_1_2012.txt
    hicoFitCoeffsFile: radiometric calibration file name
      => default: coeff_080706g10emncexp0137ls4l_it13_ftt111_smoothed_2ndcorrsm074x1_asm10_bsm20_
                  wsm12_iss_1_2012.dat
    wvlFile: wavelength file name
      => default: "HICO_Wavlns_FWHM_128_Chnls_Gao_Shift_0p84nm_v2p0_1column.txt"
    issXVVFile: ISS data file name
      => default: "ISS_MINUS_XVV.TXT"
    secOrFile: second order correction coefficients file name
      => default: "HICO_Empirical_water_scale_ftt1.11ms_1_2012.txt"
    r2rInputParamFile: processing parameters, defaults will be applied if not
      supplied
    doCal: calibration flag
      => default = true.

    """
    @TimeIt
    def __init__(self, iArgs, **kwds):
      ocdataroot = os.environ['OCDATAROOT']
      ocvarroot = os.environ['OCVARROOT']
      mainLoggerName = kwds.pop('parentLoggerName', __name__)
      self.logger = logging.getLogger('%s.HicoL0toL1B' % mainLoggerName)
      self.inpDict = {}
      self.fNamesDict = {}
      self.inpDict["l0File"] = iArgs.l0file
      self.inpDict["hdrFile"] = iArgs.hdrfile
      self.inpDict["csvFile"] = iArgs.csvfile
      if "bandScaleFactorFile" in kwds:
          self.inpDict["bandScaleFactorFile"] = kwds["bandScaleFactorFile"]
      else:
          self.inpDict["bandScaleFactorFile"] = ocdataroot + \
                       "/hico/cal/Gao_scaling_factor_Final_1_2012.txt"
      if "hicoFitCoeffsFile" in kwds:
          self.inpDict["hicoFitCoeffsFile"] = kwds["hicoFitCoeffsFile"]
      else:
          self.inpDict["hicoFitCoeffsFile"] = ocdataroot + \
              "/hico/cal/coeff_080706g10emncexp0137ls4l_" \
              + "it13_ftt111_smoothed_2ndcorrsm074x1_" + \
                "asm10_bsm20_wsm12_iss_1_2012.dat"
      if "wvlFile" in kwds:
          self.inpDict["wvlFile"] = kwds["wvlFile"]
      else:
          self.inpDict["wvlFile"] = ocdataroot + \
              "/hico/cal/HICO_Wavlns_FWHM_128_Chnls_Gao" \
              + "_Shift_0p84nm_v2p0_1column.txt"
      if "secOrFile" in kwds:
          self.inpDict["secOrFile"] = kwds["secOrFile"]
      else:
          self.inpDict["secOrFile"] = ocdataroot + \
              "/hico/cal/HICO_Empirical_water_scale_" \
              + "ftt1.11ms_1_2012.txt"
      if "issXVVFile" in kwds:
          self.inpDict["issXVVFile"] = kwds["issXVVFile"]
      else:
          self.inpDict["issXVVFile"] = ocvarroot+"/hico/ISS_MINUS_XVV.TXT"
      self.inpDict["r2rInputParam"] = kwds.pop("r2rInputParam", None)
      self.inpDict["doCal"] = kwds.pop("doCal", True)
      self.inpDict["verbose"] = kwds.pop("verbose", False)
      self.inpDict["labDark"] = kwds.pop("labDark", False)
      self.inpDict["doSmooth"] = kwds.pop("doSmooth", True)
      self.logger.info('inpDict content: ')
      for k, v in self.inpDict.items():
          self.logger.info('%s: %s' % (k, v))
      # populate HicoBil object
      self.params = dict.fromkeys(['base_dir', 'input_dir', 'output_dir',
                                   'dark_temp', 'ffTemp', 'flip_left_right',
                                   'flip_cal', 'cal_mult', 'cal_mult',
                                   'cal_shift', 'outputInterleave',
                                   'outputScaleFactor', 'frameTransferTime',
                                   'exposureTime', 'validL0Pixels', 'nb_read',
                                   'nb_all', 'ns', 'nl_predark',
                                   'nl_postdark', 'output_is_int',
                                   'jpg_r_range', 'jpg_g_range',
                                   'jpg_b_range', 'wave_um', 'fwhm',
                                   'coeff_13ms', 'bScaleFactor',
                                   'begin_date', 'end_date', 'begin_time', 'end_time'])
      self.debugDict = dict.fromkeys(['corr_smear', 'so2', 'smoother'])
      self.geomDict = dict.fromkeys(['imCenLat', 'imCenLon', 'fullLats',
                                     'fullLons', 'viewz', 'viewa',
                                     'sunposzen', 'sunposaz'])
      self.headerDict = dict.fromkeys(['ns', 'nl', 'nb', 'offset', 'inter',
                                       'data_type', 'byte_order',
                                       'bad_bands', 'file_type',
                                       'sensor_type', 'band_names',
                                       'descrip', 'wave_len', 'fwhm',
                                       'wavelength_unit', 'x_start',
                                       'ystart', 'map_info', 'def_bands',
                                       'z_range', 'imcenlat', 'imcenlon'])
      self.badDataDict = dict.fromkeys(['bad_pixs', 'sat_pixs'])
      self.outFilNamesDict = dict.fromkeys(['geomEnvi', 'geomHdr', 'flagEnvi',
                                            'flagHdr', 'ndviEnvi',
                                            'ndviHdr', 'rgbJpg', 'rgbBip',
                                            'rgbHdr', 'outRad', 'outRadHdr'])
      self.ancDict = dict.fromkeys(['n11', 'n12', 'n1t',
                                    'n21', 'n22', 'n2t',
                                    'n31', 'n32', 'n3t',
                                    'd2_ratio'])  # , 'rhot_red', 'rhot_nir'])
      self.L0 = HicoL0Reader(iArgs.l0file, self.logger.name)
      self.__FillOutputFileNames()
      self.__ReadWaveFile()
      self.__SetL02L1BParams()
      self.__GetISSOrientation()
      if self.inpDict['doCal']:
          self.__ReadCalFile()
      self.__ReadBandScaleFile()
      self.__FillAncDict()
      self.__SetPeriods()
      self.__FillGeomDict()
      n21 = self.ancDict['n21']
      n22p = self.ancDict['n22'] + 1
      self.floatData = self.L0.data['pic0'][n21:n22p].astype('f4')

    @TimeIt
    def __FillAncDict(self):
      self.logger.debug('filling ancDict')
      self.ancDict['n11'] = 3
      self.ancDict['n12'] = self.L0.data['n_predark'] - 1
      self.ancDict['n1t'] = self.ancDict['n12'] - self.ancDict['n11'] + 1
      self.ancDict['n20'] = self.ancDict['n12'] + 1  # start of HICO L1B image
      # start of non zeroed out  HICO L1B image
      self.ancDict['n21'] = self.L0.data['n_predark'] + self.ancDict['n11']
      self.ancDict['n2t'] = self.L0.data['nl'] - self.L0.data['n_predark'] -\
          self.L0.data['n_postdark'] - self.ancDict['n11']
      # end of HICO L1B image
      self.ancDict['n22'] = self.ancDict['n2t'] + self.ancDict['n21'] - 1
      self.ancDict['n31'] = self.L0.data['nl'] - self.L0.data['n_postdark']\
          + self.ancDict['n11']
      self.ancDict['n32'] = self.L0.data['nl'] - 1
      self.ancDict['n3t'] = self.ancDict['n32'] - self.ancDict['n31'] + 1
      theta_day = (2 * m.pi / 365) * (self.L0.data['yearDay'] - 1)
      d2_ratio = (1.00011 + 0.034221 * m.cos(theta_day) +
                  0.001280 * m.sin(theta_day) +
                  0.000719 * m.cos(2 * theta_day) -
                  0.000077 * m.sin(2 * theta_day))
      self.ancDict['d2_ratio'] = d2_ratio

    @TimeIt
    def __FillOutputFileNames(self):
      '''
      Fills dictionary containing output filenames
      '''
      if self.inpDict['doCal']:
          radtype = '_rad'
      else:
          radtype = '_nocal'
      outFileNameBase = os.path.basename(self.inpDict['l0File']) + radtype
      self.outFilNamesDict['geomEnvi'] = outFileNameBase + '_geom.bil'
      self.outFilNamesDict['geomHdr'] = outFileNameBase + '_geom.hdr'
      self.outFilNamesDict['outRad'] = outFileNameBase + '.bil'
      self.outFilNamesDict['outRadHdr'] = outFileNameBase + '.hdr'
      self.outFilNamesDict['flagEnvi'] = outFileNameBase + '_flag.bip'
      self.outFilNamesDict['flagHdr'] = outFileNameBase + '_flag.hdr'
      self.outFilNamesDict['ndviEnvi'] = outFileNameBase + '_ndvi.bip'
      self.outFilNamesDict['ndviHdr'] = outFileNameBase + '_ndvi.hdr'
      self.outFilNamesDict['rgbBip'] = outFileNameBase + '_rgb.bip'
      self.outFilNamesDict['rgbHdr'] = outFileNameBase + '_rgb.hdr'
      self.outFilNamesDict['rgbJpg'] = outFileNameBase + '_rgb.jpg'

    @TimeIt
    def __GetISSOrientation(self):
      """
         Determines ISS orientation
      """
      flipLR = True
      issOrientation = "+XVV"
      sceneid = int(os.path.basename(self.inpDict["l0File"]).split('.')[5])
      with open(self.inpDict["issXVVFile"]) as isf:
          for line in isf.readlines()[1:]:
              lElems = line.split(',')
              if sceneid >= int(lElems[0]) and sceneid <= int(lElems[1]):
                  flipLR = False
                  issOrientation = "-XVV"
      self.params['flip_left_right'] = flipLR
      self.params['issOrientation'] = issOrientation

    def __SetPeriods(self):
      '''Formatting time for writing to nc file'''
      datefmt = '%Y%m%d'
      timefmt = '%H%M%S'
      gnc_epoch = DT(1980, 1, 6)
      start_delta = TDel(0, self.L0.header['FFgncSec'])
      end_delta = TDel(0, self.L0.header['LFgncSec'])
      firstFrameDT = start_delta + gnc_epoch
      lastFrameDT = end_delta + gnc_epoch
      self.params['begin_date'] = DT.strftime(firstFrameDT, datefmt)
      self.params['end_date'] = DT.strftime(lastFrameDT, datefmt)
      self.params['begin_time'] = DT.strftime(firstFrameDT, timefmt)
      self.params['end_time'] = DT.strftime(lastFrameDT, timefmt)

    @TimeIt
    def __SetL02L1BParams(self):
      """
      Populates L02L1BParams hash (as a dict).
      In hico_L0_to_L1B_h5.bash this is equivalent to
          * filling out {basename}.raw_param
      """
      offset = -300
      ff_time = DT(self.L0.header['FFyearMSB'] * 100 + self.L0.header['FFyearLSB'],
                   self.L0.header['FFmonth'], self.L0.header['FFday'],
                   self.L0.header['FFhours'], self.L0.header['FFmin'], self.L0.header['FFsec'])
      timestamp = cal.timegm(ff_time.timetuple())
      dark_time = DT.utcfromtimestamp(timestamp + offset)
      ff_temp = self.__GetCamTemp(self.inpDict["csvFile"], ff_time)
      dark_temp = self.__GetCamTemp(self.inpDict["csvFile"], dark_time)
      self.params['base_dir'] = "./"
      self.params['input_dir'] = "./"
      self.params['output_dir'] = "./"
      self.params['dark_temp'] = dark_temp
      self.params['ffTemp'] = ff_temp
      self.params['flip_cal'] = True
      self.params['cal_mult'] = 1.32
      self.params['cal_shift'] = 0
      self.params['outputInterleave'] = "bil"
      self.params['outputScaleFactor'] = 50  # out_scale
      self.params['frameTransferTime'] = 0.00111
      self.params['exposureTime'] = 0.0137
      self.params['validL0Pixels'] = 508
      self.params['nb_read'] = 128  # nb
      self.params['nb_all'] = 170  # nbb
      self.params['ns'] = 512
      self.params['nl_predark'] = 200
      self.params['nl_postdark'] = 200
      self.params['out_type'] = 'uint16'
      self.__MakeFWHM()

    # %% STATICMETHOD HELPER FUNCTIONS
    @staticmethod
    def GetLineAveragedPic(tempData):
      good = np.where(tempData != -1, True, False)
      picXave = tempData[good].reshape(tempData.shape).mean(axis=0,
                                                            dtype='f4')
      return picXave

    @staticmethod
    def __LinearCongrid(arr, newdims):
      method = 'linear'
      minusone = True
      centre = False
      m1 = np.cast[int](minusone)
      ndims = len(arr.shape)
      dimlist = []
      ofs = np.cast[int](centre) * 0.5
      old = np.array(arr.shape)
      newdims = np.asarray(newdims, dtype=float)
      # calculate new dims
      for i in range(ndims):
          base = np.arange(newdims[i])
          dimlist.append((old[i] - m1) / (newdims[i] - m1) *
                         (base + ofs) - ofs)
      # specify old dims
      olddims = [np.arange(i, dtype=np.float) for i in list(arr.shape)]

      # first interpolation - for ndims = any
      mint = interpolate.interp1d(olddims[-1], arr, kind=method)
      newa = mint(dimlist[-1])
      # trorder = [ndims - 1] + range(ndims - 1)
      for i in range(ndims - 2, -1, -1):
          newa = newa.transpose()  # trorder)
          mint = interpolate.interp1d(olddims[i], newa, kind=method)
          newa = mint(dimlist[i])
      if ndims > 1:
          # need one more transpose to return to original dimensions
          newa = newa.transpose()  # trorder)
      return newa

    @staticmethod
    def __RotMatrix(roll, pitch, yaw):
      """
      Computes a rotation matrix given roll,pitch and yaw, in that order.
      Resulting rotation matrix shape is (3x3)
      """
      rotMat = np.zeros((3, 3))
      phi = m.radians(roll)
      theta = m.radians(pitch)
      psi = m.radians(yaw)
      ct = np.cos(theta)
      st = np.sin(theta)
      cp = np.cos(psi)
      sp = np.sin(psi)
      cf = np.cos(phi)
      sf = np.sin(phi)
      rotMat = [[ct * cp, ct * sp, -st],
                [sf * st * cp - cf * sp, sf * st * sp + cf * cp, ct * sf],
                [cf * st * cp + sf * sp, cf * st * sp - sf * cp, ct * cf]]

      return np.array(rotMat)

    @staticmethod
    def __RngBear2LatLonEllipsoid(latInArr, lonInArr, distArr, brngInArr):
      """
      Solution of the geodetic direct problem after T. Vincenty.
      Modified Rainsford's method with Helmert's elliptical terms,
      effective in any azimuth and at any distance short of antipodal.

      semiMajAx is the semi-major axis of the reference ellipsoid;
      flatFactor is the flattening of the reference elllipsoid;
      latitudes and longitudes in randians, +ve North and East;
      azims in radians clockwise from North.
      Geodesic distance s assumed in units of semiMajAx
      """
      # WGS84 constants (www.wgs84.com)
      flatFactor = 1/298.257223563
      semiMajAx = 6378137

      # Initial Values
      eps = 0.5e-13   # tolerance
      glon1 = np.radians(lonInArr)
      glat1 = np.radians(latInArr)
      faz = np.radians(brngInArr)
      r = 1 - flatFactor
      tu = r*np.sin(glat1) / np.cos(glat1)
      sf = np.sin(faz)
      cf = np.cos(faz)
      baz = np.zeros_like(faz)
      w = cf != 0
      baz = np.arctan2(tu[w], cf[w])*2
      cu = (np.sqrt(tu * tu + 1)) ** -1
      su = tu*cu
      sa = cu*sf
      c2a = -sa * sa + 1
      x = np.sqrt((1/(r**2) - 1) * c2a + 1) + 1
      x = (x - 2) / x
      c = 1 - x
      c = ((x*x/4) + 1) / c
      d = (0.375 * x**3 - x)
      tu = distArr / r / semiMajAx / c
      y = tu

      while True:
          sy = np.sin(y)
          cy = np.cos(y)
          cz = np.cos(baz + y)
          e = 2 * cz ** 2 - 1
          c = y
          x = e * cy
          y = 2 * e - 1
          y = ((((2*sy)**2 - 3) * y * cz * d / 6 + x) * d / 4 - cz) * sy *\
              d + tu
          if np.max(abs(y - c)) <= eps:
              break
      baz = cu * cy * cf - su * sy
      c = r * np.sqrt(sa**2 + baz ** 2)
      d = su * cy + cy * sy * cf
      glat2 = np.arctan2(d, c)
      c = cu * cy - su * sy * sf
      x = np.arctan2(sy * sf, c)
      c = ((-3 * c2a + 4) * flatFactor + 4) * c2a * flatFactor / 16
      d = ((e * cy * c + cz) * sy * c + y) * sa
      glon2 = glon1 + x - (1 - c) * d * flatFactor
      baz = np.arctan2(sa, baz) + m.pi
      lonOutArr = np.degrees(glon2)
      latOutArr = np.degrees(glat2)
      brngOutArr = np.degrees(baz)

      return latOutArr, lonOutArr, brngOutArr

    @staticmethod
    def __Quats2Euler(q, qtype=0):
      '''
      Takes ISS Quaternions (q) in their normal order.
      Returns angles phi, theta, psi.
      X,Y,Z are forward in, right of orbit and down, respectively.
      qtype0 returns Euler angle in the typical R1 R2 R3 order w/ first
      rotation about z-axis, then y-axis, then x-axis

      '''
      if qtype == 0:
          phi = m.degrees(m.atan2(2*(q[2] * q[3] + q[0] * q[1]),
                          q[3]**2 - q[2]**2 - q[1]**2 + q[0]**2))
          theta = m.degrees(-m.asin(2 * (q[1] * q[3] - q[0] * q[2])))
          psi = m.degrees(m.atan2(2 * (q[2] * q[1] + q[0] * q[3]),
                          -q[3]**2 - q[2]**2 + q[1]**2 + q[0]**2))
      else:
          phi = m.degrees(m.atan2(2*(-q[1] * q[2] + q[0] * q[3]),
                          -q[3]**2 - q[2]**2 + q[1]**2 + q[0]**2))
          theta = m.degrees(m.asin(2 * (q[1] * q[3] + q[0] * q[2])))
          psi = m.degrees(m.atan2(2 * (-q[2] * q[3] + q[0] * q[1]),
                          q[3]**2 - q[2]**2 - q[1]**2 + q[0]**2))
      return phi, theta, psi

    @staticmethod
    def __Vec2ZenithBearing(vec):
      """"
      Input: vec[x,y,z]
      Output: zenith and bearing angles (radians)

      """
      zenithRad = m.acos(vec[2])
      bearingRad = m.atan2(vec[1], vec[0])
      if bearingRad < 0:
          bearingRad += 2 * m.pi
      return zenithRad, bearingRad

    @staticmethod
    def __GetVecs(rotMat, angleRad):
      """
      Output 3x1 ref vector.
      Input: * 3x3 rotation matrix, rotMat
             * reference angle in radians, angleRad
      """
      dirVec = np.array([0, m.sin(angleRad), m.cos(angleRad)])
      refVec = np.array([np.sum(rotMat[:, i] * dirVec) for i in range(3)])
      return refVec

    @staticmethod
    def __GetDist(zenithAngleRad, height):
      """
      Computes distance along surface of ellipsoid (km)
      """
      earthRadiusKM = 6378
      dist = (m.asin(((earthRadiusKM + height) / earthRadiusKM) *
              m.sin(abs(zenithAngleRad))) -
              abs(zenithAngleRad)) * earthRadiusKM
      return dist

    @staticmethod
    def __ValueLocate(vec, vals):
      '''
      Equivalent to IDL's value_locate routine. Excerpt from exelis follows
          "The VALUE_LOCATE function finds the intervals within a given monotonic
          vector that brackets a given set of one or more search values."
      Assumes vec and val are 1D arrays
      '''
      return np.amax([np.where(v >= vec, np.arange(len(vec)), -1)
                      for v in vals], axis=1)

    def __GetCamTemp(self, csvFileName, timeData):
      """
      READS CSV FILE
      """
      use_columns = ["ISSTIMEYEAR", "ISSTIMEMONTH", "ISSTIMEDAY",
                     "ISSTIMEHOUR", "ISSTIMEMINUTE", "ISSTIMESECOND",
                     "HICOCAMERATEMP"]
      badval = -1
      try:
          df = pd.read_csv(csvFileName, usecols=use_columns)
      except ValueError:
          self.logger.warning("badly formatted csv file - attempting to resolve")
          ocdataroot = os.getenv('OCDATAROOT')
          hicocaldir = os.path.join(ocdataroot, 'hico/cal')
          col_name_file = os.path.join(hicocaldir, 'csv_columns.pkl')
          try:
              with open(col_name_file, 'rb') as f:
                  col_names = pickle.load(f)
          except FileNotFoundError:
              self.logger.error("column name file not found")
              sys.exit(status=1)
          try:
              df = pd.read_csv(csvFileName, skiprows=1, names=col_names,
                               usecols=use_columns)
              self.logger.info('csv format issue resolved')
          except ValueError:
              self.logger.error("Something is wrong with the CSV file")
              sys.exit(status=1)

      x = df.loc[(df.ISSTIMEYEAR == (timeData.year % 100)) &
                 (df.ISSTIMEMONTH == (timeData.month)) &
                 (df.ISSTIMEDAY == (timeData.day)) &
                 (df.ISSTIMEHOUR == (timeData.hour)) &
                 (df.ISSTIMEMINUTE == (timeData.minute)) &
                 (df.ISSTIMESECOND == (timeData.second))].HICOCAMERATEMP.values
      if x.size == 1:
          return x[0]
      else:
          return badval

    @TimeIt
    def __ReadWaveFile(self):
      '''Reads wavelength file'''
      if self.inpDict['wvlFile'] == 'DEFAULT':
          nb = self.L0.data['nb']
          if nb == 128:
              wc1 = 0.353428
              wc2 = 0.005728
          else:
              wc1 = 0.3515095
              wc2 = 0.0019095
          self.inpDict['wvlFile'] = 'lambda_b (0<=b<=' + str(nb) + '):' +\
              str(wc1) + '+' + str(wc2) + '*b'
          self.params['wave_um'] = wc1 + wc2 * np.ones(nb)
      else:
          self.params['wave_um'] = np.loadtxt(self.inpDict['wvlFile'])

    @TimeIt
    def __ReadCalFile(self):
      """
      READS hicoFitCoeffs FILE
      """
      if sys.byteorder == 'little':
          dt = np.dtype('<f')
      else:
          dt = np.dtype('>f')
      tempArray = np.fromfile(self.inpDict['hicoFitCoeffsFile'], dtype=dt)
      coeff_13ms = tempArray.reshape(self.params['nb_read'], self.params['ns'])
      if self.L0.data['nb'] == 128:
          exposure_time = self.L0.data['exposure_time'] / 1e6
          et_ratio = 13e-3 / exposure_time
          coeff_13ms *= et_ratio
          cal_sh = self.params['cal_shift']
          if self.params['flip_cal']:
              coeff_13ms = np.fliplr(coeff_13ms)
          if cal_sh != 0:
              hhold = coeff_13ms
              if cal_sh < 0:
                  coeff_13ms[:, 0:512 + cal_sh] = hhold[:, 0-cal_sh:512]
                  coeff_13ms[:, 512 + cal_sh:512] = 0
              else:
                  coeff_13ms[:, cal_sh:512] = hhold[0:512 - cal_sh, :]
                  coeff_13ms[:, 0:cal_sh] = 0.0
      coeff_13ms = np.expand_dims(coeff_13ms, axis=0)
      self.params['coeff_13ms'] = coeff_13ms

    @TimeIt
    def __ReadBandScaleFile(self):
      """
      READS bandScaleFacor FILE
      """
      bScaleFactor = np.genfromtxt(self.inpDict['bandScaleFactorFile'])
      bScaleFactor = np.expand_dims(bScaleFactor[:, 1], axis=0)
      bScaleFactor = np.expand_dims(bScaleFactor, axis=-1)
      self.params['bScaleFactor'] = bScaleFactor

    @TimeIt
    def __MSALGam(self):
      """
      Returns angle of slit on the ground, returned in degrees
      """
      msAlt = (self.L0.data['FFLatLonH'][2] + self.L0.data['LFLatLonH'][2]
               ) * 0.5
      msLatDeg = (self.L0.data['FFLatLonH'][0] +
                  self.L0.data['LFLatLonH'][0]) * 0.5
      msLatDeg = min([max([msLatDeg, -51.6]), 51.6])
      gammaRad = m.acos(m.cos(m.radians(51.6)) / m.cos(m.radians(msLatDeg)))
      if self.L0.data['LFLatLonH'][0] < self.L0.data['FFLatLonH'][0]:
          gammaRad *= -1
      gammaDeg = m.degrees(gammaRad)
      return msAlt, msLatDeg, gammaDeg

    @TimeIt
    def _DoDarkCorrection(self):
      """
      hico_to_rad_kwp_truncate.pro
      lines1760 - 1932
      Nested function implemented to avoid nest loop repetitions (thanks NRL!)
      No closure: nested function not returned by calling function.
      Also pic0 & others are  already in scope but note "nonlocal" declaration
      allowing changes in pic0.
      NOTE: THIS IS USED ONLY BY SMEAR CORRECTION FUNCTION
      """
      # nl, ns, nb = self.L0.data['nl'], self.L0.data['ns'], self.L0.data['nb']
      n11, n12 = self.ancDict['n11'], self.ancDict['n12']
      n21, n22 = self.ancDict['n21'], self.ancDict['n22']
      n1t, n2t = self.ancDict['n1t'], self.ancDict['n2t']
      n31, n32 = self.ancDict['n31'], self.ancDict['n32']
      # theData = self.theData
      intData = self.L0.data['pic0']
      dark1Pos, dark2Pos = self.L0.data['Dark1Pos'], self.L0.data['Dark2Pos']
      darkTemp, fftemp = self.params['dark_temp'], self.params['ffTemp']
      prePostDarkOK, pic13ave, ts, = 0, 0, 41
      fit20 = np.empty((n2t, 1, 1))
      fit20[:, 0, 0] = np.log(1 + np.arange(n2t) / ts)
      f13av = (1 + ts / n1t) * np.log(1 + n1t / ts) - 1
      dark_b_coef = 11.2
      if fftemp >= 10.0 and fftemp < 24.0:
          dark_b_coef = 11.15
      elif fftemp >= 24.0 and fftemp < 27.5:
          dark_b_coef = -0.6286 * fftemp + 26.2857
      elif fftemp >= 27.5 and fftemp < 33.0:
          dark_b_coef = 0.7273 * fftemp - 11.0
      elif fftemp >= 33.0 and fftemp < 50.0:
          dark_b_coef = 0.2143 * fftemp + 5.9286
      LoopDrkLine(intData, startLine=n11, stopLine=n12)  # modifies theData
      self.logger.debug('intData after n11-: %s' % intData.dtype)
      LoopDrkLine(intData, startLine=n31, stopLine=n32)  # modifies theData
      self.logger.debug('intData after n31-: %s' % intData.dtype)
      if dark1Pos == -395:
          pic13ave += self.GetLineAveragedPic(intData[n11:n12+1])
          prePostDarkOK += 1
      if dark2Pos == -395:
          pic13ave += self.GetLineAveragedPic(intData[n31:n32+1])
          prePostDarkOK += 2
      if prePostDarkOK == 3:
          pic13ave /= 2
      elif prePostDarkOK == 2:
          pic13ave *= (0.970 + 0.0036 * (darkTemp - 23))
      elif prePostDarkOK == 1:
          pic13ave *= (1.032 - 0.0040 * (darkTemp - 23))
      b = dark_b_coef + 0.9 * (pic13ave - 221) / (285 - 221)
      # self.floatData += pic13ave - b * f13av + 1.2
      # a2 = pic13ave - b * f13av + 1.2
      self.floatData = intData[n21:n22+1] - (pic13ave - b * f13av + 1.2 + b * fit20.astype('f4'))
      # self.floatData += b[np.newaxis, :, :] * fit20.astype('f4')
      self.logger.debug('floatData after drk corr.: %s' % self.floatData.dtype)

    @TimeIt
    def _DoSmearCorrection(self):
      nb = self.L0.data['nb']
      ns = self.L0.data['ns']
      exposureTime = self.params['exposureTime']
      frameTransferTime = self.params['frameTransferTime']
      stable_integration_time = np.array(exposureTime - frameTransferTime, dtype='f2')
      delta_t = np.array(frameTransferTime / 511, dtype='f2')
      frac_eqB6 = (frameTransferTime + delta_t) / (stable_integration_time - delta_t)
      if nb == 128:
          tot_count = (self.floatData.sum(axis=1) +
                       self.floatData[:, -1, :] * (171 - nb)) * \
                      (3/ns) * frac_eqB6
      else:
          tot_count = (self.floatData.sum(axis=1) +
                       self.floatData[:, -1, :] * (512 - nb)) * \
                      (1/ns) * frac_eqB6
      self.floatData *= (1 + frac_eqB6.astype('f2'))
      self.floatData -= tot_count[:, np.newaxis, :].astype('f2')
      self.logger.debug('floatData after smear correction: %s' % self.floatData.dtype)

    @TimeIt
    def __MakeFWHM(self):
      nb = self.L0.data['nb']
      fwhm = np.empty(nb)
      if self.inpDict['doSmooth']:
          fwhm[:3] = 0.005
          fwhm[3:69] = 0.01
          fwhm[69:nb-4] = 0.02
          fwhm[nb-4:] = 0.005
      else:
          if nb == 128:
              fwc = 0.005
          elif nb == 384:
              fwc = 0.00167
          fwhm = np.ones(nb) * fwc
      self.params['fwhm'] = fwhm * 1000

    @TimeIt
    def _DoGaussianSmoothing(self):
      ''' Two steps   1-prepare parameters for gaussian smoothing
                      2-perform gaussian smoothing
      '''
      nb = self.L0.data['nb']
      n2t = self.ancDict['n2t']
      # prepare gaussian smoothing
      fwhm_10nm = np.array([0.0001497, 0.0141444, 0.2166629, 0.5380859])
      fwhm_10nm = np.concatenate((fwhm_10nm, -np.sort(-fwhm_10nm[:-1])))
      fwhm_20nm = np.array([0.0070725, 0.0347489, 0.1083365, 0.2143261, 0.2690555])
      fwhm_20nm = np.concatenate((fwhm_20nm, -np.sort(-fwhm_20nm[:-1])))
      # normalize
      fwhm_10nm /= fwhm_10nm.sum()
      fwhm_20nm /= fwhm_20nm.sum()
      smoother = np.zeros((nb, nb))
      np.fill_diagonal(smoother, 1)
      for i in range(3, 69):
          smoother[i, i-3:i+4] = fwhm_10nm  # diff indexing order than idl
      for i in range(69, nb-4):
          smoother[i, i-4:i+5] = fwhm_20nm  # diff indexing order than idl
      for li in range(n2t):
          self.floatData[li] = smoother.dot(self.floatData[li])
      self.logger.debug("floatData: %s" % self.floatData.dtype)

    @TimeIt
    def _Do2ndOrdCorrection(self):
      """
      Populates a matrix that holds linear interp. weights combined with scale factors
      to subsequently allow interpolation and 2nd order calcs. to be done in a single
      matrix multiplication.
      """
      nb = self.L0.data['nb']
      n2t = self.ancDict['n2t']
      # n21 = self.ancDict['n21']
      # n22 = self.ancDict['n22']
      wl = self.params['wave_um']
      wave_half = wl / 2
      sc_1 = 0.131 * (0.86 - wl)
      sc_2 = 0.113 * (wl - 0.86)
      scale = np.where(sc_1 > sc_2, sc_1, sc_2)
      scale[wl < 2 * wl[0]] = 0
      maP = self.__ValueLocate(wl, wave_half)
      test = maP != -1
      s2data = np.loadtxt(self.inpDict['secOrFile'], skiprows=1)
      if nb == 128:
          scale = s2data[:, 2]
      else:
          scale = np.interp(wl, s2data[:, 0], s2data[:, 2])
      so2 = np.diagflat([1]*nb).astype('f8')
      for i in range(nb):
          if test[i]:
              weight = (wave_half[i]-wl[maP[i]])/(wl[maP[i] + 1] - wl[maP[i]])
              so2[i, maP[i]:maP[i]+2] = -np.array([(1 - weight), weight]) * scale[i]
      if so2.sum() != so2.shape[1]:
          for li in range(n2t):
              # self.floatData[li] = so2.dot(self.floatData[li])
              self.floatData[li] = so2.dot(self.floatData[li])

    @TimeIt
    def _ApplyMasks(self):
      '''Called by _DoCalMaskAndScale()'''
      ns = self.L0.data['ns']
      bad_sat_pix = self.badDataDict['bad_pixs'][3:] + self.badDataDict['sat_pixs'][3:]
      if bad_sat_pix.any():
          self.floatData[bad_sat_pix] = 0
      #if self.params['validL0Pixels'][0] > 0:
          #left_mask = self.params['validL0Pixels'][0]
          #self.floatData[:, :, :left_mask] = 0
      #if self.params['validL0Pixels'] < (ns - 1):
      right_mask = self.params['validL0Pixels']
      self.floatData[:, :, right_mask:ns] = 0

    @TimeIt
    def _FlipDataArray(self):
      ''' hack to use numpy's fliplr method on a 3D array. '''
      self.floatData = np.rollaxis(self.floatData, 0, 3)
      self.floatData = np.fliplr(self.floatData)
      self.floatData = np.rollaxis(self.floatData, 2, 0)

    @TimeIt
    def _DoCalMaskAndScale(self):
      '''Calibration changes depending on number of bands and doCal option.
      When doCal=False but nb=384, data is re-scaled.
      '''
      cal_mult = self.params['cal_mult']
      gao_scale = self.params['bScaleFactor']
      out_scale = self.params['outputScaleFactor']
      self.floatData *= self.params['coeff_13ms']
      self._ApplyMasks()
      if self.params['flip_left_right']:
          self._FlipDataArray()
      self.floatData *= cal_mult
      self.floatData *= gao_scale
      return None

    @TimeIt
    def ConvertRaw2Rad(self):
      self._DoDarkCorrection()  # do darkSubtraction
      self.logger.info("Dark correction done")
      self._StoreBadData()
      self.logger.info("Bad data indexed")
      self._DoSmearCorrection()
      self.logger.info("Smear correction done")
      self._Do2ndOrdCorrection()
      self.logger.info("Second order correction done")
      if self.inpDict['doSmooth']:
          self._DoGaussianSmoothing()
          self.logger.info("Gaussian smoothing done")
      if self.inpDict['doCal'] and self.L0.data['nb'] == 128:
          self._DoCalMaskAndScale()
          self.logger.info("Calibration/Masking/Scaling done")
      elif not self.inpDict['doCal'] and self.L0.data['nb'] == 384:
          self.floatData /= self.params['outputScaleFactor']
      self.logger.info("Raw to radiance conversion completed")

    # GEOM Helper files
    @TimeIt
    def __GetSunAngles(self):
      """
      Inputs
      -------
          iday: year's day number
          hr: GMT time of day in real hours, e.g. 4:30PM = 16.5
      Returns
      -------
      solarAngles {ndarray}, where
          solarAngles[0] = solar zenith angle in degrees
          solarAngles[1] = solar azimuth angle in degrees clockwise from North
      """
      iday, hr = self.L0.data['yearDay'], self.L0.data['floatHour']
      xlon, ylat = self.geomDict['fullLons'], self.geomDict['fullLats']
      # Compute solar declination angle
      thez = 360*(iday - 1) / 365
      rthez = m.radians(thez)
      sdec = 0.396372 - 22.91327 * m.cos(rthez) + 4.02543 * m.sin(rthez) - \
          0.387205 * m.cos(2 * rthez) + 0.051967 * m.sin(2 * rthez) - \
          0.154527 * m.cos(3 * rthez) + 0.084798 * m.sin(3 * rthez)
      rsdec = m.radians(sdec)
      # Time correction for solar hour angle, and solar hour angle
      tc = 0.004297 + 0.107029 * m.cos(rthez) - 1.837877 * m.sin(rthez) - \
          0.837378 * m.cos(2 * rthez) - 2.342824 * m.sin(2 * rthez)
      xha = (hr - 12) * 15 + xlon + tc
      xha[xha > 180] -= 360
      xha[xha < -180] += 360
      rlat = np.deg2rad(ylat)
      # rlon = np.deg2rad(xlon)
      rha = np.radians(xha)
      # Sun zenith
      costmp = np.sin(rlat) * np.sin(rsdec) + np.cos(rlat) * np.cos(rsdec) * np.cos(rha)
      # eps = abs(costmp)
      costmp[costmp > 1] = 1
      costmp[costmp < -1] = -1
      rsunz = np.arccos(costmp)
      sunz = np.rad2deg(rsunz)
      # Sun azimuth
      sna = m.cos(rsdec) * np.sin(rha)
      csa = np.sin(rlat) * np.cos(rha) * m.cos(rsdec) - m.sin(rsdec) * np.cos(rlat)
      rsuna = np.arctan2(sna, csa) + m.pi
      suna = np.rad2deg(rsuna)
      self.geomDict['sunposzen'] = sunz
      self.geomDict['sunposaz'] = suna

    @TimeIt
    def __GetFirstLastFrameLatsLons(self):
      """
      Returns
      -------
      latsFF,latsLF: first and last frame latitude
      lonsFF,lonsLF: first and last frame longitude
      """
      t_start = time.time()
      ALPHA_L = m.radians(-3.461515)
      ALPHA_R = m.radians(3.461514)
      # sensor location
      mean_sens_alt, mean_sens_lat, gammaDeg = self.__MSALGam()
      # quaternions to euler
      firstFramePhi, firstFrameTheta, firstFramePsi = self.__Quats2Euler(self.L0.header['FFquat'])
      lastFramePhi, lastFrameTheta, lastFramePsi = self.__Quats2Euler(self.L0.header['LFquat'])
      if (abs(firstFramePsi) <= 20) and (abs(lastFramePsi) <= 20):
          hofq = '+XVV'
      elif (abs(firstFramePsi) >= 160) and (abs(lastFramePsi) >= 160):
          hofq = '-XVV'

      if self.params['issOrientation'] != hofq:
          msg = 'ISS orientation in input file does not agree with the value \
                      derived from the quaternions in the L0 file...'
          self.logger.warning(msg)
      # build (3x3) rotation matrices for first & last frames
      firstFrameRotMat = self.__RotMatrix(firstFramePhi, firstFrameTheta, firstFramePsi)
      lastFrameRotMat = self.__RotMatrix(lastFramePhi, lastFrameTheta, lastFramePsi)
      centerAngleDeg = self.L0.data['ScenePointingAngle'] - 60.0  # signed angle from nadir
      if self.params['issOrientation'] == '-XVV':
          centerAngleDeg *= -1.0
      centerAngleRad = m.radians(centerAngleDeg)
      # north side is always to left of path in +XVV
      northAngleRad = centerAngleRad + ALPHA_L
      southAngleRad = centerAngleRad + ALPHA_R
      vecNorthFF = self.__GetVecs(firstFrameRotMat, northAngleRad)
      vecCenterFF = self.__GetVecs(firstFrameRotMat, centerAngleRad)
      vecSouthFF = self.__GetVecs(firstFrameRotMat, southAngleRad)
      vecNorthLF = self.__GetVecs(lastFrameRotMat, northAngleRad)
      vecCenterLF = self.__GetVecs(lastFrameRotMat, centerAngleRad)
      vecSouthLF = self.__GetVecs(lastFrameRotMat, southAngleRad)
      zenithNorthFF, bearingNorthFF = self.__Vec2ZenithBearing(vecNorthFF)
      zenithCenterFF, bearingCenterFF = self.__Vec2ZenithBearing(vecCenterFF)
      zenithSouthFF, bearingSouthFF = self.__Vec2ZenithBearing(vecSouthFF)
      zenithNorthLF, bearingNorthLF = self.__Vec2ZenithBearing(vecNorthLF)
      zenithCenterLF, bearingCenterLF = self.__Vec2ZenithBearing(vecCenterLF)
      zenithSouthLF, bearingSouthLF = self.__Vec2ZenithBearing(vecSouthLF)

      if self.params['issOrientation'] == '+XVV':
          distNorthFF = self.__GetDist(zenithNorthFF, self.L0.data['FFLatLonH'][2])
          distNorthLF = self.__GetDist(zenithNorthLF, self.L0.data['FFLatLonH'][2])
          distSouthFF = self.__GetDist(zenithSouthFF, self.L0.data['FFLatLonH'][2])
          distSouthLF = self.__GetDist(zenithSouthLF, self.L0.data['FFLatLonH'][2])
      else:
          distNorthFF = self.__GetDist(zenithSouthFF, self.L0.data['FFLatLonH'][2])
          distNorthLF = self.__GetDist(zenithSouthLF, self.L0.data['FFLatLonH'][2])
          distSouthFF = self.__GetDist(zenithSouthFF, self.L0.data['FFLatLonH'][2])
          distSouthLF = self.__GetDist(zenithSouthLF, self.L0.data['FFLatLonH'][2])
          bearingNorthFF, bearingSouthFF = bearingSouthFF, bearingNorthFF
          bearingNorthLF, bearingSouthLF = bearingSouthLF, bearingNorthLF
      distCenterFF = self.__GetDist(zenithCenterFF, self.L0.data['FFLatLonH'][2])
      distCenterLF = self.__GetDist(zenithCenterLF, self.L0.data['FFLatLonH'][2])
      bearingNorthFF += 90 - gammaDeg
      bearingSouthFF += 90 - gammaDeg
      bearingCenterFF += 90 - gammaDeg
      bearingNorthLF += 90 - gammaDeg
      bearingSouthLF += 90 - gammaDeg
      bearingCenterLF += 90 - gammaDeg
      latsFF, lonsFF, _ = self.__RngBear2LatLonEllipsoid(np.repeat(self.L0.data['FFLatLonH'][0],
                                                                   3),
                                                         np.repeat(self.L0.data['FFLatLonH'][1],
                                                                   3),
                                                         np.array([distNorthFF, distCenterFF,
                                                                  distSouthFF]),
                                                         np.array([bearingNorthFF,
                                                                   bearingCenterFF,
                                                                   bearingSouthFF]))
      lonsFF[lonsFF < -180] += 360
      lonsFF[lonsFF >= 180] -= 360
      latsLF, lonsLF, _ = self.__RngBear2LatLonEllipsoid(np.repeat(self.L0.data['LFLatLonH'][0],
                                                                   3),
                                                         np.repeat(self.L0.data['LFLatLonH'][1],
                                                                   3),
                                                         np.array([distNorthLF, distCenterLF,
                                                                   distSouthLF]),
                                                         np.array([bearingNorthLF,
                                                                   bearingCenterLF,
                                                                   bearingSouthLF]))
      lonsLF[lonsLF < -180] += 360
      lonsLF[lonsLF >= 180] -= 360
      angleDict = {'latsFF': latsFF, 'latsLF': latsLF,
                   'lonsFF': lonsFF, 'lonsLF': lonsLF,
                   'centerAngleDeg': centerAngleDeg, 'gammaDeg': gammaDeg}
      return angleDict

    @TimeIt
    def __GetLonsLats(self, latsFF, latsLF, lonsFF, lonsLF):
      nl_predark = self.L0.data['n_predark']
      nl_postdark = self.L0.data['n_postdark']
      nll = self.L0.data['nl']
      ns = self.L0.data['ns']
      nl = nll - nl_predark - nl_postdark
      n21 = nl_predark + 3
      n22 = nl + nl_predark - 1
      n2t = n22 - n21 + 1
      # Deal with lats/lons first
      imCenLat = (latsFF[1] + latsLF[1]) * 0.5
      if (np.abs(lonsLF[1] - lonsFF[1]) < 100):
          imCenLon = (lonsLF[1] - lonsFF[1]) * 0.5
      else:
          imCenLon = ((np.max([lonsFF[1], lonsLF[1]]) - 360) +
                      np.min(lonsLF[1], lonsFF[1]) * 0.5)
          if imCenLon < -180:
              imCenLon += 360
      self.geomDict['imCenLon'] = imCenLon
      self.geomDict['imCenLat'] = imCenLat
      fullLats = np.vstack((latsFF, latsLF))
      fullLats = self.__LinearCongrid(fullLats, (ns, n2t+3))

      ffeasthem, = np.where(((lonsFF > 0) & (lonsFF <= 180)))
      ffwesthem, = np.where(((lonsFF <= 0) | (lonsFF > 180)))
      lfwesthem, = np.where(((lonsLF > -180) & (lonsLF <= 0)))
      nffe, nffw, nlfw = ffeasthem.size, ffwesthem.size, lfwesthem.size
      cross_date_line = False
      if (nffe > 0) & (nlfw > 0):
          lonsLF[lfwesthem] += 360
          cross_date_line = True
          if nffw > 0:
              lonsFF[ffwesthem] += 360
      fullLons = np.vstack((lonsFF, lonsLF))
      fullLons = self.__LinearCongrid(fullLons, (ns, n2t+3))
      if cross_date_line:
          fullLons[fullLons > 180] -= 360
      if ((self.params['issOrientation'] == '+XVV') & (self.params['flip_left_right'])) | \
         ((self.params['issOrientation'] == '-XVV') & (not self.params['flip_left_right'])):
              fullLats = np.transpose(np.rot90(fullLats, 1))
              fullLons = np.transpose(np.rot90(fullLons, 1))
      self.geomDict['fullLats'] = fullLats.T
      self.geomDict['fullLons'] = fullLons.T

    @TimeIt
    def __GetViewAngles(self, centAngDeg, gammaDeg):
      '''
      Fill geomDict viewA and viewZ and viewAZ
      '''
      if centAngDeg < 0:
          viewaz = 180 - gammaDeg
      elif centAngDeg > 0:
          viewaz = 360 - gammaDeg
      else:
          viewaz = 0
      viewaz %= 360
      self.geomDict['viewaz'] = viewaz
      pixels = (np.arange(self.L0.data['ns']) - 256.5) / 255.5
      viewl0 = 3.487 * pixels - 0.0035 * pixels ** 3
      spa = self.L0.data['ScenePointingAngle']
      vl = (spa + viewl0 - 60)
      if self.params['issOrientation'] != '+XVV':
          vl = -vl
      viewl = np.abs(vl)
      azNlocs = np.where(vl <= 0, True, False)
      azSlocs = ~azNlocs
      viewa_line = np.zeros((self.L0.data['ns'], 1))
      viewa_line[azNlocs] = 180 - gammaDeg
      viewa_line[azSlocs] = 360 - gammaDeg
      viewa_line %= 360
      viewz = np.tile(np.reshape(viewl, (viewl.shape[0], 1)), (1, self.ancDict['n2t'] + 3))
      viewa = np.tile(viewa_line, (1, self.ancDict['n2t'] + 3))
      if self.params['flip_left_right']:
          viewa = np.transpose(np.rot90(viewa, 1))
          viewz = np.transpose(np.rot90(viewz, 1))
      self.geomDict['viewa'] = viewa.T
      self.geomDict['viewz'] = viewz.T

    @TimeIt
    def __FillGeomDict(self):
      '''
      Fills geomDict, calls four auxiliaries functions;
      __GetFirstLastFrameLatsLon()
      __GetLonsLats() -> fills imcenlon/lat and fullLons/lats
      __GetSunAngles() -> fills sola/solz
      __GetViewAngles() -> fills viewaz/viewa/viewz
      '''
      angDict = self.__GetFirstLastFrameLatsLons()
      self.__GetLonsLats(angDict['latsFF'], angDict['latsLF'],
                         angDict['lonsFF'], angDict['lonsLF'])
      self.__GetSunAngles()
      self.__GetViewAngles(angDict['centerAngleDeg'], angDict['gammaDeg'])

    @TimeIt
    def _StoreBadData(self):
      """Computes and stores bad, dropped and saturated pixels"""
      n21 = self.ancDict['n21']
      n22 = self.ancDict['n22']
      n2t = self.ancDict['n2t']
      ns = self.L0.data['ns']
      nb = self.L0.data['nb']
      self.badDataDict['bad_pixs'] = np.zeros((n2t+3, nb, ns), dtype=bool)
      self.badDataDict['sat_pixs'] = np.zeros((n2t+3, nb, ns), dtype=bool)
      intData = self.L0.data['pic0'][n21:n22+1]  # ref to subset of theData
      sat_pix = np.where(intData == 16383, True, False)
      bad_pix = np.where(intData == -1, True, False)
      if bad_pix.any():
          self.badDataDict['bad_pixs'][3:] = bad_pix
      if sat_pix.any():
          self.badDataDict['sat_pixs'][3:] = sat_pix

    @TimeIt
    def FillWriteFlags(self, ncgrp):
      '''
      Fills bad data flags array. Flags array is then saved to file.
      Flagbits used are 0, 1, 2, 4, 5, 6; also 7 if doCal option is on.
      Bit 3, NAVWARN, is set by default

      '''
      flags = np.ones((self.ancDict['n2t']+3, self.L0.data['ns']), dtype='int8') * 4
      bad_sunzen = np.where(self.geomDict['sunposzen'] > 88, True, False)
      if bad_sunzen.any():
          self.ancDict['rhot_nir'][bad_sunzen] *= m.pi
          self.ancDict['rhot_red'][bad_sunzen] *= m.pi
      if not bad_sunzen.all():
          goodsz = np.logical_not(bad_sunzen)
          self.ancDict['rhot_nir'][goodsz] *= (m.pi *
                                               np.cos(np.deg2rad
                                                      (self.geomDict['sunposzen'][goodsz])))
          self.ancDict['rhot_red'][goodsz] *= (m.pi *
                                               np.cos(np.deg2rad
                                                      (self.geomDict['sunposzen'][goodsz])))
      ptr = np.where(((self.geomDict['fullLats'] > 90) |
                     (self.geomDict['fullLats'] < -90) | (self.geomDict['fullLons'] > 180) |
                     (self.geomDict['fullLons'] < -180)), True, False)
      flags[ptr] = flags[ptr] | 2  # bit 2 is navfail
      ptr = np.where(self.geomDict['viewz'] > 60, True, False)
      flags[ptr] = flags[ptr] | 8  # bit 4 (?) is hisatzen
      ptr = np.where(self.geomDict['sunposzen'] > 75, True, False)
      flags[ptr] = flags[ptr] | 16  # bit 5 (??) is hisolzen
      ptr = np.where(self.badDataDict['sat_pixs'] > 0, True, False)
      flags[ptr] = flags[ptr] | 32
      ptr = np.where(self.badDataDict['bad_pixs'] > 0, True, False)
      flags[ptr] = flags[ptr] | 64
      flags.tofile(self.outFilNamesDict['flagEnvi'])
      if self.inpDict['doCal'] & self.L0.data['nb'] == 128:
          ptr = np.where((self.ancDict['rhot_nir'] > 0.02), True, False)
          flags[ptr] = flags[ptr] | 1  # bit 1 is land
          pos_red = np.where(self.ancDict['rhot_red'] > 0)
          ratio = np.zeros((self.L0.data['ns'], self.ancDict['n2t']+3))
          ratio[pos_red] = self.ancDict['rhot_nir'][pos_red] / self.ancDict['rhot_red'][pos_red]
          ptr = np.where((((self.ancDict['rhot_red'] > 0.5) & (self.ancDict['rhot_nir'] > 0.5))
                          | ((ratio > 0.8) & (ratio < 1.1))), True, False)
          flags[ptr] = flags[ptr] | 128
      ncgrp['scan_quality_flags'][:] = flags.T

    @TimeIt
    def WriteRadFile(self, prodGrp, periodGrp):
      '''
      Writes _rad.bil (and _rad.hdr files?).

      Steps: assess whether to export as bil or bip and transpose accordingly,
             extract data
             write data
      Note: The data array, initially, has dims nl x nb x ns, where
      nl = number of lines, nb = number of bands, ns = number of samples.
      '''
      nlsub = self.ancDict['n2t'] + 3
      nb = self.L0.data['nb']
      ns = self.L0.data['ns']
      outArray = np.zeros((nlsub, ns, nb))
      outArray[3:] = self.floatData.transpose(0, 2, 1)
      # outArray dim is (nl x nb x ns)
      # we want dims nl x ns x nb
      lt = prodGrp['Lt']
      lt[:] = outArray
      # fwhm = prodGrp['fwhm']
      # fwhm[:] = self.params["fwhm"]
      lt.wavelengths = self.params["wave_um"].astype('f4') * 1000  # switch to nm
      lt.fwhm = self.params['fwhm'].astype('f4')
      # wavelengths = prodGrp['wavelengths']
      # wavelengths[:] = self.params["wave_um"] * 1000  # switch to nm
      periodGrp.Beginning_Date = self.params['begin_date']
      periodGrp.Ending_Date = self.params['end_date']
      periodGrp.Beginning_Time = self.params['begin_time']
      periodGrp.Ending_Time = self.params['end_time']
