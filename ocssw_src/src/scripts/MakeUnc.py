#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 15:09:44 2015
script to run analytic uncertainty analysis of rrs data
@author: ekarakoy
"""
import numpy as np
import netCDF4 as nc
import glob
import sys
import os
import re
import argparse
import logging
import multiprocessing as mp
import pickle
import shutil
from datetime import datetime as DT


class MakeUnc(object):
    """
    Class to get Rrs uncertainties for a given L2 granule. Includes methods to:
        * get baseline from L2 granule
        * calculate uncertainties (as rmse) from corresponding perturbed L2
            files
        * save uncertainty variables in original unperturbed granule
    Inputs:
        *args:
            1- baselineFile
            2- noisyDir -- directory where noisy files are located


    """
    def __init__(self, pArgs, parent_logger_name):
        self.logger = logging.getLogger('%s.MakeUnc' % parent_logger_name)
        self.silFile = pArgs.ifile
        self.noisyDir = pArgs.npath
        # Process options
        self.doChla = pArgs.dochl
        self.doNflh = pArgs.doflh
        self.pSafe = pArgs.psafe
        self.doSaniCheck = pArgs.sanity
        self.rrsSilDict = dict.fromkeys(self.bands)
        self.attrRrsUncDict = dict.fromkeys(self.bands)
        self.dimsDict = dict.fromkeys(['rrs_unc_%s' % band for band in self.bands])
        self.dTypeDict = dict.fromkeys(['rrs_unc_%s' % band for band in self.bands])
        self.rrsUncArrDict = dict.fromkeys(self.bands)
        if self.pSafe:
            self._PlaySafe()
        if self.doSaniCheck:
            self.dimsDict.update(dict.fromkeys(['lt_unc_%s' % band for band in self.bands]))
            self.dTypeDict.update(dict.fromkeys(['lt_unc_%s' % band for band in self.bands]))
            self.ltUncArrDict = dict.fromkeys(self.bands)
            self.ltSilDict = dict.fromkeys(self.bands)
            self.attrLtUncDict = dict.fromkeys(self.bands)
        otherProdkeys = []
        if self.doChla:
            self.logger.info('chl_a_unc will be included')
            otherProdkeys.append('chlor_a')
            otherProdkeys.append('chlor_a_unc')
        if self.doNflh:
            self.logger.info('nflh_unc will be included')
            otherProdkeys.append('nflh')
            otherProdkeys.append('nflh_unc')
        if otherProdkeys:
            self.otherProdsDict = dict.fromkeys(otherProdkeys)
            attrUncKeys = [x for x in otherProdkeys if re.search('_unc', x)]
            self.attrOtherProdUncDict = dict.fromkeys(attrUncKeys)
            self.dTypeDict.update(dict.fromkeys(attrUncKeys))
            self.dimsDict.update(dict.fromkeys(attrUncKeys))

    def _PlaySafe(self):
        '''
        Method to copy backup of unprocessed silent L2
        This so as not to redo entire processing if a problem arises.
        If a copy already exists, it is assumed this is not the first
        processing attempt and the silent L2 is now tainted. It is removed and
        a clean copy is generated from the backup.
        '''
        orig = self.silFile
        cpy = self.silFile + '.cpy'
        if os.path.exists(cpy):  # not the first time - something's wrong
            os.remove(orig)  # remove "original - tainted file"
            self.logger.info('%s already exists. Removing original %s' % (cpy,
                                                                          orig)
                             )
            shutil.copy2(cpy, orig)
            self.logger.info('Copying silent from %s' % cpy)
        else:
            shutil.copy2(orig, cpy)
            self.logger.info('No copy for %s. Generating copy' % self.silFile)

    def WriteToSilent(self):
        # first create NC variables if necessary
        # save unc in corresponding NC variables.
        with nc.Dataset(self.silFile, 'r+') as dsSil:
            geoGr = dsSil.groups['geophysical_data']
            geoVar = geoGr.variables
            for band in self.bands:
                rrsUncProdName = 'Rrs_unc_%s' % band
                if rrsUncProdName not in geoVar:
                    try:
                        varRrsUnc = geoGr.createVariable(rrsUncProdName,
                                                         self.dTypeDict['rrs_unc_%s' % band],
                                                         self.dimsDict['rrs_unc_%s' % band])
                        varRrsUnc.setncatts(self.attrRrsUncDict[band])
                    except KeyError:  # skip nir bands
                        pass
                else:
                    varRrsUnc = geoVar[rrsUncProdName]
                varRrsUnc[:] = self.rrsUncArrDict[band]
                if self.doSaniCheck:
                    ltUncProdName = 'Lt_unc_%s' % band
                    if ltUncProdName not in geoVar:
                        varLtUnc = geoGr.createVariable(ltUncProdName,
                                                        self.dTypeDict['lt_unc_%s' % band],
                                                        self.dimsDict['lt_unc_%s' % band])
                        varLtUnc.setncatts(self.attrLtUncDict[band])
                    else:
                        varLtUnc = geoVar[ltUncProdName]
                    varLtUnc[:] = self.ltUncArrDict[band]
            if self.doChla:
                if 'chlor_a_unc' not in geoVar:
                    varChloraUnc = geoGr.createVariable('chlor_a_unc',
                                                        self.dTypeDict['chlor_a_unc'],
                                                        self.dimsDict['chlor_a_unc'])
                    varChloraUnc.setncatts(self.attrOtherProdUncDict['chlor_a_unc'])
                else:
                    varChloraUnc = geoVar['chlor_a_unc']
                varChloraUnc[:] = self.otherProdsDict['chlor_a_unc']

            if self.doNflh:
                if 'nflh_unc' not in geoVar:
                    self.logger.info('nflh_unc there; creating variable...')
                    varNflhUnc = geoGr.createVariable('nflh_unc',
                                                      self.otherProdsDict['nflh_unc'].dtype,
                                                      self.dimsDict['nflh_unc'])
                    varNflhUnc.setncatts(self.attrOtherProdUncDict['nflh_unc'])
                else:
                    self.logger.info('nflh_unc available, using existing variable...')
                    varNflhUnc = geoVar['nflh_unc']
                varNflhUnc[:] = self.otherProdsDict['nflh_unc']
        # self.logger.info("%s Processing Complete" % baseLineFname)
        return None

    def BuildUncs(self, noisySfx):
        """"
        Calculates rrs uncertainty as st.dev of rrs. Note that to save time
            I use unperturbed rrs as the rrs baseline for the simulation
        """
        fBaseName = self.silFile.split('/')[-1].split('.')[0].split('_')[0]
        matchFilPatt = os.path.join(self.noisyDir, '%s*%s*' % (fBaseName, noisySfx))
        self.logger.info("searching for %s..." % matchFilPatt)
        flis = glob.glob(matchFilPatt)
        lflis = len(flis)
        if lflis == 0:
            self.logger.error('No files to process with pattern %s' % matchFilPatt)
            sys.exit(1)
        else:
            self.logger.info("%d files to be processed" % lflis)
        rrsAggDataDict = {band: np.zeros_like(self.rrsSilDict[band])
                          for band in self.rrsSilDict.keys()}
        if self.doSaniCheck:
            ltAggDataDict = {band: np.zeros_like(self.ltSilDict[band])
                             for band in self.ltSilDict.keys()}
        if self.doChla:
            chlAggDataArr = np.zeros_like(self.otherProdsDict['chlor_a'])
        if self.doNflh:
            nflhAggDataArr = np.zeros_like(self.otherProdsDict['nflh'])
        # process noisy data
        for fcount, fname in enumerate(flis):
            prcDone = 100 * fcount/ (lflis - 1)
            self.logger.info("Loading and reading %s -- %.1f%%" %
                             (fname, prcDone))
            with nc.Dataset(fname) as nds:
                nGeoGr = nds.groups['geophysical_data']
                nGeoVar = nGeoGr.variables
                for band in self.bands:
                    try:
                        noisyRrs = nGeoVar['Rrs_%s' % band][:]
                        rrsAggDataDict[band] += (noisyRrs -
                                                 self.rrsSilDict[band]) ** 2
                    except KeyError as e:
                        self.logger.error(e)
                    if self.doSaniCheck:
                        noisyLt = nGeoVar['Lt_%s' % band][:]
                        ltAggDataDict[band] += (noisyLt -
                                                self.ltSilDict[band]) ** 2
                if self.doChla:
                    noisyChl = nGeoVar['chlor_a'][:]
                    chlAggDataArr += (noisyChl - self.otherProdsDict['chlor_a']) ** 2
                if self.doNflh:
                    noisyNflh = nGeoVar['nflh'][:]
                    nflhAggDataArr += (noisyNflh - self.otherProdsDict['nflh']) ** 2

        for band in self.bands:
            self.logger.debug("computing deviation for band %s" % band)
            try:
                self.rrsUncArrDict[band] = np.ma.sqrt(rrsAggDataDict[band]/ lflis)
            except KeyError:  # band not here just skipping
                pass
            if self.doSaniCheck:
                self.ltUncArrDict[band] = np.sqrt(ltAggDataDict[band]/ lflis)
                self.logger.debug('running sanity check for band %s' % band)
            if self.doChla:
                self.otherProdsDict['chlor_a_unc'] = np.ma.sqrt(chlAggDataArr
                                                               / lflis)
                self.logger.debug('computing deviation for chlor a')
            if self.doNflh:
                self.otherProdsDict['nflh_unc'] = np.ma.sqrt(nflhAggDataArr
                                                            / lflis)
                self.logger.debug('computing deviation for nflh')
        self.logger.info("\nProcessed %d files " % lflis)
        return None

    def ReadFromSilent(self):
        '''Reads Baseline file
            Flags: l2bin default flags, namely ATMFAIL(1), LAND(2), HIGLINT(8),
            HILT(16), HISATZEN(32), STRAYLIGHT(256), CLDICE(512),
            COCCOLITH(1024), HISOLZEN(4096), LOWLW(16384), CHLFAIL(32768),
            NAVWARN(65536), MAXAERITER(524288), CHLWARN(2097152),
            ATMWARN(4194304), NAVFAIL(33554432), FILTER(67108864)
            flagBits = 1 + 2 + 8 + 16 + 32 +  256 + 512 + 1024 + 4096 + 16384 +
            32768 +  65536 + 524288 + 2097152 + 4194304 + 33554432 + 67108864
            l2flags = geoVar['l2_flags'][:]
            flagMaskArr = (l2flags & flagBits > 0)
        '''
        self.logger.debug('attemping to open silent file %s' % self.silFile)
        with nc.Dataset(self.silFile, 'r') as dsSil:
            geoGr = dsSil.groups['geophysical_data']
            geoVar = geoGr.variables
            for band in self.bands:
                try:
                    rrs = geoVar['Rrs_%s' % band]
                    self.rrsSilDict[band] = rrs[:]
                    self.attrRrsUncDict[band] = {'long_name': 'Uncertainty in %s' % rrs.long_name,
                                                 '_FillValue': rrs._FillValue, 'units': rrs.units,
                                                 'scale_factor': rrs.scale_factor,
                                                 'add_offset': rrs.add_offset}
                    self.dimsDict['rrs_unc_%s' % band] = rrs.dimensions
                    self.dTypeDict['rrs_unc_%s' % band] = rrs.dtype
                except KeyError as e:
                    self.logger.info("Rrs_%s not available, skipping" % band)
                    self.rrsSilDict.pop(band)
                    self.attrRrsUncDict.pop(band)
                    self.dimsDict.pop('rrs_unc_%s' % band)
                    self.dTypeDict.pop('rrs_unc_%s' % band)
                if self.doSaniCheck:
                    self.logger.debug('setting up to run sanity check for band %s' % band)
                    lt = geoVar['Lt_%s' % band]
                    self.ltSilDict[band] = lt[:]
                    self.attrLtUncDict[band] = {'long_name': 'Uncertainty in %s' % lt.long_name,
                                                '_FillValue': lt._FillValue, 'units': lt.units}
                    self.dimsDict['lt_unc_%s' % band] = lt.dimensions
                    self.dTypeDict['lt_unc_%s' % band] = lt.dtype
            if self.doChla:
                self.logger.debug('setting up to compute chla uncertainty')
                chla = geoVar['chlor_a']
                self.otherProdsDict['chlor_a'] = chla[:]
                self.attrOtherProdUncDict['chlor_a_unc'] = {'long_name':
                                                            'Uncertainty in %s' % chla.long_name,
                                                            '_FillValue': chla._FillValue,
                                                            'units': chla.units,
                                                            'valid_min': 0}
                self.dTypeDict['chlor_a_unc'] = chla.dtype
                self.dimsDict['chlor_a_unc'] = chla.dimensions
            if self.doNflh:
                self.logger.debug('setting up to compute nflh uncertainty')
                nflh = geoVar['nflh']
                self.otherProdsDict['nflh'] = nflh[:]
                self.attrOtherProdUncDict['nflh_unc'] = {'long_name':
                                                         'Uncertainty in %s' % nflh.long_name,
                                                         '_FillValue': nflh._FillValue,
                                                         'units': nflh.units,
                                                         'scale_factor': nflh.scale_factor,
                                                         'add_offset': nflh.add_offset}
                self.dimsDict['nflh_unc'] = nflh.dimensions
                self.dTypeDict['nflh_unc'] = nflh.dtype
        return None


class MakeSwfUnc(MakeUnc):
    """Uncertainty subclass for SeaWiFS"""
    def __init__(self, *args, **kwargs):
        self.sensor = 'SeaWiFS'
        self.bands = kwargs.pop("bands",
                                ['412', '443', '490', '510', '555', '670',
                                 '765', '865'])
        self.colDict = {'412': '#001166', '443': '#004488', '490': '#116688',
                        '510': '#228844', '555': '#667722', '670': '#aa2211',
                        '765': '#770500', '865': '#440000'}
        super(MakeSwfUnc, self).__init__(*args)
        return None


class MakeHMA(MakeUnc):
    """Uncertainty engine for HMODISA"""
    def __init__(self, *args, **kwargs):
        self.sensor = 'HMODISA'
        self.bands = kwargs.pop("bands",
                                ['412', '443', '488', '531', '547', '555',
                                 '645', '667', '678', '748', '859', '869',
                                 '1240', '1640', '2130'])
        super(MakeHMA, self).__init__(*args)
        self.colDict = {'412': '#001166', '443': '#004488', '488': '#1166FF',
                        '531': '#337722', '547': '#557733', '555': '#669922',
                        '645': '#883311', '667': '#aa2211', '678': '#dd3300'}
        return None


def PathsGen(matchPattern, l2MainPath):
    # create generator of l2 directory paths
    l2PathsGen = glob.iglob(matchPattern)
    spatt = re.compile('(S[0-9]+)')
    for l2path in l2PathsGen:
        if os.path.isdir(l2path):
            basename = spatt.findall(l2path)[0]
            l2Pa = os.path.join(l2MainPath, basename)
            silFiPa = os.path.join(l2Pa, basename) + '_silent.L2'
            noiDiPa = os.path.join(l2Pa, 'Noisy/')
        else:
            continue
        yield [silFiPa, noiDiPa]


class CBatchManager():
    '''
    Class to manage complete uncertainty generation; from processing of L1As to
    creation of uncertainty from noisy L2 files, to the final packing of new
    uncertainty products into the baseline L2.
    '''

    def __init__(self, pArgs):
        '''
        Takes a directory containing L2 directories
        '''
        self.pArgs = pArgs
        self.l2MainPath = pArgs.ipath
        if self.pArgs.sensor == 'SeaWiFS':
            self.matchPattern = os.path.join(self.l2MainPath, 'S*/')
        return None

    def _BatchRun(self, sArgs):
        ifile, npath = sArgs
        uncObj = MakeSwfUnc(ifile, npath)
        uncObj.ReadFromSilent()
        uncObj.BuildUncs(self.pArgs.nsfx)
        uncObj.WriteToSilent()
        return uncObj.silFile

    def ProcessL2s(self):
        paramGen = (params for params in PathsGen(self.matchPattern,
                                                  self.l2MainPath))
        with mp.Pool() as pool:
            results = pool.map(self._BatchRun, paramGen)
        return results  # temporary: should be replaced by log entry


def ParseCommandLine(args):
    parser = argparse.ArgumentParser(prog="MakeUnc")
    parser.add_argument('-i', '--ifile', help='Initial L2 file path.',
                        type=str)
    parser.add_argument('-j', '--ipath',
                        help='Main L2 path for batch processing.', type=str)
    parser.add_argument('-n', '--npath', help='Path to noisy data directory.',
                        type=str)
    parser.add_argument('-s', '--nsfx',
                        help='Noisy file suffix for pattern matching.',
                        type=str, default='_noisy_')
    parser.add_argument('-c', '--dochl', help='Compute chloropyll uncertainty. Default is False',
                        action='store_true', default=False)
    parser.add_argument('-f', '--doflh',
                        help='Compute normalized fluorescence line height. Default is False',
                        action='store_true', default=False)
    parser.add_argument('-p', '--psafe', help='Back source file up. Default is False',
                        action='store_true', default=False)
    parser.add_argument('-a', '--sanity', help='Do sanity check. Default is False',
                        action='store_true', default=False)
    parser.add_argument('-b', '--batch', help='Batch processing option. Default is False',
                        action='store_true', default=False)
    parser.add_argument('-w', '--workers',
                        help='Number of concurrent processes. Default is 1',
                        type=int, default=1)
    parser.add_argument('-z', '--sensor',
                        help='Specify sensor data originates from. Default is SeaWiFS',
                        type=str, default='SeaWiFS')
    parser.add_argument('-e', '--debug', help='increase output verbosity',
                        action='store_true', default=False)
    parsedArgs = parser.parse_args(args)
    
    if ((not parsedArgs.ifile and not parsedArgs.batch) or 
       (parsedArgs.batch and not parsedArgs.ipath)):
           parser.print_help()
           sys.exit(1)
           
    return parsedArgs


def SetLogger(dbg=False):
    '''
    sets logger with more verbose format if dbg_lvl=True
    '''
    logger_name = 'MakeUNC_%s_T_%s' % (DT.date(DT.now()), DT.time(DT.now()))
    logfn = '%s.log' % logger_name
    logger = logging.getLogger(logger_name)
    if dbg:
        level = logging.DEBUG
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s -'
                                      + ' [%(module)s..%(funcName)s..%(lineno)d]'
                                      + ' - %(message)s')
    else:
        level = logging.INFO
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger.setLevel(level)
    fh = logging.FileHandler(logfn)
    fh.setLevel(level)
    fh.setFormatter(formatter)
    ch = logging.StreamHandler()
    ch.setLevel(logging.WARNING)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.addHandler(fh)
    logger.debug('logging')
    return logger


def Main(argv):

    pArgs = ParseCommandLine(argv)
    mainLogger = SetLogger(dbg=pArgs.debug)
    if pArgs.batch:
        # min. cmd line is ipath for main L2Path (all L2s should be in a
        # common directory. ) and -b
        mainLogger.info('Initializing batch processor')
        bRunner = CBatchManager(pArgs)
        res = bRunner.ProcessL2s()
        pickle.dump(res, open('L2BatchList.pkl', 'wb'))
    else:
        baseLineFile = pArgs.ifile
        noisyDataDir = pArgs.npath
        noisySfx = pArgs.nsfx
        baseLineFname = baseLineFile.split('/')[-1]
        if noisyDataDir[-1] != '/':
            noisyDataDir += '/'
        if baseLineFname[0] == 'S':
            mainLogger.info('processing SeaWiFS data')
            uncObj = MakeSwfUnc(pArgs, mainLogger.name)
        elif baseLineFname[0] == 'A':
            mainLogger.info('processing MODISA data')
            uncObj = MakeHMA(pArgs)
        uncObj.ReadFromSilent()
        uncObj.BuildUncs(noisySfx)
        uncObj.WriteToSilent()

if __name__ == "__main__":

    Main(sys.argv[1:])
