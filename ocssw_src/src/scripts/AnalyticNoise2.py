#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 09:55:41 2015
Data might need to be transformed before variance propagation
--> see this page for how to do it in python
    http://scipy.github.io/devdocs/generated/scipy.stats.boxcox.html
--> and this page for a discussion about it
    http://stats.stackexchange.com/questions/18844/when-and-why-to-take-the-log-of-a-distribution-of-numbers

@author: ekarakoy
"""
import numpy as np
import netCDF4 as nc
import argparse
#import matplotlib.pyplot as pl
#from matplotlib.backends.backend_pdf import PdfPages
import sys

class AnalyticNoise(object):
    """Class to manage uncertainty estimation using an analytical approach.
        The class takes an optional boolean argument, multiDim, that defaults
        to False, which signals whether the user wants to follow the contribu
        tion each band's Lt perturbation to each band's rrs noise estimation.

        The order of implementation would be to
         * instantiate a child object specific to a sensor
             - e.g. SwfNoise
         * load baseline data
             - object.GetBaseline()
         * process pertubed data files to generate variance
             - object.UpdVarFromPertDat --> Lt_band by Lt_band
             - note that obj.UpdateAllVars -->
         * plot to see all is well
         * save computed rrs uncertainties back into the baseline file

         I might write a full sequence method, but for now it's manual
    """


    def __init__(self,silFile,noisyDir,multiDimVar = False):

        # Paths
        self.silFile = silFile
        self.noisyDir = noisyDir
        # Options
        self.multiDimVar = multiDimVar
        # Ref. Dicts
        self.ltProdNameDict = dict.fromkeys(self.bands)
        self.rrsProdNameDict = dict.fromkeys(self.bands)
        self.rrsUncProdNameDict = dict.fromkeys(self.bands)
        # Input Data Dicts
        self.ltSilDict = dict.fromkeys(self.bands)
        self.rrsSilDict = dict.fromkeys(self.bands)
        self.attrDict = dict.fromkeys(self.bands)
        self.dimsDict = dict.fromkeys(self.bands)
        self.dTypeDict = dict.fromkeys(self.bands)
        # Output Data Dicts
        self.rrsUncDict = dict.fromkeys(self.bands)
        for band in self.bands:
            self.ltProdNameDict[band] = 'Lt_' + band
            self.rrsProdNameDict[band] = 'Rrs_' + band
            self.rrsUncProdNameDict[band] = 'Rrs_unc_' + band

    def ReadFromSilent(self):
        '''Reads Baseline file

            Flags: l2bin default flags, namely ATMFAIL(1), LAND(2), HIGLINT(8),
                HILT(16), HISATZEN(32), STRAYLIGHT(256), CLDICE(512), COCCOLITH(1024),
                HISOLZEN(4096), LOWLW(16384), CHLFAIL(32768), NAVWARN(65536),
                MAXAERITER(524288), CHLWARN(2097152), ATMWARN(4194304),
                NAVFAIL(33554432), FILTER(67108864)
        '''
        flagBits = 1 + 2 + 8 + 16 + 32 +  256 + 512 + 1024 + 4096 + 16384 + \
                32768 +  65536 + 524288 + 2097152 + 4194304 + 33554432 + 67108864
        with nc.Dataset(self.silFile,'r') as dsSil:
            geoGr = dsSil.groups['geophysical_data']
            geoVar = geoGr.variables
            l2flags = geoVar['l2_flags'][:]
            flagMaskArr = (l2flags & flagBits > 0)
            for band in self.bands:
                variable = geoVar[self.rrsProdNameDict[band]]
                self.rrsSilDict[band] = variable[:]
                self.rrsSilDict[band].mask += flagMaskArr # Apply l2bin flags
                self.attrDict[band] = {'long_name' : 'Uncertainty in ' +
                                                          variable.long_name,
                                             '_FillValue': variable._FillValue,
                                             'units': variable.units,
                                             'scale_factor':variable.scale_factor,
                                             'add_offset':variable.add_offset}
                self.dimsDict[band] = variable.dimensions
                self.dTypeDict[band] = variable.dtype
                self.ltSilDict[band] = geoVar[self.ltProdNameDict[band]][:]

    def BuildUncs(self,noisySfx,verbose=False):

        for band_rrs in self.bands:
            if self.multiDimVar:
                self.rrsUncDict[band_rrs]={}
                for bandLt in self.bands:
                    band_lt = 'Lt_' + bandLt
                    self.rrsUncDict[band_rrs][band_lt] = np.zeros_like(self.rrsSilDict[bandLt])
                self.rrsUncDict[band_rrs]['Aggregate'] = np.zeros_like(self.rrsSilDict[band_rrs])
            else:
                self.rrsUncDict[band_rrs] = np.zeros_like(self.rrsSilDict[band_rrs])

        fBaseName = self.silFile.split('/')[-1]
        for i,ltBand in enumerate(self.bands):
            noisyFileName = self.noisyDir + '/' + fBaseName + noisySfx + str(i).zfill(4) + '.L2'
            ltKey = 'Lt_' + ltBand
            print("searching for %s" % noisyFileName, flush=True)
            with nc.Dataset(noisyFileName) as nds:
                if verbose:
                    print("Loading and reading %s" % noisyFileName)
                ngv = nds.groups['geophysical_data'].variables
                if verbose:
                    print("Processing perturbation of Lt_%s" % ltBand)
                ltPert = ngv[self.ltProdNameDict[ltBand]][:]
                ltSigma = self._ApplyPolyFit(self.coeffsDict[ltBand],
                                             self.ltSilDict[ltBand],self.polyDeg)
                ltSigma = self._ApplyLims(ltSigma,self.sigLimDict[ltBand])
                dLt = ltPert - self.ltSilDict[ltBand]
                for band in self.bands:
                    rrsPert = ngv[self.rrsProdNameDict[band]][:]
                    rrsPert.mask += self.rrsSilDict[band].mask
                    sigTemp = ltSigma * (rrsPert - self.rrsSilDict[band]) / dLt
                    varTemp = sigTemp **2
                    if verbose:
                        print("----> estimating %s" % self.rrsUncProdNameDict[band])
                    if self.multiDimVar:
                        #update current bands's Lt perturbation contribution to rrs noise
                        self.rrsUncDict[band][ltKey] += varTemp # CAN I JUST TO THAT ON THE COMPRESSED DATA?
                        # aggregate current band's contribution to rrs noise
                        self.rrsUncDict[band]['Aggregate'] += varTemp
                    else:
                        # only aggregate current band's contribution to rrs noise
                        self.rrsUncDict[band] += varTemp

    def WriteToSilent(self):
        # first create NC variables if necessary
        # save unc in corresponding NC variables.
        with nc.Dataset(self.silFile,'r+') as dsSil:
            geoGr = dsSil.groups['geophysical_data']
            geoVar = geoGr.variables
            for band in self.bands:
                uncProdName  = self.rrsUncProdNameDict[band]
                if uncProdName not in geoVar:
                    dType = self.dTypeDict[band]
                    dIms = self.dimsDict[band]
                    variable = geoGr.createVariable(uncProdName,dType,dIms)
                    variable.setncatts(self.attrDict[band])
                else:
                    variable=geoVar[uncProdName]
                if self.multiDimVar:
                    variable[:] = self.rrsUncDict[band]['Aggregate']
                else:
                    variable[:] = self.rrsUncDict[band]

    @staticmethod
    def _ApplyLims(sigs,lims):
        """This function applies thresholds to sigmas from pre-computed
        upper and lower bounds

        """
        finiteSigsIdx = np.isfinite(sigs)
        sigsSub = sigs[finiteSigsIdx]
        sigsSub[np.where(sigsSub<lims[0])] = lims[0]
        sigsSub[np.where(sigsSub>lims[1])] = lims[1]
        sigs[finiteSigsIdx] = sigsSub
        return sigs

    @staticmethod
    def _ApplyPolyFit(cf,lt,deg):
        """This function applies a polynomial fit to the data in lt.
            cf: coefficients of the polynomial
            deg: degree of the polynomial
        """
        y = 0
        for i in range(deg):
            y += cf[i] * (lt** (deg - i))
        y += cf[deg]
        return 1.0/ y


class SeaWiFSNoise(AnalyticNoise):
    """Subclass to AnalyticNoise, specific to SeaWiFS"""

    def __init__(self,*args,**nargs):
        self.sensor='SeaWiFS'
        self.polyDeg = 4
        self.bands = ['412','443','490','510','555','670','765','865']
        self.coeffsDict={'412':np.array([-8.28726301e-11,3.85425664e-07,-9.10776926e-04,
                                     1.65881862e+00,4.54351582e-01]),
                    '443':np.array([-1.21871258e-10,5.21579320e-07,-1.14574109e-03,
                                    1.96509056e+00,4.18921861e-01]),
                    '490':np.array([-2.99068165e-10,1.05225457e-06,-1.90591166e-03,
                                    2.66343986e+00,6.67187489e-01]),
                    '510':np.array([-5.68939986e-10,1.67950509e-06,-2.56915149e-03,
                                    3.05832773e+00,9.34468454e-01]),
                    '555':np.array([-1.31635902e-09,3.09617393e-06,-3.73473556e-03,
                                    3.52394751e+00,3.54105899e-01]),
                    '670':np.array([-8.65458303e-09,1.18857306e-05,-8.37771886e-03,
                                    4.64496430e+00,4.14633422e-02]),
                    '765':np.array([-4.96827099e-08,4.50239057e-05,-2.10425126e-02,
                                    7.75862055e+00,5.18893137e-02]),
                    '865':np.array([-1.30487418e-07,9.35407901e-05,-3.40988182e-02,
                                    9.43414239e+00,7.84956550e-01])
                    }
        self.sigLimDict={'412':[0.0008826466289493833,0.058990630111850566],
                    '443':[0.00078708204284562292,0.050110810767361902],
                    '490':[0.00073633772526737848,0.036883976493020949],
                    '510':[0.00074975219339103519,0.031987200608546547],
                    '555':[0.00080870577569697015,0.055925717740595647],
                    '670':[0.0010890690698014294,0.04336828983700642],
                    '765':[0.00093810092188024833,0.026092951412422679],
                    '865':[0.0010906675048335769,0.02122474906498301]
                    }
        self.colDict = {'412':'#001166','443':'#004488','490':'#116688',
                        '510':'#228844','555':'#667722','670':'#aa2211',
                        '765':'#770500','865':'#440000'}
        super(SeaWiFSNoise,self).__init__(*args,**nargs)

def Main(argv):

    """ Expects baseline filepath/name and where pertubed files
        are stored.
    """
    __version__="2.0"
    parser = argparse.ArgumentParser(prog="AnalyticNoise")
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('baselinefile',  type=str, help=' input bseline file',metavar="BASELINEFILE")
    parser.add_argument('noiseDir',type=str, help='directory where pertubed files are stored',metavar="NOISEDIR")
    args = parser.parse_args()
    baseLineFile = args.baselinefile
    noisyDataDir = args.noiseDir
    
    noisySfx = '_wiggledBand_' #argv[2]
    swf = SeaWiFSNoise(baseLineFile,noisyDataDir,multiDimVar=True)
    swf.ReadFromSilent()
    swf.BuildUncs(noisySfx,verbose=True)
    swf.WriteToSilent()

if __name__ == "__main__":
    Main(sys.argv[1:])
