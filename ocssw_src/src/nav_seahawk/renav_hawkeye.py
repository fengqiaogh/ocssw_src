#!/usr/bin/env python

import argparse
import numpy as np
import os
import sys
import pyIGRF
from sgp4.api import jday 
from astropy.time import Time
from netCDF4 import Dataset
from hawknav.read_adcs import read_adcs_from_l1a
from hawknav.tle2orb import tle2orb
from hawknav.orb_interp import orb_interp
from hawknav.l_sun import l_sun
from hawknav.orb2lla import orb2lla
from hawknav.ned2ecr import ned2ecr
from hawknav.j2000_to_ecr import j2000_to_ecr
from hawknav.propagate import propagate
from hawknav.qinv import qinv
from hawknav.qprod import qprod
from hawknav.qmethod import qmethod
from hawknav.remake_orbit_objects import remake_orbit_objects
from hawknav.drop_orbit_objects import drop_orbit_objects
from hawknav.write_ncdf_data_object import write_ncdf_data_object

__version__ = '0.21.1_2021-12-07'

def renav_hawkeye(l1afile,l1arenav,tlefile,utcpolefile=None,magbiasfix=False,verbose=False):
    # Program to recompute Hawkeye navigation in L1A file using attitude sensor data
    # Input: 
    #   - l1afile, seahawk L1A file, in netcdf format
    #   - l1arenav, modified seahawk L1A file, in netcdf format
    #   - TLEline1 & TLEline2, TLE lines
    #   - utcpolefile, Earth motion file , 'utcpole.dat'
    # Output: new data in l1arenav file 
    # Ported from renav_hawkeye.pro by Fred Patt.
    # Liang Hong, 2/27/2020
    # Liang Hong, 3/4/2020, replaced geomag70 with pyIGRF lib. Note: original loadCoeffs.py in the package needs to be updated with a correction in LUT index
    # Liang Hong, 4/2/2020, reduced for loops to improve efficiency 
    # Liang Hong, 12/3/2020, check for a time misalignment between the Sun sensor and gyro/RW data; Sun vectors to correct RW rates; only FSSXP is used; 
    #                        initial polynomial correction to rwr using the gyro data
    # Liang Hong, 12/16/2020,  correction to the magnetometer data;  initial gyro bias to apply only to the data before the onboard gyro temperature correction was applied
    # Liang Hong, 2/19/2021, added missing group attributes to output l1arenav file
    # Liang Hong, 2/23/2021, added option for magnetometer bias correction on/off switch
    # Liang Hong, 4/21/2021, resolved image folding issue in FSS and antenna interference zone, generate FSS error and weighting arrays, avoid invalid FSS
    # Liang Hong, 8/4/2021, read image data year and date from time_coverage_start instead of image id from file name
    # Liang Hong, 12/7/2021, changed sun sensor to use calibration sun vector instead of calculated from FSS XP
    
    if not utcpolefile:
        utcpolefile = os.path.join(os.environ['OCVARROOT'],'modis','utcpole.dat')

    TLEline1 = None
    TLEline2 = None
    try:
        if verbose:
            print("reading %s" % tlefile)
        with open(tlefile,'r') as tle:
            TLEline1 = tle.readline()
            TLEline2 = tle.readline()
        TLEline1 = TLEline1.strip()
        TLEline2 = TLEline2.strip()
        
    except:
        print("Could not read TLE file.")
        sys.exit(110)
    
    nc_fid = Dataset(l1afile, 'r')

    # Read ADCS data from file
    print('Running renavigation (VER%s) on %s\n' % (__version__,l1afile))
    ptime,sunv,mag2,gyro,rw,stime,sunxp,sunzp,sunzn,atime,quat = read_adcs_from_l1a(l1afile)    # checked with .pro
    if len(ptime)<1:
        print('Error reading input L1A file.')
        sys.exit(110)

    # Get date from L1A file
    
    t_start= Time(nc_fid.getncattr('time_coverage_start'),format='isot') 
    iyr = int(t_start.yday[0:4])
    iday = int(t_start.yday[5:8])   # day of year
    
    t_start = t_start.unix - 30*60 # in UTC seconds, 30 minutes before image start
    t_end = Time(nc_fid.getncattr('time_coverage_end'),format='isot').unix + 30*60 # in UTC seconds, 30 minutes after image end

    # generate ECR orbit data from TLE
    t_interval = 60     # 1 minute interval in orbit data
    orb = tle2orb(TLEline1,TLEline2,t_start,t_end,t_interval)  # generate ECR data per 1 minunte interval
    
    # Check for time misalignment
    if (stime[0]-ptime[0] > 0.5):
        print('stime behind')
        ptime = ptime[1:]
        gyro = gyro[1:,:]
        rw = rw[1:,:]
        mag2 = mag2[1:,:]
        sunv = sunv[:-1,:]
        stime = stime[:-1]
        sunxp = sunxp[:-1,:]
        quat = quat[1:,:]
        
    if (ptime[0]-stime[0] > 0.5):
        print('ptime behind')
        stime = stime[1:]
        sunxp = sunxp[1:,:]
        ptime = ptime[:-1]
        sunv = sunv[1:,:]
        mag2 = mag2[:-1,:]
        gyro = gyro[:-1,:]
        rw = rw[:-1,:]
        quat = quat[:-1,:]
    
    nsen = np.size(ptime)
    nfss = np.size(stime)
    nsen = np.min([nsen,nfss])    

    # Interpolate orbit
    osec = (orb[:,0]-iday)*86400.0
    ptime_for_orb_interp = np.copy(ptime)
    if (osec[-1]>86400):
        # possible interpolation over a day
        ptime_for_orb_interp[np.where(ptime<osec[0])] += 86400
    pi,vi,orbfl = orb_interp(osec,orb[:,1:4],orb[:,4:7],ptime_for_orb_interp)
    lon,lat,hgt = orb2lla(pi)    # checked with .pro

    # Generate model mag field data
    hgt[np.where(hgt>600)] = 600     # geomag only works to 600 km
    nDecimalYear = Time(nc_fid.getncattr('time_coverage_start'),format='isot').decimalyear
    bmagg = []
    #nsen = np.size(ptime)
    for iRec in range(0,nsen):
        n2e = ned2ecr(lon[iRec],lat[iRec])
        geomagout = pyIGRF.igrf_value(lat[iRec],lon[iRec],hgt[iRec],nDecimalYear)[3:6]
        bmagg.append(np.dot(geomagout,n2e).tolist())
    magr = np.asarray(bmagg)
    
    # Generate model Sun vectors
    sunr,rs = l_sun(iyr,iday,ptime)
    sunr = np.transpose(sunr)
    
    # Transform model Sun and mag field vectors to ECI
    suni = np.zeros((nsen,3))
    magi = np.zeros((nsen,3))
    posi = np.zeros((nsen,3))
    veli = np.zeros((nsen,3))
    ecmat = j2000_to_ecr(iyr*np.ones(np.shape(ptime)),iday*np.ones(np.shape(ptime)),ptime,utcpolefile)
    suni = np.einsum('ij,ijk->ik',sunr,np.swapaxes(ecmat,1,2))
    magi = np.einsum('ij,ijk->ik',magr,np.swapaxes(ecmat,1,2))
    posi = np.einsum('ij,ijk->ik',pi,np.swapaxes(ecmat,1,2))
    veli = np.einsum('ij,ijk->ik',vi,np.swapaxes(ecmat,1,2))
    
    omegae = 7.29211585494e-5
    veli[:,0] = veli[:,0] - posi[:,1]*omegae
    veli[:,1] = veli[:,1] + posi[:,0]*omegae    
    
    # Determine which Sun sensor to use
    # sunb = np.copy(sunxp)
    sunb = np.copy(sunv)     # changed to use calibrated sun vector, due to FSS outage on Dec. 3, 2021
    
    # Correct magnetometer data
    # seahawk onboard ADPS implemented correction after Feb. 1, 2021 (1612137600)
    if (t_start<1612137600) or (magbiasfix==True):
        pm2 = np.array([ [0.969584,  0.022348,  0.077772],
       		[-0.024340,  0.909508, -0.003070],
       		[-0.059786,  0.082971,  0.907878] ])
        bm2 = [ -2.154e-06, -1.405e-06,  1.664e-06 ]
        mag2 = np.einsum('ij,ijk->ik',mag2,np.tile(pm2,[nsen,1,1]))
        mag2 = mag2 - np.tile(bm2,(nsen,1))
           
    # Use gyro data to correct RW rate data
    bmom = [11000,12000,3500]
    jd,fr = jday(iyr,1,iday,0,0,0)
    gbias = [0,0,0]
    if (jd < 2459191): # 2020-12-07 date of gyro temp correction
        gbias = [-0.0045,-0.0042,0.0045] 
    rwr = rw/bmom
    gyro = gyro - gbias
    fit = np.flipud(np.polyfit(np.arange(nsen),rwr-gyro,deg=2))
    rwr = rwr - np.transpose(np.polynomial.polynomial.polyval(np.arange(nsen), fit))
        
    # Use Sun sensor to correct rates
    cs = np.cross(sunb[1:,:],sunb[0:-1,:])
    cr = np.cross(sunb[0:-1,:],np.cross(rwr[0:-1,:],sunb[0:-1,:]))
    serr = np.zeros(nsen-1)+0.001
    swt = np.zeros(nsen)+100
    ind_err = np.argwhere(np.sum(np.square(cs-cr),1)>1e-5)
    serr[ind_err] = 1
    swt[ind_err] = 0
    swt[ind_err+1] = 0
    fit = np.flipud(np.polyfit(np.arange(nsen-1),cr-cs,deg=4,w=1/serr))
    rwr = rwr - np.transpose(np.polynomial.polynomial.polyval(np.arange(nsen),fit))
    
    # Propagate quaternions using RW-generated rates
    qbp = propagate(rwr)[0:-1,:]
        
    # Use propagated quaternions to transform vectors to start time
    magp = np.zeros((nsen,3))
    qv = np.zeros((nsen,4))
    qv[:,0:3] = sunb
    qi = qinv(qbp)
    qt1 = qprod(qbp,qv)
    qt2 = qprod(qt1,qi)
    sunp = qt2[:,0:3]
    qv[:,0:3] = mag2
    qt1 = qprod(qbp,qv)
    qt2 = qprod(qt1,qi)
    magp = qt2[:,0:3]
    
    # Normalize magnetometer data
    magpm = np.sqrt(magp[:,0]**2+magp[:,1]**2+magp[:,2]**2)
    magim = np.sqrt(magi[:,0]**2+magi[:,1]**2+magi[:,2]**2)
    magp = np.divide(magp.T,magpm).T
    magi = np.divide(magi.T,magim).T

    # Compute starting quaternion using q-method and propagate
    vecb = np.concatenate((sunp,magp))
    veci = np.concatenate((suni,magi))
    wgt = np.zeros(2*nsen)
    wgt[0:nsen] = swt
    wgt[nsen:2*nsen] = 1
    
    ind_invalid = np.argwhere(np.isnan(veci))
    if len(ind_invalid)>0:
        ind_invalid_row = np.unique(ind_invalid[:,0])
        vecb = np.delete(vecb,ind_invalid_row,0)
        veci = np.delete(veci,ind_invalid_row,0)
        wgt = np.delete(wgt,ind_invalid_row)
    q0 = qmethod(vecb,veci,wgt)
    qbb = np.zeros((nsen,4))
    qbb[0,:] = q0
    qbb[1:,:] = qprod(np.tile(q0,(nsen-1,1)),qbp[1:,:])
    
    # Re-order quaternion (scalar first)
    quati = np.zeros((nsen,4))
    quati[:,0] = qbb[:,3]
    quati[:,1:4] = qbb[:,0:3]
    
    # Generate extended orbit data array
    norb = nsen + 60
    otime = np.arange(norb) + ptime[0] - 30
    otime_for_orb_interp = np.copy(otime)
    if (osec[-1]>86400):
        # possible interpolation over a day
        otime_for_orb_interp[np.where(otime<osec[0])] += 86400
    pi,vi,orbfl = orb_interp(osec,orb[:,1:4],orb[:,4:7],otime_for_orb_interp)
    ecmat = j2000_to_ecr(iyr*np.ones(norb),iday*np.ones(norb),otime,utcpolefile)
    posi = np.einsum('ij,ijk->ik',pi,np.swapaxes(ecmat,1,2))
    veli = np.einsum('ij,ijk->ik',vi,np.swapaxes(ecmat,1,2))
    
    omegae = 7.29211585494e-5
    veli[:,0] = veli[:,0] - posi[:,1]*omegae
    veli[:,1] = veli[:,1] + posi[:,0]*omegae
    
    # Generate renav L1A file
    drop_orbit_objects(l1afile,l1arenav)

    nc_fid = Dataset(l1arenav, 'a')
    print('Writing updated navigation to %s\n'%l1arenav)

    ngid = nc_fid.groups['navigation_data']
    remake_orbit_objects(norb,nc_fid,ngid)
    
    # Navigation data
    write_ncdf_data_object(ngid,'att_time',ptime)
    write_ncdf_data_object(ngid,'att_quat',quati)
    write_ncdf_data_object(ngid,'orb_time',otime)
    write_ncdf_data_object(ngid,'orb_pos',posi)
    write_ncdf_data_object(ngid,'orb_vel',veli)

    # Close original NetCDF file.
    nc_fid.close()

if __name__ == "__main__":
    print("renav_hawkeye version %s"%__version__)
    # parse command line
    parser = argparse.ArgumentParser(prog='renav_hawkeye',
        description='Recompute Hawkeye navigation information for a L1A file using attitude sensor data')
    parser.add_argument('-v', '--verbose', help='print status messages',
                        action='store_true',default=False)
    parser.add_argument('ifile', help='name of L1A file to renavigate')
    parser.add_argument('ofile', help='output filename')
    parser.add_argument('tle', help='TLE file')  
    parser.add_argument('--utcpole',
                        help='Polar wander file',
                        default=None)
    parser.add_argument("--magbiasfix",
                        help='Correct magnetometer bias',
                        default=False,
                        dest='magbiasfix',
                        action='store_true')

    # parser.add_argument('--version', action='version', version="%(prog)s version " +__version__ )
    args = parser.parse_args()

    renav_hawkeye(args.ifile,args.ofile,args.tle,utcpolefile=args.utcpole,magbiasfix=args.magbiasfix,verbose=args.verbose)