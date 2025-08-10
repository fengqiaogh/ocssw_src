#!/usr/bin/env python3

import netCDF4 as nc4
import numpy as np
import numpy.ma as ma
from scipy import ndimage 
from scipy import signal
import sys
import re
import os,fnmatch
import argparse
from datetime import datetime, timezone

__version__ = '1.3 (2025-07-22)'

def find(pattern, path):
# function to find a file matching pattern in path
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

def make_striping_profile(npix):
#  function to generate an empirical profile that represents the 
#  OCI red band striping vs. pixel number.  The profile is bounded by 
#  two sinusoidal segments that are connected by a straight line.

    xs                         = (np.arange(38.) + 1.)/39.

    striping_profile           = np.empty((npix), dtype=np.float32)
    striping_profile[0:6]      = 0.
    striping_profile[6:70]     = np.sin((np.arange(64.)-30)*np.pi/60.)+1
    striping_profile[108:135]  = -0.775*(np.sin((np.arange(27.)-13.)*np.pi/28.)-1.)
    striping_profile[70:108]   = striping_profile[69]*(1.-xs) + striping_profile[108]*xs
    striping_profile[135:npix] = 0.

    striping_profile = striping_profile/2

    return striping_profile
    
def make_red_sender_band(pix1, pix2, iwave, woff, wsig, ipw, ilw, wred, rred, tri):

# IDL program to generate sender pseudo-band for OCI red band striping correction

#  Name		Type 	I/O	Description
#
#  pix1		int	 I	Pixel 1 index for receiver band
#  pix2		int	 I	Pixel 2 index for receiver band
#  iwave	int	 I	Index of receiver band wavelength
#  woff		float	 I	Offset from receiver band to center of spectral weighting function
#  wsig		float	 I	Sigma of spectral weighting function
#  ipw		int	 I	Width of pixel weighting function	
#  ilw		int	 I	Width of line weighting function	
#  wred		float	 I	Array of red band wavelengths
#  rred		float	 I	Array of red band reflectances (pixel, line, band)
#  rredr	float	 O	Array of red sender band reflectances (pixel, line)
#  tri          int      I      switch for triangular pixel weighting function
  
# Create output sender band array  
    rs     = rred.shape
    rredr  = np.zeros((rs[1], rs[2]), dtype=np.float32)
    nlines = rs[1]

# Create array of evenly spaced bands
#   Create index array
    irw    = np.concatenate((np.arange(17), 2*np.arange(29)+18, np.arange(11)+76, 2*np.arange(13)+88, np.arange(49)+114))
#  Resample input bands
    rredi  = rred[irw,:,:]
    wredi  = wred[irw]
    nbandi = irw.size

#  Generate spectral weighting function
    swgt   = np.exp(-((wredi - wred[iwave] - woff)/wsig)**2/2)
    swgt   = swgt/np.sum(swgt)
  
#  Generate pixel weighting function (square or triangular)
    if tri == 1 :
        print("triangular weighting function not implemented")
        exit(1)
    else:
        pwgt = np.tile(1.0,ipw)/ipw

    wgt = np.empty((nbandi, ilw, ipw), dtype=np.float32) 
    for i in range(ilw):
        wgt[:,i,:] = np.outer(swgt, pwgt)/ilw

# Compute weighted sender band
    tmp1    = rredi[:,:,pix1-100-ipw:pix2+ipw]
   #tmp2    = ndimage.correlate(tmp1, wgt)
    corr    = signal.correlate(tmp1, wgt, mode='valid' )
    tmp2    = np.zeros(tmp1.shape)
    xi      = tmp1.shape[0]//2
    zi      = (tmp1.shape[2]-corr.shape[2])//2
    zf      = zi + corr.shape[2]
    tmp2[xi,:,zi:zf] = corr[0,:,:]
    np.place(tmp2, tmp2 < -1, -32767)
    rredr[:,pix1-100-ipw:pix2+ipw] = tmp2[int(nbandi/2),:,:]

    return rredr
  

# program to compute and perform OCI red band HAM B striping correction

def main():
    print("correct_oci_ghosting", __version__)

    helpText = """
Return value
    0:   All OK
    1:   program error
    100: correction already done
    110: rhot_red does not exist
"""

    parser = argparse.ArgumentParser(
        description='OCI red band HAM B striping correction',
        epilog=helpText,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("l1bfile", help="l1b file to be corrected")
    parser.add_argument("--l1bfile2", dest="l1bfile2", help="second l1b file")
    parser.add_argument("--lut", dest="lutfile", help="path to non-default lut file")
    args = parser.parse_args()


    #l1bfiles=['striping_test/PACE_OCI.20240723T123140.L1B.V3.nc', 'striping_test/PACE_OCI.20240723T123140.L1B.V3.nc']
    #l1bfiles=['PACE_OCI.20240723T123140.L1B.V3.nc']
    #l1bfiles = []
    l1bfiles=[args.l1bfile]

    if args.lutfile:
        LUT_file = args.lutfile
    else:
        file = open(os.path.expandvars("$OCDATAROOT/oci/correct_oci_ghosting_defaults.par"))
        for line in file:
            if re.search("^lut", line):
                fields = line.split("=")
        try:
            LUT_file = os.path.expandvars(fields[1]).rstrip()
        except:
            print("default lut file not found")
            print("exiting")
            exit(1)


    if args.l1bfile2:
        l1bfiles.append(args.l1bfile2)


    print('Red band HAM B striping correction for ',l1bfiles[0])

    # Read striping LUT
    print('Reading LUT: ', LUT_file)

    #LUT_file = "IDL/OCI_striping_LUT_V1.nc"
    LUT_dataset = nc4.Dataset(LUT_file, mode='r', format='NETCDF4')

    send_band_off   = LUT_dataset['sender_band_offset'][:].data.item()
    send_band_sigma = LUT_dataset['sender_band_sigma'][:].data.item()
    send_pix_width  = LUT_dataset['sender_pixel_width'][:].data
    send_line_off   = LUT_dataset['sender_line_offset'][:].data.item()
    send_pix_off    = LUT_dataset['sender_pixel_offset'][:].data.item()
    str_start_pix   = LUT_dataset['striping_start_pixel'][:].data.item()
    str_end_pix     = LUT_dataset['striping_end_pixel'][:].data.item()
    str_start_band  = LUT_dataset['striping_start_band'][:].data.item()
    str_amp         = LUT_dataset['striping_amplitude'][:].data

    npixstr = str_end_pix - str_start_pix + 1

    LUT_dataset.close()

    # Read data from L1B file
    print('Reading red bands')

    l1b_dataset = nc4.Dataset(l1bfiles[0], mode='r', format='NETCDF4')


    # Check for correction already performed
    if 'striping_corrected' in l1b_dataset.groups['observation_data'].ncattrs():
        print('Striping correction already performed -- exiting')
        exit(100)

    try:
        rhot_red1 = l1b_dataset['observation_data']['rhot_red'][:]
        qred      = l1b_dataset['observation_data']['qual_red'][:]
    except Exception as e:
        print("-E- Could not read observation_data/rhot_red")
        exit(110)

    print('Reading HAM side and wavelengths')

    HAM_side = l1b_dataset['scan_line_attributes']['HAM_side'][:]
    red_wave = l1b_dataset['sensor_band_parameters']['red_wavelength'][:]


    # Get band and scan dimensions
    rbands  = l1b_dataset.dimensions['red_bands'].size
    nscans  = l1b_dataset.dimensions['scans'].size
    npixels = l1b_dataset.dimensions['pixels'].size
    k1 = np.nonzero(HAM_side)[0]
    nk1 = k1.size


    l1b_dataset.close()

    # Check for second file
    if len(l1bfiles) == 2:
        # Read additional scans and append to data
        l1b_dataset = nc4.Dataset(l1bfiles[-1], mode='r', format='NETCDF4')
        print('Reading red bands from file ',l1bfiles[-1])
        ##HERE
        num_scans = l1b_dataset['observation_data']['rhot_red'].shape[1]
        if(num_scans < 20):
            print("Accessory file has: ", num_scans, " scans.")
            print("Accessory file must have at least 20 scans")
            print("exiting.")
            exit(1)
        rhot_red2 = l1b_dataset['observation_data']['rhot_red'][0:rbands, 0:send_line_off, 0:npixels]
        rhot_red  = np.empty((rbands, nscans+send_line_off, npixels), dtype=np.float32)
        rhot_red[:,:nscans,:] = rhot_red1
        rhot_red[:,nscans:,:]    = rhot_red2
        l1b_dataset.close()
    else:
        # Reduce number of scans by line offset
        rhot_red = rhot_red1
        nscans   = nscans - send_line_off
        nk1      = int(nk1 - send_line_off/2)
        k1       = k1[0:nk1]

    # Get striping pixel profile
    str_prof  = make_striping_profile(npixstr)
    #str_prof1 = np.tile(str_prof, (nk1,1))  #never used
    str_corr  = np.zeros((rbands, nscans, npixstr), dtype=np.float32)

    # Loop through bands
    fill = -32767.0
    for iband in range(str_start_band, rbands):
        if(np.mod(iband, 10) == 0):
            print(iband)

        # Construct pseudo-band
        rhot_send = make_red_sender_band(str_start_pix, str_end_pix, iband, send_band_off, send_band_sigma, 
        send_pix_width[0], 1, red_wave, rhot_red, 0)

        # Loop through pixels
        for jpix in range(str_start_pix, str_end_pix+1):
            j = np.intersect1d(np.where(rhot_red1[iband,k1,jpix] != fill), np.where(rhot_send[k1+send_line_off,jpix+send_pix_off] != fill))
            str_corr[iband,k1[j],jpix-str_start_pix] = str_amp[iband]*str_prof[jpix-str_start_pix]*rhot_send[k1[j]+send_line_off,jpix+send_pix_off]

    # Apply correction and write data back to file
    rhot_red1[:,0:nscans,str_start_pix:str_end_pix+1] = rhot_red1[:,0:nscans,str_start_pix:str_end_pix+1] - str_corr
    print('Writing corrected data to ',l1bfiles[0])

    l1b_dataset = nc4.Dataset(l1bfiles[0], mode='a', format='NETCDF4')
    l1b_dataset['observation_data']['rhot_red'][:] = rhot_red1

    # clear bit 1
    qred[:,k1,str_start_pix:str_end_pix] = np.bitwise_and(qred[:,k1,str_start_pix:str_end_pix],253)
    l1b_dataset['observation_data']['qual_red'][:] = qred

    setattr(l1b_dataset['observation_data'] , 'striping_corrected', 'yes')

    # update the history global attribute
    l1b_dataset.history = l1b_dataset.history + "; " + datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S") + "Z: " + " ".join(sys.argv)

    l1b_dataset.close()

if __name__ == "__main__":
    main()
