#! /usr/bin/env python3

''' This is the aviris preprocessing script.
    Run this on the *img.hdr if it fails to run under l2gen
    or the aviris file converters.

    aviris_preproc <img file> <output file>

R.Healy 11/1/2016 (richard.healy@nasa.gov)

'''
__version__ = '1.0.0'

__author__ = 'rhealy'

import argparse
import sys, re
import glob, os
import numpy as np

def main():
    '''Specifies command line arguments and parses command line accordingly.'''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='''      Read an Aviris hdr file and generate a parameter file as input to the converter if needed.

            ''', add_help=True)
    parser.add_argument('-ifile', nargs=1, type=str, help=' Input Aviris hdr file ')
    parser.add_argument('-ofile', nargs=1, type=str, default=('aviris_preproc.hdr'),help=' output file ')

    #args = parser.parse_args('-ifile /accounts/rhealy/Downloads/aviris/oci/f100828t01p00r11rdn_b/f100828t01p00r11rdn_b_sc01_ort_img.hdr -ofile blah.hdr'.split())
    #args = parser.parse_args('-ifile /glusteruser/rhealy/aviris/f940921t01p02_r06c/f940921t01p02_r06c_img.hdr -ofile blah.hdr'.split())
    #args = parser.parse_args('-ifile /glusteruser/rhealy/aviris/f000411t01p03_r04c/f000411t01p03_r04c_img.hdr'.split())

    args=parser.parse_args()

    if not args.ifile:
        parser.error("you must specify an input file")

    ifile = ''.join(args.ifile)
    ofile = ''.join(args.ofile)
    print("HDRFILE =",checkHdrFile(ifile,ofile))


def checkHdrFile(ifile,ofile):
    avfields=['lines','samples','bands','interleave','rotation angle',
              'pixel size','Northing','Easting','UTM zone','wavelength','fwhm']
    musthave=['lines','samples','bands','fwhm','wavelength']

    avdict = {}
    avdir=os.path.dirname(ifile)
    avfile=os.path.basename(ifile)
    if avdir:
       try:
         os.chdir(avdir)
       except:
         print("Unable to change to ",avdir)
         sys.exit()

    for theLine in open(avfile,'r') :
        if '=' in theLine:
            fields = theLine.split('=')
            fields[0] = fields[0].strip()
            for f in avfields:
                pattern = re.compile(f)
                if re.search(pattern,fields[0]):
                    avdict[f] = fields[1]
                    #print(f, '=',avdict[f])
    if len(avdict) == len(avfields):
        return ifile

    for f in musthave:
        if f not in avdict.keys():
            print('Missing data...',f,avdict.keys())
            sys.exit()
    infofields=['pilot_st_lat','pilot_end_lat','pilot_st_long','pilot_end_long',
                'pilot_st_time','pilot_end_time']
    infodict = {}
    infofile=''.join(glob.glob("*info"))
    if not infofile:
        print("Missing a info file")
        sys.exit()
    for theLine in open(infofile,'r'):
        if '=' in theLine:
            fields = theLine.split('=')
            fields[0] = fields[0].strip()
            for f in infofields:
                pattern = re.compile(f)
                if re.search(pattern,fields[0]):
                    infodict[f] = fields[1]

    #for info,val in infodict.items():
        #print(info," = ",val)
    if len(infodict) < len(infofields):
        print("Incomplete information in info file:",infofile)
        print("Need these fields: ",infofields)
        print("But only have: ",infodict.keys())
        sys.exit()
    if len(infodict['pilot_st_time'].split()) > 1:
        time_st = parseMinSec(infodict['pilot_st_time'])
    else:
        time_st = np.float(infodict['pilot_st_time'])
    #print("time start=",time_st)
    if len(infodict['pilot_end_time'].split()) > 1:
        time_end = parseMinSec(infodict['pilot_end_time'])
    else:
         time_end = np.float(infodict['pilot_end_time'])

    #print("time end=",time_end)
    if len(infodict['pilot_st_lat'].split()) > 1:
        lat_st = parseMinSec(infodict['pilot_st_lat'])
    else:
         lat_st = np.float(infodict['pilot_st_lat'])
    #print("lat_st=",lat_st)
    if len(infodict['pilot_end_lat'].split()) > 1:
        lat_end = parseMinSec(infodict['pilot_end_lat'])
    else:
        lat_end = np.float(infodict['pilot_end_lat'])

    #print("lat_end=",lat_end)
    if len(infodict['pilot_st_long'].split()) > 1:
        lon_st = parseMinSec(infodict['pilot_st_long'])
    else:
        lon_st = np.float(infodict['pilot_st_long'])

    #print("lon_st=",lon_st)
    if len(infodict['pilot_st_long'].split()) > 1:
        lon_end = parseMinSec(infodict['pilot_end_long'])
    else:
        lon_end = np.float(infodict['pilot_end_long'])
    #print("lon_end=",lon_end)

    sign_lat = np.sign(lat_st)
    sign_lon = np.sign(lon_st)

    navfile=''.join(glob.glob("*nav"))
    #print("navfile <<",navfile,">>")
    if not navfile:
        print("Missing a navigation file")
        sys.exit()

    navout = ofile.split(".")[0]+"_nav.txt"
    fnavout = open(navout,"w")

    fnav = open(navfile,'r')
    line = fnav.readline()
    fnav.close()
    if len(line.split()) == 3:
        alt,lat,lon = line.split()
        # correct the sign in the nav file
        if np.sign(np.float(lat)) == sign_lat:
            sign_lat = 1.0
        if np.sign(np.float(lon)) == sign_lon:
            sign_lon = 1.0

        time_ar = np.arange(time_st,time_end,(time_end-time_st)/np.float(avdict['lines']))
        for theLine, avtime in zip(open(navfile,'r'), time_ar):
            alt,lat,lon = theLine.split()
            fnavout.write("{:.8f} {:.4f} {:.4f} {:.4f}\n".format(avtime,np.float(alt),sign_lat*np.float(lat),sign_lon*np.float(lon)))
    else:
        for theLine in open(navfile,'r'):
            nav_data = theLine.split()
            stime_ar = nav_data[3].split(':')
            avtime = np.float(stime_ar[1])+np.float(stime_ar[2])/60.+np.float(stime_ar[3])/3600.
            alt = nav_data[20]
            lat = nav_data[4]
            lon = nav_data[5]
            fnavout.write("{:.8f} {:.4f} {:.4f} {:.4f}\n".format(avtime,np.float(alt)/1000.,np.float(lat),np.float(lon)))

    fnavout.close()
    fhdr = open(ofile,"w")
    fhdr.write('AVIRIS_PREPROC_HDR\n')
    fhdr.write("hdrfile = {}\n".format(avfile))
    fhdr.write("navfile = {}\n".format(navout))
    gainfile=''.join(glob.glob("*gain"))
    if gainfile:
        fhdr.write("gainfile = {}\n".format(gainfile))
    imgfile=''.join(glob.glob("*img"))
    if imgfile:
        fhdr.write("imgfile = {}\n".format(imgfile))
    else:
        print("WARNING!! No img file\n")
        print("It is extremely likely the aviris reader will not work unless you add imgfile=the_name_of_the_img_file \n")
        print("into the file ",ofile)

    fhdr.close()
    return ofile

def parseMinSec(timeFld):
    time_stamp =timeFld.split()
    if len(time_stamp) < 3:
        time_stamp = timeFld.split()
        if len(time_stamp) < 3:
            print("Can't parse Time Field: ",timeFld)
            sys.exit()
    shr = ''.join(time_stamp[0].split(':'))
    smin = ''.join(time_stamp[1].split(':'))
    ssec = ''.join(time_stamp[2].split(':'))
    #print("hr={},min={},sec={}".format(shr,smin,ssec))
    hr = np.float(shr)
    return hr+np.sign(hr)*(np.float(smin)+np.float(ssec)/60)/60

def ParseCommandLine(args):
    '''Specifies command line arguments and parses command line accordingly.'''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='''      Read an Aviris hdr file and generate a parameter file as input to the converter if needed.

            ''', add_help=True)
    parser.add_argument('-ifile', nargs=1, type=str, help=' Aviris hdr file ')
    args = parser.parse_args('-ifile=/accounts/rhealy/Downloads/aviris/oci/f100828t01p00r11rdn_b/f100828t01p00r11rdn_b_sc01_ort_img.hdr'.split())

    parsedArgs = parser.parse_args(args)
    return parsedArgs


if __name__ == '__main__': main()
