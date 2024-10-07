"""

SeaDAS library for commonly used functions within other python scripts

"""
import hashlib
import os
import sys
import re
import subprocess
from tarfile import BLOCKSIZE
import time
from datetime import datetime, timedelta, date
import logging
import requests
from requests.adapters import HTTPAdapter
from pathlib import Path

from seadasutils.MetaUtils import readMetadata


#  ------------------ DANGER -------------------
#
# The next 5 functions:
#    getSession
#    isRequestAuthFailure
#    httpdl
#    uncompressFile
#    get_file_time
#
# exist in two places:
#    OCSSWROOT/src/manifest/manifest.py
#    OCSSWROOT/src/scripts/seadasutils/ProcUtils.py
#
# Make sure changes get into both files.
#

DEFAULT_CHUNK_SIZE = 131072
BLOCKSIZE = 65536

# requests session object used to keep connections around
obpgSession = None

def getSession(verbose=0, ntries=5):
    global obpgSession

    if not obpgSession:
        # turn on debug statements for requests
        if verbose > 1:
            logging.basicConfig(level=logging.DEBUG)

        obpgSession = requests.Session()
        obpgSession.mount('https://', HTTPAdapter(max_retries=ntries))

        if verbose:
            print("OBPG session started")
    else:
        if verbose > 1:
            print("reusing existing OBPG session")

    return obpgSession

#  ------------------ DANGER -------------------
# See comment above
def isRequestAuthFailure(req) :
    ctype = req.headers.get('Content-Type')
    if ctype and ctype.startswith('text/html'):
        if "<title>Earthdata Login</title>" in req.text:
            return True
    return False

#  ------------------ DANGER -------------------
# See comment above
def httpdl(server, request, localpath='.', outputfilename=None, ntries=5,
           uncompress=False, timeout=30., verbose=0, force_download=False,
           chunk_size=DEFAULT_CHUNK_SIZE):

    status = 0
    urlStr = 'https://' + server + request

    global obpgSession
    localpath = Path(localpath)
    getSession(verbose=verbose, ntries=ntries)

    modified_since = None
    headers = {}

    if not force_download:
        if outputfilename:
            ofile = localpath / outputfilename
            modified_since = get_file_time(ofile)
        else:
            rpath = Path(request.rstrip())
            if 'requested_files' in request:
                rpath = Path(request.rstrip().split('?')[0]) 
            ofile = localpath / rpath.name
            if re.search(r'(?<=\?)(\w+)', ofile.name):
                ofile = Path(ofile.name.split('?')[0])

            modified_since = get_file_time(ofile)

        if modified_since:
            headers = {"If-Modified-Since":modified_since.strftime("%a, %d %b %Y %H:%M:%S GMT")}

    with obpgSession.get(urlStr, stream=True, timeout=timeout, headers=headers) as req:

        if req.status_code != 200:
            status = req.status_code
        elif isRequestAuthFailure(req):
            status = 401
        else:
            if not Path.exists(localpath):
                os.umask(0o02)
                Path.mkdir(localpath, mode=0o2775, parents=True)

            if not outputfilename:
                cd = req.headers.get('Content-Disposition')
                if cd:
                    outputfilename = re.findall("filename=(.+)", cd)[0]
                else:
                    outputfilename = urlStr.split('/')[-1]

            ofile = localpath / outputfilename

            # This is here just in case we didn't get a 304 when we should have...
            download = True
            if 'last-modified' in req.headers:
                remote_lmt = req.headers['last-modified']
                remote_ftime = datetime.strptime(remote_lmt, "%a, %d %b %Y %H:%M:%S GMT").replace(tzinfo=None)
                if modified_since and not force_download:
                    if (remote_ftime - modified_since).total_seconds() < 0:
                        download = False
                        if verbose:
                            print("Skipping download of %s" % outputfilename)

            if download:
                total_length = req.headers.get('content-length')
                length_downloaded = 0
                total_length = int(total_length)
                if verbose >0:
                    print("Downloading %s (%8.2f MBs)" % (outputfilename,total_length /1024/1024))

                with open(ofile, 'wb') as fd:

                    for chunk in req.iter_content(chunk_size=chunk_size):
                        if chunk: # filter out keep-alive new chunks
                            length_downloaded += len(chunk)
                            fd.write(chunk)
                            if verbose > 0:
                                percent_done = int(50 * length_downloaded / total_length)
                                sys.stdout.write("\r[%s%s]" % ('=' * percent_done, ' ' * (50-percent_done)))
                                sys.stdout.flush()

                if uncompress:
                    if ofile.suffix in {'.Z', '.gz', '.bz2'}:
                        if verbose:
                            print("\nUncompressing {}".format(ofile))
                        compressStatus = uncompressFile(ofile)
                        if compressStatus:
                            status = compressStatus
                else:
                    status = 0

                if verbose:
                    print("\n...Done")

    return status


#  ------------------ DANGER -------------------
# See comment above
def uncompressFile(compressed_file):
    """
    uncompress file
    compression methods:
        bzip2
        gzip
        UNIX compress
    """

    compProg = {".gz": "gunzip -f ", ".Z": "gunzip -f ", ".bz2": "bunzip2 -f "}
    exten = Path(compressed_file).suffix
    unzip = compProg[exten]
    cmd = [unzip,str(compressed_file.resolve())]
    p = subprocess.Popen(cmd, shell=False)
    status = os.waitpid(p.pid, 0)[1]
    if status:
        print("Warning! Unable to decompress %s" % compressed_file)
        return status
    else:
        return 0

#  ------------------ DANGER -------------------
# See comment above
def get_file_time(localFile):
    ftime = None
    localFile = Path(localFile)
    if not Path.is_file(localFile):
        while localFile.suffix in {'.Z', '.gz', '.bz2'}:
            localFile = localFile.with_suffix('')

    if Path.is_file(localFile):
        ftime = datetime.fromtimestamp(localFile.stat().st_mtime)

    return ftime

def cleanList(filename, parse=None):
    """
    Parses file list from oceandata.sci.gsfc.nasa.gov through html source
    intended for update_luts.py, by may have other uses
    """
    oldfile = Path.resolve(filename)
    newlist = []
    if parse is None:
        parse = re.compile(r"(?<=(\"|\')>)\S+(\.(hdf|h5|dat|txt))")
    if not Path.exists(oldfile):
        print('Error: ' + oldfile + ' does not exist')
        sys.exit(1)
    else:
        of = open(oldfile, 'r')
        for line in of:
            if '<td><a' in line:
                try:
                    newlist.append(parse.search(line).group(0))
                except Exception:
                    pass
        of.close()
        Path.unlink(oldfile)
        return newlist


def date_convert(datetime_i, in_datetype=None, out_datetype=None):
    """
    Convert between datetime object and/or standard string formats

    Inputs:
        datetime_i   datetime object or formatted string
        in_datetype  input format code;
                     must be present if datetime_i is a string
        out_datetype output format code; if absent, return datetime object

        datetype may be one of:
        'j': Julian     YYYYDDDHHMMSS
        'g': Gregorian  YYYYMMDDHHMMSS
        't': TAI        YYYY-MM-DDTHH:MM:SS.uuuuuuZ
        'h': HDF-EOS    YYYY-MM-DD HH:MM:SS.uuuuuu
    """

    # define commonly used date formats
    date_time_format = {
        'd': "%Y%m%d",  # YYYYMMDD
        'j': "%Y%j%H%M%S",  # Julian    YYYYDDDHHMMSS
        'g': "%Y%m%d%H%M%S",  # Gregorian YYYYMMDDHHMMSS
        't': "%Y-%m-%dT%H:%M:%S.%fZ",  # TAI YYYY-MM-DDTHH:MM:SS.uuuuuuZ
        'h': "%Y-%m-%d %H:%M:%S.%f",  # HDF-EOS YYYY-MM-DD HH:MM:SS.uuuuuu
    }
    if in_datetype is None:
        dateobj = datetime_i
    else:
        dateobj = datetime.strptime(datetime_i, date_time_format[in_datetype])

    if out_datetype is None:
        return dateobj
    else:
        return dateobj.strftime(date_time_format[out_datetype])


def addsecs(datetime_i, dsec, datetype=None):
    """
    Offset datetime_i by dsec seconds.
    """
    dateobj = date_convert(datetime_i, datetype)
    delta = timedelta(seconds=dsec)
    return date_convert(dateobj + delta, out_datetype=datetype)

def diffsecs(time0, time1, datetype=None):
    """
    Return difference in seconds.
    """
    t0 = date_convert(time0, datetype)
    t1 = date_convert(time1, datetype)
    return (t1-t0).total_seconds()

def round_minutes(datetime_i, datetype=None, resolution=5, rounding=0):
    """Round to nearest "resolution" minutes, preserving format.

    Parameters
    ----------
    datetime_i : string
        String representation of datetime, in "datetype" format
    datetype : string
        Format of datetime, as strftime or date_convert() code
    resolution : integer, optional
        Number of minutes to round to (default=5)
    rounding : integer, optional
        Rounding "direction", where
            <0 = round down
             0 = round to nearest (default)
            >0 = round up
    """
    dateobj = date_convert(datetime_i, datetype)

    if rounding < 0: # round down
        new_minute = (dateobj.minute // resolution) * resolution
    elif rounding > 0: # round up
        new_minute = (dateobj.minute // resolution + 1) * resolution
    else:  # round to nearest value
        new_minute = ((dateobj.minute + resolution/2.0) // resolution) * resolution

    # truncate to current hour; add new minutes
    dateobj -= timedelta(minutes=dateobj.minute,
                                  seconds=dateobj.second,
                                  microseconds=dateobj.microsecond)
    dateobj += timedelta(minutes=new_minute)

    return date_convert(dateobj, out_datetype=datetype)


def remove(file_to_delete):
    """
    Delete a file from the system
    A simple wrapper for Path.unlink
    """
    file_to_delete = Path(file_to_delete)
    if Path.exists(file_to_delete):
        Path.unlink(file_to_delete)
        return 0

    return 1


def ctime(the_file):
    """
    returns days since file creation
    """

    today = date.today().toordinal()
    p = Path(the_file)
    utc_create = time.localtime(p.stat().st_ctime)

    return today - date(utc_create.tm_year, utc_create.tm_mon, utc_create.tm_mday).toordinal()


def mtime(the_file):
    """
    returns days since last file modification
    """

    today = date.today().toordinal()
    p = Path(the_file)
    utc_mtime = time.localtime(p.stat().st_mtime)

    return today - date(utc_mtime.tm_year, utc_mtime.tm_mon, utc_mtime.tm_mday).toordinal()


def cat(file_to_print):
    """
    Print a file to the standard output.
    """
    with open(file_to_print) as f:
        print(f.read())

def compare_checksum(filepath,checksum):
    hasher = hashlib.sha1()
    with open(filepath, 'rb') as afile:
        buf = afile.read(BLOCKSIZE)
        while len(buf) > 0:
            hasher.update(buf)
            buf = afile.read(BLOCKSIZE)

    if hasher.hexdigest() == checksum:
        return False
    else:
        return True

def check_sensor(inp_file):
    """
    Determine the satellite sensor from the file metadata
    if unable to determine the sensor, return 'X'
    """

    senlst = {'Sea-viewing Wide Field-of-view Sensor (SeaWiFS)': 'seawifs',
              'SeaWiFS': 'seawifs',
              'Coastal Zone Color Scanner (CZCS)': 'czcs',
              'Ocean Color and Temperature Scanner (OCTS)': 'octs',
              'Ocean Scanning Multi-Spectral Imager (OSMI)': 'osmi',
              'Ocean   Color   Monitor   OCM-2': 'ocm2', 
              'Second-generation Global Imager (SGLI)': 'sgli',
              'GOCI': 'goci', 'Hawkeye': 'hawkeye', 'hico': 'hico',
              'OCIS': 'ocis', 'OCI': 'oci','MERIS': 'meris','MOS': 'mos', 'TM': 'tm',
              'Aquarius': 'aquarius', 'VIIRS': 'viirs'}


    fileattr = readMetadata(inp_file)
    if not fileattr:
        # sys.stderr.write('empty fileattr found in ' + inp_file + '\n')
        return 'X'
    if 'ASSOCIATEDPLATFORMSHORTNAME' in fileattr:
        print(fileattr['ASSOCIATEDPLATFORMSHORTNAME'])
        return fileattr['ASSOCIATEDPLATFORMSHORTNAME']
    elif 'Instrument_Short_Name' in fileattr:
        print(senlst[str(fileattr['Instrument_Short_Name'])])
        return senlst[str(fileattr['Instrument_Short_Name'])]
    elif 'Sensor' in fileattr:
        print(senlst[(fileattr['Sensor']).strip()])
        return senlst[(fileattr['Sensor']).strip()]
    elif 'Sensor name' in fileattr:
        print(senlst[(fileattr['Sensor name']).strip()])
        return senlst[(fileattr['Sensor name']).strip()]
    elif 'SENSOR_ID' in fileattr and re.search('(OLI|ETM)', fileattr['SENSOR_ID']):
        if 'SPACECRAFT_ID' in fileattr and re.search('LANDSAT_9', fileattr['SPACECRAFT_ID']):
            print('L9')
            return 'L9'
        else:
            print(fileattr['SENSOR_ID'].strip())
            return fileattr['SENSOR_ID'].strip()
    elif 'PRODUCT' in fileattr and re.search('MER', fileattr['PRODUCT']):
        print(fileattr['PRODUCT'])
        return 'meris'
    elif 'instrument' in fileattr:
        print(fileattr['instrument'])
        if re.search('(OLCI|MSI|VIIRS)', fileattr['instrument']):
            if 'platform' in fileattr:
                return fileattr['platform']
        else:
            return senlst[(fileattr['instrument'])].strip()
    else:
        return 'X'
