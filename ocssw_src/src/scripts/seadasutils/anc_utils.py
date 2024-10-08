
import os
import sys
import gc
import re
import subprocess
import json
from datetime import datetime
from datetime import timedelta

import xml.etree.ElementTree as ElementTree
from operator import sub
from collections import OrderedDict

import seadasutils.MetaUtils as MetaUtils
import seadasutils.ProcUtils as ProcUtils
from seadasutils.timestamp_utils import goci_timestamp
from seadasutils.timestamp_utils import hawkeye_timestamp
from seadasutils.timestamp_utils import hico_timestamp
from seadasutils.timestamp_utils import meris_timestamp
from seadasutils.timestamp_utils import msi_timestamp
from seadasutils.timestamp_utils import ocm2_timestamp
from seadasutils.timestamp_utils import olci_timestamp
from seadasutils.timestamp_utils import etm_timestamp
from seadasutils.timestamp_utils import sgli_timestamp
from seadasutils.timestamp_utils import tm_timestamp
from seadasutils.timestamp_utils import l9_timestamp
from modis.modis_utils import modis_timestamp
from viirs.viirs_utils import viirs_timestamp
from seadasutils.aquarius_utils import aquarius_timestamp


import seadasutils.ancDB as db
#import modules.ancDBmysql as db


DEFAULT_ANC_DIR_TEXT = "$OCVARROOT"


class getanc:
    """
    utilities for ancillary file search
    """

    def __init__(self, filename=None,
                 start=None,
                 stop=None,
                 ancdir=None,
                 ancdb='ancillary_data.db',
                 anc_file=None,
                 curdir=False,
                 atteph=False,
                 sensor=None,
                 opt_flag=None,
                 verbose=0,
                 printlist=True,
                 download=True,
                 timeout=60,
                 refreshDB=False,
                 use_filename=False):
        self.filename = filename
        self.start = start
        self.stop = stop
        self.ancdir = ancdir
        self.ancdb = ancdb
        self.anc_file=anc_file
        self.curdir = curdir
        self.opt_flag = opt_flag
        self.atteph = atteph
        self.dl = download
        self.refreshDB = refreshDB
        self.sensor = sensor
        self.dirs = {}
        self.files = {}
        self.printlist = printlist
        self.verbose = verbose
        self.timeout = timeout
        self.server_status = None
        self.db_status = None
        self.proctype = None
        self.use_filename=use_filename
        if self.atteph:
            self.proctype = 'modisGEO'

        self.query_site = "oceandata.sci.gsfc.nasa.gov"
        self.data_site = "oceandata.sci.gsfc.nasa.gov"

    def chk(self):
        """
        Check validity of inputs to
        """

        if self.start is None and self.filename is None:
            print("ERROR: No L1A_or_L1B_file or start time specified!")
            sys.exit(32)

        if self.atteph:
            if self.sensor is None and self.filename is None:
                print("ERROR: No FILE or MISSION specified.")
                sys.exit(1)
            if self.sensor is not None and self.sensor != "modisa" and self.sensor != "modist" \
                    and self.sensor.lower() != "aqua" and self.sensor.lower() != "terra":
                print("ERROR: Mission must be 'aqua', 'modisa', 'terra', or 'modist' ")
                sys.exit(16)

        if self.curdir is True and self.ancdir is not None:
            print("ERROR: The '--use-current' and '--ancdir' arguments cannot be used together.")
            print("       Please use only one of these options.")
            sys.exit(1)

        if self.start is not None:
            if (len(self.start) != 13 and len(self.start) != 19) or int(self.start[0:4]) < 1978 or int(self.start[0:4]) > 2030:
                print("ERROR: Start time must be in YYYYDDDHHMMSS or YYYY-MM-DDTHH:MM:SS format and YYYY is between 1978 and 2030.")
                sys.exit(1)

        if self.stop is not None:
            if (len(self.stop) != 13 and len(self.start) != 19) or int(self.stop[0:4]) < 1978 or int(self.stop[0:4]) > 2030:
                print("ERROR: End time must be in YYYYDDDHHMMSS or YYYY-MM-DDTHH:MM:SS format and YYYY is between 1978 and 2030.")
                sys.exit(1)

    @staticmethod
    def get_start_end_info(info):
        """
        Extracts and returns the start time, start date, end time, and end date
        from info. Returns a None value for any item not found.
        """
        starttime = None
        stoptime = None
        startdate = None
        stopdate = None
        for line in info[0].decode("utf-8").splitlines():
            if line.find("Start_Time") != -1:
                starttime = line.split('=')[1]
            if line.find("End_Time") != -1:
                stoptime = line.split('=')[1]
            if line.find("Start_Date") != -1:
                startdate = line.split('=')[1]
            if line.find("End_Date") != -1:
                stopdate = line.split('=')[1]
        return starttime, startdate, stoptime, stopdate

    @staticmethod
    def get_goci_time(goci_time_str):
        month_dict = {'JAN': '01', 'FEB': '02', 'MAR': '03', 'APR': '04',
                      'MAY': '05', 'JUN': '06', 'JUL': '07', 'AUG': '08',
                      'SEP': '09', 'OCT': '10', 'NOV': '11', 'DEC': '12'}
        time_str = goci_time_str[8:12] + '-' + \
                   month_dict[goci_time_str[4:7]] + '-' + \
                   goci_time_str[1:3] + 'T' + goci_time_str[13:15] + ':' + \
                   goci_time_str[16:18] + ':' + goci_time_str[19:21] + '.' + \
                   goci_time_str[22:25] + 'Z'
        return time_str

    @staticmethod
    def get_time_coverage_xml(elem):
        data_elem = elem.find('Data')
        value_elem = data_elem.find('DataFromFile')
        return value_elem.text.strip()

    def get_start_end_info_from_xml(self, raw_xml):
        """
        Extracts and returns the start time, start date, end time, and end date
        from the XML tree in raw_xml. Returns a None value for any item not
        found.
        """

        xml_root = ElementTree.fromstring(raw_xml)

        time_start_list = xml_root.findall('.//Attribute[@Name="time_coverage_start"]')
        if len(time_start_list) > 0:
            if len(time_start_list) > 1:
                print("Encountered more than 1 time_coverage_start tag. Using 1st value.")
            start = self.get_time_coverage_xml(time_start_list[0])
        else:
            time_start_list = xml_root.findall('.//Attribute[@Name="Scene Start time"]')
            if len(time_start_list) > 1:
                print("Encountered more than 1 Scene Start time tag. Using 1st value.")
            start_str = self.get_time_coverage_xml(time_start_list[0])
            start = self.get_goci_time(start_str)

        time_end_list = xml_root.findall('.//Attribute[@Name="time_coverage_end"]')
        if len(time_end_list) > 0:
            if len(time_end_list) > 1:
                print("Encountered more than 1 time_coverage_end tag. Using 1st value.")
            stop = self.get_time_coverage_xml(time_end_list[0])
        else:
            time_end_list = xml_root.findall('.//Attribute[@Name="Scene end time"]')
            if len(time_end_list) > 1:
                print("Encountered more than 1 Scene end time tag. Using 1st value.")
            stop_str = self.get_time_coverage_xml(time_end_list[0])
            stop = self.get_goci_time(stop_str)
        return start, stop

    def setup(self):
        """
        Set up the basics
        """
        # global stopdate, stoptime, startdate, starttime

        # set l2gen parameter filename
        if self.filename is None:
            if len(self.start) == 13:
                start_obj = datetime.strptime(self.start, '%Y%j%H%M%S')
                self.start = datetime.strftime(start_obj, '%Y-%m-%dT%H:%M:%S')
                if self.stop is not None:
                    stop_obj = datetime.strptime(self.stop, '%Y%j%H%M%S')
                    self.stop = datetime.strftime(stop_obj, '%Y-%m-%dT%H:%M:%S')
            if self.anc_file is None:
                self.base = self.start
            else:
                self.base = self.anc_file.replace(".anc", "")
            self.server_file = self.base + ".anc.server"    

            if self.atteph:
                self.anc_file = self.base + ".atteph"
            else:
                self.anc_file = self.base + ".anc"
        elif self.use_filename: 
            if self.anc_file is None:
                self.base = os.path.basename(self.filename)
                self.server_file = "{0:>s}.anc.server".format('.'.join(self.base.split('.')[0:-1]))
                if self.atteph:
                    self.anc_file = "{0:>s}.atteph".format('.'.join(self.base.split('.')[0:-1]))
                else:
                    self.anc_file = "{0:>s}.anc".format('.'.join(self.base.split('.')[0:-1]))
            else:
                self.base = self.anc_file.replace(".anc", "") 
                self.server_file = self.base + ".anc.server"

                if self.atteph:
                    self.anc_file = '.'.join([self.base, 'atteph'])
                else:
                    self.anc_file = '.'.join([self.base, 'anc'])

            if self.server_file == '.anc.server':
                self.server_file = self.filename + self.server_file
                self.anc_file = self.filename + self.anc_file        
        else:
            if self.anc_file is None:
                self.base = os.path.basename(self.filename)
                self.server_file = "{0:>s}.anc.server".format('.'.join(self.base.split('.')[0:-1]))
            else:
                self.base = self.anc_file.replace(".anc", "") 
                self.server_file = self.base + ".anc.server"
            
            if self.atteph:
                self.anc_file = '.'.join([self.base, 'atteph'])
            else:
                self.anc_file = '.'.join([self.base, 'anc'])

            if self.server_file == '.anc.server':
                self.server_file = self.filename + self.server_file
                self.anc_file = self.filename + self.anc_file

            # Check if start time specified.. if not, obtain from HDF header or if the
            # file doesn't exist, obtain start time from filename
            if self.start is None:
                # Check for existence of ifile, and if it doesn't exist assume the
                # user wants to use this script without an actual input file, but
                # instead an OBPG formatted filename that indicates the start time.
                if not os.path.exists(self.filename):
                    if self.sensor:
                        print("*** WARNING: Input file doesn't exist! Parsing filename for start time and setting")
                        print("*** end time to 5 minutes later for MODIS and 15 minutes later for other sensors.")
                        if len(self.base) < 14 or int(self.base[1:5]) < 1978 or int(self.base[1:5]) > 2030:
                            print("ERROR: Filename must be in XYYYYDDDHHMMSS format where X is the")
                            print("sensor letter and YYYY is between 1978 and 2030.")
                            sys.exit(1)
                        else:
                            self.start = self.base[1:14]
                    else:
                        print("*** ERROR: Input file doesn't exist and mission not set...bailing out...")
                        sys.exit(1)

                else:
                    # Determine start/end times from HDF file
                    # for l1info subsample every 250 lines
                    if self.verbose:
                        print("Determining pass start and end times...")
                    senchk = ProcUtils.check_sensor(self.filename)

                    if re.search('(Aqua|Terra)', senchk):
                        #   if self.mission == "A" or self.mission == "T":
                        self.start, self.stop, self.sensor = modis_timestamp(self.filename)
                    elif senchk.find("viirs|JPSS-1|JPSS-2|Suomi_NPP") == 0:
                        self.start, self.stop, self.sensor = viirs_timestamp(self.filename)
                    elif senchk.find("aquarius") == 0:
                        self.start, self.stop, self.sensor = aquarius_timestamp(self.filename)
                    elif re.search('goci', senchk):
                        self.start, self.stop, self.sensor = goci_timestamp(self.filename)
                    elif re.search('hawkeye', senchk):
                        self.start, self.stop, self.sensor = hawkeye_timestamp(self.filename)
                    elif senchk.find("hico") == 0:
                        self.start, self.stop, self.sensor = hico_timestamp(self.filename)
                    elif re.search('meris', senchk):
                        self.start, self.stop, self.sensor = meris_timestamp(self.filename)  
                    elif re.search('(S2A|S2B)', senchk):
                        self.start, self.stop, self.sensor = msi_timestamp(self.filename)
                    elif re.search('ocm2', senchk):
                        self.start, self.stop, self.sensor = ocm2_timestamp(self.filename)
                    elif re.search('ETM', senchk):
                        self.start, self.stop, self.sensor = etm_timestamp(self.filename)
                    elif re.search('(3A|3B)', senchk):
                        self.start, self.stop, self.sensor = olci_timestamp(self.filename)
                    elif re.search('sgli', senchk):
                        self.start, self.stop, self.sensor = sgli_timestamp(self.filename)
                    elif re.search('tm', senchk):
                        self.start, self.stop, self.sensor = tm_timestamp(self.filename)
                    elif re.search('L9', senchk):
                        self.start, self.stop, self.sensor = l9_timestamp(self.filename)
                    else:
                        if self.sensor is None:
                            self.sensor = senchk
                        mime_data = MetaUtils.get_mime_data(self.filename)
                        if MetaUtils.is_netcdf4(mime_data):
                            metadata = MetaUtils.dump_metadata(self.filename)
                            starttime, stoptime = self.get_start_end_info_from_xml(metadata)
                            starttime = starttime.strip('"')
                            stoptime = stoptime.strip('"')
                            starttime = starttime.strip("'")
                            stoptime = stoptime.strip("'")
                            if starttime.find('T') != -1:
                                self.start = starttime[0:19]
                                # self.start = ProcUtils.date_convert(starttime, 't', 'j')
                            else:
                                self.start = starttime[0:19]
                                # self.start = ProcUtils.date_convert(starttime, 'h', 'j')
                            if stoptime.find('T') != -1:
                                self.stop = stoptime[0:19]
                                # self.stop = ProcUtils.date_convert(stoptime, 't', 'j')
                            else:
                                self.stop = stoptime[0:19]
                                # self.stop = ProcUtils.date_convert(stoptime, 'h', 'j')
                            pass
                        else:
                            infocmd = [os.path.join(self.dirs['bin'], 'l1info'), '-s', '-i 250', self.filename]
                            l1info = subprocess.Popen(infocmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                            info = l1info.communicate()
                            (starttime, startdate, stoptime, stopdate) = self.get_start_end_info(info)
                            if not starttime or not startdate or not stoptime or not stopdate:
                                err_msg = 'ERROR: For ' + self.base + ' could not determine: '
                                if not starttime:
                                    err_msg = err_msg + ' start time'
                                if not startdate:
                                    if not starttime:
                                        err_msg = err_msg + ', start date'
                                    else:
                                        err_msg = err_msg + ' start date'
                                if not stoptime:
                                    if not starttime or not startdate:
                                        err_msg = err_msg + ', stoptime'
                                    else:
                                        err_msg = err_msg + ' stop time'
                                if not stopdate:
                                    if not starttime or not startdate or not stoptime:
                                        err_msg = err_msg + ', stop date'
                                    else:
                                        err_msg = err_msg + ' stop date'
                                err_msg = err_msg + '. Exiting.'
                                print(err_msg)
                                if info[1]:
                                    print("l1info reported the following error:")
                                    print('    {0}'.format(info[1]))
                                sys.exit(1)
                            else:
                                # self.start = ProcUtils.date_convert(startdate + ' ' + starttime, 'h', 'j')
                                # self.stop = ProcUtils.date_convert(stopdate + ' ' + stoptime, 'h', 'j')
                                self.start = startdate + 'T' + starttime[0:8]
                                self.stop = stopdate + 'T' + stoptime[0:8]

        if self.filename and not self.sensor and not self.use_filename:
            # Make sure sensor is set (JIRA Ticket #1012)
            self.sensor = ProcUtils.check_sensor(self.filename)

        if self.verbose:
            print()
            print("Input file: " + str(self.filename))
            print("Sensor    : " + str(self.sensor))
            print("Start time: " + str(self.start))
            print("End time  : " + str(self.stop))
            print()

    def set_opt_flag(self, key, off=False):
        """
        set up the opt_flag for display_ancillary_data (type=anc)
        opt_flag values:
           0 - just the basics, MET/OZONE
           1 - include OISST
           2 - include NO2
           4 - include ICE
        """
        optkey = {'sst': 1, 'no2': 2, 'ice': 4}

        if off:
            self.opt_flag = self.opt_flag - optkey[key]
        else:
            self.opt_flag = self.opt_flag + optkey[key]

    def finddb(self):
        """
        Checks local db for anc files.
        """
        print(self.ancdb)
        if len(os.path.dirname(self.ancdb)):
            self.dirs['log'] = os.path.dirname(self.ancdb)
            self.ancdb = os.path.basename(self.ancdb)

        if not os.path.exists(self.dirs['log']):
            self.ancdb = os.path.basename(self.ancdb)
            print('''Directory %s does not exist.
Using current working directory for storing the ancillary database file: %s''' % (self.dirs['log'], self.ancdb))
            self.dirs['log'] = self.dirs['run']

        self.ancdb = os.path.join(self.dirs['log'], self.ancdb)

        if not os.path.exists(self.ancdb):
            return 0

        ancdatabase = db.ancDB(dbfile=self.ancdb)
        if not os.path.getsize(self.ancdb):
            if self.verbose:
                print("Creating database: %s " % self.ancdb)
            ancdatabase.openDB()
            ancdatabase.create_db()
        else:
            ancdatabase.openDB()
            if self.verbose:
                print("Searching database: %s " % self.ancdb)

        if self.filename:
            filekey = os.path.basename(self.filename)
            if (re.search('manifest.safe', self.filename) or re.search('xfdumanifest.xml', self.filename)):
                # filekey = *.SAFE/maifest.safe
                filekey = os.path.join(os.path.basename(os.path.dirname(self.filename)), filekey) 
        else:
            filekey = None
        if self.start and len(self.start) == 13:
                start_obj = datetime.strptime(self.start, '%Y%j%H%M%S')
                self.start = datetime.strftime(start_obj, '%Y-%m-%dT%H:%M:%S')
                if self.stop is not None:
                    stop_obj = datetime.strptime(self.stop, '%Y%j%H%M%S')
                    self.stop = datetime.strftime(stop_obj, '%Y-%m-%dT%H:%M:%S')
        status = ancdatabase.check_file(filekey,starttime=self.start)
        if status:
            if not self.refreshDB:
                self.files = ancdatabase.get_ancfiles(filekey, self.atteph, starttime=self.start)
                if self.curdir:
                    for anckey in list(self.files.keys()):
                        self.files[anckey] = os.path.basename(self.files[anckey])
                self.db_status = ancdatabase.get_status(filekey, self.atteph, starttime=self.start)
                self.start, self.stop = ancdatabase.get_filetime(filekey, starttime=self.start)
            else:
                ancdatabase.delete_record(filekey, starttime=self.start)

            ancdatabase.closeDB()

        if status and self.db_status:
            if not self.refreshDB:
                if self.db_status > 0:
                    print("Warning! Non-optimal data exist in local repository.")
                    print("Consider re-running with the --refreshDB option to check for optimal ancillary files")
                return 1
            else:
                return 0
        else:
            return 0

    def findweb(self):
        """
        Execute the display_ancillary_files search and populate the locate cache database
        """

        # dlstat = 0
    
        with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'missionID.json'), 'r') as msn_file:
            msn = json.load(msn_file)

        # msn = {"seawifs": "6", "modisa": "7", "aqua": "7", "modist": "8", "terra": "8", "octs": "9",
        #        "czcs": "11", "suomi-npp": "14", "aquarius": "15", "meris": "19", "ocm2": "21", "hico" : "25", 
        #        "goci": "27", "oli": "28", "3a": "29", "jpss-1": "33","hawkeye": "35", "3b": "36"}

        ProcUtils.remove(self.server_file)

        # Query the OBPG server for the ancillary file list
        opt_flag = str(self.opt_flag)
        anctype = 'anc_data_api'
        # if self.atteph:
        #     opt_flag = ''
        #     anctype = 'atteph_test'

        if self.sensor == 'aquarius':
            opt_flag = ''

        msnchar = '0'
        if str(self.sensor).lower() in msn:
            # msnchar = msn[str(self.sensor).lower()]
            msnchar = msn[str(self.sensor.lower())]
        elif self.sensor and self.sensor.isdigit():
            msnchar = self.sensor

        if self.stop is None and not self.use_filename:
            start_obj = datetime.strptime(self.start, '%Y-%m-%dT%H:%M:%S')
            time_change = timedelta(minutes=5)
            stop_obj = start_obj + time_change
            self.stop = datetime.strftime(stop_obj,'%Y-%m-%dT%H:%M:%S')
            anc_str = '?&m=' + msnchar + '&s=' + self.start + '&e=' + self.stop + '&missing_tags=1'
            dlstat = ProcUtils.httpdl(self.query_site, '/'.join(['/api', anctype, anc_str]),
                                      os.path.abspath(os.path.dirname(self.server_file)),
                                      outputfilename=self.server_file,
                                      timeout=self.timeout,
                                      verbose=self.verbose
                                      )
        elif self.use_filename:
            anc_str = '?filename=' + os.path.basename(self.filename) + '&missing_tags=1'
            dlstat = ProcUtils.httpdl(self.query_site,
                                      '/'.join(['/api', anctype, anc_str]),
                                      os.path.abspath(os.path.dirname(self.server_file)),
                                      outputfilename=self.server_file,
                                      timeout=self.timeout,
                                      verbose=self.verbose
                                      )
        else: 
            anc_str = '?&m=' + msnchar + '&s=' + self.start + '&e=' + self.stop + '&missing_tags=1'
            dlstat = ProcUtils.httpdl(self.query_site,
                                      '/'.join(['/api', anctype, anc_str]),
                                      os.path.abspath(os.path.dirname(self.server_file)),
                                      outputfilename=self.server_file,
                                      timeout=self.timeout,
                                      verbose=self.verbose
                                      )
        gc.collect()

        if dlstat:
            print("Error retrieving ancillary file list")
            sys.exit(dlstat)

        with open(self.server_file, 'r') as data_file:
            results = json.load(data_file)
            self.db_status = int(results['status'])

            if self.db_status == 128:
                ProcUtils.remove(self.anc_file)
                print("The time duration provided is longer than 2 hours. please try with a shorter duration." )
                print("No parameter file created.")
                ProcUtils.remove(self.server_file)
                sys.exit(128)

            if self.db_status == 64:
                ProcUtils.remove(self.anc_file)
                print("Stop time is earlier than start time." )
                print("No parameter file created.")
                ProcUtils.remove(self.server_file)
                sys.exit(64)

            self.db_status = 0

            for f in results['files']:
                if (len(str(f[1]))) == 0:
                    if not self.atteph:
                        if re.search('MET', f[0]):
                           self.db_status = self.db_status | 1
                        if re.search('OZONE', f[0]):
                           self.db_status = self.db_status | 2
                        if re.search('sst', f[0]):
                           self.db_status = self.db_status | 4
                        if re.search('no2', f[0]):
                           self.db_status = self.db_status | 8
                        if re.search('ice', f[0]):
                           self.db_status = self.db_status | 16
                    if self.atteph:
                        if re.search('ATT', f[0]):
                           self.db_status = self.db_status | 4
                        if re.search('EPH', f[0]):
                            self.db_status = self.db_status | 8
                else:
                    if self.atteph:
                        if re.search('ATT|EPH', f[0]):
                            self.files[str(f[0]).lower()] = str(f[1])
                            if (f[2] != 1):
                                if re.search('ATT', f[0]):
                                    self.db_status = self.db_status | 1
                                if re.search('EPH', f[0]):
                                    self.db_status = self.db_status | 2
                    else: 
                        if re.search('MET|OZONE|sst|ice|AER|GEO', f[0]):
                            self.files[str(f[0]).lower()] = str(f[1])
            # self.db_status = int(results['status'])

        # FOR MET/OZONE:
        # For each anc type, DB returns either a zero status if all optimal files are
        # found, or different error statuses if not. However, 3 MET or OZONE files can be
        # returned along with an error status meaning there were one or more missing
        # files that were then filled with the file(s) found, and so though perhaps
        # better than climatology it's still not optimal. Therefore check for cases
        # where all 3 MET/ozone files are returned but status is negative and then
        # warn the user there might be more files to come and they should consider
        # reprocessing at a later date.
        #
        # DB return status bitwise values:
        # -all bits off means all is well in the world
        # -bit 0 = 1 - missing one or more MET
        # -bit 1 = 1 - missing one or more OZONE
        # -bit 2 = 1 - missing SST
        # -bit 3 = 1 - missing NO2
        # -bit 4 = 1 - missing ICE

        # FOR ATT/EPH:
        # all bits off means all is well in the world
        # -bit 0 = 1 - predicted attitude selected
        # -bit 1 = 1 - predicted ephemeris selected
        # -bit 2 = 1 - no attitude found
        # -bit 3 = 1 - no ephemeris found

        if self.server_status == 1 or dlstat or self.db_status is None:
            print("ERROR: The display_ancillary_files.cgi script encountered an error and returned the following text:")
            print()
            ProcUtils.cat(self.server_file)
            sys.exit(99)

        if self.db_status == 31 and not self.atteph:
            ProcUtils.remove(self.anc_file)
            print("No ancillary files currently exist that correspond to the start time " + self.start)
            print("No parameter file created (l2gen defaults to the climatologies).")
            ProcUtils.remove(self.server_file)
            sys.exit(31)

        if self.db_status == 12 and self.atteph:
            ProcUtils.remove(self.anc_file)
            print("No att/eph files currently exist that correspond to the start time " + self.start)
            ProcUtils.remove(self.server_file)
            sys.exit(12)

        # extra checks
        for f in (list(self.files.keys())):
            if not len(self.files[f]):
                print("ERROR: display_ancillary_files.cgi script returned blank entry for %s. Exiting." % f)
                sys.exit(99)

        ancdatabase = db.ancDB(dbfile=self.ancdb)

        if not os.path.exists(ancdatabase.dbfile) or os.path.getsize(ancdatabase.dbfile) == 0:
            ancdatabase.openDB()
            ancdatabase.create_db()
        else:
            ancdatabase.openDB()

        missing = []

        for anctype in self.files:
            if self.files[anctype] == 'missing':
                missing.append(anctype)
                continue
            if (self.filename and self.dl) or (self.start and self.dl):
                path = self.dirs['anc']
                if not self.curdir:
                    year, day = self.yearday(self.files[anctype])
                    path = os.path.join(path, year, day)

                if self.filename:
                    filekey = os.path.basename(self.filename)
                    if (re.search('manifest.safe', self.filename) or re.search('xfdumanifest.xml', self.filename)):
                        # filekey = *.SAFE/maifest.safe
                        filekey = os.path.join(os.path.basename(os.path.dirname(self.filename)), filekey) 
                else:
                    filekey = None
                ancdatabase.insert_record(satfile=filekey, starttime=self.start, stoptime=self.stop, anctype=anctype,
                                          ancfile=self.files[anctype], ancpath=path, dbstat=self.db_status,
                                          atteph=self.atteph)

        ancdatabase.closeDB()
        # remove missing items
        for anctype in missing:
            self.files.__delitem__(anctype)

    def yearday(self, ancfile):
        if ancfile.startswith('RIM_'):
            ymd = ancfile.split('_')[2]
            dt = datetime.strptime(ymd, '%Y%m%d')
            year = dt.strftime('%Y')
            day = dt.strftime('%j')
            return year, day

        if ancfile.startswith('MERRA'):
            ymd = ancfile.split('.')[4]
            dt = datetime.strptime(ymd, '%Y%m%d')
            year = dt.strftime('%Y')
            day = dt.strftime('%j')
            return year, day      
        
        # For YYYYmmddT or YYmmddHHMMSS
        if re.search('[\d]{8}T', ancfile) or re.search('[\d]{14}', ancfile): 
            ymd = re.search('[\d]{8}', ancfile).group()
            dt = datetime.strptime(ymd, '%Y%m%d')
            year = dt.strftime('%Y')
            day = dt.strftime('%j')
            return year, day   
        # For YYYYDDDHHMMSS
        elif re.search('[\d]{13}', ancfile):
            ymd = re.search('[\d]{7}', ancfile).group()
            year = ymd[0:4]
            day = ymd[4:7]
            return year, day
        # For YYYYmmddHHMM
        elif re.search('[\d]{12}', ancfile):
            ymd = re.search('[\d]{8}', ancfile).group()
            dt = datetime.strptime(ymd, '%Y%m%d')
            year = dt.strftime('%Y')
            day = dt.strftime('%j')
            return year, day
        # For YYYYDDDHHMM
        elif re.search('[\d]{11}', ancfile):
            ymd = re.search('[\d]{7}', ancfile).group()
            year = ymd[0:4]
            day = ymd[4:7]
            return year, day
        # For YYYYmmddHH
        elif re.search('[\d]{10}', ancfile):
            ymd = re.search('[\d]{8}', ancfile).group()
            dt = datetime.strptime(ymd, '%Y%m%d')
            year = dt.strftime('%Y')
            day = dt.strftime('%j')
            return year, day
        # For YYYYDDDHH
        elif re.search('[\d]{9}', ancfile):
            ymd = re.search('[\d]{7}', ancfile).group()
            year = ymd[0:4]
            day = ymd[4:7]
            return year, day
        # For YYYYmmdd
        elif re.search('[\d]{8}', ancfile):
            ymd = re.search('[\d]{8}', ancfile).group()
            dt = datetime.strptime(ymd, '%Y%m%d')
            year = dt.strftime('%Y')
            day = dt.strftime('%j')
            return year, day
        # For YYYYDDD
        elif re.search('[\d]{7}', ancfile):
            ymd = re.search('[\d]{7}', ancfile).group()
            year = ymd[0:4]
            day = ymd[4:7]
            return year, day

        if ancfile.startswith('GMAO'):
            ymd = ancfile.split('.')[1][0:8]
            dt = datetime.strptime(ymd, '%Y%m%d')
            year = dt.strftime('%Y')
            day = dt.strftime('%j')
            return year, day     

        if ancfile[:1].isdigit():
            ymd = ancfile[0:8]
            dt = datetime.strptime(ymd, '%Y%m%d')
            year = dt.strftime('%Y')
            day = dt.strftime('%j')
            return year, day

        if ancfile.startswith('SIF'):
            yyyyddd = ancfile.split('.')[0]
            offset = 3
        elif ancfile.startswith('PERT_'):
            yyyyddd = ancfile.split('_')[2]
            offset = 0
        elif self.atteph and not re.search(".(att|eph)$", ancfile):
            yyyyddd = ancfile.split('.')[1]
            offset = 1
        else:
            yyyyddd = ancfile
            offset = 1  # skip only 1st char

        year = yyyyddd[offset:offset + 4]
        day = yyyyddd[offset + 4:offset + 7]
        return year, day

    def locate(self, forcedl=False):
        """
        Find the files on the local system or download from OBPG
        """

        FILES = []
        for f in (list(self.files.keys())):
            if self.atteph:
                if re.search('scat|atm|met|ozone|file|geo|aer', f):
                    continue
            else:
                if re.search('att|eph', f):
                    continue
            FILES.append(os.path.basename(self.files[f]))

        dl_msg = 1

        for FILE in list(OrderedDict.fromkeys(FILES)):
            year, day = self.yearday(FILE)

            # First check hard disk...unless forcedl is set

            if self.curdir:
                self.dirs['path'] = '.'  # self.dirs['anc']
                if os.path.exists(FILE) and forcedl is False:
                    download = 0
                    if self.verbose:
                        print("  Found: %s" % FILE)
                else:
                    download = 1
            else:
                ancdir = self.dirs['anc']

                self.dirs['path'] = os.path.join(ancdir, year, day)
                if os.path.exists(os.path.join(ancdir, year, day, FILE)) and forcedl is False:
                    download = 0
                    if self.verbose:
                        print("  Found: %s/%s" % (self.dirs['path'], FILE))
                else:
                    download = 1

            # Not on hard disk, download the file

            if download == 1:
                if self.dl is False:
                    if dl_msg == 1:
                        if self.verbose:
                            print("Downloads disabled. The following missing file(s) will not be downloaded:")
                        dl_msg = 0
                    if self.verbose:
                        print("  " + FILE)
                else:

                    if self.verbose:
                        print("Downloading '" + FILE + "' to " + self.dirs['path'])
                    status = ProcUtils.httpdl(self.data_site, ''.join(['/ob/getfile/', FILE]),
                            self.dirs['path'], timeout=self.timeout, uncompress=True,force_download=forcedl,
                            verbose=self.verbose)
                    gc.collect()
                    if status:
                        if status == 304:
                            if self.verbose:
                                print("%s is not newer than local copy, skipping download" % FILE)
                        elif status == 401:
                            print("*** ERROR: Authentication Failure retrieving:")
                            print("*** " + '/'.join([self.data_site, 'ob/getfile', FILE]))
                            print("*** Please check that your ~/.netrc file is setup correctly and has proper permissions.")
                            print("***")
                            print("*** see: https://oceancolor.gsfc.nasa.gov/data/download_methods/")
                            print("***\n")
                        else:
                            print("*** ERROR: The HTTP transfer failed with status code " + str(status) + ".")
                            print("*** Please check your network connection and for the existence of the remote file:")
                            print("*** " + '/'.join([self.data_site, 'ob/getfile', FILE]))
                            print("***")
                            print("*** Also check to make sure you have write permissions under the directory:")
                            print("*** " + self.dirs['path'])
                            print()

                        if status != 304:
                            ProcUtils.remove(os.path.join(self.dirs['path'], FILE))
                            ProcUtils.remove(self.server_file)
                            sys.exit(1)

            for f in (list(self.files.keys())):
                if self.atteph:
                    if re.search('met|ozone|file|aer|geo', f):
                        continue
                else:
                    if re.search('att|eph', f):
                        continue
                if FILE == self.files[f]:
                    self.files[f] = os.path.join(self.dirs['path'], FILE)

    def write_anc_par(self):
        """
        create the .anc parameter file
        """
        ProcUtils.remove(self.anc_file)

        NONOPT = ""
        if not self.atteph:

            if self.sensor == 'aquarius':
                inputs = {'MET': {'bitval': 1, 'required': ['met1', 'met2', 'atm1', 'atm2']},
                          'SST': {'bitval': 4, 'required': ['sstfile1', 'sstfile2']},
                          'SeaIce': {'bitval': 16, 'required': ['icefile1', 'icefile2']},
                          'Salinity': {'bitval': 32, 'required': ['sssfile1', 'sssfile2']},
                          'XRAY': {'bitval': 64, 'required': ['xrayfile1']},
                          # 'SCAT': {'bitval': 128, 'required': ['scat']},
                          'TEC': {'bitval': 256, 'required': ['tecfile']},
                          'SWH': {'bitval': 512, 'required': ['swhfile1', 'swhfile2']},
                          'Frozen': {'bitval': 1024, 'required': ['frozenfile1', 'frozenfile2']},
                          'GEOS': {'bitval': 2048, 'required': ['geosfile']},
                          'ARGOS': {'bitval': 4096, 'required': ['argosfile1', 'argosfile2']},
                          'SIF': {'bitval': 8192, 'required': ['sif']},
                          'PERT': {'bitval': 16384, 'required': ['pert']},
                          'Matchup': {'bitval': 32768, 'required': ['sssmatchup']},
                          'Rainfall': {'bitval': 65536, 'required': ['rim_file']}}
                for anc in inputs:
                    if self.db_status & inputs[anc]['bitval']:
                        NONOPT = " ".join([NONOPT, anc])
                    else:
                        for ancfile in (inputs[anc]['required']):
                            if ancfile not in self.files:
                                NONOPT = " ".join([NONOPT, anc])
                                print('*** WARNING: No optimal {0} files found.'.format(ancfile))
                                break

            else:  # not aquarius

                if self.db_status & 1:
                    NONOPT = " ".join([NONOPT, 'MET'])
                else:
                    for key in (['met1', 'met2', 'met3']):
                        if key not in self.files:
                            NONOPT = " ".join([NONOPT, 'MET'])
                            print("*** WARNING: No optimal MET files found.")
                            break

                if self.db_status & 2:
                    NONOPT = " ".join([NONOPT, 'OZONE'])
                else:
                    for key in (['ozone1', 'ozone2', 'ozone3']):
                        if key not in self.files:
                            NONOPT = " ".join([NONOPT, 'OZONE'])
                            print("*** WARNING: No optimal OZONE files found.")
                            break

                if self.opt_flag & 1 and ('sstfile' not in self.files or (self.db_status & 4)):
                    NONOPT = " ".join([NONOPT, 'SST'])
                    print("*** WARNING: No optimal SST files found.")

                if self.opt_flag & 2 and ('no2file' not in self.files or (self.db_status & 8)):
                    NONOPT = " ".join([NONOPT, 'NO2'])
                    print("*** WARNING: No optimal NO2 files found.")

                if self.opt_flag & 4 and ('icefile' not in self.files or (self.db_status & 16)):
                    NONOPT = " ".join([NONOPT, 'Sea Ice'])
                    print("*** WARNING: No optimal ICE files found.")

        ancpar = open(self.anc_file, 'w')

        for key in sorted(self.files.keys()):
            if self.atteph:
                if key.find('att') != -1 or key.find('eph') != -1:
                    ancpar.write('='.join([key, self.files[key]]) + '\n')
            else:
                if key.find('att') == -1 and key.find('eph') == -1:
                    ancpar.write('='.join([key, self.files[key]]) + '\n')

        ancpar.close()

        if self.verbose:
            if self.atteph:
                if self.dl:
                    print("All required attitude and ephemeris files successfully determined and downloaded.")
                else:
                    print("All required attitude and ephemeris files successfully determined.")
            else:
                print()
                print("Created '" + self.anc_file + "' l2gen parameter text file:\n")

        if self.verbose or self.printlist:
            ProcUtils.cat(self.anc_file)

        if len(NONOPT):
            if self.db_status == 31:
                print("No optimal ancillary files were found.")
                print("No parameter file was created (l2gen defaults to the climatological ancillary data).")
                print("Exiting.")
                ProcUtils.remove(self.server_file)
            else:
                print()
                print("*** WARNING: The following ancillary data types were missing or are not optimal: " + NONOPT)
                if self.db_status == 3:
                    print("*** Beware that certain MET and OZONE files just chosen by this program are not optimal.")
                    print("*** For near real-time processing the remaining files may become available soon.")
                elif self.db_status == 1:
                    print("*** Beware that certain MET files just chosen by this program are not optimal.")
                    print("*** For near real-time processing the remaining files may become available soon.")
                elif self.db_status == 2:
                    print("*** Beware that certain OZONE files just chosen by this program are not optimal.")
                    print("*** For near real-time processing the remaining files may become available soon.")
        else:
            if self.verbose:
                if self.dl:
                    print()
                    print("- All optimal ancillary data files were determined and downloaded. -")
                else:
                    print()
                    print("- All optimal ancillary data files were determined. -")

    def cleanup(self):
        """
        remove the temporary 'server' file and adjust return status - if necessary
        """

        ProcUtils.remove(self.server_file)

        # if an anc type was turned off but its db_status bit was on, turn off the
        # status bit so the user (and GUI) won't think anything's wrong
        if not self.atteph:
            if self.db_status & 4 and self.opt_flag & 1:
                self.db_status = sub(self.db_status, 4)
            if self.db_status & 8 and self.opt_flag & 2:
                self.db_status = sub(self.db_status, 8)
            if self.db_status & 16 and self.opt_flag & 4:
                self.db_status = sub(self.db_status, 16)
