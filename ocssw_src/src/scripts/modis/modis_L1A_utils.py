#! /usr/bin/env python3


import os
import sys
import subprocess
import shutil
from seadasutils.ParamUtils import ParamProcessing
import mlp.get_obpg_file_type
import mlp.obpg_data_file
import mlp.get_output_name_utils
import seadasutils.ProcUtils as ProcUtils


class modis_l1a:
    """
    This class defines the parameters and sets up the environment necessary
    to process a MODIS L0 granule to produce a L1A granule.  It is implemented
    by the modis_L1A code.  This class also contains the methods necessary to
    run the l0_write_construct and l1agen_modis binaries.
    """

    def __init__(self, filename=None,
                 parfile=None,
                 l1a=None,
                 nextgranule=None,
                 startnudge=10,
                 stopnudge=10,
                 satellite=None,
                 fix=True,
                 rounding=True,
                 lutver=None,
                 lutdir=None,
                 log=False,
                 verbose=True):
        self.filename = filename
        self.parfile = parfile
        self.l1a = l1a
        self.proctype = 'modisL1A'
        self.nextgranule = nextgranule
        self.stopnudge = stopnudge
        self.startnudge = startnudge
        self.sat_name = satellite
        self.verbose = verbose
        self.fix = fix
        self.l0fixed = False
        self.rounding = rounding
        self.lutversion = lutver
        self.lutdir = lutdir
        self.log = log
        self.ancdir = None
        self.curdir = False
        self.pcf_template = None
        self.start = None
        self.stop = None
        self.length = 0
        self.dirs = {}

        # version-specific variables
        self.collection_id = '061'
        self.pgeversion = '6.1.9'
#        self.lutversion = '0'

        if self.parfile:
            p = ParamProcessing(parfile=self.parfile)
            p.parseParFile(prog='l1agen')
            phash = p.params['l1agen']
            for param in (phash.keys()):
                if not self[param]:
                    self[param] = phash[param]
        self.l0file = os.path.basename(self.filename)

    def __setitem__(self, index, item):
        self.__dict__[index] = item

    def __getitem__(self, index):
        return self.__dict__[index]


    def chk(self):
        """
        check input parameters
        """
        if not os.path.exists(self.filename):
            print("ERROR: File: " + self.filename + " does not exist.")
            sys.exit(1)

        if self.nextgranule is not None:
            if not os.path.exists(self.nextgranule):
                print("ERROR: File " + self.nextgranule + " does not exist.")
                sys.exit(1)
            else:
                if self.verbose and self.stopnudge != 0:
                    print("* Next L0 granule is specified, therefore setting stopnudge = 0 *")
                self.stopnudge = 0

    def get_constructor(self):

        # rudimentary PCF file for logging, leapsec.dat
        self.pcf_file = self.l0file + '.pcf'
        ProcUtils.remove(self.pcf_file)

        os.environ['PGS_PC_INFO_FILE'] = self.pcf_file
        pcf = [line for line in open(self.pcf_template, 'r')]
        sed = open(self.pcf_file, 'w')
        for line in pcf:
            line = line.replace('LOGDIR', self.dirs['run'])
            line = line.replace('L1AFILE', os.path.basename(self.filename))
            line = line.replace('VARDIR', os.path.join(self.dirs['var'], 'modis'))
            sed.write(line)
        sed.close()

        # create constructor file
        if self.verbose:
            print("Determining pass start and stop time...\n")

        l0cnst_write_cmd = [os.path.join(self.dirs['bin'], 'l0cnst_write_modis'), self.l0file]

        l0cnst_write_output = subprocess.run(l0cnst_write_cmd, capture_output=True, text=True, shell=False)
        if l0cnst_write_output.returncode:
            return l0cnst_write_output.returncode
        else:
            for line in l0cnst_write_output.stdout.splitlines():
                try:
                    key, val = line.split('=')
                except ValueError:
                    continue  # line doesn't contain '='
                if 'starttime' in key:
                    self.start = val.strip()
                elif 'stoptime' in key:
                    self.stop = val.strip()
                elif 'length' in key:
                    self.length = val.strip()
            return 0


    def l1a_name(self):
    # determine output file name
        file_typer = mlp.get_obpg_file_type.ObpgFileTyper(self.filename)
        ftype, sensor = file_typer.get_file_type()
        stime, etime = file_typer.get_file_times()
        data_files_list = list([mlp.obpg_data_file.ObpgDataFile(self.filename,
                                                                ftype, sensor,
                                                                stime, etime)])
        self.l1a = mlp.get_output_name_utils.get_output_name(data_files_list, 'modis_L1A', None)

    def l0(self):
        """
        Write L0 Constructor File
        """

        # The L0 file and constructor file must reside in the same directory,
        # so create a symlink to the L0 file as needed.

        if not os.path.exists(self.l0file):
            os.symlink(self.filename, self.l0file)

        # create constructor file
        status = self.get_constructor()

        if status != 0 and status != 3:
            # bad error - exit now
            print("l0cnst_write_modis: Unrecoverable error encountered while attempting to generate constructor file.")
            sys.exit(1)

        if status == 3:
            # recoverable? try to fix l0
            if self.verbose:
                print("l0cnst_write_modis: A corrupt packet or a packet of the wrong size")
                print("                    was detected while generating the constructor file.")
            if not self.fix:
                print("Fixing of Level-0 file is currently disabled.")
                print("Please re-run without the '--disable-fix_L0' option.")
                sys.exit(1)

            # 1st call to l0fix_modis: get pass times
            if self.verbose:
                print("Attempting to fix the Level 0 file using l0fix_modis")

            l0fixcmd = [os.path.join(self.dirs['bin'], 'l0fix_modis'), self.l0file, '-1', '-1']
            if self.verbose:
                print(l0fixcmd)
            l0fix = subprocess.run(l0fixcmd, capture_output=True, text=True, shell=False)
            if l0fix.returncode:
                print("l0fix_modis: Unrecoverable error in l0fix_modis!")
                sys.exit(1)

            # 2nd call to l0fix_modis: fix packets
            for line in l0fix.stdout.splitlines():
                if "taitime_start:" in line:
                    self.taitime_start = line.split()[1]
                if "taitime_stop:" in line:
                    self.taitime_stop = line.split()[1]

            l0fixcmd = [os.path.join(self.dirs['bin'], 'l0fix_modis'), self.l0file, self.taitime_start, self.taitime_stop,
                 self.l0file + '.fixed']
            if self.verbose:
                print(l0fixcmd)
            l0fix = subprocess.run(l0fixcmd, capture_output=True, shell=False)
            if l0fix.returncode:
                print("l0fix_modis: Unrecoverable error in l0fix_modis!")
                sys.exit(1)
            if self.verbose:
                print("New Level 0 file successfully generated. Regenerating constructor file...")
            self.l0file += '.fixed'

            # try again to make constructor file
            status = self.get_constructor()
            if status:
                print("Failed to generate constructor file after running l0fix_modis.")
                print("Please examine your Level 0 file to determine if it is completely corrupt.")
                sys.exit(1)
            else:
                self.l0fixed = True

        # Determine pass start and stop times and duration of pass
        if self.l1a is not None:
            if self.verbose:
                print("Using specified output L1A filename: %s" % self.l1a)
        else:
            self.l1a_name()

        # Adjust end time to nominal 5-minute boundary to ensure processing of final scan
        if self.rounding:
            starttime = ProcUtils.round_minutes(self.start,'t',5)
            stoptime  = ProcUtils.round_minutes(self.stop, 't',5)
            if starttime == stoptime:
                if self.verbose:
                    print("Short granule: extending end time to 5-min boundary")
                stoptime  = ProcUtils.round_minutes(self.stop, 't',5,rounding=1)
            if starttime != stoptime:
                self.stop = stoptime  # don't change self.start

        # Adjust starttime, stoptime, and gransec for L0 processing
        self.start = ProcUtils.addsecs(self.start, self.startnudge, 't')
        self.stop = ProcUtils.addsecs(self.stop, -1 * self.stopnudge, 't')
        self.length = ProcUtils.diffsecs(self.start, self.stop, 't')
        self.granmin = str(float(self.length) / 60.0)

        # set output file name
        if self.verbose:
            print("Input Level 0:", self.filename)
            print("Output Level 1A:", self.l1a)
            print("Satellite:", self.sat_name)
            print("Start Time:", self.start)
            print("Stop Time:", self.stop)
            print("Granule Duration:", self.length, "seconds")
            print("")

    def run(self):
        """
        Run l1agen_modis (MOD_PR01)
        """

        # if next granule is set, create temp concatenated file
        if self.nextgranule is not None:
            shutil.move(self.l0file, self.l0file + '.tmp')
            cat = open(self.l0file, 'wb')
            shutil.copyfileobj(open(self.l0file + '.tmp', 'rb'), cat)
            shutil.copyfileobj(open(self.nextgranule, 'rb'), cat)
            cat.close()

        if self.verbose:
            print("Processing MODIS L0 file to L1A...")
        status = subprocess.run(os.path.join(self.dirs['bin'], 'l1agen_modis'), shell=False).returncode
        if self.verbose:
            print('l1agen_modis exit status:', str(status))

        # if next granule is set, move original L0 file back to original name
        if self.nextgranule is not None:
            shutil.move(self.l0file + '.tmp', self.l0file)

        # clean up log files
        if not status:
            ProcUtils.remove(self.l0file + '.constr')
            ProcUtils.remove(self.l0file + '.pcf')
            base = os.path.basename(self.l0file)
            ProcUtils.remove(os.path.join(self.dirs['run'], 'LogReport.' + base))
            ProcUtils.remove(os.path.join(self.dirs['run'], 'LogStatus.' + base))
            ProcUtils.remove(os.path.join(self.dirs['run'], 'LogUser.' + base))
            ProcUtils.remove(os.path.join(self.dirs['run'], "GetAttr.temp"))
            ProcUtils.remove(os.path.join(self.dirs['run'], "ShmMem"))
            if self.l0fixed:
                fixed_l0file = self.l0file.replace('.fixed','')
                ProcUtils.remove(fixed_l0file + '.constr')
                ProcUtils.remove(fixed_l0file + '.pcf')
                base = os.path.basename(fixed_l0file)
                ProcUtils.remove(os.path.join(self.dirs['run'], 'LogReport.' + base))
                ProcUtils.remove(os.path.join(self.dirs['run'], 'LogStatus.' + base))
                ProcUtils.remove(os.path.join(self.dirs['run'], 'LogUser.' + base))

            if os.path.dirname(self.l1a) != self.dirs['run']:
                origfile = os.path.join(self.dirs['run'],
                                            os.path.basename(self.l1a))
                if os.path.exists(origfile):
                        shutil.move(origfile, self.l1a)
                        ProcUtils.remove('.'.join([origfile, 'met']))
                else:
                    ProcUtils.remove('.'.join([self.l1a, 'met']))
            # ProcUtils.remove(self.l1a + '.met')

            if os.path.islink(self.l0file):
                os.remove(self.l0file)

            if os.path.islink(os.path.basename(self.filename)):
                os.remove(os.path.basename(self.filename))

            if self.log is False:
                ProcUtils.remove(self.pcf_file)
                base = os.path.basename(self.l1a)
                ProcUtils.remove(os.path.join(self.dirs['run'], 'LogReport.' + base))
                ProcUtils.remove(os.path.join(self.dirs['run'], 'LogStatus.' + base))
                ProcUtils.remove(os.path.join(self.dirs['run'], 'LogUser.' + base))

            if self.verbose:
                print("MODIS L1A processing complete.")
        else:
            print("modis_l1a: ERROR: MODIS L1A processing failed.")
            print("Please examine the LogStatus and LogUser files for more information.")
            sys.exit(1)
