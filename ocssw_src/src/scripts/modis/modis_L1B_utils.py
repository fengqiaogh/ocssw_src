
import os
import sys
import subprocess
import shutil
from seadasutils.ParamUtils import ParamProcessing
import seadasutils.ProcUtils as ProcUtils
import mlp.get_obpg_file_type
import mlp.get_output_name_utils


class ModisL1B:
    """

    Runs l1bgen_modis for a single Level 1A granule
    A MODIS Aqua or Terra Level 1B granule is produced from a L1A file.

    This script is a wrapper for MODIS L1A to L1B processing:

        -calls the modis_timestamp binary to determine platform & start time
        -constructs the PCF file
        -calls the l1bgen_modisa or l1bgen_modist binary
        -checks exit status from the L1B processing and cleans up if no error

    """

    def __init__(self, inp_file=None,
                 parfile=None,
                 geofile=None,
                 okm=None,
                 hkm=None,
                 qkm=None,
                 obc=None,
                 lutver=None,
                 lutdir=None,
                 delfiles=0,
                 log=False,
                 verbose=False):
        self.filename = inp_file
        self.parfile = parfile
        self.geofile = geofile
        self.okm = okm
        self.hkm = hkm
        self.qkm = qkm
        self.obc = obc
        self.proctype = 'modisL1B'
        self.lutversion = lutver
        self.lutdir = lutdir
        self.delfiles = delfiles
        self.log = log
        self.ancdir = None
        self.pcf_file = None
        self.curdir = False
        self.verbose = verbose
        self.dirs = {}
        self.sensor = None
        self.prefix = None

        # version specific variables
        self.collection_id = '061'

        if self.parfile:
            if self.verbose:
                print("Reading parameter file: %s" % self.parfile)
            param_proc = ParamProcessing(parfile=self.parfile)
            param_proc.parseParFile(prog='l1bgen')
            phash = param_proc.params['l1bgen']
            for param in (list(phash.keys())):
                if not self[param]:
                    self[param] = phash[param]

        # determine geo file name
        if self.geofile is None:
            file_typer = mlp.get_obpg_file_type.ObpgFileTyper(self.filename)
            ftype, sensor = file_typer.get_file_type()
            stime, etime = file_typer.get_file_times()
            data_files_list = list([mlp.obpg_data_file.ObpgDataFile(self.filename,
                                                                   ftype,
                                                                   sensor,
                                                                   stime,
                                                                   etime)])
            self.geofile = mlp.get_output_name_utils.get_output_name(data_files_list, 'modis_GEO', None)

            if not os.path.exists(self.geofile):
                self.geofile = os.path.join(os.path.dirname(self.filename), self.geofile)

            if not os.path.exists(self.geofile):
                self.geofile = '.'.join([self.filename.split('.')[0], "GEO"])
                if not os.path.exists(self.geofile):
                    geofile_parts = self.filename.split('.')[:-1]
                    geofile_parts.append('GEO')
                    self.geofile = '.'.join(geofile_parts)
            
            if self.verbose:
                print("Assuming GEOFILE is %s" % self.geofile)


    def __setitem__(self, index, item):
        self.__dict__[index] = item

    def __getitem__(self, index):
        return self.__dict__[index]

    def chk(self):
        """
        check parameters
        """
        if not os.path.exists(self.filename):
            print("ERROR: File", self.filename, "does not exist.")
            sys.exit(1)

        if not os.path.exists(self.geofile):
            print("ERROR: File", self.geofile, "does not exist.")
            sys.exit(1)

    def _clean_files(self):
        """
        Removes any unwanted files from the L1A -> L1B processing.
        """
        if self.delfiles & 1:
            ProcUtils.remove(self.okm)
        if self.delfiles & 2:
            ProcUtils.remove(self.hkm)
        if self.delfiles & 4:
            ProcUtils.remove(self.qkm)
        if self.delfiles & 8:
            ProcUtils.remove(self.obc)

        if self.log is False:
            ProcUtils.remove(self.pcf_file)
            base = os.path.basename(self.okm)
            ProcUtils.remove(os.path.join(self.dirs['run'],
                                          '.'.join(['LogReport', base])))
            ProcUtils.remove(os.path.join(self.dirs['run'],
                                          '.'.join(['LogStatus', base])))
            ProcUtils.remove(os.path.join(self.dirs['run'],
                                          '.'.join(['LogUser', base])))

    def run(self):
        """
        Run l1bgen_modis (MOD_PR02)
        """

        if self.verbose:
            print("Processing MODIS L1A file to L1B...")
        l1bgen = os.path.join(self.dirs['bin'],
                              ''.join(['l1bgen_', self.sensor]))
        status = subprocess.run(l1bgen, shell=False).returncode
        if self.verbose:
            print(l1bgen, "exit status:", str(status))

        if not status:
            ProcUtils.remove(os.path.join(self.dirs['run'], "GetAttr.temp"))
            ProcUtils.remove(os.path.join(self.dirs['run'], "ShmMem"))
            for l1bfile in (self.okm, self.hkm, self.qkm, self.obc):
                if os.path.dirname(l1bfile) != self.dirs['run']:
                    origfile = os.path.join(self.dirs['run'],
                                            os.path.basename(l1bfile))
                    if os.path.exists(origfile):
                        shutil.move(origfile, l1bfile)
                        ProcUtils.remove('.'.join([origfile, 'met']))
                else:
                    ProcUtils.remove('.'.join([l1bfile, 'met']))

            self._clean_files()

            if self.verbose:
                print("MODIS L1B processing complete.")
        # else .. it failed
        else:
            print("ERROR: MODIS L1B processing failed.")
            print("Please examine the LogStatus and LogUser files for more information.")
            sys.exit(1)
