import os
from pathlib import Path
import logging
import sys


def env(self):
    """
    A simple module to populate some important environment variables
    """
    if os.getenv("OCSSWROOT") is None:
        scriptdir = Path(__file__).resolve().parent
        if 'src' in scriptdir.parts:
            self.dirs['packageroot'] =  scriptdir.parents[2]
        else:
            self.dirs['packageroot'] =  scriptdir.parents[1]
    else:
        self.dirs['packageroot'] = Path(os.getenv("OCSSWROOT"))

    self.dirs['root'] = self.dirs['packageroot'] / "share"
    if os.getenv("OCDATAROOT"):
        self.dirs['root'] = Path(os.getenv("OCDATAROOT"))

    if not self.dirs['root'].exists():
        print("ERROR: The OCDATAROOT {} directory does not exist.".format(self.dirs['root']))
        print("Please make sure you have downloaded and installed OCSSW and sourced the OCSSW environment")
        sys.exit(1)

    self.dirs['scripts'] = self.dirs['packageroot'] / "scripts"
    if os.getenv("OCSSW_SCRIPTS"):
        self.dirs['scripts'] = Path(os.getenv("OCSSW_SCRIPTS"))

    self.dirs['var'] = self.dirs['packageroot'] / "var"
    if os.getenv("OCVARROOT"):
        self.dirs['var'] = Path(os.getenv("OCVARROOT"))

    self.dirs['bin'] = self.dirs['packageroot'] / "bin"
    if os.getenv("OCSSW_BIN"):
        self.dirs['bin'] = Path(os.getenv("OCSSW_BIN"))

    self.dirs['bin3'] = self.dirs['packageroot'] / "opt" / "bin"
    if os.getenv("LIB3_BIN"):
        self.dirs['bin3'] = Path(os.getenv("LIB3_BIN"))

    self.dirs['log'] = self.dirs['var'] / "log"
    Path(self.dirs['log']).mkdir(parents=True, exist_ok=True)

    if os.getenv("OCSSW_DEBUG") is not None and int(os.getenv("OCSSW_DEBUG")) > 0:
        if not os.path.exists(self.dirs['bin']):
            print("Error:  OCSSW_DEBUG set, but...\n\t%s\ndoes not exist!" % self.dirs['bin'])
            sys.exit(1)
        else:
            print("Running debug binaries...\n\t%s" % self.dirs['bin'])

    self.dirs['run'] = str(Path.cwd())
    if self.curdir:
        self.dirs['anc'] = self.dirs['run']
    else:
        if self.ancdir is None:
            if os.getenv("L2GEN_ANC") is None and os.getenv("USER_L2GEN_ANC") is None:
                if self.verbose:
                    print("Neither the L2GEN_ANC nor USER_L2GEN_ANC environment variables are set.")
                    print("...using the current working directory for ancillary file download.")
                self.curdir = True
                self.dirs['anc'] = self.dirs['run']
            else:
                if os.getenv("L2GEN_ANC") is not None:
                    self.dirs['anc'] = os.getenv("L2GEN_ANC")
                else:
                    self.dirs['anc'] = os.getenv("USER_L2GEN_ANC")
        else:
            self.dirs['anc'] = Path(self.ancdir)

def build_executable_path(prog_name):
    """
    Returns the path to the program named in prog_name as a pathlib Path object.
    None is returned if the program is not found.
    """
    packageroot = None
    if os.getenv("OCSSWROOT"):
        packageroot = Path(os.getenv("OCSSWROOT"))

    prog_path = packageroot / 'bin' / prog_name
    if not Path.exists(prog_path):
        err_msg = "Cannot find program %s" % prog_name
        logging.error(err_msg)
        sys.exit(err_msg)

    return prog_path
