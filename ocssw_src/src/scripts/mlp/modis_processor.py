
"""
Class for handling the processing of MODIS data files.
"""

import modis.modis_L1A_utils
import seadasutils.ProcUtils
import os
import subprocess
import sys

__author__ = 'melliott'

class ModisProcessor():
    """
    A class for doing MODIS Processing
    """
    def __init__(self, l0_name, l1a=None, geo=None):
        self.ocssw_root = os.getenv('OCSSWROOT')
        self.l0_file = l0_name
        if l1a is not None:
            self.l1a = l1a
        else:
            self.l1a = self.derive_l1a_basename() + '.L1A_LAC'
        if geo is not None:
            self.geo = l1a
        else:
            self.geo = self.derive_l1a_basename() + '.GEO'

    def derive_l1a_basename(self):
        """
        Determine what the default basename for the L1A file (and GEO file) should be.
        """
        starttime = None
        create_constructor_cmd = ['l0cnst_write_modis ', self.l0_file]
        print('Running: ', create_constructor_cmd)
        try:
            l0cnst_write_output = subprocess.run(create_constructor_cmd, capture_output=True, text=True, shell=False)
            if l0cnst_write_output.returncode:
                print('Error! Could not run l0const')
                sys.exit(41)
            else:
                print('l0cnst_write_modis run complete')
                for line in l0cnst_write_output.stdout.splitlines():
                    try:
                        key, val = line.split('=')
                    except ValueError:
                        continue  # line doesn't contain '='
                    if 'starttime' in key:
                        self.start = val.strip()
                        break
        except OSError as ose:
            print(ose.errno, ose.strerror)

        if starttime is None:
            print('Error!  Could not determine start time.')
            sys.exit(42)
        granule_time = seadasutils.ProcUtils.date_convert(starttime, 't', 'j')
        return ''.join(['A', granule_time])

    def run_modis_geo(self, out_file=None):
        """
        Do the MODIS Geonavigation processing.
        """
        geo_cmd = [os.path.join(self.ocssw_root, 'bin', 'modis_GEO'), self.l1a]
        if out_file:
            geo_cmd.append(['-o',out_file])
        ret_val = subprocess.run(geo_cmd, shell=False).returncode
        return ret_val

    def run_modis_l1a(self, out_file=None):
        """
        Do the MODIS L1A processing.
        """
        l1a_cmd = [os.path.join(self.ocssw_root, 'bin', 'modis_L1A'), '-v', self.l0_file]
        if out_file:
            l1a_cmd.append(['-o',out_file])
        ret_val = subprocess.run(l1a_cmd, shell=False).returncode
        return ret_val

    def run_modis_l1b(self):
        """
        Do the MODIS L1B processing.
        """
        l1b_cmd = [os.path.join(self.ocssw_root, 'bin', 'modis_L1B'), self.l1a]
        ret_val = subprocess.run(l1b_cmd, shell=False).returncode
        return ret_val
