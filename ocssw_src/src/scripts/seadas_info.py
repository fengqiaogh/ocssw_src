#!/usr/bin/env python3

"""
Print various items that may be of interest for SeaDAS and OCSSW
troubleshooting or forum posts.
"""

from __future__ import print_function

import argparse
import os
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
import platform


__version__ = '0.0.3.devel'

def get_cleaned_python_version():
    """
    Return the Python version, cleaned of extra text (especially relevant for
    Anaconda installations. Return "Not found" if the Python version is not
    found (although how that happens is a major mystery!).
    """
    py_ver = 'Not found'
    py_ver = sys.version
    if sys.version.upper().find('ANACONDA') != -1:
        py_ver = ' '.join([sys.version.split(' ')[0], '(part of an Anaconda installation)'])
    else:
        py_ver = sys.version.split(' ')[0]
    return py_ver

def get_command_output(command):

    cmd_output = None
    try:
        process = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        cmd_output = stdout.decode("utf-8") + stderr.decode("utf-8")
    except:
        cmd_output = None
    return cmd_output

def get_env():
    """
    Returns env output, which is obtained by running env
    via subprocess.
    """
    env_output = None
    cmd = ["env"]
    output = get_command_output(cmd)
    if output:
        env_output = output.splitlines()
        env_output.sort()
    return env_output

def get_l_program_version(l_program, ocssw_root):
    """
    Returns the version for various "l programs" (e.g. l2gen, l3bin)
    """
    lprog_ver = 'program not found'
    exe_dirs = ['bin', 'run/bin', 'run/bin/linux_64']
    ver_args = ['version', '-version']
    exe_fnd = False
    prog = os.path.join(ocssw_root, l_program)
    if os.path.exists(prog):
        exe_fnd = True
    else:
        for cand_dir in exe_dirs:
            prog = os.path.join(ocssw_root, cand_dir, l_program)
            if os.path.exists(prog):
                exe_fnd = True
                break
    if exe_fnd:
        for ver_arg in ver_args:
            cmd_line = [prog, ver_arg]
            out_txt = get_command_output(cmd_line)
            if out_txt:
                ver_lines = out_txt.split('\n')
                ver_lines.remove('')
                lprog_ver = ver_lines[-1].strip()
                break
    return lprog_ver

def get_installed_mission(ocdata_root):
    missionInstalled = []
    missionListPath = os.path.join(ocdata_root, 'common/missionList.xml')
    if os.path.exists(missionListPath):
        tree = ET.parse(missionListPath)
        root = tree.getroot()
        for elem in root:
            name = elem.get('name')
            for subelem in elem:
                missionDir = subelem.text
                pathAbs = os.path.join(ocdata_root, missionDir) 
                if os.path.exists(pathAbs):
                    missionInstalled.append(name)
    return missionInstalled

def get_java_version():
    """
    Returns the java version, which is obtained by running java -version
    via subprocess.
    """
    java_ver = 'Not installed'
    cmd = ['java', '-version']
    output = get_command_output(cmd)
    if output and output.upper().find('RUNTIME') != -1:
        lines = output.split('\n')
        java_ver = re.sub('(^[^"]+|(?<=")[^"]+$)', '', lines[0]).strip('"')
    return java_ver

def get_os_distribution():
    """
    Returns the distribution of the operating system.
    """
    dist = ''
    rel_lines = None
    uname_info = os.uname()
    if uname_info[0] == 'Linux':
        with open('/etc/os-release', 'rt') as os_rel_file:
            rel_lines = os_rel_file.readlines()
        if rel_lines:
            rel_info_dict = dict()
            for line in rel_lines:
                line_parts = line.split('=')
                if len(line_parts) == 2:
                    if line_parts[0].strip():
                        rel_info_dict[line_parts[0].strip()] = line_parts[1].strip().strip('"')
            if 'PRETTY_NAME' in rel_info_dict:
                dist = rel_info_dict['PRETTY_NAME']
    elif uname_info[0] == 'Darwin':
        dist = 'MacOS ' + platform.mac_ver()[0]
    else:
        dist = ' '.join([uname_info[0], uname_info[2]])
    return dist

def get_python3_path():
    """
    Returns the python3 path, which is obtained by running which python3
    via subprocess.
    """
    python3_path = ''
    cmd = ['which', 'python3']
    output = get_command_output(cmd)
    if output:
        lines = output.split('\n')
        lines.remove('')
        python3_path = lines[-1].strip()
    return python3_path

def get_seadas_version(seadas_root):
    """
    Returns the SeaDAS version as held in the VERSION.txt file in the
    directory pointed to by the SEADAS_ROOT environment variable.
    """
    seadas_ver = 'Not installed'
    ver_path = os.path.join(seadas_root, 'VERSION.txt')
    if os.path.exists(ver_path):
        with open(ver_path, 'rt') as ver_file:
            in_line = ver_file.readline()
            if in_line.find('VERSION') != -1:
                seadas_ver = in_line.split(' ', 1)[1].strip()
    return seadas_ver

def handle_command_line():
    """
    Handle help or version being requested on the command line.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--AppDir", \
                        help="Application installation directory", \
                        type=str)
    args = parser.parse_args()
    return args

def print_seadas_info(ocssw_root, ocdata_root):
    """
    Print out the items of interest for the SeaDAS installation specified
    by the seadas_root and ocssw_root parameters.
    """
    print('NASA Science Processing (OCSSW):')
    print('OCSSWROOT={0}'.format(ocssw_root))
    print('OCDATAROOT={0}'.format(ocdata_root))
    print('l2gen version: {0}'.format(get_l_program_version('l2gen', ocssw_root).strip()))
    print('l2bin version: {0}'.format(get_l_program_version('l2bin', ocssw_root).strip()))
    print('l3bin version: {0}'.format(get_l_program_version('l3bin', ocssw_root).strip()))
    print('l3mapgen version: {0}'.format(get_l_program_version('l3mapgen', ocssw_root).strip()))
    print('Installed Missions: {0}'.format(get_installed_mission(ocdata_root)))
    print('')
    return

def print_sys_info():
    """
    Print out information about the system (OS, Python version, Java version).
    """
    print('General System and Software:')
    print('Operating system: {0}'.format(get_os_distribution()))
    print('Java version: {0}'.format(get_java_version()))
    print('Python3 version: {0}'.format(get_cleaned_python_version()))
    print('Python3 Path: {0}'.format(get_python3_path()))
    env_output = get_env()
    printenvs = ['CC','CXX','FC','OCSSW_ARCH','OCSSW_BIN','OCSSW_DEBUG',
                 'OCVARROOT','ELEMENTS','EOS_LIB_PREFIX','GCC_TUNE','HDFEOS_LIB',
                 'HRPT_STATION_IDENTIFICATION_FILE','L2GEN_ANC','LIB3_BIN','LIB3_CHECK',
                 'LIB3_DIR','LIB3_INC','LIB3_LIB','NAVCTL','NAVQC','ORBCTL','PGSINC','PGSLIB',
                 'OCTS_REGISTRATION_TABLES','PROJ_DATA','PROJ_LIB','SWFTBL','SWTBL']
    if env_output:
        print('Env:')
        for line in env_output:
            key, value = line.split('=', 1)
            if key in printenvs:
                print('   {0}={1}'.format(key,value))

    return

def main():
    """
    The program's main function - gets and prints items of interest.
    """
    if 'OCSSWROOT' in os.environ:
        ocssw_root = os.environ['OCSSWROOT']
    else:
        ocssw_root = 'Not defined'
    if 'OCDATAROOT' in os.environ:
        ocdata_root = os.environ['OCDATAROOT']   
    else:
        ocdata_root = 'Not defined'
    print_seadas_info(ocssw_root, ocdata_root)
    print_sys_info()
    return 0

if __name__ == '__main__':
    sys.exit(main())
