#!/usr/bin/env python3

"""
Rewrite the seadas configuration file (usually 
$SEADAS/config/seadas.config), changing the value
of the seadas.ocssw.root line to a value passed
in via the command line call.
"""
import argparse
import shutil
import sys

def rewrite_seadas_config(seadas_config, install_dir):
    """
    Reads the input file then writes it to a temporary file,
    changing the appropriate entry.  Then it copies the temporary
    file over the original.
    """
    with open(seadas_config, 'rt') as in_file:
        in_content = in_file.readlines()

    with open('seadas_config.tmp', 'wt') as out_file:
        for line in in_content:
            if line.startswith('seadas.ocssw.root'):
                out_line = 'seadas.ocssw.root = ' + install_dir + '\n'
            else:
                out_line = line
            out_file.write(out_line)
    shutil.copy('seadas_config.tmp', seadas_config)
    return

def main():
    version = "1.0"
    parser = argparse.ArgumentParser(prog="rewrite_seadas_config")
    parser.add_argument('--version', action='version', version='%(prog)s ' + version)
    parser.add_argument("seadasConfig", 
                      help="SeaADAS config file to edit", metavar="CONFIGFILE")  
    parser.add_argument("OCSSWROOT", 
                      help="OCSSW root directory", metavar="OCSSWROOT")  

    args = parser.parse_args()
    if args.seadasConfig is None or args.OCSSWROOT is None:
        parser.print_help()
        sys.exit(1)

    rewrite_seadas_config(args.seadasConfig, args.OCSSWROOT)
    return 0

if __name__ == '__main__':
    sys.exit(main())
