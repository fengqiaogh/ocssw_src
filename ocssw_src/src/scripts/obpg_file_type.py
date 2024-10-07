#!/usr/bin/env python3

"""
A class for determining the OBPG type of a file.
"""
__version__ = '2.0'

__author__ = 'melliott'

import optparse
import os
import sys
import mlp.get_obpg_file_type

def get_usage_text():
    usage_text = \
        """usage: %prog [options] FILE_NAME [FILE_NAME ...]

  The following file types are recognized:
    Instruments: CZCS, GOCI, HICO, Landsat OLI, MODIS Aqua,
                 MODIS Terra, OCM2, OCTS, SeaWiFS, VIIRSN, VIIRSJ1, VIIRSJ2
    Processing Levels: L0 (MODIS only), L1A, L1B, L2, L3 binned,
                       L3 mapped """
    return usage_text

def process_command_line(cl_parser):
    """
    Uses optparse to get the command line options & arguments.
    """
    cl_parser.add_option('-t', '--times', action='store_true',
                         dest='times', default=False,
                         help='output start and end times for the file(s)')
    (opts, args) = cl_parser.parse_args()
    return opts, args

def main():
    """
    Main function to drive the program when invoked as a program.
    """
    use_msg = get_usage_text()
    ver_msg = ' '.join(['%prog', __version__])
    cl_parser = optparse.OptionParser(usage=use_msg, version=ver_msg)
    (opts, args) = process_command_line(cl_parser)

    if len(args) > 0:
        for arg in args:
            fname = arg
            file_typer = mlp.get_obpg_file_type.ObpgFileTyper(fname)
            (obpg_file_type, instrument) = file_typer.get_file_type()
            output = ''
            if obpg_file_type == 'unknown':
                if instrument != 'unknown':
                    output = '{0}: {1}: unknown'.format(
                                            os.path.basename(fname), instrument)
                else:
                    output = '{0}: unknown: unknown'.format(
                                os.path.basename(fname))
            else:
                if instrument != 'unknown':
                    output = '{0}: {1}: {2}'.format(os.path.basename(fname),
                                                    instrument, obpg_file_type)
                else:
                    output = '{0}: unknown: {1}'.format(
                                        os.path.basename(fname), obpg_file_type)
            if opts.times:
                if obpg_file_type != 'unknown' and instrument != 'unknown':
                    start_time, end_time = file_typer.get_file_times()
                    output += ': {0} : {1}'.format(start_time, end_time)
                else:
                    output += ': unable to determine file start and end times'
            print(output)
    else:
        print('\nError!  No file specified for type identification.\n')
        cl_parser.print_help()
    return 0

if __name__ == '__main__':
    sys.exit(main())
