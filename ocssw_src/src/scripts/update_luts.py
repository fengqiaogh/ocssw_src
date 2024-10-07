#! /usr/bin/env python3
"""
update_luts.py

Updates LUTS for the various sensors.
"""

import argparse
import seadasutils.LutUtils as Lut


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawTextHelpFormatter):
    pass


if __name__ == '__main__':

    version = '2.1'
    description = 'Retrieve latest lookup tables for specified sensor.'
    sensors = [ 'all', 'common', 'seawifs', 'hico', 'modisa', 'modist', 'viirsn', 'viirsj1', 'viirsj2', 'oci']
    platforms = ['aqua', 'terra', 'npp', 'j1', 'j2']

    # Define commandline options

    parser = argparse.ArgumentParser(prog='update_luts',formatter_class=CustomFormatter,
                                     description=description, add_help=True)

    parser.add_argument('mission', metavar='MISSION',
                        help='sensor or platform to process; one of:\n%(choices)s',
                        choices= sensors + platforms)

    parser.add_argument('-e', '--eval', action='store_true', dest='evalluts',
                        help='also download evaluation LUTs')

    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help='print status messages')

    parser.add_argument('-n', '--dry-run', action='store_true', dest='dry_run',
                        help='no action; preview files to be downloaded')

    parser.add_argument('--timeout', type=float, default=10,
                        help='network timeout in seconds')

    parser.add_argument('--version', action='version', version='%(prog)s ' + version)

    parser.add_argument('-d', '--debug', action='store_true',
                        help=argparse.SUPPRESS) # hidden option

    # Read options and take action
    args = parser.parse_args()
    if args.debug:
        import logging
        logging.basicConfig(level=logging.DEBUG,
                            format='%(levelname)s:%(message)s')

    # install all of the luts
    if args.mission == 'all':
        for tmpMission in sensors:
            if tmpMission != 'all':
                luts = Lut.LutUtils(verbose=args.verbose,
                                    mission=tmpMission,
                                    evalluts=args.evalluts,
                                    timeout=args.timeout,
                                    dry_run=args.dry_run)
                luts.get_luts()
                if luts.status != 0:
                    parser.exit(luts.status)
        parser.exit(luts.status)

    # always update the common directory
    luts = Lut.LutUtils(verbose=args.verbose,
                        mission='common',
                        evalluts=False,
                        timeout=args.timeout,
                        dry_run=args.dry_run)
    luts.get_luts()
    if luts.status != 0:
        parser.exit(luts.status)

    # stop here if only the common dir was requested
    if args.mission == 'common':
        parser.exit(luts.status)

    # on to the requested sensor
    luts = Lut.LutUtils(verbose=args.verbose,
                        mission=args.mission,
                        evalluts=args.evalluts,
                        timeout=args.timeout,
                        dry_run=args.dry_run)

    luts.get_luts()
    parser.exit(luts.status)
