#!/usr/bin/env python3
# coding: utf-8

"""
A Perl script to create and output satellite matchups from a SeaBASS file given an 
OB.DAAC L2 satellite SST file and a valid SeaBASS file containing
lat, lon, date, and time as /field entries.
written by J.Scott on 2018/07/20 (joel.scott@nasa.gov)
"""

def main():

    import argparse
    import os
    import re
    import subprocess
    from datetime import datetime, timedelta
    from copy import copy
    from math import isnan
    from collections import OrderedDict
    from seabass.SB_support import readSB

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,description='''\
      This program create and output satellite matchups from a given SeaBASS file.

      REQUIRED inputs:
          1) --sat_file        an OB.DAAC L2 satellite SST file(s)
          2) --seabass_file    a valid SeaBASS file(s) with latitude, longitude, and date-time information as field entries.

      Outputs:
          1) the original SeaBASS data
          AND
          2) collocated satellite products as additional columns appended to --seabass_file

      Example usage call:
         mk_matchup.py --sat_file=[file name].nc --seabass_file=[file name].sb

      Caveats:
        * This script is designed to work with files that have been properly
          formatted according to SeaBASS guidelines (i.e. Files that passed FCHECK).
          Some error checking is performed, but improperly formatted input files
          could cause this script to error or behave unexpectedly. Files
          downloaded from the SeaBASS database should already be properly formatted, 
          however, please email seabass@seabass.gsfc.nasa.gov and/or the contact listed
          in the metadata header if you identify problems with specific files.

        * It is always HIGHLY recommended that you check for and read any metadata
          header comments and/or documentation accompanying data files. Information 
          from those sources could impact your analysis.

        * Compatibility: This script was developed for Python 3.6.

      License:
        /*=====================================================================*/
                         NASA Goddard Space Flight Center (GSFC) 
                 Software distribution policy for Public Domain Software

         The fd_matchup.py code is in the public domain, available without fee for 
         educational, research, non-commercial and commercial purposes. Users may 
         distribute this code to third parties provided that this statement appears
         on all copies and that no charge is made for such copies.

         NASA GSFC MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THE SOFTWARE
         FOR ANY PURPOSE. IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED
         WARRANTY. NEITHER NASA GSFC NOR THE U.S. GOVERNMENT SHALL BE LIABLE FOR
         ANY DAMAGE SUFFERED BY THE USER OF THIS SOFTWARE.
        /*=====================================================================*/
      ''',add_help=True)

    parser.add_argument('--sat_file', nargs='+', type=argparse.FileType('r'), required=True, help='''\
      REQUIRED: input OB.DAAC Level-2 satellite netCDF file
      ''')

    parser.add_argument('--seabass_file', nargs='+', type=argparse.FileType('r'), required=True, help='''\
      REQUIRED: input SeaBASS file
      Must be a valid SeaBASS file, passing FHCHECK with no errors.
      Matched-up satellite variables will be appended as additional fields to the data matrix and relevant headers.
      File must contain latitude and longitude and date-time expressed as FIELD entries.
      ''')

    parser.add_argument('--box_size', nargs=1, default=([5]), type=int, help=('''\
      OPTIONAL: box size of the satellite data extract made around the in situ point
      Valid values are odd numbers between 3 and 11, default = 5
      '''))

    parser.add_argument('--min_valid_sat_pix', nargs=1, default=([50.0]), type=float, help=('''\
      OPTIONAL: percent minimum valid satellite pixels required to create an extract
      Valid value: (0.0 - 100.0), default = 50.0
      '''))

    parser.add_argument('--max_time_diff', nargs=1, default=([0.5]), type=float, help=('''\
      OPTIONAL: maximum time difference between satellite and in situ point
      Valid value: decimal number of hours (0 - 36 hours), default = 3
      '''))

    parser.add_argument('--verbose', default=False, action='store_true', help=('''\
      OPTIONAL: Displays reason for failed matchup for each in situ target called.
      '''))

    parser.add_argument('--no_header_comment', default=False, action='store_true', help=('''\
      OPTIONAL: Flag to NOT append exclusion criteria to the OFILE header. Useful when running script repeatedly. 
      '''))

    args=parser.parse_args()
    dict_args=vars(args)

    # input verification
    if ((dict_args["box_size"][0] % 2) == 0) or (dict_args["box_size"][0] > 11) or (dict_args["box_size"][0] < 3):
        parser.error("invalid --box_size specified, must be an ODD integer between 3 and 11")

    if (dict_args["min_valid_sat_pix"][0] > 100.0) or (dict_args["min_valid_sat_pix"][0] < 0.0):
        parser.error("invalid --min_valid_sat_pix specified, must be a percentage expressed as a floating point number between 0.0 and 100.0")

    if (dict_args["max_time_diff"][0] > 36) or (dict_args["max_time_diff"][0] < 0):
        parser.error("invalid --max_time_diff specified, must be a decimal number between 0 and 36")
    else:
        twin_Hmin = -1 * int(dict_args['max_time_diff'][0])
        twin_Mmin = -60 * (dict_args['max_time_diff'][0] - int(dict_args['max_time_diff'][0]))
        twin_Hmax = 1 * int(dict_args['max_time_diff'][0])
        twin_Mmax = 60 * (dict_args['max_time_diff'][0] - int(dict_args['max_time_diff'][0]))

    for filein_sb in dict_args['seabass_file']:

        # read and verify SeaBASS file and required fields
        if os.path.isfile(filein_sb.name):
            ds = readSB(filename=filein_sb.name, 
                        mask_missing=False, 
                        mask_above_detection_limit=False, 
                        mask_below_detection_limit=False, 
                        no_warn=True)
        else:
            parser.error('ERROR: invalid --seabass_file specified; does ' + filein_sb.name + ' exist?')

        ds.datetime = ds.fd_datetime()
        if not ds.datetime:
            parser.error('missing fields in SeaBASS file -- file must contain a valid FIELDS combination of date/year/month/day/sdy and time/hour/minute/second')

        print('Looking for satellite/in situ match-ups for:',filein_sb.name)

        for filein_sat in dict_args['sat_file']:

            if not re.search('\.nc', filein_sat.name.lower()) or \
               not re.search('l2', filein_sat.name.lower()) or \
               not re.search('sst', filein_sat.name.lower()):
                parser.error("invalid --sat_file specified, must be a Level-2 (L2) OB.DAAC SST netCDF (nc) file")
            else:
                #set l2_flags to check for OC/IOP versus SST/SST4 product suites
                flag_arg = ' ignore_flags=LAND\ NAVFAIL\ NAVWARN' + \
                           ' count_flags=LAND\ NAVFAIL\ NAVWARN'

            print('Checking:',filein_sat.name)
            write_flag = 0

            # loop through input SeaBASS file data rows
            for row,dt in enumerate(ds.datetime):

                # create time range of satellite obs to extract
                tim_min = dt + timedelta(hours=twin_Hmin,minutes=twin_Mmin)
                tim_max = dt + timedelta(hours=twin_Hmax,minutes=twin_Mmax)

                # verify lat/lon inputs from file
                try:
                    ds.lon = [float(i) for i in ds.data['lon']]
                    ds.lat = [float(i) for i in ds.data['lat']]
                except Exception as E:
                    print(E)
                    parser.error('Missing fields in SeaBASS file. File must contain lat and lon as fields, or specify --slat and --slon.')

                if isnan(ds.lat[row]) or isnan(ds.lon[row]):
                    continue

                if abs(ds.lon[row]) > 180.0:
                    parser.error('invalid longitude input: all longitude values in ' + filein_sb.name + ' MUST be between -180/180E deg.')
                if abs(ds.lat[row]) > 90.0:
                    parser.error('invalid latitude input: all latitude values in ' + filein_sb.name + ' MUST be between -90/90N deg.')

                # construct sys call to sstval_extract
                sys_call_str = 'sstval_extract' + \
                               ' ifile=' + filein_sat.name + \
                               ' slon=' + str(ds.lon[row]) + \
                               ' slat=' + str(ds.lat[row]) + \
                               ' qual_check=qual_sst' + \
                               ' qual_check_distance=10' + \
                               ' global_att=1' + \
                               ' variable_att=1' + \
                               ' boxsize=' + str(dict_args['box_size'][0]) + flag_arg
                # variable_att flag needed to extract units
                # global_att flag needed to extract sensor/instrument names

                pid = subprocess.run(sys_call_str, shell=True, encoding='ascii', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                if pid.returncode == 99:
                    #if dict_args['verbose']:
                        #print('No matchup: in situ target not in granule.')
                    continue #no valid matchup
                elif pid.returncode == 101 or pid.returncode == 102:
                    parser.error('sstval_extract failed -- only accepts Level-2 (L2) satellite files. ' + \
                                    filein_sat.name + ' is not a valid L2 file')
                elif pid.returncode != 0:
                    print(pid.returncode)
                    print(pid.stdout)
                    print(pid.stderr)
                    parser.error('sstval_extract failed -- verify that the sstval_extract binary is compiled and on your PATH and that ' + \
                                    filein_sat.name + ' exists')

                # define structures to keep track of sstval_extract's output files
                file_ls = OrderedDict()
                file_del = []
                var_ls = []

                upix_ct = 0
                fpix_ct = dict_args['box_size'][0]^2
                pix_ct  = 0

                tims = 0
                tim_sat = 0

                cvs = []

                # parse the extract information
                file_del.append(filein_sat.name + '.qc');
                try:
                    fileobj = open(filein_sat.name + '.qc','r')
                    lines = fileobj.readlines()
                    for line in lines:
                        newline = re.sub("[\r\n]+",'',line)
                        if 'unflagged_pixel_count' in newline:
                            upix_ct = int(newline.split('=')[1])
                        elif 'flagged_pixel_count' in newline:
                            fpix_ct = int(newline.split('=')[1])
                        elif 'pixel_count' in newline:
                            pix_ct = int(newline.split('=')[1])
                        elif 'time' in newline:
                            try:
                                tims = re.search("(\d+)-(\d+)-(\d+)\s+(\d+):(\d+):(\d+)", newline.split('=')[1]);
                                tim_sat = datetime(year=int(tims.group(1)), \
                                          month=int(tims.group(2)), \
                                          day=int(tims.group(3)), \
                                          hour=int(tims.group(4)), \
                                          minute=int(tims.group(5)), \
                                          second=int(tims.group(6)))
                            except:
                                continue
                        elif 'variables' in newline:
                            var_ls = newline.split('=')[1].split(',')
                            for var in var_ls:
                                file_del.append(filein_sat.name + '.qc.' + var);
                                if 'l2_flags'   in var or \
                                   'cntl_pt_rows'   in var or \
                                   'clat'   in var or \
                                   'clon'   in var or \
                                   'day'   in var or \
                                   'elat'   in var or \
                                   'elon'   in var or \
                                   'msec'   in var or \
                                   'slat'   in var or \
                                   'slon'   in var or \
                                   'year'   in var:
                                    continue
                                file_ls[var] = filein_sat.name + '.qc.' + var
                    fileobj.close()
                except Exception as E:
                    print(E)
                    parser.error(' unable to open and read file ' + filein_sat.name + '.qc')

                # parse the satellite nc file information
                file_del.append(filein_sat.name + '.qc.global_attrs');
                [inst, plat] = readValEglobatt(filein_sat.name + '.qc.global_attrs', parser)
                if not inst:
                    inst = 'na'
                if not plat:
                    plat = 'na'

                # apply exclusion criteria
                # compute and evaluate the max time diff test
                if tim_sat > tim_max or tim_sat < tim_min:
                    clean_file_lis(file_del)
                    if dict_args['verbose']:
                        print('No matchup: failed MAX_TIME_DIFF, required =',dict_args["max_time_diff"][0],'Exclusion level = 1, Matrix row =',row)
                    continue #no valid matchup

                # compute and evaluate the min valid sat pix test
                if (pix_ct - fpix_ct) != 0:
                    if upix_ct >= dict_args['box_size'][0]:
                        pix_thresh = 100.0 * (upix_ct / (pix_ct - fpix_ct))
                        if pix_thresh < dict_args['min_valid_sat_pix'][0]:
                            clean_file_lis(file_del)
                            if dict_args['verbose']:
                                print('No matchup: failed MIN_VALID_SAT_PIX, required =',dict_args['min_valid_sat_pix'][0],'found =',pix_thresh,'Exclusion level = 4, Matrix row =',row)
                            continue #no valid matchup
                    else:
                        clean_file_lis(file_del)
                        if dict_args['verbose']:
                            print('No matchup: failed MIN_VALID_SAT_PIX, extracted satellite pixels less than box size, required =',dict_args['box_size'][0],'found =',upix_ct,'Exclusion level = 3, Matrix row =',row)
                        continue #no valid matchup
                else:
                    clean_file_lis(file_del)
                    if dict_args['verbose']:
                        print('No matchup: failed MIN_VALID_SAT_PIX, division by zero when deriving pix_thresh due to required L2FLAG criteria, Exclusion level = 2, Data row =',row)
                    continue #no valid matchup

                write_flag = 1 #only write out (write_flag == true), if matchups found

                #save L2_fname
                L2file_varname = inst + '_' + plat + '_l2fname'
                ds.addDataToOutput(row, L2file_varname.lower(), 'none', os.path.basename(filein_sat.name), True)

                #save tdiff
                tdiff_varname = inst + '_' + plat + '_tdiff'
                tdiff = tim_sat - dt
                ds.addDataToOutput(row, tdiff_varname.lower(), 'seconds', tdiff.total_seconds(), True)

                # save extract-variables
                for var in file_ls:

                    #save extracted values for each var in file_lis
                    [fmax, fmin, fmean, fmedian, fstdev, centerval, units] = readValEfile(file_ls[var], parser)

                    var_name = inst + '_' + plat + '_' + var.lower() + '_max'
                    ds.addDataToOutput(row, var_name.lower(), units, fmax, True)

                    var_name = inst + '_' + plat + '_' + var.lower() + '_min'
                    ds.addDataToOutput(row, var_name.lower(), units, fmin, True)

                    var_name = inst + '_' + plat + '_' + var.lower() + '_mean'
                    ds.addDataToOutput(row, var_name.lower(), units, fmean, True)

                    var_name = inst + '_' + plat + '_' + var.lower() + '_median'
                    ds.addDataToOutput(row, var_name.lower(), units, fmedian, True)

                    var_name = inst + '_' + plat + '_' + var.lower() + '_stddev'
                    ds.addDataToOutput(row, var_name.lower(), units, fstdev, True)

                    var_name = inst + '_' + plat + '_' + var.lower() + '_center_pixel_value'
                    ds.addDataToOutput(row, var_name.lower(), units, centerval, True)

                clean_file_lis(file_del)

            if write_flag == 1:
        
                for line in ds.comments:
                    if 'File ammended by OCSSW match-up maker script' in line:
                        comment_flag = True
                        break
            
                    else:
                        comment_flag = False
        
                if not dict_args['no_header_comment']:
                    ds.comments.append(' ')
                    ds.comments.append(' File ammended by OCSSW match-up maker script: mk_matchup_sst.py on ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '.')
                    ds.comments.append(' WARNINGS: This script does NOT adjust in situ data to water-leaving values')
                    ds.comments.append('           This script does NOT account for potential oversampling by the in situ data in time or space.')
                    ds.comments.append('           If successive calls to this script are made for a single in situ file AND multiple-valid-overpasses exist,')
                    ds.comments.append('               only the data from the last successive call will be saved to the output file. This may NOT be the best')
                    ds.comments.append('               quality satellite data in space and time.')
                    ds.comments.append(' EXCLUSION CRITERIA applied to this satellite file:')
                    ds.comments.append('     Box size of satellite extract = ' + str(dict_args['box_size'][0]) + ' pixels by ' + str(dict_args['box_size'][0]) + ' pixels')
                    ds.comments.append('     Minimum percent valid satellite pixels = ' + str(dict_args['min_valid_sat_pix'][0]))
                    ds.comments.append('     Maximum time difference between satellite and in situ = ' + str(dict_args['max_time_diff'][0]) + ' hours')
                    ds.comments.append('     The qual_sst value varies between 0 (best) and 4 (worst).')
                    ds.comments.append('     The qual_sst_mean (qual_sst_max) is the mean (max) of the ' + \
                                       str(dict_args['box_size'][0]) + ' by ' + str(dict_args['box_size'][0]) + ' pixel satellite extract.')
                    ds.comments.append(' ')

                print('Satellite/in situ match-up(s) found.')
                
                ds.writeSBfile(filein_sb.name)
        
            else:
                print('No valid satellite match-ups found.')
            
        print(' ')

    return


def clean_file_lis(file_ls):
    import os
    for d in file_ls:
        try:
            os.remove(d)
        except Exception as E:
            print(E)
            print('WARNING: Cleanup of ',d,' failed. Verify that you have read/write priviledges in the current working directory.')
    return


def readValEglobatt(fname, parser):
    import re
    inst = ''
    plat = ''
    try:
        fileobj = open(fname,'r')
        lines = fileobj.readlines()
        for line in lines:
            newline = re.sub("[\r\n]+",'',line)
            if 'instrument=' in newline:
                inst = newline.lower().split('=')[1]
            elif 'platform=' in newline:
                plat = newline.lower().split('=')[1]
        fileobj.close()
    except Exception as E:
        print(E)
        parser.error(' unable to open and read file ' + fname)
    return(inst, plat)


def readValEfile(fname, parser):

    import re

    missing   = ''
    fmax      = ''
    fmin      = ''
    fmean     = ''
    fmedian   = ''
    fstdev    = ''
    centerval = ''
    units     = ''

    try:

        fileobj = open(fname,'r')
        lines = fileobj.readlines()
        fileobj.close()

        for line in lines:

            newline = re.sub("[\r\n]+",'',line)

            if '_FillValue' in newline:
                missing = newline.split('=')[1]

        for line in lines:

            newline = re.sub("[\r\n]+",'',line)

            if 'filtered_max' in newline:
                fmax = newline.split('=')[1]
                if fmax == missing:
                    fmax = ''

            elif 'filtered_min' in newline:
                fmin = newline.split('=')[1]
                if fmin == missing:
                    fmin = ''

            elif 'filtered_mean' in newline:
                fmean = newline.split('=')[1]
                if fmean == missing:
                    fmean = ''

            elif 'filtered_median' in newline:
                fmedian = newline.split('=')[1]
                if fmedian == missing:
                    fmedian = ''

            elif 'filtered_stddev' in newline:
                fstdev = newline.split('=')[1]
                if fstdev == missing:
                    fstdev = ''

            elif 'center_value' in newline:
                centerval = newline.split('=')[1]
                if centerval == missing:
                    centerval = ''

            elif 'units' in newline:
                units = re.sub('\s', '_', newline.split('=')[1])

    except Exception as E:

        print(E)
        parser.error(' unable to open and read file ' + fname)

    return(fmax, fmin, fmean, fmedian, fstdev, centerval, units)


if __name__ == "__main__": main()
