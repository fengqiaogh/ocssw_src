"""
Utility functions for get_output_name.
"""

import calendar
import datetime
import os
import re
import sys
import types
#from lxml.html.diff import start_tag

import mlp.get_obpg_file_type
import mlp.obpg_data_file as obpg_data_file
import seadasutils.ProcUtils as ProcUtils
import seadasutils.time_utils as time_utils

__author__ = 'byang'

__version__ = '1.0.6-2021-09-29'

DEBUG = False
#DEBUG = True

def convert_str_to_int(short_str):
    """
    Returns an integer taken from the passed in string.
    """
    try:
        int_value = int(short_str)
    except ValueError:
        err_msg = "Error! Unable to convert {0} to integer.".format(short_str)
        sys.exit(err_msg)
    return int_value

def find_extension(format_data_list, search_term):
    """
    Returns the extension from format_data_list that is indicated by
    search_term.
    """
    extension = None
    try:
        # Are we searching for a matching index number ...
        int(search_term)
        tuple_index = 0
    except ValueError:
        # ... or a matching format name?
        tuple_index = 1
    # Use a generator to find the match.
    format_index = next((i for i, t in enumerate(format_data_list) if format_data_list[i][tuple_index].lower() == search_term.lower()), None)
    if (format_index != None) and (format_index < len(format_data_list)):
        extension = format_data_list[format_index][2]
    else:
        for ext_candidate in format_data_list:
            if search_term.lower() == ext_candidate[2].lower():
                extension = ext_candidate[2]
                break
    return extension

def get_base_element(data_files, target_program, clopts):
    """
    Returns the base element from input filename.
    """
    indicator = get_platform_indicator(data_files[0])
    if len(data_files) > 1:
        for file in data_files:
            indctr = get_platform_indicator(file)
            if indctr.find(indicator) == -1:
                if indctr.find('MODIS') != -1 and indicator.find('MODIS') != -1:
                    indicator = 'CROSS_MODIS.'
                elif indctr.find('VIIRS') != -1 and indicator.find('VIIRS') != -1:
                     indicator = 'CROSS_VIIRS.'
                else:
                    indicator = 'CROSS_SENSOR.'
                    break
    if target_program.find('bin') != -1 or \
                        target_program == 'mapgen' or target_program == 'l3mapgen':
        base_element_name = indicator + get_l3_time(data_files)
    else:
        if data_files[0].file_type.find('Level 0') != -1:
            time_stamp = get_l0_timestamp(data_files[0].name) 
        else:
            time_stamp = data_files[0].start_time
        dt_obj = datetime.datetime.strptime(time_stamp, "%Y%j%H%M%S")  
    
        base_element_name = "{}{}T{}".format(indicator, dt_obj.strftime("%Y%m%d"), time_stamp[7:])
    
    return base_element_name

def get_end_doy_year(data_files):
    """
    Extract a day of year and year from an L0 file's metadata and return
    them as integer values .
    """
    if data_files[-1].end_time:
        year = convert_str_to_int(data_files[-1].end_time[0:4])
        day = convert_str_to_int(data_files[-1].end_time[4:7])
    elif data_files[-1].metadata:
        day_str = 'End Day'
        yr_str = 'End Year'
        day = convert_str_to_int(data_files[-1].metadata[day_str])
        year = convert_str_to_int(data_files[-1].metadata[yr_str])
    else:
        err_msg = 'Error! Cannot find end time for {0}'.format(
                data_files[-1].name)
        sys.exit(err_msg)
    return day, year

def _get_data_files_info(flf):
    """
    Returns a list of data files read from the specified input file.
    """
    data_file_list = []
    with open(flf, 'rt') as file_list_file:
        inp_lines = file_list_file.readlines()
    for line in inp_lines:
        filename = line.strip()
        if os.path.exists(filename):
            file_typer = mlp.get_obpg_file_type.ObpgFileTyper(filename)
            file_type, sensor = file_typer.get_file_type()
            stime, etime = file_typer.get_file_times()
            data_file = obpg_data_file.ObpgDataFile(filename, file_type,
                                                    sensor, stime, etime)
            data_file_list.append(data_file)
    data_file_list.sort()
    return data_file_list

def get_l0_timestamp(l0_file_name):
    """
    A method to get the date/time stamp from L0 files.
    """
    # Todo: Add check & handling for time stamp in metadata.
    if os.path.exists(l0_file_name + '.const'):
        with open(l0_file_name + '.const') as constructor_file:
            constructor_data = constructor_file.readlines()
        for line in constructor_data:
            if line.find('starttime=') != -1:
                start_time = line[line.find('=') + 1].strip()
                break
        time_stamp = ProcUtils.date_convert(start_time, 't', 'j')
    else:
        input_basename = os.path.basename(l0_file_name)
        matched_name = re.match(r"MOD00.?.[AP](\d\d\d\d\d\d\d).(\d\d\d\d)", input_basename)
        if matched_name is None:
            matched_name = re.match(r"[AP](\d\d\d\d\d\d\d)(\d\d\d\d)\d\d\.L0_.{3}", input_basename)
        if matched_name:
            time_stamp = matched_name.group(1)+matched_name.group(2) + '00'
        else:
            err_msg = "Unable to determine time stamp for input file {0}".\
            format(l0_file_name)
            sys.exit(err_msg)  
    return time_stamp 

def get_days_diff(day1, day2):
    """
    Returns the number of days between two days, by subtracting day2 from day1.
    """
    return (day1 - day2).days

def get_end_day_year(metadata):
    """
    Returns the end day and year for a file, determined from the contents of
    metadata as ints.
    """
    if 'End Day' in metadata:
        eday = convert_str_to_int(metadata['End Day'])
    elif 'Period End Day' in metadata:
        eday = convert_str_to_int(metadata['Period End Day'])
    elif 'time_coverage_end' in metadata:
        eday = time_utils.convert_month_day_to_doy(
            metadata['time_coverage_end'][5:7],
            metadata['time_coverage_end'][8:10],
            metadata['time_coverage_end'][0:4])
    else:
        err_msg = 'Error! Cannot determine end day.'
        sys.exit(err_msg)
    if 'End Year' in metadata:
        eyear = convert_str_to_int(metadata['End Year'])
    elif 'Period End Year' in metadata:
        eyear = convert_str_to_int(metadata['Period End Year'])
    elif 'time_coverage_end' in metadata:
        eyear = convert_str_to_int(metadata['time_coverage_end'][0:4])
    else:
        err_msg = 'Error! Cannot determine end year.'
        sys.exit(err_msg)
    return eday, eyear

def get_extension(program, clopts):
    """
    Returns the extension appropriate for the program.
    """
    extension_dict = {'level 1a': '.nc',
                'modis_L1A': '.hdf',
                'geo': '.nc',
                'modis_GEO': '.hdf',
                'geolocate_hawkeye': '.nc',
                'geolocate_viirs': '.nc',
                'l1aextract': '.nc',
                'l1aextract_modis': '.hdf',
                'l1aextract_viirs': '.nc',
                'l1aextract_seawifs': '.hdf',
                'l1brsgen': '.hdf',
                'l1mapgen': '.png',
                'level 1b': '.nc',
                'modis_L1B': '.hdf',
                'calibrate_viirs': '.nc',
                'l1bgen': '.nc',
                'l2gen': '.nc',
                'l2extract': '.nc',
                'l2brsgen': '.hdf',
                # 'l2mapgen': '.ppm',
                'l2bin': '.nc',
                'l3bin': '.nc',
                'l3mapgen': '.nc',
                'mapgen': '.png'}
    extension_allowed_dict = {'level 1a': '.nc',
                'modis_L1A': '.hdf',
                'geo': {'.hdf','.nc'},
                'modis_GEO': '.hdf',
                'geolocate_hawkeye': '.nc',
                'geolocate_viirs': '.nc',
                'l1aextract': '.nc',
                'l1aextract_modis': '.hdf',
                'l1aextract_viirs': '.nc',
                'l1aextract_seawifs': '.hdf',
                'l1brsgen': {'.hdf', '.bin', '.png', '.ppm'},
                'l1mapgen': {'.ppm', '.png', '.tiff'},
                'level 1b': {'.hdf','nc'},
                'modis_L1B': '.hdf',
                'calibrate_viirs': '.nc',
                'l1bgen': {'.nc', '.hdf'},
                'l2gen': {'.nc', '.hdf'},
                'l2extract': {'.nc', '.hdf'},
                'l2brsgen': {'.hdf', '.png', '.ppm'},
                # 'l2mapgen': {'ppm', 'png', 'tiff'},
                'l2bin': '.nc',
                'l3bin': '.nc',
                'l3mapgen': {'.nc', '.ppm', '.png', '.tiff'},
                'mapgen': {'.nc', '.ppm', '.png', '.tiff'}}
    if program in list(extension_dict.keys()):
        if clopts and 'oformat' in clopts and clopts['oformat'] != None :
            file_formats = read_fileformats()
            format_ext = '.' + find_extension(file_formats, clopts['oformat'])
            if format_ext == '.':
                ext = '.hdf'
            elif format_ext in extension_allowed_dict[program]:
                ext = format_ext
            else: 
                err_msg = 'Error! The oformat {0} is not supported by {1}.'.format(
                    clopts['oformat'], program)
                sys.exit(err_msg)
        else:
            ext = extension_dict[program]
    return ext

def get_extra_bits(data_files, target_program, clopts):
    """
    A method to get the extra bits for l2bin, l3mapgen.
    """
    extra_bits =''
    if target_program.find('bin') != -1 or \
                        target_program == 'l3mapgen' or target_program == 'mapgen':
        sday, syear = get_start_doy_year(data_files)
        eday, eyear = get_end_doy_year(data_files)
        if sday and syear and sday > 0 and syear > 0:
            sdate = datetime.datetime.strptime(str(syear) + '-' + str(sday),
                                                '%Y-%j')
        else:
            err_msg = 'Error! Cannot process start date data: year = ' \
                '{0}, doy = {1}'.format(syear, sday)
            sys.exit(err_msg)
        if eday and eyear and eday > 0 and eyear > 0:
            edate = datetime.datetime.strptime(str(eyear) + '-' + str(eday),
                                                '%Y-%j')
        else:
            err_msg = 'Error! Cannot process end date data: year = {0},' \
                'doy = {1}'.format(eyear, eday)
            sys.exit(err_msg)
        days_diff = get_days_diff(edate, sdate)
        if clopts and 'suite' in clopts and clopts['suite']!= None:
            suite = '.' + clopts['suite']
        # elif data_files[0].metadata != None and 'suite' in data_files[0].metadata and \
        #                                 data_files[0].metadata['suite'].strip() != '':
        #     suite = '.' + data_files[0].metadata['suite'].strip()
        else:
            suite = ''
        if suite == None:
            suite = ''
        if days_diff == 0:
            extra_bits = '.DAY' + suite
        else:
            if days_diff == 7:
                extra_bits = '.8D' + suite
            else:
                extra_bits = '.CU' + suite
        if (target_program.find('l3mapgen') != -1 or target_program.find('mapgen') != -1)\
                and clopts and 'resolution' in clopts and clopts['resolution'] != None:
            extra_bits += '.' + clopts['resolution']
    elif target_program.find('l2gen') != -1:
        if clopts and 'suite' in clopts and clopts['suite'] != None:
            extra_bits = '.' + clopts['suite']   
    if data_files[0].name.find('sub') != -1:  
        extra_bits += '.sub'
    return extra_bits

def get_l3_time(data_files):
        """
        An internal method to return the L3bin time from an L2 or
        L3bin file.
        """
        l3_time = ''
        if len(data_files) == 1 and (re.search("\.\d\d\d\d\d\d\d\dT\d\d\d\d\d\d\.", data_files[0].name) or 
                                    re.search("\d\d\d\d\d\d\d\d\d\d\d\d\d\.L2", data_files[0].name)):
            time_stamp = data_files[0].start_time
            dt_obj = datetime.datetime.strptime(time_stamp, "%Y%j%H%M%S")
            l3_time = "{}T{}".format(dt_obj.strftime("%Y%m%d"), time_stamp[7:]) 
        else:     
            sday, syear = get_start_doy_year(data_files)
            eday, eyear = get_end_doy_year(data_files)
            if sday and syear and sday > 0 and syear > 0:
                sdate = datetime.datetime.strptime(str(syear) + '-' + str(sday),
                                                '%Y-%j')
            else:
                err_msg = 'Error! Cannot process start date data: year = {0}' \
                ', doy = {1}'.format(syear, sday)
                sys.exit(err_msg)
            if eday and eyear and eday > 0 and eyear > 0:
                edate = datetime.datetime.strptime(str(eyear) + '-' + str(eday),
                                                '%Y-%j')
            else:
                err_msg = 'Error! Cannot process end date data: year = {0},' \
                'doy = {1}'.format(eyear, eday)
                sys.exit(err_msg)
            days_diff = (edate, sdate)
            if days_diff == 0:
                l3_time = '%d%02d%02d' % (syear, sdate.month, sdate.day)
            else:
                l3_time = '%d%02d%02d%d%02d%02d' % (syear, sdate.month, sdate.day, eyear, edate.month, edate.day)
        return l3_time

def get_level(program, data_files):
    """
    Returns the level element for the target_program.
    """
    level_dict = {'level 1a': '.L1A',
                'modis_L1A': '.L1A',
                'geo':  '.GEO', 
                'modis_GEO': '.GEO',
                'geolocate_hawkeye': '.GEO',
                'geolocate_viirs': '.GEO',
                'l1aextract': '.L1A.sub',
                'l1aextract_modis': '.L1A.sub',
                'l1aextract_viirs': '.L1A.sub',
                'l1aextract_seawifs': '.L1A.sub',
                'l1brsgen': '.L1BRS',
                'level 1b': '.L1B',
                'modis_L1B': '.L1B',
                'calibrate_viirs': '.L1B',
                'l1bgen': '.L1B',
                'l2gen': '.L2',
                'l2extract': '.L2.sub',
                'l2brsgen': '.L2BRS',
                # 'l2mapgen': '.L2',
                'l2bin': '.L3b',
                'l3bin': '.L3b',
                'l3mapgen': '.L3m',
                'mapgen': '.L3m'}
    if program == 'geo' and data_files[0].sensor.find('VIIRS') != -1:
        program = 'geolocate_viirs'
    if program in list(level_dict.keys()):
        level = level_dict[program]
    elif program == 'l1mapgen':
        if data_files[0].file_type.find('Level 1A') != -1:
            level = 'L1A_MAP'
        elif data_files[0].file_type.find('Level 1B') != -1:
            level = '.L1B_MAP'
    return level

def get_output_name(data_files, target_program, clopts):
    """
    Returns the file name derived from the input file name, target program name and oformat .
    """
    if clopts and not isinstance(clopts, dict):
        # Assuming the clopts passed in is a group of options from optparse.
        clopts = vars(clopts)
    if target_program == 'mapgen':
        if data_files[0].file_type.find('Level 1') != -1 and len(data_files) ==1:
            target_program = 'l1mapgen'
    output_name = get_base_element(data_files, target_program, clopts) + get_level(target_program, data_files)\
                + get_extra_bits(data_files, target_program, clopts) + get_extension(target_program, clopts)
    return output_name

def get_platform_indicator(data_file):
        """
        Returns a character which indicates what platform (instrument) the
        data in the file is from.  
        """
        indicator_dict = {'Aquarius': 'SACD', 
                          'CZCS': 'NIMBUS7', 
                          'GOCI': 'COMS',
                          'HICO': 'ISS',
                          'MERIS': 'ENVISAT',
                          'MOS': 'IRSP3', 
                          'HAWKEYE': 'SEAHAWK1',
                          'OCI': 'PACE',
                          'OCIS': 'PACE',
                          'OCM2': 'OCEANSAT2', 
                          'OCTS': 'ADEOS',
                          'OLCI S3A': 'S3A',
                          'OLCI S3B': 'S3B',
                        #   'OLI L8': 'LANDSAT8', 
                        #   'OLI L9': 'LANDSAT9', 
                          'OSMI': 'KOMSAT1',
                          'SeaWiFS': 'SEASTAR',
                          'SGLI': 'GC1'}
        data_type = ''
        if data_file.name.find('CROSS_SENSOR') != -1:
            indicator = 'CROSS_SENSOR.'
            return indicator
        elif data_file.name.find('CROSS_MODIS') != -1:
            indicator = 'CROSS_MODIS.'
            return indicator
        elif data_file.name.find('CROSS_VIIRS') != -1:
            indicator = 'CROSS_VIIRS.'
            return indicator
        if data_file.sensor in list(indicator_dict.keys()):
            sensor = data_file.sensor.upper()
            indicator = indicator_dict[data_file.sensor]
            if sensor.find('OCTS') != -1:
                data_type = 'GAC'
            elif sensor.find('MERIS') != -1:
                if 'FRS' in data_file.name:
                    data_type = 'FRS'
                elif 'RR' in data_file.name:
                    data_type = 'RR'
            elif sensor.find('SEAWIFS') != -1:
                if 'GAC' in data_file.name:
                    data_type = 'GAC'
                elif 'infile' in data_file.metadata and 'GAC' in data_file.metadata['infile']:
                    data_type = "GAC"
                else:
                    data_type = "LAC"   
            elif sensor.find('OLCI') != -1:
                sensor = 'OLCI'
                if 'data_type' in data_file.metadata:
                    data_type = data_file.metadata['data_type']
                elif 'EFR' in data_file.name:
                    data_type = 'EFR'
                elif 'ERR' in data_file.name:
                    data_type = 'ERR'
        elif data_file.sensor.find('MODIS') != -1:
            sensor = 'MODIS'
            if data_file.sensor.find('Aqua') != -1:
                indicator = 'AQUA'
            elif data_file.sensor.find('Terra') != -1:
                indicator = 'TERRA'
            else:
                err_msg = 'Error!  Could not determine platform indicator for MODIS file {0}.'.\
                          format(data_file.name)
                sys.exit(err_msg)
        elif data_file.sensor.find('VIIRS') != -1:
            sensor = 'VIIRS'
            if data_file.sensor.find('J1') != -1:
                indicator = 'JPSS1'
            elif data_file.sensor.find('J2')!= -1:
                indicator = 'JPSS2'
            elif data_file.sensor.find('NPP')!= -1:
                indicator = 'SNPP'
            else:
                err_msg = 'Error!  Could not determine platform indicator for VIIRS file {0}.'.\
                          format(data_file.name)
                sys.exit(err_msg)
        elif data_file.sensor.find('MSI') != -1:
            sensor = 'MSI'
            if data_file.sensor == 'MSI S2A':
                indicator = 'S2A'
            elif data_file.sensor == 'MSI S2B':
                indicator = 'S2B'
            else:
                err_msg = 'Error!  Could not determine platform indicator for MSI file {0}.'.\
                          format(data_file.name)
                sys.exit(err_msg)
        elif data_file.sensor.find('OLI') != -1:
            sensor = 'OLI'
            if data_file.sensor == 'OLI L8':
                indicator = 'LANDSAT8'
            elif data_file.sensor == 'OLI L9':
                indicator = 'LANDSAT9'
            else:
                err_msg = 'Error!  Could not determine platform indicator for MSI file {0}.'.\
                          format(data_file.name)
                sys.exit(err_msg)
        else:
            err_msg = 'Error!  Platform indicator, {0}, for {1} is not known.'.\
                      format(data_file.sensor, data_file.name)
            sys.exit(err_msg)
        # for dfile in self.data_files[1:]:
        #     if dfile.sensor in list(indicator_dict.keys()):
        #         if indicator != indicator_dict[dfile.sensor]:
        #             indicator = 'X'
        #             break
        #     else:
        #         indicator = 'X'
        #         break
        if data_type:
            indicator += '_' + sensor + '_' + data_type + '.'
        else:
            indicator += '_' + sensor + '.'
        return indicator

def get_start_doy_year(data_files):
        """
        Extract a day of year and year from a file's metadata and return
        them as integer values .
        """
        if data_files[0].end_time:
            year = convert_str_to_int(data_files[0].start_time[0:4])
            day = convert_str_to_int(data_files[0].start_time[4:7])
        elif data_files[0].metadata:
            day_str = 'Start Day'
            yr_str = 'Start Year'
            day = convert_str_to_int(data_files[0].metadata[day_str])
            year = convert_str_to_int(data_files[0].metadata[yr_str])
        else:
            err_msg = 'Error! Cannot find end time for {0}'.format(
                data_files[0].name)
            sys.exit(err_msg)
        return day, year

def get_start_day_year(metadata):
    """
    Returns the start day and year for a file, determined from the contents of
    metadata as ints.
    """
    if 'Start Day' in metadata:
        sday = convert_str_to_int(metadata['Start Day'])
    elif 'Period Start Day' in metadata:
        sday = convert_str_to_int(metadata['Period Start Day'])
    elif 'time_coverage_start' in metadata:
        sday = time_utils.convert_month_day_to_doy(
            metadata['time_coverage_start'][5:7],
            metadata['time_coverage_start'][8:10],
            metadata['time_coverage_start'][0:4])
    else:
        err_msg = 'Error! Cannot determine start day.'
        sys.exit(err_msg)
    if 'Start Year' in metadata:
        syear = convert_str_to_int(metadata['Start Year'])
    elif 'Period Start Year' in metadata:
        syear = convert_str_to_int(metadata['Period Start Year'])
    elif 'time_coverage_start' in metadata:
        syear = convert_str_to_int(metadata['time_coverage_start'][0:4])
    else:
        err_msg = 'Error! Cannot determine start year.'
        sys.exit(err_msg)
    return sday, syear

def get_time_period_extension(start_date_str, end_date_str):
    """
    Return the part of the file extension based on the time period within the
    start and end dates.
    """
    first_date = datetime.datetime.strptime(start_date_str, '%Y%j%H%M%S')
    last_date = datetime.datetime.strptime(end_date_str, '%Y%j%H%M%S')
    date_diff = last_date - first_date
    if date_diff.days == 0:
        time_ext = '.DAY'
    elif date_diff.days == 7:
        time_ext = '.8D'
    elif is_month(first_date, last_date):
        time_ext = '.MO'
    elif is_year(first_date, last_date):
        time_ext = '.YR'
    else:
        time_ext = '.CU'
    return time_ext

def is_month(day1, day2):
    """
    Returns True if the days are the endpoints of a month; False otherwise.
    """
    return day1.month == day2.month and day1.day == 1 and\
           day2.day == calendar.monthrange(day1.year, day1.month)[1]

def is_year(day1, day2):
    """
    Returns True if the days are the endpoints of a year; False otherwise.
    """
    return day1.year == day2.year and day1.month == 1 and day1.day == 1 and\
           day2.month == 12 and day2.day == 31

def read_fileformats():
    """
    Returns a tuple containing the file formats.
    """

    format_file_path = os.path.join(os.getenv('OCDATAROOT'), 'common',
                                    'file_formats.txt')
    if os.path.exists(format_file_path):
        file_formats = []
        format_file_hndl = open(format_file_path)
        inp_lines = format_file_hndl.readlines()
        format_file_hndl.close()
        for line in inp_lines:
            cleaned_line = line.strip()
            if cleaned_line[0] != '#':
                #format = get_format(cleaned_line)
                file_format = tuple(cleaned_line.split(':'))

                file_formats.append(file_format)

        return file_formats
    else:
        err_msg = 'Error! Cannot find file {0}.'.format(format_file_path)
        sys.exit(err_msg)