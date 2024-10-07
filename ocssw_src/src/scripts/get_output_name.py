#!/usr/bin/env python3

"""
Program to return the name of the next level file that would be created from
the input file name.
"""

import optparse
import os
import re
import sys
import traceback

import datetime
import mlp.get_obpg_file_type
import seadasutils.MetaUtils
import seadasutils.ProcUtils as ProcUtils
#import namer_constants
import mlp.name_finder_utils
import mlp.next_level_name_finder
import mlp.get_output_name_utils
import mlp.obpg_data_file
import mlp.viirs_next_level_name_finder

__version__ = '1.0.6-2021-08-01'
__author__ = 'byang'

#2345678901234567890123456789012345678901234567890123456789012345678901234567890

#########################################

FILE_TYPES_CONVERTER = {'Level 0': 'Level 0',
                        'Level 1A': 'Level 1A',
                        'l1agen' : 'Level 1A',
                        'Level 1B': 'Level 1B',
                        'Level 2': 'Level 2',
                        'L3bin': 'L3bin',
                        'Level 3 Binned': 'L3bin'}

#def get_1_file_name(input_name, file_typer, file_type, sensor,
#                    target_program, clopts):
# def convert_str_to_int(short_str):
#     """
#     Returns an integer taken from the passed in string.
#     """
#     try:
#         int_value = int(short_str)
#     except ValueError:
#         err_msg = "Error! Unable to convert {0} to integer.".format(short_str)
#         sys.exit(err_msg)
#     return int_value

# def find_extension(format_data_list, search_term):
#     """
#     Returns the extension from format_data_list that is indicated by
#     search_term.
#     """
#     extension = None
#     try:
#         # Are we searching for a matching index number ...
#         int(search_term)
#         tuple_index = 0
#     except ValueError:
#         # ... or a matching format name?
#         tuple_index = 1
#     # Use a generator to find the match.
#     format_index = next((i for i, t in enumerate(format_data_list) if format_data_list[i][tuple_index].lower() == search_term.lower()), None)
#     if (format_index != None) and (format_index < len(format_data_list)):
#         extension = format_data_list[format_index][2]
#     else:
#         for ext_candidate in format_data_list:
#             if search_term.lower() == ext_candidate[2].lower():
#                 extension = ext_candidate[2]
#                 break
#     return extension

# def get_1_file_name(data_file, target_program, clopts):
#     """
#     Return the next level name for a single file.
#     """
#     level_finder = mlp.name_finder_utils.get_level_finder([data_file],
#                                                       target_program,
#                                                       clopts)
#     next_level_name = level_finder.get_next_level_name()
#     return next_level_name

def get_multifile_next_level_name(data_files_list_info, target_program, clopts):
    """
    Return the next level name for a set of files.
    """
    for file_info in data_files_list_info:
        if file_info.file_type == 'unknown':
            err_msg = 'Error!  File {0} is of unknown type.'.format(
                file_info.name)
            sys.exit(err_msg)
    next_level_name = get_multifile_output_name(data_files_list_info,
                                                target_program, clopts)
    return next_level_name

def get_multifile_output_name(data_files_list_info, target_program, clopts):
    """
    Returns the file name derived from a group of files names.
    """
    list_file_type = data_files_list_info[0].file_type
    for data_file in data_files_list_info[1:]:
        if data_file.file_type != list_file_type:
            err_msg = 'Error!  File types do not match for {0} and {1}'.\
                      format(data_files_list_info[0].name, data_file.name)
            sys.exit(err_msg)
    output_name = mlp.get_output_name_utils.get_output_name(data_files_list_info, target_program, clopts)
    return output_name

# def read_fileformats():
#     """
#     Returns a tuple containing the file formats.
#     """

#     format_file_path = os.path.join(os.getenv('OCDATAROOT'), 'common',
#                                     'file_formats.txt')
#     if os.path.exists(format_file_path):
#         file_formats = []
#         format_file_hndl = open(format_file_path)
#         inp_lines = format_file_hndl.readlines()
#         format_file_hndl.close()
#         for line in inp_lines:
#             cleaned_line = line.strip()
#             if cleaned_line[0] != '#':
#                 #format = get_format(cleaned_line)
#                 file_format = tuple(cleaned_line.split(':'))

#                 file_formats.append(file_format)

#         return file_formats
#     else:
#         err_msg = 'Error! Cannot find file {0}.'.format(format_file_path)
#         sys.exit(err_msg)

#########################################

def get_command_line_data():
    """
    Returns the options and arguments from a command line call.
    """
    ver_msg = ' '.join(['%prog', __version__])
    use_msg = 'usage: %prog INPUT_FILE TARGET_PROGRAM'
    cl_parser = optparse.OptionParser(usage=use_msg, version=ver_msg)
    cl_parser.add_option('--oformat', dest='oformat', action='store',
                         type='string', help='output format')
    cl_parser.add_option('--odir', dest='odir', action='store',
                         type='string', help='output directory')
    cl_parser.add_option('--resolution', dest='resolution', action='store',
                         type='string',
                         help='resolution for l3mapgen')
    cl_parser.add_option('--suite', dest='suite', action='store',
                         type='string', help='data type suite')
    # cl_parser.add_option('--product', dest='product', action='store',
    #                      type='string', help='product type (for smigen)')
    (clopts, clargs) = cl_parser.parse_args()
    if len(clargs) == 0:
        print ("\nError! No input file or target program specified.\n")
        cl_parser.print_help()
        sys.exit(0)
    elif len(clargs) == 1:
        print ("\nError! No target program specified.\n")
        cl_parser.print_help()
        sys.exit(0)
    elif len(clargs) > 2:
        print ('\nError!  Too many arguments specified on the command line.')
        cl_parser.print_help()
        sys.exit(0)
    else:
        return clopts, clargs[0], clargs[1]

def get_data_files_info(file_list_file):
    """
    Returns a list of of data files.
    """
    file_info = []
    with open(file_list_file, 'rt') as file_list_file:
        inp_lines = file_list_file.readlines()
    for line in inp_lines:
        filename = line.strip()
        if os.path.exists(filename):
            file_typer = mlp.get_obpg_file_type.ObpgFileTyper(filename)
            file_type, sensor = file_typer.get_file_type()
            if file_type != 'unknown':
                stime, etime = file_typer.get_file_times()
                file_metadata = file_typer.attributes
                data_file = mlp.obpg_data_file.ObpgDataFile(filename, file_type,
                                                        sensor, stime, etime, file_metadata)
                file_info.append(data_file)
            else:
                err_msg = 'Error!  {0} is not an OBPG file.'.\
                          format(filename)
                sys.exit(err_msg)
        else:
            err_msg = 'Error! File {0} could not be found.'.format(filename)
            sys.exit(err_msg)
    file_info.sort(key=myfunc)
    return file_info

def myfunc(n):
    return n.start_time

def handle_unexpected_exception(exc_info):
    """
    Builds and prints an error message from the exception information,
    then exits.
    """
    exc_parts = exc_info
    err_type = str(exc_parts[0]).split('.')[1][0:-2]
    err_msg = 'Error!  Encountered {0}:'.format(str(err_type))
    print (err_msg)
    if DEBUG:
        traceback.print_exc()
    sys.exit(1)

def main():
    """
    Main function for when this module is called as a program.
    """
    ret_status = 0
    clopts, inp_name, targ_prog = get_command_line_data()

    #if not targ_prog in namer_constants.PROCESSABLE_PROGRAMS:
    if not targ_prog in PROCESSABLE_PROGRAMS:
        err_msg = 'Error!  The target program, "{0}", is not known.'.\
                  format(targ_prog)
        sys.exit(err_msg)
    if os.path.exists(inp_name):
        try:
            file_typer = mlp.get_obpg_file_type.ObpgFileTyper(inp_name)
            ftype, sensor = file_typer.get_file_type()
            if ftype == 'unknown':
                if seadasutils.MetaUtils.is_ascii_file(inp_name):
                    # Try treating the input file as a file list file.
                    data_files_info = get_data_files_info(inp_name)
                    if len(data_files_info) > 0:
                        next_level_name = get_multifile_next_level_name(
                            data_files_info, targ_prog, clopts)
                    else:
                        err_msg = "Error!  No OBPG files found in {0}".\
                                  format(inp_name)
                        sys.exit(err_msg)
                else:
                    # The input file wasn't a file list file.
                    err_msg = "File {0} is not an OBPG file.".format(inp_name)
                    sys.exit(err_msg)
            else:
                # The file is an OBPG file
                stime, etime = file_typer.get_file_times()
                file_metadata = file_typer.attributes
                data_file = mlp.obpg_data_file.ObpgDataFile(inp_name, ftype, sensor,
                                                        stime, etime,
                                                        file_metadata)
                next_level_name = mlp.get_output_name_utils.get_output_name([data_file], targ_prog, clopts)
            print ('Output Name: ' + next_level_name)
        except SystemExit as sys_ex:
            # The intention here is to catch exit exceptions we throw in other
            # parts of the program and continue with the exit, outputting
            # whatever error message was created for the exit.
            sys.exit(sys_ex)
        except:
            handle_unexpected_exception(sys.exc_info())
    else:
        err_msg = "Error!  File {0} was not found.".format(inp_name)
        sys.exit(err_msg)
    return ret_status

##########################################

#global DEBUG
DEBUG = False
DEBUG = True  # Comment out for production use
PROCESSABLE_PROGRAMS = \
    set(list(mlp.next_level_name_finder.NextLevelNameFinder.PROCESSING_LEVELS.keys()) +\
        list(mlp.next_level_name_finder.HawkeyeNextLevelNameFinder.PROCESSING_LEVELS.keys()) +\
        list(mlp.next_level_name_finder.ModisNextLevelNameFinder.PROCESSING_LEVELS.keys()) +\
        list(mlp.next_level_name_finder.SeawifsNextLevelNameFinder.PROCESSING_LEVELS.keys()) +\
        list(mlp.viirs_next_level_name_finder.ViirsNextLevelNameFinder.PROCESSING_LEVELS.keys()))

if __name__ == '__main__':
    sys.exit(main())
