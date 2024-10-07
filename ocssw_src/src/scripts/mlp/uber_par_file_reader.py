
"""
Module providing the ParReader class for reading parameter files for the
multilevel_processor.py program.
"""

import os
import re
import sys

__author__ = 'melliott'

SECTION_HEADER_TEXT = 'section'
VALID_KEYS = ['deletefiles', 'overwrite', 'use_existing', 'combine_files', 'use_ancillary', 'odir', 'ifile']

# class DuplicateEntry(Exception):
#     """
#     Exception class for duplicate entries in a section dictionary.
#     """
#     def __init__(self, value):
#         super(DuplicateEntry, self).__init__()
#         self.value = value

#     def __str__(self):
#         return repr(self.value)

# class InvalidEntry(Exception):
#     """
#     Exception class for invalid entries in a main section dictionary.
#     """
#     def __init__(self, value):
#         super(InvalidEntry, self).__init__()
#         self.value = value

#     def __str__(self):
#         return repr(self.value)

def add_par_entry(the_dict, the_key, val_to_add):
    """
    Adds an entry to the par file part of the passed in section dictionary.
    """
    if the_key in the_dict:
        the_dict[the_key].append(val_to_add)
    else:
        the_dict[the_key] = [val_to_add]

def add_sect_entry(filename, sect_key, the_dict, the_key, val_to_add):
    """
    Adds an entry to a section dictionary.
    """
    if sect_key != 'main':
        if the_key not in the_dict:
            the_dict[the_key] = val_to_add
        else:
            err_msg = 'Duplicate entry found for {0} in {1}'.format(the_key, filename)
            sys.exit(err_msg)
    else: 
        if is_key_valid(the_key):
            if the_key not in the_dict:
                the_dict[the_key] = val_to_add
            else:
                err_msg = 'Duplicate entry found for {0} in {1}'.format(the_key, filename)
                sys.exit(err_msg)
        else:
            err_msg = 'Invalid entry {0} found in {1}'.format(the_key, filename)
            sys.exit(err_msg)
        

def get_sect_key(line):
    """
    Returns the section name from a line of text.
    The line is expected to be of the form '# section SECTION_NAME' (without the quotes).
    """
    sect_key = ''
    sect_key = re.sub(r'^\s*\[\s*(.*?)\s*\]\s*$', '\\1', line)
#    sect_key = re.sub('\s*\[\s*', '', line)
#    sect_key = re.sub('\s*\]\s*', '', sect_key)
    return sect_key

def is_section_header(line):
    """
    Returns True if a line is the header for a new section; returns False otherwise.
    """
    section_pattern = r'\s*\[\s*\S+.*\]\s*'
    compiled_pattern = re.compile(section_pattern, re.IGNORECASE)
    if re.search(compiled_pattern, line.strip()):
        return True
    else:
        return False

def is_key_valid(key_str):
    """
    Returns True if opt is one of the valid options for Main section.
    """
    key_valid = False
    if key_str in VALID_KEYS:
        key_valid = True
    return key_valid

def is_whole_line_comment(line):
    """
    Returns True if an entire line is a comment; returns False otherwise.
    """
    if line.lstrip()[0:1] == '#':
        return True
    else:
        return False

class ParReader(object):
    """
    A class to perform reading of a seadas_processor.py parameter file.
    """
    def __init__(self, fname, acceptable_single_keys, section_converter=None):
        """
        Initializes the ParReader object: Confirms that the passed in file
        exists, is a regular file, and is readable.  If those are all true,
        sets the object filename value to the passed in value.
        """
        if os.path.exists(fname):
            if os.path.isfile(fname):
                if os.access(fname, os.R_OK):
                    self.filename = fname
                    self.acceptable_single_keys = acceptable_single_keys
                    self.section_converter = section_converter
                    self.par_results = {}
                else:
                    err_msg = "Error!  Unable to read {0}.".format(fname)
                    sys.exit(err_msg)
            else:
                err_msg = "Error!  {0} is not a regular file.".format(fname)
                sys.exit(err_msg)
        else:
            err_msg = "Error!  File {0} could not be found.".format(fname)
            sys.exit(err_msg)

    def _start_new_section(self, hdr_line):
        """
        Initializes a new section of the dictionary to be returned.
        """
        sect_dict = {}
        sect_key = get_sect_key(hdr_line)
        self.par_results[sect_key] = sect_dict
        if self.section_converter:
            if not (sect_key in list(self.section_converter.keys())):
                err_msg = 'Error! Section name "{0}" is not recognized.'.\
                          format(sect_key)
                sys.exit(err_msg)
        return sect_dict, sect_key

    def read_par_file(self):
        """
        Parses a parameter file, returning the contents in a dictionary of
        dictionaries.  The "outer" dictionary contains each section.  The "inner"
        dictionaries contain the parameter/value pairs.
        """
        sect_key = None
        sect_dict = {}
        with open(self.filename, 'rt') as par_file:
            for line in par_file.readlines():
                line = line.strip()
                if line != '' and not is_whole_line_comment(line):
                    if is_section_header(line):
                        sect_dict, sect_key = self._start_new_section(line)
                    else:
                        if sect_key != None:
                            if line.find('=') != -1:
                                key, val = line.split('=', 2)
                                if key == 'par':
                                    add_par_entry(sect_dict, 'par',
                                                   val.strip().strip('"').strip("'"))
                                else:
                                    
                                    add_sect_entry(self.filename, sect_key, sect_dict, key,
                                                       val.strip().strip('"').strip("'"))
                                        # if sect_key != 'main':
                                        #     add_sect_entry(sect_dict, key,
                                        #                val.strip().strip('"').strip("'"))
                                        # else:
                                        #     if is_key_valid(key):
                                        #        add_sect_entry(sect_dict, key,
                                        #                val.strip().strip('"').strip("'")) 
                                        #     else:
                                        #         err_msg = '{0} is not a valid key for main section'.format(key)
                                        #         sys.exit(err_msg)
                                    # except DuplicateEntry as dup_exc:
                                    #     err_msg = 'Duplicate entry found for {0} in {1}'.format(str(dup_exc), self.filename)
                                    #     sys.exit(err_msg)
                                    # except InvalidEntry as invalid_exc:
                                    #     err_msg = 'Invalid entry found for {0} in {1}'.format(str(invalid_exc), self.filename)
                                    #     sys.exit(err_msg)
                            elif line.strip in self.acceptable_single_keys:
                                sect_dict[key] = 'True'
                            else:
                                if line.startswith('--'):
                                    add_sect_entry(self.filename, sect_key, sect_dict, line, None)
                                else:
                                    err_msg = 'Found entry "{0}" with no value in {1}'.format(line, self.filename)
                                    sys.exit(err_msg)
                        else:
                            err_msg = 'Error in {0}, no section header found!'.\
                                      format(self.filename)
                            sys.exit(err_msg)
            if sect_key != '':
                self.par_results[sect_key] = sect_dict
        return self.par_results
