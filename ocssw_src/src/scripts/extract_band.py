#!/usr/bin/env python3

"""

Program to extract a band (Variable) from one OBPG netCDF4 L2 file and create
a new file containing just that band.

The program is also meant to serve as an example of how OBPG files in netCDF4
format can be accessed and manipulated.

Helpful link:
    http://schubert.atmos.colostate.edu/~cslocum/netcdf_example.html

"""

import argparse
import os
import sys
import time
import traceback
import netCDF4 as nc

__version__ = '1.0.1_2016-10-18'

def copy_attribute(attr, in_group, out_group):
    """
    Copies the attribute named in attr from in_group to out_group.

    The checks for unicode data are needed because netCDF stores its data in
    Unicode. However, setncattr will prepend  "string" to the attributes if
    they are unicode values.  A long discussion of this appears at:
         https://github.com/Unidata/netcdf4-python/pull/389
    """
    attr_val = in_group.getncattr(attr)
    if isinstance(attr, unicode):
        attr = str(attr)
    if isinstance(attr_val, unicode):
        attr_val = str(attr_val)
    out_group.setncattr(attr, attr_val)

def copy_variable(src_var, src_grp, dest_grp):
    """
    Copies the netCDF4 Variable held in src_var from the src_grp Group to the
    dest_grp Group.
    """
    var_name = src_var[0]
    var_type = src_var[1].datatype
    new_var = dest_grp.createVariable(var_name, var_type,
                                      dimensions=src_var[1].dimensions)
    var_attrs = src_var[1].ncattrs()
    for var_attr in var_attrs:
        copy_attribute(var_attr, src_var[1], new_var)
    dest_grp.variables[var_name][:] =\
        src_grp.variables[var_name][:]

def create_subgroup(in_parent_grp, in_grp_name, out_parent_grp,
                    select_var=None):
    """
    Creates a Group in out_grp using cur_grp_name as the name of the created
    Group. Then recursively does the same for any subgroups of in_grp. If
    select_var is not None, then only that Variable will be copied into the new
    Group; otherwise, all Variables are copied.
    """
    #Create the new group
    new_grp = out_parent_grp.createGroup(in_grp_name)
    # Copy Variables
    grp_variables = in_parent_grp.groups[in_grp_name].variables
    if select_var:
        var_to_copy = grp_variables[select_var]
        copy_variable((select_var, var_to_copy),
                      in_parent_grp.groups[in_grp_name], new_grp)
    else:
        for var_item in grp_variables.items():
            if len(var_item) > 1:
                copy_variable(var_item, in_parent_grp.groups[in_grp_name],
                              new_grp)

    # Copy Attributes
    attrs = in_parent_grp[in_grp_name].ncattrs()
    for attr in attrs:
        copy_attribute(attr, in_parent_grp[in_grp_name], new_grp)
    # Copy subgroups
    subgroups = in_parent_grp[in_grp_name].groups
    if len(subgroups) > 0:
        for grp in subgroups:
            try:
                create_subgroup(in_parent_grp.groups[in_grp_name], grp, new_grp)
            except AttributeError:
                print (traceback.format_exc())

def get_cla():
    """
    Returns the command line arguments: band to extract, input filename, and
    output filename.
    """
    parser = argparse.ArgumentParser(description=\
        'Extracts a band from a netCDF4 OBPG file.')
    parser.add_argument('--version', action='version',
                        version='%(prog)s ' + __version__)
    parser.add_argument('band', type=str, help='band to be extracted')
    parser.add_argument('input_file', type=str, help='path to the input file')
    parser.add_argument('output_file', type=str, help='path to the output file')
    args = parser.parse_args()
    arg_list = [args.band, args.input_file, args.output_file]
    return arg_list

def main():
    """
    Primary driver of the program; get command line arguments, check the files
    specified and kick off the processing
    """
    args = get_cla()
    band = args[0]
    in_file = args[1]
    out_file = args[2]
    # './test/data/A2015086003500.L2_LAC_OC.nc'
    if os.path.exists(out_file):
        err_msg = 'Error! A file named {0} already exists!'.format(out_file)
        sys.exit(err_msg)

    if os.path.exists(in_file):
        with nc.Dataset(in_file) as in_dataset:
            if band in in_dataset.groups['geophysical_data'].variables:
                with nc.Dataset(out_file, 'w') as out_dataset:
                    process_main_dataset(in_dataset, band, out_dataset,
                                         ['date_created', 'history',
                                          'product_name'])
                    # Set attribute values that need to be customized.
                    cmd_line = ' '.join(sys.argv)
                    orig_history = in_dataset.getncattr('history')
                    if isinstance(orig_history, unicode):
                        orig_history = str(orig_history)
                    history_str = ' ; '.join([orig_history, cmd_line])
                    creation_date_str = time.strftime('%Y-%m-%dT%H:%M:%S.000Z',
                                                      time.gmtime())
                    out_dataset.setncattr('history', history_str)
                    out_dataset.setncattr('product_name', sys.argv[3])
                    out_dataset.setncattr('date_created', creation_date_str)
            else:
                err_msg = 'Error! Cannot locate band {0} in {1}.'.\
                    format(band, in_file)
                sys.exit(err_msg)
    else:
        err_msg = 'Could NOT find {0}.'.format(in_file)
        sys.exit(err_msg)

def process_main_dataset(in_dataset, band, out_dataset, attrs_to_skip=None):
    """
    Copy the top level dimensions, Variables, and Attributes from in_dataset
    to out_dataset. Copy Groups other than the 'geophysical_data' group. In the
    'geophysical_data' Group, call copy_band to copy only the Variable named
    in band.
    """
    for dimension in in_dataset.dimensions:
        out_dataset.createDimension(dimension,
                                    len(in_dataset.dimensions[dimension]))
    for var in in_dataset.variables:
        copy_variable(var, in_dataset, out_dataset)
    global_attrs = in_dataset.ncattrs()
    if attrs_to_skip:
        for attr in global_attrs:
            if not attr in attrs_to_skip:
                copy_attribute(attr, in_dataset, out_dataset)
    else:
        for attr in global_attrs:
            copy_attribute(attr, in_dataset, out_dataset)
    top_groups = in_dataset.groups
    for grp in top_groups:
        if grp == 'geophysical_data':
            create_subgroup(in_dataset, grp, out_dataset, band)
        else:
            create_subgroup(in_dataset, grp, out_dataset)

# The following allows the file to be imported without immediately executing.
if __name__ == '__main__':
    sys.exit(main())
