#! /usr/bin/env python3

'''
This L2 aggregator script will:
1. Take in L2 files
2. Merge data
'''

import argparse
import pathlib
import sys
import netCDF4
import os

__version__ = "1.0"
verbose = False


# Copy groups from topB to topA

def copyGroups(topA, topB):
    for group_name, group in topB.groups.items():
        topA.createGroup(group_name)

        # Copy group attributes
        for attr in group.ncattrs():
            topA[group_name].setncattr(attr, group.getncattr(attr))
            
        # bail before copying variables from the geophysical_data group
        if group_name == "geophysical_data":
            continue

        if verbose:
            print("copying variables for group", group_name, ":")

        # Iterate through variables and copy all variables and their attributes
        for name, variable in group.variables.items():
            if verbose:
                print("   ", name)

            if "_FillValue" in variable.ncattrs():
                fill_value = variable.getncattr("_FillValue")
            else:
                fill_value = None

            topA[group_name].createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value)

            for attr in variable.ncattrs():
                if attr != "_FillValue":
                    topA[group_name].variables[name].setncattr(attr, variable.getncattr(attr)) 

            topA[group_name].variables[name][:] = variable[:]
  
        if len(group.groups) != 0:
            copyGroups(topA[group_name], group)
    

def merge(idatasets, mergefile, products):
    
    # delete output file if it exists
    try:
        os.remove(mergefile)
    except OSError as e:
        pass

    # Create output file 
    orootgrp = netCDF4.Dataset(pathlib.Path(mergefile), "w", format="NETCDF4")

    # Create dimensions for output file. Note all L2 files will/should have the same dimensions.
    for name, dim in idatasets[0].dimensions.items():
        orootgrp.createDimension(name, len(dim) if not dim.isunlimited() else None)

    # copy global attributes from first L2 file into mergefile 
    for attr in idatasets[0].ncattrs():
        if attr == "product_name":
            orootgrp.setncattr(attr, mergefile)
        else:
            orootgrp.setncattr(attr, idatasets[0].getncattr(attr))

    # recursively copy first dataset's contents into orootgrp except for geophysical_data
    copyGroups(orootgrp, idatasets[0])

    if verbose:
        print("copying variables for group geophysical_data:")

    # merge geophysical_data 
    # combining geophysical data from each dataset in idatasets
    geophysical_data = orootgrp["/geophysical_data"]
    for ds in idatasets: 
        for name, ivariable in ds['/geophysical_data'].variables.items(): 
            if name in products: 
                if (name not in geophysical_data.variables):
                    if verbose:
                        print("   ", name)

                    if "_FillValue" in ivariable.ncattrs():
                        fill_value = ivariable.getncattr("_FillValue")
                    else:
                        fill_value = None
                    ovariable = geophysical_data.createVariable(name, ivariable.datatype, ivariable.dimensions, fill_value=fill_value)

                    for attr in ivariable.ncattrs():
                        #create a new attribute name/value pair 
                        if attr != "_FillValue":
                            ovariable.setncattr(attr, ivariable.getncattr(attr)) 

                    ovariable[:] = ivariable[:]
    
    orootgrp.close()


if __name__ == "__main__":

    # parse command line
    parser = argparse.ArgumentParser(
        description='Merge Level 2 files')
    parser.add_argument('--version', action='store_true', help='print program version')
    parser.add_argument('-v', '--verbose', help='print status messages', action='store_true')
    parser.add_argument('ifile', help='Input text file specifying L2 files to be merged')
    parser.add_argument('ofile', help='Output netcdf file which will contain merged L2 file data')
    parser.add_argument('--product', type=str, help='product(s) to map; comma separated', default = "None")
    
    print("l2merge", __version__)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    if(sys.argv[1] == "--version"):
        sys.exit(0)
        
    args = parser.parse_args()
    verbose = args.verbose

    # Read and store each line in the .txt infile as a netcdf dataset
    idatasets = []
    try:
        with open(args.ifile) as infile:
            for line in infile:
                idatasets.append(netCDF4.Dataset(line.strip(), 'r'))
    except IOError as e:
        print("Error when trying to open input file: ", e)
    infile.close()

    products = args.product.split(",")
    products.append("l2_flags")

    print("Merging datasets...")
    try: 
        merge(idatasets, args.ofile, products)
        print("Merging was successful")
    except Exception as e:
        print("Unable to merge L2 files: ", e)