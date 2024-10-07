#! /usr/bin/env python3

'''
This L3 aggregator script will:
1. Take in L3 mapped files
2. Merge data
'''

import argparse

import sys
import netCDF4 as nc
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
    


def merge(idatasets,ofile, products):
        # delete output file if it exists
    try:
        os.remove(ofile)
    except OSError as e:
        pass
    pass
        # Create output file 
    orootgrp = nc.Dataset(ofile, "w", format="NETCDF4")
    # Create dimensions for output file. Note all L2 files will/should have the same dimensions.
    for name, dim in idatasets[0].dimensions.items():
        orootgrp.createDimension(name, len(dim) if not dim.isunlimited() else None)

    # copy global attributes from first L3 file into mergefile 
    for attr in idatasets[0].ncattrs():
        if attr == "product_name":
            orootgrp.setncattr(attr, ofile)
        else:
            orootgrp.setncattr(attr, idatasets[0].getncattr(attr))

    # recursively copy first dataset's contents into orootgrp
    copyGroups(orootgrp, idatasets[0])
    dims = orootgrp.dimensions
    nlat : int = dims["lat"].size
    nlon : int = dims["lon"].size
    if verbose:
        print("copying variables from the file roots:")
    for ds in idatasets:
        dims = ds.dimensions
        if dims["lat"].size != nlat:
            print(f"--Error--: latitude dimenstions mismatch for {ds.product_name}")
            exit(1)
        if dims["lon"].size != nlon:
            print(f"--Error--: longitude dimenstions mismatch for {ds.product_name}")
            exit(1)
        for name, ivariable in ds.variables.items(): 
            if name in products: 
                if (name not in orootgrp.variables):
                    if verbose:
                        print("   ", name)

                    if "_FillValue" in ivariable.ncattrs():
                        fill_value = ivariable.getncattr("_FillValue")
                    else:
                        fill_value = None
                    ovariable = orootgrp.createVariable(name, ivariable.datatype, ivariable.dimensions, fill_value=fill_value)

                    for attr in ivariable.ncattrs():
                        #create a new attribute name/value pair 
                        if attr != "_FillValue":
                            ovariable.setncattr(attr, ivariable.getncattr(attr)) 

                    ovariable[:] = ivariable[:]
    
    orootgrp.close()


if __name__ == "__main__":
        # parse command line
    parser = argparse.ArgumentParser(
        description='Merge Level 3 mapped files')
    parser.add_argument('--version', action='store_true', help='print program version')
    parser.add_argument('-v', '--verbose', help='print status messages', action='store_true')
    parser.add_argument('ifile', help='Input files (comma separated list of filenames) or a ASCII file which contains paths of L3 netcdf4 files to be merged')
    parser.add_argument('ofile', help='Output netcdf file which will contain merged L3 file data')
    parser.add_argument('--product', type=str, help='product(s) to map; comma separated', default = None)

    print("l3mapmerge", __version__)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    if(sys.argv[1] == "--version"):
        sys.exit(0)

    args = parser.parse_args()
    verbose = args.verbose
    ifile : str = args.ifile
    ifiles = ifile.split(sep=',')
    # Read and store each line in the .txt infile as a netcdf dataset
    idatasets = []
    if len(ifiles) == 1: # supplied txt file
        try:
            with open(args.ifile) as infile:
                for line in infile:
                    idatasets.append(nc.Dataset(line.strip(), 'r'))
        except IOError as e:
            print("--Error-- when trying to open input file: ", e)
            sys.exit(1)  
    else:
        try:
            for ncfile in ifiles:
                idatasets.append(nc.Dataset(ncfile.strip(), 'r'))     
        except IOError as e:
            print("--Error-- when trying to open input file: ", e)
            sys.exit(1)
    if not args.product:
        print("No products were specified. Exiting ...")
        exit(1) 
    products : list[str] = args.product.split(",") + ['lat', 'lon']
     
    print("Merging datasets...")
    merge(idatasets, args.ofile, products)