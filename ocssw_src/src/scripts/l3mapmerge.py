#! /usr/bin/env python3

'''
This L3 aggregator script will:
1. Take in L3 mapped files
2. Merge data
'''
import datetime
import argparse
from typing import List, Dict,Set
import sys
import netCDF4 as nc
import os

__version__ = "1.2"
verbose = False

l3_merge_attributes_exclude: Set[str] = {"l2_flag_names", "suggested_image_scaling_minimum",
                                         "suggested_image_scaling_maximum", "suggested_image_scaling_type","_lastModified",
                                         "suggested_image_scaling_applied", "data_minimum", "data_maximum", "data_bins",
                                         "id", "identifier_product_doi", "identifier_product_doi_authority", "keywords", 
                                         "processing_version","date_created"}
l3_merge_input_attrs = ["ifile","ofile","pversion","doi","product"]
l3_merge_processing_attrs = ["software_name","software_version", "l3map_files"]

def get_compression(var):
    # Set default compression parameters:
    comp_zlib = False
    deflate_level = None
    chunk_sizes = None
    # check for compression:
    if 'zlib' in var.filters():
    # Get the deflate level
        comp_zlib = True
        deflate_level = var.filters()['zlib']
        if var.chunking() != "contiguous":
            chunk_sizes = tuple(var.chunking())
    return comp_zlib,deflate_level,chunk_sizes

    
def get_variables(nc_handle):
    var_names = nc_handle.variables
    all_variables : List[str] = list()
    for var in var_names:
        var_nc: nc.Variable = nc_handle[var]
        # l3 var condition
        if len(var_nc.dimensions) == 2:
            if var_nc.dimensions[0] == "lat" and var_nc.dimensions[1] == "lon":
                all_variables.append(var)
    groups = nc_handle.groups
    for grp in groups:
        grp_nc: nc.Group = nc_handle[grp]
        all_variables += get_variables(grp_nc)
    return all_variables


def merge(idatasets,ofile, products : List[str],l3mapmerge_attrs : Dict[str,str] ):
        # delete output file if it exists
    try:
        os.remove(ofile)
    except OSError as e:
        pass
    pass
        # Create output file 
    orootgrp :  nc.Dataset = nc.Dataset(ofile, "w", format="NETCDF4")
    # Create dimensions for output file. Note all L2 files will/should have the same dimensions.
    for name, dim in idatasets[0].dimensions.items():
        orootgrp.createDimension(name, len(dim) if not dim.isunlimited() else None)

    # copy global attributes from first L3 file into mergefile 
    # create processing_control group, input_parameters group
    process_control_grp : nc.Group =  orootgrp.createGroup("processing_control")
    input_parameters_grp : nc.Group = process_control_grp.createGroup("input_parameters")

    for attr in idatasets[0].ncattrs():
        if attr in l3_merge_attributes_exclude:
            continue
        if attr == "product_name":
            orootgrp.setncattr(attr, ofile)
        else:
            orootgrp.setncattr(attr, idatasets[0].getncattr(attr))

    for attr_merge, val in l3mapmerge_attrs.items():
        if attr_merge in l3_merge_input_attrs:
            input_parameters_grp.setncattr(attr_merge, val)
            continue
        if attr_merge in l3_merge_processing_attrs:
            process_control_grp.setncattr(attr_merge, val)
            continue       
        orootgrp.setncattr(attr_merge, val)



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
                    # Set default compression parameters:
                    comp_zlib, deflate_level, chunk_sizes = get_compression(ivariable)
                    
                    ovariable = orootgrp.createVariable(name, ivariable.datatype, ivariable.dimensions, fill_value=fill_value, zlib = comp_zlib,complevel = deflate_level,  chunksizes = chunk_sizes)

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
    parser.add_argument('--product', type=str, help='product(s) to merge; comma separated', default = None)
    parser.add_argument('--doi', type=str, help='DOI', required=False)
    parser.add_argument('--pversion', type=str, help='processing version', required=False)
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
    input_files : List[str] = list()
    keywords : Set[str] = set()
    if len(ifiles) == 1: # supplied txt file 
        try:
            with open(args.ifile) as infile:
                for line in infile:
                    idatasets.append(nc.Dataset(line.strip(), 'r'))
                    input_files.append(line.strip())
        except IOError as e:
            print("--Error-- when trying to open input file: ", e)
            sys.exit(1)  
    else:
        try:
            for ncfile in ifiles:
                idatasets.append(nc.Dataset(ncfile.strip(), 'r'))
                input_files.append(ncfile.strip())     
        except IOError as e:
            print("--Error-- when trying to open input file: ", e)
            sys.exit(1)
    products : List[str]  = list()
    # needed to track  files the products come from
    found_products : Dict[str,str] = dict()
    # default looking for all l3 products
    if not args.product:
        for inp_file, idataset in zip(input_files,idatasets):
            file_prods = get_variables(idataset)
            if "keywords" in idataset.ncattrs():
                keywords.add(idataset.getncattr("keywords"))
            for fprod in file_prods:
                if fprod in found_products:
                    print(f"Warning: product {fprod} are found in {found_products[fprod]} and {inp_file}. Ignoring the variable from {inp_file}")
                else:
                    found_products[fprod] = inp_file
        for fprod in found_products:
            products.append(fprod)
    else: 
        products = args.product.split(",")
    products = ['lat', 'lon'] + products
    # l3merge attributes
    # history
    history = " ".join(sys.argv)
    l3mapmerge_attrs : Dict[str,str] = dict()
    l3mapmerge_attrs["history"] = history
    #input_sources
    l3map_files = ",".join(input_files)
    l3mapmerge_attrs["l3map_files"] = l3map_files
    #doi
    if args.doi:
        l3mapmerge_attrs["doi"] = args.doi
        l3mapmerge_attrs["identifier_product_doi"] = args.doi
        l3mapmerge_attrs["identifier_product_doi_authority"] = "http://dx.doi.org"
    else:
        l3mapmerge_attrs["doi"] = ""
    #pversion
    if args.pversion:
        l3mapmerge_attrs["pversion"] = args.pversion
        l3mapmerge_attrs["processing_version"] = args.pversion
    else:
        l3mapmerge_attrs["pversion"] = ""
    #ifile
    l3mapmerge_attrs["ifile"] = ifile
    #ofile
    l3mapmerge_attrs["ofile"] = args.ofile
    # product
    l3mapmerge_attrs["product"] = ",".join(found_products.keys())
    # date_created
    date_created = datetime.datetime.now()  
    l3mapmerge_attrs["date_created"] = date_created.strftime("%Y-%m-%dT%H:%M:%S.%fZ") 
    # software_name
    software_name = "l3mapmerge"
    l3mapmerge_attrs["software_name"] = software_name
    # software_version
    software_version = __version__
    l3mapmerge_attrs["software_version"] = software_version
    # keywords
    l3mapmerge_attrs["keywords"] = "; ".join(keywords)
    # keywords
    l3mapmerge_attrs["id"] = f"L3/{args.ofile}"
    print("Merging L3 map datasets...")
    merge(idatasets, args.ofile, products,l3mapmerge_attrs)