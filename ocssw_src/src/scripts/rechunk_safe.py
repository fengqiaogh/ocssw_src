#!/usr/bin/env python3
from typing import Dict, List, Tuple,Union
import netCDF4 as nc
import sys
import argparse
import zipfile
import tarfile
import os
from shutil import move

PROGRAM_NAME = "rechunk_safe"
VERSION = "1.01"


def split_str(string: str, sep: str) -> Tuple:
    temp: List[int] = []
    for x in string.split(sep=sep):
        temp.append(int(x))
    return tuple(temp)


def move_files(extract_path: str):
    for root, _, files in os.walk(extract_path):
        for file in files:
            old_path = os.path.join(root, file)
            new_path = os.path.join(extract_path, os.path.basename(file))
            move(old_path, new_path)


def extract_file(archive: str, extract_path: str):
    print(f"Extracting {archive}")
    try:
        if archive.endswith(".zip") or archive.endswith(".ZIP"):
            if not os.path.exists(extract_path):
                os.makedirs(extract_path)
            with zipfile.ZipFile(archive, 'r') as zip_ref:
                zip_ref.extractall(extract_path)
            return 0
        elif archive.endswith(".tar"):
            tar_ext = "r"
        elif archive.endswith(".tar.gz") or archive.endswith(".tgz") or archive.endswith(".tar.Z"):
            tar_ext = 'r:gz'
        elif archive.endswith(".tar.bz2") or archive.endswith(".tbz"):
            tar_ext = 'r:bz2'
        else:
            raise (Exception("Invalid zip format"))
        try:
            with tarfile.open(archive, tar_ext) as tar:
                tar.extractall(extract_path)
            print(f"Successfully extracted '{archive}' to '{extract_path}'")
            return 0
        except FileNotFoundError:
            print(f"Error: File not found: '{archive}'")
        except tarfile.ReadError as e:
            print(f"Error: Could not open '{archive}' as a tar file: {e}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
    except Exception as e:
        print("Error while extracting %s (%s)" % (archive, sys.exc_info()[0]))
        print("Exception thrown: ",e)
    return 1


def select_files(dir: str, regex: str):
    outs = list()
    for root, _, files in os.walk(dir):
        for file in files:
            if regex in os.path.basename(file):
                inp_path = os.path.join(root, file)
                outs.append(inp_path)
    return outs


def converter(inp: Union[nc.Dataset, nc.Group], out_nc: Union[nc.Dataset, nc.Group], **kwargs):
    dims: Dict[str, nc.Dimension] = inp.dimensions
    for name, dim in dims.items():
        out_nc.createDimension(dim.name, dim.size)
    out_nc.setncatts(inp.__dict__)
    vars_nc: Dict[str, nc.Variable] = inp.variables
    for name, var in vars_nc.items():
        attr_var = var.ncattrs()
        chunking = tuple(var.chunking())
        if "chunk_size" in kwargs:
            if isinstance(kwargs["chunk_size"], dict):
                if name in kwargs["chunk_size"]:
                    chunking = kwargs["chunk_sizes"][name]
            else:
                chunking = kwargs["chunk_size"]
        dims: Tuple[nc.Dimension] = var.get_dims()
        dim_names = tuple([x.name for x in dims])
        var_out = out_nc.createVariable(varname=name, datatype=var.dtype, dimensions=dim_names, zlib=True,
                                        chunksizes=chunking, fill_value=var.getncattr("_FillValue"))
        for name_attr in attr_var:
            if name_attr != "_FillValue":
                var_out.setncattr(name_attr, var.getncattr(name_attr))
        var_out[:] = var[:]
        grps = inp.groups
        for grp in grps:
            inp_grp: nc.Group = inp[grp]
            out_grp: nc.Group = out_nc.createGroup(grp)
            converter(inp_grp, out_grp, **kwargs)


def full_converter(file: str, out_file: str, **kwargs):
    with nc.Dataset(file) as inp, nc.Dataset(out_file, mode="w") as out:
        inp: nc.Dataset = inp
        out: nc.Dataset = out
        converter(inp, out, **kwargs)


if __name__ == "__main__":
    print(PROGRAM_NAME, VERSION)
    parser = argparse.ArgumentParser(
        prog=PROGRAM_NAME,
        description='''\
                        This program unpacks an archive file (zip, tar), places all the archived files in a specified directory, and optionally rechunks variable with user specified chunk sizes 
                        ifile - one file or list of archive files, space separated. Must be supplied
                        odir - output directory. Default is "./"
                        chunk_size - chunk size. Default is none, optional
                        regex - regex for filename. Default is _radiance.nc , optional
                        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''\
                Supported ifile format: zip, tar.gz, tar.bzip2, tar 
                chunk_size must be specified as follows: "-c 123,568" or "-c var1:242,453 var2:111,565,878". In the former case the chunk size is applied uniformly to all variables 
                regex is a full file name or a substring of a full filename. Allows to select netCDF files for chunking. Default set to _radiance.nc for MERIS/OLCI
         ''')
    parser.add_argument("-v", "--version", action="version", version=VERSION)
    parser.add_argument(
        "--ifile", help="input archive file(s)", nargs='+', type=str, required=True)
    parser.add_argument(
        "--odir", help="output directory", type=str, required=False, default="./")
    parser.add_argument('-c', '--chunk_size', nargs='*', help='chunk sizes', required=False)
    parser.add_argument('-r', '--regex', type=str, help='regex for filename', required=False, default="_radiance.nc")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()
    odir = args.odir
    regex = args.regex
    ifiles = args.ifile
    for ifile in ifiles:
        # unzip
        status = extract_file(ifile, odir)
        if status != 0:
            print(f"Extraction failed for {ifile}")
            sys.exit(1)
        move_files(odir)
        chunk_to_pass = None
        if "chunk_size" in args and args.chunk_size:
            chunks: List[str] = args.chunk_size
            if len(chunks) == 1 and not (":" in chunks):
                chunk_to_pass = split_str(chunks[0], ',')
            else:
                chunk_to_pass = dict()
                for chunk in chunks:
                    out = chunk.split(sep=':')
                    if len(out) != 2:
                        print(f"The format for the chuck size {out} is not correct")
                        sys.exit(1)
                    var_name = out[0]
                    chunk_to_pass[var_name] = split_str(out[1], ',')
        if chunk_to_pass:
            input_files = select_files(odir, regex)
            for inp_file in input_files:
                out_temp = "temp.nc"
                print(f"Re-chunking {inp_file}")
                full_converter(inp_file, out_temp, chunk_size=chunk_to_pass)
                move(out_temp, inp_file)
    print("Done")
