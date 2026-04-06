"""
Module containing utilities to manipulate netCDF4 files.
"""

__author__ = "gfireman"

import re
import sys
import time
from os.path import basename

import netCDF4
import numpy as np


def sanitize(name):
    name = re.sub(r"/+", r"#", name)
    name = re.sub(r"^\.+", r"", name)
    name = re.sub(r"^\%+", r"Percent", name)
    return name


def nccopy_var(srcvar, dstgrp, indices=None, verbose=False, debug=False):
    """
    Copy a netCDF4 variable, optionally subsetting some dimensions.

    Function to copy a single netCDF4 variable and associated attributes.
    Optionally subset specified dimensions.

    Parameters
    ----------
    srcvar : netCDF4.Variable
        Open variable to be copied
    dstgrp : netCDF4.Group
        Open Group or Dataset destination object to copy stuff to
    indices : dict, optional
        Dict of dimname:[indexarr] to subset a dimension
    verbose : boolean, optional
        Print extra info

    Side Effects
    ------------
    Strings are written as H5T_CSET_ASCII, not H5T_CSET_UTF8
    Empty attributes are written as scalar "" instead of NULL
    """

    grpname = srcvar.group().path
    if grpname == "/":
        grpname = ""
    varname = sanitize(srcvar.name)
    varpath = "/".join([grpname, varname])
    # if verbose:
    #    print(f"var: {varpath}")

    # get variable definition
    zlib = srcvar.filters().get("zlib", False)
    shuffle = srcvar.filters().get("shuffle", False)
    complevel = srcvar.filters().get("complevel", 0)

    # ensure dimensions are correctly handled
    outputDims = {}
    if indices:
        for dimname in srcvar.dimensions:
            if dimname in indices:
                outputDims[dimname] = len(indices[dimname])

    # define chunk sizes
    if srcvar.chunking() == "contiguous":  # srcvar is not chunked
        dstChunks = None
    else:
        srcChunks = list(srcvar.chunking())
        dstChunks = []
        for idx, dimname in enumerate(srcvar.dimensions):
            if dimname in outputDims:
                dstChunks.append(min(srcChunks[idx], outputDims[dimname]))
            else:
                dstChunks.append(srcChunks[idx])

    # define compound datatypes as needed
    if srcvar._iscompound:
        name = sanitize(srcvar._cmptype.name)
        if debug:
            print(f"compound type: {name}")
        dtype = dstgrp.createCompoundType(srcvar.dtype, name)
        fillvalue = None  # no fill value for compound types
    else:  # _isenum, _isprimitive, _isvlen
        dtype = srcvar.dtype
        fillvalue = srcvar._FillValue if "_FillValue" in srcvar.ncattrs() else None

    # create variable with same storage format
    dstvar = dstgrp.createVariable(
        sanitize(srcvar.name),
        dtype,
        srcvar.dimensions,
        fill_value=fillvalue,
        zlib=zlib,
        shuffle=shuffle,
        complevel=complevel,
        chunksizes=dstChunks,
    )

    # set variable attributes
    for name, value in srcvar.__dict__.items():
        if name == "_FillValue":
            continue
        name = sanitize(name)
        if debug:
            print(f"\t{name} = {value}")
        dstvar.setncattr(name, value)

    # if no dimension changes, copy all
    if not indices or not any(k in indices for k in srcvar.dimensions):
        if verbose:
            print(f"\tcopying {srcvar.name}")
        if srcvar._isvlen:
            dstvar[0] = srcvar[0]
        else:
            dstvar[:] = srcvar[:]

    # otherwise, copy only the subset
    else:
        if verbose:
            print(f"\tsubsetting {srcvar.name}")
        tmpvar = srcvar[:]
        for dimname in indices:
            try:
                axis = srcvar.dimensions.index(dimname)
            except ValueError:
                continue
            tmpvar = np.take(tmpvar, indices[dimname], axis=axis)
        dstvar[:] = tmpvar

    # make sure it's written out
    dstgrp.sync()


def nccopy_grp(
    srcgrp,
    dstgrp,
    indices=None,
    verbose=False,
    debug=False,
    copyvars=True,
    copygrps=True,
):
    """
    Recursively copy a netCDF4 group, optionally subsetting some dimensions.

    Function to recursively copy a netCDF4 group,
    with associated attributes, dimensions and variables.
    Optionally subset specified dimensions.

    Parameters
    ----------
    srcgrp : netCDF4.Group
        Open Group or Dataset source object containing stuff to be copied
    dstgrp : netCDF4.Group
        Open Group or Dataset destination object to copy stuff to
    indices : dict, optional
        Dict of dimname:[indexarr] to subset a dimension
    verbose : boolean, optional
        Print extra info
    """

    if verbose:
        print(f"grp: {srcgrp.path}")

    # copy all group attributes
    for name, value in srcgrp.__dict__.items():
        name = sanitize(name)
        if debug:
            print(f"att: {name} = {value}")
        dstgrp.setncattr(name, value)

    # define each dimension
    for dimname, dim in srcgrp.dimensions.items():
        if dim.isunlimited():
            dimsize = None
        elif indices and dimname in indices:
            dimsize = len(indices[dimname])
        else:
            dimsize = len(dim)
        if debug:
            print(f"dim: {dimname} = {dimsize}")
        dstgrp.createDimension(dimname, dimsize)

    # define each variable
    if copyvars:
        for varname, srcvar in srcgrp.variables.items():
            nccopy_var(srcvar, dstgrp, indices=indices, verbose=verbose, debug=debug)

    # define each subgroup
    if copygrps:
        for grpname, srcsubgrp in srcgrp.groups.items():
            dstsubgrp = dstgrp.createGroup(grpname)
            nccopy_grp(
                srcsubgrp, dstsubgrp, indices=indices, verbose=verbose, debug=debug
            )


def nccopy(srcfile, dstfile, verbose=False, debug=False):
    """
    Copy a netCDF4 file.

    Function to copy a netCDF4 file to a new file.
    Intended mostly as a demonstration.

    Parameters
    ----------
    srcfile : str
        Path to source file; must be netCDF4 format.
    dstfile : str
        Path to destination file; directory must exist.
    verbose : boolean, optional
        Print extra info
    """

    with netCDF4.Dataset(srcfile, "r") as src, netCDF4.Dataset(dstfile, "w") as dst:
        if verbose:
            print(f"src: {src.filepath()}")
            print(f"dst: {dst.filepath()}")

        # copy values without masking
        src.set_auto_mask(False)
        dst.set_auto_mask(False)
        nccopy_grp(src, dst, verbose=verbose, debug=debug)


def ncsubset_vars(srcfile, dstfile, subset, verbose=False, debug=False, **kwargs):
    """
    Copy a netCDF4 file, with some dimensions subsetted.

    Parameters
    ----------
    srcfile : str
        Path to source file; must be netCDF4 format.
    dstfile : str
        Path to destination file; directory must exist.
    subset : dict, optional
        Dict of dimname:[startindex,endindex] to subset a dimension
    verbose : boolean, optional
        Print extra info

    Side Effects
    ------------
    Strings are written as H5T_CSET_ASCII, not H5T_CSET_UTF8
    Empty attributes are written as scalar "" instead of NULL
    """

    if verbose:
        print(f"opening {srcfile}")
    with netCDF4.Dataset(srcfile, "r") as src:
        src.set_auto_mask(False)

        # validate input
        for dimname in subset:
            if subset[dimname][0] > subset[dimname][1]:
                print(f'Invalid indices for dimension "{dimname}"; exiting.')
                return
        for dimname, dim in src.dimensions.items():
            if (dimname in subset) and any(
                (0 > d or d > len(dim) - 1) for d in subset[dimname]
            ):
                oldsubset = subset.copy()
                subset[dimname] = np.clip(
                    subset[dimname], a_min=0, a_max=len(dim) - 1
                ).tolist()
                print(
                    f'Clipping "{dimname}" dimension indices to match input file: {oldsubset[dimname]} -> {subset[dimname]}'
                )

        # construct index arrays
        indices = {k: np.arange(subset[k][0], subset[k][1] + 1) for k in subset}

        # copy source file
        if verbose:
            print(f"opening {dstfile}")
        with netCDF4.Dataset(dstfile, "w") as dst:
            dst.set_auto_mask(False)
            nccopy_grp(src, dst, indices=indices, verbose=verbose, debug=debug)
            update_history(dst, **kwargs)
            # dstfile closes automatically

        # srcfile closes automatically


def update_history(dataset, timestamp=None, cmdline=None):
    """
    Update 'date_created' and 'history' attributes

    Function to add or update 'date_created' and 'history'
    attributes for specified dataset (usually root).

    Parameters
    ----------
    dataset : netCDF4.Group
        Open Group or Dataset destination object to update
    timestamp : time.struct_time, optional
        Timestamp to add to history attribute
        Defaults to current time
    cmdline : string, optional
        Description to add to history attribute
    """

    if not timestamp:
        timestamp = time.gmtime()
    fmt = "%Y-%m-%dT%H:%M:%SZ"  # ISO 8601 extended date format
    date_created = time.strftime(fmt, timestamp)

    if not cmdline:
        cmdline = " ".join([basename(sys.argv[0])] + sys.argv[1:])
    cmdline = "".join([date_created, ": ", cmdline])
    if "history" in dataset.ncattrs():
        history = "".join([dataset.history.strip(), "; ", cmdline])
    else:
        history = cmdline

    dataset.setncattr("date_created", date_created)
    dataset.setncattr("history", history)
