# L3Bin 5.14

## Overview

## Usage
```l3bin ifile=input-file ofile=output-file prod=prodlist```

or

```l3bin par=parfile.par```

The input file is a list of L3 binned files.
The argument-list is a set of keyword=value pairs. The arguments can
be specified on the commandline, or put into a parameter file, or the
two methods can be used together, with commandline over-riding.

  return value: 0=OK, 1=error, 110=no pixels binned. 

The list of valid keywords follows:

   ```help (boolean) (alias=h) (default=false) = print usage information```

   ```version (boolean) (default=false) = print the version information```

   ```verbose (boolean) (default=off) = Allow more verbose screen messages```

   ```dump_options (boolean) (default=false) = print information about each option```

   ```dump_options_paramfile (ofile) = print information about each option to paramfile```

   ```dump_options_xmlfile (ofile) = print information about each option to XML file```

   ```par (ifile) (alias=parfile) = input parameter file```

   ```pversion (string) (default=unspecified) = production version```

   ```ifile (ifile) (alias=in,infile) = input file name with list of L3 files```

   ```ofile (ofile) (alias=out) (default=output) = output bin file name```

   ```oformat (string) (default=netCDF4) = output file format```
          netCDF4: output a netCDF4 file
          hdf5:    output a HDF5 file
   ```merged (ofile) = merged file name```

   ```latnorth (float) (default=+90) = northern most latitude```

   ```latsouth (float) (default=-90) = southern most latitude```

   ```loneast (float) (default=+180) = eastern most longitude```

   ```lonwest (float) (default=-180) = western most longitude```

   ```sday (int) (default=1970001) = start datadate (YYYYDDD) ```

   ```eday (int) (default=2038018) = end datadate (YYYYDDD)```

   ```deflate (int) (default=5) = deflate level```

   ```orbit1 (int) (default=-1) = sorbit```

   ```orbit2 (int) (default=-1) = eorbit```

   ```median (int) (default=0) = median```

   ```noext (boolean) (default=off) = set to 1 to suppress generation of external files```

   ```unit_wgt (int) (default=0) = unit_wgt```

   ```composite_scheme (string) = composite scheme (min/max)```

   ```composite_prod (string) = composite product fieldname```

   ```reduce_fac (int) (default=1) = scale reduction factor (power of 2)```

   ```resolve (string) = bin resolution, overrides reduce_frac if defined```

          H: 0.5km
          Q: 250m
          HQ: 100m
          HH: 50m
          1: 1.1km
          2: 2.3km
          4: 4.6km
          9: 9.2km
          18: 18.5km
          36: 36km
          1D: 1 degree
          HD: 0.5 degree
          QD: 0.25 degree
   ```prod (string) (alias=out_parm) (default=DEFAULT) = bin products```
        
        [default=all products in L3 file]
   ```doi (string) = Digital Object Identifier (DOI) string```
