# l2gen 

## Overview
L2gen creates Level 2 files out of any Level 1 file (A, B, or C). Examples of Level 2 files can be found [here.](https://oceancolor.gsfc.nasa.gov/cgi/browse.pl?sen=amod)
A major dependency of l2gen is NetCDF, which is used for its capabilities in storing scientific data in a compact and easy-to-use format.

## Usage
l2gen argument-list

The argument-list is a set of keyword=value pairs. The arguments can be specified on the command line, or put into a parameter file, or the two methods can be used together, with command line taking precedence. Parameter files are defined with a file extension of .par, and individual arguments are separated by newlines, with variadic options being defined on the same line. An example of a .par file, foo.par:
```
ifile=OCI_EFLA_V11_20220321T183042_L1B.nc
ofile=OCI_EFLA_V11_20220321T183042_L2.nc
l2prod=rhos ndvi evi evi2 
proc_land=1
```
Note that these arguments are not exhaustive. Run ``` l2gen -h ``` for a comprehensive listing of options and how to use them.
For ease of use, the user may specify a parameter file, in which the parameters to L2Gen are specified. Parameter files are typically defined using the file extension ```.par```. This can make for a cleaner-looking command line for multiple runs, and increase reproducability between runs. Individual arguments are separated by newlines. To employ a parameter file, specify the ```par``` argument on the command line. Note that if an argument is specified on the command line as well as in the parameter file, the command line argument will override.
