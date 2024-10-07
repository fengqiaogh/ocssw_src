# l3mapgen

## Overview
This program takes a product (or products if netCDF output) from an L3 bin or SMI file, reprojects the data using proj.4 and writes a mapped file in the requested output format.  

## Usage
```l3mapgen <arguments>```

The argument list is a set of keyword=value pairs.  Arguments can be specified on the command line, or put into a parameter (.par) file, or the two methods can be used together, with command line taking precedence. Parameter files are typically defined using the file extension ```.par```, and contain arguments (parameters) for the call to l3mapgen. This can make for a cleaner-looking command line for multiple runs, and increase reproducability between runs. Individual arguments are separated by newlines. The content of a par file (ex. foo.par) might look ike this: 
```
ifile=OCI_FOO_20220321T191111_L3bin.nc
ofile=OCI_FOO_20220321T191111_L3map.nc
```
With a corresponding execution of ```l1mapgen par=foo.par```
Return values:
    0 = No error
    1 = Error
    110 = No valid data to map
