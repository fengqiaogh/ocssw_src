# L2Extract

## Overview
L2Extract takes in a product or products if the user wants netCDF output from an L2 file and outputs a subset of the input, based on the parameters.

## Usage
```l2extract <arguments>```
```l2extract ifile spix epix sline eline pix_sub sc_sub ofile <prodlist>```
L2Extract's arguments are keyword=value pairs, and can either be specified on the command line or in a parameter file. Parameter files are typically defined using the file extension ```.par```. This can make for a cleaner-looking command line for multiple runs, and increase reproducability between runs. Individual arguments are separated by newlines. To employ a parameter file, specify the ```par``` argument on the command line. Note that if an argument is specified on the command line as well as in the parameter file, the command line argument will override.
The user must specify an input file ```ifile``` and output file ```ofile```, and four boundaries that specify the desired subset of the input file. These are ```spix```, ```epix```, ```sline```, and ```eline```.