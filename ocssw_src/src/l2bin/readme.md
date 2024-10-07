# L2Bin

## Overview
L2Bin takes in a Level 2 file and outputs a binned version of the same file.

## Usage
```l2bin <arguments>```

L2Bin's arguments are keyword=value pairs, and can either be specified on the command line or in a parameter file. Parameter files are typically defined using the file extension ```.par```. This can make for a cleaner-looking command line for multiple runs, and increase reproducability between runs. Individual arguments are separated by newlines. To employ a parameter file, specify the ```par``` argument on the command line. Note that if an argument is specified on the command line as well as in the parameter file, the command line argument will override.