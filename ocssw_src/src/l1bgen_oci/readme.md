# L1BGen_OCI
Generates l1b files for OCI using l1a files

### Usage 
```l1bgen_oci <arguments>```

Available arguments:
```version (boolean) (default=false) = print the version information
dump_options (boolean) (default=false) = print information about each option
dump_options_paramfile (ofile) = print information about each option to paramfile
dump_options_xmlfile (ofile) = print information about each option to XML file
par (ifile) = input parameter file
ifile (ifile) = Input L1A file
ofile (ofile) = Output L1B file
cal_lut (ifile) = CAL LUT file
geo_lut (ifile) = GEO LUT file
doi (string) = Digital Object Identifier (DOI) string
pversion (string) (default=Unspecified) = processing version string
demfile (string) (default=$OCDATAROOT/common/gebco_ocssw_v2020.nc) = Digital elevation map file
radiance (boolean) (default=false) = Generate radiances
``` 

Some notable arguments are ```ifile```, ```ofile```, and ```par```. These specify the input file, the output file, and the parameter file, respectively.

Parameter files are typically defined using the file extension ```.par```, and contain arguments (parameters) for the call to l1bgen_oci. This can make for a cleaner-looking command line for multiple runs, and increase reproducability between runs. Individual arguments are separated by newlines. The content of a par file (ex. foo.par) might look ike this: 
```
ifile=OCI_FOO_20220321T191111_L1A.nc
ofile=OCI_FOO_20220321T191111_L1B.nc
```
With a corresponding execution of ```l1bgen_oci par=foo.par```
To employ a parameter file, specify the ```par``` argument on the command line. Note that if an argument is specified on the command line as well as in the parameter file, the command line argument will override.