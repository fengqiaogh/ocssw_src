# L1AGen_OCI

## Overview
L1agen_oci takes in Level 0 files (typically ending in .oci), and creates a Level 1A file out of them. Examples of Level 1 files can be found [here.](https://oceancolor.gsfc.nasa.gov/cgi/browse.pl?sen=amod) The foremost dependency of l1agen_oci is NetCDF, which is used for its capabilities in storing scientific data in a compact and easy-to-use format. The purpose of going from Level 0 to Level 1 in the context of l1agen_oci is to take the raw data sent from the instrument and reformat it into something more usable. 

## Usage
```
l1agen_oci OCI_packet_file granule_len

    [-g | --maxgap maximum missing scans allowed in an output file]
    [-k | --hktlist SC_HKT_input_list]
    [-s | --swir_loff_set SWIR_LOFF_config, list of 9 comma separated integers]
    [-t | --start_time YYYYmmddTHHMMSS or YYYY-mm-ddTHH:MM:SS]
    [-o | --outlist output_list_file]
    [-f | --outfile output_file_name]
    [-d | --doi doi_string]
    [-v | --pversion processing_version]
    [-n | --noSPW]

Return   1 Fatal error
       110 No ancillary packets found in file
       120 No L1A file was generated
```
For ease of use, the user may specify a parameter file, in which the parameters to L1AGen_OCI are specified. Parameter files are typically defined using the file extension ```.par```. This can make for a cleaner-looking command line for multiple runs, and increase reproducability between runs. Individual arguments are separated by newlines. To employ a parameter file, specify the ```par``` argument on the command line. Note that if an argument is specified on the command line as well as in the parameter file, the command line argument will override.
