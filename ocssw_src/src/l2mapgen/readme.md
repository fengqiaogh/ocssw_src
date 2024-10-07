# l2mapgen 1.0

## Overview
This program takes a product from an L2 file, maps it using a Plate Carree cylindrical projection, and produces a gray scale PGM or color PPM file.

## Usage 
```l2mapgen <argument-list>```

The argument list is a set of keyword=value pairs.  Arguments can be specified on the command line, or put into a parameter (.par) file, or the two methods can be used together, with command line taking precedence. Parameter files are typically defined using the file extension ```.par```, and contain arguments (parameters) for the call to l3mapgen. This can make for a cleaner-looking command line for multiple runs, and increase reproducability between runs. 

The list of valid keywords follows:
   
   ```help (boolean) (alias=h) (default=false) = print usage information```
   
   ```version (boolean) (default=false) = print the version information```
   
   ```dump_options (boolean) (default=false) = print information about each option```
   
   ```dump_options_paramfile (ofile) = print information about each option to paramfile```
   
   ```dump_options_xmlfile (ofile) = print information about each option to XML file```
   
   ```par (ifile) = input parameter file```
   
   ```ifile (ifile) = input L2 file name or file with a list of files names```
   
   ```ofile (ofile) = output map filename (NULL=STDOUT)```
   
   ```prod (string) = product name```
   
   ```apply_pal (boolean) (default=false) = apply color palette, false = grayscale```
   
   ```palfile (ifile) (default=default) = palette filename```
   
   ```palette_dir (ifile) (default=$OCDATAROOT/common/palette) = palette directory```
   
   ```product_table (ifile) (default=$OCDATAROOT/common/smigen_product_table.dat) = product table```
   
   ```flaguse (string) = flags to be masked```
   
   ```quality (int) (default=2) = minimum allowable quality level for SST.  Valid only for SST  and only if ```qual_sst or qual_sst4 SDS exist
   
   ```mask (boolean) (default=no) = apply mask to land, cloud and glint (see below)```
   
   
   ```datamin (float) (default=0.0) = minimum value for data scaling```
        (default see SMI product table)
   
   ```datamax (float) (default=0.0) = maximum value for data scaling```
        (default see SMI product table)
   
   ```stype (int) (default=0) = scaling type (default see SMI product table)```
        1: LINEAR
        2: LOG
   
   ```east (float) (default=0.0) = Map East longitude```
        (default=scene(s) Easternmost Longitude)
   
   ```west (float) (default=0.0) = Map West longitude```
        (default=scene(s) Westernmost Longitude)
   
   ```north (float) (default=0.0) = Map North latitude```
        (default=scene(s) Northernmost Longitude)
   
   ```south (float) (default=0.0) = Map South latitude```
        (default=scene(s) Southernmost Longitude)
   
   ```width (int) (default=800) = width of the output image```
   
   ```threshold (float) (default=5) = minimum percentage of the area of interest that must receive valid pixel data before an image is generated```
   
   ```outmode (string) (default=ppm) = format of the output file```
        ppm: PPM or PGM image file (alias 1)
        png: PNG color or grayscale image file (alias 2)
        tiff: TIFF color or grayscale geo tiff image file (alias 3)

   
   If the "mask" option is set, the output PGM image will be masked for
   flags defined in the flaguse parameter. The "no data" pixel value will
   change from 0 to 255, and pixel values 252, 253, and 254 will represent the
   sunglint, land, and all other (e.g. clouds/ice,hilt,atmfail,navfail,chlfail)
   masks, respectively. NOTE: sunglint is NOT masked by default, but if it is
   added to the flaguse parameter, it will be distinguished in the masking as
   medium gray.  If a palette is applied and the mask option is set, the
   palette values will be modified:
                  Value   R       G       B
                  252     128     128     128
                  253     160     82      45
                  254     255     255     255
                  255     0       0       0

   By default, this program sends its results to standard output as a
   PGM-formatted binary data stream.  Save it to a file via ">" or pipe it
   to your favorite image display program.  The output image is rendered in
   a Plate Carree projection.
