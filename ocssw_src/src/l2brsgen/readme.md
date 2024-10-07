# l2brsgen 2.0

## Overview
  This program takes a product from a L2 file, subsamples the file and writes a browse file

## Usage: 
```l2brsgen <arguments>```


The argument list is a set of keyword=value pairs.  Arguments can be specified on the command line, or put into a parameter (.par) file, or the two methods can be used together, with command line taking precedence. Parameter files are typically defined using the file extension ```.par```, and contain arguments (parameters) for the call to l3mapgen. This can make for a cleaner-looking command line for multiple runs, and increase reproducability between runs. 

The list of valid keywords follows:   
   ```help (boolean) (alias=h) (default=false) = print usage information```

   ```version (boolean) (default=false)      = print version information```
   
   ```dump_options (boolean) (default=false) = print information about each option```
   
   ```dump_options_paramfile (ofile)         = print information about each option to paramfile```
   
   ```dump_options_xmlfile (ofile)           = print information about each option to XML file```
   
   ```par       (ifile) = input parameter file```
   
   ```ifile     (ifile) = input L2 file name```
   
   ```ofile     (ofile) (default=output)    = output filename```
   
   ```prod      (string) (default=chlor_a)  = product name```
   
   ```quality   (int) (default=999)         = highest quality value acceptable```
   
   ```rflag     (string) (default=ORIGINAL) = replacement flag```
   
   ```flaguse   (string)                    = Flags used to mask data```
   
   ```chl_flags (string) (default=ATMFAIL,HILT,STRAYLIGHT,CLDICE,LOWLW,CHLWARN,CHLFAIL,NAVWARN,MAXAERITER,NAVFAIL,FILTER,HIGLINT) = Flags used to mask data for chl product if flaguse not set```
   
   ```sst_flags (string) (default=SSTFAIL) = Flags used to mask data for sst product if flaguse not set```
   
   ```spixl (int) (default=1) = start pixel number```
   
   ```epixl (int) (default=-1) = end pixel number (-1=the last pixel)```
   
   ```dpixl (int) (default=1) = pixel subsampling interval```
   
   ```sline (int) (default=1) = start line number```
   
   ```eline (int) (default=-1) = end line number (-1=the last line)```
   
   ```dline (int) (default=1) = line subsampling interval```
   
   ```apply_pal (boolean) (default=no) = apply color palette, false = grayscale```
   
   ```palfile (ifile) (default=default) = palette filename.  "default" means the```
        
        palette is chosen using the product table file
   
   ```palette_dir (ifile) (default=$OCDATAROOT/common/palette) = directory```
        
        containing the palette files
   
   ```product_table (ifile) (default=$OCDATAROOT/common/l2brsgen_product_table.dat) = product table```
   
   ```datamin (float) (default=0.0) = minimum value for data scaling```
        
        default see product_table
   
   ```datamax (float) (default=0.0) = maximum value for data scaling```
        
        default see product_table
   
   ```stype (int) (default=0) = scaling type (default see product_table)```
        
          1: LINEAR
        
          2: LOG
   
   ```oformat (string) (default=HDF4) = format of the output file```
        
        hdf4: (1) HDF browse file
        
        png:  (5) PNG color or grayscale image file
        
        ppm:  (7) PPM color or PGM grayscale image file

