make_L3 README.txt

Compiling
---------

1) set environment variables for HDBLIB and HDFINC to the appropriate HDF4.1r3 
/lib and /include directories respectively:

% setenv HDFLIB /usr/local/HDF4.1r3/lib
% setenv HDFLIB /usr/local/HDF4.1r3/include

2) go to the LIB directory and run the shell script 'arch'

% cd ./LIB 
% arch

This step will create a library that  is subsequently used by the make_L3 
program.

3) edit path.h to reflect the appropriate directories

path.h defines several standard directories used by make_L3

4) compile make_L3_v1.0.1.c

% cc -O2 -fullwarn -woff 1209,1210 -64 make_L3_v1.0.1.c -o make_L3 /
./LIB/lib.a -I./LIB -I$HDFINC -L$HDFLIB -lmfhdf -ldf -lz -ljpeg /
-lm -lmalloc

Running
-------

make_L3 takes input from standard input, so the SeaWIFS level 2 files
should be piped to the program via the command line as such:

% ls /data6/seawifs/L2/S1998018* | make_L3 -pixsz=4633.22 -proj=1 /
-ndvi -evi -angles -method=minsenz -TOA -maxmem=500 -bufstep=900

Some command line parameters:

make_L3 <parameters>

where <parameters> are a combination of:

-pixsz           -> the pixel size in meters
-proj=<type>     -> the projection, where <type>:
	            1 -> MODIS Sinusoidal
-ndvi            -> add the NDVI layer
-evi             -> add the EVI layer
-angles          -> add all of the angle layers (sensor and sun)
-method=<type>   -> compositing method where <type>:
	            minsenz -> minimum sensor zenith angle
-TOA             -> add the Top Of Atmosphere layer
-maxmem=<size>   -> the maximum memory, in Megabytes, to be 
                    allocated by the program to buffer the 
                    compositing.  The bigger, the better.  Ideally
                    it should be 1700 to fit the entire global
                    image (with all the layers) at once.  Otherwise
                    the process will be buffered on disk.
-bufstep=<lines> -> the number of lines to move "down" after 
                    each buffered horizontal segment has been
                    processed.  This is essentially an overlap
                    between subsequent horizontal strips of the
                    buffer, required because of the skewness 
                    of the map projection.  If the number of 
                    lines is too small, then the compositing 
                    process may be slow.  But if the number of
                    lines is too big, then there may be blank
                    horizontal gaps on the final image.








