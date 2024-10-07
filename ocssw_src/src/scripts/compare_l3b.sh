#!/bin/bash
#
# compare two L3 bin files since nccmp does not (yet) compare user defined types
#

# look for the leave files option
cleanup=1
if [ $# -eq 3 ]; then
    if [ $1 = "-l" ]; then
	cleanup=0
	shift
    fi
fi


if [ $# -ne 2 ]; then
    progName=`basename $0`
    echo "usage: $progName [-l] <binFile1> <binFile2>"
    echo
    echo "    convert bin files to raw projection map files of the mean for all"
    echo "    variables using l3mapgen then use nccmp to compare the map files."
    echo
    echo "    -l  leave temporary map files, no cleanup."
    echo
    exit 1
fi

mapFile1=$1.map.nc
mapFile2=$2.map.nc

l3mapgen "ifile=$1" "ofile=$mapFile1" "oformat=netcdf4" "projection=raw"
if [ $? -ne 0 ]; then
    echo "l3mapgen failed for $1"
    if [ $cleanup -eq 1 ]; then
        rm -f $mapFile1
    fi
    exit 1
fi

l3mapgen "ifile=$2" "ofile=$mapFile2" "oformat=netcdf4" "projection=raw"
if [ $? -ne 0 ]; then
    echo "l3mapgen failed for $2"
    if [ $cleanup -eq 1 ]; then
        rm -f $mapFile1 $mapFile2
    fi
    exit 1
fi

nccmp -m -g -d -f -C 10 -G date_created,software_version,_lastModified,ifile,ofile $OCTEST_TOLERANCE $mapFile1 $mapFile2

exitCode=$?

if [ $cleanup -eq 1 ]; then
    rm -f $mapFile1 $mapFile2
fi

exit $exitCode

