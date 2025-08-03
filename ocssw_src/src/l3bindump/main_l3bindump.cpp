/*
 This program extracts data from level-3 HDF bin files and
 writes them out as an ASCII table.

 Regions of interest can be specified by latitude and longitude
 boundaries or by a central coordinate and a radius in kilometers.

 Norman Kuring		14-Feb-2003	Original development
 Norman Kuring		11-Dec-2007	Fix memory-overrun bug and add a
 couple of calls to free().
 Norman Kuring		21-Dec-2011	Give a precision when printing out
 bit strings to avoid unwanted printing
 of uninitialized memory.
 Norman Kuring		21-Mar-2013	Change the latbin array from 32-bit
 floats to 64-bit to reduce rounding-
 error effects on bin numbers at smaller
 bin sizes.  Thanks to George White for
 pointing this one out.
 Sean Bailey        25-Mar-2014 Redesigned to use libbin++ to be able
 to work with new netCDF4 formatted files, as well as the Aquarius HDF5
 versions, modified output to include mean and stdev but exclude selcat,
 flag and time trend bits - Joel Gales made numerous mods to libbin++
 for this to work...
 Don Shea           28-Luly-2015 Started using L3File class to read bin files. 
 Changed the looping which increased the speed quite a bit.
 */
#include "l3bindump.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <clo.h>
#include <strings.h>
#include <L3File.h>

#include "compareObj.h"

#define EARTH_RADIUS 6371.229
#define BIN_NOT_FOUND 110

using namespace l3;

void printBin(FILE *outfp, int32_t row, int32_t col, L3Bin *l3Bin, L3Shape* shape) {
    int64_t bin = l3Bin->getBinNum();
    float lat, lon;
    shape->rowcol2latlon(row, col, lat, lon);
    float north, south, east, west;
    shape->rowcol2bounds(row, col, north, south, east, west);

    fprintf(outfp, "%07ld %9.5f %10.5f %9.5f %9.5f %10.5f %10.5f ", (long)bin,
            lat, lon, north, south, west, east);
    fprintf(outfp, "%4d %3d ", l3Bin->getNobs(), l3Bin->getNscenes());
    fprintf(outfp, "% .8e % .8e % .8e ", l3Bin->getSum(), l3Bin->getSumSquares(), l3Bin->getWeights());
    fprintf(outfp, "%10.5f %10.5f", l3Bin->getMean(), l3Bin->getStdev());
    fprintf(outfp, "\n");
}

void printWholeFile(FILE *outfp, L3File &l3File) {
    L3Shape* shape = l3File.getShape();
    int32_t row, col;
    L3Bin *l3Bin;
    bool binFound = false;
    for (row = 0; row < shape->getNumRows(); row++) {
        for (col = 0; col < shape->getNumCols(row); col++) {
            l3Bin = l3File.getBin(row, col);
            if (l3Bin) {
                binFound = true;
                printBin(outfp, row, col, l3Bin, shape);
            }
        }
    }
    if (!binFound) {
        fprintf(stderr, "No bins found in file.\n");
        exit(BIN_NOT_FOUND);
    }
    exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[]) {
    want_verbose = 0;

    int32_t i;
    instr *input; // input parameters structure

    FILE *outfp = stdout;
    char outfile[FILENAME_MAX] = "\0";

    L3File l3File;

    /* hold all of the command line options */
    clo_optionList_t* list;

    if (argc == 1) {
        l3bindump_usage();
        return 0;
    }

    for (i = 0; i < argc; i++) {
        // see if help on command line
        if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "-help") == 0)) {
            l3bindump_usage();
            return 0;
        }
    }

    input = (instr*) malloc(sizeof (instr));
    if (!input) {
        fprintf(stderr, "-E- %s %d: Error allocating input structure.\n", __FILE__,
                __LINE__);
        exit(EXIT_FAILURE);
    }

    l3bindump_input_init(input);

    // create an option list
    list = clo_createList();

    /* initialize the option list with descriptions and default values */
    l3bindump_init_options(list);

    // Parse input parameters

    if (l3bindump_read_options(list, argc, argv) != 0) {
        fprintf(stderr, "-E- %s: Error reading program options.\n", __FILE__);
        return (EXIT_FAILURE);
    }

    /*                                                                  */
    /* Now, loop through command arguments again and update input struct*/
    /*                                                                  */
    if (l3bindump_load_input(list, input) != 0) {
        fprintf(stderr, "-E- %s: Error loading options into input structure.\n",
                __FILE__);
        return (EXIT_FAILURE);
    }

    if (input->verbose == 1) {
        want_verbose = 1;
    }

    // open the file
    if (!l3File.open(input->ifile)) {
        fprintf(stderr, "%s: Error: Unable to open %s for reading.\n",
                argv[0], input->ifile);
        exit(EXIT_FAILURE);
    }

    // set the product to read
    if (!l3File.setActiveProductList(input->l3bprod)) {
        fprintf(stderr, "%s: Error: Product %s not found in file.\n",
                argv[0], input->l3bprod);
        exit(EXIT_FAILURE);
    }

    // Open output file if requested
    if (strlen(input->ofile) > 0) {
        strcpy(outfile, input->ofile);
        if ((outfp = fopen(outfile, "w")) == NULL) {
            fprintf(stderr, "%s: Error: Unable to open %s for writing.\n",
                    argv[0], outfile);
            exit(EXIT_FAILURE);
        }
    }

    //
    // write out the header info
    //
    if (strcmp(input->oformat, "SeaBASS") == 0) {
        fprintf(outfp, "/begin_header\n");
        fprintf(outfp, "/delimiter=space\n");
        fprintf(outfp, "/missing=-999\n");
        if (input->north > -999) {
            fprintf(outfp, "/north_latitude=%10.5f\n", input->north);
            fprintf(outfp, "/south_latitude=%10.5f\n", input->south);
            fprintf(outfp, "/east_longitude=%10.5f\n", input->east);
            fprintf(outfp, "/west_longitude=%10.5f\n", input->west);
        }
        fprintf(outfp, "/fields=bin,lat,lon,north,south,west,east,nobs,nscenes,sum,sum_squared,weight,mean,stdev\n");
        fprintf(outfp, "/end_header\n");

    } else {
        fprintf(outfp, "%80s%15.15s %15.15s\n", " ", input->l3bprod, input->l3bprod);
        fprintf(outfp, "    bin centerlat  centerlon");
        fprintf(outfp, "     north     south       west       east");
        fprintf(outfp, "    n   N             sum     sum_squared          weight");
        fprintf(outfp, "       mean      stdev\n");
        fprintf(outfp, "------- --------- ----------");
        fprintf(outfp, " --------- --------- ---------- ----------");
        fprintf(outfp, " ---- --- --------------- --------------- ---------------");
        fprintf(outfp, " ---------- ----------\n");
    }

    if (want_verbose)
        fprintf(outfp, "! Input file: %s\n", input->ifile);

    CompareObj *compareObj = NULL;
    L3Bin *l3Bin;
    L3Shape *shape = l3File.getShape();
    int32_t row, col;
    int32_t startRow=0, stopRow=0;
    if (want_verbose)
        fprintf(outfp, "! Number of Rows: %d\n", shape->getNumRows());

    // handle lat/lon without radius
    if(input->radius == -999 && input->lat != -999 && input->lon != -999) {
        input->bin_number = shape->latlon2bin(input->lat, input->lon);
    }

    if (input->bin_number > 0) {
        /* Input arguments are: main_file_path parameter bin_number. */
        shape->bin2rowcol(input->bin_number, row, col);

        l3Bin = l3File.getBin(row, col);
        if (l3Bin) {
            printBin(outfp, row, col, l3Bin, shape);
            exit(EXIT_SUCCESS);
        } else {
            fprintf(stderr, "Requested bin number = %ld was not in the file.\n", (long)input->bin_number);
            exit(BIN_NOT_FOUND);
        }


    } else if (input->radius > 0) {
        /* Input arguments are: main_file_path parameter lat lon radius. */

        input->lat = constrainLat(input->lat);
        input->lon = constrainLon(input->lon);

        float radius_degrees = (input->radius / EARTH_RADIUS) * OEL_RADEG;
        if (radius_degrees > 180) {

            /* The entire globe has been selected. */
            printWholeFile(outfp, l3File);

        } else {
            double north, south;

            north = input->lat + radius_degrees;
            south = input->lat - radius_degrees;
            if (north >= 90.0) {
                stopRow = shape->getNumRows();
            } else {
                stopRow = shape->lat2row(north) + 1; // add an extra row
                if (stopRow >= shape->getNumRows())
                    stopRow = shape->getNumRows();
            }

            if (south <= -90.0) {
                startRow = 0;
            } else {
                startRow = shape->lat2row(south) - 1; // add an extra row
                if (startRow < 0)
                    startRow = 0;
            }
            compareObj = new CompareObjCircle(input->lat, input->lon, input->radius);
        }


    } else if (input->north != -999) {
        /* Input arguments are: main_file_path parameter north south west east. */

        input->north = constrainLat(input->north);
        input->south = constrainLat(input->south);
        input->west = constrainLon(input->west);
        input->east = constrainLon(input->east);

        if (input->south > input->north) {
            double tmp;

            fprintf(stderr, "-W- %s line %d: ", __FILE__, __LINE__);
            fprintf(stderr,
                    "Specified south latitude is greater than specified north latitude.\n");
            fprintf(stderr, "I will swap the two.\n");
            tmp = input->north;
            input->north = input->south;
            input->south = tmp;
        }

        if (input->north == 90.0 && input->south == -90.0 &&
                input->east == 180.0 && input->west == -180.0) {

            /* The entire globe has been selected. */
            printWholeFile(outfp, l3File);

        }

        stopRow = shape->lat2row(input->north);
        startRow = shape->lat2row(input->south);
        compareObj = new CompareObjLonOnly(input->east, input->west);

    } else {
        l3bindump_usage();
        return (EXIT_FAILURE);
    }

    float lat, lon;
    bool binFound = false;
    for (row = startRow; row <= stopRow; row++) {
        for (col = 0; col < shape->getNumCols(row); col++) {
            shape->rowcol2latlon(row, col, lat, lon);
            if (compareObj->isInside(lat, lon)) {
                l3Bin = l3File.getBin(row, col);
                if (l3Bin) {
                    binFound = true;
                    printBin(outfp, row, col, l3Bin, shape);
                } // bin found in file
            } // if inside
        } // col
    } // row
    if (!binFound) {
        fprintf(stderr, "No bins found in file.\n");
        exit(BIN_NOT_FOUND);
    }

    return (EXIT_SUCCESS);
} /* End of main() function */
