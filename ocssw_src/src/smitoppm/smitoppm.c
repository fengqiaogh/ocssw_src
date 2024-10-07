/* ========================================================================
 * smitoppm - converts any SeaWiFS SMI map file to a ppm file
 * 
 * Synopsis:
 *
 *   smitoppm input-smi-filename > output-ppm-filename
 *
 * Description:
 * 
 *   Uses the SeaWiFS L3_Mapfiles interface to read the byte-scaled image and
 *   palette.  Uses the bytescaled image to index into the palette, and writes
 *   a three-color image to standard output. The filename, product type,
 *   binning period, and other relevant meta-data are stored in the ppm file
 *   header as ascii comments.
 *
 * Modification history:
 *
 *     Programmer     Organization      Date      Description of change
 *   --------------   ------------    --------    ---------------------
 *   Bryan A. Franz   GSC             10/06/97    Original development
 *   Joel M. Gales    Futuretech      09/04/01    Support non-standard 
 *                                                smi maps
 *   Joel M. Gales    Futuretech      05/14/10    Don't get quality for
 *                                                HDF5 files & use datasize
 *                                                instead of datatype for 
 *                                                non-byte processing.
 *   Joel M. Gales    Futuretech      09/21/11    Add support for arctan
 *                                                scaling (SSS)
 *   Joel M. Gales    Futuretech      06/14/13    Add support for NETCDF4
 *
 *   Joel M. Gales    Futuretech      09/27/13    Removed netcdf code
 *                                                HDF5 calls can be used
 *                                                for NETCDF4 files
 * ======================================================================== */

#include <stdio.h>
#include <unistd.h>

typedef unsigned char byte;
#include "netcdf.h"
#include "map.h"
#include "mapproto_o.h"
#include "hdf5.h"
#include "mapproto.h"

#include <math.h>

#define NX          4096                 /* Map grid x-dimension           */
#define NY          2048                 /* Map grid y-dimension           */
#define MAXSTRNGLEN 255                  /* Max generic string length      */

#define CMD_ARGS    "b:w:s"   /* Valid commandline options */

#define VERSION "3.4"

/* ------------------------------------------------------------------------- *
 * usage - display proper calling sequence for this program                  *
 * ------------------------------------------------------------------------- */
void usage(char *progname) {
    fprintf(stderr, "%s %s (%s %s)\n", progname, VERSION, __DATE__, __TIME__);
    fprintf(stderr, "\nUsage: %s input-smi-filename > output-ppm-filename",
            progname);
    fprintf(stderr, "\n");

}

/* ------------------------------------------------------------------------ *
 *                              main                                        *
 * ------------------------------------------------------------------------ */
int main(int argc, char *argv[]) {

    int32 ix; /* Longitudinal bin index               */
    int32 iy; /* Latitudinal bin index                */
    int32 ic; /* Color bin index                      */
    /*    byte        rgb[NY][NX][3];*/ /* Three-color image array              */
    byte *rgb; /* Three-color image array              */

    int16 syear;
    int16 sday;
    int32 smsec;
    int16 eyear;
    int16 eday;
    int32 emsec;
    int32 status;
    int32 nrows;
    int32 ncols;
    size_t datasize;

    char prod_type[MAXSTRNGLEN];
    char l3m_name[MAXSTRNGLEN];

    /*    byte        image[NY][NX];*/
    byte *image = NULL;
    byte *quality;
    byte palette[3 * 256];
    meta_struct meta_l3m;

    uint16 ui16;
    uint8 best = 0, worst = 254;
    uint8 bio = 0;
    float32 flt32;

    uint16 fill_val_int16 = 65535;
    float32 fill_val_float = -32767.0;

    int maxval = 254;
    int islog = 0;
    int isatan = 0;
    double slope;
    double intercept;

    /* Parameters for getopt() */
    extern int opterr;
    extern int optind;
    extern char *optarg;
    int c;

    if (argc == 1) {
        usage("smitoppm");
        return 0;
    }

    /* Process optional command-line arguments */
    while ((c = getopt(argc, argv, CMD_ARGS)) != EOF) {
        switch (c) {
        case 'b':
            best = atoi(optarg);
            break;
        case 'w':
            worst = atoi(optarg);
            break;
        case 's':
            bio = 1;
            break;
        }
    }

    /* 
     * Read SMI palette and meta-data
     */
    status = get_l3m(argv[optind],
            &syear, &sday, &smsec,
            &eyear, &eday, &emsec,
            prod_type, NULL,
            NULL, palette, &meta_l3m);

    if (status < 0) {
        fprintf(stderr, "%s: Error accessing %s as L3 SMI file\n",
                argv[0], argv[optind]);
        exit(1);
    }

    nrows = meta_l3m.nrows;
    ncols = meta_l3m.ncols;


    if (H5Fis_hdf5(argv[optind])) {
        if (H5Tget_size(meta_l3m.bintype) == 1) {
            image = (byte *) calloc(nrows*ncols, sizeof (byte));
            datasize = 1;
        } else if (H5Tget_size(meta_l3m.bintype) == 2) {
            image = (byte *) calloc(nrows*ncols, 2 * sizeof (byte));
            datasize = 2;
        } else if (H5Tget_size(meta_l3m.bintype) == 4) {
            image = (byte *) calloc(nrows*ncols, 4 * sizeof (byte));
            datasize = 4;
            /*
                fprintf(stderr,"PPM file cannot be generated from FLOAT L3m file.\n");
                exit(-1);
             */
        }
    } else {
        if (meta_l3m.bintype == DFNT_UINT8) {
            image = (byte *) calloc(nrows*ncols, sizeof (byte));
            datasize = 1;
        } else if (meta_l3m.bintype == DFNT_UINT16) {
            image = (byte *) calloc(nrows*ncols, 2 * sizeof (byte));
            datasize = 2;
        } else if (meta_l3m.bintype == DFNT_FLOAT32) {
            image = (byte *) calloc(nrows*ncols, 4 * sizeof (byte));
            datasize = 4;
            /*
            printf("PPM file cannot be generated from FLOAT L3m file.\n");
            exit(-1);
             */
        }
    }

    rgb = (byte *) calloc(nrows * ncols * 3, sizeof (byte));
    quality = (byte *) calloc(nrows*ncols, sizeof (byte));

    /* 
     * Read SMI image
     */
    status = get_l3m(argv[optind],
            &syear, &sday, &smsec,
            &eyear, &eday, &emsec,
            prod_type, "l3m_name",
            &image[0], NULL, NULL);

    if (status < 0) {
        fprintf(stderr, "%s: Error accessing %s as L3 SMI file\n",
                argv[0], argv[optind]);
        exit(1);
    }

    status = get_l3m(argv[optind],
            &syear, &sday, &smsec,
            &eyear, &eday, &emsec,
            prod_type, "l3m_qual",
            &quality[0], NULL, NULL);

    if (status < 0) {
        best = 0;
        worst = 254;
    }

    if (datasize == 2) {

        for (iy = 0; iy < nrows; iy++) {
            for (ix = 0; ix < ncols; ix++) {

                memcpy(&ui16, &image[2 * (iy * ncols + ix)], 2);

                if (ui16 == fill_val_int16)
                    image[iy * ncols + ix] = 255;
                else {
                    if (bio) {
                        image[iy * ncols + ix] = (byte) (ui16 / 256);
                        if (ui16 == 32767) image[iy * ncols + ix] = 128;
                    } else if (strcmp(l3m_name, "Number of pixels per bin") == 0) {
                        image[iy * ncols + ix] =
                                (byte) (255 * (ui16 - meta_l3m.scaled_data_min) /
                                (meta_l3m.scaled_data_max - meta_l3m.scaled_data_min));
                    } else {
                        image[iy * ncols + ix] = (byte) (250 * ui16 / 65534);
                    }
                }
            } // col loop
        } // row loop

    } else if (datasize == 4) {

        if (strcmp(meta_l3m.scaled_data_type, "LOG") == 0) {
            islog = 1;
            intercept = log10(meta_l3m.scaled_data_min);
            slope = (log10(meta_l3m.scaled_data_max) - intercept) / maxval;
        } else if (strcmp(meta_l3m.scaled_data_type, "ATAN") == 0) {
            isatan = 1;
            intercept = 17.5;
            slope = 125.0;
        } else {
            intercept = meta_l3m.scaled_data_min;
            slope = (meta_l3m.scaled_data_max - intercept) / maxval;
        }
        for (iy = 0; iy < nrows; iy++) {
            for (ix = 0; ix < ncols; ix++) {

                memcpy(&flt32, &image[4 * (iy * ncols + ix)], 4);

                if (flt32 < (fill_val_float + 1E-5))
                    image[iy * ncols + ix] = 255;
                else if (flt32 <= meta_l3m.scaled_data_min)
                    image[iy * ncols + ix] = (uint8) 0;
                else if (flt32 >= meta_l3m.scaled_data_max)
                    image[iy * ncols + ix] = (uint8) maxval;
                else {
                    if (islog)
                        image[iy * ncols + ix] = (uint8) ((log10(flt32) - intercept) / slope);
                    else if (isatan)
                        image[iy * ncols + ix] = (uint8) slope * ((atan(0.5 * flt32 - intercept) / atan(intercept)) + 1);
                    else
                        image[iy * ncols + ix] = (uint8) ((flt32 - intercept) / slope);
                }
            }
        }
    }


    /* 
     * Convert from greyscale and palette to rgb array
     */
    for (iy = 0; iy < nrows; iy++)
        for (ix = 0; ix < ncols; ix++)
            for (ic = 0; ic < 3; ic++)
                rgb[iy * (ncols * 3) + ix * 3 + ic] =
                    palette[image[iy * ncols + ix] * 3 + ic ];


    /* Quality Filtering */
    /* ----------------- */
    for (iy = 0; iy < nrows; iy++)
        for (ix = 0; ix < ncols; ix++)
            if ((quality[iy * ncols + ix] < best) ||
                    (quality[iy * ncols + ix] > worst))
                for (ic = 0; ic < 3; ic++)
                    rgb[iy * (ncols * 3) + ix * 3 + ic] =
                        palette[255 * 3 + ic ];


    /*                      
     * Write output ppm file
     */
    printf("P6\n");
    printf("# File:       %s\n", argv[optind]);
    printf("# Dataset:    %s\n", l3m_name);
    printf("# Bin Period: %s\n", prod_type);
    printf("# Start Time: %d %d %d\n", syear, sday, smsec);
    printf("# End Time:   %d %d %d\n", eyear, eday, emsec);
    printf("%d\n", ncols);
    printf("%d\n", nrows);
    printf("255\n");
    fwrite(&rgb[0], 1, ncols * nrows * 3, stdout);

    free(image);
    free(quality);
    free(rgb);

    exit(0);
}




