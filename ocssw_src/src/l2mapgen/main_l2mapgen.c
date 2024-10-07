/*
Inputs:
        see usage...
Output:
        Eight-bit PGM data stream.

Code draws heavily on swplatecar, smigen, and l2bin

Original Development: 
        Sean Bailey -- 30 June 2006
	
Modifications:
        SWB, 17 Jul 2008, fixed the threshold parameter - long/float discrepancy

 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

#include <netcdf.h>
#include <png.h>
#include "miscstruct.h"
#include "miscanfill.h"
#include "mipoly.h"
#include "mfhdf.h"
#include <clo.h>
#include "l2mapgen.h"
#include "l2mapgen_input.h"
#include <genutils.h>
#include <readL2scan.h>

#include <stdio.h>
#include <X11/X.h>
#include <X11/Xlib.h>
#include <unistd.h>

#include <geotiffio.h>
#include <geo_normalize.h>
#include <geo_tiffp.h>
#include <geo_keyp.h>
#include <xtiffio.h>
#include <geokeys.h>


#define ROOT2 1.4142135623730950488016887242096981
#define PI 3.14159265358979323846
#define BINBELOWTHRESH  110

#define CALLOC(ptr,typ,num) {      \
  (ptr) = (typ *)calloc((num) , sizeof(typ));    \
  if((ptr) == NULL){       \
    fprintf(stderr,"-E- %s line %d: Memory allocation failure.\n", \
    __FILE__,__LINE__);       \
    exit(EXIT_FAILURE);       \
  }         \
}

/* function prototypes */
int scan_convert(XPoint *ptsIn);
int collect_bins(int number_of_initial_points, XPoint *initial_point,
        int *span_width);
int miCreateETandAET();
int miInsertionSort();

static int *binsoverlain = NULL;
static int binsperpixel = 0;
static int numoverlain;
static int width, height;
static double scale;

int32_t l3m_params = 0;
char **parmname_list;
char **parmname_short;
char **unit_list;
char **scaling_list;
float *maximum_list;
float *minimum_list;
char **palette_list;
char **precision_list;

int main(int argc, char *argv[]) {

    static instr input;
    static l2_prod l2_str[MAXNFILES];
    static meta_l2Type meta_l2;
    static char buf[65535];
    char qual[9];

    int32_t npix, nscans;
    int32_t nfiles;
    int status;
    int32_t l2_flags;
    byte quality;
    byte goodpix;
    int prodidx = -1;
    int qualidx = -1;
    double latmax = -90.0;
    double latmin = 90.0;
    double lonmin = 180.0;
    double lonmax = -180.0;

    double *sum;
    unsigned char *mean;
    unsigned char *mask;
    int *count;

    uint32_t numbins, numoutpix;
    XPoint corners[8]; /* actually, 4 pixel corners + 4 side midpoints */
    float64 threshold;
    float nrad, srad, wrad, erad;
    int i, j, k, m, n, ii;
    int progress;
    uint32_t mask_252, mask_253, all_masks;
    uint32_t mask_glint = 0;
    uint32_t required_mask;

    char scale_type[50];
    uint8_t default_palfile = 0;
    char *cptr;

    static uint16_t maxval = 250;
    float val;
    float slope = 1.0;
    float intercept = 0.0;
    int32_t prod_num = -1;
    byte *rgb;
    short *r, *g, *b;

    FILE *fp = NULL;
    FILE *outfp = NULL;


    /* hold all of the command line options */
    clo_optionList_t* list = clo_createList();

    if (argc == 1 || strcmp(argv[1], "-version") == 0 || strcmp(argv[1], "--version") == 0) {
        want_verbose = 0;
        l2mapgen_init_options(list);
        l2mapgen_read_options(list, argc, argv);
        clo_printUsage(list);
        exit(EXIT_FAILURE);
    }

    if (l2mapgen_input(argc, argv, list, &input) != 0) {
        clo_printUsage(list);
        exit(EXIT_FAILURE);
    }

    if (strcasecmp(input.palfile, "default") == 0) {
        default_palfile = 1;
    }

    /* Read product table */
    fp = fopen(input.product_table, "r");
    if (fp == 0x0) {
        fprintf(stderr, "SMIGEN product table \"%s\" not found.\n", input.product_table);
        exit(EXIT_FAILURE);
    }

    l3m_params = 0;
    while (fgets(buf, 128, fp) != NULL) {
        if ((buf[0] >= 0x41) && (buf[0] <= 0x5a)) l3m_params++;
    }
    fseek(fp, 0, SEEK_SET);

    parmname_list = (char**) calloc(l3m_params, sizeof (char *));
    parmname_short = (char**) calloc(l3m_params, sizeof (char *));
    unit_list = (char**) calloc(l3m_params, sizeof (char *));
    scaling_list = (char**) calloc(l3m_params, sizeof (char *));
    palette_list = (char**) calloc(l3m_params, sizeof (char *));
    maximum_list = (float *) calloc(l3m_params, sizeof (float));
    minimum_list = (float *) calloc(l3m_params, sizeof (float));
    precision_list = (char **) calloc(l3m_params, sizeof (char *));

    i = 0;
    while (fgets(buf, 128, fp) != NULL) {
        if ((buf[0] >= 0x41) && (buf[0] <= 0x5a)) {

            cptr = strtok(buf, ":");
            parmname_list[i] = (char*) malloc(strlen(cptr) + 1);
            strcpy(parmname_list[i], cptr);

            cptr = strtok(NULL, ":");
            parmname_short[i] = (char*) malloc(strlen(cptr) + 1);
            strcpy(parmname_short[i], cptr);

            cptr = strtok(NULL, ":");
            unit_list[i] = (char*) malloc(strlen(cptr) + 1);
            strcpy(unit_list[i], cptr);

            cptr = strtok(NULL, ":");
            scaling_list[i] = (char*) malloc(strlen(cptr) + 1);
            strcpy(scaling_list[i], cptr);

            cptr = strtok(NULL, ":");
            minimum_list[i] = (float) atof(cptr);

            cptr = strtok(NULL, ":");
            maximum_list[i] = (float) atof(cptr);

            cptr = strtok(NULL, ":");
            precision_list[i] = (char*) malloc(strlen(cptr) + 1);
            strcpy(precision_list[i], cptr);

            cptr = strtok(NULL, "\n");
            palette_list[i] = (char*) malloc(strlen(cptr) + 1);
            strcpy(palette_list[i], cptr);

            i++;
        }
    }
    fclose(fp);


    for (i = 0; i < l3m_params; i++) {

        if (strcmp(parmname_short[i], input.prod) == 0) {

            prod_num = i;

            /* define scaling */
            strcpy(scale_type, scaling_list[i]);
            if (input.stype == 1) strcpy(scale_type, "linear");
            if (input.stype == 2) strcpy(scale_type, "logarithmic");

            if (input.datamin == 0.0) input.datamin = minimum_list[i];
            if (input.datamax == 0.0) input.datamax = maximum_list[i];


            if (strcmp(scale_type, "linear") == 0) {

                strcpy(scale_type, "LINEAR");
                intercept = input.datamin;
                slope = (input.datamax - intercept) / maxval;
            } else if (strcmp(scale_type, "logarithmic") == 0) {

                strcpy(scale_type, "LOG");

                intercept = log10(input.datamin);
                slope = (log10(input.datamax) - intercept) / maxval;
            }

            /* Read palette file */
            if (default_palfile) {
                strcpy(input.palfile, input.palette_dir);
                strcat(input.palfile, "/");
                strcat(input.palfile, palette_list[i]);
                strcat(input.palfile, ".pal");
            }

            if (!(r = (short *) calloc(256, sizeof (short)))) {
                fprintf(stderr, "smigen: Error allocating space for red.\n");
                return -1;
            };
            if (!(g = (short *) calloc(256, sizeof (short)))) {
                fprintf(stderr, "smigen: Error allocating space for green.\n");
                return -1;
            };
            if (!(b = (short *) calloc(256, sizeof (short)))) {
                fprintf(stderr, "smigen: Error allocating space for blue.\n");
                return -1;
            };

            if (getlut_file(input.palfile, r, g, b)) {
                fprintf(stderr, "Error reading palette file %s\n", input.palfile);
                free(r);
                free(g);
                free(b);
                return -1;
            }
            if (input.mask) {
                r[252] = 128;
                g[252] = 128;
                b[252] = 128;
                r[253] = 160;
                g[253] = 82;
                b[253] = 45;
                r[254] = 255;
                g[254] = 255;
                b[254] = 255;
                r[255] = 0;
                g[255] = 0;
                b[255] = 0;
            }
            for (i = 0; i < 256; i++) {
                input.palette[i * 3] = r[i];
                input.palette[i * 3 + 1] = g[i];
                input.palette[i * 3 + 2] = b[i];
            }

            break;
        } // if prod matches
    } // for i

    if (prod_num == -1) {
        /* define scaling */
        strcpy(scale_type, "LINEAR");
        if (input.stype == 1) strcpy(scale_type, "linear");
        if (input.stype == 2) strcpy(scale_type, "logarithmic");

        if (input.datamin == 0.0) input.datamin = 0.001;
        if (input.datamax == 0.0) input.datamax = 1.0;


        if (strcmp(scale_type, "linear") == 0) {

            strcpy(scale_type, "LINEAR");
            intercept = input.datamin;
            slope = (input.datamax - intercept) / maxval;
        }

        if (strcmp(scale_type, "logarithmic") == 0) {

            strcpy(scale_type, "LOG");

            intercept = log10(input.datamin);
            slope = (log10(input.datamax) - intercept) / maxval;
        }

        /* Read default.pal palette file */
        strcpy(input.palfile, input.palette_dir);
        strcat(input.palfile, "/default.pal");

        if (!(r = (short *) calloc(256, sizeof (short)))) {
            fprintf(stderr, "smigen: Error allocating space for red.\n");
            return -1;
        };
        if (!(g = (short *) calloc(256, sizeof (short)))) {
            fprintf(stderr, "smigen: Error allocating space for green.\n");
            return -1;
        };
        if (!(b = (short *) calloc(256, sizeof (short)))) {
            fprintf(stderr, "smigen: Error allocating space for blue.\n");
            return -1;
        };

        if (getlut_file(input.palfile, r, g, b)) {
            fprintf(stderr, "Error reading palette file %s\n", input.palfile);
            free(r);
            free(g);
            free(b);
            return -1;
        }
        if (input.mask) {
            r[252] = 128;
            g[252] = 128;
            b[252] = 128;
            r[253] = 160;
            g[253] = 82;
            b[253] = 45;
            r[254] = 255;
            g[254] = 255;
            b[254] = 255;
            r[255] = 0;
            g[255] = 0;
            b[255] = 0;
        }
        for (i = 0; i < 256; i++) {
            input.palette[i * 3] = r[i];
            input.palette[i * 3 + 1] = g[i];
            input.palette[i * 3 + 2] = b[i];
        }
    }

    /* Single HDF input */
    /* ---------------- */
    int singleInputFile = 0;
    if (Hishdf(input.ifile) == TRUE) {
        singleInputFile = 1;
    } else {
        int fid;
        status = nc_open(input.ifile, NC_NOWRITE, &fid);
        if (status == NC_NOERR) {
            singleInputFile = 1;
            nc_close(fid);
        }
    }
    if (singleInputFile) {
        nfiles = 1;
        status = openL2(input.ifile, 0x0, &l2_str[0]);

        status = readL2meta(&meta_l2, 0);
        if (meta_l2.northlat >= latmax) latmax = meta_l2.northlat;
        if (meta_l2.southlat <= latmin) latmin = meta_l2.southlat;
        if (meta_l2.westlon <= lonmin) lonmin = meta_l2.westlon;
        if (meta_l2.eastlon >= lonmax) lonmax = meta_l2.eastlon;

        closeL2(&l2_str[0], 0);
        fprintf(stderr, "Single HDF input\n");
    } else {

        /* Filelist input - Determine number of input files */
        /* ------------------------------------------------ */
        nfiles = 0;
        fp = fopen(input.ifile, "r");
        if (fp == NULL) {
            fprintf(stderr, "Input listing file: \"%s\" not found.\n", input.ifile);
            return -1;
        }
        while (fgets(buf, 256, fp) != NULL) nfiles++;
        fclose(fp);
        fprintf(stderr, "%d input files\n", nfiles);


        /* Open L2 input files  */
        /* ------------------- */
        fp = fopen(input.ifile, "r");
        for (i = 0; i < nfiles; i++) {
            fgets(buf, 256, fp);
            buf[strlen(buf) - 1] = 0;
            status = openL2(buf, 0x0, &l2_str[i]);

            status = readL2meta(&meta_l2, i);
            if (meta_l2.northlat >= latmax) latmax = meta_l2.northlat;
            if (meta_l2.southlat <= latmin) latmin = meta_l2.southlat;
            if (meta_l2.westlon <= lonmin) lonmin = meta_l2.westlon;
            if (meta_l2.eastlon >= lonmax) lonmax = meta_l2.eastlon;

            closeL2(&l2_str[i], i);


        } /* ifile loop */
        fclose(fp);
    }
    /* Setup flag masking */
    /* --------------- */
    if (input.flaguse[0] == 0) {
        strcpy(input.flaguse, DEF_FLAG);
    } else {
        input.mask = 1;
    }
    strcpy(buf, l2_str[0].flagnames);
    setupflags(buf, "HIGLINT", &mask_252, &required_mask, &status,l2_str[0].l2_bits);
    setupflags(buf, "LAND", &mask_253, &required_mask, &status,l2_str[0].l2_bits);
    setupflags(buf, input.flaguse, &all_masks, &required_mask, &status,l2_str[0].l2_bits);

    if ((all_masks & mask_252) != 0)
        mask_glint = 1;

    /* Make sure L2 product exists for every input L2 file */
    /* ---------------------------------------------------------------- */
    status = 0;
    for (j = 0; j < nfiles; j++) {
        for (i = 0; i < l2_str[j].nprod; i++) {
            if (strcmp(l2_str[j].prodname[i], input.prod) == 0) {
                status++;
            }
        }
    }

    if (status != nfiles) {
        fprintf(stderr, "Product %s not found in all L2 file(s)\n", input.prod);
        exit(EXIT_FAILURE);
    }

    /* Get the box boundaries. */
    if (input.north == input.south && input.west == input.east) {
        input.north = latmax;
        input.south = latmin;
        input.west = lonmin;
        input.east = lonmax;
    }

    nrad = input.north * PI / 180;
    srad = input.south * PI / 180;
    wrad = input.west * PI / 180;
    erad = input.east * PI / 180;

    fprintf(stderr, "Mapping data to:\n N: %8.5f\n S: %8.5f\n W: %8.5f\n E: %8.5f\n",
            input.north, input.south, input.west, input.east);

    fprintf(stderr, "Scale Type      : %s\n", scale_type);
    fprintf(stderr, "Data Min (abs)  : %8.4f\n", input.datamin);
    fprintf(stderr, "Data Max (abs)  : %8.4f\n", input.datamax);
    fprintf(stderr, "Scale Slope     : %8.4f\n", slope);
    fprintf(stderr, "Scale Intercept : %8.4f\n", intercept);
    if (input.apply_pal >= 1 || default_palfile == 0)
        fprintf(stderr, "Applying palette: %s\n", input.palfile);
    if (input.mask)
        fprintf(stderr, "Applying masking to flagged pixels\n");

    if (nrad <= srad) {
        fprintf(stderr, "The northernmost boundary must be greater than the ");
        fprintf(stderr, "southernmost boundary.\n");
        exit(EXIT_FAILURE);
    }

    /* Get the size of the output image. */
    width = input.width;
    if (width < 32) {
        fprintf(stderr, "Width (%d) is too small to produce a useful image.\n",
                input.width);
        exit(EXIT_FAILURE);
    }
    if (wrad < erad)
        height = rint((nrad - srad) * width / (erad - wrad));
    else
        height = rint((nrad - srad) * width / (erad - wrad + 2 * PI));

    numbins = width * height;
    scale = height / (nrad - srad);

    /* Allocate memory for the output data. */
    CALLOC(sum, double, numbins);
    CALLOC(count, int, numbins);
    CALLOC(mean, unsigned char, numbins);
    CALLOC(mask, unsigned char, numbins);
    CALLOC(rgb, byte, numbins * 3);
    /* Get the threshold percentage. */
    threshold = (double) input.threshold;

    /* For each input image... */
    for (k = 0; k < nfiles; k++) {
        fprintf(stderr, "Opening HDF file, %s ...\n", l2_str[k].filename);
        status = reopenL2(k, &l2_str[k]);

        nscans = l2_str[k].nrec;
        npix = l2_str[k].nsamp;
        fprintf(stderr, "pix: %d scan: %d\n", npix, nscans);

        if (strcmp(input.prod, "sst") == 0 || strcmp(input.prod, "sst4") == 0) {
            strcpy(qual, "qual_");
            strcat(qual, input.prod);
            for (i = 0; i < l2_str[k].nprod; i++) {
                if (strcmp(qual, l2_str[k].prodname[i]) == 0) {
                    qualidx = i;
                    break;
                }
            }
        }
        for (i = 0; i < l2_str[k].nprod; i++)
            if (strcmp(input.prod, l2_str[k].prodname[i]) == 0) {
                prodidx = i;
                break;
            }

        if (npix <= 0 || nscans <= 0) {
            fprintf(stderr,
                    "-E- %s line %d: Bad scene dimension: npix=%d nscans=%d in file, %s.\n",
                    __FILE__, __LINE__, npix, nscans, l2_str[k].filename);
            exit(EXIT_FAILURE);
        }

        /* Use the following variable to show this program's progress. */
        progress = (int) ceil(((double) nscans / 78));


        /* For each scan line ... */
        fprintf(stderr, "Mapping swath pixels to Plate Carree projection...\n");
        for (i = 0; i < nscans; i++) {

            /* Read swath record from L2 */
            /* ------------------------- */
            status = readL2(&l2_str[k], k, i, -1, NULL);

            /* For each pixel ... */
            for (j = 0; j < npix; j++) {

                double plat, plon; /* pixel center */
                double sinplat, cosplat;
                double c[2]; /* edge/corner distance */
                double sinc[2], cosc[2];
                double sinaz[8], cosaz[8];
                int out = 0; /* Number of pixel "corners"   */
                /* that fall outside the image */

                plat = l2_str[k].latitude[j] * PI / 180.0;
                plon = l2_str[k].longitude[j] * PI / 180.0;

                sinplat = sin(plat);
                cosplat = cos(plat);

                if (j < npix - 1) {

                    double nlat, nlon; /* next pixel center */
                    double dlat, dlon; /* delta lat and lon */
                    double sindlat2, sindlon2, cosnlat;
                    double az; /* azimuth to next pixel */
                    double phw; /* pixel half width */


                    nlat = l2_str[k].latitude[j + 1] * PI / 180.0;
                    nlon = l2_str[k].longitude[j + 1] * PI / 180.0;

                    cosnlat = cos(nlat);

                    dlat = nlat - plat;
                    dlon = nlon - plon;

                    sindlat2 = sin(dlat / 2);
                    sindlon2 = sin(dlon / 2);

                    phw = asin(sqrt(sindlat2 * sindlat2 + /* Page 30, */
                            sindlon2 * sindlon2 * /* Equation */
                            cosplat * cosnlat)); /*     5-3a */

                    /*
                    Make the pixel's coverage slightly broader to avoid
                    black pin holes in the output image where the estimated
                    pixel boundaries do not quite line up.  Don't increase
                    the size too much, however, because that will cause
                    unwanted image blurring.
                     */
                    phw *= 1.3;

                    c[0] = phw; /* center-to-edge   distance */
                    c[1] = phw * ROOT2; /* center-to-corner distance */

                    sinc[0] = sin(c[0]);
                    sinc[1] = sin(c[1]);
                    cosc[0] = cos(c[0]);
                    cosc[1] = cos(c[1]);

                    az = /* Page 30, Equation 5-4b */
                            atan2(
                            cosnlat * sin(dlon),
                            cosplat * sin(nlat) - sinplat * cosnlat * cos(dlon)
                            );

                    sinaz[0] = sin(az);
                    cosaz[0] = cos(az);
                    for (ii = 1; ii < 8; ii++) {
                        az += PI / 4;
                        sinaz[ii] = sin(az);
                        cosaz[ii] = cos(az);
                    }
                }

                for (ii = 0; ii < 8; ii++) {
                    double lat, lon;
                    int ec; /* edge = 0, corner = 1 */

                    ec = ii & 1;

                    /* Page 31, Equation 5-5 */
                    lat = asin(sinplat * cosc[ec]
                            + cosplat * sinc[ec] * cosaz[ii]);

                    lon = plon /* Page 31, Equation 5-6 */
                            + atan2(sinc[ec] * sinaz[ii],
                            cosplat * cosc[ec]
                            - sinplat * sinc[ec] * cosaz[ii]);

                    if (wrad > erad) { /* 180-degree meridian included */
                        if (lon >= wrad && lon > erad) {
                            corners[ii].x = (lon - wrad) * scale;
                        } else if (lon < wrad && lon <= erad) {
                            corners[ii].x = (lon - wrad + 2 * PI) * scale;
                        } else if (lon - erad < wrad - lon) {
                            corners[ii].x = (lon - wrad + 2 * PI) * scale;
                        } else {
                            corners[ii].x = (lon - wrad) * scale;
                        }
                    } else { /* 180-degree meridian not included */
                        corners[ii].x = (lon - wrad) * scale;
                    }
                    corners[ii].y = (nrad - lat) * scale;

                    if (corners[ii].x < 0 || corners[ii].x >= width
                            || corners[ii].y < 0 || corners[ii].y >= height) out++;
                }

                /*
                If out == 8, then the entire pixel is outside the
                image area, so I skip to the next one.
                 */
                if (out == 8) continue;

                l2_flags = l2_str[k].l2_flags[j];
                if (qualidx >= 0) {
                    quality = l2_str[k].l2_data[qualidx][j];
                }


                /* Use scan conversion to determine
                 *  which bins are overlain by this pixel. */
                numoverlain = 0;
                scan_convert(corners);

                for (m = 0; m < numoverlain; m++) {
                    n = binsoverlain[m];

                    /*
                    Keep track of the masks if so requested.
                    I am basing the mask logic on the code in put_l2brs.c
                    which is part of the level-2 browse file generator (l2brsgen).
                    At the time I write this, pixel values 251 and 255 are unused
                    in the SeaWiFS level-2 browse files.	Norman Kuring	10-Jul-2001
                     */

                    if (input.mask && qualidx < 0) {
                        if (l2_flags & mask_252 && mask_glint &&
                                strcmp(input.prod, "sst") != 0 &&
                                strcmp(input.prod, "sst4") != 0) {
                            mask[n] = 252;
                        } else if (l2_flags & mask_253) {
                            if (mask[n] != 252) mask[n] = 253;
                        } else if (l2_flags & all_masks) {
                            if (mask[n] != 252 && mask[n] != 253) mask[n] = 254;
                        }
                    } else if (input.mask && qualidx >= 0) {
                        if (l2_flags & mask_253) {
                            mask[n] = 253;
                        } else if (quality > input.quality) {
                            if (mask[n] != 253) mask[n] = 254;
                        }
                    }

                    /* Accumulate the sums for each bin. Only count unmasked pixels. */
                    goodpix = 1;
                    if (qualidx >= 0 && quality > input.quality) goodpix = 0;
                    else if (qualidx < 0 && (l2_flags & all_masks) != 0) goodpix = 0;

                    if (goodpix) {
                        val = l2_str[k].l2_data[prodidx][j];
                        if (val > input.datamax) val = input.datamax;
                        if (val < input.datamin) val = input.datamin;
                        sum[n] += val;
                        count[n]++;
                    }
                }
            }
            if (i != 0 && i % progress == 0) fprintf(stderr, ".");
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "Finished Plate Carree projection mapping\n");


    /* Calculate the means for each bin. */
    fprintf(stderr, "Computing means...\n");
    numoutpix = 0;
    for (n = 0; n < numbins; n++) {
        if (count[n] == 0) {
            if (input.mask) {
                if (mask[n] > 0) {
                    mean[n] = mask[n];
                    numoutpix++;
                } else {
                    mean[n] = 255;
                }
            }
        } else {
            if (strcmp(scale_type, "LINEAR") == 0) {
                mean[n] = (sum[n] / count[n] - intercept) / slope;
            } else if (strcmp(scale_type, "LOG") == 0) {
                mean[n] = (log10(sum[n] / count[n]) - intercept) / slope;
            }
            numoutpix++;
        }
    }

    if ((double) 100 * numoutpix / numbins < threshold) {
        fprintf(stderr, "Number of output pixels (%u of a possible %u) ",
                numoutpix, numbins);
        fprintf(stderr, "< threshold (%.2f%%). Output image not generated.\n",
                threshold);
        exit(BINBELOWTHRESH);
    }

    if (input.outmode == 1) { // ppm

        // setup output file pointer
        if (strlen(input.ofile) > 0) {
            fprintf(stderr, "Writing out file: %s ...\n", input.ofile);
            outfp = fopen(input.ofile, "wb");
        } else {
            fprintf(stderr, "Writing out data to STDOUT...\n");
            outfp = stdout;
        }

        if (input.apply_pal == 0 && default_palfile == 1) {
            fprintf(outfp, "P5\n%d %d\n255\n", width, height);
            for (n = 0; n < numbins; n++) {
                fputc((int) mean[n], outfp);
            }
        } else {
            /*
             * Convert from greyscale and palette to rgb array
             */
            for (i = 0; i < numbins; i++) {
                rgb[3 * i] = r[mean[i]];
                rgb[3 * i + 1] = g[mean[i]];
                rgb[3 * i + 2] = b[mean[i]];
            }

            fprintf(outfp, "P6\n%d %d\n255\n", width, height);
            fwrite(&rgb[0], 1, width * height * 3, outfp);
        }

        // close file if not using stdout
        if (stdout != outfp)
            fclose(outfp);

        // end ppm

    } else if (input.outmode == 2) { // png

        // setup output file pointer
        if (strlen(input.ofile) > 0) {
            fprintf(stderr, "Writing out file: %s ...\n", input.ofile);
            outfp = fopen(input.ofile, "wb");
        } else {
            fprintf(stderr, "Writing out data to STDOUT...\n");
            outfp = stdout;
        }

        png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                NULL, NULL, NULL);
        if (!png_ptr) {
            fprintf(stderr, "l2mapgen: Error: Unable to create PNG write structure.\n");
            exit(EXIT_FAILURE);
        }

        png_infop info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr) {
            png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
            fprintf(stderr, "l2mapgen: Error: Unable to create PNG info structure.\n");
            exit(EXIT_FAILURE);
        }
        if (setjmp(png_jmpbuf(png_ptr))) {
            png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
            fprintf(stderr, "l2mapgen: Error: Unable to call PNG setjmp().\n");
            exit(EXIT_FAILURE);
        }
        png_init_io(png_ptr, outfp);

        uint8_t * row_pointers[height];
        for (i = 0; i < height; i++)
            row_pointers[i] = mean + (i * width);
        png_set_rows(png_ptr, info_ptr, row_pointers);

        if (input.apply_pal == 0 && default_palfile == 1) {

            // Grayscale
            png_set_IHDR(png_ptr, info_ptr, width, height,
                    8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
                    PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
        } else {

            // color palette
            png_set_IHDR(png_ptr, info_ptr, width, height,
                    8, PNG_COLOR_TYPE_PALETTE, PNG_INTERLACE_NONE,
                    PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
            png_set_PLTE(png_ptr, info_ptr, (png_const_colorp) input.palette, 256);
        }

        png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

        /* clean up after the write, and free any memory allocated */
        png_destroy_write_struct(&png_ptr, (png_infopp) NULL);

        // close file if not using stdout
        if (stdout != outfp)
            fclose(outfp);

        // end png

    } else if (input.outmode == 3) { // geotiff
        uint16_t *rr, *gg, *bb;
        TIFF *tiff;
        GTIF *gtif;
        float *fmean;

        if (strlen(input.ofile) == 0) {
            fprintf(stderr, "Can not write TIFF file to stdout.\n");
            exit(EXIT_FAILURE);
        }

        tiff = XTIFFOpen(input.ofile, "w");
        if (tiff == NULL) {
            fprintf(stderr, "Could not open outgoing image\n");
            exit(EXIT_FAILURE);
        }
        gtif = GTIFNew(tiff);
        if (gtif == NULL) {
            fprintf(stderr, "Could not create geoTIFF structure\n");
            exit(EXIT_FAILURE);
        }

        // calc geo TIFF  tags
        double tiepoints[6] = {0, 0, 0, 0, 0, 0};
        double pixscale[3] = {0, 0, 0};

        // pixel width
        pixscale[0] = (input.east - input.west) / width;

        // pixel height
        pixscale[1] = (input.north - input.south) / height;

        // set top left corner pixel lat, lon
        tiepoints[3] = input.west;
        tiepoints[4] = input.north;

        TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, width);
        TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, height);
        TIFFSetField(tiff, GTIFF_PIXELSCALE, 3, pixscale);
        TIFFSetField(tiff, GTIFF_TIEPOINTS, 6, tiepoints);

        if (input.apply_pal == 0 && default_palfile == 1) {
            // Grayscale

            TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
            TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
            TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 32);
            TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 1);
            TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);

            // write geo TIFF keys
            GTIFKeySet(gtif, GTModelTypeGeoKey, TYPE_SHORT, 1, ModelGeographic);
            GTIFKeySet(gtif, GTRasterTypeGeoKey, TYPE_SHORT, 1, RasterPixelIsArea);
            GTIFKeySet(gtif, GeographicTypeGeoKey, TYPE_SHORT, 1, GCS_WGS_84);

            GTIFWriteKeys(gtif);

            fmean = (float*) malloc(numbins * sizeof (float));
            if (fmean == NULL) {
                fprintf(stderr, "Could not allocate memory for TIFF grayscale mean\n");
                exit(EXIT_FAILURE);
            }

            // calc the mean as a float
            for (i = 0; i < numbins; i++) {
                if (count[i] > 0) {
                    fmean[i] = sum[i] / count[i];
                } else {
                    fmean[i] = -32767.0;
                }
            }

            // Actually write the image
            if (TIFFWriteEncodedStrip(tiff, 0, fmean, numbins * sizeof (float)) == 0) {
                fprintf(stderr, "Could not write TIFF image\n");
                exit(EXIT_FAILURE);
            }

            free(fmean);

        } else {
            // Colormap
            rr = (uint16_t*) malloc(256 * sizeof (uint16_t));
            gg = (uint16_t*) malloc(256 * sizeof (uint16_t));
            bb = (uint16_t*) malloc(256 * sizeof (uint16_t));
            if (rr == NULL || gg == NULL || bb == NULL) {
                fprintf(stderr, "Could not allocate memory for TIFF color map\n");
                exit(EXIT_FAILURE);
            }

            // need a colormap of shorts not bytes
            for (i = 0; i < 256; i++) {
                rr[i] = r[i] << 8;
                gg[i] = g[i] << 8;
                bb[i] = b[i] << 8;
            }

            TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
            TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_PALETTE);
            TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 8);
            TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 1);
            TIFFSetField(tiff, TIFFTAG_COLORMAP, rr, gg, bb);

            // write geo TIFF keys
            GTIFKeySet(gtif, GTModelTypeGeoKey, TYPE_SHORT, 1, ModelGeographic);
            GTIFKeySet(gtif, GTRasterTypeGeoKey, TYPE_SHORT, 1, RasterPixelIsArea);
            GTIFKeySet(gtif, GeographicTypeGeoKey, TYPE_SHORT, 1, GCS_WGS_84);

            GTIFWriteKeys(gtif);

            // Actually write the image
            if (TIFFWriteEncodedStrip(tiff, 0, mean, width * height) == 0) {
                fprintf(stderr, "Could not write TIFF image\n");
                exit(EXIT_FAILURE);
            }

            free(rr);
            free(gg);
            free(bb);

        } // color map

        GTIFFree(gtif);
        XTIFFClose(tiff);

        // end geotiff

    } else {

        fprintf(stderr, "l2mapgen: Error: invalid value for outmode - %d.\n", input.outmode);
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < l3m_params; i++) {
        free(parmname_list[i]);
        free(parmname_short[i]);
        free(unit_list[i]);
        free(scaling_list[i]);
        free(palette_list[i]);
        free(precision_list[i]);
    }
    free(parmname_list);
    free(parmname_short);
    free(unit_list);
    free(scaling_list);
    free(palette_list);
    free(maximum_list);
    free(minimum_list);
    free(precision_list);

    free(sum);
    free(count);
    free(mean);
    free(mask);
    free(rgb);
    free(r);
    free(g);
    free(b);
    fprintf(stderr, "Done.\n");

    return (EXIT_SUCCESS);

}



/*------------------------------------------------------------------------------

   Function: scan_convert

   Returns type: int

   Description:	This function is derived from the scan-conversion
                code used by an X-server for filling in polygons.

   Parameters: (in calling order)
      Type		Name		I/O	Description
      ----		----		---	-----------
      struct _DDXPoint*	ptsIn		In	vertex specifications

   Modification history:
      Programmer	Date		Description of change
      ----------	----		---------------------
      Norman Kuring	15-Sep-1992	Adaptation from X-server code
                                        for my own nefarious purposes.

------------------------------------------------------------------------------*/

/***********************************************************
Copyright 1987 by Digital Equipment Corporation, Maynard, Massachusetts,
and the Massachusetts Institute of Technology, Cambridge, Massachusetts.

                        All Rights Reserved

Permission to use, copy, modify, and distribute this software and its 
documentation for any purpose and without fee is hereby granted, 
provided that the above copyright notice appear in all copies and that
both that copyright notice and this permission notice appear in 
supporting documentation, and that the names of Digital or MIT not be
used in advertising or publicity pertaining to distribution of the
software without specific, written prior permission.  

DIGITAL DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING
ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL
DIGITAL BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR
ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
SOFTWARE.

 ******************************************************************/
/* $XConsortium: mipolygen.c,v 5.0 89/06/09 15:08:35 keith Exp $ */

extern void miloadAET(), micomputeWAET(), miFreeStorage();

/*
 *
 *     Written by Brian Kelleher;  Oct. 1985
 *
 *     Routine to fill a polygon.  Two fill rules are
 *     supported: frWINDING and frEVENODD.
 *
 *     See fillpoly.h for a complete description of the algorithm.
 *
 *     Modified by Norman Kuring to fill a bin-number array instead of
 *     filling shapes on the server.
 */

int scan_convert(XPoint *ptsIn) {
    int count = 8;
    register EdgeTableEntry *pAET; /* the Active Edge Table   */
    register int y; /* the current scanline    */
    register int nPts = 0; /* number of pts in buffer */
    register EdgeTableEntry *pWETE; /* Winding Edge Table      */
    register ScanLineList *pSLL; /* Current ScanLineList    */
    register XPoint *ptsOut; /* ptr to output buffers   */
    int *wdth;
    XPoint FirstPoint[NUMPTSTOBUFFER]; /* the output buffers */
    int FirstWidth[NUMPTSTOBUFFER];
    EdgeTableEntry *pPrevAET; /* previous AET entry      */
    EdgeTable ET; /* Edge Table header node  */
    EdgeTableEntry AET; /* Active ET header node   */
    EdgeTableEntry *pETEs; /* Edge Table Entries buff */
    ScanLineListBlock SLLBlock; /* header for ScanLineList */
    int fixWAET = 0;

    if (!(pETEs = (EdgeTableEntry *) malloc(sizeof (EdgeTableEntry) * count))) {
        return (0);
    }
    ptsOut = FirstPoint;
    wdth = FirstWidth;
    if (!miCreateETandAET(count, ptsIn, &ET, &AET, pETEs, &SLLBlock)) {
        free(pETEs);
        return (0);
    }
    pSLL = ET.scanlines.next;

    /*
     *  for each scanline
     */
    for (y = ET.ymin; y < ET.ymax; y++) {
        /*
         *  Add a new edge to the active edge table when we
         *  get to the next edge.
         */
        if (pSLL && y == pSLL->scanline) {
            miloadAET(&AET, pSLL->edgelist);
            micomputeWAET(&AET);
            pSLL = pSLL->next;
        }
        pPrevAET = &AET;
        pAET = AET.next;
        pWETE = pAET;

        /*
         *  for each active edge
         */
        while (pAET) {
            /*
             *  if the next edge in the active edge table is
             *  also the next edge in the winding active edge
             *  table.
             */
            if (pWETE == pAET) {
                ptsOut->x = pAET->bres.minor;
                ptsOut++->y = y;
                *wdth++ = pAET->nextWETE->bres.minor - pAET->bres.minor;
                nPts++;

                /*
                 *  send out the buffer
                 */
                if (nPts == NUMPTSTOBUFFER) {
                    collect_bins(nPts, FirstPoint, FirstWidth);
                    ptsOut = FirstPoint;
                    wdth = FirstWidth;
                    nPts = 0;
                }

                pWETE = pWETE->nextWETE;
                while (pWETE != pAET)
                    EVALUATEEDGEWINDING(pAET, pPrevAET, y, fixWAET);
                pWETE = pWETE->nextWETE;
            }
            EVALUATEEDGEWINDING(pAET, pPrevAET, y, fixWAET);
        }

        /*
         *  reevaluate the Winding active edge table if we
         *  just had to resort it or if we just exited an edge.
         */
        if (miInsertionSort(&AET) || fixWAET) {
            micomputeWAET(&AET);
            fixWAET = 0;
        }
    }

    /*
     *     Get any spans that we missed by buffering
     */
    collect_bins(nPts, FirstPoint, FirstWidth);
    free(pETEs);
    miFreeStorage(SLLBlock.next);
    return (1);
}

int collect_bins(
        int number_of_initial_points,
        XPoint *initial_point,
        int *span_width
        ) {

    int i;
    short x, y;
    int end;

    /* The variables, width, height, binsoverlain, and numoverlain, are
     *  static global varibles. */

    for (i = 0; i < number_of_initial_points; i++) {

        y = initial_point[i].y;
        if (y >= height || y < 0) continue;

        for (
                x = initial_point[i].x,
                end = initial_point[i].x + span_width[i];
                x < end;
                x++
                ) {
            if (x >= width || x < 0) continue;
            if (numoverlain >= binsperpixel) {
                binsperpixel += 1024;
                binsoverlain = (int *) realloc(binsoverlain, binsperpixel * sizeof (int));
                if (binsoverlain == NULL) {
                    fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
                    fprintf(stderr,
                            "Memory allocation error (numoverlain=%d binsperpixel=%d\n",
                            numoverlain, binsperpixel);
                    exit(EXIT_FAILURE);
                }
                return (0);
            }
            binsoverlain[numoverlain++] = y * width + x;
        }
    }
    return (1);
}
