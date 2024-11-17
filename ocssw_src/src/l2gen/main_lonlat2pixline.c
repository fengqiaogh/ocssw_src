/* ========================================================================
 * MSscanpixlimits - Determines bounding scans and pixels for given lon/lat box
 *
 * Synopsis:
 *
 *   MSscanpixlimits input-filename [geo file] SW-lon SW-lat NE-lon NE-lat
 *
 * Description:
 *
 * Modification history:
 *
 *     Programmer     Organization      Date       Description of change
 *   --------------   ------------    --------     ---------------------
 *   Joel M. Gales    Futuretech      20 Sept 1999 Original Development
 *   Joel M. Gales    Futuretech      29 Sept 2000 Perform single pixel
 *                                                 search within initial
 *                                                 scan loop.
 *   Joel M. Gales    Futuretech      03 Oct  2000 Change uv, uvpix, maxcos
 *                                                 and dot to float64
 *
 *   Joel M. Gales    Futuretech      08 Nov  2006 Processes both MODIS &
 *                                                 non-MODIS L1 & L2 granules
 *                                                 (OCTS, CZCS, SEAWIFS)
 *
 *   Joel M. Gales    Futuretech      08 Feb  2007 Remove obsolete argc=3,5
 *                                                 cases.
 *                                                 Update usage statement
 *                                                 Set ver# to 6.00
 *
 *   Joel M. Gales    Futuretech      12 Feb  2007 Put back argc=3,5
 *                                                 cases.
 *                                                 Update usage statement
 *                                                 Set ver# to 6.01
 *
 *   Joel M. Gales    Futuretech      22 Feb  2007 Subtract extract pixel
 *                                                 offset if present
 *                                                 Set ver# to 6.03
 *
 *   Joel M. Gales    Futuretech      13 Mar  2007 Subtract extract pixel
 *                                                 offset for pixsrch also
 *                                                 Set ver# to 6.04
 *
 *   Joel M. Gales    Futuretech      15 Mar  2007 Move executable code to
 *                                                 after declaration statements
 *
 *   Joel M. Gales    Futuretech      16 May  2007 Reset min scan & pix to -1
 *                                                 if region not found
 *
 *   Joel M. Gales    Futuretech      15 Oct  2007 Add support for box
 *                                                 extraction
 * * ======================================================================== */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "lonlat2pixline.h"

#define CMD_ARGS    "+x:y:r:o:Fv"   /* Valid commandline options */
#define VERSION_SPL "6.5.2"

void usage(const char *prog) {
    printf("%s %s (%s %s)\n",
            "lonlat2pixline", VERSION_SPL, __DATE__, __TIME__);
    printf("Usage (box): %s [-x box-halfwidth] [-y box-halfheight] [-r res]\n", prog);
    printf("             [-v] infile [geofile] SWlon SWlat NElon NElat\n");
    printf("Usage (pix): %s [-x box-halfwidth] [-y box-halfheight] [-r res]\n", prog);
    printf("             [-v] [-F] infile [geofile] lon lat\n");
    printf("where:\n");
    printf("    infile  is either a L1B, L2, or GEO filename\n");
    printf("    geofile is the GEO filename used for MODIS L1B or VIIRS L1 files\n");
    printf("    -r res  resolution for MODIS Files (1000,500,250)\n");
    printf("    -x val  min number of pixels on either side of location to include\n");
    printf("    -y val  min number of scan lines on top and bottom of location to include\n");
    printf("    -F      return 110 if full-box can't be extracted and\n");
    printf("              120 if full file would be extracted\n");
    printf("    -v      print more verbose messages\n");
    printf("    -o file output to file instead of stdout\n");
    printf("    return    0   everything OK\n");
    printf("            100   location not found in scene\n");
    printf("            110   full box not extracted using the -F option\n");
    printf("            120   the whole file is inside the given coordinates, using -F\n");
    printf("            other error number\n\n");

    exit(FATAL_ERROR);
}

/* -------------------------------------------------------------------- *
 *                              main                                    *
 * -------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
    int result;
    int32_t want_FAIL = 0;
    char* ofileName = NULL;
    FILE* ofile = stdout;

    lonlat2pixline_t params;

    /* Parameters for getopt() */
    extern int opterr;
    extern int optind;
    extern char *optarg;
    int c;

    // cdata_();

    // set some normal defaults
    params.want_pixbox = 0;
    params.pix_srch = 0;
    params.resolution = -1;
    params.xbox = 0;
    params.ybox = 0;
    params.geo_filename[0] = 0;

    want_verbose = 0;

    while ((c = getopt(argc, argv, CMD_ARGS)) != EOF) {
        switch (c) {
        case 'x':
            params.xbox = atoi(optarg);
            if (params.xbox != 0 && params.ybox == 0)
                params.ybox = params.xbox;
            params.want_pixbox = 1;
            break;
        case 'y':
            params.ybox = atoi(optarg);
            if (params.xbox == 0 && params.ybox != 0)
                params.xbox = params.ybox;
            params.want_pixbox = 1;
            break;
        case 'r':
            params.resolution = atoi(optarg);
            break;
        case 'F':
            want_FAIL = 1;
            break;
        case 'v':
            want_verbose = 1;
            break;
        case 'o':
            ofileName = optarg;
            break;
        default:
            usage(argv[0]);
            break;
        }
    }

    if (ofileName != NULL) {
        ofile = fopen(ofileName, "w");
        if (ofile == NULL) {
            printf("-E- Can not open %s as output file\n", ofileName);
            exit(EXIT_FAILURE);
        }
    }

    /*                                                                  */
    /* Process required command-line arguments                          */
    /*                                                                  */
    switch (argc - optind + 1) {
    case 4:
        strcpy(params.input_filename, argv[1 + optind - 1]);
        params.SWlon = atof(argv[2 + optind - 1]);
        params.SWlat = atof(argv[3 + optind - 1]);
        params.pix_srch = 1;
        break;
    case 5:
        strcpy(params.input_filename, argv[1 + optind - 1]);
        strcpy(params.geo_filename, argv[2 + optind - 1]);
        params.SWlon = atof(argv[3 + optind - 1]);
        params.SWlat = atof(argv[4 + optind - 1]);
        params.pix_srch = 1;
        break;
    case 6:
        strcpy(params.input_filename, argv[1 + optind - 1]);
        params.SWlon = atof(argv[2 + optind - 1]);
        params.SWlat = atof(argv[3 + optind - 1]);
        params.NElon = atof(argv[4 + optind - 1]);
        params.NElat = atof(argv[5 + optind - 1]);
        if ((params.SWlon == params.NElon) && (params.SWlat == params.NElat)) {
            params.pix_srch = 1;
        } else {
            params.pix_srch = 0;
        }
        break;
    case 7:
        strcpy(params.input_filename, argv[1 + optind - 1]);
        strcpy(params.geo_filename, argv[2 + optind - 1]);
        params.SWlon = atof(argv[3 + optind - 1]);
        params.SWlat = atof(argv[4 + optind - 1]);
        params.NElon = atof(argv[5 + optind - 1]);
        params.NElat = atof(argv[6 + optind - 1]);
        if ((params.SWlon == params.NElon) && (params.SWlat == params.NElat)) {
            params.pix_srch = 1;
        } else {
            params.pix_srch = 0;
        }
        break;
    default:
        usage(argv[0]);
        break;
    }

    result = lonlat2pixline(&params);

    if ((result == 110 || result == 120) && !want_FAIL) {
        result = 0;
    }

    if (result && (result != 110) && (result != 120)) {
        exit(result);
    }

    if (params.pix_srch) {
        fprintf(ofile, "#Lon=%f\n", params.pixLon);
        fprintf(ofile, "#Lat=%f\n", params.pixLat);
    } else {
        fprintf(ofile, "# %d %d %d %d\n", params.spixl, params.epixl,
                params.sline, params.eline);
    }

    fprintf(ofile, "sline=%d\n", params.sline);
    fprintf(ofile, "eline=%d\n", params.eline);
    fprintf(ofile, "spixl=%d\n", params.spixl);
    fprintf(ofile, "epixl=%d\n", params.epixl);
    if (ofileName != NULL) {
        fprintf(ofile, "status=%d\n", result);
        fclose(ofile);
    }

    exit(result);
}
