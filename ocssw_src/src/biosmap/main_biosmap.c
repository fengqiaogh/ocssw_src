/* ========================================================================
 * SWbiosmap - generates a biosphere map from a set of SeaWiFS SMI NDVI and
 *            chlor_a maps.
 * 
 * Synopsis:
 *
 *   SWbiosmap -l input-filename-list | (-c chlorophyll-file -n ndvi-file)
 *            [-o output-filename] [-r replacement-filename]
 *
 *   output-filename      : output map file, def="map.hdf"
 *   chlorophyll-file     : name of input Chlorophyll SMI file
 *   ndvi-file-file       : name of input NDVI SMI file
 *   input-filename-list  : file containing list of input SMI filenames
 *   replacement-filename : name of file which this new output supercedes
 *
 * Description:
 * 
 *   SWbiosmap reads through a set of SeaWiFS SMI NDVI and Chlor-a files
 *   files, and outputs a merged SMI map file containing NDVI over land and 
 *   Chl_a over ocean. Multiple observations of the same grid point are 
 *   resolved by "last value in". The input is the name of a file which 
 *   contains a <CR> separated list of NDVI and/or Chl_a SMI files. If only
 *   one Chlorophyll file and one NDVI file are to be combined, the -c and
 *   -n options can be used instead of -l.  The -r option allows 
 *   specification of a filename which the new output product is intended 
 *   to replace.
 *
 * Modification history:
 *
 *     Programmer     Organization      Date      Description of change
 *   --------------   ------------    --------    ---------------------
 *   Bryan A. Franz   GSC             09/12/97    Original development
 *
 * ======================================================================== */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <libgen.h>

typedef unsigned char byte;
#include "map.h"
#include "seaproto.h"
#include "meta_l3b.h"
#include "mapproto_o.h"
#include <timeutils.h>
#include "mapproto.h"

#define VERSION    "2.0"
#define NX          4320                 /* Map grid x-dimension           */
#define NY          2160                 /* Map grid y-dimension           */
#define FILENAMELEN 255                  /* Max filename length            */
#define MAXSTRNGLEN 255                  /* Max generic string length      */
#define MAXINFILES  50                   /* Max number of input files      */
#define NODATA8     255                  /* Missing data value (byte)      */
#define NODATA16    65535                /* Missing data value (short int) */
#define CMD_ARGS    "o:l:r:c:n:"         /* Valid commandline options      */

/* Macro definitions */
#define max(A,B)    ((A) > (B) ? (A) : (B))  /* Returns greater of A and B */ 
#define min(A,B)    ((A) < (B) ? (A) : (B))  /* Returns lesser of A and B  */

/* Function prototypes */
void usage(char* progname);
double ydmsec2jul(int16 year, int16 day, int32 msec);
void jul2ydmsec(double jul, int16 *year, int16 *day, int32 *msec);

/* ------------------------------------------------------------------------ *
 *                              main                                        *
 * ------------------------------------------------------------------------ */
int main(int argc, char *argv[]) {

    static char infileList[MAXINFILES][FILENAMELEN]; /* Input file list      */
    static uint8 map08[NY][NX]; /* Equiangular grid of map data          */
    static uint16 map16[NY][NX]; /* Equiangular grid of map data          */
    int32 ix; /* Longitudinal bin index                */
    int32 iy; /* Latitudinal bin index                 */
    FILE *fp; /* Pointer to list file                  */
    int16 nfiles = 0; /* Number of input files in infileList   */
    int16 iarg; /* Command-line argument index          */
    int16 ifile; /* Input file index                     */
    double jstart = 9e7; /* Julian start time of mapped data     */
    double jstop = -9e7; /* Julian stop time of mapped data      */
    double sjul; /* Julian start time of input file      */
    double ejul; /* Julian Stop time of input file       */
    double deljul; /* Julian time difference               */
    int16 syear;
    int16 sday;
    int32 smsec;
    int16 eyear;
    int16 eday;
    int32 emsec;
    int32 status;
    uint16 offset;

    /* Parameters for getopt() */
    extern int opterr;
    extern int optind;
    extern char *optarg;
    int c;

    /* Parameters for L3 map interface */
    char outfile[FILENAMELEN] = "map.hdf";
    char replaces[MAXSTRNGLEN] = "ORIGINAL";
    int16 bin_syear;
    int16 bin_sday;
    int16 bin_eyear;
    int16 bin_eday;
    int16 map_syear;
    int16 map_sday;
    int32 map_smsec;
    int16 map_eyear;
    int16 map_eday;
    int32 map_emsec;
    float map_latrange[2] = {90., -90.};
    float map_lonrange[2] = {-180., 180.};
    int32 lines = NY;
    int32 columns = NX;
    char flag_use[MAXSTRNGLEN] = "";
    byte eng_q_use[4] = {0, 0, 0, 0};
    char ptime[17];
    static char infiles[(FILENAMELEN + 1) * MAXINFILES];
    char prod_type[MAXSTRNGLEN];
    int32 nbins = 0;
    char l3m_name[MAXSTRNGLEN];
    char measure[MAXSTRNGLEN] = "mean";
    char l3_flag_names[MAXSTRNGLEN];
    char proc_con[MAXSTRNGLEN];
    char proc_log[MAXSTRNGLEN];

    char precision[2];
    void *map;
    int32 sdfid;
    int32 attr_index;
    char buffer[1024];
    char *cptr;
    char last_prec[2];

    meta_l3bType meta_l3b;

    /* Additional parameters for l3m_read() */
    static byte image08[NY][NX];
    static uint16 image16[NY][NX];
    static byte palette[3 * 256];
    void *image;
    meta_struct meta_l3m;


    /* Log command-line sequence to process-control string */
    strcpy(proc_con, argv[0]);
    for (iarg = 1; iarg < argc; iarg++) {
        strcat(proc_con, " ");
        strcat(proc_con, argv[iarg]);
    }

    /* Process command-line arguments */
    while ((c = getopt(argc, argv, CMD_ARGS)) != EOF) {
        switch (c) {
        case 'c': /* Input chlorophyll file */
        case 'n': /* Input NDVI file */
            strncpy(&infileList[nfiles][0], optarg, FILENAMELEN);
            nfiles++;
            break;
        case 'o': /* Output biosphere file */
            strncpy(outfile, optarg, FILENAMELEN);
            break;
        case 'l': /* Input file list */
            if ((fp = fopen(optarg, "r")) == NULL) {
                fprintf(stderr, "%s: Error opening %s for reading.\n",
                        argv[0], optarg);
                exit(1);
            }
            while (fscanf(fp, "%s", &infileList[nfiles][0]) != EOF)
                nfiles++;
            break;
        case 'r': /* File to replace at DACC */
            strncpy(replaces, optarg, FILENAMELEN);
            break;
        default:
            printf("Invalid argument: %d\n", c);
            usage(argv[0]);
            exit(1);
            break;
        }
    }
    if (nfiles == 0) {
        printf("No input files specified.\n");
        usage(argv[0]);
        exit(1);
    }


    /* Get process start time */
    get_time(ptime);

    strcpy(precision, "B");
    for (ifile = 0; ifile < nfiles; ifile++) {
        sdfid = SDstart(infileList[ifile], DFACC_RDONLY);
        if (sdfid != -1) {
            attr_index = SDfindattr(sdfid, "Input Parameters");
            status = SDreadattr(sdfid, attr_index, buffer);
            cptr = strstr(buffer, "PRECISION=");
            if (cptr != NULL) {
                cptr += 10;
                if (*cptr == 'I') strcpy(precision, "I");
                if (ifile > 0 && strcmp(last_prec, precision) != 0) {
                    printf("Input files are of varying precision\n");
                    exit(2);
                }
                strncpy(last_prec, precision, 1);
            }
        }
    }

    /* Initialize output map array to no data */
    if (precision[0] == 'I') {
        map = map16;
        image = image16;
        for (iy = 0; iy < NY; iy++)
            for (ix = 0; ix < NX; ix++)
                map16[iy][ix] = NODATA16;
    } else {
        map = map08;
        image = image08;
        for (iy = 0; iy < NY; iy++)
            for (ix = 0; ix < NX; ix++)
                map08[iy][ix] = NODATA8;
    }

    /*                                                      */
    /* Loop through each input file and process accordingly */
    /* ---------------------------------------------------- */
    /*                                                      */
    for (ifile = 0; ifile < nfiles; ifile++) {

        printf("Processing L3 Map file: %s\n", infileList[ifile]);

        status = get_l3m(infileList[ifile],
                &syear, &sday, &smsec,
                &eyear, &eday, &emsec,
                prod_type, l3m_name,
                image, palette, &meta_l3m);

        printf("syear = %d, sday = %d, smsec = %d\n", syear, sday, smsec);
        printf("eyear = %d, eday = %d, emsec = %d\n", eyear, eday, emsec);

        printf("name = %s\n", l3m_name);
        if (strncmp(l3m_name, "Chlorophyll", 11) == 0)
            offset = 0;
        else
            if (precision[0] == 'I') offset = 32768 - 1;
        else offset = 128;


        /* Update data time range */
        sjul = ydmsec2jul(syear, sday, smsec);
        ejul = ydmsec2jul(eyear, eday, emsec);
        jstart = min(jstart, sjul);
        jstop = max(jstop, ejul);


        /* 
         * Loop through each bin of the input image.  If bin has valid
         * data, store value in merged map.
         */

        if (precision[0] == 'I') {
            for (iy = 0; iy < NY; iy++) {
                for (ix = 0; ix < NX; ix++) {
                    if (image16[iy][ix] != NODATA16)
                        map16[iy][ix] = image16[iy][ix] / 2 + offset;

                }
            }
        } else {
            for (iy = 0; iy < NY; iy++) {
                for (ix = 0; ix < NX; ix++) {
                    if (image08[iy][ix] != NODATA8)
                        map08[iy][ix] = image08[iy][ix] / 2 + offset;

                }
            }
        }


    } /* End loop over input files */


    /*                                          */
    /* Create remaining meta data for Map file  */
    /* ---------------------------------------- */
    /*                                          */

    strcpy(l3m_name, "biosphere");
    strcpy(l3_flag_names, "");
    strcpy(flag_use, "");
    strcpy(proc_log, "");

    jul2ydmsec(jstart, &bin_syear, &bin_sday, &smsec);
    jul2ydmsec(jstop + 1, &bin_eyear, &bin_eday, &emsec);
    jul2ydmsec(jstart, &map_syear, &map_sday, &map_smsec);
    jul2ydmsec(jstop, &map_eyear, &map_eday, &map_emsec);

    /* Create file input list */
    strcpy(infiles, basename(infileList[0]));
    for (ifile = 0; ifile < nfiles; ifile++) {
        strcat(infiles, ", ");
        strcat(infiles, basename(infileList[ifile]));
    }

    /* Compute number of filled bins */
    if (precision[0] == 'I') {
        for (iy = 0; iy < NY; iy++)
            for (ix = 0; ix < NX; ix++)
                if (map16[iy][ix] != NODATA16)
                    nbins++;
    } else {
        for (iy = 0; iy < NY; iy++)
            for (ix = 0; ix < NX; ix++)
                if (map08[iy][ix] != NODATA8)
                    nbins++;
    }

    /* Identify product binning period */
    deljul = jstop - jstart;
    if (deljul <= 0.1)
        strcpy(prod_type, "scene");
    else if (deljul > 0.1 && deljul <= 1.0)
        strcpy(prod_type, "day");
    else if (deljul > 1.0 && deljul <= 8.0)
        strcpy(prod_type, "8-day");
    else if (deljul > 8.0 && deljul <= 16.0)
        strcpy(prod_type, "16-day");
    else if (deljul > 16.0 && deljul <= 31.0)
        strcpy(prod_type, "month");
    else if (deljul > 31.0 && deljul <= 367.0)
        strcpy(prod_type, "year");
    else
        strcpy(prod_type, "multi-year");

    /* Set the L3_bin metadata block to zeros */
    memset(&meta_l3b, '\0', sizeof (meta_l3bType));


    /*                       */
    /* Write output map file */
    /* --------------------- */
    /*                       */
    printf("Writing L3 Map file: %s\n", outfile);
    printf("syear = %d, sday = %d, smsec = %d\n", map_syear, map_sday, map_smsec);
    printf("eyear = %d, eday = %d, emsec = %d\n", map_eyear, map_eday, map_emsec);
    printf("Number of filled bins: %d of %d\n", nbins, NX * NY);

    /*
    To avoid changes the calling sequence to "put_l3m", we pass
    info about the precision through the "nbins" parameter.
     */
    if (precision[0] == 'I') nbins = -nbins;

    status = put_l3m(outfile, replaces,
            bin_syear, bin_sday, bin_eyear, bin_eday,
            map_syear, map_sday, map_smsec,
            map_eyear, map_eday, map_emsec,
            map_latrange, map_lonrange, lines, columns,
            l3_flag_names, flag_use, eng_q_use, ptime, infiles,
            prod_type, nbins, l3m_name, map, measure,
            proc_con, proc_log, &meta_l3b);
    if (status < 0) {
        printf("%s: Error opening or writing %s as L3 SMI file\n",
                argv[0], outfile);
        exit(1);
    }

    return 0;
}

/* ------------------------------------------------------------------------- *
 * usage - display proper calling sequence for this program                  *
 * ------------------------------------------------------------------------- */
void usage(char *progname) {

    printf("%s %s (%s %s)\n", progname, VERSION, __DATE__, __TIME__);

    fprintf(stderr, "\nUsage: %s ", progname);
    fprintf(stderr, "-l input-filename-list | ");
    fprintf(stderr, "(-c chlorophyll-file -n ndvi-file)\n");
    fprintf(stderr, "                [-o output-filename] ");
    fprintf(stderr, "[-r replacement-filename]\n");
    fprintf(stderr, "\n output-filename      : output map file, def=\"map.hdf\"\n");
    fprintf(stderr, " input-filename-list  : file containing list of "
            "L3 SMI CHLO or NDVI filenames\n");
    fprintf(stderr, " chlorophyll-file     : input L3 SMI CHLO filename\n");
    fprintf(stderr, " input-filename       : input L3 SMI NDVI filename\n");
    fprintf(stderr, " replacement-filename : filename which this product is intended to replace\n");
}

/* ------------------------------------------------------------------------- *
 * ydmsec2jul - converts year, day-of-year, millisecs-of-day to Julian       * 
 * ------------------------------------------------------------------------- */
double ydmsec2jul(int16 year, int16 day, int32 msec) {
    double jul;

    jul = (367 * year - (7 * year) / 4 + day + 1721044) + msec / 8.64e7;

    return ( jul);
}

/* ------------------------------------------------------------------------- *
 * jul2ydmsec - converts Julian to year, day-of-year, millisecs-of-day       *
 * ------------------------------------------------------------------------- */
void jul2ydmsec(double jul, int16 *year, int16 *day, int32 *msec) {
    int32 days_since;
    int32 years_since;

    /* Compute days since January 0, 1900 */
    days_since = (int32) jul - 2415020;

    /* Compute years since 1900 */
    years_since = 4 * days_since / 1461;

    /* Compute year, day-of-year, msecs of day */
    *year = years_since + 1900;
    *day = days_since - 1461 * (years_since - 1) / 4 - 365;
    *msec = (int32) (fmod(jul, 1.0)*8.64e7);
}


