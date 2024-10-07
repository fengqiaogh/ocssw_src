/*-----------------------------------------------------------------------------
    Program:   anc_seaice

    Description:  convert binary file of sea ice from wgrib2 into a
      stabndard file of sea ice

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       argc          I    count of command line args - 3
        char *[]  argv          I    command line arguments:
                                     see check_usage for all options

    Modification history:

    W. Robinson, SAIC   9 Dec 2009  initial routine
    W. Robinson, SAIC  14 Jan 2011  Add smoothed icemask

----------------------------------------------------------------------------*/
#include "ancil.h"
#include <time.h>
#include <limits.h>
#include <unistd.h>
#include "ancnrt_proto.h"

/*
 * Met data specific settings
 */

#define VGROUPNAME   "Geophysical Data"

#define VSIZE         0.0833
#define HSIZE         0.0833
#define MAX_NORTH    90.0
#define MAX_SOUTH   -90.0
#define MAX_WEST   -180.0
#define MAX_EAST    180.0
#define GRIBLATSZ    2160   /* GRIB file size in lat */
#define GRIBLONSZ    4320   /* GRIB file size in lon */
#define OUTLATSZ   2160    /* output # lines  */
#define OUTLONSZ   4320    /* output # pixels */
#define GRIB_MODE_1 1  /*  the GRIB modes, GRIB 1 old, GRIB 2 post jan 08 */
#define GRIB_MODE_2 2

/* #define VSIZE         2.5 */
/* #define HSIZE         2.5 */
/* #define MAX_EAST    177.5 */
/* #define WPHLATSZ     73 */ /* latitude  of GRIB files */
/* #define WPHLONSZ    144 */ /* long of MODIFIED GRIB files */

/* fortran prototype */
void mk_smooth_ice_map_(char* char_ice, float* frac_ice_smoothed);

int main(int argc, char *argv[]) {
    int i, j, k, l, anctyp, int_char_min, int_char_max;
    int rank, SDSref = 0;
    int result = 0;
    int array_size = 0;
    int numannarr = 0;
    int one = 1, rgmode = 1;
    int year = 0, month = 0, day = 0, hour = 0,
            julval = 0, pyear = 0;
    int32 shape[2], shape_smooth[2];

    float32 *datarr1, dat;

    static float32 latscl[GRIBLATSZ], lonscl[GRIBLONSZ], dlat, dlon;
    char message[MAXNAMELNG], inz_wind[MAXNAMELNG], inm_wind[MAXNAMELNG];
    char inpress[MAXNAMELNG], inp_water[MAXNAMELNG], inrel_hum[MAXNAMELNG];
    char in_oz[MAXNAMELNG], longname[MAXNAMELNG], *sdsname;
    char *units, *dataattr, *longnamelabel, *longnamevalue;
    char *unitslabel, *unitsvalue, *datafmt, outfile[MAXNAMELNG];
    char outfilename[MAXNAMELNG], outdir[MAXNAMELNG];
    char annotfile[MAXNAMELNG], source_name[100];
    int n_opt_arg;
    float latstep, lonstep; /* lat/lon incr. */
    int latsz, lonsz, np, nl, fcst; /* lat/lon size and fcst info holders */
    float nmostlat, smostlat, wmostlon, emostlon; /* lat/lon corners */
    time_t t; /* processing time */
    struct tm *localtm; /* processing time */
    struct annotation *xannot;
    int grib_mode, npix = 4320, nlin = 2160;
    char grib2_t_str[300];

    /*
     * HDF datafile variables
     */

    int32 sdfid, fid, gridid, dimid, sdsid, sdsref;
    int32 datatype;

    /*
     * data type array pointers
     */

    float32 *float_SDSz_wind, *float_SDSm_wind, *float_SDSpress;
    float32 *float_SDSp_water, *float_SDSrel_hum;
    int16 *float_SDS_oz, *oz_newsiz;
    int8 *int8_SDSdataQC, *seaice8;

    /* 
     * external functions used 
     */

    int readgrib(), julian_(), wrtattr(), startHDF();
    int addAttr(), setSDSref();
    int32 wrtsds();
    void deattachHDFGrid(), pexit();
    int8 check_usage();
    struct annotation * fillenv_annot();

    /*
     * ------- check command line arguments and set args  ------------------
     */

    if ((check_usage(argc, argv, &anctyp, &n_opt_arg, source_name,
            &grib_mode, grib2_t_str)) != 0)
        exit(-1);
    strcpy(annotfile, argv[ n_opt_arg ]);
    strcpy(in_oz, argv[ 1 + n_opt_arg ]);
    if (argc - n_opt_arg < 3) strcpy(outdir, "./");
    else strcpy(outdir, argv[ 2 + n_opt_arg ]);

    if ((datarr1 = malloc(GRIBLATSZ * GRIBLONSZ * sizeof ( float)))
            == NULL) {
        printf("%s - %d: Unable to allocate space for datarr1 array\n", __FILE__,
                __LINE__);
        return 1;
    }
    /*
     * read the data
     */
    if (grib_mode == GRIB_MODE_2) {
        result = readgrib2(in_oz, npix, nlin, rgmode, datarr1);
        grib2_t(grib2_t_str, &year, &julval, &hour, &np, &nl, &fcst);
    } else {
        result = readgrib(in_oz, npix, nlin, datarr1, &year, &month, &day, &hour);
        julval = julian_(&day, &month, &year);
    }

    if (result == 1) {
        strcpy(message, "readgrib  file: ");
        strcat(message, in_oz);
        pexit(message);
    }

    /*
     * Create outfile name
     * Note that for grib 2 files, whole year is available already
     */
    if (grib_mode == GRIB_MODE_1) {
        if (year < 90) year = year + 2000; /* 2 digit 19xx or 20xx yrs to 4 */
        else year = year + 1900;
    }

    sprintf(outfile, "%s/N%04d%03d00_SEAICE_%s_24h.hdf",
            outdir, year, julval, source_name);

    sprintf(outfilename, "N%04d%03d00_SEAICE_%s_24h.hdf",
            year, julval, source_name);
    /* For SDPS' benefit! */
    printf("%s+%04d%03d000000000+%04d%03d235959999\n", outfile, year, julval,
            year, julval);


    /*
     * Write HDF annotations, lat and lon SDSs
     */

    /**** check for existing HDF  ***/

    /*  WDR/4/21/2015 removing the stupid failure if output is there already
     if (!access(outfile, F_OK)) pexit ("....output file exists");
     */

    if ((numannarr = count_annot(annotfile)) == 0)
        pexit("no annotations found");

    if ((xannot = (struct annotation *)
            malloc(sizeof (struct annotation) * numannarr))
            == NULL) pexit("malloc Annotation");

    xannot = fillenv_annot(annotfile);

    /*
     * Determine processing time
     */

    (void) time(&t);
    localtm = localtime(&t);
    pyear = localtm->tm_year;
    if (pyear < 90) pyear = pyear + 2000; /* 2 digit 19xx or 20xx yrs to 4 */
    else pyear = pyear + 1900;

    /*
     * -------  assign metadata values to local descriptive variables ----
     * -------  insert dates and other values to metadata array  ---------
     * -------  follow order of fillenv data file  -----------------------
     */

    for (i = 0; i < numannarr; i++) {

        /* Insert values to descr array element */

        if (!strcmp(xannot[i].label, "Product Name"))
            sprintf(xannot[i].descr, "%s", outfilename);

        else if (!strcmp(xannot[i].label, "Processing Time"))
            sprintf(xannot[i].descr, "%04d%03d%02d%02d%02d%03d",
                pyear, (localtm->tm_yday + 1), localtm->tm_hour, localtm->tm_min,
                localtm->tm_sec, 0);

        else if (!strcmp(xannot[i].label, "Input Files"))
            sprintf(xannot[i].descr, "%s", in_oz);

        else if (!strcmp(xannot[i].label, "Processing Control"))
            sprintf(xannot[i].descr, "%s %s %s %s",
                argv[0], annotfile, in_oz, outdir);

        else if (!strcmp(xannot[i].label, "Start Time"))
            sprintf(xannot[i].descr, "%04d%03d000000000", year, julval);

        else if (!strcmp(xannot[i].label, "End Time"))
            sprintf(xannot[i].descr, "%04d%03d230000000", year, julval);

        else if (!strcmp(xannot[i].label, "Start Year"))
            sprintf(xannot[i].descr, "%04d", year);

        else if (!strcmp(xannot[i].label, "Start Day"))
            sprintf(xannot[i].descr, "%03d", julval);

        else if (!strcmp(xannot[i].label, "Start Millisec"))
            sprintf(xannot[i].descr, "00000000");

        else if (!strcmp(xannot[i].label, "End Year"))
            sprintf(xannot[i].descr, "%04d", year);

        else if (!strcmp(xannot[i].label, "End Day"))
            sprintf(xannot[i].descr, "%03d", julval);

        else if (!strcmp(xannot[i].label, "End Millisec"))
            sprintf(xannot[i].descr, "82800000");

            /* Extract values from descr array element for program use */

        else if (!strcmp(xannot[i].label, "Northernmost Latitude"))
            sscanf(xannot[i].descr, "%f", &nmostlat);

        else if (!strcmp(xannot[i].label, "Southernmost Latitude"))
            sscanf(xannot[i].descr, "%f", &smostlat);

        else if (!strcmp(xannot[i].label, "Westernmost Longitude"))
            sscanf(xannot[i].descr, "%f", &wmostlon);

        else if (!strcmp(xannot[i].label, "Easternmost Longitude"))
            sscanf(xannot[i].descr, "%f", &emostlon);

        else if (!strcmp(xannot[i].label, "Latitude Step"))
            sscanf(xannot[i].descr, "%f", &latstep);
    }

    /*
     * Create HDF file
     */

    if ((result = startHDF(outfile, &sdfid, &fid, DFACC_CREATE)) != 0)
        pexit("Fatal error starting HDF file");

    /*
     * Write attribute array to HDF file
     */

    if ((result = wrtattr(sdfid, xannot, numannarr)) != 0) pexit("wrtattr");
    free(xannot);

    /*
     * -------- Allocate space for 2D data arrays  -----------------------
     */

    rank = 2;

    shape[0] = OUTLATSZ; /* lat */
    shape[1] = OUTLONSZ; /* lon */
    array_size = shape[0] * shape[1];

    if ((int8_SDSdataQC =
            (int8 *) calloc(array_size, sizeof (int8))) == NULL)
        pexit("calloc int8_SDSdataQC");

    if ((float_SDS_oz =
            (int16 *) malloc(sizeof (int16) * GRIBLATSZ * GRIBLONSZ)) == NULL)
        pexit("malloc float_SDS_oz");

    if ((oz_newsiz =
            (int16 *) malloc(sizeof (int16) * array_size)) == NULL)
        pexit("malloc oz_newsiz");
    if ((seaice8 = (int8 *) malloc(sizeof ( int8) * array_size)) == NULL)\
    pexit("malloc seaice8");

    /*
     *  Process ozone the same way here 
     */
    l = 0;
    for (j = 0; j < GRIBLATSZ; j++) {
        for (k = (GRIBLONSZ / 2); k < GRIBLONSZ; k++) {
            dat = *(datarr1 + j * GRIBLONSZ + k);
            if (dat < 0 || dat > SHRT_MAX)
                float_SDS_oz[l] = 0;
            else
                float_SDS_oz[l] = (int16) (100. * dat);
            l++; /* index counter */
        } /* for k */

        for (k = 0; k < (GRIBLONSZ / 2); k++) { /* lons */
            dat = *(datarr1 + j * GRIBLONSZ + k);
            if (dat < 0 || dat > SHRT_MAX)
                float_SDS_oz[l] = 0;
            else
                float_SDS_oz[l] = (int16) (100. * dat);
            l++; /* index counter */
        } /* for k */
    } /* for j */

    /*
     *  interpolate the array to the traditional TOVS ozone size
     */
    resize_oz((short *) float_SDS_oz, GRIBLONSZ, GRIBLATSZ,
            OUTLONSZ, OUTLATSZ, (short *) oz_newsiz);
    int_char_min = 0;
    int_char_max = 100;
    for (j = 0; j < array_size; j++) {
        if (*(oz_newsiz + j) < int_char_min) {
            *(int8_SDSdataQC + j) = 1;
            *(seaice8 + j) = -1;
        } else if (*(oz_newsiz + j) > int_char_max) {
            *(int8_SDSdataQC + j) = 1;
            *(seaice8 + j) = 101;
        } else
            *(seaice8 + j) = *(oz_newsiz + j);
    }
    /* 
     * ----------------- Create Vgroup  ----------------------------
     */

    gridid = setupGrid(fid, VGROUPNAME);

    /*
     * ----------------- Write ozone SDS ----------------------------
     */

    sdsname = "seaice";
    strcpy(longname, "Sea ice concentration");
    units = "percent coverage";
    datafmt = "int8";
    datatype = DFNT_INT8;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            seaice8, gridid)) != 0)
        pexit("SDSinFile ozone");


    // Smoothed 1/2 degree icemask for Aquarius processing
    float frac_ice_smooth[720 * 360];
    mk_smooth_ice_map_((char *) seaice8, frac_ice_smooth);

    sdsname = "seaice_smoothed";
    strcpy(longname, "Sea ice concentration (smoothed)");
    units = "percent coverage";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;
    shape_smooth[0] = 360; /* lat */
    shape_smooth[1] = 720; /* lon */
    SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape_smooth,
            frac_ice_smooth, gridid);

    free(float_SDS_oz);

    /*
     * ---- Write QC flag SDS, units attribute, and add to grid Vgroup -----
     */

    sdsname = "seaice_QC";
    strcpy(longname, "Sea ice Q/C flag");
    units = "";
    datafmt = "int8";
    datatype = DFNT_INT8;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            int8_SDSdataQC, gridid)) != 0)
        pexit("SDSinFile ozone_QC");

    /*
     * free storage for QC array (used for all products)
     * also free the input array
     */

    free(int8_SDSdataQC);
    free(datarr1);
    free(seaice8);
    free(oz_newsiz);
    /*
     * deattach HDF grid (used for all products)
     */

    deattachHDFgrid(gridid);

    /*
     * close HDF structures
     */

    if ((result = closeHDFstructs(sdfid, fid)) != 0) pexit("closeHDFstructs");

    return 0;
} /* main */

int8 check_usage(int argc, char *argv[], int *anctyp, int *n_opt_arg,
        char *source_name, int *grib_mode, char *grib2_t_str)
/*****************************************************************
   File:        check_usage
  
   Purpose:     check command line arguments
  
   Description: check command line arguements and report proper usage
                on argument count error or no argumnts.
  
   Returns:     success (0) or failure(-1)

      Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               argc             I      arg count
      char *[]          argv             I      input arg list
      int *             anctyp           O      type of anc data, either
                                                ANC_MET or ANC_OZONE
      int *             n_opt_arg        O      # option args
      char *            source_name      O      An alternate source (as in gbl 
                                              attrib source) name
      int *             grib_mode        O      Type of grib file the binary
                                              file came from: GRIB_MODE_1 for
                                              a GRIB 1 file, GRIB_MODE_2 for a
                                              GRIB 2 file (NCEP std as of
                                              Jan 2008)
      char *            grib2_t_str       O      string with time of the GRIB 2
                                                data - used for setting file
                                                name and attributes
  
   Note:        This routine is custom for each program that uses it since
                the arguments are different of course.
  
   Author:      Brian D. Schieber, GSC, 10/93
  
   Mods:           10/4/93 BDS, new HDF inputs handle all 3 NCEP parms
            W. Robinson, GSC, 4 Feb 97  changes for 5th param:rel_hum
            W. Robinson, GSC, 18 Jun 98  add flexability to either accept
                7 or 8 inputs for met data or 3 or 4 for ozone data
            W. Robinson, SAIC, 12 May 2005  have optional argument -s
                for an optional source name
            W. Robinson, SAIC, 24 Jan 2008  add a grib mode -2 for grib 2 
               format files
            J. Gales, Futuretech, 20 Jan 2010  Add version number
 *****************************************************************/ {
    extern char *optarg;
    extern int optind, optopt;
    int c, opt_source_found, errflg;

    /*
     *  get the options first, one for the source name
     */
    *grib_mode = GRIB_MODE_1;
    opt_source_found = 0;
    errflg = 0;
    *n_opt_arg = 0;
    while ((c = getopt(argc, argv, ":s:2:")) != -1) {
        switch (c) {
        case 's':
            strcpy(source_name, optarg);
            opt_source_found = 1;
            break;
        case '2':
            strcpy(grib2_t_str, optarg);
            *grib_mode = GRIB_MODE_2;
            break;
        case ':': /* -s without operand */
            fprintf(stderr,
                    "Option -%c requires an operand\n", optopt);
            errflg++;
            break;
        case '?':
            fprintf(stderr, "Unrecognized option: -%c\n", optopt);
            errflg++;
        }
    }
    if (!errflg) {
        *n_opt_arg = optind;
        /*
         *  for regular arguments
         */
        if (argc - optind == 3 || argc - optind == 4) {
            if (!opt_source_found) strcpy(source_name, "NCEP");
        } else {
            errflg++;
        }
    }
    if (errflg) {
        printf("\nVersion 1.0 (%s  %s)", __DATE__, __TIME__);
        printf("\n\nUsage:\n");
        printf("\t%s [-s<sname>] <metadata> <file> [outdir]\n", argv[0]);
        printf("\nWhere:\n");
        printf("\tmetadata:   DAAC-style HDF metadata (ASCII file)\n");
        printf("\tfile:       GRIB converted binary file to process\n");
        printf("\t(this can be either from a GRIB or a GRIB2 file\n");
        printf("\toutdir:     Optional output directory ('.' default)\n");
        printf("\t   Options:\n");
        printf("\t-s <sname>   Optional source name to output file with:\n");
        printf("\t   file name form: NYYYYDDD00_SEAICE_<sname>_24h.hdf\n");
        printf("\t   default is NCEP\n");
        printf(
                "\t-2 <GRIB str>   Specify the -2 with the <GRIB str> to process\n");
        printf(
                "\t        binary files derived from a GRIB 2 format file.  The\n");
        printf(
                "\t        binary files are derived using wgrib2 with the -bin \n");
        printf(
                "\t        option.  The <GRIB str> is a required date / time\n");
        printf(
                "\t        identification string gotten using wgrib2 -t\n");
        printf("\n");
        exit(-1);
    }
    return (0);
}
