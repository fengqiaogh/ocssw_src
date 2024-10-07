/***************************************************************** 
 * File:        ancnrt_aquarius.c
 *
 * Purpose:     create HDF ancillary real time datafile from 
 *              NCEP (formerly NMC) (unpacked GRIB format) W, P, PW, H input
 *		or a TOVS ozone unpacked GRIB file from NMC.
 *
 * Description: this program will create a HDF with the components:
 *
 *          Met data:
 *              HDF DAAC compliant Annotation object  
 *              1-D Latitude SDS
 *              1-D Longitude SDS
 *              2-D data parameter SDS (geohgt)
 *              2-D data parameter SDS (tmp_sfc)
 *              2-D data parameter SDS (press)
 *              2-D data parameter SDS (c_water)
 *              2-D data parameter SDS (soil_w)
 *
 *
 * Input parms: 
 *    char *annotfile - input of DAAC compliant annotations ("fillenv.met")
 *    char *ingeohgt  - input to process ("geohgt.dat")
 *    char *intmp_sfc  - input to process ("tmp_sfc.dat")
 *    char *inpress   - input to process ("press.dat")
 *    char *inc_water - input to process ("c_water.dat")
 *    char *inc_soil  - input to process ("soil_w.dat")
 *
 *    char *outdir    - output directory to write file
 *
 * Output parms:none            
 *
 *       output name is derived from GRIB converted file header. 
 *       Example name: "Q199312300_NCEP.MET"
 *
 * Returns:     run status    
 *
 * Local vars:  numerous variables for storing HDF annotations and SDSs.
 *
 * Subs called:         
 *              wrtattr    - write HDF annotation to file
 *              fill_annot - read HDF annotations from file
 *              wrtsds     - write SDS from array to file
 *              pexit      - print fatal errors
 *              check_usage- check command line args
 *              readgrib  - FORTRAN routine to read binary 
 *                              (GRIB extracted) files
 *
 * History:     none            
 *
 * Note:        This program is meant to generate SeaWiFS style HDF 
 *              ancillary real-time datafiles.
 *
 * Author:      Brian D. Schieber, GSC, 4/93
 *
 * Modification history:        
 *       BDS, 4/19/93 - HDF missing flag removed after discussions
 *                      with Jim F. and M.Darzi.
 *                    - NCEP data in Pascals, divided by 100 => millibars.
 *       BDS, 4/30/93 - Modified GRIB read WPHLONSZ to remove "world wrapped"
 *                      extra Lon. Will skip last longitude in HDF write. 
 *       BDS, 10/4/93 - Modified for new HDF design and to produce and 
 *                      output file which holds all three met fields (WPH).
 *       BDS, 11/8/93 - Modified to create z/m wind SDS rather than windspeed.
 *       BDS, 11/28/93- Modified MAX_EAST to 177.5 (from 180.0).
 *       BDS, 12/6/93 - Assuming no missing values, 0.0's passed through.
 *       BDS, 5/19/95 - Support new SeaWiFS v 2.7 product specs.
 *       BDS, 8/21/96 - Renamed 'perror' to 'pexit' to avoid conflict with
 *                      HDF4.0.
 *       BDS, 9/18/96 - Changed annot structure.  Modify to support new 1x1
 *                      GDAS1 files.
 *       BDS, 9/19/96 - Precipitable water has units kg/m^2
 *       W. Robinson, GSC, 4 Feb 97 put rh back in
 *       W. Robinson, GSC, 18 Jun 98  rename to ancnrt to include added 
 *                    role of processing the TOVS ozone data
 *       W. Robinson, GSC, 18 Dec 98 for the ozone output, write the
 *                    start and end times to be compatible with o3nrt
 *       W. Robinson, SAIC, 23 Jan 2008  to add ability to read the binary 
 *             files derived from the GRIB 2 format files (just a few 
 *             differences, see readgrib2 in readgrib.c)
 *       W. Robinson, SAIC, modify for added argument to readgrib2
 *       W. Robinson, SAIC, 16 Feb 2010  modify readgrib2, readgrib2_3d 
 *              call args npix, nlin
 *       W. Robinson, SAIC, 6 Apr 2015  new readgrib2 call and added grib2_t
 *              call - also name RH wit h2m AGL instead of entire atmosphere
 *       
 *****************************************************************/

#include "ancil.h"
#include <time.h>
#include <limits.h>
#include <unistd.h>
#include "ancnrt_proto.h"

/*
 * Met data specific settings
 */

#define VGROUPNAME   "Geophysical Data"

#define BIN_METH      2
#define REGISTRATION  CENTER
#define VSIZE         1.0
#define HSIZE         1.0
#define MAX_NORTH    90.0
#define MAX_SOUTH   -90.0
#define MAX_WEST   -180.0
#define MAX_EAST    180.0
#define WPHLATSZ    181                /* latitude  of GRIB files */
#define WPHLONSZ    360                /* long of MODIFIED GRIB files */
#define WPHLATSZ_GFS 361               /* latitude  of GRIB GFS files */
#define WPHLONSZ_GFS 720               /* long of MODIFIED GRIB GFS files */
#define NPRFL        21
#define NPRFL_HGT    26

#define GRIB_MODE_1 1  /*  the GRIB modes, GRIB 1 old, GRIB 2 post jan 08 */
#define GRIB_MODE_2 2

#define VERSION "1.10"

/* #define VSIZE         2.5 */
/* #define HSIZE         2.5 */
/* #define MAX_EAST    177.5 */
/* #define WPHLATSZ     73 */ /* latitude  of GRIB files */
/* #define WPHLONSZ    144 */ /* long of MODIFIED GRIB files */

int main(int argc, char *argv[]) {
    int i, j, k, l, anctyp;
    int rank, rank_3d, SDSref = 0;
    int result = 0;
    int array_size = 0, array_size_gfs = 0;
    int numannarr = 0;
    int one = 1, rgmode = 0;
    int year = 0, month = 0, day = 0, hour = 0, msec = 0,
            julval = 0, pyear = 0;
    int32 shape[2], shape_3d[3], shape_gfs[2];
    int np, nl, fcst;
    int npix = WPHLONSZ, nlin = WPHLATSZ;
    int npix_gfs = WPHLONSZ_GFS, nlin_gfs = WPHLATSZ_GFS;

    float32 datarr1[WPHLATSZ][WPHLONSZ];
    float32 datarr2[WPHLATSZ][WPHLONSZ];
    float32 datarr3[WPHLATSZ][WPHLONSZ];
    float32 datarr4[WPHLATSZ][WPHLONSZ];
    float32 datarr5[WPHLATSZ][WPHLONSZ];
    float32 datarr6[WPHLATSZ][WPHLONSZ];

    float32 latscl[WPHLATSZ], lonscl[WPHLONSZ],
            dlat, dlon;
    char message[MAXNAMELNG], ingeohgt[MAXNAMELNG], intmp_sfc[MAXNAMELNG];
    char inpress[MAXNAMELNG], inp_water[MAXNAMELNG], inc_water[MAXNAMELNG];
    char insoil_w[MAXNAMELNG], intmp_bgr[MAXNAMELNG];
    ;
    char inu_wind[MAXNAMELNG], inv_wind[MAXNAMELNG];
    char inrh[MAXNAMELNG];
    char ingeohgt_prl[MAXNAMELNG], intmp_prl[MAXNAMELNG];
    char inrh_prl[MAXNAMELNG], inclwmr_prl[MAXNAMELNG];
    char insnow_gfs[MAXNAMELNG];
    char longname[MAXNAMELNG], *sdsname;
    char *units, *dataattr, *longnamelabel, *longnamevalue;
    char *unitslabel, *unitsvalue, *datafmt, outfile[MAXNAMELNG];
    char outfilename[MAXNAMELNG], outdir[MAXNAMELNG];
    char annotfile[MAXNAMELNG], source_name[100];
    int n_opt_arg;
    float latstep, lonstep; /* lat/lon incr. */
    int latsz, lonsz; /* lat/lon size */
    float nmostlat, smostlat, wmostlon, emostlon; /* lat/lon corners */
    time_t t; /* processing time */
    struct tm *localtm; /* processing time */
    struct annotation *xannot;
    int grib_mode;
    char grib2_t_str[300];

    /*
     * HDF datafile variables
     */

    int32 sdfid, fid, gridid, dimid, sdsid, sdsref;
    int32 datatype;

    /*
     * data type array pointers
     */

    float32 *float_SDSgeohgt, *float_SDStmp_sfc, *float_SDSpress;
    float32 *float_SDSrh, *float_SDSp_water, *float_SDSc_water;
    float32 *float_SDSsoil_w, *float_SDStmp_bgr;
    float32 *float_SDSgeohgt_prl, *float_SDStmp_prl, *float_SDSrh_prl;
    float32 *float_SDSclwmr_prl;
    float32 *float_SDSu_wind, *float_SDSv_wind;
    float32 *float_SDSsnow_gfs;
    //  int8      *int8_SDSdataQC;

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

    if ((check_usage(argc, argv, &n_opt_arg, source_name,
            &grib_mode, grib2_t_str)) != 0)
        exit(-1);

    if (argc == 1) {
        printf("%s %s (%s %s)\n\n", "ancnrt_aquarius", VERSION, __DATE__, __TIME__);
        printf("ancnrt_aquarius junk geohgt_grib2.dat tmp_sfc_grib2.dat\n");
        printf("press_grib2.dat rh_grib2.dat p_water_grib2.dat\n");
        printf("c_water_grib2.dat soil_w_grib2.dat\n");
        printf("u_wind_grib2.dat v_wind_grib2.dat\n");
        printf("geohgt_prl_grib2.dat tmp_prl_grib2.dat rh_prl_grib2.dat\n");
        printf("clwmr_prl_grib2.dat snow_gfs_grib2.dat tmp_bgr_grib2.dat -2 gtime\n");
        exit(0);
    }

    printf("%s %s (%s %s)\n\n", "ancnrt_aquarius", VERSION, __DATE__, __TIME__);

    strcpy(annotfile, argv[ n_opt_arg ]);

    strcpy(ingeohgt, argv[ 1 + n_opt_arg ]);
    strcpy(intmp_sfc, argv[ 2 + n_opt_arg ]);
    strcpy(inpress, argv[ 3 + n_opt_arg ]);
    strcpy(inrh, argv[ 4 + n_opt_arg ]);
    strcpy(inp_water, argv[ 5 + n_opt_arg ]);
    strcpy(inc_water, argv[ 6 + n_opt_arg ]);
    strcpy(insoil_w, argv[ 7 + n_opt_arg ]);
    strcpy(inu_wind, argv[ 8 + n_opt_arg ]);
    strcpy(inv_wind, argv[ 9 + n_opt_arg ]);

    strcpy(ingeohgt_prl, argv[10 + n_opt_arg ]);
    strcpy(intmp_prl, argv[11 + n_opt_arg ]);
    strcpy(inrh_prl, argv[12 + n_opt_arg ]);
    strcpy(inclwmr_prl, argv[13 + n_opt_arg ]);

    strcpy(insnow_gfs, argv[14 + n_opt_arg ]);
    strcpy(intmp_bgr, argv[15 + n_opt_arg ]);

    strcpy(outdir, argv[16 + n_opt_arg ]);

    /*
     * ------- Read each binary file, extract data array and time  ----------
     */

    /*
     * -------- Allocate space for 2D data arrays  -----------------------
     */

    rank = 2;

    shape[0] = WPHLATSZ; /* lat */
    shape[1] = WPHLONSZ; /* lon */
    array_size = shape[0] * shape[1];

    shape_gfs[0] = WPHLATSZ_GFS; /* lat */
    shape_gfs[1] = WPHLONSZ_GFS; /* lon */
    array_size_gfs = shape_gfs[0] * shape_gfs[1];

    if ((float_SDSgeohgt =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSgeohgt");

    if ((float_SDStmp_sfc =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDStmp_sfc");

    if ((float_SDSpress =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSpress");

    if ((float_SDSrh =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSrh");

    if ((float_SDSp_water =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSp_water");

    if ((float_SDSc_water =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSc_water");

    if ((float_SDSsoil_w =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSsoil_w");

    if ((float_SDSu_wind =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSu_wind");

    if ((float_SDSv_wind =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSv_wind");

    if (strcmp(insnow_gfs, "NONE") != 0) {
        if ((float_SDSsnow_gfs =
                (float32 *) malloc(sizeof (float32) * array_size_gfs)) == NULL)
            pexit("malloc float_SDSsnow_gfs");
    }

    if ((float_SDStmp_bgr =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDStmp_bgr");

    if (grib_mode == GRIB_MODE_2) {
        result = readgrib2(ingeohgt, npix, nlin, rgmode, float_SDSgeohgt);
        grib2_t(grib2_t_str, &year, &julval, &hour, &np, &nl, &fcst);
    } else
        result = readgrib(ingeohgt, npix, nlin, datarr1, &year, &month,
            &day, &hour);

    if (result == 1) {
        strcpy(message, "readgrib  file: ");
        strcat(message, ingeohgt);
        pexit(message);
    }

    result = readgrib2(intmp_sfc, npix, nlin, rgmode, float_SDStmp_sfc);

    if (result == 1) {
        strcpy(message, "readgrib  file: ");
        strcat(message, intmp_sfc);
        pexit(message);
    }

    result = readgrib2(inpress, npix, nlin, rgmode, float_SDSpress);

    if (result == 1) {
        strcpy(message, "readgrib  file: ");
        strcat(message, inpress);
        pexit(message);
    }

    result = readgrib2(inrh, npix, nlin, rgmode, float_SDSrh);

    if (result == 1) {
        strcpy(message, "readgrib  file: ");
        strcat(message, inrh);
        pexit(message);
    }

    result = readgrib2(inp_water, npix, nlin, rgmode, float_SDSp_water);

    if (result == 1) {
        strcpy(message, "readgrib  file: ");
        strcat(message, inp_water);
        pexit(message);
    }

    result = readgrib2(inc_water, npix, nlin, rgmode, float_SDSc_water);

    if (result == 1) {
        strcpy(message, "readgrib  file: ");
        strcat(message, inc_water);
        pexit(message);
    }

    result = readgrib2(insoil_w, npix, nlin, rgmode, float_SDSsoil_w);

    if (result == 1) {
        strcpy(message, "readgrib  file: ");
        strcat(message, insoil_w);
        pexit(message);
    }

    result = readgrib2(inu_wind, npix, nlin, rgmode, float_SDSu_wind);

    if (result == 1) {
        strcpy(message, "readgrib  file: ");
        strcat(message, inu_wind);
        pexit(message);
    }

    result = readgrib2(inv_wind, npix, nlin, rgmode, float_SDSv_wind);

    if (result == 1) {
        strcpy(message, "readgrib  file: ");
        strcat(message, inv_wind);
        pexit(message);
    }

    if (strcmp(insnow_gfs, "NONE") != 0) {
        result = readgrib2(insnow_gfs, npix_gfs, nlin_gfs, rgmode,
                float_SDSsnow_gfs);

        if (result == 1) {
            strcpy(message, "readgrib  file: ");
            strcat(message, insnow_gfs);
            pexit(message);
        }
    }

    result = readgrib2(intmp_bgr, npix, nlin, rgmode, float_SDStmp_bgr);

    if (result == 1) {
        strcpy(message, "readgrib  file: ");
        strcat(message, intmp_bgr);
        pexit(message);
    }

    // Profile data

    rank_3d = 3;

    shape_3d[0] = NPRFL_HGT; /* prfl */
    shape_3d[1] = WPHLATSZ; /* lat */
    shape_3d[2] = WPHLONSZ; /* lon */
    array_size = shape_3d[0] * shape_3d[1] * shape_3d[2];

    if ((float_SDSgeohgt_prl =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSgeohgt_prl");

    if ((float_SDStmp_prl =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDStmp_prl");

    if (grib_mode == GRIB_MODE_2)
        result = readgrib2_3d(ingeohgt_prl, grib2_t_str, npix, nlin,
            float_SDSgeohgt_prl, NPRFL_HGT,
            &year, &month, &day, &hour);

    if (result == 1) {
        strcpy(message, "readgrib  file: ");
        strcat(message, ingeohgt_prl);
        pexit(message);
    }


    if (grib_mode == GRIB_MODE_2)
        result = readgrib2_3d(intmp_prl, grib2_t_str, npix, nlin,
            float_SDStmp_prl, NPRFL_HGT,
            &year, &month, &day, &hour);

    if (result == 1) {
        strcpy(message, "readgrib  file: ");
        strcat(message, intmp_prl);
        pexit(message);
    }


    shape_3d[0] = NPRFL; /* prfl */
    shape_3d[1] = WPHLATSZ; /* lat */
    shape_3d[2] = WPHLONSZ; /* lon */
    array_size = shape_3d[0] * shape_3d[1] * shape_3d[2];

    if ((float_SDSrh_prl =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSrh_prl");

    if ((float_SDSclwmr_prl =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSclwmr_prl");

    if (grib_mode == GRIB_MODE_2)
        result = readgrib2_3d(inrh_prl, grib2_t_str, npix, nlin, float_SDSrh_prl,
            NPRFL,
            &year, &month,
            &day, &hour);

    if (result == 1) {
        strcpy(message, "readgrib  file: ");
        strcat(message, inrh_prl);
        pexit(message);
    }

    if (grib_mode == GRIB_MODE_2)
        result = readgrib2_3d(inclwmr_prl, grib2_t_str, npix, nlin,
            float_SDSclwmr_prl, NPRFL,
            &year, &month,
            &day, &hour);

    if (result == 1) {
        strcpy(message, "readgrib  file: ");
        strcat(message, inclwmr_prl);
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

    if (grib_mode == GRIB_MODE_1)
        julval = julian_(&day, &month, &year);
    //  sprintf(outfile, "%s/Q%04d%03d%02d_%s.MET",
    //	  outdir, year, julval, hour, source_name );

    //sprintf(outfilename, "Q%04d%03d%02d_%s.MET",
    //	  year, julval, hour, source_name );

    sprintf(outfile, "%s/N%04d%03d%02d_QMET_NCEP_6h",
            outdir, year, julval, hour);

    sprintf(outfilename, "N%04d%03d%02d_QMET_NCEP_6h",
            year, julval, hour);
    printf("%s\n", outfile);



    /*
     * Write HDF annotations, lat and lon SDSs
     */

    /**** check for existing HDF  ***/
    /*
    // JMG
      if (!access(outfile, F_OK)) pexit ("....output file exists");

      if ((numannarr = count_annot(annotfile)) == 0) 
          pexit ("no annotations found");

      if ((xannot = (struct annotation *) 
                    malloc (sizeof(struct annotation) * numannarr)) 
                    == NULL) pexit  ("malloc Annotation");

      xannot = fillenv_annot(annotfile);
     */

    /*
     * Determine processing time
     */

    (void) time(&t);
    localtm = localtime(&t);
    pyear = localtm->tm_year;
    if (pyear < 90) pyear = pyear + 2000; /* 2 digit 19xx or 20xx yrs to 4 */
    else pyear = pyear + 1900;

    if (hour > 0) msec = hour * 60 * 60 * 1000;
    else msec = 0;

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
            sprintf(xannot[i].descr, "%s %s %s %s %s",
                ingeohgt, intmp_sfc, inpress, inc_water, insoil_w);

        else if (!strcmp(xannot[i].label, "Processing Control"))
            sprintf(xannot[i].descr, "%s %s %s %s %s %s %s %s",
                argv[0], annotfile, ingeohgt, intmp_sfc, inpress,
                inc_water, insoil_w, outdir);

        else if (!strcmp(xannot[i].label, "Start Time"))
            sprintf(xannot[i].descr, "%04d%03d%02d0000000", year, julval, hour);

        else if (!strcmp(xannot[i].label, "End Time"))
            sprintf(xannot[i].descr, "%04d%03d%02d0000000", year, julval, hour);

        else if (!strcmp(xannot[i].label, "Start Year"))
            sprintf(xannot[i].descr, "%04d", year);

        else if (!strcmp(xannot[i].label, "Start Day"))
            sprintf(xannot[i].descr, "%03d", julval);

        else if (!strcmp(xannot[i].label, "Start Millisec"))
            sprintf(xannot[i].descr, "%08d", msec);

        else if (!strcmp(xannot[i].label, "End Year"))
            sprintf(xannot[i].descr, "%04d", year);

        else if (!strcmp(xannot[i].label, "End Day"))
            sprintf(xannot[i].descr, "%03d", julval);

        else if (!strcmp(xannot[i].label, "End Millisec"))
            sprintf(xannot[i].descr, "%08d", msec);

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

            /*
             *  WDR now a variable longitude step, # rows and columns
             */
        else if (!strcmp(xannot[i].label, "Longitude Step")) {
            lonstep = 1.;
            sprintf(xannot[i].descr, "%f", lonstep);
        }
        else if (!strcmp(xannot[i].label, "Number of Rows")) {
            latsz = 181;
            sprintf(xannot[i].descr, "%d", latsz);
        }
        else if (!strcmp(xannot[i].label, "Number of Columns")) {
            lonsz = 360;
            sprintf(xannot[i].descr, "%d", lonsz);
        }
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
    //  free(xannot);


    /* 
     * ----------------- Create Vgroup  ----------------------------
     */

    gridid = setupGrid(fid, VGROUPNAME);

    /* 
     * ----------------- Write geohgt SDS ----------------------------
     */
    sdsname = "geohgt";
    strcpy(longname, "Geopotential height, surface");
    units = "gpm";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            float_SDSgeohgt, gridid)) != 0)
        pexit("SDSinFile geohgt");

    free(float_SDSgeohgt);

    /* 
     * ----------------- Write tmp_sfc SDS ----------------------------
     */
    sdsname = "tmp_sfc";
    strcpy(longname, "Temperature at surface");
    units = "degrees K";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            float_SDStmp_sfc, gridid)) != 0)
        pexit("SDSinFile tmp_sfc");

    free(float_SDStmp_sfc);


    /* 
     * ----------------- Write pressure SDS ----------------------------
     */
    sdsname = "press";
    strcpy(longname, "Atmospheric pressure at surface");
    units = "Pascals";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            float_SDSpress, gridid)) != 0)
        pexit("SDSinFile press");

    free(float_SDSpress);

    /* 
     * ----------------- Write relative humidity 2 m AGL SDS ---------
     */
    sdsname = "rh";
    strcpy(longname, "Relative Humidity, 2 m above ground level");
    units = "percent";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            float_SDSrh, gridid)) != 0)
        pexit("SDSinFile rh");

    free(float_SDSrh);


    /*
     * ----------------- Write Precipitable water SDS -------------------------
     */
    sdsname = "p_water";
    strcpy(longname, "Precipitable water");
    units = "kg m^-2";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            float_SDSp_water, gridid)) != 0)
        pexit("SDSinFile p_water");

    free(float_SDSp_water);


    /*
     * ----------------- Write cloud water SDS ----------------------------
     */
    sdsname = "c_water";
    strcpy(longname, "Cloud water");
    units = "kg m^-2";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            float_SDSc_water, gridid)) != 0)
        pexit("SDSinFile c_water");

    free(float_SDSc_water);


    /*
     * ----------------- Write soil moisture SDS ----------------------------
     */
    sdsname = "soil_w";
    strcpy(longname, "Soil moisture");
    units = "";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            float_SDSsoil_w, gridid)) != 0)
        pexit("SDSinFile soil_w");

    free(float_SDSsoil_w);


    /*
     * ----------------- Write u wind SDS ----------------------------
     */
    sdsname = "u_wind";
    strcpy(longname, "u wind 10 m");
    units = "m sec^-1";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            float_SDSu_wind, gridid)) != 0)
        pexit("SDSinFile u_wind");

    free(float_SDSu_wind);


    /*
     * ----------------- Write v wind SDS ----------------------------
     */
    sdsname = "v_wind";
    strcpy(longname, "v wind 10 m");
    units = "m sec^-1";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            float_SDSv_wind, gridid)) != 0)
        pexit("SDSinFile v_wind");

    free(float_SDSv_wind);


    /*
     * ----------------- Write GFS snow SDS ----------------------------
     */
    sdsname = "snow_gfs";
    strcpy(longname, "Snow (GFS)");
    units = "kg/m^2";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;

    if (strcmp(insnow_gfs, "NONE") != 0) {
        if ((SDSinFile(sdsname, longname, units, datafmt,
                datatype, sdfid, rank, shape_gfs,
                float_SDSsnow_gfs, gridid)) != 0)
            pexit("SDSinFile snow_gfs");

        free(float_SDSsnow_gfs);
    }


    /* 
     * ----------------- Write geohgt profile SDS -------------------------
     */
    sdsname = "geohgt_prfl";
    strcpy(longname, "Geopotential height, profile");
    units = "m sec^-1";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;

    shape_3d[0] = NPRFL_HGT;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank_3d, shape_3d,
            float_SDSgeohgt_prl, gridid)) != 0)
        pexit("SDSinFile geohgt");

    free(float_SDSgeohgt_prl);

    /*
     * ----------------- Write tmp_bgr SDS ----------------------------
     */
    sdsname = "tmp_bgr";
    strcpy(longname, "Temperature 0-0.1 m below ground");
    units = "degrees K";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            float_SDStmp_bgr, gridid)) != 0)
        pexit("SDSinFile tmp_bgr");

    free(float_SDStmp_bgr);



    /* 
     * ----------------- Write tmp profile SDS -------------------------
     */
    sdsname = "tmp_prfl";
    strcpy(longname, "Temperature, profile");
    units = "degrees K";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank_3d, shape_3d,
            float_SDStmp_prl, gridid)) != 0)
        pexit("SDSinFile tmp");

    free(float_SDStmp_prl);


    /* 
     * ----------------- Write relative humidity profile SDS -----------------
     */
    sdsname = "rh_prfl";
    strcpy(longname, "Relative Humidity, profile");
    units = "percent";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;

    shape_3d[0] = NPRFL;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank_3d, shape_3d,
            float_SDSrh_prl, gridid)) != 0)
        pexit("SDSinFile rh");

    free(float_SDSrh_prl);


    /* 
     * ----------------- Write Cloud Mixing Ratio profile SDS -----------------
     */
    sdsname = "clwmr_prfl";
    strcpy(longname, "Cloud Water Mixing Ratio, profile");
    units = "";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank_3d, shape_3d,
            float_SDSclwmr_prl, gridid)) != 0)
        pexit("SDSinFile clwmr");

    free(float_SDSclwmr_prl);


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

int8 check_usage(int argc, char *argv[], int *n_opt_arg,
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
        strcpy(source_name, "NCEP");

        if (errflg) {
            printf("\n\nUsage:\n");
            printf("\t%s <metadata> <f1> <f2> <f3> <f4> <file> [outdir] [options]\n",
                    argv[0]);
            printf("\t%s [-s<sname>] <metadata> [outdir]\n", argv[0]);
            printf("\nWhere:\n");
            printf("\tmetadata:   DAAC-style HDF metadata (ASCII file)\n");
            printf("\tfiles:      GRIB converted binary file to process\n");
            printf("\tFor met data, 6 input files (fX):\n");
            printf(
                    "\t (geohgt.dat, tmp_sfc.dat, press.dat, c_water.dat soil_w.dat)\n");
            printf("\toutdir:     Optional output directory ('.' default)\n");
            printf("\t   Options:\n");
            printf("\t-s <sname>   Optional source name to output file with:\n");
            printf("\t   for MET data: QYYYYDDD00_<sname>.MET, default sname = NCEP\n");
            printf(
                    "\t-2 <GRIB str>   Specify the -2 with the <GRIB str> to process\n");
            printf(
                    "\t        binary files derived from a GRIB 2 format file.  The\n");
            printf(
                    "\t        binary files are derived using wgrib2 with the -bin \n");
            printf(
                    "\t        option.  The <GRIB str> is a required date / time\n");
            printf(
                    "\t        identification string gotten using wgrib -t\n");
            printf(
                    "\t         (see $SWFAPP/scripts/met_ingest2.scr for full example)\n");
            printf("\n");
            exit(-1);
        }
    }
    return (0);
}
