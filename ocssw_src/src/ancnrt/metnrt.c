/***************************************************************** 
 * File:        metnrt.c
 *
 * Purpose:     create HDF ancillary real time datafile from 
 *              NCEP (formerly NMC) (GRIB format) W, P, PW, H input.
 *
 * Description: this program will create a HDF with the components:
 *
 *              HDF DAAC compliant Annotation object  
 *              1-D Latitude SDS
 *              1-D Longitude SDS
 *              2-D data parameter SDS (z_wind)
 *              2-D Q/C flag           (z_wind_QC)
 *              2-D data parameter SDS (m_wind)
 *              2-D Q/C flag           (m_wind_QC)
 *              2-D data parameter SDS (press)
 *              2-D Q/C flag           (press_QC)
 *              2-D data parameter SDS (rel_hum)
 *              2-D Q/C flag           (rel_hum_QC)
 *              2-D data parameter SDS (p_water)
 *              2-D Q/C flag           (p_water_QC)
 *
 * Input parms: 
 *    char *annotfile - input of DAAC compliant annotations ("fillenv.met")
 *    char *inz_wind  - input to process ("z_wind.dat")
 *    char *inm_wind  - input to process ("m_wind.dat")
 *    char *inpress   - input to process ("press.dat")
 *    char *inp_water - input to process ("p_water.dat")
 *    char *inrel_hum   - input to process ("rel_hum.dat")
 *    char *outdir    - output directory to write file
 *
 * Output parms:none            
 *
 *       output name is derived from GRIB converted file header. 
 *       Example name: "S199312300_NCEP.MET"
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
 *              readgrid_  - FORTRAN routine to read binary 
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
 *****************************************************************/

#include "ancil.h"
#include <time.h>

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

/* #define VSIZE         2.5 */
/* #define HSIZE         2.5 */
/* #define MAX_EAST    177.5 */
/* #define WPHLATSZ     73 */ /* latitude  of GRIB files */
/* #define WPHLONSZ    144 */ /* long of MODIFIED GRIB files */

main(int argc, char *argv[]) {
    int i, j, k, l;
    int rank, SDSref = 0;
    int result = 0;
    int array_size = 0;
    int numannarr = 0;
    int one = 1;
    int year = 0, month = 0, day = 0, hour = 0, msec = 0,
            julval = 0, pyear = 0;
    int32 shape[2];
    /* float32   datarr1[WPHLATSZ][WPHLONSZ + 1];     actual GRID has extra LON  */
    /* float32   datarr2[WPHLATSZ][WPHLONSZ + 1];     */
    /* float32   datarr3[WPHLATSZ][WPHLONSZ + 1];     */
    /* float32   datarr4[WPHLATSZ][WPHLONSZ + 1];     */

    float32 datarr1[WPHLATSZ][WPHLONSZ];
    float32 datarr2[WPHLATSZ][WPHLONSZ];
    float32 datarr3[WPHLATSZ][WPHLONSZ];
    float32 datarr4[WPHLATSZ][WPHLONSZ];
    float32 datarr5[WPHLATSZ][WPHLONSZ];

    float32 latscl[WPHLATSZ], lonscl[WPHLONSZ],
            dlat, dlon;
    char message[MAXNAMELNG];
    char inz_wind[MAXNAMELNG];
    char inm_wind[MAXNAMELNG];
    char inpress[MAXNAMELNG];
    char inp_water[MAXNAMELNG];
    char inrel_hum[MAXNAMELNG];
    char longname[MAXNAMELNG];
    char *sdsname;
    char *units;
    char *dataattr;
    char *longnamelabel;
    char *longnamevalue;
    char *unitslabel;
    char *unitsvalue;
    char *datafmt;
    char outfile[MAXNAMELNG];
    char outfilename[MAXNAMELNG];
    char outdir[MAXNAMELNG];
    char annotfile[MAXNAMELNG];
    float latstep, lonstep; /* lat/lon incr. */
    int latsz, lonsz; /* lat/lon size */
    float nmostlat, smostlat, wmostlon, emostlon; /* lat/lon corners */
    time_t t; /* processing time */
    struct tm *localtm; /* processing time */
    struct annotation *xannot;


    /*
     * HDF datafile variables
     */

    int32 sdfid, fid, gridid, dimid, sdsid, sdsref;
    int32 datatype;

    /*
     * data type array pointers
     */

    float32 *float_SDSz_wind;
    float32 *float_SDSm_wind;
    float32 *float_SDSpress;
    float32 *float_SDSp_water;
    float32 *float_SDSrel_hum;
    int8 *int8_SDSdataQC;

    /* 
     * external functions used 
     */

    int readgrib_();
    void julday_();
    int wrtattr();
    int32 wrtsds();
    int startHDF();
    int setupGrid();
    void deattachHDFGrid();
    int writeGeom();
    int addAttr();
    int setSDSref();
    void pexit();
    int8 check_usage();
    struct annotation * fillenv_annot();

    /*
     * ------- check command line arguments and set args  ------------------
     */

    if ((check_usage(argc, argv)) != 0) exit(-1);
    strcpy(annotfile, argv[1]);
    strcpy(inz_wind, argv[2]);
    strcpy(inm_wind, argv[3]);
    strcpy(inpress, argv[4]);
    strcpy(inp_water, argv[5]);
    strcpy(inrel_hum, argv[6]);
    if (argc < 8) strcpy(outdir, "./");
    else strcpy(outdir, argv[7]);

    /*
     * ------- Read each binary file, extract data array and time  ----------
     */

    result = readgrib_(inz_wind, datarr1, &year, &month, &day, &hour, MAXNAMELNG);
    if (result == 1) {
        strcpy(message, "readgrib_  file: ");
        strcat(message, inz_wind);
        pexit(message);
    }

    result = readgrib_(inm_wind, datarr2, &year, &month, &day, &hour, MAXNAMELNG);
    if (result == 1) {
        strcpy(message, "readgrib_  file: ");
        strcat(message, inm_wind);
        pexit(message);
    }

    result = readgrib_(inpress, datarr3, &year, &month, &day, &hour, MAXNAMELNG);
    if (result == 1) {
        strcpy(message, "readgrib_  file: ");
        strcat(message, inpress);
        pexit(message);
    }

    result = readgrib_(inp_water, datarr4, &year, &month, &day, &hour, MAXNAMELNG);
    if (result == 1) {
        strcpy(message, "readgrib_  file: ");
        strcat(message, inp_water);
        pexit(message);
    }

    result = readgrib_(inrel_hum, datarr5, &year, &month, &day, &hour, MAXNAMELNG);
    if (result == 1) {
        strcpy(message, "readgrib_  file: ");
        strcat(message, inrel_hum);
        pexit(message);
    }

    /* 
     *** This used to dump the GRIB array to stdout for a data check ***
      dlat = 2.5;
      dlon = 2.5;
      for (i = 0; i < WPHLATSZ; i++) {
         for (j = 0; j < WPHLONSZ; j++) {
            printf("\nLat: %f", 90.0 - (i * dlat));
            printf(" Lon: %f", (-180.0 + (j * dlon)));
            printf(" %f ", datarr1[i][j]/100.0);
         }
      }
     */

    /*
     * Create outfile name
     */

    if (year < 90) year = year + 2000; /* 2 digit 19xx or 20xx yrs to 4 */
    else year = year + 1900;

    julval = julian_(&day, &month, &year);
    sprintf(outfile, "%s/S%04d%03d%02d_NCEP.MET",
            outdir, year, julval, hour);

    sprintf(outfilename, "S%04d%03d%02d_NCEP.MET",
            year, julval, hour);

    /* For SDPS' benefit! */
    printf("%s\n", outfile);

    /*
     * Write HDF annotations, lat and lon SDSs
     */

    /**** check for existing HDF  ***/

    if (!access(outfile, F_OK)) pexit("....output file exists");

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
                inz_wind, inm_wind, inpress, inrel_hum, inp_water);

        else if (!strcmp(xannot[i].label, "Processing Control"))
            sprintf(xannot[i].descr, "%s %s %s %s %s %s %s %s",
                argv[0], annotfile, inz_wind, inm_wind, inpress,
                inrel_hum, inp_water, outdir);

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

        else if (!strcmp(xannot[i].label, "Longitude Step"))
            sscanf(xannot[i].descr, "%f", &lonstep);

        else if (!strcmp(xannot[i].label, "Number of Rows"))
            sscanf(xannot[i].descr, "%d", &latsz);

        else if (!strcmp(xannot[i].label, "Number of Columns"))
            sscanf(xannot[i].descr, "%d", &lonsz);

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
    shape[0] = WPHLATSZ; /* lat */
    shape[1] = WPHLONSZ; /* lon */
    array_size = shape[0] * shape[1];

    if ((float_SDSz_wind =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSz_wind");

    if ((float_SDSm_wind =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSm_wind");

    if ((float_SDSpress =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSpress");

    if ((float_SDSp_water =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSp_water");

    if ((float_SDSrel_hum =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSrel_hum");

    if ((int8_SDSdataQC =
            (int8 *) calloc(array_size, sizeof (int8))) == NULL)
        pexit("calloc int8_SDSdataQC");


    /* 
     * Process z/m_winds, press, rel_hum, p_water
     *  Shift all arrays so that 0 lon is centered, divide pressure by 100.
     */

    l = 0;
    for (j = 0; j < shape[0]; j++) { /* lats */
        for (k = (shape[1] / 2); k < shape[1]; k++) { /* lons */
            float_SDSz_wind[l] = datarr1[j][k];
            float_SDSm_wind[l] = datarr2[j][k];
            float_SDSpress[l] = datarr3[j][k] / 100.0;
            float_SDSp_water[l] = datarr4[j][k];
            float_SDSrel_hum[l] = datarr5[j][k];
            l++; /* index counter */
        } /* for k */

        for (k = 0; k < (shape[1] / 2); k++) { /* lons */
            float_SDSz_wind[l] = datarr1[j][k];
            float_SDSm_wind[l] = datarr2[j][k];
            float_SDSpress[l] = datarr3[j][k] / 100.0;
            float_SDSp_water[l] = datarr4[j][k];
            float_SDSrel_hum[l] = datarr5[j][k];
            l++; /* index counter */
        } /* for k */
    } /* for j */

    /* 
     * ----------------- Create Vgroup  ----------------------------
     */

    gridid = setupGrid(fid, VGROUPNAME);

    /* 
     * ----------------- Write z_wind SDS ----------------------------
     */


    sdsname = "z_wind";
    strcpy(longname, "Zonal wind at 10 m");
    units = "m sec^-1";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            float_SDSz_wind, gridid)) != 0)
        pexit("SDSinFile z_wind");

    free(float_SDSz_wind);

    /*
     * ---- Write QC flag SDS, units attribute, and add to grid Vgroup -----
     */

    sdsname = "z_wind_QC";
    strcpy(longname, "Zonal wind at 10 m Q/C flag");
    units = "";
    datafmt = "int8";
    datatype = DFNT_INT8;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            int8_SDSdataQC, gridid)) != 0)
        pexit("SDSinFile z_wind_QC");


    /* 
     * ----------------- Write m_wind SDS ----------------------------
     */

    sdsname = "m_wind";
    strcpy(longname, "Meridional wind at 10 m");
    units = "m sec^-1";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            float_SDSm_wind, gridid)) != 0)
        pexit("SDSinFile m_wind");

    free(float_SDSm_wind);

    /*
     * ---- Write QC flag SDS, units attribute, and add to grid Vgroup -----
     */

    sdsname = "m_wind_QC";
    strcpy(longname, "Meridional wind at 10 m Q/C flag");
    units = "";
    datafmt = "int8";
    datatype = DFNT_INT8;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            int8_SDSdataQC, gridid)) != 0)
        pexit("SDSinFile m_wind_QC");

    /* 
     * ----------------- Write pressure SDS ----------------------------
     */


    sdsname = "press";
    strcpy(longname, "Atmospheric pressure at mean sea level");
    units = "millibars";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            float_SDSpress, gridid)) != 0)
        pexit("SDSinFile press");

    free(float_SDSpress);

    /*
     * ---- Write QC flag SDS, units attribute, and add to grid Vgroup -----
     */

    sdsname = "press_QC";
    strcpy(longname, "Atmospheric pressure at mean sea level Q/C flag");
    units = "";
    datafmt = "int8";
    datatype = DFNT_INT8;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            int8_SDSdataQC, gridid)) != 0)
        pexit("SDSinFile press_QC");

    /*
     * ----------------- Write humidity SDS ----------------------------
     */

    sdsname = "rel_hum";
    strcpy(longname, "Relative humidity at 1000 mb");
    units = "percent";
    datafmt = "float32";
    datatype = DFNT_FLOAT32;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            float_SDSrel_hum, gridid)) != 0)
        pexit("SDSinFile rel_hum");

    free(float_SDSrel_hum);

    /*
     * ---- Write QC flag SDS, units attribute, and add to grid Vgroup -----
     */

    sdsname = "rel_hum_QC";
    strcpy(longname, "Relative humidity at 1000 mb Q/C flag");
    units = "";
    datafmt = "int8";
    datatype = DFNT_INT8;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            int8_SDSdataQC, gridid)) != 0)
        pexit("SDSinFile rel_hum_QC");

    /*
     * ----------------- Write precipitable water SDS ----------------------------
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
     * ---- Write QC flag SDS, units attribute, and add to grid Vgroup -----
     */

    sdsname = "p_water_QC";
    strcpy(longname, "Precipitable water Q/C flag");
    units = "";
    datafmt = "int8";
    datatype = DFNT_INT8;

    if ((SDSinFile(sdsname, longname, units, datafmt,
            datatype, sdfid, rank, shape,
            int8_SDSdataQC, gridid)) != 0)
        pexit("SDSinFile p_water_QC");

    /*
     * free storage for QC array (used for all products)
     */

    free(int8_SDSdataQC);

    /*
     * deattach HDF grid (used for all products)
     */

    deattachHDFgrid(gridid);

    /*
     * close HDF structures
     */

    if ((result = closeHDFstructs(sdfid, fid)) != 0) pexit("closeHDFstructs");

} /* main */

/*****************************************************************
 * File:        check_usage
 *
 * Purpose:     check command line arguments
 *
 * Description: check command line arguements and report proper usage
 *              on argument count error or no argumnts.
 *
 * Returns:     success (0) or failure(-1)
 *
 * Subs called: none
 *
 * History:     none
 *
 * Note:        This routine is custom for each program that uses it since
 *              the arguments are different of course.
 *
 * Author:      Brian D. Schieber, GSC, 10/93
 *
 * Mods:           10/4/93 BDS, new HDF inputs handle all 3 NCEP parms
 *          W. Robinson, GSC, 4 Feb 97  changes for 5th param:rel_hum
 *****************************************************************/

int8 check_usage(int argc, char *argv[]) {
    if (argc < 7) {
        printf("\n\nUsage:\n");
        printf("\t%s <metadata> <file> <file> <file> <file> <file> [outdir]\n",
                argv[0]);
        printf("\nWhere:\n");
        printf("\tmetadata:   DAAC-style HDF metadata (ASCII file)\n");
        printf("\tfiles:      GRIB converted binary file to process\n");
        printf(
                "\t (z_wind.dat, m_wind.dat, press.dat, p_water.dat, rel_hum.dat)\n");
        printf("\toutdir:     Optional output directory ('.' default)\n");
        printf("\n");
        exit(-1);
    }
    return (0);
}
