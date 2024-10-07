/*******************************************************************
 * File:        o3qc.c
 *
 * Purpose:     QC HDF ancillary real time OZONE datafile against
 *              climatology for the same month.
 *
 * Description: This program will create a HDF with:
 *              2-D of standard devs of real to clim (parameter = OZONE)
 *              A mailable (stdout) text file of %points in each STD category.
 *
 * The real time parameter is read into array RT. The climatological average
 * and standard deviation for the same time period (month) is read into
 * CLMAVG and CLMSTD respectively.  The CLMAVG array is subtracted from the
 * RT array and divided by the CLMSTD array to produce an array of the
 * number of std deviations from the mean for the real-time data.
 * As the data is processed a tally of points falling + and - 1, 2, and 3
 * STDs from the mean will be stored for later writing of percentages in
 * each category.
 *
 *  For OZONE data:
 *    The real time parameters are read into array ozoneRT
 *
 *    The climatological averages and standard deviations for the same time
 *    period (month) are read into:
 *       ozoneCLMAVG and ozoneCLMSTD respectively.
 *
 *    Real-time and climatology files are dimensioned 180x288.
 *    S19891991_TOMS.OZONE does not require regidding like the MET file did.
 *
 *
 * Input parms (command line arguments):
 *  o3qc <nrtfile> <climfile> <outfile> <monthstr> ozone <stdlimit1>
 *    <stdlimit2> <stdlimit3> <maxmissing> <maxmissing> [loval] [hival]
 *    [outmax] [min_lat] [max_lat] [gridflag] [ifileflag]
 *
 *  <nrtfile> (char *) - real-time input file to process
 *            ("a199310000_TOVS.OZONE")
 *  <climfile> (char *) - clim input file to process
 *            ("S19891991_TOMS.OZONE")
 *  <outfile> (char *) - output QC filename
 *            ("S199310000_TOVS.OZONE.qc")
 *  <monthstr> (char *) - month being processed
 *  <stdlimit1> (float32) - % w/in 1 STD
 *  <stdlimit2> (float32) - % w/in 2 STD
 *  <stdlimit3> (float32) - % w/in 3 STD
 *  <maxmissing> (int) - max num missing pts.
 *
 *  Optional inputs: (defaults are available, described below)
 *  [loval] (int) - lowest acceptable value in NRT file
 *  [hival] (int) - highest acceptable value in NRT file
 *  [outmax] (int) - # of points that can be outside lo/hival before
 *                an error flag is set
 *  [min_lat] (float) - minimum latitude above which to check, default -90.
 *  [max_lat] (float) - maximum latitude below which to check, default 90.
 *        NOTE that if min_lat = 95, the limits will be determined using
 *        the julian day in the equation:
 *        zen = -23. * cos( ( jday + 10 ) * 360 / 365 )
 *        min_lat = -90. + zen + mmoffset
 *        max_lat =  90. - zen - mmoffset
 *        where mmoffset = max_lat and defaults to 10. if nothing there
 *  [gridflag] (char) - 's' or 'd' for std dev or direct numeric difference
 *                     only the 's' option is valid for statistical reporting.
 *  [ifileflag] (int) - a 1 if the input realtime file name is actually a 
 *              climatology file.  0 or blank if not.
 *
 * Output parms:
 *   Example name: "S199312300_TOVS.OZONE.qc"
 *
 * Returns:     SUCCESS or FAILURE depending on threshold checks.
 *
 * Local vars:  numerous variables for storing HDFs.
 *
 * Subs called: int rdsds       - reads an HDF SDS data grid
 *              int8 check_usage  - confirm args
 *              int startHDF    - open HDF file for output
 *              int setupGrid    - setup HDF geometry struct
 *              int writeGeom    - write HDF geometry struct
 *              int wrtsds       - write HDF SDS array
 *              int addAttr      - add HDF SDS attributes
 *              int setSDSref    - set HDF SDS reference
 *              void deattachHDFgrid - detach HDF grid struct
 *              int closeHDFstructs - close HDF files
 *
 * History:     none
 *
 * Note:        This program is meant to generate HDF QC SDS data sets.
 *
 * Author:      Brian D. Schieber, GSC, 4/93
 *
 * Modification history:
 *   10/8/93 BDS, work with new HDF design format
 *   03/94 BDS, mods to run from script and work in CVDB
 *   05/94 BDS, outfile is passed in rather than derived.
 *   08/95 BDS, mod to reflect NRT and CLM HDF format changes (spec 2.7).
 *   11/95 BDS, add ability to sum number of values in image < or >
 *              two specified range arguments and return FAIL if > outmax.
 *           Values for LO(lower range), HI (upper range),
 *           and OUTMAX (MAX. Pts allowed outside range) hardcoded within.
 *           Also, print min and max value of data.
 *    8/96 BDS, renamed 'perror' to 'pexit' to avoid HDF4.0 conflict.
 *    W. Robinson, GSC, 13 Mar 97  add a latitude range to perform check 
 *         within (for use with TOMS data which misses the pole(s))
 *    6/97 KJS. Problems with SUN port due to strings not being
 *         long enough - gridflg, sdsname. Changed malloc to calloc-
 *         QC array outside lats was not being initialized. Hdiff issue.
 *    W. Robinson, GSC, 16 jun 97  add feature for day-of-year controlled
 *         min, max latitude
 *    W. Robinson, SAIC, 28 Sep 07  for count of high and low O3, weight by
 *         cos( lat )
 **************************************************************************/

#include "ancil.h"
#include "l1io.h"
#include <genutils.h>

/*
 * Ozone specific settings
 */

#define VGROUPCLASS  "PlanetaryGrid"
#define BIN_METH     2
#define REGISTRATION CENTER
#define VSIZE         1.25
#define HSIZE         1.00
#define MAX_NORTH    89.5
#define MAX_SOUTH   -89.5
#define MAX_WEST   -179.375
#define MAX_EAST    179.375
#define TOMSLATSZ   180
#define TOMSLONSZ   288
#define CLMfile 1

int main(int argc, char *argv[]) {
    int i, j, l;
    int rank;
    int result = 0;
    int missing = 0;
    int minval = 1000; /* initial minimum value of all points */
    int maxval = -1; /* initial maximum value of all points */
    int rtmiss = 0;
    int neg1cnt = 0;
    int neg2cnt = 0;
    int neg3cnt = 0;
    int pos1cnt = 0;
    int pos2cnt = 0;
    int pos3cnt = 0;
    int failed = 0;
    int totalpts = 0;
    int array_size = 0;
    int32 inShape[2];
    int32 shape[2];
    char infileRT[MAXNAMELNG];
    char infileCLM[MAXNAMELNG];
    char monthstr[10];
    char sdsname[MAXNAMELNG];
    char vgroupname[MAXNAMELNG];
    char gridflag[2];
    char datalabel[MAXNAMELNG];
    char *dataattr;
    char *dataunit;
    char outfile[MAXNAMELNG];
    int ifileflag = 0;
    float32 stdlimit1, stdlimit2, stdlimit3;
    int maxmissing;
    int loval = 100; /* default lo range threshold */
    int hival = 500; /* default hi range threshold */
    int outmax = 25; /* max num pts allowed to fall outside range */
    int locnt; /* pts under lo limit */
    int hicnt; /* pts above hi limit */
    float flocnt = 0., fhicnt = 0., coslat;

    int min_lin = (-20), max_lin = (TOMSLATSZ + 20), pts_used;
    float min_lat = -90., max_lat = 90., mmoffset = 10., zen, lat;
    float pi;
    int jday;

    int lo_cnts[5] = {0, 0, 0, 0, 0}, hi_cnts[5] = {0, 0, 0, 0, 0};
    int z_st[4] = {-1, -1, -1, -1}, z_en[4] = {-1, -1, -1, -1};
    int lo_vals[5] = {50, 100, 150, 200, 250},
    hi_vals[5] = {400, 450, 500, 550, 600};
    int in_z_run = 0, zcur = 0, zline, ipx, iv; /* for extra info on data */

    /*
     * HDF datafile variables
     */

    int32 sdfid, fid;
    int32 gridid, sdsid, geomid;
    int32 datatype;
    char vgname[MAXNAMELNG];

    /*
     * data type array pointers
     */

    int16 *int_SDSdataRT;
    int16 *int_SDSdataCLMAVG;
    int16 *int_SDSdataCLMSTD;
    float32 *float_SDSdataQC;

    /* external functions used */

    int8 check_usage();

    int anc_daymon(char *, int*, char*);
    /*
     * ------- check command line arguments and set args  ------------------
     */
    pi = acos(-1.);

    strcpy(gridflag, "s");
    if ((check_usage(argc, argv)) != 0) pexit("insufficient args provided");
    strcpy(infileRT, argv[1]);
    strcpy(infileCLM, argv[2]);
    strcpy(outfile, argv[3]);
    strcpy(monthstr, argv[4]);

    if (!strcmp(lowcase(argv[5]), "ozone"))
        strcpy(vgroupname, "Geophysical Data");
    else if (!strcmp(lowcase(argv[5]), "uwind"))
        pexit("Only OZONE parameter supported by this program.");
    else if (!strcmp(lowcase(argv[5]), "vwind"))
        pexit("Only OZONE parameter supported by this program.");
    else if (!strcmp(lowcase(argv[5]), "pres"))
        pexit("Only OZONE parameter supported by this program.");
    else if (!strcmp(lowcase(argv[5]), "rhum"))
        pexit("Only OZONE parameter supported by this program.");
    else pexit("Invalid Parameter string specified");

    sscanf(argv[6], "%f", &stdlimit1);
    sscanf(argv[7], "%f", &stdlimit2);
    sscanf(argv[8], "%f", &stdlimit3);
    sscanf(argv[9], "%d", &maxmissing);
    if (argc >= 11) sscanf(argv[10], "%d", &loval);
    if (argc >= 12) sscanf(argv[11], "%d", &hival);
    if (argc >= 13) sscanf(argv[12], "%d", &outmax);
    if (argc >= 14) sscanf(argv[13], "%f", &min_lat);
    if (argc >= 15) sscanf(argv[14], "%f", &max_lat);

    if (argc >= 16) {
        strncpy(gridflag, argv[15], 1);
        gridflag[1] = '\0';
    }
    if (argc == 17) ifileflag = atoi(argv[16]);

    /*
     * -------- Allocate space for 2D data arrays  -----------------------
     */

    rank = 2;
    shape[0] = TOMSLATSZ; /* lat */
    shape[1] = TOMSLONSZ; /* lon */
    array_size = shape[0] * shape[1];

    if ((int_SDSdataRT =
            /*      (int16 *) malloc (sizeof(int16) * array_size)) == NULL) */
            (int16 *) calloc(array_size, sizeof (int16))) == NULL)
        pexit("calloc int_SDSdataRT");

    if ((int_SDSdataCLMAVG =
            /*      (int16 *) malloc (sizeof(int16) * array_size)) == NULL) */
            (int16 *) calloc(array_size, sizeof (int16))) == NULL)
        pexit("calloc int_SDSdataCLMAVG");

    if ((int_SDSdataCLMSTD =
            /*      (int16 *) malloc (sizeof(int16) * array_size)) == NULL) */
            (int16 *) calloc(array_size, sizeof (int16))) == NULL)
        pexit("calloc int_SDSdataCLMSTD");

    if ((float_SDSdataQC =
            /*      (float32 *) malloc (sizeof(float32) * array_size)) == NULL) */
            (float32 *) calloc(array_size, sizeof (float32))) == NULL)
        pexit("calloc float_SDSdataQC");

    /*
     * ------- Determine appropriate real-time (SDS) to match clim. file -----
     *         (The "real-time" file could also be a clim. itself.)
     *         In OZONE NRT there are two SDS's in the realtime file:
     *          "ozone" and "ozone_QC"
     *
     *         In OZONE CLM there are 3 SDSs X 12 months or 36 SDS arrays.
     *         The below matching must consider or determine:
     *         - The NRT month desired
     *         - The corresponding CLM "Statistics" SDS
     */

    /*
     *  Read CLM type file if ifileflag == CLMfile.
     *  Read NRT type file otherwise and by default.
     */

    if (ifileflag == CLMfile) {
        strcpy(vgname, monthstr);
        strcpy(sdsname, "ozone_mean");
    } else {
        strcpy(vgname, vgroupname);
        strcpy(sdsname, "ozone");
        /*
         *  The month will be derived from the file start day rather than the 
         *  input month.
         *  Read start day attribute
         */
        if (anc_daymon(infileRT, &jday, monthstr) != 0)
            pexit("in getting day, month");
        /*
         *  add logic to get proper min, max lat
         */
        if (min_lat == 95.) {
            if (argc >= 15) mmoffset = max_lat;
            zen = -23. * cos((jday + 10.) * 2 * pi / 365.);
            min_lat = -90. + zen + mmoffset;
            if (min_lat < -90.) min_lat = -90.;
            max_lat = 90. + zen - mmoffset;
            if (max_lat > 90.) max_lat = 90.;
            printf("day dependent min, max latitude chosen:\n");
            printf("day = %d, mmoffset = %f, sun lat = %f\n",
                    jday, mmoffset, zen);
        }
    }

    if ((rdsds(infileRT, vgname, sdsname, inShape,
            int_SDSdataRT)) != 0) pexit("rdsds RT");

    if ((inShape[0] != shape[0]) || (inShape[1] != shape[1]))
        pexit("real-time dimensions not matching expected");

    /*
     * ------- Determine appropriate climatology month (SDS) to match ------
     *         real time file.  Get the average and std.dev. SDSs.
     */

    strcpy(vgname, monthstr);
    strcpy(sdsname, "ozone_mean");

    if ((rdsds(infileCLM, vgname, sdsname, inShape,
            int_SDSdataCLMAVG)) != 0) pexit("rdsds CLMAVG");

    if ((inShape[0] != shape[0]) || (inShape[1] != shape[1]))
        pexit("clim avg dimensions not matching expected");

    strcpy(sdsname, "ozone_std_dev");

    if ((rdsds(infileCLM, vgname, sdsname, inShape,
            int_SDSdataCLMSTD)) != 0) pexit("rdsds CLMSTD");

    if ((inShape[0] != shape[0]) || (inShape[1] != shape[1]))
        pexit("clim std dimensions not matching expected");

    /*  Out for now
       l = 0;
       for (i = 0; i < TOMSLATSZ; i++) {
          printf ("\n");
          for (j = 0; j < TOMSLONSZ; j++) {
             printf ("%d ", int_SDSdataRT[l]);
             l++;
          }
       }
     */

    /*
     * ------- Compute variance of real-time from climatology -----------
     *         and tally missing points, #+ & - 1,2,3 std dev
     *         for writing to file for mailing.
     */

    l = 0;
    missing = 0;
    locnt = 0;
    hicnt = 0;
    rtmiss = 0;
    neg1cnt = 0;
    neg2cnt = 0;
    neg3cnt = 0;
    pos1cnt = 0;
    pos2cnt = 0;
    pos3cnt = 0;
    pts_used = 0;

    for (i = 0; i < TOMSLATSZ; i++) {
        /*
         *  this code will find the 1st 4 runs of lines with all 0 values
         */
        zline = 1; /* say this line is a zero line */
        if (zcur < 4) /* only can do first 4 0 runs (should be max of 2 */ {
            for (ipx = 0; ipx < TOMSLONSZ; ipx++) /* look for != 0 pix i line */ {
                if (int_SDSdataRT[l + ipx] > 0) {
                    zline = 0; /* this is not a zero line */
                    if (in_z_run == 1) {
                        z_en[zcur] = i - 1; /*  if we find a non-zero line while in a */
                        zcur++; /*  run of zero lines, note the end of the */
                        in_z_run = 0; /*  run and say we are out of a zero run */
                    }
                    break;
                }
            }
            /*
             *  if a line is found = 0, we either start a run if we are not in one
             *  or continue in current zero run
             */
            if (zline == 1) {
                if (in_z_run == 0) {
                    z_st[ zcur ] = i;
                    in_z_run = 1;
                }
                if (i == TOMSLATSZ - 1) z_en[ zcur ] = i;
            }
        }
        /*** end run location code  ***/
        /*
         *  only do data within the min, max latitude range
         */
        lat = 90. - (i * 180. / TOMSLATSZ);
        coslat = cos(lat);
        if (lat >= min_lat && lat <= max_lat) {
            pts_used += TOMSLONSZ;
            if (i > min_lin) min_lin = i;
            if (i < max_lin) max_lin = i;

            for (j = 0; j < TOMSLONSZ; j++) {
                float_SDSdataQC[l] = 0.0; /* initialize each point */

                /* pt. is missing, NRT missing, a MIN val, and counted as out of
                 lower range */
                if ((int_SDSdataRT[l] <= 0) || (int_SDSdataCLMAVG[l] <= 0)) {
                    missing++;
                    if (int_SDSdataRT[l] <= 0) rtmiss++;
                    if (int_SDSdataRT[l] < minval) minval = int_SDSdataRT[l];
                    if (int_SDSdataRT[l] < loval) locnt++;
                    /*
                     *  get lo_cnts
                     */
                    for (iv = 0; iv < 5; iv++) {
                        if (int_SDSdataRT[l] < lo_vals[iv]) lo_cnts[iv]++;
                    }
                } else { /* non-missing point */

                    if ((!strcmp(gridflag, "d")) || (!strcmp(gridflag, "D"))) {
                        float_SDSdataQC[l] =
                                (float32) int_SDSdataRT[l] - int_SDSdataCLMAVG[l];
                    } else {
                        if (int_SDSdataCLMSTD[l] > 0)
                            float_SDSdataQC[l] = (float32)
                            (int_SDSdataRT[l] - int_SDSdataCLMAVG[l]) /
                            int_SDSdataCLMSTD[l];
                    }

                    /* find min and max of all points */
                    if (int_SDSdataRT[l] < minval) minval = int_SDSdataRT[l];
                    if (int_SDSdataRT[l] > maxval) maxval = int_SDSdataRT[l];

                    /* sum up points outside lo/hi range (LO adds to missing counts)*/
                    if (int_SDSdataRT[l] < loval) flocnt += coslat;
                    if (int_SDSdataRT[l] > hival) fhicnt += coslat;
                    /*
                     *  get lo_cnts and hi_cnts          */
                    for (iv = 0; iv < 5; iv++) {
                        if (int_SDSdataRT[l] < lo_vals[iv]) lo_cnts[iv]++;
                        if (int_SDSdataRT[l] > hi_vals[iv]) hi_cnts[iv]++;
                    }

                    if (float_SDSdataQC[l] < -1.0) neg1cnt++;
                    if (float_SDSdataQC[l] < -2.0) neg2cnt++;
                    if (float_SDSdataQC[l] < -3.0) neg3cnt++;

                    if (float_SDSdataQC[l] > 1.0) pos1cnt++;
                    if (float_SDSdataQC[l] > 2.0) pos2cnt++;
                    if (float_SDSdataQC[l] > 3.0) pos3cnt++;
                }

                /* Not in use
                      if ((j > 100) && (j < 150)) {
                         printf("QC  [%d][%d]: %f\n", i,j,float_SDSdataQC[l]);
                         printf("RT  [%d][%d]: %d\n", i,j,int_SDSdataRT[l]);
                         printf("AVG [%d][%d]: %d\n", i,j,int_SDSdataCLMAVG[l]);
                         printf("STD [%d][%d]: %d\n", i,j,int_SDSdataCLMSTD[l]);
                         }
                         printf("[%d][%d]: %d %d %d\n", i,j,int_SDSdataRT[l],
                          int_SDSdataCLMAVG[l], int_SDSdataCLMSTD[l]);
                 */
                l++;
            } /* for j */
        }
        else {
            l += TOMSLONSZ; /* we still have to move along thru the array */
        }
    } /* for i */

    totalpts = pts_used - missing;

    free(int_SDSdataRT);
    free(int_SDSdataCLMAVG);
    free(int_SDSdataCLMSTD);

    /*
     * ----------------- Write QC SDS ----------------------------
     */

    /*
     * Create HDF file
     */

    if ((startHDF(outfile, &sdfid, &fid, DFACC_CREATE)) != 0)
        pexit("Fatal error starting HDF file");

    /*
     * Create grid structure
     */

    /* gridid = setupGrid(fid, VGROUPCLASS, "QC difference data"); */
    gridid = setupGrid(fid, "QC difference data");

    /*
     * Write geometry Vdata, get ref, and add to grid Vgroup (geomid not used)
     */

    if ((geomid = writeGeom(fid, gridid, GEOMNAME, BIN_METH, REGISTRATION,
            VSIZE, HSIZE, MAX_NORTH, MAX_SOUTH,
            MAX_WEST, MAX_EAST)) == ERROR)
        pexit("Fatal error writing geometry");

    /*
     * Write SDS grid
     */

    strcpy(datalabel, "QC array");
    datatype = DFNT_FLOAT32;
    dataattr = "UNITS";
    dataunit = "std dev diff";

    /* datalabel, dataunit, datafmt, */
    if ((sdsid = wrtsds(sdfid, rank, shape, datatype,
            datalabel,
            float_SDSdataQC)) < 0) pexit("main wrtsds");

    free(float_SDSdataQC);

    /*
     * set SDS attribute
     */

    if ((result = addAttr(sdsid, dataattr, DFNT_CHAR, dataunit)) != 0)
        pexit("addAttr");

    /*
     * add SDS to Vgroup and deattach SDS
     */

    if ((result = setSDSref(sdsid, gridid)) != 0)
        pexit("setSDSref");

    /*
     * deattach HDF grid
     */

    deattachHDFgrid(gridid);

    /*
     * close HDF structures
     */

    if ((result = closeHDFstructs(sdfid, fid)) != 0) pexit("closeHDFstructs");

    /*
     *  make integer hi and low counts
     */
    locnt = (int) flocnt;
    hicnt = (int) fhicnt;
    /*
     * ----------------- Print stats to standard output -----------
     */


    if (!strcmp(gridflag, "s")) {
        printf("\n");
        printf(
                "---------------------------------------------------------------\n");
        printf("Results of comparison of real-time and climatological files:\n");
        printf("%s\t%s\n", infileRT, infileCLM);
        printf("Month: %s\tParameter: %s\n", monthstr, "Ozone");
        printf("Thresholds: %f %f %f Max Missings: %d\n",
                stdlimit1, stdlimit2, stdlimit3, maxmissing);
        printf("Lines considered range from %d (lat %f) to %d (lat %f)\n",
                max_lin, max_lat, min_lin, min_lat);

        printf("\n");
        printf("Minimum value: %d Maximum: %d\n", minval, maxval);
        printf("Total # non-missing values: %6d / considered: %d\n",
                totalpts, pts_used);
        printf("Total # values <%d: %d  >%d: %d  allowed: %d",
                loval, locnt, hival, hicnt, outmax);
        if (locnt >= outmax || hicnt >= outmax) {
            failed = 1;
            printf("  ***\n\n");
        } else {
            printf("\n\n");
        }

        /* if (locnt >= outmax) failed = 1;  lo points out of range exceeded */
        /* if (hicnt >= outmax) failed = 1;  hi points out of range exceeded */

        printf("Total # points/percentage  < -1 STD: %6d (%5.2f percent)\n",
                neg1cnt, (float) neg1cnt / totalpts * 100.0);
        printf("Total # points/percentage +/- 1 STD: %6d (%5.2f percent)",
                totalpts - (neg1cnt + pos1cnt),
                (float) (totalpts - (pos1cnt + neg1cnt)) / totalpts * 100.0);

        if (((float) (totalpts - (pos1cnt + neg1cnt)) / totalpts * 100.0) <
                stdlimit1) {
            failed = 1;
            printf("  ***\n");
        } else {
            printf("\n");
        }

        printf("Total # points/percentage  >  1 STD: %6d (%5.2f percent)\n",
                pos1cnt, (float) pos1cnt / totalpts * 100.0);
        printf("\n");

        printf("Total # points/percentage  < -2 STD: %6d (%5.2f percent)\n",
                neg2cnt, (float) neg2cnt / totalpts * 100.0);
        printf("Total # points/percentage +/- 2 STD: %6d (%5.2f percent)",
                totalpts - (neg2cnt + pos2cnt),
                (float) (totalpts - (pos2cnt + neg2cnt)) / totalpts * 100.0);

        if (((float) (totalpts - (pos2cnt + neg2cnt)) / totalpts * 100.0) <
                stdlimit2) {
            failed = 1;
            printf("  ***\n");
        } else {
            printf("\n");
        }

        printf("Total # points/percentage  >  2 STD: %6d (%5.2f percent)\n",
                pos2cnt, (float) pos2cnt / totalpts * 100.0);
        printf("\n");

        printf("Total # points/percentage  < -3 STD: %6d (%5.2f percent)\n",
                neg3cnt, (float) neg3cnt / totalpts * 100.0);
        printf("Total # points/percentage +/- 3 STD: %6d (%5.2f percent)",
                totalpts - (neg3cnt + pos3cnt),
                (float) (totalpts - (pos3cnt + neg3cnt)) / totalpts * 100.0);

        if (((float) (totalpts - (pos3cnt + neg3cnt)) / totalpts * 100.0) <
                stdlimit3) {
            failed = 1;
            printf("  ***\n");
        } else {
            printf("\n");
        }

        printf("Total # points/percentage  >  3 STD: %6d (%5.2f percent)\n",
                pos3cnt, (float) pos3cnt / totalpts * 100.0);
        printf("\n");

        printf("Total # missing values:       %6d\n", missing);
        printf("Total # missing real-time:    %6d, max allowed: %6d",
                rtmiss, maxmissing);
        if (rtmiss > maxmissing) {
            failed = 1;
            printf("  ***\n");
        } else {
            printf("\n");
        }
        printf(
                "---------------------------------------------------------------\n");
        /*
         *  do 2 line terse report of the info for easy extraction
         */
        printf("++%6d+%6d+%6d+%6d+%6d+%6d+%6d+%6d+%6d+%6d+\n", lo_cnts[0],
                lo_cnts[1], lo_cnts[2], lo_cnts[3], lo_cnts[4], hi_cnts[0],
                hi_cnts[1], hi_cnts[2], hi_cnts[3], hi_cnts[4]);
        printf("++%5.2f+%5.2f+%5.2f+%6d+%5d+%5d+%5d+%5d+%5d+%5d+%5d+%5d+\n",
                (float) (totalpts - (pos1cnt + neg1cnt)) / totalpts * 100.0,
                (float) (totalpts - (pos2cnt + neg2cnt)) / totalpts * 100.0,
                (float) (totalpts - (pos3cnt + neg3cnt)) / totalpts * 100.0,
                rtmiss, z_st[0], z_en[0], z_st[1], z_en[1], z_st[2], z_en[2],
                z_st[3], z_en[3]);
        printf("++%5.2f+%5.2f+%5.2f+%5.2f+%5.2f+%5.2f+%5d+%5d+%6.2f+%6.2f+\n",
                (float) neg1cnt / totalpts * 100.0, (float) neg2cnt / totalpts * 100.0,
                (float) neg3cnt / totalpts * 100.0, (float) pos1cnt / totalpts * 100.0,
                (float) pos2cnt / totalpts * 100.0, (float) pos3cnt / totalpts * 100.0,
                max_lin, min_lin, max_lat, min_lat);
        printf(
                "---------------------------------------------------------------\n");
        printf("\n\n");

    }

    /* set exit value to reflect threshold checks */

    if (failed) {
        printf("Threshold Status: FAILED\n");
        return (FAILURE);
    } else {
        printf("Threshold Status: SUCCESS\n");
        return (SUCCESS);
    }

} /* main */

/*****************************************************************
 *
 * check user parameters and show example if not all parms given
 *
 *****************************************************************/

int8 check_usage(int argc, char *argv[]) {
    if (argc < 10) {
        printf("\n\nUsage:\n");
        printf("\t%s <nrtfile><clmfile><outfile><monthstr>\n", argv[0]);
        printf("\t <param><t1><t2><t3><maxmiss><loval><hival><outmax>\n");
        printf("\t <min_lat><max_lat>[diff][type]\n");
        printf("\nWhere:\n");
        printf("\tnrtfile:   Real-time file to process\n");
        printf("\tclmfile:   Climatology file\n");
        printf("\toutfile:   output QC file\n");
        printf("\tmonthstr:  Month in string form (i.e., January)\n");
        printf("\t           (only used for examining climatology file,\n\t           otherwise derived from the input file name\n");
        printf("\tparam:     Parameter to run: OZONE ...only\n");
        printf(
                "\tt1:        Threshold percent of pts w/in 1 STD DEV of climatology\n");
        printf(
                "\tt2:        Threshold percent of pts w/in 2 STD DEV of climatology\n");
        printf(
                "\tt3:        Threshold percent of pts w/in 3 STD DEV of climatology\n");
        printf("\tmaxmiss:   Max num of missing NRT points permitted.\n");
        printf("\tloval:     lowest acceptable value in NRT file\n");
        printf("\thival:     highest acceptable value in NRT file\n");
        printf("\toutmax:    # of points, weighted by cos( latitude )\n");
        printf("\t           that can be outside lo/hival before error set\n");
        printf("\tmin_lat:   minimum latitude to consider (default -90.).\n");
        printf("\tmax_lat:   maximum latitude to consider (default 90.).\n");
        printf("\t           NOTE that if min_lat = 95, the limits will be\n");
        printf("\t           determined using the julian day in the equation:\n");
        printf("\t           min_lat = -90. + zen + max_lat\n");
        printf("\t           max_lat =  90. - zen - max_lat\n");
        printf("\t           where max_lat = limit adjustment, default, 10.\n");
        printf("\t           and zen = -23. * cos( ( jday + 10 ) * 360 / 365 )\n");
        printf("\t        ** [optional] parameters follow.  Preceeding args \n");
        printf(
                "\t           are required for subsequently ones used (fill spaces).\n");
        printf(
                "\tdiff:      [optional] Enter 'd' to see a simple difference output\n");
        printf("\t           grid rather than the default STD variance 's'.\n");
        printf("\t           Difference calculations do not product reports.\n");
        printf(
                "\ttype       [optional] Enter 1 if NRT file is actually a CLM file.\n");
        printf("\t           This will allow the program to read the CLM\n");
        printf("\t           as of test of CLM to itself.\n");
        printf("\n");
        printf("Example: \n\n");

        printf("\t o3qc $SDSDEMO/S199407100_TOVS.OZONE.hdf \n");
        printf("\t $SDSDATA/S19891991_TOMS.OZONE.hdf output.hdf \n");
        printf("\t March OZONE 25 50 75 50\n\n");
        return (ERROR);
    }
    return 0;
}
