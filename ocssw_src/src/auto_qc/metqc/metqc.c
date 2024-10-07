/*****************************************************************
 * FILE:        metqc.c
 *
 * PURPOSE:     QC HDF ancillary real time datafile against climatology
 *              for one of four parameters.
 *
 * DESCRIPTION: This program will create a HDF with:
 *
 *              HDF Geometry Vdata of Lat/Lon limits and steps.
 *              2-D std dev of NRT vs CLM for z_wind or,
 *              2-D std dev of NRT vs CLM for m_wind or,
 *              2-D std dev of NRT vs CLM for press  or,
 *              2-D std dev of NRT vs CLM for rel_hum.
 *
 *              A ASCII text file of %points in each STD category.
 *
 *  For MET data:
 *    The real time parameter is read into array
 *       z_windRT, m_windRT, pressRT or rel_humRT
 *
 *    Real-time files are dimensioned 73x144.  The COADS based
 *    climatology is 90x180.  A regridded COADS climatology
 *    (S19461990_COADS.REGRID.MET) is therefore used by this program.
 *
 *    The climatological averages and standard deviations for the same time
 *    period (month) are read into:
 *       z_windCLMAVG, m_windCLMAVG, pressCLMAVG or rel_humCLMAVG and
 *       z_windCLMSTD, m_windCLMAVG, pressCLMSTD or rel_humCLMSTD respectively.

 *
 * PROCEDURE:
 *
 *    The CLMAVG array is subtracted from the RT array and
 *    divided by the CLMSTD array to produce an array of the number
 *    of std deviations from the mean for the real-time data.
 *
 *    An option for an array of just the differences can be run by envoking
 *    the 'd' (difference) flag on the command line.
 *
 *    As the data is processed a tally of points falling + and - 1, 2, and 3
 *    STDs from the mean will be stored for later writing of percentages in
 *    each category (if 'difference' option not used).
 *
 * INPUT PARMS:  command line arguments (argv[?])
 *   [1] char *nrtfile  - real-time input file to process
 *                  ("S199310200_NMC.MET")
 *   [2] char *climfile - regridded climatology input file to process
 *                  ("S19461990_COADS.REGRID.MET")
 *   [3] char *outfile - output QC filename.
 *       Note: FNOC could be used in place of NMC.
 *   [4] char monthstr  - month being processed
 *   [5] - field to compare: 'z_wind', 'm_wind', 'press', 'rel_hum' or 'p_water'
 *   [6] float32 stdlimit1 - % w/in 1 STD
 *   [7] float32 stdlimit2 - % w/in 2 STD
 *   [8] float32 stdlimit3 - % w/in 3 STD
 *   [9] int maxmissing    - max num missing pts.
 *
 *  Optional inputs:  on command line
 *   (args are read in order so preceeding optionals are required in order
 *    to use succeeding ones.)
 *
 *   [10] float32 loval  - lowest allowed value
 *   [11] float32 hival  - highest allowed value
 *   [12] int outmax     - max num allowed outliers before error flag set.
 *   [13] char gridflag  - 'S'tandard deviation [default] or 'D'ifference calc.
 *                        Note, difference calculations do not product reports.
 *   [14] int ifileflag  -  0 = file is NRT file [default],
 *                      1 = file in CLM file (stats on one CLM vs another)
 *
 * OUTPUT PARMS:
 *    prints report to standard out.
 *    returns a SUCCESS (0) or FAILURE(5) value as Standard Deviation
 *    threshold exit status.  This value can be checked with a C-shell
 *    $status flag.
 *
 *    output name is similar to the input real-time file.
 *    Example name: "S199310200_NMC.MET.'parm'.qc"
 *                 ('parm' = z_wind, m_wind, press, rel_hum or p_water)
 *
 * RETURNS:     Program returns success or failure codes to a calling
 *              'shell' script.
 *
 * LOCAL VARS:  variables for storing HDFs.
 *
 * Subs called: int rdsds        - reads an HDF SDS data grid
 *              int8 check_usage - confirm args
 *              int startHDF     - open HDF file for output
 *              int setupGrid    - setup HDF geometry struct
 *              int writeGeom    - write HDF geometry struct
 *              int wrtsds       - write HDF SDS array
 *              int addAttr      - add HDF SDS attributes
 *              int setSDSref    - set HDF SDS reference
 *              void deattachHDFgrid - detach HDF grid struct
 *              int closeHDFstructs - close HDF files
 *              int rd_size - read field size
 *              int resize_2d - resize climatology field to match nrt
 *
 * HISTORY:     none
 *
 * NOTE:        This program is meant to generate HDF QC SDS data sets.
 *
 * AUTHOR:      Brian D. Schieber, GSC, 4/93
 *
 * MODIFICATION HISTORY:
 *
 * 10/14/93 BDS - Subtantial mods to process all three MET parameters
 *                at once and to handle the NEW HDF format specs.
 * 12/06/93 BDS - Make that all 4 parms (z_wind,m_wind,press,rel_hum).
 *                Verifying the exclusion of missings from CLM in
 *                these calcs (in NRT, 0.0's passed through).
 *                Defined VALIDMIN to -99.1 (catches CDFMINVAL (-99.9)
 *                in CLM) as min threshold for QC valid points to
 *                calculate.
 * 3/94 BDS     - prep for call from script which used SYBASE calls.
 *                Append parm type to 'qc' file name.
 *                Return pass or fail STD DEV thresholds as file status.
 * 5/94 BDS     - allow specification of output directory
 * 8/14/95 BDS  - revised for OAPS spec 2.7
 *   11/95 BDS, add ability to sum number of values in image < or >
 *              two specified range arguments and return FAIL if > outmax.
 *              Values for LO(lower range), HI (upper range),
 *              and OUTMAX (MAX. Pts allowed outside range) hardcoded within
 *              for each of the four parameters.
 *              Also, print min and max value of data.
 * 8/21/96 BDS  - renamed 'perror' to 'pexit' to avoid HDF4.0 conflict.
 * W. Robinson, GSC, 6 Feb 97  add p_water and make code work independently
 *              of the nrt and climatology sizes
 *****************************************************************/
#include <mfhdf.h>

#include "ancil.h"
#include "l1io.h"
#include <genutils.h>

/*
 *  Note on all data sizes:  old (pre march 97) data sizes
 *  met nrt: 73 lines, 144 pixels
 *  met climatology: 90 lines, 180 pixels
 *  met climatology gridded: 73 lines, 144 pixels - same as nrt
 */
/*
 * MET specific settings
 */

#define VGROUPCLASS "PlanetaryGrid" /* all inner Vgroups use same name */

#define BIN_METH     2                 /* geometry structure settings    */
#define REGISTRATION CENTER            /*             "                  */

/* MET file settings */
#define VSIZE         2.5              /*             "                  */
#define HSIZE         2.5              /*             "                  */
#define MAX_NORTH    89.0              /*             "                  */
#define MAX_SOUTH   -89.0              /*             "                  */
#define MAX_WEST   -179.0              /*             "                  */
#define MAX_EAST    179.0              /*             "                  */

#define VALIDMIN    -99.1              /* using '.1' to compensate float */
#define CLMfile       1

struct annotation *annot;

int main(int argc, char *argv[]) {
    int i, j, p;
    int rank;
    int result = 0;
    int missing = 0;

    float32 minval = 1200.0; /* initial minimum value */
    float32 maxval = -100.0; /* initial maximum value */

    float32 loZ_WIND = -30.0; /* default lo range threshold */
    float32 hiZ_WIND = 30.0; /* default hi range threshold */
    float32 loM_WIND = -30.0; /* default lo range threshold */
    float32 hiM_WIND = 30.0; /* default hi range threshold */
    float32 loPRESS = 850.0; /* default lo range threshold */
    float32 hiPRESS = 1100.0; /* default hi range threshold */
    float32 loREL_HUM = 5.0; /* default lo range threshold */
    float32 hiREL_HUM = 100.0; /* default hi range threshold */
    float32 loP_WATER = 0.; /* as above precip. water */
    float32 hiP_WATER = 200.;
    float32 loval; /* default actual range threshold */
    float32 hival; /* default actual range threshold */
    int outmax = 25; /* max num pts allowed to fall
                                       outside range */

    int locnt; /* pts under lo range limit */
    int hicnt; /* pts above hi range limit */
    int rtmiss = 0;
    int totalpts = 0;
    int neg1cnt = 0;
    int neg2cnt = 0;
    int neg3cnt = 0;
    int pos1cnt = 0;
    int pos2cnt = 0;
    int pos3cnt = 0;
    int failed = 0;
    int array_size = 0;
    int32 shape_nrt[2], shape_clm[2];
    int32 inShape[2];
    char infileRT[MAXNAMELNG];
    char infileCLM[MAXNAMELNG];
    char outfile[MAXNAMELNG];
    char gridflag[2];
    char datalabel[MAXNAMELNG];
    int32 datatype;
    char *dataattr;
    char *dataunit;
    char vgname[MAXNAMELNG];
    char monthstr[13];
    char sdsname[MAXNAMELNG];
    char vgroupname[21];
    int32 sdsid, sdfid, fid, gridid, geomid;
    int ifileflag = 0;
    float32 stdlimit1, stdlimit2, stdlimit3;
    int maxmissing;
    /*
     * data type array pointers
     */

    float32 *float_SDSdataRT;
    float32 *float_SDSdataCLMAVG;
    float32 *float_SDSdataCLMSTD;
    float32 *float_SDSdataQC;
    float32 *c_avg, *c_std; /* new avg, std dev for pre-resized climatology */

    /* functions used from this file */

    int8 check_usage();

    int resize_2d(float *, int, int, int, int, float *);
    int rd_size(char *, int *, int *);

    /*
     * ------- check command line arguments and set args  ------------------
     */

    strcpy(gridflag, "s"); /* says: use std as default output QC field */
    if ((check_usage(argc, argv)) != 0) pexit("insufficient args provided");

    strcpy(infileRT, argv[1]);
    strcpy(infileCLM, argv[2]);
    strcpy(outfile, argv[3]);
    strcpy(monthstr, argv[4]);

    if (!strcmp(lowcase(argv[5]), "ozone"))
        pexit("OZONE parameter type not supported by this program.");

    if (!strcmp(lowcase(argv[5]), "z_wind")) {
        loval = loZ_WIND;
        hival = hiZ_WIND;
    } else if (!strcmp(lowcase(argv[5]), "m_wind")) {
        loval = loM_WIND;
        hival = hiM_WIND;
    } else if (!strcmp(lowcase(argv[5]), "press")) {
        loval = loPRESS;
        hival = hiPRESS;
    } else if (!strcmp(lowcase(argv[5]), "rel_hum")) {
        loval = loREL_HUM;
        hival = hiREL_HUM;
    } else if (!strcmp(lowcase(argv[5]), "p_water")) {
        loval = loP_WATER;
        hival = hiP_WATER;
    } else
        pexit(
            "parameter must be 'z_wind', 'm_wind', 'press', 'rel_hum' or 'p_water'.");

    strcpy(vgroupname, argv[5]);

    sscanf(argv[6], "%f", &stdlimit1);
    sscanf(argv[7], "%f", &stdlimit2);
    sscanf(argv[8], "%f", &stdlimit3);
    sscanf(argv[9], "%d", &maxmissing);
    if (argc >= 11) sscanf(argv[10], "%f", &loval);
    if (argc >= 12) sscanf(argv[11], "%f", &hival);
    if (argc >= 13) sscanf(argv[12], "%d", &outmax);
    if (argc >= 14) {
        strncpy(gridflag, argv[14], 1);
        gridflag[1] = '\0';
    }
    if (argc >= 15) ifileflag = atoi(argv[15]);


    /*
     * ------- Read Z/M_WIND, PRESS, REL_HUM or OZONE real-time array files ---------
     */

    if (ifileflag == CLMfile) {
        strcpy(vgname, monthstr);
        strcpy(sdsname, vgroupname);
        strcat(sdsname, "_mean");
    } else {
        strcpy(vgname, "Geophysical Data");
        strcpy(sdsname, vgroupname);
    }
    /*
     *  Read the size of the NRT field and set the size here
     */
    if (rd_size(infileRT, (int*) &shape_nrt[1], (int*) &shape_nrt[0]) != 0) {
        pexit("metqc: Unable to read the data field size");
    }
    /*
     *  allocate the nrt-size arrays here
     * RT - met field
     * CLMAVG - mean climatology
     * CLMSTD - Standard deviation of observations around the mean
     * QC - output QC field
     */
    rank = 2;
    array_size = shape_nrt[0] * shape_nrt[1];

    if ((float_SDSdataRT =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSdataRT");

    if ((float_SDSdataCLMAVG =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSdataCLMAVG");

    if ((float_SDSdataCLMSTD =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSdataCLMSTD");

    if ((float_SDSdataQC =
            (float32 *) malloc(sizeof (float32) * array_size)) == NULL)
        pexit("malloc float_SDSdataQC");
    /*
     *  Read Z/M_WIND, PRESS, REL_HUM, P_WATER or OZONE real-time array files 
     */

    if ((rdsds(infileRT, vgname, sdsname, inShape,
            float_SDSdataRT)) != 0) pexit("rdsds RT");

    if ((inShape[0] != shape_nrt[0]) || (inShape[1] != shape_nrt[1])) {
        printf("inShape[0] %d shape_nrt[0] %d inShape[1] %d shape_nrt[1] %d\n",
                inShape[0], shape_nrt[0], inShape[1], shape_nrt[1]);
        pexit("real-time dimensions not matching expected");
    }

    /*
     * ------- Determine appropriate climatology month (SDS) to match ------
     *         real time file.  Get the average and std.dev. SDSs.
     */

    /*
     *  Read the size of the CLM field and set the size here
     */
    if (rd_size(infileCLM, (int*) &shape_clm[1], (int*) &shape_clm[0]) != 0) {
        pexit("metqc: Unable to read the climatology data field size");
    }
    /*
     *  allocate the space
     */
    if ((c_avg =
            (float32 *) malloc(sizeof (float32) * shape_clm[0] * shape_clm[1]))
            == NULL) {
        pexit("malloc float_SDSdataCLMAVG");
    }
    if ((c_std =
            (float32 *) malloc(sizeof (float32) * shape_clm[0] * shape_clm[1]))
            == NULL) {
        pexit("malloc float_SDSdataCLMAVG");
    }

    strcpy(vgname, monthstr);
    strcpy(sdsname, vgroupname);
    strcat(sdsname, "_mean");

    if ((rdsds(infileCLM, vgname, sdsname, inShape,
            c_avg)) != 0) pexit("rdsds CLMAVG");

    if ((inShape[0] != shape_clm[0]) || (inShape[1] != shape_clm[1])) {
        printf("inShape[0] %d shape_clm[0] %d inShape[1] %d shape_clm[1] %d\n",
                inShape[0], shape_clm[0], inShape[1], shape_clm[1]);
        pexit("clim avg dimensions not matching expected");
    }

    strcpy(sdsname, vgroupname);
    strcat(sdsname, "_std_dev");

    if ((rdsds(infileCLM, vgname, sdsname, inShape,
            c_std)) != 0) pexit("rdsds CLMSTD");

    if ((inShape[0] != shape_clm[0]) || (inShape[1] != shape_clm[1])) {
        printf("inShape[0] %d shape_clm[0] %d inShape[1] %d shape_clm[1] %d\n",
                inShape[0], shape_clm[0], inShape[1], shape_clm[1]);
        pexit("clim std dimensions not matching expected");
    }

    /*
       l = 0;
       for (i = 0; i < shape_nrt[0]; i++) {
          printf ("\n");
          for (j = 0; j < shape_nrt[1]; j++) {
             printf ("%f ", float_SDSdataRT[l]);
             l++;
          }
       }
     */

    /*
     *  reconcile the climatology array size to that of the met field
     */
    printf("NRT field size is - pixels: %d, lines: %d\n",
            shape_nrt[1], shape_nrt[0]);
    if (shape_clm[0] != shape_nrt[0] || shape_clm[1] != shape_nrt[1]) {
        printf("NOTE: Climatology field size is different\n");
        printf("      Pixels: %d, lines: %d\n", shape_clm[1], shape_clm[0]);
        printf("      Resizing climatology to match NRT\n");
    }

    if ((shape_clm[0] <= 0) || (shape_clm[1] <= 0) ||
            (shape_nrt[0] <= 0) || (shape_nrt[1] <= 0)) {
        pexit("metqc: one of shape_clm, shape_nrt dimensions <= 0!");
    }
    if (resize_2d(c_avg, shape_clm[1], shape_clm[0],
            shape_nrt[1], shape_nrt[0], float_SDSdataCLMAVG) != 0) {
        pexit("metqc: resize_2d problem for CLMAVG, exiting");
    }

    if (resize_2d(c_std, shape_clm[1], shape_clm[0],
            shape_nrt[1], shape_nrt[0], float_SDSdataCLMSTD) != 0) {
        pexit("metqc: resize_2d problem for CLMSTD, exiting");
    }
    free(c_avg);
    free(c_std);

    /*
     * ------- Compute variance of real-time from climatology -----------
     *         and tally #0s (and/or missing), #+ & - 1,2,3 std dev
     *         for writing to file for mailing.
     */

    p = 0;
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

    for (i = 0; i < shape_nrt[0]; i++) {
        for (j = 0; j < shape_nrt[1]; j++) {
            float_SDSdataQC[p] = 0.0; /* initialize each point */

            /* pt. is missing, NRT missing, a MIN val, and counted as out of
               lower range */
            if ((float_SDSdataRT[p] <= VALIDMIN) ||
                    (float_SDSdataCLMAVG[p] <= VALIDMIN)) {
                missing++;
                if (float_SDSdataRT[p] <= VALIDMIN) rtmiss++;
                if (float_SDSdataRT[p] < minval) minval = float_SDSdataRT[p];
                if (float_SDSdataRT[p] < loval) locnt++;
            } else {
                if ((!strcmp(gridflag, "d")) || (!strcmp(gridflag, "D"))) {
                    float_SDSdataQC[p] =
                            (float32) float_SDSdataRT[p] - float_SDSdataCLMAVG[p];
                } else {
                    if (float_SDSdataCLMSTD[p] > 0)
                        float_SDSdataQC[p] = (float32)
                        (float_SDSdataRT[p] - float_SDSdataCLMAVG[p]) /
                        float_SDSdataCLMSTD[p];
                }

                /* find min and max of all points */
                if (float_SDSdataRT[p] < minval) minval = float_SDSdataRT[p];
                if (float_SDSdataRT[p] > maxval) maxval = float_SDSdataRT[p];

                /* sum up points outside lo/hi range (LO adds to missing counts)*/
                if (float_SDSdataRT[p] < loval) locnt++;
                if (float_SDSdataRT[p] > hival) hicnt++;

                if (float_SDSdataQC[p] < -1.0) neg1cnt++;
                if (float_SDSdataQC[p] < -2.0) neg2cnt++;
                if (float_SDSdataQC[p] < -3.0) neg3cnt++;

                if (float_SDSdataQC[p] > 1.0) pos1cnt++;
                if (float_SDSdataQC[p] > 2.0) pos2cnt++;
                if (float_SDSdataQC[p] > 3.0) pos3cnt++;
            }
            /*
                     if ((j > 100) && (j < 150)) {
                     printf("QC  [%d][%d]: %f\n", i,j,float_SDSdataQC[p]);
                     printf("RT  [%d][%d]: %f\n", i,j,float_SDSdataRT[p]);
                     printf("AVG [%d][%d]: %f\n", i,j,float_SDSdataCLMAVG[p]);
                     printf("STD [%d][%d]: %f\n", i,j,float_SDSdataCLMSTD[p]);
                     }
                     printf("[%d][%d]: %f %f %f\n", i,j,float_SDSdataRT[p],
                           float_SDSdataCLMAVG[p], float_SDSdataCLMSTD[p]);
             */
            p++;
        } /* for j */
    } /* for i */

    totalpts = array_size - missing;

    free(float_SDSdataRT);
    free(float_SDSdataCLMAVG);
    free(float_SDSdataCLMSTD);

#if 0
    p = 0;
    for (i = 0; i < shape_nrt[0]; i++) {
        for (j = 0; j < shape_nrt[1]; j++) {
            float_SDSdataQC[p] = 0.0; /* initialize each point */
            p++;
        } /* for j */
    } /* for i */
#endif

    /*
     * ----------------- Write QC SDS ----------------------------
     */

    /*
     * Create HDF file
     */


    if ((startHDF(outfile, &sdfid, &fid, DFACC_CREATE)) != 0)
        pexit("Fatal error starting output HDF file");

    /*
     * Create grid structure
     */

    /* gridid = setupGrid(fid, VGROUPCLASS, vgroupname); */
    gridid = setupGrid(fid, vgroupname);

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
    datatype = DFNT_FLOAT;
    dataattr = "UNITS";
    dataunit = "std dev diff";

    if ((sdsid = wrtsds(sdfid, rank, shape_nrt, datatype, datalabel,
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

    if ((result = setSDSref(sdsid, gridid)) != 0) pexit("setSDSref");

    /*
     * deattach HDF grid
     */

    deattachHDFgrid(gridid);

    /*
     * close HDF structures
     */

    if ((result = closeHDFstructs(sdfid, fid)) != 0) pexit("closeHDFstructs");

    /*
     * ----------------- Print stats to standard output -----------
     */

    /* printf("gridflag [%s]\n",gridflag); */

    if (!strcmp(gridflag, "s")) {
        printf("\n");
        printf("-------------------------------------------------------------\n");
        printf("Results of comparison of real-time and climatological files:\n");
        printf("%s\n%s\n", infileRT, infileCLM);
        printf("Month: %s\t\tParameter: %s\n", monthstr, argv[5]);
        printf("Thresholds: %8.3f %8.3f %8.3f Max Missings: %d",
                stdlimit1, stdlimit2, stdlimit3, maxmissing);
        printf("\n");
        printf("Minimum value: %8.3f Maximum: %8.3f\n", minval, maxval);

        printf("Total # non-missing values: %6d\n", totalpts);
        printf("Total # values <%8.3f: %d  >%8.3f: %d  allowed: %d",
                loval, locnt, hival, hicnt, outmax);

        /* if (locnt >= outmax) failed = 1;  lo points out of range exceeded */
        /* if (hicnt >= outmax) failed = 1;  hi points out of range exceeded */
        /*
         *  WDR if fail in either case, append a set of astrisks
         */
        if (locnt >= outmax || hicnt >= outmax) {
            failed = 1;
            printf("  ***\n");
        } else {
            printf("\n");
        }

        printf("\n");
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

        printf("Total # missing values:             %6d\n", missing);
        printf("Total # missing real-time:          %6d", rtmiss);
        if (rtmiss > maxmissing) {
            failed = 1;
            printf("  ***\n");
        } else {
            printf("\n");
        }
        printf("-------------------------------------------------------------\n");
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
    if (argc < 9) {
        printf("\n\nUsage:\n");
        printf("\t%s <file><file><outfile><monthstr><param><t1><t2><t2><maxmiss><lo><hi><outmax>[diff][type] \n", argv[0]);
        printf("\nWhere:\n");
        printf("\tfile:      Real-time file to process\n");
        printf("\tfile:      Climatology file\n");
        printf("\toutfile:   Output QC file\n");
        printf("\tmonthstr:  Month in string form (i.e., January)\n");
        printf("\tparam:     Parameter to run: z_wind m_wind press rel_hum p_water\n");
        printf("\tt1:        Threshold percent of pts w/in 1 STD DEV of clim\n");
        printf("\tt2:        Threshold percent of pts w/in 2 STD DEV of clim\n");
        printf("\tt3:        Threshold percent of pts w/in 3 STD DEV of clim\n");
        printf("\tmaxmiss:   Max num of missing NRT points permitted\n");
        printf("\tlo:        Lowest allowed value in NRT file\n");
        printf("\thi:        Highest allowed value in NRT file\n");
        printf("\tlo:        Maximum # of NRT points allowed outside of lo/hi before returning an error flag\n");
        printf("\t           [optional] parameters follow.  Preceeding args \n");
        printf("\t           are required for subsequent ones used\n");
        printf("\tdiff:      [optional] Enter 'd' to see a simple difference output\n");
        printf("\t           grid rather than the default STD variance 's'\n");
        printf("\t           Difference calculations do not product reports\n");
        printf("\ttype       [optional] Enter 1 if NRT file is actually a CLM file\n");
        printf("\t           This will allow the program to read the CLM\n");
        printf("\t           as of test of CLM to itself\n");
        printf("\n");
        return (ERROR);
    }
    return 0;
}
