/*****************************************************************
 * File:        o3nrt_omi.c
 *
 * Purpose:     create HDF ancillary real time datafile from
 *              OMI TOMS ascii input only - from o3nrt
 *
 * NOTE:        Due to the non-ovelapping nature of the EPTOMS product
 *              2 days worth of data are required to fill all pixels in
 *              an image.   All other products cover the globe with one
 *              days's worth of data.
 *
 * Description: this program will create a HDF with the components:
 *
 *              HDF DAAC compliant global attributes
 *              Geometry Vdata for storage of LAT/LON endpoints and step
 *              2-D data parameter OZONE SDS
 *              2-D OZONE QC flag (padded all zeroes for TOMS and TOVS,
 *                with EARLIERVAL to EPTOMS missing points filled with
 *                prior day's data)
 *
 * Input parms:
 * char *annotfile - input file of DAAC compliant annots ("fillenv.eptomsnrt")
 * char *infile    - input ozone file to process ("ga960903.ept")
 * char *directory - output file directory (".")
 * char *itype     - "OMITOMS"
 * char *prevfile  - input ozone file for previous day ("OMITOMS" only)
 *
 * Returns:     false on error
 *
 * Subs called:
 *              fillenv_annot - fill an array with metadata read from
 *                              an ASCII file.
 *              count_annot- determine number of entries in ASCII annot file.
 *              rdgrid_    - read TOMS/TOVS ASCII grid (FORTRAN routine)
 *              gregor_    - calc gregorian date from julian number.
 *              startHDF   - open HDF file structures
 *              setupGrid  - open HDF grid struct.
 *              wrtattr    - write HDF attributes to file
 *              wrtsds     - write SDS from array to file
 *              pexit      - print fatal errors
 *              check_usage- check command line args
 *              closeHDFstructs- close HDF file structs
 *              writeGeom  - write HDF file geometry (lat/lon)
 *              addAttr    - set file local or global attributes
 *              setSDSref  - get SDS ref and set to Vgroup
 *
 * History:     none
 *
 * Note:        This program is meant to generate SeaWiFS style HDF
 *              ancillary real-time datafiles.
 *              The string parsing used by rdgrid may need revised to
 *              match any new TOMS/TOVS/EPTOMS file formats.
 *
 * Author:      Brian D. Schieber, GSC, 3/93
 *
 *
 * Modification history:
 *   BDS, 4/22/93 - mod int32 to int16 (per L. Kumar, Jimf, Darzi)
 *   BDS, 6/7/93  - program now generates output filenames consistent
 *                  with SeaWiFS specs.
 *   BDS, 8/30/93 - redesign for new HDF specs.  Refer to
 *                  "SeaWiFS ANCILLARY DATA FORMAT (HDF) UPDATE"
 *                  paper (SeaWiFS tech. memo) for details.
 *   BDS, 5/17/95 - Comply with mods in "Operational Archive Prod. Specs."
 *                      Version 2.7 (2 May 1995) Many changes.
 *   BDS, 10/1/95 - Comply with 2.8 specs:
 *                   With this modification, TOVS and TOMS input ASCII
 *                   files diverge:
 *                    1. TOVS input now contains an extra line of header to
 *                      describe the starting and ending times of the data.
 *                    2. Filenames of HDF files made from TOVS input files
 *                      now list the starting AND ending times of the input.
 *   BDS, 12/8/95 - change metadata "Process Control" to "Processing Control"
 *                   - append terminators to some metadata strings.
 *   BDS, 8/21/96 - renamed 'perror' to 'pexit' to avoid HDF4.0 conflict.
 *   BDS, 9/4/96 - modifications to support new EPTOMS files.
 *   BDS, 9/12/96 - mod to add Ascention time metadata.  Also support all
 *                  three ozone sources: TOMS, TOVS, EPTOMS
 *   BDS, 9/19/96 - also support ADTOMS (not [optionally] using 2nd file as in EPTOMS)
 *                  Add M.Darzi's new flags for type of EPTOMS data fill:
 *                    EPTOMS values for points in "QC grid" are:
 *                    0  - no change to original value.
 *                    1  - value changed by QC routine (not this program)
 *                    10 - value from earlier source product.
 *                    11 - value from nearest non-zero grid point.
 *  W. Robinson, GSC, 3 Jul 97  call a smoothing routine: world_avg for the 
 *                  EP TOMS data.  Also, set start time of EP TOMS to be start
 *                  time of the primary file, do not use start of second
 *                  file used for filler.
 *  W. Robinson, SAIC, 15 Jun 09  Quick mod to read the OMI with 360 pix, 
 *                    180 lines
 *****************************************************************/
#include <mfhdf.h>
#include "ancil.h"
#include <time.h>
#include <unistd.h>
#include "ancnrt_proto.h"

/*
 * Ozone specific settings
 */

#define VGROUPNAME   "Geophysical Data"
#define BIN_METH     2
#define REGISTRATION CENTER
#define EARLIERVAL    10  /* QC grid value for EPTOMS pts filled from prev day */
#define NEARESTVAL    11  /* QC grid value for EPTOMS pts filled from nearest
                             (file 1 or 2) non-zero point */

int main(int argc, char *argv[]) {
    char *datafmt;
    char *datalabel;
    char dataunit[MAXNAMELNG];
    char dataattr[MAXNAMELNG];
    char directory[MAXNAMELNG];
    char itype[MAXNAMELNG];
    char annotfile[MAXNAMELNG];
    char outfile[MAXNAMELNG];
    char outfilename[MAXNAMELNG];
    char infile[MAXNAMELNG];
    char prevfile[MAXNAMELNG];
    char start_time[17];
    char end_time[17];
    char node_time[17];
    char acntime[9];
    int i, j, k, l, kk;
    int rank;
    int type;
    int result;
    int array_size = 0;
    int numannarr = 0;
    int modified_pts = 0; /* EPTOMS, track pts modded */
    int i_skip; /* integer dummy var */
    char tmpvar[MAXNAMELNG]; /* temp time variable */
    int year, julval, hour, min; /* data times */
    short startyear, startday, starthour, /* start data times */
            startmin, startsec;
    int startmsec;
    short endyear, endday, endhour, /* end data times */
            endmin, endsec;
    int endmsec;
    int startmsecday; /* millisec of day */
    int endmsecday; /* millisec of day */
    int pyear = 0; /* local proc times */
    float latstep, lonstep; /* lat/lon incr. */
    float nmostlat, smostlat, wmostlon, emostlon; /* lat/lon corners */
    int latsz, lonsz; /* lat/lon size */
    int16 datarr[180][360];
    int16 prevdatarr[180][360];
    int32 shape[2];
    time_t t; /* processing time */
    struct tm *localtm; /* processing time */
    struct annotation *xannot;

    /*
     * HDF datafile variables
     */

    int32 sdfid, fid, gridid;
    int32 datatype;

    /*
     * data type array pointers
     */

    int8 *int8_SDSdataQC;
    int16 *int16_SDSdataO3;

    /*
     * functions used
     */

    int wrtattr();
    int32 wrtsds();
    int startHDF();
    int setSDSref();
    int closeHDFstructs();
    void pexit();
    int8 check_usage();
    struct annotation * fillenv_annot();
    char *upcase();
    void world_avg(int16 *, int32, int32, int, int, int16 *);

    /*
     * ------- check command line arguments and set args  ------------------
     */

    if ((check_usage(argc, argv)) != 0) exit(-1);

    strcpy(annotfile, argv[1]);
    strcpy(infile, argv[2]);
    strcpy(directory, argv[3]);
    strcpy(itype, argv[4]);

    /*
     * ------- Read ozone ASCII file and extract data grid ----------
     */

    if (!strcmp(upcase(itype), "OMITOMS")) {
        type = 2;
        if (argc == 6) strcpy(prevfile, argv[5]);
        else strcpy(prevfile, ""); /* no 2nd file */
    } else {
        pexit("ERROR: OMITOMS not specified");
    }

    /*
     * Read ASCII OZONE grid
     */

    rdgrid_(infile, &type, (int16*) datarr, &julval, &year, start_time, end_time,
            acntime, &result, 255, 17, 17, 9);
    if (result == 1) pexit("rdgrid_");

    /*
     * ------- Create component date/time ----------
     *   TOVS input ASCII data provides start and end times in an additional
     *   header line.
     *   TOMS, EPTOMS and ADTOMS input ASCII does not provide start and end time,
     *   only one header line for year, Gregorian day, and Julian day.
     */

    if (type == 0) { /* TOMS data */

        gregor_(&julval, &year, &i_skip, &i_skip); /* get year from julian */

        hour = 12;
        sprintf(outfile, "%s/S%04d%03d%02d_%s.OZONE", directory, year,
                julval, hour, "TOMS");
        sprintf(outfilename, "S%04d%03d%02d_%s.OZONE", year, julval,
                hour, "TOMS");

        /* Make hour, min, sec, ms all zeros for TOMS */

        sprintf(start_time, "%04d%03d%02d0000000", year, julval, hour);
        sprintf(end_time, "%04d%03d%02d0000000", year, julval, hour);

        startyear = year;
        endyear = year;
        startday = julval;
        endday = julval;
        startmsecday = endmsecday = 0;
    } else if (type == 1) { /* TOVS */

        strcpy(&start_time[16], "\0");
        strcpy(&end_time[16], "\0");

        strncpy(tmpvar, &start_time[0], 4);
        strcpy(&tmpvar[4], "\0");
        startyear = (short) atoi(tmpvar);

        strncpy(tmpvar, &start_time[4], 3);
        strcpy(&tmpvar[3], "\0");
        startday = (short) atoi(tmpvar);

        strncpy(tmpvar, &start_time[7], 2);
        strcpy(&tmpvar[2], "\0");
        starthour = (short) atoi(tmpvar);

        strncpy(tmpvar, &start_time[9], 2);
        strcpy(&tmpvar[2], "\0");
        startmin = (short) atoi(tmpvar);

        strncpy(tmpvar, &start_time[11], 2);
        strcpy(&tmpvar[2], "\0");
        startsec = (short) atoi(tmpvar);

        strncpy(tmpvar, &start_time[13], 3);
        strcpy(&tmpvar[3], "\0");
        startmsec = atoi(tmpvar);

        /* ******************************** */

        strncpy(tmpvar, &end_time[0], 4);
        strcpy(&tmpvar[4], "\0");
        endyear = (short) atoi(tmpvar);

        strncpy(tmpvar, &end_time[4], 3);
        strcpy(&tmpvar[3], "\0");
        endday = (short) atoi(tmpvar);

        strncpy(tmpvar, &end_time[7], 2);
        strcpy(&tmpvar[2], "\0");
        endhour = (short) atoi(tmpvar);

        strncpy(tmpvar, &end_time[9], 2);
        strcpy(&tmpvar[2], "\0");
        endmin = (short) atoi(tmpvar);

        strncpy(tmpvar, &end_time[11], 2);
        strcpy(&tmpvar[2], "\0");
        endsec = (short) atoi(tmpvar);

        strncpy(tmpvar, &end_time[13], 3);
        strcpy(&tmpvar[3], "\0");
        endmsec = atoi(tmpvar);

#ifdef DEBUG
        /* Print SeaWiFS Project required output to STDOUT */

        printf("startyear: %04d\n", startyear);
        printf("startday:  %03d\n", startday);
        printf("starthour: %02d\n", starthour);

        printf("endyear:   %04d\n", endyear);
        printf("endday:    %03d\n", endday);
        printf("endhour:   %02d\n", endhour);

        printf("%s/S%04d%03d%02d%03d%02d_%s.OZONE\n", directory,
                startyear, startday, starthour, endday, endhour, "TOVS");
#endif

        sprintf(outfile, "%s/S%04d%03d%02d%03d%02d_%s.OZONE", directory,
                startyear, startday, starthour, endday, endhour, "TOVS");
        sprintf(outfilename, "S%04d%03d%02d%03d%02d_%s.OZONE",
                startyear, startday, starthour, endday, endhour, "TOVS");

        startmsecday =
                (starthour * 60 * 60 * 1000) +
                (startmin * 60 * 1000) +
                (startsec * 1000) +
                startmsec;

        endmsecday =
                (endhour * 60 * 60 * 1000) +
                (endmin * 60 * 1000) +
                (endsec * 1000) +
                endmsec;
    } else if (type == 2) { /* EPTOMS */

        hour = 00;
        sprintf(outfile, "%s/N%04d%03d%02d_O3_TOMSOMI_24h.hdf",
                directory, year, julval, hour);

        sprintf(outfilename, "N%04d%03d%02d_O3_TOMSOMI_24h.hdf",
                year, julval, hour);

        /* ending times based on 1st file in EPTOMS case */
        sprintf(end_time, "%04d%03d235959999", year, julval);

        endyear = year;
        endday = julval;
        endmsecday = 86399999;

        /* set the starting times as 0Z on day of primary file */
        sprintf(start_time, "%04d%03d000000000", year, julval);
        startyear = year;
        startday = julval;
        startmsecday = 0;

        /* Determine metadata "Node Crossing Time" value from hh:mm [AP]M form */

        strcpy(&acntime[8], "\0");

        strncpy(tmpvar, &acntime[0], 2);
        strcpy(&tmpvar[2], "\0");
        hour = (short) atoi(tmpvar);

        strncpy(tmpvar, &acntime[6], 2); /* AM or PM? */
        if (!strcmp(tmpvar, "pm")) hour += 12;

        strncpy(tmpvar, &acntime[3], 2);
        strcpy(&tmpvar[2], "\0");
        min = (short) atoi(tmpvar);

        sprintf(node_time, "%04d%03d%02d%02d00000",
                year, julval, hour, min);

        /* Read previous day */

        rdgrid_(prevfile, &type, (int16*) prevdatarr, &julval, &year,
                start_time, end_time, acntime, &result, 255, 17, 17, 9);
        if (result == 1) pexit("rdgrid_");

        if (type != 2) pexit("2nd EPTOMS file not of EPTOMS type.  Exit.");

    } else if (type == 3) { /* ADTOMS */

        hour = 12; /* change to 12 from 0 */
        sprintf(outfile, "%s/S%04d%03d%02d_%s.OZONE",
                directory, year, julval, hour, "ADTOMS");

        sprintf(outfilename, "S%04d%03d%02d_%s.OZONE",
                year, julval, hour, "ADTOMS");

        sprintf(end_time, "%04d%03d235959999", year, julval);

        endyear = year;
        endday = julval;
        endmsecday = 86399999;

        /* Determine metadata "Node Crossing Time" value from hh:mm [AP]M form */

        strcpy(&acntime[8], "\0");

        strncpy(tmpvar, &acntime[0], 2);
        strcpy(&tmpvar[2], "\0");
        hour = (short) atoi(tmpvar);

        strncpy(tmpvar, &acntime[6], 2); /* AM or PM? */
        /*
         *  Note that we want the day-side crossing node.  For ADEOS, it is
         *  descending during daylight, so take ascending time and correct by 12 
         *  hours
         */
        if (!strcmp(tmpvar, "AM")) hour += 12;

        strncpy(tmpvar, &acntime[3], 2);
        strcpy(&tmpvar[2], "\0");
        min = (short) atoi(tmpvar);

        sprintf(node_time, "%04d%03d%02d%02d00000",
                year, julval, hour, min);

        sprintf(start_time, "%04d%03d000000000", year, julval);
        startyear = year;
        startday = julval;
        startmsecday = 0;
    } else pexit("ERROR: Neither TOMS, TOVS, EPTOMS or ADTOMS specified");

    if (!access(outfile, F_OK))
        pexit("....output file exists.  Won't overwrite.");

    /*
     * --- read annotation ASCII file for metadata placement in new HDF -----
     */

    if ((numannarr = count_annot(annotfile)) == 0)
        pexit("no annotations found");

    if ((xannot = (struct annotation *)
            malloc(sizeof (struct annotation) * numannarr))
            == NULL) pexit("malloc Annotation");

    xannot = fillenv_annot(annotfile);

    /*
     * Extract values from descr array element for program use
     */

    for (i = 0; i < numannarr; i++) {

        if (!strcmp(xannot[i].label, "Northernmost Latitude"))
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

    } /* for i */

    /*
     * Lat and Lon limits adjusted by 1/2 step size
     */

    nmostlat = nmostlat - (latstep / 0.5); /*   89.500 =  90 - (1.00 / 0.5) */
    smostlat = smostlat + (latstep / 0.5); /*  -89.500 = -90 + (1.00 / 0.5) */
    wmostlon = wmostlon + (lonstep / 0.5); /* -179.5 = 180 + (1.0 / 0.5) */
    emostlon = emostlon - (lonstep / 0.5); /*  179.5 = 180 - (1.0 / 0.5) */

    /*
     * Create HDF file
     */

    if ((result = startHDF(outfile, &sdfid, &fid, DFACC_CREATE)) != 0)
        pexit("Fatal error starting HDF file");

    /*
     * Allocate space for 2D data arrays
     */

    rank = 2;
    shape[0] = latsz; /* lat */
    shape[1] = lonsz; /* lon */
    array_size = shape[0] * shape[1];

    if ((int16_SDSdataO3 =
#ifdef USE_MALLOC
            (int16 *) malloc(sizeof (int16) * array_size)) == NULL)
        pexit("malloc int16_SDSdataO3");
#else
            (int16 *) calloc(sizeof (int16), array_size)) == NULL)
        pexit("calloc int16_SDSdataO3");
#endif

    if ((int8_SDSdataQC =
#ifdef USE_MALLOC
            (int8 *) malloc(sizeof (int8) * array_size)) == NULL)
        pexit("malloc int8_SDSdataQC");
#else
            (int8 *) calloc(sizeof (int8), array_size)) == NULL)
        pexit("calloc int8_SDSdataQC");
#endif
    for (i = 0; i < array_size; i++) int8_SDSdataQC[i] = 0; /* init */

    /*
     * Create SDS array for writing
     */

    l = 0;
    for (j = 0; j < shape[0]; j++) { /* lats */
        for (k = 0; k < shape[1]; k++) { /* lons */

            /* For EPTOMS data, overwrite missings (0's) with matching
                  position from previous day's data.  Mark QC array (as EARLIERVAL) */

            if ((type == 2) && (datarr[j][k] == 0)) {
                if (prevdatarr[j][k] != 0) {
                    datarr[j][k] = prevdatarr[j][k];
                    int8_SDSdataQC[l] = EARLIERVAL;
                } else { /* even 2nd file is missing for this point */

                    /* Scan lons in either longitudinal direction from current point */

                    for (kk = 1; kk <= 50; kk++) {
                        if ((k - kk >= 0) && (datarr[j][k - kk] > 0)) {
                            datarr[j][k] = datarr[j][k - kk];
                            int8_SDSdataQC[l] = NEARESTVAL;
                            break;
                        } else if ((k + kk < shape[1]) && (datarr[j][k + kk] > 0)) {
                            datarr[j][k] = datarr[j][k + kk];
                            int8_SDSdataQC[l] = NEARESTVAL;
                            break;
                        }
                    }
                }
                modified_pts++;
            }

            if (datarr[j][k] > 0) {
                int16_SDSdataO3[l] = (int16) datarr[j][k];
            } else int16_SDSdataO3[l] = 0;
            l++; /* index counter */
        } /* for k */
    } /* for j */

    /*
     *  WDR for EP-TOMS, apply an average to the data
     */
    if (type == 2) world_avg((int16 *) datarr, shape[0], shape[1], 5, 5,
            int16_SDSdataO3);

    /*
     * -------  Assign metadata values to local descriptive variables ----
     * -------   insert dates and other values to metadata array  ---------
     */

    for (i = 0; i < numannarr; i++) {

        if (!strcmp(xannot[i].label, "Product Name"))
            sprintf(xannot[i].descr, "%s", outfilename);

        else if (!strcmp(xannot[i].label, "Processing Time")) {

            /* Determine processing time */
            (void) time(&t);
            localtm = localtime(&t);
            pyear = localtm->tm_year;
            if (pyear < 90) pyear = pyear + 2000; /* 2 digit to 4 */
            else pyear = pyear + 1900;

            sprintf(xannot[i].descr, "%04d%03d%02d%02d%02d%03d",
                    pyear, localtm->tm_yday + 1, localtm->tm_hour, localtm->tm_min,
                    localtm->tm_sec, 0);
        }
        else if (!strcmp(xannot[i].label, "Input Files")) {
            if (type == 2)
                sprintf(xannot[i].descr, "%s,%s", infile, prevfile);
            else
                sprintf(xannot[i].descr, "%s", infile);
        }
        else if (!strcmp(xannot[i].label, "Processing Control"))
            if (type == 2)
                sprintf(xannot[i].descr, "%s %s %s %s %s %s",
                    argv[0], annotfile, infile, directory, itype, prevfile);
            else
                sprintf(xannot[i].descr, "%s %s %s %s %s",
                    argv[0], annotfile, infile, directory, itype);

        else if (!strcmp(xannot[i].label, "Start Time"))
            sprintf(xannot[i].descr, "%s", start_time);

        else if (!strcmp(xannot[i].label, "End Time"))
            sprintf(xannot[i].descr, "%s", end_time);

        else if (!strcmp(xannot[i].label, "Start Year"))
            sprintf(xannot[i].descr, "%04d", startyear);

        else if (!strcmp(xannot[i].label, "Start Day"))
            sprintf(xannot[i].descr, "%03d", startday);

        else if (!strcmp(xannot[i].label, "Start Millisec"))
            sprintf(xannot[i].descr, "%08d", startmsecday);

        else if (!strcmp(xannot[i].label, "End Year"))
            sprintf(xannot[i].descr, "%04d", endyear);

        else if (!strcmp(xannot[i].label, "End Day"))
            sprintf(xannot[i].descr, "%03d", endday);

        else if (!strcmp(xannot[i].label, "End Millisec"))
            sprintf(xannot[i].descr, "%08d", endmsecday);

        else if (!strcmp(xannot[i].label, "Node Crossing Time"))
            sprintf(xannot[i].descr, "%s", node_time);

        else if (!strcmp(xannot[i].label, "Points Modified"))
            sprintf(xannot[i].descr, "%08d", modified_pts);

    } /* for i */


    /*
     * Write attribute array to HDF file
     */

    if ((result = wrtattr(sdfid, xannot, numannarr)) != 0) pexit("wrtattr");
    free(xannot);


    /*
     * --- Write OZONE SDS, units attribute, and add to grid Vgroup ----------
     */

    datalabel = "ozone";
    strcpy(dataunit, "Total ozone");
    strcpy(dataattr, "Dobson units");
    datafmt = "int16";
    datatype = DFNT_INT16;

    gridid = setupGrid(fid, VGROUPNAME);

    if ((SDSinFile(datalabel, dataunit, dataattr, datafmt,
            datatype, sdfid, rank, shape, int16_SDSdataO3, gridid)) != 0)
        pexit("SDSinFile ozone wrtsds");

    free(int16_SDSdataO3);

    /*
     * --- Write OZONE QC flag SDS, units attribute, and add to grid Vgroup ---
     */

    datalabel = "ozone_QC";
    strcpy(dataunit, "Total ozone Q/C flag");
    strcpy(dataattr, "");
    datafmt = "int8";
    datatype = DFNT_INT8;

    if ((SDSinFile(datalabel, dataunit, dataattr, datafmt,
            datatype, sdfid, rank, shape, int8_SDSdataQC, gridid)) != 0)
        pexit("SDSinFile ozone QC wrtsds");

    free(int8_SDSdataQC);

    deattachHDFgrid(gridid); /* deattach HDF grid */

    /*
     * close HDF structures
     */

    result = 0;
    if ((result = closeHDFstructs(sdfid, fid)) != 0) pexit("closeHDFstructs");


    /*
     * Print info to STDOUT for SeaWiFS Project
     */

    printf("%s+%s+%s\n", outfile, start_time, end_time);

    exit(0);

} /* main */

/*****************************************************************
 * File:        check_usage
 *
 * Purpose:     check command line arguments
 *
 * Description: check command line arguements and report proper usage
 *              on argument count error or no argumnts.
 *
 * Input parms:
 * char *annotfile - input file of DAAC compliant annots ("fillenv.tomsnrt")
 * char *infile    - input ozone file to process ("nimoz100.y93")
 * char *directory - output file directory (".")
 *
 * Output parms:none
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
 * Author:      Brian D. Schieber, GSC, 9/93
 *
 * Modification history:
 *****************************************************************/


int8 check_usage(int argc, char *argv[]) {
    if (argc < 5) {
        printf("\n\nUsage:\n");
        printf("\t%s <metafile> <file> <directory> <itype> [<prevfile>]\n", argv[0]);
        printf("\nWhere:\n");
        printf("\tmetafile:   HDF metadata file (in ASCII)\n");
        printf("\tfile:       file to process\n");
        printf("\tdirectory:  output directory (no ending slash)\n");
        printf("\titype:      input data type 'TOMS', 'TOVS', 'EPTOMS' or 'ADTOMS'\n");
        printf("\tprevfile:   [optional] previous day if type 'EPTOMS'\n");
        printf("\n\n");
        printf("\tExample:\n");
        printf("\t\t o3nrt $SDSDEMO/fillenv.eptomsnrt $SDSDEMO/ga960903.ept\n");
        printf("\t\t ./ 'EPTOMS' $SDSDEMO/ga960902.ept\n\n");
        return -1;
    }
    return 0;
}
