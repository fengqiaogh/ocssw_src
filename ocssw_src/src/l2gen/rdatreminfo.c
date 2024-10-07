/* ========================================================== */
/* rdatreminfo() - reads parameters from ATREM sensor files   */
/*                                                            */
/* Inputs:                                                    */
/*     sensorID - sensor ID as defined in l12_parms.h         */
/*     pname    - parameter name from sensor table            */
/*                                                            */
/* Outputs:                                                   */
/*     pval     - pointer to scalar or array containing param */
/*                                                            */
/* Written By:                                                */
/*     R. Healy, March 2015                                   */
/*     Based on rdsensorinfo.c                                */
/*                                                            */
/* ========================================================== */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "sensorInfo.h"
#include "genutils.h"

// this is really defined in l2gen/l12_params.h
// but I need it here also.
#define NEWSENSINFO     32   /* use test sensor info file             */
void parseline(char *line, char *name, char *value);

/**
 * read in atrem initialization information.
 *
 * @param sensorID id of the sensor to lookup
 * @param evalmask set to 32 to get info out of the "eval" directory
 * @param pname name of the parameter to retrieve
 * @param pval pointer to the requested data
 * @return the number of bands, -1 if error
 */
int32_t rdatreminfo(int32_t sensorID, int32_t evalmask, const char *pname, void **pval) {
    FILE *fp;
    char *filedir;
    char filename[FILENAME_MAX];
    char line[80];
    char name[80];
    char value[80];
    char param[80];
    char *p;
    int i;
    int status;

    static int32_t sensorID_s = -999;
    static int32_t nbands = 0;

    static float *fwhm, *fwave;
    static int32_t h2o, co2, co, ch4, no2, n2o, o3, o2;
    static float vrto3, sno2;
    static float window1, window2, window3, window4, wp94c, w1p14c;
    static int32_t nb1, nb2, nb3, nb4, nbp94, nb1p14, full_calc, dogeom;

    want_verbose = 0;
    if (sensorID != sensorID_s) {

        if (want_verbose)
            printf("\nLoading atrem (atmospheric removal - water vapor,etc.) information for %s\n",
                sensorId2SensorName(sensorID));

        /*								*/
        /* Locate atmocor datafile using environment variable           */
        /*								*/
        if ((filedir = getenv("OCDATAROOT")) == NULL) {
            printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
            return (-1);
        }

        strcpy(filename, filedir);
        if ((evalmask & NEWSENSINFO) != 0)
            strcat(filename, "/eval/");
        else
            strcat(filename, "/");
        strcat(filename, sensorId2SensorDir(sensorID));

        strcat(filename, "/msl12_atrem_info.dat");

        if (want_verbose)
            printf("Opening atrem information file %s\n Looking for %s\n", filename, pname);

        if ((fp = fopen(filename, "r")) == NULL) {
            fprintf(stderr, "-E- %s line %d: unable to open %s for reading\n",
                    __FILE__, __LINE__, filename);
            return (-1);
        }

        // read file to determine number of bands

        while (fgets(line, 80, fp)) {

            memset(name, '\0', sizeof (name));
            memset(value, '\0', sizeof (value));

            // skip comment lines, empty lines, and lines without a name = value pair.

            if (line[0] == '#' || line[0] == '\n')
                continue;
            if (!(p = strchr(line, '=')))
                continue;

            // parse parameter name and value, looking for Nbands

            parseline(line, name, value);

            if (strcmp(name, "Nbands") == 0) {
                nbands = (int32_t) atoi(value);
                break;
            }
        }

        if (nbands <= 0) {
            printf("-E- %s: error find Nbands in %s.\n", __FILE__, filename);
            return (-1);
        }

        // allocate space for static data

        status = 0;
        if ((fwhm = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((fwave = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if (status == 1) {
            printf(
                    "-E- %s: error allocating space for %d bands in sensor info.\n",
                    __FILE__, nbands);
            return (-1);
        }


        // Loop through each line of file looking for identifiers 

        while (fgets(line, 80, fp)) {

            memset(name, '\0', sizeof (name));
            memset(value, '\0', sizeof (value));

            // skip comment lines, empty lines, and lines without a name = value pair.

            if (line[0] == '#' || line[0] == '\n')
                continue;
            if (!(p = strchr(line, '=')))
                continue;

            // parse parameter name and value

            parseline(line, name, value);

            // Copy value to appropriate variable

            if (strcmp(name, "Nbands") == 0)
                nbands = (int32_t) atoi(value);
            else if (strcmp(name, "h2o") == 0)
                h2o = (int32_t) atoi(value);
            else if (strcmp(name, "o2") == 0)
                o2 = (int32_t) atoi(value);
            else if (strcmp(name, "o3") == 0)
                o3 = (int32_t) atoi(value);
            else if (strcmp(name, "ch4") == 0)
                ch4 = (int32_t) atoi(value);
            else if (strcmp(name, "co") == 0)
                co = (int32_t) atoi(value);
            else if (strcmp(name, "n2o") == 0)
                n2o = (int32_t) atoi(value);
            else if (strcmp(name, "no2") == 0)
                no2 = (int32_t) atoi(value);
            else if (strcmp(name, "co2") == 0)
                co2 = (int32_t) atoi(value);
            else if (strcmp(name, "vrto3") == 0)
                vrto3 = (float) atof(value);
            else if (strcmp(name, "sno2") == 0)
                sno2 = (float) atof(value);
            else if (strcmp(name, "window1") == 0)
                window1 = (float) atof(value);
            else if (strcmp(name, "window2") == 0)
                window2 = (float) atof(value);
            else if (strcmp(name, "window3") == 0)
                window3 = (float) atof(value);
            else if (strcmp(name, "window4") == 0)
                window4 = (float) atof(value);
            else if (strcmp(name, "w1p14c") == 0)
                w1p14c = (float) atof(value);
            else if (strcmp(name, "wp94c") == 0)
                wp94c = (float) atof(value);
            else if (strcmp(name, "nb1") == 0)
                nb1 = (int32_t) atoi(value);
            else if (strcmp(name, "nb2") == 0)
                nb2 = (int32_t) atoi(value);
            else if (strcmp(name, "nb3") == 0)
                nb3 = (int32_t) atoi(value);
            else if (strcmp(name, "nb4") == 0)
                nb4 = (int32_t) atoi(value);
            else if (strcmp(name, "nbp94") == 0)
                nbp94 = (int32_t) atoi(value);
            else if (strcmp(name, "nb1p14") == 0)
                nb1p14 = (int32_t) atoi(value);
            else if (strcmp(name, "full_calc") == 0)
                full_calc = (int32_t) atoi(value);
            else if (strcmp(name, "dogeom") == 0)
                dogeom = (int32_t) atoi(value);
            else {
                for (i = 0; i < nbands; i++) {
                    sprintf(param, "fwhm(%d)", i + 1);
                    if (strcmp(name, param) == 0) {
                        fwhm[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "Lambda(%d)", i + 1);
                    if (strcmp(name, param) == 0) {
                        fwave[i] = (float) atof(value);
                        break;
                    }

                }
            }

        }


        /*                                                            */
        /* Write-out what was read, for informational purposes        */
        /*                                                            */
        if (want_verbose) {

            printf("window1=%f\n", window1);
            printf("  %3s %8s\n", "Bnd",
                    "Fwhm");
            for (i = 0; i < nbands; i++)
                printf(
                    "  %3d %8.3f\n",
                    i, fwhm[i]);
            printf("\n");

            printf("\n");
        } // want_verbose
    }

    /*                                              */
    /* Return a pointer to the requested parameter  */
    /*                                              */
    if (pval != NULL) {
        if (strcmp(pname, "Nbands") == 0)
            *pval = (void *) &nbands;
        else if (strcmp(pname, "fwhm") == 0)
            *pval = (void *) fwhm;
        else if (strcmp(pname, "Lambda") == 0)
            *pval = (void *) fwave;
        else if (strcmp(pname, "o2") == 0)
            *pval = (void *) &o2;
        else if (strcmp(pname, "o3") == 0)
            *pval = (void *) &o3;
        else if (strcmp(pname, "n2o") == 0)
            *pval = (void *) &n2o;
        else if (strcmp(pname, "ch4") == 0)
            *pval = (void *) &ch4;
        else if (strcmp(pname, "co") == 0)
            *pval = (void *) &co;
        else if (strcmp(pname, "no2") == 0)
            *pval = (void *) &no2;
        else if (strcmp(pname, "h2o") == 0)
            *pval = (void *) &h2o;
        else if (strcmp(pname, "co2") == 0)
            *pval = (void *) &co2;
        else if (strcmp(pname, "vrto3") == 0)
            *pval = (void *) &vrto3;
        else if (strcmp(pname, "sno2") == 0)
            *pval = (void *) &sno2;
        else if (strcmp(pname, "window1") == 0)
            *pval = (void *) &window1;
        else if (strcmp(pname, "window2") == 0)
            *pval = (void *) &window2;
        else if (strcmp(pname, "window3") == 0)
            *pval = (void *) &window3;
        else if (strcmp(pname, "window4") == 0)
            *pval = (void *) &window4;
        else if (strcmp(pname, "wp94c") == 0)
            *pval = (void *) &wp94c;
        else if (strcmp(pname, "w1p14c") == 0)
            *pval = (void *) &w1p14c;
        else if (strcmp(pname, "nb1") == 0)
            *pval = (void *) &nb1;
        else if (strcmp(pname, "nb2") == 0)
            *pval = (void *) &nb2;
        else if (strcmp(pname, "nb3") == 0)
            *pval = (void *) &nb3;
        else if (strcmp(pname, "nb4") == 0)
            *pval = (void *) &nb4;
        else if (strcmp(pname, "nbp94") == 0)
            *pval = (void *) &nbp94;
        else if (strcmp(pname, "nb1p14") == 0)
            *pval = (void *) &nb1p14;
        else if (strcmp(pname, "full_calc") == 0)
            *pval = (void *) &full_calc;
        else if (strcmp(pname, "dogeom") == 0)
            *pval = (void *) &dogeom;

        else
            return (-1);
    }

    sensorID_s = sensorID;

    if (pname != NULL)
        return (nbands);

    return (0);
}

