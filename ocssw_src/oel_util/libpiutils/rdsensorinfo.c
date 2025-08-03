/* ========================================================== */
/* rdsensorinfo() - reads parameters from MSL12 sensor files  */
/*                                                            */
/* Inputs:                                                    */
/*     sensorID - sensor ID as defined in l12_parms.h         */
/*     pname    - parameter name from sensor table            */
/*                                                            */
/* Outputs:                                                   */
/*     pval     - pointer to scalar or array containing param */
/*                                                            */
/* Written By:                                                */
/*     B. Franz, October 2002                                 */
/*                                                            */
/* Converted to dynamic allocation, February 2013, BAF        */
/*                                                            */
/* ========================================================== */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <sensorInfo.h>
#include <genutils.h>

#define SENSOR_INFO_FILENAME "msl12_sensor_info.dat"

// this is really defined in l2gen/l12_params.h
// but I need it here also.
#define NEWSENSINFO     32   /* use test sensor info file             */

void parseline(char *line, char *name, char *value) {
    char *p;
    char *p1;
    char *p2;

    // Parse parameter name string

    p = strchr(line, '=');

    p1 = line;
    while (isspace(*p1))
        p1++;
    p2 = p - 1;
    while (isspace(*p2))
        p2--;
    strncpy(name, p1, p2 - p1 + 1);

    // Parse parameter value string

    p1 = p + 1;
    while (isspace(*p1))
        p1++;
    p2 = p1;
    while (!isspace(*p2))
        p2++;
    strncpy(value, p1, p2 - p1);
}

/**
 * lookup information about a sensor.
 *
 * @param sensorID id of the sensor to lookup
 * @param evalmask set to 32 to get info out of the "eval" directory
 * @param pname name of the parameter to retrieve
 * @param pval pointer to the requested data
 * @return the number of bands, -1 if error
 */
int32_t rdsensorinfo(int32_t sensorID, int32_t evalmask, const char *pname,
        void **pval) {
    FILE *fp;
    char *filedir;
    char filename[FILENAME_MAX];
    char prefix[FILENAME_MAX];
    const char* sensorDir;
    const char* subsensorDir;
    char line[80];
    char name[80];
    char value[80];
    char param[80];
    char wparam[80];
    char *p;
    int i;
    int status;

    static int32_t sensorID_s = -999;
    static int32_t nbands = 0;
    static int32_t nbandsir = 0;
    static int32_t nbandsvis = 0;

    static int32_t *bindx;
    static int32_t *iwave;
    static float *fwave;
    static float *Fo;
    static float *Tau_r;
    static float *k_oz;
    static float *t_co2;
    static float *k_no2;
    static float *a_h2o;
    static float *b_h2o;
    static float *c_h2o;
    static float *d_h2o;
    static float *e_h2o;
    static float *f_h2o;
    static float *g_h2o;
    static float *awhite;
    static float *aw;
    static float *bbw;
    static float *wed;
    static float *waph;

    static float *ooblw01;
    static float *ooblw02;
    static float *ooblw03;

    static float *oobwv01;
    static float *oobwv02;
    static float *oobwv03;
    static float *oobwv04;
    static float *oobwv05;
    static float *oobwv06;
    static float *oobwv07;
    static float *oobwv08;
    static float *oobwv09;
    static float *oobwv10;
    static float *oobwv11;
    static float *oobwv12;

    if (sensorID != sensorID_s) {

        if (want_verbose)
            printf("\nLoading characteristics for %s\n", sensorId2SensorName(sensorID));

        /*								*/
        /* Locate atmocor datafile using environment variable           */
        /*								*/
        if ((filedir = getenv("OCDATAROOT")) == NULL) {
            printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
            return (-1);
        }

        strcpy(prefix, filedir);
        if ((evalmask & NEWSENSINFO) != 0)
            strcat(prefix, "/eval/");
        else
            strcat(prefix, "/");

        sensorDir = sensorId2SensorDir(sensorID);
        if (sensorDir)
            strcat(prefix, sensorDir);
        else {
            printf("\n -E- Can't find sensor directory for ID %d\n", sensorID);
            exit(EXIT_FAILURE);
        }
        strcat(prefix, "/");

        fp = NULL;
        
        // first try in subsensor directory
        subsensorDir = subsensorId2SubsensorDir(sensorId2SubsensorId(sensorID));
        if(subsensorDir) {
            strcpy(filename, prefix);
            strcat(filename, subsensorDir);
            strcat(filename, "/");
            strcat(filename, SENSOR_INFO_FILENAME);
            fp = fopen(filename, "r");
        }

        // if needed look in the sensor directory
        if(fp == NULL) {
            strcpy(filename, prefix);
            strcat(filename, SENSOR_INFO_FILENAME);
            if ((fp = fopen(filename, "r")) == NULL) {
                fprintf(stderr, "-E- %s line %d: unable to open %s for reading\n",
                        __FILE__, __LINE__, filename);
                return (-1);
            }
        }
        
        if (want_verbose)
            printf("Opening sensor information file %s\n", filename);

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
        if ((bindx = (int32_t*) calloc(nbands, sizeof (int32_t))) == NULL)
            status = 1;
        if ((iwave = (int32_t*) calloc(nbands, sizeof (int32_t))) == NULL)
            status = 1;
        if ((fwave = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((Fo = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((Tau_r = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((k_oz = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((awhite = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((aw = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((bbw = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((wed = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((waph = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((a_h2o = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((b_h2o = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((c_h2o = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((d_h2o = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((e_h2o = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((f_h2o = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((g_h2o = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((t_co2 = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((k_no2 = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((ooblw01 = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((ooblw02 = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((ooblw03 = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((oobwv01 = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((oobwv02 = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((oobwv03 = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((oobwv04 = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((oobwv05 = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((oobwv06 = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((oobwv07 = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((oobwv08 = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((oobwv09 = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((oobwv10 = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((oobwv11 = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if ((oobwv12 = (float*) calloc(nbands, sizeof (float))) == NULL)
            status = 1;
        if (status == 1) {
            printf(
                    "-E- %s: error allocating space for %d bands in sensor info.\n",
                    __FILE__, nbands);
            return (-1);
        }

        // initialize non-zero output data

        for (i = 0; i < nbands; i++) {
            bindx[i] = i; // don't delete yet
            t_co2[i] = 1.0;
            k_no2[i] = 1.0;
            a_h2o[i] = 1.0;
            // initialize water vap. correction coefficients if not specified in file
            oobwv01[i] = 1;
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
            else {
                for (i = 0; i < nbands; i++) {
                    sprintf(param, "Lambda(%d)", i + 1);
                    if (strcmp(name, param) == 0) {
                        fwave[i] = (float) atof(value);
                        iwave[i] = (int32_t) roundf(fwave[i]);
                        break;
                    }
                }
                for (i = 0; i < nbands; i++) {
                    sprintf(param, "F0(%d)", i + 1);
                    sprintf(wparam, "F0(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        Fo[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "Tau_r(%d)", i + 1);
                    sprintf(wparam, "Tau_r(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        Tau_r[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "k_oz(%d)", i + 1);
                    sprintf(wparam, "k_oz(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        k_oz[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "t_co2(%d)", i + 1);
                    sprintf(wparam, "t_co2(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        t_co2[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "k_no2(%d)", i + 1);
                    sprintf(wparam, "k_no2(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        k_no2[i] = (float) atof(value);
                        break;
                    }

                    sprintf(param, "a_h2o(%d)", i + 1);
                    sprintf(wparam, "a_h2o(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        a_h2o[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "b_h2o(%d)", i + 1);
                    sprintf(wparam, "b_h2o(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        b_h2o[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "c_h2o(%d)", i + 1);
                    sprintf(wparam, "c_h2o(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        c_h2o[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "d_h2o(%d)", i + 1);
                    sprintf(wparam, "d_h2o(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        d_h2o[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "e_h2o(%d)", i + 1);
                    sprintf(wparam, "e_h2o(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        e_h2o[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "f_h2o(%d)", i + 1);
                    sprintf(wparam, "f_h2o(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        f_h2o[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "g_h2o(%d)", i + 1);
                    sprintf(wparam, "g_h2o(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        g_h2o[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "awhite(%d)", i + 1);
                    sprintf(wparam, "awhite(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        awhite[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "aw(%d)", i + 1);
                    sprintf(wparam, "aw(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        aw[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "bbw(%d)", i + 1);
                    sprintf(wparam, "bbw(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        bbw[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "wed(%d)", i + 1);
                    sprintf(wparam, "wed(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        wed[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "waph(%d)", i + 1);
                    sprintf(wparam, "waph(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        waph[i] = (float) atof(value);
                        break;
                    }

                    sprintf(param, "ooblw01(%d)", i + 1);
                    sprintf(wparam, "ooblw01(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        ooblw01[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "ooblw02(%d)", i + 1);
                    sprintf(wparam, "ooblw02(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        ooblw02[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "ooblw03(%d)", i + 1);
                    sprintf(wparam, "ooblw03(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        ooblw03[i] = (float) atof(value);
                        break;
                    }

                    sprintf(param, "oobwv01(%d)", i + 1);
                    sprintf(wparam, "oobwv01(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        oobwv01[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "oobwv02(%d)", i + 1);
                    sprintf(wparam, "oobwv02(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        oobwv02[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "oobwv03(%d)", i + 1);
                    sprintf(wparam, "oobwv03(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        oobwv03[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "oobwv04(%d)", i + 1);
                    sprintf(wparam, "oobwv04(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        oobwv04[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "oobwv05(%d)", i + 1);
                    sprintf(wparam, "oobwv05(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        oobwv05[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "oobwv06(%d)", i + 1);
                    sprintf(wparam, "oobwv06(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        oobwv06[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "oobwv07(%d)", i + 1);
                    sprintf(wparam, "oobwv07(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        oobwv07[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "oobwv08(%d)", i + 1);
                    sprintf(wparam, "oobwv08(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        oobwv08[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "oobwv09(%d)", i + 1);
                    sprintf(wparam, "oobwv09(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        oobwv09[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "oobwv10(%d)", i + 1);
                    sprintf(wparam, "oobwv10(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        oobwv10[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "oobwv11(%d)", i + 1);
                    sprintf(wparam, "oobwv11(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        oobwv11[i] = (float) atof(value);
                        break;
                    }
                    sprintf(param, "oobwv12(%d)", i + 1);
                    sprintf(wparam, "oobwv12(%d)", iwave[i]);
                    if (strcmp(name, param) == 0 || strcmp(name, wparam) == 0) {
                        oobwv12[i] = (float) atof(value);
                        break;
                    }

                }
            }

        }

        /*                                                            */
        /* Separate bands into reflected and emmissive and visible    */
        /*                                                            */
        for (i = 0; i < nbands; i++)
            if (iwave[i] > MINWAVE_IR)
                break;

        nbandsir = nbands - i;
        nbands = i;

        for (i = 0; i < nbands; i++)
            if (iwave[i] > MAXWAVE_VIS)
                break;

        nbandsvis = i;

        /*                                                            */
        /* Write-out what was read, for informational purposes        */
        /*                                                            */
        if (want_verbose) {

            printf("  %3s %5s %8s %8s %8s %8s %8s %8s %8s %8s\n", "Bnd", "Lam",
                    "Fo", "Tau_r", "k_oz", "k_no2", "t_co2", "awhite", "aw",
                    "bbw");
            for (i = 0; i < nbands; i++)
                printf(
                    "  %3d %8.3f %8.3f %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e %8.3e\n",
                    i, fwave[i], Fo[i], Tau_r[i], k_oz[i], k_no2[i],
                    t_co2[i], awhite[i], aw[i], bbw[i]);
            printf("\n");

            if (nbandsir > 0) {
                printf("  %3s %5s\n", "Bnd", "Lam");
                for (i = 0; i < nbandsir; i++)
                    printf("  %3d %8.3f\n", i + nbands, fwave[i + nbands]);
            }
            printf("\n");
        } // want_verbose
    }

    /*                                              */
    /* Return a pointer to the requested parameter  */
    /*                                              */
    if (pval != NULL) {
        if (strcmp(pname, "Nbands") == 0)
            *pval = (void *) &nbands;
        else if (strcmp(pname, "NbandsIR") == 0)
            *pval = (void *) &nbandsir;
        else if (strcmp(pname, "NbandsVIS") == 0)
            *pval = (void *) &nbandsvis;
        else if (strcmp(pname, "Bindx") == 0)
            *pval = (void *) bindx;
        else if (strcmp(pname, "fwave") == 0)
            *pval = (void *) fwave;
        else if (strcmp(pname, "iwave") == 0)
            *pval = (void *) iwave;
        else if (strcmp(pname, "Lambda") == 0)
            *pval = (void *) iwave;
        else if (strcmp(pname, "Fobar") == 0)
            *pval = (void *) Fo;
        else if (strcmp(pname, "Tau_r") == 0)
            *pval = (void *) Tau_r;
        else if (strcmp(pname, "k_oz") == 0)
            *pval = (void *) k_oz;
        else if (strcmp(pname, "t_co2") == 0)
            *pval = (void *) t_co2;
        else if (strcmp(pname, "k_no2") == 0)
            *pval = (void *) k_no2;

        else if (strcmp(pname, "a_h2o") == 0)
            *pval = (void *) a_h2o;
        else if (strcmp(pname, "b_h2o") == 0)
            *pval = (void *) b_h2o;
        else if (strcmp(pname, "c_h2o") == 0)
            *pval = (void *) c_h2o;
        else if (strcmp(pname, "d_h2o") == 0)
            *pval = (void *) d_h2o;
        else if (strcmp(pname, "e_h2o") == 0)
            *pval = (void *) e_h2o;
        else if (strcmp(pname, "f_h2o") == 0)
            *pval = (void *) f_h2o;
        else if (strcmp(pname, "g_h2o") == 0)
            *pval = (void *) g_h2o;

        else if (strcmp(pname, "awhite") == 0)
            *pval = (void *) awhite;

        else if (strcmp(pname, "aw") == 0)
            *pval = (void *) aw;
        else if (strcmp(pname, "bbw") == 0)
            *pval = (void *) bbw;

        else if (strcmp(pname, "wed") == 0)
            *pval = (void *) wed;
        else if (strcmp(pname, "waph") == 0)
            *pval = (void *) waph;

        else if (strcmp(pname, "ooblw01") == 0)
            *pval = (void *) ooblw01;
        else if (strcmp(pname, "ooblw02") == 0)
            *pval = (void *) ooblw02;
        else if (strcmp(pname, "ooblw03") == 0)
            *pval = (void *) ooblw03;

        else if (strcmp(pname, "oobwv01") == 0)
            *pval = (void *) oobwv01;
        else if (strcmp(pname, "oobwv02") == 0)
            *pval = (void *) oobwv02;
        else if (strcmp(pname, "oobwv03") == 0)
            *pval = (void *) oobwv03;
        else if (strcmp(pname, "oobwv04") == 0)
            *pval = (void *) oobwv04;
        else if (strcmp(pname, "oobwv05") == 0)
            *pval = (void *) oobwv05;
        else if (strcmp(pname, "oobwv06") == 0)
            *pval = (void *) oobwv06;
        else if (strcmp(pname, "oobwv07") == 0)
            *pval = (void *) oobwv07;
        else if (strcmp(pname, "oobwv08") == 0)
            *pval = (void *) oobwv08;
        else if (strcmp(pname, "oobwv09") == 0)
            *pval = (void *) oobwv09;
        else if (strcmp(pname, "oobwv10") == 0)
            *pval = (void *) oobwv10;
        else if (strcmp(pname, "oobwv11") == 0)
            *pval = (void *) oobwv11;
        else if (strcmp(pname, "oobwv12") == 0)
            *pval = (void *) oobwv12;

        else
            return (-1);
    }

    sensorID_s = sensorID;

    if (pname != NULL) {
        if (strcmp(pname, "NbandsIR") == 0)
            return (nbandsir);
        else if (strcmp(pname, "NbandsVIS") == 0)
            return (nbandsvis);
        else
            return (nbands);
    } else
        return (nbands);
}

