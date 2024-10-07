/* ====================================================== */
/* Module l1_mos_hdf.c                                    */
/*                                                        */
/* Functions to open and read a MOS HDF l1b file.         */
/*                                                        */
/* Written By:                                            */
/*    JT                                                  */
/*    NASA/SIMBIOS Project                                */
/*    11/98                                               */
/*                                                        */
/* ====================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h> 
#include "l1.h"
#include "filehdr_struc.h"
#include <hdf4utils.h>
#include "l1_mos.h"

#include <hdf.h>
#include <mfhdf.h>

#define NP 384       /* # pixels for MOS files               */
#define NS 384       /* # scans for MOS files                */
#define NB 13        /* Number of channels                   */
#define MAXLEN 200   /* length of file description array     */
#define PSIZE 2      /* Array size for the tabulated values  */
#define SPACE 32     /* ASCII Blank                          */
#define NUM 63       /* # elements in id[] and desc[] arrays */
#define NBANDS_MOS 8

static float gscan[NBANDS_MOS][4];



/* Function prototypes */
/* ------------------- */
int rd_gscan(void);
int ReadSDS(filehandle *l1file);
int32_t GetTime(char arr[MAXLEN]);
int32_t TimeToSec(int32_t yr, int32_t mo, int32_t dm, int32_t hr, int32_t mn, int32_t sec);
int GetYearDayMsec(double dtime[NS]);
int LeapChk(int32_t yr);
int RemoveWhitespace(char arr[MAXLEN]);
int bilinear(float p[PSIZE][PSIZE], float x[NS][NP], float y[NS][NP], float z[NS][NP]);

/* Static globals */
/* -------------- */
static char *desc_buffer = NULL; /* buffer for the file description */

typedef uint16_t array_data_t[NS];
static array_data_t* array_data;
//static uint16_t array_data[NP*14][NS]; /* full image data array */

typedef float data_2d_t[NP];
static data_2d_t *lon;
static data_2d_t *lat;
static data_2d_t *solz;
static data_2d_t *sola;
static data_2d_t *senz;
static data_2d_t *sena;
//static float lon     [NS][NP];
//static float lat     [NS][NP];
//static float solz    [NS][NP];
//static float sola    [NS][NP];
//static float senz    [NS][NP];
//static float sena    [NS][NP];

typedef float image_t[NP][NBANDS_MOS];
static image_t* image;
//static float image   [NS][NP][NBANDS_MOS];

static int32_t* yeararr;
static int32_t* dayarr;
static int32_t* msecarr;
//static int32_t  yeararr [NS];
//static int32_t  dayarr  [NS];
//static int32_t  msecarr [NS];

static int32_t begdm;
static int32_t enddm;
static int32_t begmo;
static int32_t endmo;
static int32_t begyr;
static int32_t endyr;

static int bands[NBANDS_MOS] = {0, 1, 2, 3, 4, 7, 8, 10};
static float acal[NB];


/* ------------------------------------------------------ */
/* openl1_read_mos() - opens a MOS HDF level-1 file   */
/*                    for reading, if not already opened. */
/*                    Stores the file data for lon, lat,  */
/*                    etc. in 384 x 384 arrays.           */
/*                                                        */

/* ------------------------------------------------------ */
int openl1_read_mos(filehandle *l1file) {
    data_2d_t *x;
    data_2d_t *y;
    //float x[NS][NP];
    //float y[NS][NP];

    float ul[4], ur[4], ll[4], lr[4];
    float p[PSIZE][PSIZE];
    int i, j, k, m, n, len;
    int status;
    char desc[MAXLEN], id[40], calstr[20];
    char *loc1, *loc2, *idend, *cal1, *cal2;
    int32_t hr, mn, sec;
    int32_t begMsec, endMsec, endTime, begTime;
    double eTime, bTime, dtime[NS];

    // allocate global data
    array_data = (array_data_t*) malloc((NP * 14) * sizeof (array_data_t));
    lon = (data_2d_t*) malloc(NS * sizeof (data_2d_t));
    lat = (data_2d_t*) malloc(NS * sizeof (data_2d_t));
    solz = (data_2d_t*) malloc(NS * sizeof (data_2d_t));
    sola = (data_2d_t*) malloc(NS * sizeof (data_2d_t));
    senz = (data_2d_t*) malloc(NS * sizeof (data_2d_t));
    sena = (data_2d_t*) malloc(NS * sizeof (data_2d_t));
    image = (image_t*) malloc(NS * sizeof (image_t));
    yeararr = (int32_t*) malloc(NS * sizeof (int32_t));
    dayarr = (int32_t*) malloc(NS * sizeof (int32_t));
    msecarr = (int32_t*) malloc(NS * sizeof (int32_t));

    // local data
    x = (data_2d_t*) malloc(NS * sizeof (data_2d_t));
    y = (data_2d_t*) malloc(NS * sizeof (data_2d_t));

    /* Get relative cal gains */
    if (rd_gscan() != 0)
        return (1);

    /* Set header info */
    l1file->npix = NP;
    l1file->nscan = NS;

    /* Get file description from the HDF file */
    if ((desc_buffer = GetFileDesc(l1file->name)) == NULL)
        return (1);

    /* Load the id[] and desc[] arrays */
    loc1 = desc_buffer;
    for (i = 0; i < NUM; i++) {
        loc2 = strchr(loc1, ';');
        if (loc2 == NULL) break;
        len = loc2 - loc1;
        desc[0] = '\0';
        strncat(desc, loc1, len);

        /* Remove whitespace characters */
        status = RemoveWhitespace(desc);

        /* Break into id[] and desc[] arrays */
        idend = strchr(desc, '=');
        len = idend - desc;
        id[0] = '\0';
        strncat(id, desc, len);
        for (j = 0; j < MAXLEN; j++) {
            desc[j] = desc[j + len];
        }

        /* Remove whitespace characters */
        status = RemoveWhitespace(desc);

        loc1 = loc2 + 1;

        /* Convert header information to usable format */
        for (k = 0; k < 1; k++) {
            if (strcmp(id, "LONG_UL") == 0) {
                ul[0] = atof(desc);
                break;
            }
            if (strcmp(id, "LATI_UL") == 0) {
                ul[1] = atof(desc);
                break;
            }
            if (strcmp(id, "SUNZ_UL") == 0) {
                ul[2] = atof(desc);
                break;
            }
            if (strcmp(id, "SUNA_UL") == 0) {
                ul[3] = atof(desc);
                break;
            }
            if (strcmp(id, "LONG_UR") == 0) {
                ur[0] = atof(desc);
                break;
            }
            if (strcmp(id, "LATI_UR") == 0) {
                ur[1] = atof(desc);
                break;
            }
            if (strcmp(id, "SUNZ_UR") == 0) {
                ur[2] = atof(desc);
                break;
            }
            if (strcmp(id, "SUNA_UR") == 0) {
                ur[3] = atof(desc);
                break;
            }
            if (strcmp(id, "LONG_LL") == 0) {
                ll[0] = atof(desc);
                break;
            }
            if (strcmp(id, "LATI_LL") == 0) {
                ll[1] = atof(desc);
                break;
            }
            if (strcmp(id, "SUNZ_LL") == 0) {
                ll[2] = atof(desc);
                break;
            }
            if (strcmp(id, "SUNA_LL") == 0) {
                ll[3] = atof(desc);
                break;
            }
            if (strcmp(id, "LONG_LR") == 0) {
                lr[0] = atof(desc);
                break;
            }
            if (strcmp(id, "LATI_LR") == 0) {
                lr[1] = atof(desc);
                break;
            }
            if (strcmp(id, "SUNZ_LR") == 0) {
                lr[2] = atof(desc);
                break;
            }
            if (strcmp(id, "SUNA_LR") == 0) {
                lr[3] = atof(desc);
                break;
            }
            if (strcmp(id, "BCAL_VAL") == 0) {
                len = 0;
                for (m = 0; m < NB; m++) {
                    for (n = 0; n < MAXLEN; n++) {
                        desc[n] = desc[n + len];
                    }
                    status = RemoveWhitespace(desc);
                    cal1 = desc;
                    cal2 = strchr(cal1, SPACE);
                    len = cal2 - cal1;
                    calstr[0] = '\0';
                    strncat(calstr, cal1, len);
                    acal[m] = atof(calstr);
                    cal1 = cal2 + 1;
                }
                break;
            }
            if (strcmp(id, "OP_BEG_DATE") == 0) {
                begdm = GetTime(desc);
                begmo = GetTime(desc);
                begyr = GetTime(desc);
                if (begyr > 50) {
                    begyr = begyr + 1900;
                } else {
                    begyr = begyr + 2000;
                }
                hr = GetTime(desc);
                mn = GetTime(desc);
                sec = GetTime(desc);
                begTime = TimeToSec(begyr, begmo, begdm, hr, mn, sec);
                break;
            }
            if (strcmp(id, "OP_END_DATE") == 0) {
                enddm = GetTime(desc);
                endmo = GetTime(desc);
                endyr = GetTime(desc);
                if (endyr > 50) {
                    endyr = endyr + 1900;
                } else {
                    endyr = endyr + 2000;
                }
                hr = GetTime(desc);
                mn = GetTime(desc);
                sec = GetTime(desc);
                endTime = TimeToSec(endyr, endmo, enddm, hr, mn, sec);
                break;
            }
            if (strcmp(id, "OP_BEG_MS") == 0) {
                begMsec = atol(desc);
                break;
            }
            if (strcmp(id, "OP_END_MS") == 0) {
                endMsec = atol(desc);
                break;
            }
            //         if (strcmp(id,"CLOUD_PERC") == 0) {
            //             l1file->percent_cloud = atoi(desc);
            //             break;
            //         }
            //         if (strcmp(id,"LAND_PERC") == 0) {
            //             l1file->percent_land = atoi(desc);
            //             break;
            //         }
            //         if (strcmp(id,"WATER_PERC") == 0) {
            //             l1file->percent_water = atoi(desc);
            //             break;
            //         }
        } /* k loop */
    } /* i loop */

    bTime = (double) begTime + begMsec / 1000.0;
    eTime = (double) endTime + endMsec / 1000.0;

    /* set yeararr, dayarr & msecarr arrays */
    for (i = 0; i < NS; i++) {
        dtime[i] = bTime + (eTime - bTime) / (NS - 1) * i;
    }
    status = GetYearDayMsec(dtime);

    /* set x and y matrices, using i for rows, j for columns */
    for (i = 0; i < NS; i++) {
        for (j = 0; j < NP; j++) {
            x[i][j] = (1.0 * j) / (NP - 1);
            y[i][j] = (1.0 * i) / (NS - 1);
        }
    }

    p[0][0] = ul[0];
    p[0][1] = ur[0];
    p[1][0] = ll[0];
    p[1][1] = lr[0];
    status = bilinear(p, x, y, lon);
    p[0][0] = ul[1];
    p[0][1] = ur[1];
    p[1][0] = ll[1];
    p[1][1] = lr[1];
    status = bilinear(p, x, y, lat);
    p[0][0] = ul[2];
    p[0][1] = ur[2];
    p[1][0] = ll[2];
    p[1][1] = lr[2];
    status = bilinear(p, x, y, solz);
    p[0][0] = ul[3];
    p[0][1] = ur[3];
    p[1][0] = ll[3];
    p[1][1] = lr[3];
    status = bilinear(p, x, y, sola);

    /* keep the longitude in the -180 to 180 range      */
    for (i = 0; i < NS; i++) {
        for (j = 0; j < NP; j++) {
            if (lon[i][j] > 180) {
                lon[i][j] = lon[i][j] - 360;
            } else {
                if (lon[i][j] < -180)
                    lon[i][j] = lon[i][j] + 360;
            }
        }
    }


    /* No sensor geometry is available.  Assume zenith is zero */
    /* at center, 7-deg at edge.  Azimuth is 90-deg relative to */
    /* solar azimuth */
    for (i = 0; i < NS; i++) {
        for (j = 0; j < NP; j++) {
            senz[i][j] = fabs(-7.0 + x[i][j] * 14.0);
            sena[i][j] = sola[i][j] + 270;
            if (sena[i][j] > 360.0) sena[i][j] -= 360.0;
        }
    }

    printf("\nMOS Scan-Dependent Calibration Coefficients (Band C0 C1 C2 C3)\n");
    for (i = 0; i < NBANDS_MOS; i++) {
        printf("  %d: %e, %e, %e, %e\n", i + 1,
                gscan[i][0], gscan[i][1], gscan[i][2], gscan[i][3]);
    }

    free(x);
    free(y);

    return (status);
}


/* ------------------------------------------------------ */
/* readl1_mos() - reads a MOS HDF level-1 record.     */
/*                                                        */

/* ------------------------------------------------------ */
int readl1_mos(filehandle *l1file, int32_t recnum, l1str *l1rec) {
    int i, j, k;
    static char current_file[FILENAME_MAX] = "";
    int ibnd;
    float gain, x;
    float Lt7;

    int32_t ip, ib, iw;
    int32_t nwave = l1rec->l1file->nbands;
    int32_t *bindx = l1rec->l1file->bindx;

    /* check current recnum */
    if (recnum >= NS) {
        fprintf(stderr, "-W- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "readl1_mos() called with record number (%d) ", recnum);
        fprintf(stderr, "that is inappropriate to the number of scanlines (%d) ", NS);
        fprintf(stderr, "in the current HDF file, %s . Call ignored.\n", l1file->name);
        return (1);
    }


    /* On first access, load entire image array into memory */
    if (strcmp(current_file, l1file->name) != 0) {

        printf("Reading image data block.\n");

        strcpy(current_file, l1file->name);

        /* read in the image data */
        if (ReadSDS(l1file) != 0)
            return (1);

        /* create the image array */
        for (i = 0; i < NBANDS_MOS; i++) {
            for (j = 0; j < NS; j++) {
                ibnd = j * 14 + bands[i];
                for (k = 0; k < NP; k++) {
                    x = k + 1;
                    gain = gscan[i][0] + x * (gscan[i][1] + x * (gscan[i][2] + x * gscan[i][3]));
                    image[j][k][i] = array_data[ibnd][k] * acal[bands[i]] * gain;
                }
            }
        }
    }
    l1rec->scantime = yds2unix((int16_t) yeararr[recnum], (int16_t) dayarr[recnum],
            (double) (msecarr[recnum] / 1000.0));

    /* copy the scan arrays */
    memcpy(l1rec->lon, &lon [recnum][0], sizeof (float)*NP);
    memcpy(l1rec->lat, &lat [recnum][0], sizeof (float)*NP);
    memcpy(l1rec->solz, &solz[recnum][0], sizeof (float)*NP);
    memcpy(l1rec->sola, &sola[recnum][0], sizeof (float)*NP);
    memcpy(l1rec->senz, &senz[recnum][0], sizeof (float)*NP);
    memcpy(l1rec->sena, &sena[recnum][0], sizeof (float)*NP);

    /* copy the image data for this line */
    /* The input files do not contain the view angles per band, so we 
       must replicate them here from the nominal view angles */
    for (ip = 0; ip < NP; ip++) {
        for (iw = 0; iw < nwave; iw++) {
            ib = bindx[iw];
            l1rec->Lt[ip * nwave + ib] = image[recnum][ip][iw];
        }
    }

    /* Correction for 750 channel non-linearity in MOS */
    for (j = 0; j < l1file->npix; j++) {
        Lt7 = l1rec->Lt[j * nwave + bindx[6]];
        if (Lt7 > 0.0)
            l1rec->Lt[j * nwave + bindx[6]] = (Lt7 * Lt7 + 0.04) / Lt7;
    }

    return (0);
}


/* ------------------------------------------------------ */
/* closel1_mos() - closes the level 1 HDF file        */
/*                                                        */

/* ------------------------------------------------------ */
int closel1_mos(filehandle *l1file) {
    free(array_data);
    free(lon);
    free(lat);
    free(solz);
    free(sola);
    free(senz);
    free(sena);
    free(image);
    free(yeararr);
    free(dayarr);
    free(msecarr);

    return (0);
}


/* ------------------------------------------------------ */
/* ReadSDS() - reads "Data-Set-4". If it is empty, reads  */
/*             "Data-Set-6".                              */
/*                                                        */

/* ------------------------------------------------------ */
int ReadSDS(filehandle *l1file) {
    int32_t sd_id, sds_id, status;
    int32_t sds_index, rank, nt, dims[HDF4_UTILS_MAX_DIM], nattrs;
    int32_t start[2], edges[2];
    char name[HDF4_UTILS_MAX_NAME];

    /* Open the file and initiate the SD interface */
    sd_id = SDstart(l1file->name, DFACC_RDONLY);

    /* Get the "Data-Set-4" SDS index */
    sds_index = SDnametoindex(sd_id, "Data-Set-4");

    /* Check that the "Data-Set-4" SDS exists; if not,
    open the "Data-Set-6" SDS */
    if (sds_index == -1) {
        sds_index = SDnametoindex(sd_id, "Data-Set-6");
        if (sds_index == -1) {
            printf("-E- %s:  Error seeking SDS\n", __FILE__);
            return (1);
        }
    }

    /* Select the SDS */
    sds_id = SDselect(sd_id, sds_index);

    /* Verify the characteristics of the array */
    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);

    /* Define the location, pattern and size of the data to read */
    start[0] = start[1] = 0;
    edges[0] = dims[0];
    edges[1] = dims[1];

    /* Read the array */
    status = SDreaddata(sds_id, start, NULL, edges, (VOIDP) array_data);

    /* Terminate access to the array */
    status = SDendaccess(sds_id);

    /* Terminate access to the SD interface and close the file */
    status = SDend(sd_id);

    return (status);
}

/* ------------------------------------------------------ */
/* GetTime - gets the time from a string in the file      */
/*           description                                  */
/*                                                        */

/* ------------------------------------------------------ */
int32_t GetTime(char arr[MAXLEN]) {
    int i;
    int length, ch;
    int32_t itime;
    char tstr[5];

    RemoveWhitespace(arr);

    length = 0;
    ch = arr[length];
    while (isalnum(ch)) {
        length = length + 1;
        ch = arr[length];
    }
    tstr[0] = '\0';
    strncat(tstr, arr, length);
    itime = atol(tstr);
    for (i = 0; i < MAXLEN - length - 2; i++) {
        arr[i] = arr[i + length + 1];
    }
    arr[i] = '\0';
    return (itime);
}

/* ------------------------------------------------------ */
/* TimeToSec - calculates number of seconds since         */
/*             Jan 1,1968                                 */
/*                                                        */

/* ------------------------------------------------------ */
int32_t TimeToSec(int32_t yr, int32_t mo, int32_t dm, int32_t hr, int32_t mn, int32_t sec) {
    int i;
    int32_t totalSec = 0;
    int mdays[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243,
        273, 304, 334};

    /* get T68 seconds to 1st day of year */
    /* loop thru years from 1968 */
    for (i = 1968; i < yr; i++) {
        if (LeapChk(i) == -1) {
            totalSec = totalSec + 3600 * 24 * 366; /* leap year */
        } else {
            totalSec = totalSec + 3600 * 24 * 365; /* non-leap year */
        }
    }

    /* add additional seconds */
    totalSec = totalSec + mdays[mo - 1]*86400 +
            (dm - 1)*86400 + hr * 3600 + mn * 60 + sec;

    /* adjust for leap years */
    if ((LeapChk(yr) == -1) && (mo > 2)) {
        totalSec = totalSec + 86400;
    }
    return (totalSec);
}

/* ------------------------------------------------------ */
/* GetYearDayMsec() - sets the yeararr, dayarr            */
/*                    msecarr arrays                      */
/*                                                        */

/* ------------------------------------------------------ */
int GetYearDayMsec(double dtime[NS]) {
    int i;
    int32_t yr, sec;
    int mdays[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243,
        273, 304, 334};

    /* if the year doesn't change, set all yeararr elements the same */
    if ((endyr - begyr) == 0) {
        for (i = 0; i < NS; i++) {
            yeararr[i] = begyr;
        }
    } else {
        /* loop through years from 1968 */
        for (i = 0; i < NS; i++) {
            yr = 1968;
            while (TimeToSec(yr, 1, 0, 0, 0, 0) < (int) dtime[i]) {
                yr++;
            }
            yeararr[i] = yr - 1;
        }
    }

    /* if the day doesn't change, set all dayarr elements the same */
    if (((endmo - begmo) == 0) && ((enddm - begdm) == 0)) {
        for (i = 0; i < NS; i++) {
            dayarr[i] = mdays[begmo - 1] + begdm;
        }
    } else {
        /* get seconds up to this year */
        for (i = 0; i < NS; i++) {
            sec = TimeToSec(yeararr[i], 1, 0, 0, 0, 0);
            /* seconds to convert into day is the difference dtime[i] - sec */
            dayarr[i] = (int32_t) (dtime[i] - sec) / (86400);
        }
    }

    /* set msecarr */
    for (i = 0; i < NS; i++) {
        /* find seconds to beginning of the day */
        sec = TimeToSec(yeararr[i], 1, dayarr[i], 0, 0, 0);
        msecarr[i] = ((int) dtime[i] - sec)*1000 +
                (dtime[i] - (double) ((int) dtime[i]))* 1000;
    }
    return (0);
}

/* ------------------------------------------------------ */
/* LeapChk() - checks if the year is a leap year:         */
/*             returns -1 if leap year                    */
/*             returns 0 if non-leap year                 */
/*                                                        */

/* ------------------------------------------------------ */
int LeapChk(int32_t yr) {
    int year, leap;

    year = yr;
    leap = 0;
    if ((year % 4) == 0) {
        leap = -1;
    }
    if ((year % 100) == 0) {
        leap = 0;
    }
    if ((year % 400) == 0) {
        leap = -1;
    }
    return (leap);
}

/* ------------------------------------------------------ */
/* RemoveWhitespace() - removes leading whitespace        */
/*                      characters from an array          */
/*                                                        */

/* ------------------------------------------------------ */
int RemoveWhitespace(char arr[MAXLEN]) {
    int i, j, ch, length;

    i = 0;
    ch = (int) arr[0];
    while (!isalnum(ch) && i++ < MAXLEN) {
        ch = (int) arr[i];
    }
    length = MAXLEN - i;
    for (j = 0; j < length; j++) {
        arr[j] = arr[j + i];
    }

    return (0);
}


/* ------------------------------------------------------ */
/* bilinear - bilinearly interpolates a set of reference  */
/*            points.                                     */
/*                                                        */
/* NAME  I/O   DESCRIPTION                                */
/* ----  ---   -----------                                */
/*  p     I    2 x 2 array holding the 4 tabulated values */
/*  x     I    x coordinates to interpolate at            */
/*  y     I    y coordinates to interpolate at            */
/*  z     O    array of interpolated values               */
/*                                                        */

/* ------------------------------------------------------ */
int bilinear(float p[PSIZE][PSIZE], float x[NS][NP],
        float y[NS][NP], float z[NS][NP]) {
    /* Cleanup and fix bilinear routine JMG  06/08/01 */

    float dx1[NS][NP], dy1[NS][NP];

    int i, j, a, b, c, d;

    /* create the arrays with elements 1 - (element from dx)
    or 1- (element from dy) */
    for (i = 0; i < NS; i++) {
        for (j = 0; j < NP; j++) {
            dx1[i][j] = 1. - x[i][j];
            dy1[i][j] = 1. - y[i][j];
        }
    }

    /* calculate interpolated values */
    for (i = 0; i < NS; i++) {
        for (j = 0; j < NP; j++) {

            a = 0;
            b = 0;
            c = 1;
            d = 1;

            z[j][i] = p[a][b] * dx1[i][j] * dy1[i][j] +
                    p[a][c] * dx1[i][j] * y[i][j] +
                    p[d][b] * x[i][j] * dy1[i][j] +
                    p[d][c] * x[i][j] * y[i][j];
        }
    }
    return (0);
}

int rd_gscan(void) {
    FILE *fp;
    char line[255];
    int iband = 0;
    char *tmp_str;
    char file[1024];
    int n;

    if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
        printf("OCDATAROOT environment variable is not defined.\n");
        return (1);
    }
    strcpy(file, tmp_str);
    strcat(file, "/mos/cal/relgain");
    if ((fp = fopen(file, "r")) == NULL) {
        printf("Error on opening the MOS relative cal file - %s\n", file);
        return (1);
    }

    while ((fgets(line, 255, fp) != NULL) && (iband < 8)) {

        /* skip the comment or blank line */
        if (line[0] == '#' || line[0] == ' ' || line[0] == '\0' || line[0] == '\n')
            continue;

        n = sscanf(line, "%f, %f, %f, %f",
                &gscan[iband][0], &gscan[iband][1],
                &gscan[iband][2], &gscan[iband][3]);

        if (n != 4) {
            printf("Error in format of %s\n", file);
            return (1);
        }

        iband++;

    }

    return (0);
}
