/* ============================================================================ */
/* module l1_aci_hdf.c - functions to read Pathfinder Class data for MSL12      */
/*                                                                              */
/* Modified By: S. Walsh, RSMAS, January 2008 to read Pathfinder                */
/*                                                                              */
/* ============================================================================ */

#include "l1_aci.h"
#include <math.h>
#include "l1.h"
#include "rawcal.h"
#include "runcal.h"
#include <hdf4utils.h>

#include <hdf.h>
#include <mfhdf.h>

#define NAVBANDS  5
#define NVIRBANDS 2
#define NTIRBANDS 3

#define TEMP_REF 273.15  /* convert Bt's from Deg K to Deg C */

#define NRADTAB 100     /* Radiance to temperature table size */

#define SDNAME "Channel_"
//#define SDNAME ""

static int32_t sd_id;
static int32_t sd_id_g;
static int32_t spix = 0;
//static int32_t sscan  = 0;  /* start scan in original granule */

static float T_LOW = 185.; /* Limit of 180K in etbpsub.f */
static float T_HIGH = 370.; /* Limit of 375K in etbpsub.f */

static int pdbg1 = -1;
static int ldbg1 = -1;
//static int pdbg1 = 174;	/* 2003131211434.N16 */
//static int ldbg1 = 1952;
//static int pdbg1 = 141; /* 2009200225057.N16 */
//static int ldbg1 = 122;

float etinvert_(int32_t *ich, float *radi);
float etintegrate_(int32_t *ich, float *etemp);
void etloadresp_(int32_t *lin, char *cal);
int32_t avconsh_(int32_t *lunin, int32_t *npix, int32_t *jday);
int32_t avlooph_(float *satz, float *solz, float *delphi, float *rayly,
        float *aersol, float *aglint);

/* ----------------------------------------------------------------------------------- */
/* openl1_aci_hdf() - opens a ACI L1B file for reading.                            */
/*                                                                                     */
/* B. Franz, SAIC, February 2003.                                                      */

/* ----------------------------------------------------------------------------------- */
int openl1_aci(filehandle *file) {
    int32_t npix;
    int32_t nscan;
    const char xsatida[10];

    int32_t sds_id;
    int32_t rank;
    int32_t dims[3];
    int32_t type;
    int32_t numattr;

    /* Open the HDF input file */
    sd_id = SDstart(file->name, DFACC_RDONLY);
    if (sd_id == FAIL) {
        fprintf(stderr, "-E- %s line %d: SDstart(%s, %d) failed.\n",
                __FILE__, __LINE__, file->name, DFACC_RDONLY);
        return (HDF_FUNCTION_ERROR);
    }

    /* Get pixel and scan dimensions */
    sds_id = SDselect(sd_id, SDnametoindex(sd_id, (SDNAME"1")));
    if (SDgetinfo(sds_id, NULL, &rank, dims, &type, &numattr) == -1) {
        fprintf(stderr, "-E- %s line %d: error getting dimension info.\n",
                __FILE__, __LINE__);
        return (HDF_FUNCTION_ERROR);
    }
    npix = dims[1];
    /* can't use dims for nscan because avhrr file is created before the */
    /*   actual number of scans is known - because missing scans need to */
    /*   be filled in */
    //    nscan = dims[0];
    if (getHDFattr(sd_id, "Number of scans", "", (VOIDP) & nscan) != 0) {
        printf("-E- %s line %d: Error reading Number of scans attribute.\n",
                __FILE__, __LINE__);
        return (1);
    }

    if (getHDFattr(sd_id, "Orbit Number", "", (VOIDP) & file->orbit_number)
            != 0) {
        printf("-E- %s line %d: Error reading Orbit attribute.\n",
                __FILE__, __LINE__);
        return (1);
    }

    /* Open the HDF geolocation input file.*/

    sd_id_g = SDstart(file->geofile, DFACC_RDONLY);
    if (sd_id_g == FAIL) {
        printf("Error opening geolocation file.\n");
        fprintf(stderr, "-E- %s line %d: SDstart(%s, %d) failed.\n",
                __FILE__, __LINE__, file->geofile, DFACC_RDONLY);
        return (HDF_FUNCTION_ERROR);
    }

    if (getHDFattr(sd_id, "Satellite", "", (VOIDP) xsatida) != 0) {
        printf("-E- %s line %d: Error reading Satellite attribute.\n",
                __FILE__, __LINE__);
        return (1);
    }

    if (getHDFattr(sd_id, "spatial_resolution", "", (VOIDP) & file->spatialResolution) != 0) {
        printf("-W- %s line %d: Error reading spatial_resolution attribute.\n",
                __FILE__, __LINE__);
    }
    file->subsensorID = satname2xsatid(xsatida);

    file->npix = npix;
    file->ndets = 1;
    file->nscan = nscan;
    file->sd_id = sd_id;

    return (LIFE_IS_GOOD);
}

/* ----------------------------------------------------------------------------------- */
/* readl1_aci_hdf() - reads 1 line (scan) from a ACI ingested CLASS file, loads l1rec.        */
/*                                                                                     */

/* ----------------------------------------------------------------------------------- */
int readl1_aci(filehandle *file, int32_t scan, l1str *l1rec) {
    static int firstCall = 1;
    static int firstScan = 1;

    static int16_t *ch1 = NULL;
    static int16_t *ch2 = NULL;
    static int16_t *ch3 = NULL;
    static int16_t *ch4 = NULL;
    static int16_t *ch5 = NULL;

    static float m[NAVBANDS];
    static float b[NAVBANDS];

    static float rayly[2];
    static float aersol[2];
    static int year, jday, RelDay, LY, leap, month, day, hh, mm, ss;
    static int32_t ms = 0;
    static char starttime[17];

    static int startOfMonth[2][12] = {
        { 0, 31, 59, 90, 120, 151, 181, 212, 243,
            273, 304, 334},
        { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305,
            335}
    };

    /* Jon Mittaz's calibration for ch4 and ch5 */
    static double ch4_coeff[5] = {0.0351618, 1.24004, 0.979493, 8.6071114e-06,
        8.6071114e-06};
    static double ch5_coeff[5] = {0.0326272, 1.06256, 0.986256, 4.0256298e-06,
        4.0256298e-06};

    static char ingester[5];

    int16_t *data = NULL;
    int32_t npix = (int32_t) file->npix;
    int32_t nbands = (int32_t) file->nbands;
    int32_t nbandsir = (int32_t) file->nbandsir;
    int32_t detnum = 1;
    int32_t ich, iw, ip, ipb;
    int32_t ii, jj;
    int32_t dims[3];
    /* lunin is used to open the noaa*.cal file.  It's an almost random
     * number chosen by searching for 'open' in the other fortran routines
     * and selecting an unused number */
    int32_t lunin = 23;

    static float radinc[3];
    static float radbas[3];
    static float radinv[3][NRADTAB];
    static float teminc[3];
    static float tembas[3];
    static float teminv[3][NRADTAB];

    static int vbcksc[30];
    static int vspace[50];
    float rspc[3];
    float radprt[3];
    float avgbck[3];
    float avgspc[3];
    float64 sumbck[2][3];
    float64 sumspc[2][5];

    char scaletype[128];
    char startdate[15];
    char satnam[12];

    static RAW_CAL *rawrec;
    static RUN_CAL *runrec;

    static int32_t *scantime;

    char tmpstr[80];

    int32_t status;
    int32_t irun, jrun, idif, iraw, jraw, ni;
    float ai;
    float ascan;
    float percnt;
    float trad, ttem, tnlc;
    float tcor;
    float prtemp;

    if (firstCall) {

        /* initialize index arrays */

        firstCall = 0;

        jj = 0;
        for (ii = 0; ii < 50; ii++) {
            vspace[ii] = jj;
            jj += 5;
            if (jj >= 50) {
                jj = jj - 50 + 1;
            }
        }
        jj = 0;
        for (ii = 0; ii < 30; ii++) {
            vbcksc[ii] = jj;
            jj += 3;
            if (jj >= 30) {
                jj = jj - 30 + 1;
            }
        }

        if ((ch1 = (int16_t *) calloc(npix, sizeof (int16_t))) == NULL) {
            printf("-E- %s line %d: Error allocating data space.\n",
                    __FILE__, __LINE__);
            return (1);
        }

        if (getHDFattr(sd_id, "Slope", SDNAME"1", (VOIDP) & m[0]) != 0) {
            printf("-E- %s line %d: Error reading Channel_1 Slope attribute.\n",
                    __FILE__, __LINE__);
            return (1);
        }

        if (getHDFattr(sd_id, "Intercept", SDNAME"1", (VOIDP) & b[0]) != 0) {
            printf(
                    "-E- %s line %d: Error reading Channel_1 Intercept attribute.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        if (getHDFattr(sd_id, "Scale_type", SDNAME"1", (VOIDP) & scaletype)
                != 0) {
            printf(
                    "-E- %s line %d: Error reading Channel_1 Scale_type attribute.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        if (strcmp(scaletype, "y= Slope * x + Intercept;") != 0) {
            printf("-E- %s line %d: Channel_1 Scale_type must be linear.\n",
                    __FILE__, __LINE__);
            return (1);
        }

        if ((ch2 = (int16_t *) calloc(npix, sizeof (int16_t))) == NULL) {
            printf("-E- %s line %d: Error allocating data space.\n",
                    __FILE__, __LINE__);
            return (1);
        }

        if (getHDFattr(sd_id, "Slope", SDNAME"2", (VOIDP) & m[1]) != 0) {
            printf("-E- %s line %d: Error reading Channel_2 Slope attribute.\n",
                    __FILE__, __LINE__);
            return (1);
        }

        if (getHDFattr(sd_id, "Intercept", SDNAME"2", (VOIDP) & b[1]) != 0) {
            printf(
                    "-E- %s line %d: Error reading Channel_2 Intercept attribute.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        if (getHDFattr(sd_id, "Scale_type", SDNAME"2", (VOIDP) & scaletype)
                != 0) {
            printf(
                    "-E- %s line %d: Error reading Channel_2 Scale_type attribute.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        if (strcmp(scaletype, "y= Slope * x + Intercept;") != 0) {
            printf("-E- %s line %d: Channel_2 Scale_type must be linear.\n",
                    __FILE__, __LINE__);
            return (1);
        }

        if ((ch3 = (int16_t *) calloc(npix, sizeof (int16_t))) == NULL) {
            printf("-E- %s line %d: Error allocating data space.\n",
                    __FILE__, __LINE__);
            return (1);
        }

        if (getHDFattr(sd_id, "Slope", SDNAME"3", (VOIDP) & m[2]) != 0) {
            printf("-E- %s line %d: Error reading Channel_3 Slope attribute.\n",
                    __FILE__, __LINE__);
            return (1);
        }

        if (getHDFattr(sd_id, "Intercept", SDNAME"3", (VOIDP) & b[2]) != 0) {
            printf(
                    "-E- %s line %d: Error reading Channel_3 Intercept attribute.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        if (getHDFattr(sd_id, "Scale_type", SDNAME"3", (VOIDP) & scaletype)
                != 0) {
            printf(
                    "-E- %s line %d: Error reading Channel_3 Scale_type attribute.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        if (strcmp(scaletype, "y= Slope * x + Intercept;") != 0) {
            printf("-E- %s line %d: Channel_3 Scale_type must be linear.\n",
                    __FILE__, __LINE__);
            return (1);
        }

        if ((ch4 = (int16_t *) calloc(npix, sizeof (int16_t))) == NULL) {
            printf("-E- %s line %d: Error allocating data space.\n",
                    __FILE__, __LINE__);
            return (1);
        }

        if (getHDFattr(sd_id, "Slope", SDNAME"4", (VOIDP) & m[3]) != 0) {
            printf("-E- %s line %d: Error reading Channel_4 Slope attribute.\n",
                    __FILE__, __LINE__);
            return (1);
        }

        if (getHDFattr(sd_id, "Intercept", SDNAME"4", (VOIDP) & b[3]) != 0) {
            printf(
                    "-E- %s line %d: Error reading Channel_4 Intercept attribute.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        if (getHDFattr(sd_id, "Scale_type", SDNAME"4", (VOIDP) & scaletype)
                != 0) {
            printf(
                    "-E- %s line %d: Error reading Channel_4 Scale_type attribute.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        if (strcmp(scaletype, "y= Slope * x + Intercept;") != 0) {
            printf("-E- %s line %d: Channel_4 Scale_type must be linear.\n",
                    __FILE__, __LINE__);
            return (1);
        }

        if ((ch5 = (int16_t *) calloc(npix, sizeof (int16_t))) == NULL) {
            printf("-E- %s line %d: Error allocating data space.\n",
                    __FILE__, __LINE__);
            return (1);
        }

        if (getHDFattr(sd_id, "Slope", SDNAME"5", (VOIDP) & m[4]) != 0) {
            printf("-E- %s line %d: Error reading Channel_5 Slope attribute.\n",
                    __FILE__, __LINE__);
            return (1);
        }

        if (getHDFattr(sd_id, "Intercept", SDNAME"5", (VOIDP) & b[4]) != 0) {
            printf(
                    "-E- %s line %d: Error reading Channel_5 Intercept attribute.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        if (getHDFattr(sd_id, "Scale_type", SDNAME"5", (VOIDP) & scaletype)
                != 0) {
            printf(
                    "-E- %s line %d: Error reading Channel_5 Scale_type attribute.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        if (strcmp(scaletype, "y= Slope * x + Intercept;") != 0) {
            printf("-E- %s line %d: Channel_5 Scale_type must be linear.\n",
                    __FILE__, __LINE__);
            return (1);
        }

        if (getHDFattr(sd_id, "start date", "", (VOIDP) & startdate) != 0) {
            printf("-E- %s line %d: Error reading start date attribute.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        /* convert yyyy-mm-dd to year, day */
        sscanf(startdate, "%04d-%02d-%02d", &year, &month, &day);
        if (((year % 400) == 0) || (((year % 4) == 0) && ((year % 100) != 0))) {
            leap = 1;
            LY = 366;
        } else {
            leap = 0;
            LY = 365;
        }

        jday = startOfMonth[leap][month - 1] + day;

        if (getHDFattr(sd_id, "start time", "", (VOIDP) & starttime) != 0) {
            printf("-E- %s line %d: Error reading start time attribute.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        /* starttime is a float with a value in this form: hhmmss.ss */
        //	hh = starttime/10000.;
        //	mm = ((int)(starttime/100.) % 100);
        //	ss = ((int)(starttime-(100.0 * (mm+ (100.*hh)))) % 100);
        //	ms = (int32_t)(1000. * (starttime- (int)starttime));
        /* starttime is "hh:mm:ss.sss UTC" */
        sscanf(starttime, "%02d:%02d:%02d.%03d", &hh, &mm, &ss, &ms);

        if (getHDFattr(sd_id, "Satellite Name", "", (VOIDP) satnam) != 0) {
            printf("-E- %s line %d: Error reading Satellite Name attribute.\n",
                    __FILE__, __LINE__);
            return (1);
        }

        if (getHDFattr(sd_id, "PRTEMP", "", (VOIDP) & prtemp) != 0) {
            printf("-E- %s line %d: Error reading PRTEMP attribute.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        //printf(" I Mean baseplate temperature =%6.1fK\n",prtemp);
        //        Sue Walsh said to comment the ACI name out...
        //        strcpy(ingester, "ACI ");
        //        if (getDims(sd_id, "ACI Raw Running Calibration", dims) != 0) {
        //            printf(
        //                    "-E- %s line %d: Error reading ACI Raw Running Calibration dim.\n",
        //                    __FILE__, __LINE__);
        //            return (1);
        //        }
        /* allocate number of bytes needed for data in this file */
        if ((rawrec = (RAW_CAL *) malloc(dims[0])) == NULL) {
            printf("-E- %s line %d: Error allocating rawcal data space.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        strcpy(tmpstr, ingester);
        strcat(tmpstr, "Raw Running Calibration");
        status = rdSDS(sd_id, tmpstr, 0, 0, 0, 0, (VOIDP) rawrec);
        if (status != 0) {
            printf("-E- %s line %d: Error reading %s.\n",
                    __FILE__, __LINE__, tmpstr);
            return (1);
        }

        strcpy(tmpstr, ingester);
        strcat(tmpstr, "Running Calibration");
        if (getDims(sd_id, tmpstr, dims) != 0) {
            printf("-E- %s line %d: Error reading %s dim.\n",
                    __FILE__, __LINE__, tmpstr);
            return (1);
        }
        /* allocate the number of bytes needed for the data in this file */
        if ((runrec = (RUN_CAL *) malloc(dims[0])) == NULL) {
            printf("-E- %s line %d: Error allocating runcal data space.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        status = rdSDS(sd_id, tmpstr, 0, 0, 0, 0, (VOIDP) runrec);

        /* Work-around for single-scan time errors */
        //        if (TAIsec < 0.0) {
        //            printf("-W- %s: bad time in geolocation file at frame %d, using previous.\n",
        //                   __FILE__,(int)frame);
        //            TAIsec = lastTAIsec;
        //        } else
        //            lastTAIsec = TAIsec;
        //
        //        usec = TAIsec + 725846400.0;
        //        unix2yds(usec,&year,&day,&dsec);
        if (getDims(sd_id, "Scan start time", dims) != 0) {
            printf("-E- %s line %d: Error reading Scan start time dim.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        if (dims[0] != l1rec->l1file->nscan) {
            printf(
                    "-E- %s line %d: Error: Number of Scan start times (%d) must match Number of Scans (%d).\n",
                    __FILE__, __LINE__, dims[0], l1rec->l1file->nscan);
            return (1);
        }
        if ((scantime = (int32_t *) malloc(dims[0] * sizeof (int32_t))) == NULL) {
            printf("-E- %s line %d: Error allocating scantime data space.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        status = rdSDS(sd_id, "Scan start time", 0, 0, 0, 0, (VOIDP) scantime);
        if (status != 0) {
            printf("-E- %s line %d: Error reading scan start time\n",
                    __FILE__, __LINE__);
            return (1);
        }

        //	/* degc option is now ONLY for sstref in nlsst* equations */
        //	if (file->input->degc == 1) {
        //	    /* Convert to Deg C for calc */
        //	    rabias = TEMP_REF;
        //	} else {
        //	    /* Do calc in Deg K */
        //	    rabias = 0.;
        //	}

        switch (file->subsensorID) {
        case NO07:
            printf("I Applying NOAA-7 visible degradation.\n");
            m[0] = m[0] / (0.916 - 0.049 * ((float) year + ((float) jday / (float) LY) - 1981.5)
                    + 0.0050 * pow((float) year + ((float) jday / (float) LY) - 1981.5, 2));
            b[0] = b[0] / (0.916 - 0.049 * ((float) year + ((float) jday / (float) LY) - 1981.5)
                    + 0.0050 * pow((float) year + ((float) jday / (float) LY) - 1981.5, 2));
            m[1] = m[1] / (0.882 - 0.080 * ((float) year + ((float) jday / (float) LY) - 1981.5)
                    + 0.0110 * pow((float) year + ((float) jday / (float) LY) - 1981.5, 2));
            b[1] = b[1] / (0.882 - 0.080 * ((float) year + ((float) jday / (float) LY) - 1981.5)
                    + 0.0110 * pow((float) year + ((float) jday / (float) LY) - 1981.5, 2));
            break;
        case NO09:
            printf("I Applying NOAA-9 visible degradation.\n");
            m[0] = m[0] / (0.953 - 0.051 * ((float) year + ((float) jday / (float) LY) - 1985.0));
            b[0] = b[0] / (0.953 - 0.051 * ((float) year + ((float) jday / (float) LY) - 1985.0));
            m[1] = m[1] / (0.866 - 0.026 * ((float) year + ((float) jday / (float) LY) - 1985.0));
            b[1] = b[1] / (0.866 - 0.026 * ((float) year + ((float) jday / (float) LY) - 1985.0));
            break;
        case NO11:
            printf("I Applying NOAA-11 visible degradation.\n");
            m[0] = m[0] / (0.797 - 0.010 * ((float) year + ((float) jday / (float) LY) - 1989.0));
            b[0] = b[0] / (0.797 - 0.010 * ((float) year + ((float) jday / (float) LY) - 1989.0));
            m[1] = m[1] / (0.683 - 0.020 * ((float) year + ((float) jday / (float) LY) - 1989.0));
            b[1] = b[1] / (0.683 - 0.020 * ((float) year + ((float) jday / (float) LY) - 1989.0));
            break;
        case NO14:
            /* Dec 30, 1994 is base date */
            RelDay = ((year - 1994) * 365 + (jday - 364));
            printf("I Applying NOAA-14 visible degradation.\n");
            /* Post-Launch Calibration for Noaa-14 */
            /* C R Nagaraja Rao and Jianhua Chen   */
            /* NOAA/NESDIS                         */
            /* (from web page, March 1997)         */
            m[0] = 0.0000232 * (float) RelDay + 0.109;
            b[0] = -41.0 * m[0];
            m[1] = 0.0000373 * (float) RelDay + 0.129;
            b[1] = -41.0 * m[1];
            break;
        default:
            printf("I No visible degradation available for NOAA-%2.0d.\n",
                    file->subsensorID);
            break;

        }

        /* start running calibration */

        etloadresp_(&lunin, satnam);
        /* for debug / info
         printf(" t_low  %10.4f\n",T_LOW);
         printf(" t_high %10.4f\n",T_HIGH);
         */
        for (iw = 0; iw <= 2; iw++) {
            ich = iw + 3;
            radinc[iw] = (float) ((log10(etintegrate_(&ich, &T_HIGH))
                    - log10(etintegrate_(&ich, &T_LOW)))
                    / (float) (NRADTAB - 1));
            radbas[iw] = (float) (log10(etintegrate_(&ich, &T_LOW)));
            /* for debug / info
             printf(" Ch %1d radbas %14.6E radinc %14.6E\n",ich,
             radbas[iw],radinc[iw]);
             printf(" Ch Idx radinv (LOG)\n");
             */
            radinv[iw][0] = T_LOW;
            /* for debug / info
             printf(" %2d %3d %14.6E\n",ich, 1, radinv[iw][0]);
             */
            for (ip = 1; ip < NRADTAB; ip++) {
                ttem = pow(10., radbas[iw] + radinc[iw] * (float) (ip));
                radinv[iw][ip] = etinvert_(&ich, &ttem);
                /* for debug / info
                 printf(" %2d %3d %14.6E\n",ich, ip+1, radinv[iw][ip]);
                 */
            }
            teminc[iw] = (T_HIGH - T_LOW) / (float) (NRADTAB - 1);
            tembas[iw] = T_LOW;
            /* for debug / info
             printf(" Ch %1d tembas %14.6E teminc %14.6E\n",ich,
             tembas[iw],teminc[iw]);
             printf(" Ch Idx teminv\n");
             */
            teminv[iw][0] = etintegrate_(&ich, &radinv[iw][0]);
            /* for debug / info
             printf(" %2d %3d %14.6E\n",ich, 1, teminv[iw][0]);
             */
            for (ip = 1; ip < NRADTAB; ip++) {
                ttem = tembas[iw] + (teminc[iw] * (float) (ip));
                teminv[iw][ip] = etintegrate_(&ich, &ttem);
                /* for debug / info
                 printf(" %2d %3d %14.6E\n",ich, ip+1, teminv[iw][ip]);
                 */
            }
        }

        status = avconsh_(&lunin, &npix, &jday);

        /* end first */
    }

    /*    if (scan == ldbg1) {			*/
    /*	printf(" got to scan 1780\n");		*/
    /*    }						*/
    /* check to see if this scan line is a different day than the previous */
    if (scan > 0 && (scantime[scan] < scantime[scan - 1])) {
        /* crossed to the next day */
        jday = jday + 1;
        if (jday > LY) {
            /* crossed to the next year */
            jday -= LY;
            year = year + 1;
        }
    }

    l1rec->scantime = yds2unix(year, jday, (double) (scantime[scan] / 1.e3));

    /* Get position and path geometry */
    READ_SDS_ID(sd_id_g, "Longitude", l1rec->lon, scan, spix, 0, 0, 1, npix, 1,
            1);
    READ_SDS_ID(sd_id_g, "Latitude", l1rec->lat, scan, spix, 0, 0, 1, npix, 1,
            1);
    READ_SDS_ID(sd_id_g, "Solar Zenith Angle", l1rec->solz, scan, spix, 0, 0, 1,
            npix, 1, 1);
    READ_SDS_ID(sd_id_g, "Sensor Zenith Angle", l1rec->senz, scan, spix, 0, 0,
            1, npix, 1, 1);
    READ_SDS_ID(sd_id_g, "Relative Azimuth Angle", l1rec->delphi, scan, spix, 0,
            0, 1, npix, 1, 1);

    for (ip = 0; ip < npix; ip++) {

        l1rec->pixnum[ip] = spix + ip;

        if (l1rec->lon[ip] < -181.0 || l1rec->lon[ip] > 181.0
                || l1rec->lat[ip] < -91.0 || l1rec->lat[ip] > 91.0)
            l1rec->navfail[ip] = 1;
    }

    /* Read L1B data, scale to radiance, and copy relevant bands to L1 record  */
    READ_SDS_ID(sd_id, SDNAME"1", ch1, scan, spix, 0, 0, 1, npix, 1, 1);
    READ_SDS_ID(sd_id, SDNAME"2", ch2, scan, spix, 0, 0, 1, npix, 1, 1);
    READ_SDS_ID(sd_id, SDNAME"3", ch3, scan, spix, 0, 0, 1, npix, 1, 1);
    READ_SDS_ID(sd_id, SDNAME"4", ch4, scan, spix, 0, 0, 1, npix, 1, 1);
    READ_SDS_ID(sd_id, SDNAME"5", ch5, scan, spix, 0, 0, 1, npix, 1, 1);

    /* This check should somehow check the scan line number to see if it's
     * the first or last boxsiz/2 lines at the beginning or end of the file.
     *
     * I'm assuming the boxsiz is 3 so only the first line doesn't get
     * processed and that this routine is called for each scan line in order
     * from first to last.
     *
     */

    if (firstScan) {
        firstScan = 0;
        /* what else happens on first line? */
        return (LIFE_IS_GOOD);
    }

    ascan = (float) scan + 1.0; /* ascan is one based */
    irun = -1;
    idif = 99999;

    /* numcal and numraw should always be the same */
    for (ip = 0; ip < runrec->numcal; ip++) {
        if (abs((int)(runrec->runcal[ip].cenlin - ascan)) < idif) {
            idif = abs((int)(runrec->runcal[ip].cenlin - ascan));
            irun = ip;
        }
    }

    if (irun == -1 || idif > runrec->intrvl) {
        printf("-W- %s line %d: Hey, I can't find the RUNCAL entry!\n",
                __FILE__, __LINE__);
        printf("-W- scan: %f;  first center scan: %f; last: %f", ascan,
                runrec->runcal[0].cenlin,
                runrec->runcal[runrec->numcal - 1].cenlin);
        if (irun == -1) {
            if (ascan < runrec->runcal[1].cenlin) {
                irun = 0;
            } else {
                irun = runrec->numcal - 1;
            }
        }
    }
    if (ascan > runrec->runcal[irun].cenlin) {
        if (irun < runrec->numcal - 1) {
            jrun = irun + 1; /* next available entry */
        } else {
            jrun = irun; /* duplicate */
        }
    } else {
        jrun = irun; /* next entry */
        if (irun > 0) {
            irun = irun - 1; /* previous entry */
        }
    }

    iraw = irun; /* should be the same */
    jraw = jrun;

    if (runrec->runcal[irun].cenlin != runrec->runcal[jrun].cenlin) {
        percnt = (ascan - runrec->runcal[irun].cenlin)
                / (runrec->runcal[jrun].cenlin - runrec->runcal[irun].cenlin);
        if (percnt < 0.0) {
            percnt = 0.0; /* earlier than first set */
        } else if (percnt > 1.0) {
            percnt = 1.0; /* later than last set */
        }
    } else {
        percnt = 0.0; /* use first set, duplicates */
    }

    /* Decommutate raw telemetry data for sensor calibration */

    iw = 0; /* Force initialization for group */
    for (ii = 1; ii <= 30; ii++) {
        if (ii > iw * 10) {
            iw++; /* Initialize for next group */
            for (jj = 0; jj < 2; jj++) {
                sumbck[jj][iw - 1] = 0.0;
            }
        }
        sumbck[0][iw - 1] += rawrec->rawcal[iraw].bckscn[0][vbcksc[ii - 1]]
                * (1.0 - percnt)
                + rawrec->rawcal[jraw].bckscn[0][vbcksc[ii - 1]] * (percnt);
        sumbck[1][iw - 1] += rawrec->rawcal[iraw].bckscn[1][vbcksc[ii - 1]]
                * (1.0 - percnt)
                + rawrec->rawcal[jraw].bckscn[1][vbcksc[ii - 1]] * (percnt);
    }
    for (iw = 0; iw < 3; iw++) {
        if (sumbck[0][iw] > 0.0) {
            avgbck[iw] = sumbck[1][iw] / sumbck[0][iw];
        } else {
            avgbck[iw] = 0.0;
        }
    }

    iw = 0; /* Force initialization for group 1 */
    for (ii = 1; ii <= 50; ii++) {
        if (ii > iw * 10) {
            iw++; /* Initialize for next group */
            for (jj = 0; jj < 2; jj++) {
                sumspc[jj][iw - 1] = 0.0;
            }
        }
        sumspc[0][iw - 1] += rawrec->rawcal[iraw].space[0][vspace[ii - 1]]
                * (1.0 - percnt)
                + rawrec->rawcal[jraw].space[0][vspace[ii - 1]] * (percnt);
        sumspc[1][iw - 1] += rawrec->rawcal[iraw].space[1][vspace[ii - 1]]
                * (1.0 - percnt)
                + rawrec->rawcal[jraw].space[1][vspace[ii - 1]] * (percnt);
    }
    for (iw = 0; iw < 3; iw++) {
        if (sumspc[0][iw + 2] > 0.0) {
            avgspc[iw] = sumspc[1][iw + 2] / sumspc[0][iw + 2];

        } else {
            avgspc[iw] = 0.0;
        }
    }

    prtemp = (1.0 - percnt) * runrec->runcal[irun].prtemp
            + (percnt) * runrec->runcal[jrun].prtemp;

    for (iw = 0; iw < 3; iw++) {
        /* Convert temperature to radiance */
        ttem = prtemp;
        if (ttem < tembas[iw])
            ttem = tembas[iw];
        ai = (ttem - tembas[iw]) / teminc[iw];
        if (ai < 0)
            ai = 0;
        else if (ai > NRADTAB)
            ai = NRADTAB - 2;
        ni = ai; /* whole part */
        ai = ai - ni; /* fractional part */
        radprt[iw] = (1.0 - ai) * teminv[iw][ni] + ai * teminv[iw][ni + 1];
    }

    rspc[0] = 0.0;
    switch (file->subsensorID) {
    case NO07:
        rspc[1] = -5.86;
        rspc[2] = -4.95;
        break;
    case NO09:
        rspc[1] = -5.53; /* Based on Table 2A; edited data, Apr 93 */
        rspc[2] = -3.06;
        break;
    case NO10:
        rspc[1] = -7.29; /* Based on Rao memo, Jan 4, 94 */
        rspc[2] = -7.29;
        break;
    case NO11:
        rspc[1] = -8.05; /* Based on Rao memo, Jan 4, 94 */
        rspc[2] = -3.51;
        break;
    case NO12:
        rspc[1] = -5.51; /* Based on Rao memo, Jan 4, 94 */
        rspc[2] = -2.51;
        break;
    case NO14:
        rspc[1] = -4.05; /* Based on S. Brown memo, Nov 22, 94 */
        rspc[2] = -2.29;
        break;
    case NO15:
        rspc[1] = -4.50; /* Based on Draft KLM, app D.1-14 */
        rspc[2] = -3.61;
        break;
    case NO16:
        rspc[1] = -2.467; /* Based on Draft KLM, app D.2-15 */
        rspc[2] = -2.009;
        break;
    case NO17:
        rspc[1] = -8.55; /* Based on Draft KLM, app D.3-2 */
        rspc[2] = -3.97;
        break;
    case NO18:
        rspc[1] = -5.53; /* Based on Draft KLM, app D.4-2 */
        rspc[2] = -2.22;
        break;
    case NO19:
        rspc[1] = -5.49; /* Based on Draft KLM, app D.6-2 */
        rspc[2] = -3.39;
        break;
    default:
        rspc[1] = -99.;
        rspc[2] = -99.;
    }

    /* noaa-16,17,18,19 have channel 3a and 3b */
    /*	3a is 1.6 um during the day */
    /*  3b is 3.75 um at night */
    /* earlier avhrr satellites just have the 3.75 um (day and night) */
    /* so process channel 3 as visible during the day, and IR at night */
    for (iw = 0; iw < nbands + 1; iw++) {

        switch (iw) {
        case 0:
            data = ch1;
            break;
        case 1:
            data = ch2;
            break;
        case 2:
            data = ch3;
            break;
        }

        for (ip = 0; ip < npix; ip++) {

            ipb = ip * nbands + iw;
            l1rec->Lt[ipb] = 0.0;

            /* check for sentinel values and flag as appropriate */
            if (data[ip] >= 32000) {

                switch (data[ip]) {

                default:
                    l1rec->hilt[ip] = 1;
                    l1rec->Lt[ipb] = 1000.0;
                    break;
                }

            } else if (l1rec->solz[ip] < SOLZNIGHTA) {
                /* only do channels 1 and 2 for daytime pixels */
                /* only do channel 3 for daytime for noaa-16,17,18,19 */
                if (iw == 2 && file->subsensorID != NO16 && file->subsensorID != NO17 &&
                        file->subsensorID != NO18 && file->subsensorID != NO19) {
                    /* this channel doesn't make the pixel bad */
                    //		  l1rec->hilt[ip] = 1;
                    l1rec->Lt[ipb] = 1000.0;
                } else {
                    status = avlooph_(&l1rec->senz[ip], &l1rec->solz[ip],
                            &l1rec->delphi[ip], rayly, aersol,
                            &l1rec->glint_coef[ip]);

                    trad = (data[ip] * m[iw]) + b[iw];
                    trad -= rayly[iw];

                    if (aersol[iw] != 0.0) {
                        trad /= aersol[iw];
                    }
                    /* aermlt in pathnlc is 1.0 so we don't need it here */

                    l1rec->Lt[ipb] = trad;
                }
            }
        }
    }

    for (iw = 0; iw < nbandsir; iw++) {

        switch (iw) {
        case 0:
            data = ch3;
            break;
        case 1:
            data = ch4;
            break;
        case 2:
            data = ch5;
            break;
        }

        for (ip = 0; ip < npix; ip++) {

            ipb = ip * NBANDSIR + iw;
            l1rec->Ltir[ipb] = 0.0;
            l1rec->Bt[ipb] = BT_LO;

            /* don't do channel 3 for daytime for noaa-16,17,18,19 */
            if (iw == 0 && l1rec->solz[ip] < SOLZNIGHTA
                    && (file->subsensorID == NO16
                    || file->subsensorID == NO17
                    || file->subsensorID == NO18
                    || file->subsensorID == NO19)) {
                continue;
            } else {

                /* check for sentinel values and flag as appropriate */
                if (data[ip] >= 32000) {

                    switch (data[ip]) {

                    default:
                        l1rec->hilt[ip] = 1;
                        l1rec->Lt[ipb] = 1000.0;
                        break;
                    }

                } else {
                    //	      if (ip == pdbg1 && scan == ldbg1) {
                    //		 printf(" got to 320,40\n");
                    //	      }
                    /* Miami calibration */
                    if (avgbck[iw] != avgspc[iw]) {
                        if (l1_input->newavhrrcal == 1
                                && file->subsensorID == NO16
                                && iw > 0) {
                            /* Jon Mittaz's calibration for ch4 and ch5 */
                            /* T_inst = prtemp,  R_ICT = radprt[iw] */
                            /* C_S = avgspc[iw], C_BB = avgbck[iw] */
                            /* C = data[ip] */
                            if (iw == 1) {
                                /* avhrr channel 4 Bt11  */
                                if (ip == pdbg1 && scan == ldbg1) {
                                    printf(" ch4_coeffs=%f %f %f %f %f\n",
                                            ch4_coeff[0], ch4_coeff[1],
                                            ch4_coeff[2], ch4_coeff[3],
                                            ch4_coeff[4]);
                                    printf(" T_inst, prtemp = %f \n", prtemp);
                                    printf(" R_ICT, radprt = %f \n",
                                            radprt[iw]);
                                    printf(" C_S, avgspc = %f \n", avgspc[iw]);
                                    printf(" C_BB, avgbck = %f \n", avgbck[iw]);
                                    printf(" C, data = %d \n", data[ip]);
                                }
                                trad = ch4_coeff[0] * (prtemp - 288.) + ch4_coeff[1]
                                        + ((((ch4_coeff[2] * radprt[iw]) - (ch4_coeff[3]
                                        * (pow(avgspc[iw] - avgbck[iw], 2.0))))
                                        / (avgspc[iw] - avgbck[iw])) * (avgspc[iw] - data[ip]))
                                        + ch4_coeff[4] * (pow(avgspc[iw] - data[ip], 2.0));
                            } else if (iw == 2) {
                                /* avhrr channel 5  Bt12 */
                                if (ip == pdbg1 && scan == ldbg1) {
                                    printf(" ch5_coeffs=%f %f %f %f %f\n",
                                            ch5_coeff[0], ch5_coeff[1],
                                            ch5_coeff[2], ch5_coeff[3],
                                            ch5_coeff[4]);
                                    printf(" T_inst, prtemp = %f \n", prtemp);
                                    printf(" R_ICT, radprt = %f \n",
                                            radprt[iw]);
                                    printf(" C_S, avgspc = %f \n", avgspc[iw]);
                                    printf(" C_BB, avgbck = %f \n", avgbck[iw]);
                                    printf(" C, data = %d \n", data[ip]);
                                }
                                trad = ch5_coeff[0] * (prtemp - 288.) + ch5_coeff[1]
                                        + ((((ch5_coeff[2] * radprt[iw]) - (ch5_coeff[3]
                                        * (pow(avgspc[iw] - avgbck[iw], 2.0))))
                                        / (avgspc[iw] - avgbck[iw])) * (avgspc[iw] - data[ip]))
                                        + ch5_coeff[4] * (pow(avgspc[iw] - data[ip], 2.0));
                            }

                        } else {
                            trad = radprt[iw] * (avgspc[iw] - data[ip])
                                    / (avgspc[iw] - avgbck[iw])
                                    + rspc[iw] * (avgbck[iw] - data[ip])
                                    / (avgbck[iw] - avgspc[iw]);
                            if (iw == 1) {
                                /* apply non-linearity correction for channel 4 */
                                switch (file->subsensorID) {
                                case NO07:
                                    tnlc = 5.7843 - 1.0754e-1 * trad + 4.8042e-4 * pow(trad, 2);
                                    break;
                                case NO09: /* Based on Table 2A; edited data, Apr 93 */
                                    tnlc = 5.24 + (.88643 - 1.) * trad + 6.033e-4 * pow(trad, 2);
                                    break;
                                case NO10: /* Based on Rao memo, Jan 4, 94 */
                                    tnlc = 5.76 + (.88428 - 1.) * trad + 5.882e-4 * pow(trad, 2);
                                    break;
                                case NO11: /* Based on Rao memo, Jan 4, 94 */
                                    tnlc = 7.21 + (.8412 - 1.) * trad + 8.739e-4 * pow(trad, 2);
                                    break;
                                case NO12: /* Based on Rao memo, Jan 4, 94 */
                                    tnlc = 5.11 + (.88929 - 1.) * trad + 5.968e-4 * pow(trad, 2);
                                    break;
                                case NO14: /* Based on S. Brown memo, Nov 22, 94 */
                                    tnlc = 3.72 + (.92378 - 1.) * trad
                                            + 3.822e-4 * pow(trad, 2);
                                    break;
                                case NO15: /* Based on Draft KLM, app D.1-14 */
                                    tnlc = 4.76 + (-0.0932) * trad + 4.524e-4 * pow(trad, 2);
                                    break;
                                case NO16: /* Based on Draft KLM, app D.2-15 */
                                    tnlc = 2.96 + (-0.05411) * trad + 2.4532e-4 * pow(trad, 2);
                                    break;
                                case NO17: /* Based on Draft KLM, app D.3-2 */
                                    tnlc = 8.22 + (-0.15795) * trad + 7.5579e-4 * pow(trad, 2);
                                    break;
                                case NO18: /* Based on Draft KLM, app D.4-2 */
                                    tnlc = 5.82 + (-0.11069) * trad + 5.2337e-4 * pow(trad, 2);
                                    break;
                                case NO19: /* Based on Draft KLM, app D.6-2 */
                                    tnlc = 5.70 + (-0.11187) * trad + 5.4668e-4 * pow(trad, 2);
                                    break;
                                default:
                                    tnlc = -99.;
                                }
                                trad += tnlc;
                            } else if (iw == 2) {
                                /* apply non-linearity correction for channel 5 */
                                switch (file->subsensorID) {
                                case NO07:
                                    tnlc = 4.4035 - 6.9038e-2 * trad + 2.5741e-4 * pow(trad, 2);
                                    break;
                                case NO09: /* Based on Table 2A; edited data, Apr 93 */
                                    tnlc = 2.42 + (.95311 - 1.) * trad + 2.198e-4 * pow(trad, 2);
                                    break;
                                case NO10: /* Based on Rao memo, Jan 4, 94 */
                                    tnlc = 5.76 + (.88428 - 1.) * trad + 5.882e-4 * pow(trad, 2);
                                    break;
                                case NO11: /* Based on Rao memo, Jan 4, 94 */
                                    tnlc = 2.92 + (.94598 - 1.) * trad + 2.504e-4 * pow(trad, 2);
                                    break;
                                case NO12: /* Based on Rao memo, Jan 4, 94 */
                                    tnlc = 1.91 + (.96299 - 1.) * trad + 1.775e-4 * pow(trad, 2);
                                    break;
                                case NO14: /* Based on S. Brown memo, Nov 22, 94 */
                                    tnlc = 2.00 + (.96194 - 1.) * trad + 1.742e-4 * pow(trad, 2);
                                    break;
                                case NO15: /* Based on Draft KLM, app D.1-14 */
                                    tnlc = 3.83 + (-0.0659) * trad + 2.811e-4 * pow(trad, 2);
                                    break;
                                case NO16: /* Based on Draft KLM, app D.2-15 */
                                    tnlc = 2.25 + (-0.03665) * trad + 1.4854e-4 * pow(trad, 2);
                                    break;
                                case NO17: /* Based on Draft KLM, app D.3-2 */
                                    tnlc = 4.31 + (-0.07318) * trad + 3.0976e-4 * pow(trad, 2);
                                    break;
                                case NO18: /* Based on Draft KLM, app D.4-2 */
                                    tnlc = 2.67 + (-0.04360) * trad + 1.7715e-4 * pow(trad, 2);
                                    break;
                                case NO19: /* Based on Draft KLM, app D.6-2 */
                                    tnlc = 3.58 + (-0.05991) * trad + 2.4985e-4 * pow(trad, 2);
                                    break;
                                default:
                                    tnlc = -99.;
                                }
                                trad += tnlc;
                            }
                        }

                        if (ip == pdbg1 && scan == ldbg1) {
                            printf(" ch %d trad = %f\n", iw + 3, trad);
                            printf(" newavhrrcal = %d\n",
                                    l1_input->newavhrrcal);
                            printf(" xsatid = %d NO16=%d\n",
                                    file->subsensorID, NO16);
                        }
                        /* convert radiance to temperature */
                        if (trad < 1.e-10)
                            trad = 1.e-10;
                        if (radinc[iw] != 0.0)
                            ai = (log10(trad) - radbas[iw]) / radinc[iw];
                        else
                            ai = NRADTAB - 1;
                        if (ai < 0.0)
                            ai = 0.0;
                        else if (ai >= NRADTAB - 1)
                            ai = NRADTAB - 2;
                        ni = ai; /* whole part */
                        ai = ai - ni; /* fractional part */
                        ttem = (1.0 - ai) * radinv[iw][ni] + ai * radinv[iw][ni + 1];
                        if (ip == pdbg1 && scan == ldbg1) {
                            printf(" ch %d ttem = %f\n", iw + 3, ttem);
                        }
                        if (l1_input->newavhrrcal
                                == 1 && file->subsensorID == NO16) {
                            /* convert Mittaz Bt's to old so we can test sst with old coeffs */
                            //		if (iw == 1) {
                            //		    /* convert Mittaz Bt11 to old so we can test sst with old coeffs */
                            //		    tcor = -10.306526 + 0.068776865*ttem - 0.00010926801*pow(ttem,2);
                            //		    /* tcor is Walton-Mittaz so add Mittaz to get Walton */
                            //		    ttem = ttem + tcor;
                            //		} else if (iw == 2) {}
                            if (iw == 2) {
                                /* use ch5 (Bt12) correction for ch4 (Bt11) also */
                                /* convert Mittaz Bt12 to old so we can test sst with old coeffs */
                                tcor = -6.9176245 + 0.038554339 * ttem - 4.7661371e-05 * pow(ttem, 2);
                                /* tcor is Walton-Mittaz so add Mittaz to get Walton */
                                ttem = ttem + tcor;
                                /* add ch5 correction to ch4 also */
                                /* ipb is ip*NBANDSIR+iw, so ch4 should be ipb-1 */
                                l1rec->Bt[ipb - 1] = l1rec->Bt[ipb - 1] + tcor;
                                if (ip == pdbg1 && scan == ldbg1) {
                                    printf(" real ch = %d, ttem = %f, tcor = %f\n", iw + 3 - 1, ttem, tcor);
                                }
                            }
                            if (ip == pdbg1 && scan == ldbg1) {
                                printf(" ch = %d, ttem = %f, tcor = %f\n", iw + 3, ttem, tcor);
                            }
                        }
                        ttem = ttem - TEMP_REF; /* convert from Deg K to Deg C */
                        l1rec->Bt[ipb] = ttem;
                    }
                    /* put top of the atmosphere counts in Ltir */
                    l1rec->Ltir[ipb] = data[ip];
                    // l1rec->prtemp[ip] = prtemp - TEMP_REF;

                }
            }
        }
    }

    l1rec->npix = file->npix;
    l1rec->detnum = (int32_t) detnum;
    l1rec->mside = 0;

    /* Convert IR bands to brightness temperature */
    //    radiance2bt(l1rec,-1);    
    return (LIFE_IS_GOOD);
}

int closel1_aci(filehandle *file) {
    if ((sd_id_g != sd_id) && SDend(sd_id_g)) {
        fprintf(stderr, "-E- %s line %d: SDend(%d) failed for file, %s.\n",
                __FILE__, __LINE__, sd_id, file->geofile);
        return (HDF_FUNCTION_ERROR);
    }
    if (SDend(sd_id)) {
        fprintf(stderr, "-E- %s line %d: SDend(%d) failed for file, %s.\n",
                __FILE__, __LINE__, sd_id, file->name);
        return (HDF_FUNCTION_ERROR);
    }

    return (LIFE_IS_GOOD);
}

const char* xsatid2name(int xsatid) {
    switch (xsatid) {
    case NO06:
        return "NO06";
    case NO07:
        return "NO07";
    case NO08:
        return "NO08";
    case NO09:
        return "NO09";
    case NO10:
        return "NO10";
    case NO11:
        return "NO11";
    case NO12:
        return "NO12";
    case NO14:
        return "NO14";
    case NO15:
        return "NO15";
    case NO16:
        return "NO16";
    case NO17:
        return "NO17";
    case NO18:
        return "NO18";
    case NO19:
        return "NO19";
    default:
        return "UNKNOWN";
    }
}

int satname2xsatid(const char* satname) {
    if (strcmp(satname, "NOA6") == 0)
        return NO06;
    else if (strcmp(satname, "NOA7") == 0)
        return NO07;
    else if (strcmp(satname, "NOA8") == 0)
        return NO08;
    else if (strcmp(satname, "NOA9") == 0)
        return NO09;
    else if (strcmp(satname, "NO10") == 0)
        return NO10;
    else if (strcmp(satname, "NO11") == 0)
        return NO11;
    else if (strcmp(satname, "NO12") == 0)
        return NO12;
    else if (strcmp(satname, "NO14") == 0)
        return NO14;
    else if (strcmp(satname, "NO15") == 0)
        return NO15;
    else if (strcmp(satname, "NO16") == 0)
        return NO16;
    else if (strcmp(satname, "NO17") == 0)
        return NO17;
    else if (strcmp(satname, "NO18") == 0)
        return NO18;
    else if (strcmp(satname, "NO19") == 0)
        return NO19;
    else
        return -999;
}

