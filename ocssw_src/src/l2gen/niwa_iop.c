
/* ---------------------------------------------------------------------------
 * niwa_iop.c - IOP algorithm: Gerald Moore, Sam Lavender, Matt Pinkerton
 *
 * Products:
 *	niwa_a:  total absorption coefficient - water absorption
 *	niwa_bb: total backscatter coefficient - water backscatter
 *
 * Algorithm citation: Pinkerton et al., 2006. 'A method for estimating
 *	inherent optical properties of New Zealand continental shelf waters
 *	from satellite ocean colour measurements', New Zealand Journal of
 *	Marine and Freshwater Research 40, pp 227-247.
 *
 * Implementation of IOP only: J.N. Schwarz, Simon Wood, NIWA, Sept. 2008
 *
 * --------------------------------------------------------------------------- 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "l12_proto.h"
#include "l2prod.h"

#include "niwa_iop.h"

/* macros */
#define radians(degrees) ((degrees) * OEL_PI / 180.0)
#define degrees(radians) ((radians) * 180.0 / OEL_PI)

/* constant sizes */
#define LBANDS  8                               /* 8 bands in LUT */
#define VBANDS  6                               /* results returned for 6 visible bands, regardless which sensor */
#define SEAWIFS_NBANDS 8
#define MODIS_NBANDS   9

/* Look-up table dimensions */
#define TH_NS   10
#define DPHI_NS 19
#define AP_NS   16
#define BP_NS   16
#define ABSIND_NS  3891200                 /* 8 bands (10.10.19.16.16) from 172800 */

#define IOPF_NIWA_BADGEOM       0x0020          /* was SPARE2 */

#define IOPmodel 1                              /* uses only 490 and 510 bands */

/* Switch on/off VERBose mode for error checking */
#define VERB 0
//#define VERB 1

/* value to return when we cannot generate an IOP value for whatever reason
 * -- make sure this matches 'l2prodstr.badData' */
const float NO_DATA = BAD_FLT;


/* LUT - IOP variables */
static int th_s_n, th_v_n, dphi_n;
static float th_s_lev[TH_NS], th_v_lev[TH_NS], dphi_lev[DPHI_NS];
static float th_s_ind[TH_NS], th_v_ind[TH_NS], dphi_ind[DPHI_NS];
static int ap_n, bp_n;
static float band_lev[LBANDS], ap_lev[AP_NS], bp_lev[BP_NS];
static float a_ind[AP_NS], b_ind[BP_NS];

/* IOP base table variables */
static float abst[AP_NS][BP_NS][LBANDS];
static float *absen; /* the pointer for absen */
static int absind;

/* mapping from our VBAND index to sensor's l2rec bands */
static const int seawifs_map[VBANDS] = {0, 1, 2, 3, 4, 5};
static const int modis_map[VBANDS] = {0, 1, 2, 3, 4, 5};
static int nbands;
static const int *band_map;

/* sensor / band specific data values */
static float aw[VBANDS], bbw[VBANDS], lambda[VBANDS];

/* structure for passing iop retrievals per pixel */
typedef struct abs_res_ {
    float a[VBANDS];
    float bb[VBANDS];
} abs_res_t;


/* ---------------------------------------------------------------------------  */

/* ---------------------------------------------------------------------------
 * interpol()
 *  Converted from IDL interpol function, but now only does linear 1-D interpolation
 *
 *  inputs:  
 *  outputs: 
 *  returns: 
 * ---------------------------------------------------------------------------
 */
static float interpol(float *v, float *x, float u, int size) {
    static const int rsize = 1;

    int s1, i, ix, m2;
    float d, r;

    /* last subs in v and x */
    m2 = size - 2;

    /* Incr or Decr X */
    s1 = ((x[1] - x[0]) >= 0) ? 1 : -1;

    /* current point */
    ix = 0;

    /* point loop */
    for (i = 0; i <= rsize; i++) {
        /* difference */
        d = s1 * (u - x[ix]);
        if (d == 0.0) {
            r = v[ix];
        } else {
            if (d > 0.0) {
                while ((s1 * (u - x[ix + 1]) > 0.0) && (ix < m2)) {
                    ix = ix + 1;
                }
            } else {
                while ((s1 * (u - x[ix]) < 0.0) && (ix > 0)) {
                    ix = ix - 1;
                }
            }
            r = v[ix] + (u - x[ix]) * (v[ix + 1] - v[ix]) / (x[ix + 1] - x[ix]);
        }
    }
    return r;
}

/* ---------------------------------------------------------------------------
 * setgeom()
 *
 *  inputs:  
 *  outputs: 
 *  returns: TRUE if viewing geometry within accepatble limits; FALSE otherwise
 * ---------------------------------------------------------------------------
 */
static int setgeom(float sun_theta, float sen_theta, float dphi) {
    /* angles passed from niwa_iop() are in radians */
    int tabind;
    int x, y, i, pos;
    float th_s_ent, th_v_ent, dphi_ent;

#if (VERB) 
    printf("setgeom starting: [%f %f %f]\n", sun_theta * 180. / OEL_PI,
            sen_theta * 180. / OEL_PI, dphi * 180. / OEL_PI);
#endif

    /* Page in IOP table */
    /* check that the angle matches index */
    th_s_ent = floor(interpol(th_s_ind, th_s_lev, sun_theta, TH_NS) + 0.5);
    th_v_ent = floor(interpol(th_v_ind, th_v_lev, sen_theta, TH_NS) + 0.5);
    if (dphi < 0.0)
        dphi = dphi + OEL_PI * 2.0;
    if (dphi > OEL_PI * 2.0)
        dphi = dphi - OEL_PI * 2.0;
    dphi_ent = floor(interpol(dphi_ind, dphi_lev, dphi, DPHI_NS) + 0.5);

    /* LUT index */
    tabind = (int) (th_s_ent * TH_NS * DPHI_NS * AP_NS * BP_NS * LBANDS)
            + (th_v_ent * DPHI_NS * AP_NS * BP_NS * LBANDS) + (dphi_ent * AP_NS * BP_NS * LBANDS);

#if (VERB)
    printf("[th_s_ent th_v_ent dphi_ent tabind]=[%f %f %f %d]\n",
            th_s_ent, th_v_ent, dphi_ent, tabind);
#endif

    /* Check geometry within table limits */
    if (tabind > ABSIND_NS) {
        printf("\nWARNING: NIWA-IOP: Table limits exceeded.\n");
        return FALSE;
    }

    /* Load data */
    for (i = 0; i < LBANDS; i++) {
        for (y = 0; y < BP_NS; y++) {
            pos = (int) tabind + (i * BP_NS * AP_NS) + (y * AP_NS);
            for (x = 0; x < AP_NS; x++) {
                abst[x][y][i] = absen[pos + x];
            }
        }
    }

#if (VERB)
    printf("LUT loaded: abst=%d\n"), abst;
    printf("Check: endpos %d endpos %d\n", pos + AP_NS, tabind + (AP_NS * BP_NS * LBANDS));

    printf("Bands first:");
    for (x = 0; x < LBANDS; x++)
        printf(" %f", abst[0][0][x]);
    printf("\n");
    printf("Bands last:");
    for (x = 0; x < LBANDS; x++)
        printf(" %f", abst[AP_NS - 1][BP_NS - 1][x]);
    printf("\n");
#endif

    return TRUE;
}

/* ---------------------------------------------------------------------------
 * load_work_tab()
 *  Load the look-up tables from disk
 *
 *  inputs:  none
 *  outputs: global table data initialised from file values
 *  returns: size of 'absen' array
 * ---------------------------------------------------------------------------
 */
static int load_work_tab(void) {
    char *tmp_str;
    char lutdir[FILENAME_MAX];
    char fname[FILENAME_MAX];
    FILE *fp;
    int i, lut_nbands;

    if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
        printf("\nERROR: NIWA-IOP: OCDATAROOT environment variable is not defined.\n");
        exit(EXIT_FAILURE);
    }

    strcpy(lutdir, tmp_str);
    strcat(lutdir, "/common");
#if (VERB) 
    printf("LUT dir = %s\n", lutdir);
#endif

    /* Open header table */
    sprintf(fname, "%s/niwa_iop_h.txt", lutdir);

    printf("NIWA-IOP: Reading LUT header file %s\n", fname);

    fp = fopen(fname, "r");
    if (fp == NULL) {
        printf("\nERROR: NIWA-IOP: Cannot open %s\n", fname);
        exit(EXIT_FAILURE);
    }

    fscanf(fp, "%d", &lut_nbands);
    if (lut_nbands != LBANDS) {
        printf("\nERROR: NIWA-IOP: Band mismatch opening %s\n", fname);
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < LBANDS; i++)
        fscanf(fp, "%f", &band_lev[i]);

    fscanf(fp, "%d", &th_s_n);
    for (i = 0; i < TH_NS; i++)
        fscanf(fp, "%f", &th_s_lev[i]);
    for (i = 0; i < TH_NS; i++)
        th_s_ind[i] = i;
#if (VERB)
    printf("th_s_n %d th_s_lev %f %f\n", th_s_n, th_s_lev[0], th_s_lev[1]);
#endif

    fscanf(fp, "%d", &th_v_n);
    for (i = 0; i < TH_NS; i++)
        fscanf(fp, "%f", &th_v_lev[i]);
    for (i = 0; i < TH_NS; i++)
        th_v_ind[i] = i;
#if (VERB)
    printf("th_v_n %d th_v_lev %f %f\n", th_v_n, th_v_lev[0], th_v_lev[1]);
#endif

    fscanf(fp, "%d", &dphi_n);
    for (i = 0; i < DPHI_NS; i++)
        fscanf(fp, "%f", &dphi_lev[i]);
    for (i = 0; i < DPHI_NS; i++)
        dphi_ind[i] = i;
#if (VERB)
    printf("dphi_n %d dphi_lev %f %f\n", dphi_n, dphi_lev[0], dphi_lev[1]);
#endif

    fscanf(fp, "%d", &ap_n);
    for (i = 0; i < AP_NS; i++)
        fscanf(fp, "%f", &ap_lev[i]);
#if (VERB)
    printf("ap_n %d ap_lev %f %f\n", ap_n, ap_lev[0], ap_lev[1]);
#endif

    fscanf(fp, "%d", &bp_n);
    for (i = 0; i < BP_NS; i++)
        fscanf(fp, "%f", &bp_lev[i]);
#if (VERB)
    printf("bp_n %d bp_lev %f %f\n", bp_n, bp_lev[0], bp_lev[1]);
#endif

    fclose(fp);
#if (VERB)
    printf("NIWA-IOP: Loaded OK %s\n", fname);
#endif
    /* ---End of header--- */

    /*--- Setup the table indexes ---*/
    for (i = 0; i < BP_NS; i++)
        b_ind[i] = i;
    for (i = 0; i < AP_NS; i++)
        a_ind[i] = i;

    /*--- Open data table ---*/
    sprintf(fname, "%s/niwa_iop.dat", lutdir); /* LUT data file */
    printf("NIWA-IOP: Loading LUT data from %s\n", fname);

    fp = fopen(fname, "r");
    if (fp == NULL) {
        printf("\nERROR: NIWA-IOP: Cannot open %s\n", fname);
        exit(EXIT_FAILURE);
    }

    /* Allocate memory for LUT array */
    absind = th_s_n * th_v_n * dphi_n * ap_n * bp_n * LBANDS;
#if (VERB)
    printf("S_F_T: absind= %d\n", absind);
#endif

    if ((absen = (float *) calloc((absind), sizeof (float))) == 0) {
        printf("\nERROR: NIWA-IOP: memory allocation failure, absen\n");
        exit(EXIT_FAILURE);
    }

    fread(absen, sizeof (float), absind, fp);
    fclose(fp);

#if (VERB)
    printf("NIWA-IOP: %s Loaded OK\n", fname);
    printf("absen first %f %f %f %f %f %f %f %f\n", absen[0], absen[1], absen[2], absen[3],
            absen[4], absen[5], absen[6], absen[7]);
    printf("absen last %f %f %f %f %f %f %f %f\n", absen[absind - 8], absen[absind - 7],
            absen[absind - 6], absen[absind - 5], absen[absind - 4], absen[absind - 3],
            absen[absind - 2], absen[absind - 1]);
#endif

    return absind;
}

/* ---------------------------------------------------------------------------
 * mbcucof()
 *
 *  inputs:  
 *  outputs: 
 *  returns: none
 * ---------------------------------------------------------------------------
 */
static void mbcucof(float y[], float y1[], float y2[], float y12[], float d1, float d2, float *cl) {
    static const int wt[16][16] = {
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {-3, 0, 0, 3, 0, 0, 0, 0, -2, 0, 0, -1, 0, 0, 0, 0},
        {2, 0, 0, -2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, -3, 0, 0, 3, 0, 0, 0, 0, -2, 0, 0, -1},
        {0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0, 0, 1, 0, 0, 1},
        {-3, 3, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0},
        {9, -9, 9, -9, 6, 3, -3, -6, 6, -6, -3, 3, 4, 2, 1, 2},
        {-6, 6, -6, 6, -4, -2, 2, 4, -3, 3, 3, -3, -2, -1, -1, -2},
        {2, -2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0},
        {-6, 6, -6, 6, -3, -3, 3, 3, -4, 4, 2, -2, -2, -2, -1, -1},
        {4, -4, 4, -4, 2, 2, -2, -2, 2, -2, -2, 2, 1, 1, 1, 1}
    };

    int k, i;
    float xx, d1d2, x[16];

    d1d2 = d1 * d2;
    for (i = 0; i < 4; i++) {
        x[i] = y[i];
        x[i + 4] = y1[i] * d1;
        x[i + 8] = y2[i] * d2;
        x[i + 12] = y12[i] * d1d2;
    }
    for (i = 0; i < 16; i++) {
        xx = 0.0;
        for (k = 0; k < 16; k++)
            xx += wt[i][k] * x[k];
        cl[i] = xx;
    }
}

/* ---------------------------------------------------------------------------
 * f_abx()
 *
 *  inputs:  
 *  outputs: 
 *  returns: 
 * ---------------------------------------------------------------------------
 */
static float f_abx(float a, float b, int band) {
    int x, y, i;
    float ain, bin;
    float res;
    float tarr[AP_NS][BP_NS];
    int al, bl, ah, bh, xl, xh, yl, yh;
    float t, u, d1, d2, c[4][4], cl[16];
    float yy[4], y1[4], y2[4], y12[4];

    /* don't accept bad a, b inputs */
    if ((a < 0.0) || (b < 0.0))
        return 0.0;

    /* Interpolate loaded LUT values to input a and b */
    ain = interpol(a_ind, ap_lev, a, AP_NS);
    bin = interpol(b_ind, bp_lev, b, BP_NS);

    /* set array indices */
    for (y = 0; y < BP_NS; y++)
        for (x = 0; x < AP_NS; x++)
            tarr[x][y] = abst[x][y][band];

    /* set lower and upper range in table a, b */
    al = (ain > 0) ? (int) floor(ain) : 0;
    bl = (bin > 0) ? (int) floor(bin) : 0;

    bh = bl + 1;
    ah = al + 1;

    /* catch inconsistent upper and lower table a,b values */
    if ((bh > BP_NS - 1) || (ah > AP_NS - 1))
        return 0.0;
    if ((al == ah) || (bl == bh))
        return 0.0;

    /* set indices for smoothed interpolation over 3x3 table coefficient array */
    for (i = 0; i < 4; i++) {
        x = (i == 0 || i == 3) ? al : ah;
        y = (i < 2) ? bl : bh;
        xl = (x == 0) ? x : x - 1;
        xh = (x == AP_NS - 1) ? x : x + 1;
        yl = (y == 0) ? y : y - 1;
        yh = (y == BP_NS - 1) ? y : y + 1;

        yy[i] = tarr[x][y];
        y1[i] = (tarr[xh][y] - tarr[xl][y]) / (ap_lev[xh] - ap_lev[xl]);
        y2[i] = (tarr[x][yh] - tarr[x][yl]) / (bp_lev[yh] - bp_lev[yl]);
        y12[i] = (tarr[xh][yh] - tarr[xh][yl] - tarr[xl][yh] + tarr[xl][yl])
                / ((ap_lev[xh] - ap_lev[xl]) * (bp_lev[yh] - bp_lev[yl]));
    }
    d1 = ap_lev[ah] - ap_lev[al];
    d2 = bp_lev[bh] - bp_lev[bl];

    /* call function to retrieve 3x3 interpolated table look-up value of F */
    mbcucof(yy, y1, y2, y12, d1, d2, cl);
    i = 0;
    for (x = 0; x < 4; x++)
        for (y = 0; y < 4; y++)
            c[x][y] = cl[i++];
    t = (a - ap_lev[al]) / d1;
    u = (b - bp_lev[bl]) / d2;
    res = 0.0;
    for (i = 3; i >= 0; i--)
        res = t * res + ((c[i][3] * u + c[i][2]) * u + c[i][1]) * u + c[i][0];

    return res;
}

/* ---------------------------------------------------------------------------
 * f_ab()
 *  Sensor specific version of f_abx() which corrects for MODIS vs SeaWIFs bands
 *
 *  inputs:  
 *  outputs: 
 *  returns: 
 * ---------------------------------------------------------------------------
 */
static float f_ab(float a, float b, int band, int sensor_id) {
    float res, res1;

    res = f_abx(a, b, band);
    if (band == 3 && sensor_id != SEAWIFS) {
        if (sensor_id == MODIST || sensor_id == MODISA) {
            res1 = f_abx(a, b, 4);
            res += ((530 - 510) / (555 - 510)) * (res1 - res);
        } else res = 0; /* probably better to assert here? */
    }
    return res;
}

/* ---------------------------------------------------------------------------
 * iop_model1() -- calculate a and bb values for one pixel
 *
 *  inputs:  input values for this pixel
 *  outputs: a and bb values (placed in caller's result struct)
 *           (all values will be NO_DATA if not success)
 *  returns: TRUE if success (ie some good data) else FALSE
 *
 * Note: currently only works for SeaWiFS
 * ---------------------------------------------------------------------------
 */
static int iop_model1(float Rrs[], int sensor_id, abs_res_t *result) {
    static const float b_tilde = 0.01756; /* Mobley phase function */
    static const float init_a[] = {0.05, 0.06, 0.05, 0.04, 0.02, 0.03};
    static const int maxIts = 20; /* local limit on no. iterations for F, a and bb */
    static const float tol = 0.001;

    float minEpsa, maxEpsa; /* limits on epsilon_a, Pinkerton et al., 2006 */
    int iter, iw;
    float new_err, old_err;
    float b490, hib490, lob490, bstep, temp1; /* iteration parameters for bb */
    float bestb[2], besta[2];
    float F[VBANDS], a[VBANDS], temp[VBANDS]; /* local variables */
    float b[VBANDS], bn[VBANDS], rho[VBANDS];
    int success;

    /* set min max Epsa according to sensor */
    if (sensor_id == SEAWIFS) {
        minEpsa = 1.22;
        maxEpsa = 1.42;
    }
    else if (sensor_id == MODIST || sensor_id == MODISA) {
        //minEpsa = 1.80;                         /* NOMAD 20th percentile */
        //maxEpsa = 2.34;                         /* NOMAD 80th percentile */
        minEpsa = 1.65; /* NOMAD 5th percentile */
        maxEpsa = 2.77; /* NOMAD 95th percentile */
    }

    /* initialise */
    iter = 0;
    lob490 = 0.001;
    hib490 = 60.0;
    bestb[0] = -1.0;
    bestb[1] = -1.0;
    besta[0] = -1.0;
    besta[1] = -1.0;

    /* convert Rrs to rho */
    for (iw = 0; iw < VBANDS; iw++) {
        rho[iw] = OEL_PI * Rrs[iw];
    }

    /* initialise outputs */
    success = FALSE; /* until we produce some good data */
    for (iw = 0; iw < VBANDS; iw++) {
        result->a[iw] = NO_DATA;
        result->bb[iw] = NO_DATA;
        if (rho[iw] < 0.0)
            rho[iw] = 0.0;

        // b[iw] = 0.24 * ( pow (490.0 / lambda[iw], 2.22)); /* NOMAD data for eps_b spectral slope */
        // b[iw] = (1.62517 - 0.00113 * lambda[iw])/1.00367; */
        b[iw] = (lambda[iw] / lambda[2]) * (0.841 * (lambda[iw] / lambda[2]) - 2.806) + 2.965; /* EBoP */
        a[iw] = init_a[iw];
    }

    /* normalise to b to 489 nm, bn remains fixed from here on */
    for (iw = 0; iw < VBANDS; iw++) {
        bn[iw] = b[iw] / b[2];
    }

#if (VERB)
    printf("--- START iop_model1 ---\n");
    printf("rho=[%f %f %f %f %f %f]\n", rho[0], rho[1], rho[2], rho[3], rho[4], rho[5]);
    printf("b=[%f %f %f %f %f %f]\n", b[0], b[1], b[2], b[3], b[4], b[5]);
    printf("bn=[%f %f %f %f %f %f]\n", bn[0], bn[1], bn[2], bn[3], bn[4], bn[5]);
#endif

    /*** Begin calculations ***/
    b[2] = 0.05; /* start search value of b(490) */
    for (iw = 3; iw < 5; iw++)
        b[iw] = b[2] * bn[iw]; /* start values of b(510) and b(555) */
    for (iw = 2; iw < 5; iw++)
        temp[iw] = b[iw]; /* initialise */
    iter = 0;
    old_err = 100.; /* high starting value */
    new_err = old_err - 1;

    while ((old_err - new_err) > tol && iter < maxIts) { /* stop iterating */
        old_err = new_err;
        new_err = 0.0;
        for (iw = 2; iw < 5; iw++) {
            F[iw] = f_ab(0.0, b[iw], iw, sensor_id); /* 1st estimate of f(lambda, geom, iop) */
            if (F[iw] > 0.0) {
                temp[iw] = (rho[iw] * aw[iw] / F[iw] - bbw[iw]) / b_tilde; /* estimate better b(lambda) */
#if (VERB)
                printf("b[%d]=%f\n", iw, temp[iw]);
#endif
                if (temp[iw] < 0.0)
                    temp[iw] = 0.001; /* if negative, set to low value */
            } else
                temp[iw] = 0.001; /* catch negative F case */

            new_err += fabsf(temp[iw] - b[iw]) / b[iw]; /* calculate proportional change in b */
            b[iw] = temp[iw];

        } /* end loop over 490:555 for this iteration */
        temp1 = (b[2] + b[3] / bn[3] + b[4] / bn[4]) / 3.; /* average value of b2 */
        for (iw = 2; iw < 5; iw++)
            b[iw] = temp1 * bn[iw]; /* reset */
#if (VERB)
        printf("%d getting bF [b2,F2] =[%f %f]  -> %f\n", iter, b[2], F[2], new_err);
#endif
        iter++;
    }

    /* now we have a stable value for F */
#if (VERB)
    printf("got a+F: iter=%d new_err=%f\n", iter, new_err);
    printf("F=[%f %f %f]\n", F[2], F[3], F[4]);
    printf("a=[%f %f %f]\n", a[2], a[3], a[4]);
    printf("b=[%f %f %f]\n", b[2], b[3], b[4]);
    printf("lob490=%f\n", lob490);
#endif

    /* select highest b-value calculated in 490:555 range for conservative limit in eps_a iterations */
    /*for (iw=2; iw<5; iw++) if (b[iw]*bn[iw] > lob490) lob490=b[iw]*bn[iw]; */
    lob490 = b[2];
    for (iw = 3; iw < 5; iw++)
        if (b[iw] / bn[iw] > lob490)
            lob490 = b[iw] / bn[iw];

#if (VERB)
    printf("lob490=%f\n", lob490);
    //printf("lob490=%f  F[2]=%f\n",lob490,F[2]);
#endif

    /*** Take 15 steps in b490, seeking highest and lowest b490 which yield acceptable eps_a values ***/
    /* calculate stepsize */
    bstep = pow((hib490 / lob490), (1.0 / 15.0)) - 1.0;
#if (VERB)
    printf("bstep=%f lob490\n", bstep, lob490);
#endif

    /* initialise starting values */
    b490 = lob490;
    do {
        iter = 0; /* reset iteration counter for this step */

        for (iw = 2; iw < 4; iw++) { /* re-initialise a and b for this step */
            b[iw] = b490 * bn[iw];
            a[iw] = F[2] * (b[iw] * b_tilde + bbw[iw]) / rho[iw] - aw[iw]; /* use a F[2] to init a */
            if (a[iw] < 0.0)
                a[iw] = 0.001;
        }

        /* Iterate for a and F */
        old_err = 100.; /* reset error for this step */
        new_err = old_err - 1;
        while ((old_err - new_err) > tol && iter < maxIts) { /* end iteration for a at this b490 step */
            old_err = new_err;
            new_err = 0.0; /* reset error within loop */

            for (iw = 2; iw < 4; iw++) { /* loop through 490:510 */

                F[iw] = f_ab(a[iw], b[iw], iw, sensor_id); /* get next estimate of F */
                if (F[iw] > 0.0)
                    temp1 = F[iw] * (b[iw] * b_tilde + bbw[iw]) / rho[iw] - aw[iw]; /* use new F to get better a */
                else
                    temp1 = 0.001; /* if F unrealistic, set a=0.001 */
                if (temp1 < 0.0)
                    temp1 = 0.001;
                new_err += fabsf(temp1 - a[iw]) / a[iw]; /* calculate proportional change in estimate */
                a[iw] = temp1;
            }
            iter++;
#if (VERB)
            printf("%d 1getting aF [a2,F2] =[%f %f]  -> %f\n", iter, a[2], F[2], new_err);
            printf("%d 1getting aF [a3,F3] =[%f %f]  -> %f\n", iter, a[3], F[3], new_err);
#endif
        }
#if (VERB)
        printf("finished iters\n");
#endif

        /* Check b and eps_a are still in range (negative values indicate first run) */
        /* on first run through, set lower b limit to lob490 */
        if (bestb[0] < 0.0 || a[2] / a[3] > maxEpsa) {
            bestb[0] = b490;
            besta[0] = a[2];
        }            /* If b and eps_a within range, store upper and lower limits of b and a */
            /* on second run through, set upper b to lob490+bstep */
        else {
            if (bestb[1] < 0.0) {
                bestb[1] = b490;
                besta[1] = a[2];
                //bestb[0]=b490;
                //besta[0]=a[2];
            } else if (a[2] / a[3] > minEpsa) {
                bestb[1] = b490; /* store this upper b-value */
                besta[1] = a[2]; /* store this upper a-value */
            }
        }

#if (VERB)
        printf("besta=[%f %f]\n", besta[0], besta[1]);
        printf("bestb=[%f %f]\n", bestb[0], bestb[1]);
#endif

        b490 += b490 * bstep;
#if (VERB)
        printf(" [a,b,epsa]=[%f %f %f]\n", a[2], b490, a[2] / a[3]);
        printf("---finished while\n");
#endif

    } while (b490 < 1.1 * hib490 && a[2] / a[3] > minEpsa); /* end loop through 15 steps in b490 */

    /* Quality checks */
    if (bestb[0] > 0.0 && bestb[1] > 0.0) {

        /* choose best value for b490 */
        temp1 = pow(besta[0] * besta[1], 0.5);

        if (b490 > 1.1 * hib490)
            b490 = bestb[0];
        else if (bestb[1] / bestb[0] > 20.)
            b490 = pow(bestb[0] * bestb[1], 0.5);
        else if (temp1 > 1.0)
            b490 = bestb[1];
        else if (temp1 < 0.01)
            b490 = bestb[0];
        else
            b490 = 0.5 * (bestb[0] + bestb[1]);

#if (VERB)
        printf("FINAL1 b490=%f\n", b490);
#endif

        /* keep b within range used for look-up tables */
        if (b490 > hib490)
            b490 = hib490;
        else if (b490 < lob490)
            b490 = lob490;
#if (VERB)
        printf("FINAL2 b490=%f\n", b490);
#endif

        /* calculate full spectra for a and bb */
        for (iw = 0; iw < VBANDS; iw++) {
            /* final b is best b490 adjusted to required spectral shape */
            b[iw] = b490 * bn[iw];

            /* initialise a using final best b */
            a[iw] = 0.2 * (b[iw] * b_tilde + bbw[iw]) / rho[iw] - aw[iw];
            if (a[iw] < 0.0)
                a[iw] = 0.001;

            /* iterate for a and F */
            iter = 0; /* reset iteration counter */
            old_err = 100.;
            new_err = new_err - 1.;
            while ((old_err - new_err) > tol && iter < maxIts) {
                old_err = new_err;
                new_err = 0.0;
                F[iw] = f_ab(a[iw], b[iw], iw, sensor_id);
                if (F[iw] > 0.0)
                    temp1 = F[iw] * (b[iw] * b_tilde + bbw[iw]) / rho[iw] - aw[iw];
                else
                    temp1 = 0.001; /* set a to 0.001 if F unrealistic */
                if (temp1 < 0.0)
                    temp1 = 0.001; /* set a to 0.001 if derived a was unrealistic */
                new_err += fabsf(temp1 - a[iw]) / a[iw];
                a[iw] = temp1;
                iter++;
#if (VERB)
                printf("%d 2getting aF [a2,F2] =[%f %f]  -> %f\n", iter, a[2], F[2], new_err);
#endif
            }
#if (VERB)
            printf("finished iters\n");
#endif

            /* transfer to main result */
            result->bb[iw] = b[iw] * b_tilde; /* all 'bb' values considered good at this point */
            if (a[iw] > 0.001001) {
                result->a[iw] = a[iw];
                success = TRUE; /* need just one good 'a' band to consider success initially  */
            }
#if (VERB)
            printf("%d main.result1 [a,bb]=[%f %f]\n", iw, result->a[iw], result->bb[iw]);
#endif
        } /* end wavelength loop for final a,bb */

        /* final quality checks -- both a_490 and a_510 must be good */
        if (result->a[2] < 0.001 || result->a[3] < 0.001) {
            success = FALSE;
        }
    } /* end if upper and lower b490 values were positive */

#if (VERB)
    printf("Returning\n");
    (success) ? printf("success\n") : printf("fail\n");
    printf("result a=%f %f %f %f %f %f\n", result->a[0], result->a[1], result->a[2],
            result->a[3], result->a[4], result->a[5]);
    printf("result bb=%f %f %f %f %f %f\n", result->bb[0], result->bb[1], result->bb[2],
            result->bb[3], result->bb[4], result->bb[5]);
#endif

    return success;
}

/* ---------------------------------------------------------------------------
 * init()
 *  perform first time initialisations and load LUTs etc
 *
 *  inputs:  ptr to the l2 record structure (of first scan line)
 *  outputs: none
 *  returns: none
 * ---------------------------------------------------------------------------
 */
static void init(l2str *l2rec) {
    int iw, ib;

    filehandle *l1file = l2rec->l1rec->l1file;

    printf("\nUsing NIWA / UoP / Moore IOP algorithm.\n");

    if (load_work_tab() != ABSIND_NS) {
        printf("\nERROR: NIWA-IOP: failed to load look-up tables.\n");
        exit(EXIT_FAILURE);
    }
#if (VERB)
    printf("loaded work table OK\n");
#endif

    switch (l1file->sensorID) {
    case SEAWIFS:
        nbands = SEAWIFS_NBANDS;
        band_map = seawifs_map;
        break;
    case MODIST:
    case MODISA:
        nbands = MODIS_NBANDS;
        band_map = modis_map;
        break;
    default:
        printf("\nERROR: NIWA-IOP: The NIWA IOP algorithm requires SeaWiFS or MODIS LAC input data.\n");
        printf("                 [sensor ID = %d]\n", l1file->sensorID);
        exit(EXIT_FAILURE);
    }

    /* make local copy of aw, bbw and lambda for bands we use */
    for (iw = 0; iw < VBANDS; iw++) {
        ib = band_map[iw];
        aw[iw] = l1file->aw[ib];
        bbw[iw] = l1file->bbw[ib];
        lambda[iw] = l1file->fwave[ib];
    }

    /* The look-up tables assume that no brdf correction has been applied 
     * we require either brdf_opt=0 or valid numbers in the brdf array so can reverse out correction */
    if (input->brdf_opt > 0) {
        printf("\nWARNING: NIWA-IOP: brdf correction detected, attempting to reverse (use brdf_opt=0 in future).\n\n");
    }
}

/* --------------------------------------------------------------------------- 
 * niwa_iop() 
 *  this is the main iop function which must be called once for each scan line
 *
 *  inputs:  the L2 data structure for a single scan line
 *  outputs: a & bb values for VBANDS channels in supplied arrays, flags for each pixel
 *           (a & bb arrays expected size = NBANDS x npix)
 *  returns: none
 * ---------------------------------------------------------------------------
 */
void niwa_iop(l2str *l2rec, float niwa_a[], float niwa_bb[], int16 niwa_iopf[]) {
    static int first_time = 1;

    int iw, ib, ip, ibp, badrrs;
    float dphiRad, dphi, theta0Rad, thetaRad, brdf[VBANDS], Rrs[VBANDS];
    abs_res_t result;

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;

    /* do inits first time only */
    if (first_time) {
        init(l2rec);
        first_time = 0;
    }

    /* --- loop across scan --- */
    for (ip = 0; ip < l1rec->npix; ip++) {

        /* start by writing 'no data' value to all output bands 
         * (note this includes bands we do not use) */
        for (ib = 0; ib < nbands; ib++) {
            ibp = ip * l1file->nbands + ib;
            niwa_a[ibp] = NO_DATA;
            niwa_bb[ibp] = NO_DATA;
        }
        niwa_iopf[ip] = 0;

        /* skip if this pixel already masked out */
        if (l1rec->mask[ip]) {
            niwa_iopf[ip] |= IOPF_ISMASKED;
            continue;
        }

        /* take copy of Rrs */
        badrrs = FALSE;
        for (iw = 0; iw < VBANDS; iw++) {

            /* get index to 1-d array l2rec variables */
            ibp = ip * l1file->nbands + band_map[iw];

            /* make local copy of Rrs(lambda) at this pixel */
            Rrs[iw] = l2rec->Rrs[ibp];

            /* if brdf has been applied, we need to back-correct */
            if (input->brdf_opt > 0) {
                brdf[iw] = l2rec->brdf[ibp];
                Rrs[iw] = (brdf[iw] > 0) ? Rrs[iw] / brdf[iw] : -1;
            }

            /* quality check -- all must be >= 0 to procede */
            if (Rrs[iw] < 0)
                badrrs = TRUE;
        }
        if (badrrs) {
            niwa_iopf[ip] |= IOPF_BADRRS;
            continue;
        }
        /* Rrs is now non-brdf-corrected */

        /* set geometry */
        theta0Rad = radians(l1rec->solz[ip]);
        thetaRad = radians(l1rec->senz[ip]);

        /* relative azimuth in 0 - 180.0 (BFranz) */
        dphi = l1rec->sena[ip] - l1rec->sola[ip];
        if (dphi > 180.0)
            dphi -= 360.0;
        if (dphi < -180.0)
            dphi += 360.0;
        dphiRad = radians(dphi);

        if (!setgeom(theta0Rad, thetaRad, dphiRad)) {
            printf("WARNING: NIWA-IOP: Viewing geometry is out of look-up table bounds for NIWA IOP algorithm.\n");
            niwa_iopf[ip] |= IOPF_NIWA_BADGEOM;
            continue;
        }

        /* all looks good -- solve for IOPs at this pixel */
        if (iop_model1(Rrs, l1file->sensorID, &result)) {

            /* success -- copy data to output arrays */
            for (iw = 0; iw < VBANDS; iw++) {
                ibp = ip * l1file->nbands + band_map[iw];
                niwa_a[ibp] = result.a[iw];
                niwa_bb[ibp] = result.bb[iw];
            }
        } else {
            /* record that we attempted to calculate IOPs but failed */
            niwa_iopf[ip] |= IOPF_FAILED;

            //### debug for Matt
            //printf("WARNING: NIWA-IOP: iop_model failed at line %d, pixel %d:\n", l2rec->iscan, ip);
            //for (iw = 0; iw < VBANDS; ++iw) {
            //    printf("band %d: lambda = %f, Rrs = %f, aw = %f, bbw = %f, a = %f, bb = %f\n", 
            //           iw, lambda[iw], Rrs[iw], aw[iw], bbw[iw], result.a[iw], result.bb[iw]);
            //}
        }
    }
    /* report fail for all pixels where any flags were set */
    for (ip = 0; ip < l1rec->npix; ip++)
        if (niwa_iopf[ip] != 0)
            l1rec->flags[ip] |= PRODFAIL;
}
