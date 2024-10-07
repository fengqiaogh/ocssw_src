/* =========================================================================
 *
 * Calculates aph675 and ag400 algebraically from Rrs model equations.
 * chl is then calculated from aph675.
 *
 * Original Authors:
 *   Kendall Carder, F. R. Chen, and Steve Hawes
 *   University of South Florida
 *   Department of Marine Science
 *   140 7th Avenue South
 *   St. Petersburg, FL 33701
 *   kcarder@monty.marine.usf.edu
 *
 * References:
 * - Carder et al., Semi-Analytic MODIS Algorithms for Chlorophyll a and
 *   Absorption with Bio-Optical Domains Based on Nitrate-Depletion
 *   Temperatures, JGR MODIS special issue, 1998.
 * - Carder et al., Reflectance model for quantifying chlorophyll _a_
 *   in the presence of productivity degradation products, JGR, 96(C11),
 *   20599-20611, 1991.
 * - Lee et al., Model for the interpretation of hyperspectral remote-
 *   sensing reflectance, Appl. Opt, 33(24), 5721, 1994.
 *
 * Modifications:
 *   18 May 2004, B. Franz
 *     Generalized for use in MSl12.  Pass in sensor wavelengths.  Compute 
 *     aw and bbw at input wavelengths. Add use of NDT for packaging switch
 *     based on modis_chl-1.2.c code. Code simplification. Enhanced flag
 *     tracking.  General reorganization for efficiency.
 *
 * ========================================================================= */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "l12_proto.h"

#define N_PRM_HEADERS  35  /* leading comment lines in PARAM_FILE            */
#define NX             32  /* elements minus one in aph675 array.            */
#define N_ITER          5  /* iterations in bisection loop.  2^(N_ITER) = NX */
#define APH_MIN    0.0001  /* minimum aph675 value                           */
#define CARDER_FAIL    -1  /* value of products returned if failure          */

/* flag bit settings */
#define CFLAG_FAIL    0x0001           /* failure (bad inputs)               */
#define CFLAG_DEFAULT 0x0002           /* using default (empirical) algorithm*/
#define CFLAG_BLEND   0x0004           /* blending with default algorithm    */

#define CFLAG_UNPKG   0x0010           /* unpackaged case                    */
#define CFLAG_PKG     0x0020           /* packaged case                      */
#define CFLAG_HIPKG   0x0040           /* highly packaged case               */
#define CFLAG_GLOBAL  0x0080           /* global case                        */

#define CFLAG_LO412   0x0100           /* low input 412 reflectance          */
#define CFLAG_LO555   0x0200           /* low input 55x reflectance          */
#define CFLAG_CHLINC  0x0400           /* large inconsistency with default   */

static int pk; /* packaging parameter                */
static int band1, band2, band3, band4; /* band indices to make Rrs ratios    */
static int lam[6]; /* seawifs wavelengths                */
static int bb_denom; /* toggle bb into denom of Rrs eqn    */
static double bbw[6]; /* water backscattering coefs.        */
static double aw[6]; /* water absorption coefs             */
static double a0[6], a1[6], a2[6], a3[6]; /* aph model prms                  */
static double ga0[6], ga3[6]; /* aph model prms global              */
static double upa0[6], upa3[6]; /* aph model prms unpackaged          */
static double pa0[6], pa3[6]; /* aph model prms packaged            */
static double hpa0[6], hpa1[6], hpa2[6], hpa3[6]; /* aph model prms hipackaged*/
static double s; /* ag exponential slope prm           */
static double x0, x1; /* regr. coefs for x                  */
static double y_0, y_1; /* regr. coefs for y                  */
static double c0, c1, c2, c3; /* coefs for 490:555 chl alg          */
static double gc0, gc1, gc2, gc3; /* coefs for 490:555 chl alg global   */
static double upc0, upc1, upc2, upc3; /* coefs for 490:555 chl alg unpackaged*/
static double pc0, pc1, pc2, pc3; /* coefs for 490:555 chl alg packaged */
static double hpc0, hpc1, hpc2, hpc3; /* coefs for 490:555 chl alg hipackaged*/
static double p0, p1, p2; /* coefs for chl=fn(aph675)           */
static double gp0, gp1, gp2; /* coefs for chl=fn(aph675) global    */
static double upp0, upp1, upp2; /* coefs for chl=fn(aph675) unpackaged*/
static double pp0, pp1, pp2; /* coefs for chl=fn(aph675) packaged  */
static double hpp0, hpp1, hpp2; /* coefs for chl=fn(aph675) hipackaged*/
static double d1, d2, d3; /* filter table parameters            */
static double low_412_thresh, low_555_thresh, chl_inconsistent_thresh;
static double aph_hi, aph_lo; /* boundaries of aph675 array         */
static double delta;
static double aphb1, aphb2, aphb3, aphb4, fmid;
static double tx[NX + 1];
static double arr0, arr1;

/* function declarations */

static void get_params(void);
static double get_aph(int band, double aph675);


/* ========================================================================= */

/* ========================================================================= */
int carder_model(int32_t nbands, double *rrs, double *rrswl, float sst, float ndt,
        double *atot, double *aph,
        double *adg, double *bb, double *bbp,
        double *aphmPtr, double *adgmPtr, double *chlPtr,
        int16 *flags) {

    double br15, br25, br35, exponent, laph;
    double adg_def, aph_def, chl_def;
    double abr15, abr25, abr35;
    double x, y, bb1, bb2, bb3, bb4; /* backscattering prms */
    double r12, r34; /* Rrs/bb ratios */
    double g12, g34, f0, f1, f2, f3, f4, wph;
    double funk_lo, funk_hi, flo, fhi;
    double wt;
    double chl_mod, aph_mod, adg_mod; /* modeled chl,aph675,and ag400  */
    double tchl_mod, taph_mod, tadg_mod; /* modeled chl,aph675,and ag400 (temp) */

    double weit; /* package switching */
    double uptemp, gtemp, ptemp;
    double sw0lo, sw1lo, sw2lo, sw3lo;
    double dtsw1, dtsw2, dtsw3;
    int it, pktran[2];

    int i, xlo, xmid, xhi;
    static double cpa1[6], cpa2[6]; /* copies of a1 and a2 */
    static int firstCall = 1;

    if (firstCall) {

        firstCall = 0;

        /* read model parameters */
        get_params();

        /* copy the a1 and a2 arrays */
        for (i = 0; i < 6; i++) cpa1[i] = a1[i];
        for (i = 0; i < 6; i++) cpa2[i] = a2[i];

        /* create aph675 array, 0.0001-aph_hi (1/m) */
        arr0 = log10(aph_lo);
        arr1 = log10(aph_hi) - arr0;
        for (i = 0; i <= NX; i++)
            tx[i] = pow(10.0, arr0 + arr1 * (double) i / (double) NX);

        /* over-ride parameter inputs with sensor-band-specific aw, bbw  */
        /* using 10-nm square bandpass average to be consistent with Rrs */
        for (i = 0; i < 6; i++) {
            if (rrswl[i] > 0.0) {
                aw [i] = (double) aw_spectra(rrswl[i], BANDW);
                bbw[i] = (double) bbw_spectra(rrswl[i], BANDW);
            } else {
                aw [i] = (double) 1.0;
                bbw[i] = (double) 1.0;
            }
        }
    }

    /* initialize */
    for (i = 0; i < nbands; i++) {
        bb [i] = CARDER_FAIL;
        atot[i] = CARDER_FAIL;
        aph [i] = CARDER_FAIL;
        adg [i] = CARDER_FAIL;
    }
    *chlPtr = CARDER_FAIL;
    *adgmPtr = CARDER_FAIL;
    *aphmPtr = CARDER_FAIL;
    *flags = 0x0000;

    /* if we have any negative reflectances, just bail */
    if (rrs[0] <= 0. || rrs[1] <= 0. || rrs[2] <= 0. || rrs[4] <= 0.) {
        *flags |= CFLAG_FAIL;
        return 1;
    }

    /* set some quality flags */
    if (rrs[0] < low_412_thresh) {
        *flags |= CFLAG_LO412;
    }
    if (rrs[4] < low_555_thresh) {
        *flags |= CFLAG_LO555;
    }

    /* set up the limits based on NDT */
    uptemp = ndt + 2.;
    gtemp = ndt + 0.8;
    ptemp = ndt - 1.;

    sw0lo = uptemp + 1.;
    sw1lo = gtemp + (uptemp - gtemp) / 2.;
    sw2lo = ptemp + (gtemp - ptemp) / 2.;
    sw3lo = ptemp - 1.;

    dtsw1 = sw0lo - sw1lo;
    dtsw2 = sw1lo - sw2lo;
    dtsw3 = sw2lo - sw3lo;

    /* sw = 0 */
    if (sst >= sw0lo) {
        weit = 1.;
        pktran[0] = 1;
        pktran[1] = 1;
    }

    /* sw = 1 */
    if (sst >= sw1lo && sst < sw0lo) {
        weit = (sst - sw1lo) / dtsw1;
        pktran[0] = 0;
        pktran[1] = 1;
    }

    /* sw = 2 */
    if (sst >= sw2lo && sst < sw1lo) {
        weit = (sst - sw2lo) / dtsw2;
        pktran[0] = 2;
        pktran[1] = 0;
    }

    /* sw = 3 */
    if (sst >= sw3lo && sst < sw2lo) {
        weit = (sst - sw3lo) / dtsw3;
        pktran[0] = 3;
        pktran[1] = 2;
    }

    /* sw = 4 */
    if (sst < sw3lo) {
        weit = 1.;
        pktran[0] = 3;
        pktran[1] = 3;
    }


    /* calculate some bandratios */
    br15 = rrs[0] / rrs[4];
    br25 = rrs[1] / rrs[4];
    br35 = rrs[2] / rrs[4];
    abr15 = log10(br15);
    abr25 = log10(br25);
    abr35 = log10(br35);

    /* calculate default ag400 and aph675 */
    exponent = MAX(MIN(-1.147 - 1.963 * abr15 - 1.01 * abr15 * abr15 + 0.856 * abr25 + 1.702 * abr25*abr25, 10.0), -10.0);
    adg_def = 1.5 * pow(10.0, exponent);
    exponent = MAX(MIN(-0.919 + 1.037 * abr25 - 0.407 * abr25 * abr25 - 3.531 * abr35 + 1.579 * abr35*abr35, 10.0), -10.0);
    aph_def = MAX((pow(10.0, exponent) - 0.008) / 3.05, APH_MIN);

    /* calculate bb */
    if (rrs[4] > 0.001)
        x = MAX(x0 + x1 * rrs[4], 0.0);
    else
        x = MAX(x0 + x1 * 0.001, 0.0);

    y = MAX(y_0 + y_1 * rrs[1] / rrs[2], 0.0);

    for (i = 0; i < 6; i++) {
        bbp[i] = x * pow(551. / rrswl[i], y);
        bb [i] = bbw[i] + bbp[i];
    }

    /* calc. f0, f1, f2, f3, f4; coefs. for funk() */
    bb1 = bb[band1];
    bb2 = bb[band2];
    bb3 = bb[band3];
    bb4 = bb[band4];
    r12 = (rrs[band1] / bb1) / (rrs[band2] / bb2);
    r34 = (rrs[band3] / bb3) / (rrs[band4] / bb4);
    g12 = r12 * exp(-s * (rrswl[band1] - 400.)) - exp(-s * (rrswl[band2] - 400.));
    g34 = r34 * exp(-s * (rrswl[band3] - 400.)) - exp(-s * (rrswl[band4] - 400.));
    if (bb_denom) {
        f0 = g12 * (aw[band4] + bb4 - r34 * (aw[band3] + bb3))
                - g34 * (aw[band2] + bb2 - r12 * (aw[band1] + bb1));
    } else {
        f0 = g12 * (aw[band4] - r34 * aw[band3])
                - g34 * (aw[band2] - r12 * aw[band1]);
    }
    f1 = g34*r12;
    f2 = -g34;
    f3 = -g12*r34;
    f4 = g12;


    /* Loop through twice, to blend packaging effects */
    for (it = 0; it <= 1; it++) {

        pk = pktran[it];

        switch (pk) {
        case 1: /* unpackaged */
            a0[0] = upa0[0];
            a0[1] = upa0[1];
            a0[2] = upa0[2];
            a0[3] = upa0[3];
            a0[4] = upa0[4];
            a0[5] = upa0[5];
            a1[0] = cpa1[0];
            a1[1] = cpa1[1];
            a1[2] = cpa1[2];
            a1[3] = cpa1[3];
            a1[4] = cpa1[4];
            a1[5] = cpa1[5];
            a2[0] = cpa2[0];
            a2[1] = cpa2[1];
            a2[2] = cpa2[2];
            a2[3] = cpa2[3];
            a2[4] = cpa2[4];
            a2[5] = cpa2[5];
            a3[0] = upa3[0];
            a3[1] = upa3[1];
            a3[2] = upa3[2];
            a3[3] = upa3[3];
            a3[4] = upa3[4];
            a3[5] = upa3[5];
            c0 = upc0;
            c1 = upc1;
            c2 = upc2;
            c3 = upc3;
            p0 = upp0;
            p1 = upp1;
            p2 = upp2;
            *flags |= CFLAG_UNPKG;
            break;
        case 2: /* packaged */
            a0[0] = pa0[0];
            a0[1] = pa0[1];
            a0[2] = pa0[2];
            a0[3] = pa0[3];
            a0[4] = pa0[4];
            a0[5] = pa0[5];
            a1[0] = cpa1[0];
            a1[1] = cpa1[1];
            a1[2] = cpa1[2];
            a1[3] = cpa1[3];
            a1[4] = cpa1[4];
            a1[5] = cpa1[5];
            a2[0] = cpa2[0];
            a2[1] = cpa2[1];
            a2[2] = cpa2[2];
            a2[3] = cpa2[3];
            a2[4] = cpa2[4];
            a2[5] = cpa2[5];
            a3[0] = pa3[0];
            a3[1] = pa3[1];
            a3[2] = pa3[2];
            a3[3] = pa3[3];
            a3[4] = pa3[4];
            a3[5] = pa3[5];
            c0 = pc0;
            c1 = pc1;
            c2 = pc2;
            c3 = pc3;
            p0 = pp0;
            p1 = pp1;
            p2 = pp2;
            *flags |= CFLAG_PKG;
            break;
        case 3: /*hipackaged */
            a0[0] = hpa0[0];
            a0[1] = hpa0[1];
            a0[2] = hpa0[2];
            a0[3] = hpa0[3];
            a0[4] = hpa0[4];
            a0[5] = hpa0[5];
            a1[0] = hpa1[0];
            a1[1] = hpa1[1];
            a1[2] = hpa1[2];
            a1[3] = hpa1[3];
            a1[4] = hpa1[4];
            a1[5] = hpa1[5];
            a2[0] = hpa2[0];
            a2[1] = hpa2[1];
            a2[2] = hpa2[2];
            a2[3] = hpa2[3];
            a2[4] = hpa2[4];
            a2[5] = hpa2[5];
            a3[0] = hpa3[0];
            a3[1] = hpa3[1];
            a3[2] = hpa3[2];
            a3[3] = hpa3[3];
            a3[4] = hpa3[4];
            a3[5] = hpa3[5];
            c0 = hpc0;
            c1 = hpc1;
            c2 = hpc2;
            c3 = hpc3;
            p0 = hpp0;
            p1 = hpp1;
            p2 = hpp2;
            *flags |= CFLAG_HIPKG;
            break;
        default: /* global */
            a0[0] = ga0[0];
            a0[1] = ga0[1];
            a0[2] = ga0[2];
            a0[3] = ga0[3];
            a0[4] = ga0[4];
            a0[5] = ga0[5];
            a1[0] = cpa1[0];
            a1[1] = cpa1[1];
            a1[2] = cpa1[2];
            a1[3] = cpa1[3];
            a1[4] = cpa1[4];
            a1[5] = cpa1[5];
            a2[0] = cpa2[0];
            a2[1] = cpa2[1];
            a2[2] = cpa2[2];
            a2[3] = cpa2[3];
            a2[4] = cpa2[4];
            a2[5] = cpa2[5];
            a3[0] = ga3[0];
            a3[1] = ga3[1];
            a3[2] = ga3[2];
            a3[3] = ga3[3];
            a3[4] = ga3[4];
            a3[5] = ga3[5];
            c0 = gc0;
            c1 = gc1;
            c2 = gc2;
            c3 = gc3;
            p0 = gp0;
            p1 = gp1;
            p2 = gp2;
            *flags |= CFLAG_GLOBAL;
            break;
        }


        /* calculate default chl, note coeffs depend on packaging */
        exponent = MAX(MIN(c0 + c1 * abr35 + c2 * abr35 * abr35 + c3 * abr35 * abr35*abr35, 10.0), -10.0);
        chl_def = pow(10.0, exponent);


        /* solve model */
        /* ----------- */

        aphb1 = get_aph(band1, tx[NX]);
        aphb2 = get_aph(band2, tx[NX]);
        aphb3 = get_aph(band3, tx[NX]);
        aphb4 = get_aph(band4, tx[NX]);
        funk_hi = f0 + f1 * aphb1 + f2 * aphb2 + f3 * aphb3 + f4*aphb4;

        aphb1 = get_aph(band1, tx[0]);
        aphb2 = get_aph(band2, tx[0]);
        aphb3 = get_aph(band3, tx[0]);
        aphb4 = get_aph(band4, tx[0]);
        funk_lo = f0 + f1 * aphb1 + f2 * aphb2 + f3 * aphb3 + f4*aphb4;

        if (funk_lo > 0.0 && funk_hi < 0.0) {

            /* solve for aph675 via bisection */

            xlo = 0;
            xhi = NX;

            for (i = 1; i <= N_ITER; i++) {

                xmid = (xlo + xhi) / 2;
                aphb1 = get_aph(band1, tx[xmid]);
                aphb2 = get_aph(band2, tx[xmid]);
                aphb3 = get_aph(band3, tx[xmid]);
                aphb4 = get_aph(band4, tx[xmid]);
                fmid = f0 + f1 * aphb1 + f2 * aphb2 + f3 * aphb3 + f4*aphb4;
                if (fmid < 0.0)
                    xhi = xmid;
                else
                    xlo = xmid;
            }

            aphb1 = get_aph(band1, tx[xlo]);
            aphb2 = get_aph(band2, tx[xlo]);
            aphb3 = get_aph(band3, tx[xlo]);
            aphb4 = get_aph(band4, tx[xlo]);
            flo = f0 + f1 * aphb1 + f2 * aphb2 + f3 * aphb3 + f4*aphb4;

            aphb1 = get_aph(band1, tx[xhi]);
            aphb2 = get_aph(band2, tx[xhi]);
            aphb3 = get_aph(band3, tx[xhi]);
            aphb4 = get_aph(band4, tx[xhi]);
            fhi = f0 + f1 * aphb1 + f2 * aphb2 + f3 * aphb3 + f4*aphb4;

            /* aph675; linear interp across funk to find zero */
            taph_mod = tx[xlo] + (tx[xhi] - tx[xlo]) * flo / (flo - fhi);

            /* ag400 = fn(aph675) */
            if (bb_denom) {
                wph = aw[band4] + get_aph(band4, taph_mod) + bb4
                        - r34 * (aw[band3] + get_aph(band3, taph_mod) + bb3);
            } else {
                wph = aw[band4] + get_aph(band4, taph_mod)
                        - r34 * (aw[band3] + get_aph(band3, taph_mod));
            }
            tadg_mod = wph / g34;

            /* chlorophyll */
            laph = log10(taph_mod);
            exponent = p0 + p1 * laph + p2 * laph * laph;
            tchl_mod = pow(10.0, exponent);

            if (fabs(log10(tchl_mod / chl_def)) > chl_inconsistent_thresh) {
                *flags |= CFLAG_CHLINC;
            }

            /* if aph_mod > aph_hi/2, blend with default alg. */
            if (taph_mod > aph_hi / 2.) {
                wt = (tx[NX] - taph_mod) / (tx[NX] - aph_hi / 2.);
                tchl_mod = wt * tchl_mod + (1.0 - wt) * chl_def;
                tadg_mod = wt * tadg_mod + (1.0 - wt) * adg_def;
                taph_mod = wt * taph_mod + (1.0 - wt) * aph_def;
                *flags |= CFLAG_BLEND;
            }

        } else {

            /* if funk doesn't bracket zero, return defaults */

            tchl_mod = chl_def;
            tadg_mod = adg_def;
            taph_mod = aph_def;
            *flags |= CFLAG_DEFAULT;
        }

        /* weight results */
        if (it == 0) {
            chl_mod = tchl_mod;
            adg_mod = tadg_mod;
            aph_mod = taph_mod;
            for (i = 0; i < 6; i++)
                aph[i] = get_aph(i, aph_mod);
        } else {
            chl_mod = chl_mod * (1. - weit) + tchl_mod * weit;
            adg_mod = adg_mod * (1. - weit) + tadg_mod * weit;
            aph_mod = aph_mod * (1. - weit) + taph_mod * weit;
            for (i = 0; i < 6; i++)
                aph[i] = aph[i]*(1. - weit) + get_aph(i, aph_mod) * weit;
        }
    }

    /* Store results for output */
    for (i = 0; i < 6; i++) {
        adg [i] = adg_mod * exp(-s * (rrswl[i] - 400.0));
        atot[i] = aw[i] + aph[i] + adg[i];
    }

    *aphmPtr = aph_mod; /* model value aph(675) */
    *adgmPtr = adg_mod; /* model value adg(400) */
    *chlPtr = chl_mod;

    return 0;

}

/********************  FUNCTION get_params  ******************************
 *
 *  Reads in the the seawifs wavebands and Rrs model parameters to use
 *  in the algorithm.
 *
 *************************************************************************/

static void get_params(void) {
    int i;
    char dummy[100], full_file[2048];
    FILE *fin;
    char *tmp_str;

    if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
        printf("OCDATAROOT environment variable is not defined.\n");
        exit(1);
    }

    strcpy(full_file, tmp_str);
    strcat(full_file, "/common/carder.par");
    fin = fopen(full_file, "r");
    if (!fin) {
        printf("-E- %s line %d: error opening (%s) file", __FILE__, __LINE__, full_file);
        exit(1);
    }

    printf("\nLoading Carder parameters from file %s\n", full_file);

    for (i = 1; i <= N_PRM_HEADERS; i++)
        fgets(dummy, 100 - 1, fin);
    fscanf(fin, "%d\n", &pk);
    fscanf(fin, "%lf %lf %lf %lf\n", &d1, &d2, &d3, &delta); /* filter table param*/
    fscanf(fin, "%d %d %d %d\n", &band1, &band2, &band3, &band4);
    for (i = 0; i < 6; i++)
        fscanf(fin, "%d %lf %lf %lf %lf %lf %lf\n",
            &lam[i], &bbw[i], &aw[i], &ga0[i], &a1[i], &a2[i], &ga3[i]);
    fscanf(fin, "%lf\n", &s); /* ag exp. slope        */
    fscanf(fin, "%lf %lf\n", &x0, &x1); /* regr. coefs for x    */
    fscanf(fin, "%lf %lf\n", &y_0, &y_1); /* regr. coefs for y    */
    fscanf(fin, "%lf %lf %lf %lf\n", &gc0, &gc1, &gc2, &gc3); /* coefs for r35 chl*/
    fscanf(fin, "%lf %lf %lf\n", &gp0, &gp1, &gp2); /* coefs for chl-aph675 */
    fscanf(fin, "%lf %lf %lf\n", &low_412_thresh, &low_555_thresh,
            &chl_inconsistent_thresh);
    fscanf(fin, "%lf %lf\n", &aph_lo, &aph_hi);
    fscanf(fin, "%d\n", &bb_denom);
    fscanf(fin, "%lf %lf %lf %lf %lf %lf\n", &upa0[0], &upa0[1],
            &upa0[2], &upa0[3], &upa0[4], &upa0[5]);
    fscanf(fin, "%lf %lf %lf %lf %lf %lf\n", &upa3[0], &upa3[1],
            &upa3[2], &upa3[3], &upa3[4], &upa3[5]);
    fscanf(fin, "%lf %lf %lf %lf %lf %lf\n", &pa0[0], &pa0[1],
            &pa0[2], &pa0[3], &pa0[4], &pa0[5]);
    fscanf(fin, "%lf %lf %lf %lf %lf %lf\n", &pa3[0], &pa3[1],
            &pa3[2], &pa3[3], &pa3[4], &pa3[5]);
    fscanf(fin, "%lf %lf %lf %lf %lf %lf\n", &hpa0[0], &hpa0[1],
            &hpa0[2], &hpa0[3], &hpa0[4], &hpa0[5]);
    fscanf(fin, "%lf %lf %lf %lf %lf %lf\n", &hpa1[0], &hpa1[1],
            &hpa1[2], &hpa1[3], &hpa1[4], &hpa1[5]);
    fscanf(fin, "%lf %lf %lf %lf %lf %lf\n", &hpa2[0], &hpa2[1],
            &hpa2[2], &hpa2[3], &hpa2[4], &hpa2[5]);
    fscanf(fin, "%lf %lf %lf %lf %lf %lf\n", &hpa3[0], &hpa3[1],
            &hpa3[2], &hpa3[3], &hpa3[4], &hpa3[5]);
    fscanf(fin, "%lf %lf %lf %lf\n", &upc0, &upc1, &upc2, &upc3);
    fscanf(fin, "%lf %lf %lf %lf\n", &pc0, &pc1, &pc2, &pc3);
    fscanf(fin, "%lf %lf %lf %lf\n", &hpc0, &hpc1, &hpc2, &hpc3);
    fscanf(fin, "%lf %lf %lf\n", &upp0, &upp1, &upp2);
    fscanf(fin, "%lf %lf %lf\n", &pp0, &pp1, &pp2);
    fscanf(fin, "%lf %lf %lf\n", &hpp0, &hpp1, &hpp2);
    fclose(fin);
}

/********************  FUNCTION aph  *************************************
 *
 *  Returns a_phi at a given waveband as a function of aph675.
 *
 *************************************************************************/

static double get_aph(int band, double aph675) {
    if (aph675 == -1.0) {
        return (-1.0);
    } else {
        return a0[band] * exp(a1[band] * tanh(a2[band] * log(aph675 / a3[band]))) * aph675;
    }
}




/* ----------------------------------------------------------------------- */
/*  get_ndt() - read and interpolate NDT map                               */

/* ----------------------------------------------------------------------- */
static float get_ndt(float lon, float lat) {

#define NDTNX 2048
#define NDTNY 1024

    static int firstCall = 1;

    typedef float map_t[NDTNX];
    static map_t* map;
    static float slope;
    static float offset;

    int i, j;

    if (firstCall) {

        char *tmp_str;
        int32 sd_id;
        int32 sds_id;
        int32 status;
        int32 sds_index;
        int32 rank;
        int32 nt;
        int32 dims[H4_MAX_VAR_DIMS];
        int32 nattrs;
        int32 start[2];
        int32 edges[2];
        char name[H4_MAX_NC_NAME];
        char sdsname[] = "ndt";
        char file[FILENAME_MAX];

        firstCall = 0;

        // allocate map array
        map = (map_t*) allocateMemory(NDTNY * sizeof (map_t), "ndt map");


        if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
            printf("OCDATAROOT environment variable is not defined.\n");
            exit(1);
        }
        strcpy(file, tmp_str);
        strcat(file, "/common/ndt.hdf");

        printf("\nLoading nitrate depletion temperature map from %s\n", file);

        /* Open the file and initiate the SD interface */
        sd_id = SDstart(file, DFACC_RDONLY);
        if (sd_id == -1) {
            printf("-E- %s:  Error opening file %s.\n",
                    __FILE__, file);
            exit(1);
        }

        /* Get the SDS index */
        sds_index = SDnametoindex(sd_id, sdsname);
        if (sds_index == -1) {
            printf("-E- %s:  Error seeking %s SDS from %s.\n",
                    __FILE__, sdsname, file);
            exit(1);
        }

        /* Select the SDS */
        sds_id = SDselect(sd_id, sds_index);

        /* Verify the characteristics of the array */
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS info for %s from %s.\n",
                    __FILE__, sdsname, file);
            exit(1);
        }
        if (dims[0] != NDTNY || dims[1] != NDTNX) {
            printf("-E- %s:  Dimension mis-match on %s array from %s.\n",
                    __FILE__, sdsname, file);
            printf("  Expecting %d x %d\n", NDTNX, NDTNY);
            printf("  Reading   %d x %d\n", dims[1], dims[0]);
            exit(1);
        }

        start[0] = 0; /* row offset */
        start[1] = 0; /* col offset */
        edges[0] = NDTNY; /* row count  */
        edges[1] = NDTNX; /* col count  */

        status = SDreaddata(sds_id, start, NULL, edges, (VOIDP) map);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n",
                    __FILE__, sdsname, file);
            exit(1);
        }

        /* Read the scaling info */
        status = SDreadattr(sds_id, SDfindattr(sds_id, "slope"), (VOIDP) & slope);
        if (status != 0) {
            printf("-E- %s:  Error reading slope attribute for SDS %s from %s.\n",
                    __FILE__, sdsname, file);
            exit(1);
        }
        status = SDreadattr(sds_id, SDfindattr(sds_id, "intercept"), (VOIDP) & offset);
        if (status != 0) {
            printf("-E- %s:  Error reading intercept attribute for SDS %s from %s.\n",
                    __FILE__, sdsname, file);
            exit(1);
        }

        /* Terminate access to the array */
        status = SDendaccess(sds_id);

        /* Terminate access to the SD interface and close the file */
        status = SDend(sd_id);
    }

    /* No reason for perfection here */

    i = MIN(MAX((int) ((lon + 180.) / 360.0 * NDTNX), 0), NDTNX - 1);
    j = MIN(MAX((int) ((lat + 90.) / 180.0 * NDTNY), 0), NDTNY - 1);

    return (map[j][i] * slope + offset);
}


/* ----------------------------------------------------------------------- */
/*  MSl12 interface to Carder model (from carder_iop.c)                    */
/* ----------------------------------------------------------------------- */

#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"

/* These hold 2-D arrays of BIP for each of the products calculated by carder. */

static double *atot;
static double *aph, *aphm;
static double *adg, *adgm;
static double *bbp, *bb;
static double *chl;
static int16 *flagm;
static int CarderRecNum = -99;

/*  Allocate private arrays for a single scan line */

void alloc_carder(int32_t npix, int32_t nbands) {
    atot = (double*) malloc(npix * nbands * sizeof (double));
    aph = (double*) malloc(npix * nbands * sizeof (double));
    adg = (double*) malloc(npix * nbands * sizeof (double));
    bb = (double*) malloc(npix * nbands * sizeof (double));
    bbp = (double*) malloc(npix * nbands * sizeof (double));
    chl = (double*) malloc(npix * sizeof (double));
    aphm = (double*) malloc(npix * sizeof (double));
    adgm = (double*) malloc(npix * sizeof (double));
    flagm = (int16*) malloc(npix * sizeof (int16));
}

/* Determine if we have already run carder's algorithm for this scan line */

int carder_ran(int recnum) {

    if (recnum == CarderRecNum)
        return 1;
    else
        return 0;

}

/* Set-up and run the Carder model for each pixel in this scan */

void carder(l2str *l2rec) {
    static double *rrswl;
    static int *indx;
    static int firstCall = 1;

    double *rrs;
    double *aph1;
    double *adg1;
    double *atot1;

    double *bb1;
    double *bbp1;
    double aph_mod; /* aph(675) */
    double adg_mod; /* adg(400) */
    double chl_mod;
    int16 flags;
    int ip, ic;
    int ib;
    int32_t ipb;
    float ndt;
    float sst;
    int status, nbands;

    l1str *l1rec = l2rec->l1rec;
    nbands = l1rec->l1file->nbands;


    if (firstCall) {
        firstCall = 0;
        if ((rrswl = (double *) calloc(nbands, sizeof (double))) == NULL) {
            printf("-E- : Error allocating memory to rrswl\n");
            exit(FATAL_ERROR);
        }
        if ((indx = (int *) calloc(nbands, sizeof (int))) == NULL) {
            printf("-E- : Error allocating memory to indx\n");
            exit(FATAL_ERROR);
        }

        alloc_carder(l1rec->npix, nbands);
        indx[0] = bindex_get(412);
        indx[1] = bindex_get(443);
        indx[2] = bindex_get(490);
        indx[3] = bindex_get(510);
        if (indx[3] < 0) indx[3] = bindex_get(520);
        if (indx[3] < 0) indx[3] = bindex_get(530);
        indx[4] = bindex_get(551);
        indx[5] = bindex_get(670);
        for (ib = 0; ib < 6; ib++) {
            rrswl[ib] = (double) l1rec->l1file->fwave[indx[ib]];
            printf("Carder wavelength %d is %f\n", ib, rrswl[ib]);
        }
    }

    if ((rrs = (double *) calloc(nbands, sizeof (double))) == NULL) {
        printf("-E- : Error allocating memory to rrs\n");
        exit(FATAL_ERROR);
    }
    if ((aph1 = (double *) calloc(nbands, sizeof (double))) == NULL) {
        printf("-E- : Error allocating memory to aph1\n");
        exit(FATAL_ERROR);
    }
    if ((adg1 = (double *) calloc(nbands, sizeof (double))) == NULL) {
        printf("-E- : Error allocating memory to adg1\n");
        exit(FATAL_ERROR);
    }
    if ((atot1 = (double *) calloc(nbands, sizeof (double))) == NULL) {
        printf("-E- : Error allocating memory to atot1\n");
        exit(FATAL_ERROR);
    }
    if ((bb1 = (double *) calloc(nbands, sizeof (double))) == NULL) {
        printf("-E- : Error allocating memory to bb1\n");
        exit(FATAL_ERROR);
    }
    if ((bbp1 = (double *) calloc(nbands, sizeof (double))) == NULL) {
        printf("-E- : Error allocating memory to bbp1\n");
        exit(FATAL_ERROR);
    }

    for (ip = 0; ip < l1rec->npix; ip++) {

        ipb = ip*nbands;

        /* grab the remote-sensing reflectances for this pixel */
        for (ic = 0; ic < 6; ic++)
            rrs[ic] = (double) (l2rec->Rrs[ipb + indx[ic]]);

        /* Use retrieved sst, OR reference sst if no retrieval */
        if (l2rec->sst != NULL && l2rec->sst[ip] > -2.0)
            sst = l2rec->sst[ip];
        else
            sst = l1rec->sstref[ip];

        /* grab ndt for this lon/lat */
        ndt = get_ndt(l1rec->lon[ip], l1rec->lat[ip]);

        /* run the model */
        status = carder_model(nbands, rrs, rrswl, sst, ndt, atot1, aph1, adg1, bb1, bbp1, &aph_mod, &adg_mod, &chl_mod, &flags);

        /* store results of this pixel into static global arrays */
        if (status == 0) {
            for (ic = 0; ic < 6; ic++) {
                atot[ipb + indx[ic]] = atot1[ic];
                aph [ipb + indx[ic]] = aph1 [ic];
                adg [ipb + indx[ic]] = adg1 [ic];
                bb [ipb + indx[ic]] = bb1 [ic];
                bbp [ipb + indx[ic]] = bbp1 [ic];
            }
            flagm[ip] = flags;
            adgm [ip] = adg_mod;
            aphm [ip] = aph_mod;
            chl [ip] = chl_mod;
        } else {
            for (ic = 0; ic < 6; ic++) {
                atot[ipb + indx[ic]] = CARDER_FAIL;
                aph [ipb + indx[ic]] = CARDER_FAIL;
                adg [ipb + indx[ic]] = CARDER_FAIL;
                bb [ipb + indx[ic]] = CARDER_FAIL;
                bbp [ipb + indx[ic]] = CARDER_FAIL;
            }
            flagm[ip] = flags;
            adgm [ip] = CARDER_FAIL;
            aphm [ip] = CARDER_FAIL;
            chl [ip] = CARDER_FAIL;
            l1rec->flags[ip] |= PRODFAIL;
        }
    }

    CarderRecNum = l1rec->iscan;

    free(rrs);
    free(aph1);
    free(adg1);
    free(atot1);
    free(bb1);
    free(bbp1);
}

/* Interface to l2_hdf_generic() to return Carder flags */

int16 *get_flags_carder(l2str *l2rec) {
    if (!carder_ran(l2rec->l1rec->iscan))
        carder(l2rec);

    return (flagm);
}

/* Interface to l2_hdf_generic() to return Carder products */

void get_carder(l2str *l2rec, l2prodstr *p, float prod[]) {
    int band = p->prod_ix;
    int prodID = p->cat_ix;
    int ip, ipb;

    if (!carder_ran(l2rec->l1rec->iscan))
        carder(l2rec);

    for (ip = 0; ip < l2rec->l1rec->npix; ip++) {

        ipb = ip * l2rec->l1rec->l1file->nbands + band;

        switch (prodID) {

        case CAT_chl_carder:
            prod[ip] = (float) chl[ip];
            break;

        case CAT_adg_carder:
            if (band < 0)
                prod[ip] = (float) adgm[ip];
            else
                prod[ip] = (float) adg[ipb];
            break;

        case CAT_aph_carder:
            if (band < 0)
                prod[ip] = (float) aphm[ip];
            else
                prod[ip] = (float) aph[ipb];
            break;

        case CAT_a_carder:
            prod[ip] = (float) atot[ipb];
            break;

        case CAT_bbp_carder:
            prod[ip] = (float) bbp[ipb];
            break;

        case CAT_bb_carder:
            prod[ip] = (float) bb[ipb];
            break;

        default:
            printf("-E- %s line %d : erroneous product ID %d passed to get_carder().\n",
                    __FILE__, __LINE__, prodID);
            exit(1);
        }
    }

    return;
}

/* Interface to convl12() to return Carder iops */

void iops_carder(l2str *l2rec) {
    int32_t ib, ip, ipb;
    int32_t nbands = l2rec->l1rec->l1file->nbands;

    if (!carder_ran(l2rec->l1rec->iscan))
        carder(l2rec);

    for (ip = 0; ip < l2rec->l1rec->npix; ip++) for (ib = 0; ib < nbands; ib++) {
            ipb = ip * nbands + ib;
            l2rec->a [ipb] = (float) atot[ipb];
            l2rec->bb[ipb] = (float) bb [ipb];
        }

    return;
}



/* ========================================================================= */

/* ========================================================================= */
int carder_empirical(float *rrs, float sst, float ndt, float *chlPtr) {
    static int firstCall = 1;

    float exponent;
    float abr35;
    float tchl;
    float chl;
    float weit; /* package switching */

    float uptemp, gtemp, ptemp;
    float sw0lo, sw1lo, sw2lo, sw3lo;
    float dtsw1, dtsw2, dtsw3;
    int it, pktran[2];

    if (firstCall) {

        firstCall = 0;

        /* read model parameters */
        get_params();
    }

    /* initialize */
    *chlPtr = CARDER_FAIL;

    /* if we have any negative reflectances, just bail */
    if (rrs[2] <= 0. || rrs[4] <= 0.) {
        return 1;
    }

    /* calculate bandratio */
    abr35 = log10(rrs[2] / rrs[4]);

    /* set up the limits based on NDT */
    uptemp = ndt + 2.;
    gtemp = ndt + 0.8;
    ptemp = ndt - 1.;

    sw0lo = uptemp + 1.;
    sw1lo = gtemp + (uptemp - gtemp) / 2.;
    sw2lo = ptemp + (gtemp - ptemp) / 2.;
    sw3lo = ptemp - 1.;

    dtsw1 = sw0lo - sw1lo;
    dtsw2 = sw1lo - sw2lo;
    dtsw3 = sw2lo - sw3lo;

    /* sw = 0 */
    if (sst >= sw0lo) {
        weit = 1.;
        pktran[0] = 1;
        pktran[1] = 1;
    }

    /* sw = 1 */
    if (sst >= sw1lo && sst < sw0lo) {
        weit = (sst - sw1lo) / dtsw1;
        pktran[0] = 0;
        pktran[1] = 1;
    }

    /* sw = 2 */
    if (sst >= sw2lo && sst < sw1lo) {
        weit = (sst - sw2lo) / dtsw2;
        pktran[0] = 2;
        pktran[1] = 0;
    }

    /* sw = 3 */
    if (sst >= sw3lo && sst < sw2lo) {
        weit = (sst - sw3lo) / dtsw3;
        pktran[0] = 3;
        pktran[1] = 2;
    }

    /* sw = 4 */
    if (sst < sw3lo) {
        weit = 1.;
        pktran[0] = 3;
        pktran[1] = 3;
    }


    /* Loop through twice, to blend packaging effects */
    for (it = 0; it <= 1; it++) {

        pk = pktran[it];

        switch (pk) {
        case 1: /* unpackaged */
            c0 = upc0;
            c1 = upc1;
            c2 = upc2;
            c3 = upc3;
            break;
        case 2: /* packaged */
            c0 = pc0;
            c1 = pc1;
            c2 = pc2;
            c3 = pc3;
            break;
        case 3: /*hipackaged */
            c0 = hpc0;
            c1 = hpc1;
            c2 = hpc2;
            c3 = hpc3;
            break;
        default: /* global */
            c0 = gc0;
            c1 = gc1;
            c2 = gc2;
            c3 = gc3;
            break;
        }


        /* calculate default chl, note coeffs depend on packaging */
        exponent = MAX(MIN(c0 + c1 * abr35 + c2 * abr35 * abr35 + c3 * abr35 * abr35*abr35, 10.0), -10.0);
        tchl = pow(10.0, exponent);

        /* weight results */
        if (it == 0) {
            chl = tchl;
        } else {
            chl = chl * (1. - weit) + tchl * weit;
        }
    }

    *chlPtr = chl;

    return 0;

}

void chl_carder_empirical(l2str *l2rec, float prod[]) {
    static int *indx;
    static int firstCall = 1;

    int ip, ic;
    int32_t ipb;
    float ndt;
    float sst;
    int status;
    float rrs[6];
    float chl;

    l1str *l1rec = l2rec->l1rec;
    int32_t nbands = l1rec->l1file->nbands;

    if (firstCall) {
        firstCall = 0;
        if ((indx = (int *) calloc(nbands, sizeof (int))) == NULL) {
            printf("-E- : Error allocating memory to indx\n");
            exit(FATAL_ERROR);
        }
        indx[0] = bindex_get(412);
        indx[1] = bindex_get(443);
        indx[2] = bindex_get(490);
        indx[3] = bindex_get(510);
        if (indx[3] < 0) indx[3] = bindex_get(520);
        if (indx[3] < 0) indx[3] = bindex_get(530);
        indx[4] = bindex_get(551);
        indx[5] = bindex_get(670);
    }

    for (ip = 0; ip < l1rec->npix; ip++) {

        ipb = ip*nbands;

        /* grab the remote-sensing reflectances for this pixel */
        for (ic = 0; ic < 6; ic++)
            rrs[ic] = (l2rec->Rrs[ipb + indx[ic]]);

        /* Use retrieved sst, OR reference sst if no retrieval */
        if (l2rec->sst != NULL && l2rec->sst[ip] > -2.0)
            sst = l2rec->sst[ip];
        else
            sst = l1rec->sstref[ip];

        /* grab ndt for this lon/lat */
        ndt = get_ndt(l1rec->lon[ip], l1rec->lat[ip]);

        /* run the model */
        status = carder_empirical(rrs, sst, ndt, &chl);

        /* store results of this pixel into static global arrays */
        if (status == 0) {
            prod[ip] = chl;
        } else {
            prod[ip] = BAD_FLT;
            l1rec->flags[ip] |= PRODFAIL;
        }
    }
}


