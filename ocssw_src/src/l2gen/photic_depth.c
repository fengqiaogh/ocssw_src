#include "l12_proto.h"

/* =================================================================== */
/* Algorithms for computing photic depths                              */
/*                                                                     */
/* B. Franz, NASA Ocean Biology Processing Group, SAIC, March 2005.    */
/* =================================================================== */

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_poly.h>
#include "l12_proto.h"

#define Z_CHL_MAX 15.0

/* ------------------------------------------------------------------- */
/* Zphotic_lee() - photic depth (ZP Lee)                               */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/*     p     - product catalog entry                                   */
/*                                                                     */
/* Outputs:                                                            */
/*    Zp - euphtic depth, 1 value per pixel.                           */
/*                                                                     */
/* Description:                                                        */
/*                                                                     */
/* Reference:                                                          */
/* Lee, at al (2007) "Euphotic zone depth:  Its derivation and         */
/* implication to ocean-color remote sensing", JGR, 112, C03009,       */
/* pp 1-12.                                                            */
/*                                                                     */
/* Original Implementation: B. Franz, August 2007                      */

/*---------------------------------------------------------------------*/
void Zphotic_lee(l2str *l2rec, l2prodstr *p, float Zp[]) {
    static int ib490 = -1;

    float stheta = 0.0; /* use solz = 0.0 */
    float ctheta = 1.0; /* use solz = 0.0 */

    const float c1 [] = {-0.057, 0.482, 4.221};
    const float c2 [] = {0.183, 0.702, -2.567};
    const float alpha[] = {0.090, 1.465, -0.667};

    double k1, k2;
    double y1, y2, y3;
    double z1, z2, z3;
    double tau;
    float a490, bb490;
    int ip, ipb;

    l1str *l1rec = l2rec->l1rec;

    switch (p->prod_ix) {
    case -1: /*Zsd*/
        tau = -log(0.1);
        break;
    case -2: /*Zeu*/
        tau = -log(0.01);
        break;
    case -3: /* 1st optical depth */
        tau = 1.0;
        break;
    default:
        if (p->prod_ix <= 0 || p->prod_ix > 100) {
            printf("Zeu_lee: percent light should be between 1 and 100.\n");
            exit(1);
        }
        tau = -log(p->prod_ix / 100.0);
        break;
    }

    if (input->iop_opt == IOPNONE) {
        printf("IOP-based Z*_lee product requires iop model selection (iop_opt).  ");
        printf("Using default model.\n");
        input->iop_opt = IOPDEFAULT;
        get_iops(l2rec, input->iop_opt);
    }

    if (ib490 < 0) {
        ib490 = bindex_get(490);
        if (ib490 < 0) {
            printf("Zeu_lee: a 490nm channel is required for this algorithm\n");
            exit(1);
        }
    }

    for (ip = 0; ip < l1rec->npix; ip++) {

        ipb = ip * l1rec->l1file->nbands + ib490;

        a490 = l2rec->a[ipb];
        bb490 = l2rec->bb[ipb];

        Zp[ip] = p->badData;

        if (l1rec->mask[ip] || a490 <= 0.0 || bb490 <= 0.0) {
            l1rec->flags[ip] |= PRODFAIL;
        } else {
            stheta = sin(l1rec->solz[ip] / RADEG);
            ctheta = l1rec->csolz[ip];
            k1 = (c1[0] + c1[1] * sqrt(a490) + c1[2] * bb490)*(1 + alpha[0] * stheta);
            k2 = (c2[0] + c2[1]* (a490) + c2[2] * bb490)*(alpha[1] + alpha[2] * ctheta);
            y1 = (k1 * k1 - k2 * k2 - 2.0 * tau * k1) / (k1 * k1);
            y2 = (tau * tau - 2.0 * tau * k1) / (k1 * k1);
            y3 = (tau * tau) / (k1 * k1);
            gsl_poly_solve_cubic(y1, y2, y3, &z1, &z2, &z3);
            if (z2 > 0.0 && z3 > 0.0)
                Zp[ip] = MIN(z2, z3);
            else if (z1 > 0.0 && z2 > 0.0)
                Zp[ip] = MIN(z1, z2);
            else if (z1 > 0.0 && z3 > 0.0)
                Zp[ip] = MIN(z1, z3);
            else
                l1rec->flags[ip] |= PRODFAIL;
        }
    }

    return;
}


/* ------------------------------------------------------------------- */
/* zeu_morel() - depth of the euphotic layer (1%) by Morel             */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/*     p     - product catalog entry                                   */
/*                                                                     */
/* Outputs:                                                            */
/*    Zeu - euphtic depth, 1 value per pixel.                          */
/*                                                                     */
/* Description:                                                        */
/*                                                                     */
/* Reference:                                                          */
/*                                                                     */
/* Morel, A., Y. Huot, B. Gentili, P.J. Werdell, S.B. Hooker (2007).   */
/* Consistency of products derived from various ocean color sensors:   */
/* An examination before merging these products and extending their    */
/* applications, Remote Sensing of Environment, to be submitted.       */
/*                                                                     */
/* Original Implementation: B. Franz, November 2006                    */

/*---------------------------------------------------------------------*/
void Zeu_morel(l2str *l2rec, l2prodstr *p, float *Zeu) {
    float chl, x;
    int32_t ip;

    l1str *l1rec = l2rec->l1rec;

    for (ip = 0; ip < l1rec->npix; ip++) {

        chl = l2rec->chl[ip];

        if (chl > Z_CHL_MAX) {
            l1rec->flags[ip] |= PRODWARN;
        }

        if (l1rec->mask[ip] || chl <= 0.0) {
            Zeu[ip] = p->badData;
            l1rec->flags[ip] |= PRODFAIL;
        } else {
            x = log10(chl);
            Zeu[ip] = pow(10.0, 1.524 + x * (-0.460 + x * (-0.00051 + x * 0.0282)));
        }
    }

    return;
}


/* ------------------------------------------------------------------- */
/* Zsd_morel() - secchi depth by Morel                                 */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/*     p     - product catalog entry                                   */
/*                                                                     */
/* Outputs:                                                            */
/*    Zsd - secchi depth depth, 1 value per pixel.                     */
/*                                                                     */
/* Description:                                                        */
/*                                                                     */
/* Reference:                                                          */
/*                                                                     */
/* Morel, A., Y. Huot, B. Gentili, P.J. Werdell, S.B. Hooker (2007).   */
/* Consistency of products derived from various ocean color sensors:   */
/* An examination before merging these products and extending their    */
/* applications, Remote Sensing of Environment, to be submitted.       */
/*                                                                     */
/* Implementation: B. Franz, November 2006                             */

/*---------------------------------------------------------------------*/
void Zsd_morel(l2str *l2rec, l2prodstr *p, float *Zsd) {
    float chl, x;
    int32_t ip;

    l1str *l1rec = l2rec->l1rec;

    for (ip = 0; ip < l1rec->npix; ip++) {

        chl = l2rec->chl[ip];

        if (chl > Z_CHL_MAX) {
            l1rec->flags[ip] |= PRODWARN;
        }

        if (l1rec->mask[ip] || chl <= 0.0 || chl > Z_CHL_MAX) {
            Zsd[ip] = p->badData;
            l1rec->flags[ip] |= PRODFAIL;
        } else {
            x = log10(chl);
            Zsd[ip] = 8.50 + x * (-12.6 + x * (7.36 - x * 1.43));
        }
    }

    return;
}


/* ------------------------------------------------------------------- */
/* intchl_morel() - quasi depth integrated chl using morel             */

/*---------------------------------------------------------------------*/
void intchl_morel(l2str *l2rec, l2prodstr *p, float *intchl) {
    float Zeu, chl, x;
    int32_t ip;

    l1str *l1rec = l2rec->l1rec;

    for (ip = 0; ip < l1rec->npix; ip++) {

        chl = l2rec->chl[ip];

        if (chl > Z_CHL_MAX) {
            l1rec->flags[ip] |= PRODWARN;
        }

        if (l1rec->mask[ip] || chl <= 0.0) {
            intchl[ip] = p->badData;
            l1rec->flags[ip] |= PRODFAIL;
        } else {
            x = log10(chl);
            Zeu = pow(10.0, 1.524 + x * (-0.460 + x * (-0.00051 + x * 0.0282)));
            intchl[ip] = Zeu*chl;
        }
    }

    return;
}


/* ------------------------------------------------------------------- */

/*---------------------------------------------------------------------*/
void Zsd_gbr(l2str *l2rec, l2prodstr *p, float *Zsd) {
    static int firstCall = 1;
    static float a0, a1;

    int32_t ip;
    l1str *l1rec = l2rec->l1rec;

    if (firstCall) {
        firstCall = 0;
        switch (l1rec->l1file->sensorID) {
        case MODIST:
        case MODISA:
            a0 = 0.580;
            a1 = 0.742;
            break;
        default:
            a0 = 0.518;
            a1 = 0.811;
            break;
        }
    }

    p->prod_ix = -1; // 10% light
    Zphotic_lee(l2rec, p, Zsd);

    for (ip = 0; ip < l1rec->npix; ip++) {

        if (Zsd[ip] > 0.0)
            Zsd[ip] = pow(10.0, (log10(Zsd[ip]) - a0) / a1);
        else {
            Zsd[ip] = p->badData;
            l1rec->flags[ip] |= PRODFAIL;
        }
    }

    return;
}


/* ------------------------------------------------------------------- */
/* get_photic_depth() - l2_hdf_generic interface for photic depth      */

/* ------------------------------------------------------------------- */
void get_photic_depth(l2str *l2rec, l2prodstr *p, float prod[]) {

    switch (p->cat_ix) {
    case CAT_Zeu_morel:
        Zeu_morel(l2rec, p, prod);
        break;
    case CAT_Zsd_morel:
        Zsd_morel(l2rec, p, prod);
        break;
    case CAT_Zphotic_lee:
        Zphotic_lee(l2rec, p, prod);
        break;
    case CAT_Zsd_gbr:
        Zsd_gbr(l2rec, p, prod);
        break;
    }

    return;
}
