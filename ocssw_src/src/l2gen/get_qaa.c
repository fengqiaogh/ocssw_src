/*
 *  get_qaa.c
 *
 *        MSl12 wrapper for Quasi-Analytic Algorithm
 *
 *  Naval Research Laboratory
 *  Stennis Space Center, MS
 */

#include <stdlib.h>
#include <math.h>
#include "qaa.h"
#include "l12_proto.h"
#include "l2prod.h"

static int QaaRecNum = -1;

static int nbandsVis; /* number of bands computed for line       */

/* one scan line of data */

static float *atot; /* total absorption coefficient            */
static float *adg; /* detrital absorption coefficient         */
static float *aph; /* phytoplankton absorption coefficient    */
static float *bb; /* backscatter coefficient                 */
static unsigned char *flags; /* per-pixel flags                         */

/* pixel data set during initialization */

static float *bbw; /* pure-water backscattering               */
static float *aw; /* pure-water total absorption             */
static float *fwave; /* wavelength (nm)                         */

/* pixel data */

static float *Rrs; /* above-water remote sensing reflectance  */
static float *rrs; /* below-water remote sensing reflectance  */
static float *u;
static float *a_qaa; /* total absorption coefficient            */
static float *bb_qaa; /* backscatter coefficient                 */
static float *bbp_qaa; /* backscatter coefficient                 */
static float *adg_qaa; /* detrital absorption coefficient         */
static float *aph_qaa; /* phytoplankton absorption coefficient    */

static int ib410, ib440, ib490, ib555, ib670;
static int do_decomp = 1;

/* have we run for this scan line? */

static int qaa_ran(int recnum) {
    if (recnum == QaaRecNum)
        return 1;
    else
        return 0;
}

/* 
 * allocate private arrays for a single scan line (412 ... 555) 
 * check to see if the number of pixels in the line got bigger
 */
static void qaa_alloc(int32_t npix, int32_t nbands) {
    static int32_t currentNumPix = -1;

    if (npix > currentNumPix) {
        if (currentNumPix != -1) {
            free(atot);
            free(aph);
            free(adg);
            free(bb);
            free(flags);
        }
        currentNumPix = npix;
        atot = (float*) calloc(npix*nbands, sizeof (float));
        aph = (float*) calloc(npix*nbands, sizeof (float));
        adg = (float*) calloc(npix*nbands, sizeof (float));
        bb = (float*) calloc(npix*nbands, sizeof (float));
        flags = (unsigned char*) calloc(npix, sizeof (unsigned char));
    }
}

/* allocate private arrays for a single pixel (may include a 640nm) */

static void qaa_pixel_alloc(int nbands) {

    fwave = (float*) calloc(nbands, sizeof (float));
    bbw = (float*) calloc(nbands, sizeof (float));
    aw = (float*) calloc(nbands, sizeof (float));
    Rrs = (float*) calloc(nbands, sizeof (float));
    rrs = (float*) calloc(nbands, sizeof (float));
    u = (float*) calloc(nbands, sizeof (float));
    a_qaa = (float*) calloc(nbands, sizeof (float));
    bb_qaa = (float*) calloc(nbands, sizeof (float));
    bbp_qaa = (float*) calloc(nbands, sizeof (float));
    adg_qaa = (float*) calloc(nbands, sizeof (float));
    aph_qaa = (float*) calloc(nbands, sizeof (float));

}

/* NOTE: nbandsVis and l2rec->l1rec->nbands are not always equal */

static void run_qaa(l2str *l2rec) {
    static int firstCall = 1;

    l1str *l1rec = l2rec->l1rec;

    int ip, ib, ipb;
    int i;
    unsigned char flags_qaa;


    if (firstCall) {

        int *w;

        firstCall = 0;

        /* limit to visible wavelengths */
        nbandsVis = rdsensorinfo(l1rec->l1file->sensorID, l1_input->evalmask, "NbandsVIS", NULL);
        printf("QAA v6 processing for %d bands\n", nbandsVis);

        qaa_pixel_alloc(nbandsVis);

        for (ib = 0; ib < nbandsVis; ib++)
            fwave[ib] = l1rec->l1file->fwave[ib];
        get_aw_bbw(l2rec, fwave, nbandsVis, aw, bbw);

        w = input->qaa_wave;
        if (w[1] < 0 || w[2] < 0 || w[3] < 0 || w[4] < 0) {
            printf("qaa: algorithm coefficients not provided for this sensor.\n");
            exit(1);
        }

        // qaaf_v6 needs these bands at the very minimum

        ib440 = bindex_get(w[1]);
        ib490 = bindex_get(w[2]);
        ib555 = bindex_get(w[3]);
        ib670 = bindex_get(w[4]);
        if (ib440 < 0 || ib490 < 0 || ib555 < 0 || ib670 < 0) {
            printf("get_qaa: missing minimum required wavelengths "
                    "(need 440,490,555,670).\n");
            printf("get_qaa: qaa_wave[1] =%d, qaa_wave[2] =%d, "
                    "qaa_wave[3] = %d, qaa_wave[4] = %d.\n",
                    w[1], w[2], w[3], w[4]);
            exit(1);
        }

        // if we want compute the aph/adg (call qaaf_decomp) we need this additional band

        ib410 = bindex_get(w[0]);
        if (ib410 < 0) {
            printf("get_qaa: incompatible sensor wavelengths for aph/adg (need 410).\n");
            do_decomp = 0;
        }
        printf("QAA v6 bands: (%d) %d nm, (%d) %d nm, (%d) %d nm, (%d) %d nm, (%d) %d nm\n",
                ib410, w[0], ib440, w[1], ib490, w[2], ib555, w[3], ib670, w[4]);
        printf("QAA v6 wav  :");
        for (ib = 0; ib < nbandsVis; ib++)
            printf(" %10.6f", fwave[ib]);
        printf("\n");

        printf("QAA v6 aw   :");
        for (ib = 0; ib < nbandsVis; ib++)
            printf(" %10.6f", aw[ib]);
        printf("\n");

        printf("QAA v6 bbw  :");
        for (ib = 0; ib < nbandsVis; ib++)
            printf(" %10.6f", bbw[ib]);
        printf("\n");

        qaa_init(ib410, ib440, ib490, ib555, ib670);
        qaa_set_param(QAA_S_PARAM, input->qaa_adg_s);
    }

    // note that for l3gen the number of pixels changes
    qaa_alloc(l1rec->npix, nbandsVis);

    for (ip = 0; ip < l1rec->npix; ip++) {

        /* clear static globals */
        for (ib = 0; ib < nbandsVis; ib++) {
            ipb = ip * nbandsVis + ib;
            bb [ipb] = -0.1;
            atot[ipb] = -0.1;
            adg [ipb] = -0.1;
            aph [ipb] = -0.1;
        }
        flags[ip] = 0;

        if (!l1rec->mask[ip]) {

            flags_qaa = 0;
            for (i = 0; i < nbandsVis; i++)
                Rrs[i] = l2rec->Rrs[ip * l1rec->l1file->nbands + i];

            // Version 6
            qaaf_v6(nbandsVis, fwave, Rrs, aw, bbw, rrs, u, a_qaa, bb_qaa, &flags_qaa);
            if (do_decomp)
                qaaf_decomp(nbandsVis, fwave, rrs, a_qaa, aw, adg_qaa, aph_qaa,
                    &flags_qaa);

            /* store results for this pixel in static globals */

            for (ib = 0; ib < nbandsVis; ib++) {

                ipb = ip * nbandsVis + ib;

                if (isfinite(bb_qaa[ib]))
                    bb[ipb] = bb_qaa[ib];
                else {
                    bb[ipb] = -0.1;
                    l1rec->flags[ip] |= PRODFAIL;
                }

                if (isfinite(a_qaa[ib]))
                    atot[ipb] = a_qaa[ib];
                else {
                    atot[ipb] = -0.1;
                    l1rec->flags[ip] |= PRODFAIL;
                }

                if (isfinite(adg_qaa[ib]))
                    adg[ipb] = adg_qaa[ib];
                else {
                    adg[ipb] = -0.1;
                    l1rec->flags[ip] |= PRODFAIL;
                }

                if (isfinite(aph_qaa[ib]))
                    aph[ipb] = aph_qaa[ib];
                else {
                    aph[ipb] = -0.1;
                    l1rec->flags[ip] |= PRODFAIL;
                }

                /* ZP Lee, 17 August 2007 */
                if (atot[ipb] > 0.0) atot[ipb] = MAX(atot[ipb], aw [ib]*1.05);
                if (bb [ipb] > 0.0) bb [ipb] = MAX(bb [ipb], bbw[ib]*1.05);
            }
            flags[ip] = flags_qaa;

        }
    }

    QaaRecNum = l1rec->iscan;

    return;
}

/* interface to l2_hdf_generic() to return QAA flags */

unsigned char *get_flags_qaa(l2str *l2rec) {
    if (!qaa_ran(l2rec->l1rec->iscan))
        run_qaa(l2rec);

    return (flags);
}

void get_qaa(l2str *l2rec, l2prodstr *p, float prod[]) {
    int prodID = p->cat_ix;
    int ib = p->prod_ix;
    int ip, ipb;

    if (!qaa_ran(l2rec->l1rec->iscan))
        run_qaa(l2rec);

    for (ip = 0; ip < l2rec->l1rec->npix; ip++) {

        ipb = ip * nbandsVis + ib;

        switch (prodID) {

        case CAT_a_qaa:
            if (atot[ipb] > 0.0)
                prod[ip] = atot[ipb];
            else
                prod[ip] = p->badData;
            break;

        case CAT_adg_qaa:
            if (adg[ipb] > 0.0)
                prod[ip] = adg[ipb];
            else
                prod[ip] = p->badData;
            break;

        case CAT_aph_qaa:
            if (aph[ipb] > 0.0)
                prod[ip] = aph[ipb];
            else
                prod[ip] = p->badData;
            break;

        case CAT_bb_qaa:
            if (bb[ipb] > 0.0)
                prod[ip] = bb[ipb];
            else
                prod[ip] = p->badData;
            break;

        case CAT_bbp_qaa:
            if ((bb[ipb] - bbw[ib]) > 0.0)
                prod[ip] = bb[ipb] - bbw[ib];
            else
                prod[ip] = p->badData;
            break;

        default:
            printf("-E- %s line %d : erroneous product ID %d passed to get_qaa().\n",
                    __FILE__, __LINE__, prodID);
            exit(1);
        }
    }

    return;
}

/* interface to convl12() to return QAA iops */

void iops_qaa(l2str *l2rec) {
    int ib, ip;

    if (!qaa_ran(l2rec->l1rec->iscan))
        run_qaa(l2rec);

    for (ip = 0; ip < l2rec->l1rec->npix; ip++) {
        for (ib = 0; ib < nbandsVis; ib++) {
            l2rec->a [ip * l2rec->l1rec->l1file->nbands + ib] = atot[ip * nbandsVis + ib];
            l2rec->bb[ip * l2rec->l1rec->l1file->nbands + ib] = bb [ip * nbandsVis + ib];
        }
    }

    return;
}

void qaa_iops_4_bshift(l2str *l2rec, float *adg_ref, float *bbp_ref) {
    int ipb, ip;

    if (!qaa_ran(l2rec->l1rec->iscan))
        run_qaa(l2rec);
    for (ip = 0; ip < l2rec->l1rec->npix; ip++) {
        ipb = nbandsVis * ip + ib440;
        adg_ref[ip] = adg[ipb];
        bbp_ref[ip] = bb[ipb] - bbw[ib440];
    }
}

/* ------------------------------------------------------------------------------- */
/* interface to giop()                                                             */

/* ------------------------------------------------------------------------------- */
int get_bbp_qaa(l2str *l2rec, int ip, float tab_wave[], float tab_bbp[], int tab_nwave) {
    int ipb, ib;

    if (!qaa_ran(l2rec->l1rec->iscan))
        run_qaa(l2rec);

    for (ib = 0; ib < nbandsVis; ib++) {
        ipb = ip * nbandsVis + ib;
        bbp_qaa[ib] = bb[ipb] - bbw[ib];
        if (bbp_qaa[ib] < 0)
            return (0);
    }

    for (ib = 0; ib < tab_nwave; ib++) {
        tab_bbp[ib] = linterp(fwave, bbp_qaa, nbandsVis, tab_wave[ib]);
    }

    return (1);
}
