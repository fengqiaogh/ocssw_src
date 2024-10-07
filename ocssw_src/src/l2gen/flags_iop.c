#include "l12_proto.h"

/* --------------------------------------------------------------------- */

/* --------------------------------------------------------------------- */
void init_iop_flag(float32 wave[], int32 nwave, iopfstr *f) {
    float32 *aw;
    float32 *bbw;
    float32 a_hi = 5.0;
    float32 b_hi = 0.05;
    float32 w_hi = 600;

    int iw;

    if ((aw = (float32 *) calloc(nwave, sizeof (float32))) == NULL) {
        printf("-E- %s line %d : error allocating memory for aw in flags_iop.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((bbw = (float32 *) calloc(nwave, sizeof (float32))) == NULL) {
        printf("-E- %s line %d : error allocating memory for bbw in flags_iop.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    for (iw = 0; iw < nwave; iw++) {
        if (wave[iw] > w_hi) break;
        aw [iw] = aw_spectra(wave[iw], BANDW);
        bbw[iw] = bbw_spectra(wave[iw], BANDW);
    }
    f->nwave = iw;

    if ((f->a_lo = (float32 *) calloc(nwave, sizeof (float32))) == NULL) {
        printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((f->a_hi = (float32 *) calloc(nwave, sizeof (float32))) == NULL) {
        printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((f->a_on = (float32 *) calloc(nwave, sizeof (float32))) == NULL) {
        printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    if ((f->aph_lo = (float32 *) calloc(nwave, sizeof (float32))) == NULL) {
        printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((f->aph_hi = (float32 *) calloc(nwave, sizeof (float32))) == NULL) {
        printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((f->aph_on = (float32 *) calloc(nwave, sizeof (float32))) == NULL) {
        printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    if ((f->adg_lo = (float32 *) calloc(nwave, sizeof (float32))) == NULL) {
        printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((f->adg_hi = (float32 *) calloc(nwave, sizeof (float32))) == NULL) {
        printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((f->adg_on = (float32 *) calloc(nwave, sizeof (float32))) == NULL) {
        printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    if ((f->bb_lo = (float32 *) calloc(nwave, sizeof (float32))) == NULL) {
        printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((f->bb_hi = (float32 *) calloc(nwave, sizeof (float32))) == NULL) {
        printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((f->bb_on = (float32 *) calloc(nwave, sizeof (float32))) == NULL) {
        printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    if ((f->bbp_lo = (float32 *) calloc(nwave, sizeof (float32))) == NULL) {
        printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((f->bbp_hi = (float32 *) calloc(nwave, sizeof (float32))) == NULL) {
        printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((f->bbp_on = (float32 *) calloc(nwave, sizeof (float32))) == NULL) {
        printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                __FILE__, __LINE__);
        exit(1);
    }


    for (iw = 0; iw < f->nwave; iw++) {

        f->a_lo[iw] = 0.95 * aw[iw];
        f->a_hi[iw] = a_hi;
        f->a_on[iw] = 1;

        f->aph_lo[iw] = -0.05 * aw[iw];
        f->aph_hi[iw] = a_hi;
        f->aph_on[iw] = 1;

        f->adg_lo[iw] = -0.05 * aw[iw];
        f->adg_hi[iw] = a_hi;
        f->adg_on[iw] = 1;

        f->bb_lo[iw] = 0.95 * bbw[iw];
        f->bb_hi[iw] = b_hi;
        f->bb_on[iw] = 1;

        f->bbp_lo[iw] = -0.05 * bbw[iw];
        f->bbp_hi[iw] = b_hi;
        f->bbp_on[iw] = 1;
    }

    for (iw = f->nwave; iw < nwave; iw++) {
        f->a_on [iw] = 0;
        f->aph_on [iw] = 0;
        f->adg_on [iw] = 0;
        f->bb_on [iw] = 0;
        f->bbp_on [iw] = 0;
    }
    free(aw);
    free(bbw);
}


/* --------------------------------------------------------------------- */

/* --------------------------------------------------------------------- */
void set_iop_flag(float wave[], int32_t nwave,
        float a[], float aph[], float adg[],
        float bb[], float32 bbp[], int16_t *flag) {
    static int firstCall = 1;
    static iopfstr iopf;

    int iw;

    if (firstCall) {
        firstCall = 0;
        init_iop_flag(wave, nwave, &iopf);
    }

    for (iw = 0; iw < iopf.nwave; iw++) {

        if (iopf.a_on[iw] && a != NULL) {
            if (a[iw] < iopf.a_lo[iw])
                *flag |= IOPF_ALO;
            if (a[iw] > iopf.a_hi[iw])
                *flag |= IOPF_AHI;
        }

        if (iopf.aph_on[iw] && aph != NULL) {
            if (aph[iw] < iopf.aph_lo[iw])
                *flag |= IOPF_APHLO;
            if (aph[iw] > iopf.aph_hi[iw])
                *flag |= IOPF_APHHI;
        }

        if (iopf.adg_on[iw] && adg != NULL) {
            if (adg[iw] < iopf.adg_lo[iw])
                *flag |= IOPF_ADGLO;
            if (adg[iw] > iopf.adg_hi[iw])
                *flag |= IOPF_ADGHI;
        }

        if (iopf.bb_on[iw] && bb != NULL) {
            if (bb[iw] < iopf.bb_lo[iw])
                *flag |= IOPF_BBLO;
            if (bb[iw] > iopf.bb_hi[iw])
                *flag |= IOPF_BBHI;
        }

        if (iopf.bbp_on[iw] && bbp != NULL) {
            if (bbp[iw] < iopf.bbp_lo[iw])
                *flag |= IOPF_BBPLO;
            if (bbp[iw] > iopf.bbp_hi[iw])
                *flag |= IOPF_BBPHI;
        }
    }
}
