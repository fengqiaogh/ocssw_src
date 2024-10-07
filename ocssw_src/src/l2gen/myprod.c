#include "l12_proto.h"

/* ==================================================================== */
/* This module provides stubs for implementing user-defined algorithms. */
/* The products can than be output by MSL12 by specifying the name on   */
/* the l2prod list (e.g., l2prod=chlor_a,myprod1,myprod2). Users can    */
/* maintain this file and drop it into future MSL12 updates.   BAF      */
/* ==================================================================== */


/* --------------------------------------------------------------------- */
/* dummy program computes wavelength weighted integral of Rrs            */

/* --------------------------------------------------------------------- */

void myprod1(l2str *l2rec, float prod[]) {
    static int firstCall = 1;
    static int32_t mask = ATMFAIL | LAND | HIGLINT | CLOUD; /* l2_flags.h */

    l1str *l1rec = l2rec->l1rec;

    int32_t npix = l1rec->npix;
    int32_t nbands = l1rec->l1file->nbands;
    int32_t ip, ib, ipb;
    float sumwave;
    float wave;

    /* arrays of pixel interlaced by band (see l2_struc.h, alloc_l2.c) */

    //float *nLw = l2rec->nLw; /* normalized water leaving radiances  */
    float *Rrs = l2rec->Rrs; /* remote sensiing reflectances        */


    if (firstCall) {

        /* do initializations and allocations here */

        firstCall = 0;
    }

    for (ip = 0; ip < npix; ip++) {

        prod[ip] = 0.0;
        sumwave = 0.0;

        /* skip masked pixels */
        if ((l1rec->flags[ip] & mask) != 0)
            continue;

        for (ib = 0; ib < nbands; ib++) {
            ipb = ip * nbands + ib; /* index of band at pixel */
            if (Rrs[ipb] > 0.0) {
                wave = l1rec->l1file->fwave[ib];
                prod[ip] += (Rrs[ipb] * wave);
                sumwave += wave;
            }
        }

        printf("%f %f\n", prod[ip], sumwave);
        if (sumwave > 0.0)
            prod[ip] /= sumwave;
        else
            /* if you want to bin, you need a flag to indicate failure */
            /* we use the PRODFAIL bit for general products            */
            l1rec->flags[ip] |= PRODFAIL;

    }
    return;
}



/* --------------------------------------------------------------------- */
/*                      put your algorithm here                          */

/* --------------------------------------------------------------------- */
void myprod2(l2str *l2rec, float prod[]) {
}

/* --------------------------------------------------------------------- */
/*                      put your algorithm here                          */

/* --------------------------------------------------------------------- */
void myprod3(l2str *l2rec, float prod[]) {
}

/* --------------------------------------------------------------------- */
/*                      put your algorithm here                          */

/* --------------------------------------------------------------------- */
void myprod4(l2str *l2rec, float prod[]) {
}

/* --------------------------------------------------------------------- */
/*                      put your algorithm here                          */

/* --------------------------------------------------------------------- */
void myprod5(l2str *l2rec, float prod[]) {
}

/* --------------------------------------------------------------------- */
/*                      put your algorithm here                          */

/* --------------------------------------------------------------------- */
void myprod6(l2str *l2rec, float prod[]) {
}

/* --------------------------------------------------------------------- */
/*                      put your algorithm here                          */

/* --------------------------------------------------------------------- */
void myprod7(l2str *l2rec, float prod[]) {
}

/* --------------------------------------------------------------------- */
/*                      put your algorithm here                          */

/* --------------------------------------------------------------------- */
void myprod8(l2str *l2rec, float prod[]) {
}

/* --------------------------------------------------------------------- */
/*                      put your algorithm here                          */

/* --------------------------------------------------------------------- */
void myprod9(l2str *l2rec, float prod[]) {
}

/* --------------------------------------------------------------------- */
/*                      put your algorithm here                          */

/* --------------------------------------------------------------------- */
void myprod10(l2str *l2rec, float prod[]) {
}

