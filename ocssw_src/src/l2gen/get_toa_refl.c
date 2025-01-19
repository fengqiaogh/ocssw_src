/* --------------------------------------------------------------- */
/* get_toa_refl.c - compute top-of-atmosphere reflectance.         */
/*                                                                 */
/* Inputs:                                                         */
/*     l2rec - level-2 structure containing one complete scan      */
/*             after atmospheric correction.                       */
/*     band  = band number 0-7                                     */
/* Outputs:                                                        */
/*     rhot  - toa reflectance                                     */
/*                                                                 */
/* Written By: B. Franz, SAIC GSC, SIMBIOS Project, 11 April 2000  */
/*                                                                 */
/* --------------------------------------------------------------- */

#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"

void get_toa_refl(l2str *l2rec, int band, float rhot[]) {
    static float pi = 3.141592654;

    float mu0;
    int32_t ip, ipb;

    l1str *l1rec = l2rec->l1rec;
    int32_t nbands = l1rec->l1file->nbands;

    for (ip = 0; ip < l1rec->npix; ip++) {
        ipb = ip * nbands + band;
        mu0 = l1rec->csolz[ip];
        rhot[ip] = pi * l2rec->l1rec->Lt[ipb] / l1rec->Fo[band] / mu0;
    }

    return;
}
