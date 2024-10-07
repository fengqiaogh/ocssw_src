/* --------------------------------------------------------------- */
/* get_bsi.c - Biogenic Silica                                     */
/*                                                                 */
/* Inputs:                                                         */
/*     l2rec - level-2 structure containing one complete scan      */
/*             after atmospheric correction.                       */
/* Outputs:                                                        */
/*     BSi in micromol/L                                           */
/*                                                                 */
/* Algorithm Provided By: W. Balch                                 */
/* Reference: Balch, W.M., Bowler, B.C., Drapeau, D.T., Poulton,   */
/* A. and Holligan, P., 2010. Biominerals and the vertical flux of */
/* particulate organic carbon from the surface ocean. Geochemical  */
/* Research Letters, 37(L22605): 1-6.                              */
/*                                                                 */
/* --------------------------------------------------------------- */

#include "l12_proto.h"

void get_bsi(l2str *l2rec, float *BSi) {
    static float c[3] = {0.415, -0.0052, 0.1073};

    int32_t ip;
    float sst;
    float chl;

    l1str *l1rec = l2rec->l1rec;

    for (ip = 0; ip < l1rec->npix; ip++) {

        BSi[ip] = BAD_FLT;

        if (l1rec->mask[ip]) {
            l1rec->flags[ip] |= PRODFAIL;
            continue;
        }

        sst = l1rec->sstref[ip];
        chl = l2rec->chl[ip];

        if (sst < -2 || chl < 0) {
            l1rec->flags[ip] |= PRODFAIL;
            continue;
        }

        BSi[ip] = c[0] * chl + c[1] * sst + c[2];
    }
}
