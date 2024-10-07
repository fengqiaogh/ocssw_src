/* --------------------------------------------------------------- */
/* get_tricho.c - tricho classification function for MSl12.        */
/*                                                                 */
/* Inputs:                                                         */
/*     l2rec - level-2 structure containing one complete scan      */
/*             after atmospheric correction.                       */
/* Outputs:                                                        */
/*     tricho - flag 0=no, 1=yes                                   */
/*                                                                 */
/* Algorithm Provided By: Ajit Subramaniam.                        */
/* Written By: B. A. Franz, SAIC GSC, SIMBIOS Project, Sept 1999   */
/*                                                                 */
/* --------------------------------------------------------------- */

#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"

void get_tricho(l2str *l2rec, uint8 tricho[]) {

    int32_t ip, ipb;
    float de, nu, tdn;
    l1str *l1rec = l2rec->l1rec;

    for (ip = 0; ip < l1rec->npix; ip++) {

        ipb = l1rec->l1file->nbands*ip;

        tricho[ip] = 0;

        /* Skip shallow water or negative nLw conditions, */
        /* high winds, or cold water                      */
        if ((l1rec->flags && COASTZ) > 0 ||
                l1rec->ws[ip] > 8.0 ||
                l1rec->sstref[ip] < 25.0 ||
                l2rec->nLw[ipb + 1] <= 0.0 || /* 443 */
                l2rec->nLw[ipb + 2] <= 0.0 || /* 490 */
                l2rec->nLw[ipb + 4] <= 0.0) { /* 555 */
            continue;
        }

        de = l2rec->nLw[ipb + 2] - l2rec->nLw[ipb + 1];
        nu = l2rec->nLw[ipb + 2] - l2rec->nLw[ipb + 4];

        if (de > 0 && nu > 0) {
            tdn = de / nu;
            if ((l2rec->nLw[ipb + 2] > 1.4 && l2rec->nLw[ipb + 2] < 3.0) &&
                    (tdn > 0.4 && tdn < 0.6))
                tricho[ip] = 1;
        }

    }
}
