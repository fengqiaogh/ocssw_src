/* --------------------------------------------------------------- */
/* get_es.c - compute surface irradiance (Ed(0+))                  */
/*                                                                 */
/* Inputs:                                                         */
/*     l2rec - level-2 structure containing one complete scan      */
/*             after atmospheric correction.                       */
/*     band  - band number (0-7)                                   */
/*                                                                 */
/* Outputs:                                                        */
/*     Es for specified band                                       */
/*                                                                 */
/* Algorithm Provided By: M. Wang                                  */
/* Written By: B. A. Franz, SAIC GSC, SIMBIOS Project, 4 Aug 1999  */
/*                                                                 */
/* --------------------------------------------------------------- */

#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"

/**
 * now returns in W/m2/um
 */
void get_es(l2str *l2rec, int band, float Es[]) {
    int32_t ip;
    int32_t ipb;

    l1str *l1rec = l2rec->l1rec;
    int32_t nbands = l1rec->l1file->nbands;

    for (ip = 0; ip < l1rec->npix; ip++) {
        ipb = ip * nbands + band;
        if (l2rec->La[ipb] > 0.0) {
            Es[ip] = l1rec->Fo[band]
                    * l1rec->tg_sol[ipb]
                    * l1rec->t_sol[ipb]
                    * l1rec->csolz[ip]
                    * 10.0;
        } else
            Es[ip] = 0.0;
    }

    return;
}

