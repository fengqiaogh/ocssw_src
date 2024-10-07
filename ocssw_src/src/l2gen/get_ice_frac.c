/* --------------------------------------------------------------- */
/* get_ice_frac.c - dump the ice fraction ancillary data           */
/*                                                                 */
/* Inputs:                                                         */
/*     l2rec - level-2 structure containing one complete scan      */
/*             after atmospheric correction.                       */
/*                                                                 */
/* Outputs:                                                        */
/*     ice - array where ice data is written                       */
/*                                                                 */
/* Written By: Don Shea, OBPG,  23 Apr 2009                        */
/*                                                                 */
/* --------------------------------------------------------------- */

#include "l12_proto.h"

void get_ice_frac(l2str *l2rec, float ice[]) {
    int32_t ip;
    l1str *l1rec = l2rec->l1rec;

/*  WDR make current set-up so that the proc_cloud ice will not be blocked over land */
    for (ip = 0; ip < l1rec->npix; ip++) {
        if( input->proc_cloud ) {
            ice[ip] = l1rec->icefr[ip];
        } else {
            if (l1rec->flags[ip] & LAND) {
                ice[ip] = 0.0;
            } else {
                ice[ip] = ice_fraction(l1rec->lon[ip], l1rec->lat[ip]);
            }
        }
    }

}

