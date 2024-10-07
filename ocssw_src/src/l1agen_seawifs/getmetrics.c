/* ---------------------------------------------------------------- */
/* Module getmetrics - computes seawifs L1A image metrics info      */
/*                                                                  */
/* B. A. Franz, GSC, 11/97                                          */
/* ---------------------------------------------------------------- */

#include "swl0_proto.h"

static l1met totmet; /* Total metrics */
static l1met avgmet; /* Mean  metrics */

/* ---------------------------------------------------------------- */
/* addL1Metrics() - adds metrics from one L1A scan to total         */

/* ---------------------------------------------------------------- */
void addL1Metrics(INT32 scanNum, swl1rec *l1rec) {
    int ip, ib;

    if (scanNum <= 0)
        memset(&totmet, 0, sizeof (totmet));

    for (ip = 0; ip < l1rec->npix; ip++) {

        for (ib = 0; ib < NBANDS; ib++) {

            if ((l1rec->data[ip][ib] - l1rec->darkpix[ib]) < 3)
                totmet.zeropix[ib]++;

            if (l1rec->gain[ib] == 0) { /* Earth Gain 1 */

                if (l1rec->data[ip][ib] >= 1023)
                    totmet.gain1_satpix[ib]++;
                else
                    totmet.gain1_nonsatpix[ib]++;

                totmet.gain1_mean_rad[ib] += l1rec->data[ip][ib];

            } else if (l1rec->gain[ib] == 2) { /* Earth Gain 2 */

                if (l1rec->data[ip][ib] >= 1023)
                    totmet.gain2_satpix[ib]++;
                else
                    totmet.gain2_nonsatpix[ib]++;

                totmet.gain2_mean_rad[ib] += l1rec->data[ip][ib];

            }
        }

    }
    return;
}



/* ---------------------------------------------------------------- */
/* getL1metrics() - returns pointer to mean metrics structure       */

/* ---------------------------------------------------------------- */
l1met *getL1Metrics(void) {
    int ib;
    INT32 gain1Cnt;
    INT32 gain2Cnt;

    memcpy(&avgmet, &totmet, sizeof (l1met));

    for (ib = 0; ib < NBANDS; ib++) {

        gain1Cnt = totmet.gain1_satpix[ib]
                + totmet.gain1_nonsatpix[ib];

        gain2Cnt = totmet.gain2_satpix[ib]
                + totmet.gain2_nonsatpix[ib];

        if (gain1Cnt > 0)
            avgmet.gain1_mean_rad[ib] /= gain1Cnt;

        if (gain2Cnt > 0)
            avgmet.gain2_mean_rad[ib] /= gain2Cnt;
    }

    /*
    for (ib=0; ib<NBANDS; ib++) 
        printf("%2d %7d %7d %7d %7d %7d %10.2f %10.2f\n",ib,
            avgmet.gain1_satpix[ib],
            avgmet.gain2_satpix[ib],
            avgmet.gain1_nonsatpix[ib],
            avgmet.gain2_nonsatpix[ib],
            avgmet.zeropix[ib],
            avgmet.gain1_mean_rad[ib],
            avgmet.gain2_mean_rad[ib]);
     */

    return (&avgmet);
}



