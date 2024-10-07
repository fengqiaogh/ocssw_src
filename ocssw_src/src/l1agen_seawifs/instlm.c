#include <stdio.h>
#include <math.h>
#include "instlm.h"
#include "swl0_proto.h"

/* -------------------------------------------------------------- */
/* valid_instlm() - returns 1 if scan contains valid instrument   */
/*                  telemetry.                                    */
/*                                                                */
/* INT16 mnftype - minor frame type (0=LAC, 15=GAC)               */
/* INT16 mnfnume - minor frame number 1,2, or 3                   */
/* INT16 scanNum - scan number in frame (LAC=0, GAC=0,1,..,4)     */

/* -------------------------------------------------------------- */
INT16 valid_instlm(INT16 mnftype, INT16 mnfnum, INT16 scanNum) {
    if (mnftype != 15)
        if (mnfnum != 2)
            return (1);
        else
            return 0;

    else
        return ( instlm_list[mnfnum - 1][scanNum]);
}


/* ----------------------------------------------------------------*/
/* getEngQual() - sets engineering quality flags based on analog   */
/*                instrument telemetry. Returns number of values   */
/*                that exceed red limits.                          */

/* ----------------------------------------------------------------*/
INT16 getEngQual(FLOAT32 ins_ana[], BYTE eng_qual[]) {
    INT16 i;
    INT16 cnt = 0;
    INT32 *qualp = (INT32 *) eng_qual;

    uint32_t mask = 1u << 31;

    *qualp = 0;

    for (i = 0; i < 32; i++) {

        if ((ins_ana[i] < instlm_limits[i][2]) ||
                (ins_ana[i] > instlm_limits[i][3])) {

            *qualp |= mask;

            /*
            printf("INSTLM limit exceeded for word %d (%f %f %f)\n",
                i,instlm_limits[i][2],ins_ana[i],instlm_limits[i][3]);
             */

            cnt++;
        }

        mask = mask >> 1;
    }

    if (endianess() == 1)
        swapc_bytes((char *) qualp, 4, 1);

    return (cnt);
}


