/* -------------------------------------------------------------------- */
/* getorbdata() - loads telemetry structure with raw GPS data           */
/* -------------------------------------------------------------------- */

#include <stdio.h>
#include "swl0_proto.h"
#include <timeutils.h>

INT32 getorbdata(swl0indx *indx, input_sType gps[]) {
    BYTE mnf[L0LEN];
    INT32 irec;
    FILE *fp;
    INT32 ngps = 0;
    INT16 ttag[4];
    FLOAT64 secs;
    int endian = 0;

    if (endianess() == 1)
        endian = 1;

    /*									*/
    /* Open file for reading						*/
    /*									*/
    if ((fp = fopen(indx->l0file, "r")) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: unable to open %s for reading\n",
                __FILE__, __LINE__, indx->l0file);
        return (1);
    }


    /*									*/
    /* Loop through L0 index, find frames with valid GPS data, and load	*/
    /* the gps data structure with converted telemetry. We'll take any  */
    /* frame with a valid time for the first frame, as we need that     */
    /* start-time to buffer the orbit vectors.        			*/
    /*									*/
    for (irec = 0; irec < indx->nrecs; irec++) {

        /*
         * GPS data is only stored in minor frame 1
         */
        if (
                ((ngps == 0) &&
                !(indx->rec[irec].tRanError ||
                indx->rec[irec].tSeqError ||
                indx->rec[irec].tDifError)) ||

                ((indx->rec[irec].mnfnum == 1) &&
                !(indx->rec[irec].scidError ||
                indx->rec[irec].sgaError ||
                indx->rec[irec].tRanError ||
                indx->rec[irec].tSeqError ||
                indx->rec[irec].tDifError ||
                indx->rec[irec].bitError))) {

            if (fseek(fp, sizeof (swl0hdr) + irec * (L0LEN), SEEK_SET) != 0) {
                fprintf(stderr,
                        "-E- %s line %d: error reading %s\n",
                        __FILE__, __LINE__, indx->l0file);
                return (1);
            }
            if (fread(mnf, L0LEN, 1, fp) != 1) {
                fprintf(stderr,
                        "-E- %s line %d: error reading %s\n",
                        __FILE__, __LINE__, indx->l0file);
                return (1);
            }

            gps[ngps].flag = 0;
            memcpy(gps[ngps].sc_id, &mnf[O_SCID], 4);
            if (endian)
                swapc_bytes((char *) gps[ngps].sc_id, 2, 2);
            memcpy(ttag, &mnf[O_TIME], 8);
            if (endian)
                swapc_bytes((char *) ttag, 2, 4);
            ttag2ydmsec(ttag,
                    &gps[ngps].iyear, &gps[ngps].iday, (INT32 *) & gps[ngps].msec);
            conv_soh_(&mnf[O_SOH], gps[ngps].sc_ana, gps[ngps].sc_dis);

            ngps++;

            if (ngps >= MAXFRAMES - 1) {
                fprintf(stderr,
                        "-W- %s line %d: limit on MAXFRAMES exceeded\n",
                        __FILE__, __LINE__);
                break;
            }

        }

    }

    /* add dummy record with time chosen to span data period */
    unix2yds(indx->timeLast, &gps[ngps].iyear, &gps[ngps].iday, &secs);
    gps[ngps].msec = (INT32) secs * 1000;
    memset(gps[ngps].sc_ana, 0, sizeof (gps[ngps].sc_ana));
    memset(gps[ngps].sc_dis, 0, sizeof (gps[ngps].sc_dis));
    ngps++;

    fclose(fp);

    return (ngps);

}




