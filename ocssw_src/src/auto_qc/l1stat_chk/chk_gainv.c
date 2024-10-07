#include "l1stat.h"
#include "l1stat_proto.h"
extern int32 stat_status;
extern char bad_stat_str[320];

void chk_gainv(int32 sdfid, int16 dtynum, int32 nscans, thr_ctl_def thr_ctl)
/*******************************************************************

   chk_gainv

   purpose: Check the gain setting in some or all of the 8 bands
            for the GAC, LAC, SOL and KUN data which have only 1
            expected gain per channel.  The SDS gain is used for this

   Returns type: void - none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             sdfid            I      SD interface ID
      int16             dtynum           I      data type #
      int32             nscans           I      # lines in dataset
      thr_ctl_def       thr_ctl          I      structure with checking 
                                                thresholds

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       2-Nov-2001     Original development

 *******************************************************************/ {
    int32 start[3], edge[3], count_bad, band, gptr, iscan, icode;
    int16 *gain;
    /*  gval_dtyp will give the gain settings for GAC, LAC, SOL and LUN */
    int32 gval_dtyp[4][8] = {
        { 0, 0, 0, 0, 0, 0, 0, 0},
        { 0, 0, 0, 0, 0, 0, 0, 0},
        { 1, 0, 1, 1, 1, 1, 1, 1},
        { 3, 1, 1, 1, 1, 1, 1, 1}
    };
    float32 pct;
    char str[12];
    /*
     *  according to the data type, point to the correct expected gains
     */
    if (dtynum == GAC) {
        gptr = 0;
    } else if (dtynum == LAC) {
        gptr = 1;
    } else if (dtynum == SOL) {
        gptr = 2;
    } else if (dtynum == LUN) {
        gptr = 3;
    } else
        return;
    /*
     *  read in the SDS
     */
    if ((gain = (int16 *) malloc(8 * nscans * sizeof ( int16))) == NULL) {
        printf("\n*****chk_gainv: program error, unable to allocate gain space\n");
        stat_status = stat_status | 1;
        return;
    }

    start[0] = start[1] = start[2] = 0;
    edge[0] = nscans;
    edge[1] = 8;
    edge[2] = 0;

    if (rdslice(sdfid, "gain", start, edge, (void *) gain) < 0) {
        stat_status = stat_status | 1;
        printf("\n*****chk_gainv: program error, unable to read gain\n");
        return;
    }
    /*
     *  loop through the 8 bands and check if they conform
     */
    printf("\n\nGain value conformance check\n");


    /*  for lining it up...
      printf( "\n      Name code   #bad   %bad  error %
      printf( "\nGAINV_CHKX vvvv vvvvvv vvvvvvv vvvvvvv
      printf(   "---------- ---- ------ ------- -------
     */

    printf("\n      Name code   #bad   %%bad  error %%\n");
    printf("---------- ---- ------ ------- -------\n");

    for (band = 0; band < 8; band++) {
        if (thr_ctl.gainv_chk_do[band] == 1) {
            count_bad = 0;
            for (iscan = 0; iscan < nscans; iscan++) {
                if (*(gain + band + 8 * iscan) != gval_dtyp[gptr][band])
                    count_bad++;
            }
            /*
             *  report the news for this item
             */
            pct = (float32) count_bad / nscans * 100.;
            icode = 0;
            if (pct > thr_ctl.gainv_chk_pct[band]) {
                icode = 1;
                if (thr_ctl.rpt_gainv_chk == 1) {
                    stat_status = stat_status | 2;
                    sprintf(str, "GAINV_CHK%1d ", (band + 1));
                    if (strlen(bad_stat_str) <= 300)
                        strcat(bad_stat_str, str);
                }
            }
            printf("GAINV_CHK%1d %4d %6d %7.2f %7.2f\n", (band + 1),
                    icode, count_bad, pct, thr_ctl.gainv_chk_pct[band]);
        }
    }
    printf("\n\n");
    /*
     *  that's all
     */
    return;
}
