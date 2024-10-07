#include "l1stat.h"
#include "l1stat_proto.h"
extern int32 stat_status;
extern char bad_stat_str[320];

void chk_tdiv(int32 sdfid, int16 dtynum, int32 nscans, thr_ctl_def thr_ctl)
/*******************************************************************

   chk_tdiv

   purpose: Check the tdi setting in some or all of the 8 bands
            for the GAC, LAC, SOL and LUN data which have only 1
            expected tdi per channel.  The SDS tdi is used for this

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
    int32 start[3], edge[3], count_bad, band, iscan, icode;
    int16 *tdi;
    float32 pct;
    int32 count_mul[6], i;
    float32 pct_mul[6];
    char str[12];

    /*
     *  the tdi for GAC, LAC, SOL and LUN should always be 0 -
     *  4:1 TDI or 1 contribution from each of the 4 detectors
     *
     *  for TDI and IGC, at least 8% have one of the 5 TDI settings of:
     *  0 - 4:1, 37 - use det #1 for all, 88 - use det #2 for all,
     *  82 - use det #3 for all, 133 - use det #4 for all
     */
    /*
     *  read in the SDS
     */
    if ((tdi = (int16 *) malloc(8 * nscans * sizeof ( int16))) == NULL) {
        printf("\n*****chk_tdiv: program error, unable to allocate tdi space\n");
        stat_status = stat_status | 1;
        return;
    }

    start[0] = start[1] = start[2] = 0;
    edge[0] = nscans;
    edge[1] = 8;
    edge[2] = 0;

    if (rdslice(sdfid, "tdi", start, edge, (void *) tdi) < 0) {
        stat_status = stat_status | 1;
        printf("\n*****chk_tdiv: program error, unable to read tdi\n");
        return;
    }
    if ((dtynum == GAC) || (dtynum == LAC) || (dtynum == SOL) ||
            (dtynum == LUN)) {
        /*
         *  loop through the 8 bands and check if they conform
         */
        printf("\n\nTDI value conformance check for GAC, LAC, SOL, LUN\n");


        /*  for lining it up...
          printf( "\n      Name code   #bad   %bad  error %
          printf( "\nTDIV_CHKX vvvv vvvvvv vvvvvvv vvvvvvv
          printf(   "--------- ---- ------ ------- -------
         */

        printf("\n     Name code   #bad   %%bad  error %%\n");
        printf("--------- ---- ------ ------- -------\n");

        for (band = 0; band < 8; band++) {
            if (thr_ctl.tdiv_chk_do[band] == 1) {
                count_bad = 0;
                for (iscan = 0; iscan < nscans; iscan++) {
                    if (*(tdi + band + 8 * iscan) != 0)
                        count_bad++;
                }
                /*
                 *  report the news for this item
                 */
                pct = (float32) count_bad / nscans * 100.;
                icode = 0;
                if (pct > thr_ctl.tdiv_chk_pct[band]) {
                    icode = 1;
                    if (thr_ctl.rpt_tdi_vchk == 1) {
                        stat_status = stat_status | 2;
                        sprintf(str, "TDIV_CHK%1d ", (band + 1));
                        if (strlen(bad_stat_str) <= 300)
                            strcat(bad_stat_str, str);
                    }
                }
                printf("TDIV_CHK%1d %4d %6d %7.2f %7.2f\n", (band + 1),
                        icode, count_bad, pct, thr_ctl.tdiv_chk_pct[band]);
            }
        }
    } else if ((dtynum == IGC) || (dtynum == TDI)) {
        /*
         *  loop through the 8 bands and check if they conform
         */
        printf("\n\nTDI value conformance check for IGC, TDI\n");

        /*  for lining it up...
          printf( "\n     Name code    %4:1 % det 1 % det 2 % det 3 % det 4 % other error %
          printf( "\nTDIV_CHKX vvvv vvvvvvv vvvvvvv vvvvvvv vvvvvvv vvvvvvv vvvvvvv vvvvvvv
          printf(   "--------- ---- ------- ------- ------- ------- ------- ------- -------
         */

        printf("\n     Name code    %%4:1 %% det 1 %% det 2 %% det 3 %% det 4 %% other error %%\n");
        printf("--------- ---- ------- ------- ------- ------- ------- ------- -------\n");

        for (band = 0; band < 8; band++) {
            if (thr_ctl.tdiv_chk_do[band] == 1) {
                for (i = 0; i < 6; i++)
                    count_mul[i] = 0;
                for (iscan = 0; iscan < nscans; iscan++) {
                    if (*(tdi + band + 8 * iscan) == 0)
                        count_mul[0]++;
                    else if (*(tdi + band + 8 * iscan) == 37)
                        count_mul[1]++;
                    else if (*(tdi + band + 8 * iscan) == 88)
                        count_mul[2]++;
                    else if (*(tdi + band + 8 * iscan) == 82)
                        count_mul[3]++;
                    else if (*(tdi + band + 8 * iscan) == 133)
                        count_mul[4]++;
                    else
                        count_mul[5]++;
                }
                /*
                 *  report the news for this item
                 */
                icode = 0;
                for (i = 0; i < 6; i++) {
                    pct_mul[i] = (float32) count_mul[i] / nscans * 100.;
                    if (i < 5) {
                        if (pct_mul[i] < 8.) {
                            icode = 1;
                            if (thr_ctl.rpt_tdi_vchk == 1) {
                                stat_status = stat_status | 2;
                                sprintf(str, "TDIV_CHK%1d ", (band + 1));
                                if (strlen(bad_stat_str) <= 300)
                                    strcat(bad_stat_str, str);
                            }
                        }
                    } else {
                        if (pct_mul[i] > thr_ctl.tdiv_chk_pct[band]) {
                            icode = 1;
                            if (thr_ctl.rpt_tdi_vchk == 1) {
                                stat_status = stat_status | 2;
                                sprintf(str, "TDIV_CHK%1d ", (band + 1));
                                if (strlen(bad_stat_str) <= 300)
                                    strcat(bad_stat_str, str);
                            }
                        }
                    }
                }
                printf("TDIV_CHK%1d %4d %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n",
                        (band + 1), icode, pct_mul[0], pct_mul[1], pct_mul[2], pct_mul[3],
                        pct_mul[4], pct_mul[5], thr_ctl.tdiv_chk_pct[band]);
            }
        }
    }
    printf("\n\n");
    /*
     *  that's all
     */
    return;
}
