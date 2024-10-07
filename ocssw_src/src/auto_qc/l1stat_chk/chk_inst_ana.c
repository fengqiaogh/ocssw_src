#include <stdlib.h>
#include "l1stat.h"
#include "l1stat_proto.h"
extern int32 stat_status;
extern char bad_stat_str[320];

void chk_inst_ana(int32 sdfid, int32 nscans, thr_ctl_def thr_ctl)
/*******************************************************************

   chk_inst_ana

   purpose: Check some or all of the 32 instrument analog values found
            in the level-1 dataset SDS: inst_ana.  Note that this is only 
            good for GAC due to the selection criteria for good inst_ana
            lines

   Returns type: void - none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             sdfid            I      SD interface ID
      int32             nscans           I      # lines in dataset
      thr_ctl_def       thr_ctl          I      structure with checking 
                                                thresholds

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       31-Oct-2001     Original development

 *******************************************************************/ {
    int32 start[3], edge[3], ct_hi, ct_lo, item, iscan, icode, index, good_scn,
            i5, tot5, bad5;
    int32 ix275[] = {0, 1, 1, 0, 1}, ix147[] = {1, 0, 1, 1, 0},
    ix403[] = {1, 1, 0, 1, 1};
    int16 *sc_id;
    float32 pct_hi, pct_lo, pct, *inst_ana, *ana_copy, pct_loc[7];
    char str[12];
    int val_comp(const void *, const void *);
    /*
     *  read in the inst_ana SDS
     */
    if (
            (inst_ana = (float32 *) malloc(40 * nscans * sizeof ( float32)))
            == NULL ||
            (ana_copy = (float32 *) malloc(nscans * sizeof ( float32))) == NULL ||
            (sc_id = (int16 *) malloc(nscans * sizeof ( int16))) == NULL) {
        printf("\n*****chk_inst_ana: program error, unable to allocate inst_ana space\n");
        stat_status = stat_status | 1;
        return;
    }

    start[0] = start[1] = start[2] = 0;
    edge[0] = nscans;
    edge[1] = 40;
    edge[2] = 0;

    if (rdslice(sdfid, "inst_ana", start, edge, (void *) inst_ana) < 0) {
        stat_status = stat_status | 1;
        printf("\n*****chk_inst_ana: program error, unable to read inst_ana\n");
        return;
    }

    start[0] = start[1] = start[2] = 0;
    edge[0] = nscans;
    edge[1] = 1;
    edge[2] = 0;
    if (rdslice(sdfid, "sc_id", start, edge, (void *) sc_id) < 0) {
        stat_status = stat_status | 1;
        printf("\n*****chk_inst_ana: program error, unable to read sc_id\n");
        return;
    }
    /*
     *  isolate only the lines that have good telemetry.  5 lines will all have 
     *  same sc_id.  the good tlm lines depend on the sc_id value.  we'll
     *  select only the good values and re-compose inst_ana in-place
     */
    tot5 = 0;
    bad5 = 0;
    for (iscan = 0, good_scn = 0; iscan < nscans; iscan = iscan + 5) {
        if (*(sc_id + iscan) == 275 && *(sc_id + iscan + 1) == 275 &&
                *(sc_id + iscan + 2) == 275 && *(sc_id + iscan + 3) == 275 &&
                *(sc_id + iscan + 4) == 275) {
            for (i5 = 0; i5 < 5; i5++) {
                if (ix275[i5] == 1) {
                    for (item = 0; item < 40; item++) {
                        *(inst_ana + item + 40 * good_scn) =
                                *(inst_ana + item + 40 * (iscan + i5));
                    }
                    good_scn++;
                }
            }
        } else if (*(sc_id + iscan) == 147 && *(sc_id + iscan + 1) == 147 &&
                *(sc_id + iscan + 2) == 147 && *(sc_id + iscan + 3) == 147 &&
                *(sc_id + iscan + 4) == 147) {
            for (i5 = 0; i5 < 5; i5++) {
                if (ix147[i5] == 1) {
                    for (item = 0; item < 40; item++) {
                        *(inst_ana + item + 40 * good_scn) =
                                *(inst_ana + item + 40 * (iscan + i5));
                    }
                    good_scn++;
                }
            }
        } else if (*(sc_id + iscan) == 403 && *(sc_id + iscan + 1) == 403 &&
                *(sc_id + iscan + 2) == 403 && *(sc_id + iscan + 3) == 403 &&
                *(sc_id + iscan + 4) == 403) {
            for (i5 = 0; i5 < 5; i5++) {
                if (ix403[i5] == 1) {
                    for (item = 0; item < 40; item++) {
                        *(inst_ana + item + 40 * good_scn) =
                                *(inst_ana + item + 40 * (iscan + i5));
                    }
                    good_scn++;
                }
            }
        } else {
            bad5++;
        }
        tot5++;
    }
    /*
     *  only work with data that has good scan lines remaining
     */
    if (good_scn > 0) {
        /*
         *  loop through the 32 telemetry items and check if any are 
         *  outside the thresholds
         */
        printf("\n\nInstrument analog item bound threshold check\n");


        /*  for lining it up...
          printf( "\n      Name     code %outside  lo thresh  hi thresh  error % (    lo %     hi %)
                     INST_ANAvv vvvvvvvv vvvvvvvv vvvvvvvvvv vvvvvvvvvv vvvvvvvv (vvvvvvvv vvvvvvvv)
          printf(   "---------- -------- -------- ---------- ---------- --------  -------- --------
         */


        printf("\n      Name code %%outside  lo thresh   hi thresh error %% (   lo %%    hi %%)\n");
        printf("---------- ---- ------- ----------- ----------- -------  ------- -------\n");

        for (item = 0; item < 32; item++) {
            if (thr_ctl.inst_ana_do[item] == 1) {
                ct_hi = 0;
                ct_lo = 0;
                for (iscan = 0; iscan < good_scn; iscan++) {
                    if (*(inst_ana + item + 40 * iscan) > thr_ctl.inst_ana_hi[item])
                        ct_hi++;
                    if (*(inst_ana + item + 40 * iscan) < thr_ctl.inst_ana_lo[item])
                        ct_lo++;
                }
                /*
                 *  report the news for this item
                 */
                pct_hi = (float32) ct_hi / good_scn * 100.;
                pct_lo = (float32) ct_lo / good_scn * 100.;
                pct = pct_hi + pct_lo;
                icode = 0;
                if (pct > thr_ctl.inst_ana_pct[item]) {
                    icode = 1;
                    if (thr_ctl.rpt_inst_ana == 1) {
                        stat_status = stat_status | 2;
                        sprintf(str, "INST_ANA%02d ", (item + 1));
                        if (strlen(bad_stat_str) <= 300)
                            strcat(bad_stat_str, str);
                    }
                }
                printf("INST_ANA%02d %4d %7.2f %11.3f %11.3f %7.2f (%7.2f %7.2f)\n",
                        (item + 1), icode, pct,
                        thr_ctl.inst_ana_lo[item], thr_ctl.inst_ana_hi[item],
                        thr_ctl.inst_ana_pct[item], pct_lo, pct_hi);
            }
        }
        printf("\n\n");
        /*
         *  list out the values that mark % of data as a guide
         */
        printf("\nValues where following %% of data are at\n\n");
        printf("\n      Name  Low (0%%)        1%%       10%%       50%%       90%%      99%% High(100%%)\n");
        printf("---------- --------- --------- --------- --------- --------- --------- ---------\n");
        for (item = 0; item < 32; item++) {
            if (thr_ctl.inst_ana_do[item] == 1) {
                for (iscan = 0; iscan < good_scn; iscan++) {
                    *(ana_copy + iscan) = *(inst_ana + item + 40 * iscan);
                }
                /*
                 *  sort the data and find where different % of data are
                 */
                qsort((void *) ana_copy, (size_t) good_scn, (size_t) 4, val_comp);
                /*
                 * record the data values where the specified % of values are
                 */
                pct_loc[0] = ana_copy[0];
                index = (int32) (0.01 * good_scn);
                pct_loc[1] = ana_copy[index];
                index = (int32) (0.1 * good_scn);
                pct_loc[2] = ana_copy[index];
                index = (int32) (0.5 * good_scn);
                pct_loc[3] = ana_copy[index];
                index = (int32) (0.9 * good_scn);
                pct_loc[4] = ana_copy[index];
                index = (int32) (0.99 * good_scn);
                pct_loc[5] = ana_copy[index];
                pct_loc[6] = ana_copy[good_scn - 1];

                printf("INST-ANA%02d %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
                        (item + 1), pct_loc[0], pct_loc[1], pct_loc[2], pct_loc[3],
                        pct_loc[4], pct_loc[5], pct_loc[6]);
            }
        }
    }
    /*
     *  report the good scan groups and bad scan groups
     */
    icode = 0;
    if (bad5 > 0) {
        icode = 1;
        if (thr_ctl.rpt_inst_ana == 1) {
            stat_status = stat_status | 2;
            if (strlen(bad_stat_str) <= 300)
                strcat(bad_stat_str, "BAD_SCID ");
        }
    }
    pct = 100. * (float) bad5 / tot5;
    printf("\n\nAmount of 5-line GAC groups with bad sc_id:\n\n");
    printf("    Name  code  %% bad groups  # bad groups  Total # groups\n");
    printf("--------  ----  ------------  ------------  --------------\n");
    printf("BAD_SCID  %4d  %12.3f  %12d  %14d\n", icode, pct, bad5, tot5);
    /*
     *  that's all
     */
    return;
}
/*  compare routine  */

/*
static int val_comp( const void *v1, const void *v2 )
 */
int val_comp(const void *v1, const void *v2) {
    float *lv1 = (float *) v1;
    float *lv2 = (float *) v2;
    if (*lv1 == *lv2) return 0;
    return ( *lv1 < *lv2) ? -1 : 1;
}
