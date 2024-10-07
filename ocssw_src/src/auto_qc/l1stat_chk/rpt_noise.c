#include <string.h>
#include <libgen.h>
#include "l1stat.h"
#include "l1stat_proto.h"
extern int32 stat_status;
extern char bad_stat_str[320];

#define DETAIL_PRT 0  /*  1 to print detail, 0 to not print detail */

void rpt_noise(int32 sdfid, int16 dtynum, int32 nscans, int32 nsamp,
        int *spike_cnt, float *line_sd, float pct_noise_thresh,
        float pct_encrypt_thresh)
/*******************************************************************

   rpt_noise

   purpose: get information from 5 other SDSes that flag noise in
            the data, combine it with the spike count and std deviation
            computed from the chk_count routine and decide if this dataset
            has too much noise to just pass

   Returns type: void - none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             sdfid            I      SD interface ID
      int16             dtynum           I      data type number (not used now)
      int32             nscans           I      # lines in dataset
      int32             nsamp            I      # samples / line
      int *             spike_cnt        I      size nscans by 8 array of 
                                                spike counts found
      float *           line_sd          I      size nscans by 8 array of 
                                                std deviations found 
      float             pct_noise_thresh I      reporting limit for % of
                                                noisey lines
      float             pct_encrypt_thresh I    reporting limit for % of
                                                encrypted lines

      access spike_cnt and line_sd to get band b, line l with
      index = b * nrec + l

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       31-Jul-1998     Original development
      W. Robinson       24-Oct-2001     have input noise report thresholds

 *******************************************************************/ {
    int irec, ibnd, icode;
    int16 *s_satp, *s_zerop, *dark_rest, *start_syn, *stop_syn;
    int32 *msec;
    int16 code_start_syn[8] = {1023, 0, 1023, 0, 1023, 0, 1023, 0}, c_start_syn[8];
    int16 code_stop_syn[8] = {1023, 1023, 0, 0, 1023, 1023, 0, 0}, c_stop_syn[8];
    float s_zerop_thresh[8] = {0., 0., 0., 0., 0.0016, 0.008, 0.016, 0.008};
    /* this is fraction of # pixels with 0 counts and corresponds
       to 0, 0, 0, 0, 2, 10, 20, 10 for 1285 pixels (in a HRPT) */
    int32 start[3], edge[3], encrypt_cnt, ds_noise_cnt, lin_prob_cnt,
            s_satp_bad, s_zerop_bad, dark_rest_bad, start_syn_bad,
            stop_syn_bad, sum_spikes;
    /* normal noisey line threshold = 80, encrypted: 70.  now have as inputs */
    float pct_noise, pct_encrypt;
    double pct_s_satp = 0., pct_s_zerop = 0., pct_dark_rest = 0.,
            pct_start_syn = 0., pct_stop_syn = 0.;
    char str[12];

    /*
     *  set up the start and stop sync to check using the code and data type
     */
    if ((dtynum == IGC) || (dtynum == TDI) || (dtynum == SOL)) {
        for (ibnd = 0; ibnd < 7; ibnd++) {
            c_start_syn[ibnd] = 1023 - code_start_syn[ibnd];
            c_stop_syn[ibnd] = 1023 - code_stop_syn[ibnd];
        }
    } else {
        for (ibnd = 0; ibnd < 7; ibnd++) {
            c_start_syn[ibnd] = code_start_syn[ibnd];
            c_stop_syn[ibnd] = code_stop_syn[ibnd];
        }
    }

    /*
     *  allocate the space needed for the 5 arrays
     */
    if ((s_satp = (int16 *) malloc(8 * nscans * sizeof ( int16)))
            == NULL) {
        printf("\n*****rpt_noise: program error, unable to allocate s_satp space\n");
        stat_status = stat_status | 1;
        return;
    }
    if ((s_zerop = (int16 *) malloc(8 * nscans * sizeof ( int16)))
            == NULL) {
        printf("\n*****rpt_noise: program error, unable to allocate s_zerop space\n");
        stat_status = stat_status | 1;
        return;
    }
    if ((dark_rest = (int16 *) malloc(8 * nscans * sizeof ( int16)))
            == NULL) {
        printf("\n*****rpt_noise: program error, unable to allocate dark_rest space\n");
        stat_status = stat_status | 1;
        return;
    }
    if ((start_syn = (int16 *) malloc(8 * nscans * sizeof ( int16)))
            == NULL) {
        printf("\n*****rpt_noise: program error, unable to allocate start_syn space\n");
        stat_status = stat_status | 1;
        return;
    }
    if ((stop_syn = (int16 *) malloc(8 * nscans * sizeof ( int16)))
            == NULL) {
        printf("\n*****rpt_noise: program error, unable to allocate stop_syn space\n");
        stat_status = stat_status | 1;
        return;
    }

    /*
     *  for other tests, get msec also (for line id)
     */
    if ((msec = (int32 *) malloc(nscans * sizeof ( int32)))
            == NULL) {
        printf("\n*****rpt_noise: program error, unable to allocate msec space\n");
        stat_status = stat_status | 1;
        return;
    }

    /*
     *  read the arrays from the L1 dataset
     */
    start[0] = start[1] = start[2] = 0;
    edge[0] = nscans;
    edge[1] = 8;
    edge[2] = 0;

    if (rdslice(sdfid, "s_satp", start, edge, (void *) s_satp) < 0) {
        stat_status = stat_status | 1;
        printf("\n*****rpt_noise: program error, unable to read s_satp\n");
        return;
    }

    if (rdslice(sdfid, "s_zerop", start, edge, (void *) s_zerop) < 0) {
        stat_status = stat_status | 1;
        printf("\n*****rpt_noise: program error, unable to read s_zerop\n");
        return;
    }

    if (rdslice(sdfid, "dark_rest", start, edge, (void *) dark_rest) < 0) {
        stat_status = stat_status | 1;
        printf("\n*****rpt_noise: program error, unable to read dark_rest\n");
        return;
    }

    if (rdslice(sdfid, "start_syn", start, edge, (void *) start_syn) < 0) {
        stat_status = stat_status | 1;
        printf("\n*****rpt_noise: program error, unable to read start_syn\n");
        return;
    }

    if (rdslice(sdfid, "stop_syn", start, edge, (void *) stop_syn) < 0) {
        stat_status = stat_status | 1;
        printf("\n*****rpt_noise: program error, unable to read stop_syn\n");
        return;
    }

    edge[0] = nscans;
    edge[1] = 0;
    edge[2] = 0;
    if (rdslice(sdfid, "msec", start, edge, (void *) msec) < 0) {
        stat_status = stat_status | 1;
        printf("\n*****rpt_noise: program error, unable to read stop_syn\n");
        return;
    }

    /*
     *  initialize some values for the entire dataset
     */
    encrypt_cnt = 0;
    ds_noise_cnt = 0;

    /*
     *  print heading for line-by-line report (if enabled)
     */

    if (DETAIL_PRT) {
        printf("\n\nLine-by-line report of noise and noise type found\n");
        printf("line  # bands   # bands  # bands w  # bad strt  # bad stop  std dev  # spikes   msec    line\n");
        printf(" #   saturated    zero   high dark  sync wds/7  sync wds/7   band 1  all bnds   tag     type\n");
        printf("---- ---------  -------  ---------  ----------  ----------  -------  --------  -------- ----\n");
    }

    /*
     *  loop through the records and get onfo on each one
     */
    for (irec = 0; irec < nscans; irec++) {
        /*
         *  initialize counts for each line 
         */
        lin_prob_cnt = 0;
        s_satp_bad = 0;
        s_zerop_bad = 0;
        dark_rest_bad = 0;
        start_syn_bad = 0;
        stop_syn_bad = 0;
        sum_spikes = 0;

        /*
         *  loop through each band or sync word (only 7 checked for sync words )
         */
        for (ibnd = 0; ibnd < 8; ibnd++) {
            if (*(s_satp + ibnd + 8 * irec) != 0) {
                s_satp_bad++;
                lin_prob_cnt++;
            }
            if ((float) *(s_zerop + ibnd + 8 * irec) / nsamp
                    > *(s_zerop_thresh + ibnd)) {
                s_zerop_bad++;
                lin_prob_cnt++;
            }
            if (*(dark_rest + ibnd + 8 * irec) > 30) {
                dark_rest_bad++;
                lin_prob_cnt++;
            }
            if (ibnd < 7) /* only check the 1st 7 words for the sync */ {
                if (*(start_syn + ibnd + 8 * irec) != *(c_start_syn + ibnd)) {
                    start_syn_bad++;
                    lin_prob_cnt++;
                }
                if (*(stop_syn + ibnd + 8 * irec) != *(c_stop_syn + ibnd)) {
                    stop_syn_bad++;
                    lin_prob_cnt++;
                }
            }
            sum_spikes += *(spike_cnt + ibnd + 8 * irec);
        } /* end band (word) loop, report the line results */

        if (DETAIL_PRT) {
            printf("%4d%10d%9d%11d%12d%12d%9.2f%10d%9d", irec, s_satp_bad,
                    s_zerop_bad, dark_rest_bad, start_syn_bad, stop_syn_bad,
                    *(line_sd + 8 * irec), sum_spikes, *(msec + irec));
        }

        /*
         *  need to add the line type at the end
         */
        if (*(line_sd + 8 * irec) > 220.) {
            if (DETAIL_PRT) {
                printf("  **** Encrypted line\n");
            }
            encrypt_cnt++;
            ds_noise_cnt++;
        } else if (lin_prob_cnt > 0) {
            if (DETAIL_PRT) {
                printf("    ** Noisey line\n");
            }
            ds_noise_cnt++;
        } else {
            if (DETAIL_PRT) {
                printf("       Good line\n");
            }
        }

        /*
         *  count the occurences of each of the 5 SDS violations
         *  (make into percentages at the end)
         */
        if (s_satp_bad > 0) pct_s_satp++;
        if (s_zerop_bad > 0) pct_s_zerop++;
        if (dark_rest_bad > 0) pct_dark_rest++;
        if (start_syn_bad > 0) pct_start_syn++;
        if (stop_syn_bad > 0) pct_stop_syn++;

    } /* end record loop, time for the final reconing */

    /*
     *  report percentages of the 5 sdses violated
     */
    printf("\n\n#Noise report summary\n\n");
    printf("# --- Percentage bad lines ---\n");
    printf("#        bands     bands    bands w    bad strt    bad stop  \n");
    printf("#     saturated    zero   high dark  sync wds/7  sync wds/7  \n");
    printf("#     ---------  -------  ---------  ----------  ----------  \n");
    printf("#    %10.2f%9.2f%11.2f%12.2f%12.2f\n\n",
            (100. * pct_s_satp / nscans), (100. * pct_s_zerop / nscans),
            (100. * pct_dark_rest / nscans), (100. * pct_start_syn / nscans),
            (100. * pct_stop_syn / nscans));

    printf("\n\n#total # lines: %9d  # with noise: %9d  # encrypted: %d\n",
            nscans, ds_noise_cnt, encrypt_cnt);
    pct_noise = 100. * (float) ds_noise_cnt / nscans;
    pct_encrypt = 100. * (float) encrypt_cnt / nscans;
    printf
            ("                                        %9.2f %%       %9.2f %%\n\n\n",
            pct_noise, pct_encrypt);
    /*
     *  add the encryption and noise reports also
     */
    printf("\n\n#Mnemonic   Code    %% lines   %% threshold\n");
    printf(" --------   ----    -------   -----------\n");
    if ((dtynum == GAC) || (dtynum == LAC) || (dtynum == HRPT)) {
        icode = 0;
        if (pct_encrypt > pct_encrypt_thresh) {
            printf("#*** dataset failed due to > %7.2f %% encryption noise\n",
                    pct_encrypt_thresh);
            stat_status = stat_status | 2;
            icode = 1;
            sprintf(str, "ENCRYPT ");
            if (strlen(bad_stat_str) <= 300)
                strcat(bad_stat_str, str);
        }
        printf("ENCRYPT       %2d    %7.2f       %7.2f\n", icode, pct_encrypt,
                pct_encrypt_thresh);
        icode = 0;
        if (pct_noise > pct_noise_thresh) {
            printf("#*** dataset failed due to > %7.2f %% general noise\n",
                    pct_noise_thresh);
            stat_status = stat_status | 2;
            icode = 1;
            sprintf(str, "NOISE ");
            if (strlen(bad_stat_str) <= 300)
                strcat(bad_stat_str, str);
        }
        printf("NOISE         %2d    %7.2f       %7.2f\n\n", icode,
                pct_noise, pct_noise_thresh);
    } else {
        printf("ENCRYPT      N/A    %7.2f       %7.2f\n", pct_encrypt,
                pct_encrypt_thresh);
        printf("NOISE        N/A    %7.2f       %7.2f\n\n",
                pct_noise, pct_noise_thresh);
    }
    /*
     *  make a 1 line mnemonic version for easy auto-pull-out
     */
    printf("\n\n#1 line version for auto-identification\n");
    printf("\n#Mnemonic Code noise  encrypt  sat   zero  dark  start  stop\n");
    printf("#----------------------------------------------------------------------------\n");
    if ((dtynum == GAC) || (dtynum == LAC) || (dtynum == HRPT)) {
        printf("L1NOISE    %2d%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f\n",
                icode, pct_noise, pct_encrypt, (100. * pct_s_satp / nscans),
                (100. * pct_s_zerop / nscans), (100. * pct_dark_rest / nscans),
                (100. * pct_start_syn / nscans), (100. * pct_stop_syn / nscans));
    } else {
        printf("L1NOISE    N/A%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f%7.2f\n",
                pct_noise, pct_encrypt, (100. * pct_s_satp / nscans),
                (100. * pct_s_zerop / nscans), (100. * pct_dark_rest / nscans),
                (100. * pct_start_syn / nscans), (100. * pct_stop_syn / nscans));
    }

    free(s_satp);
    free(s_zerop);
    free(dark_rest);
    free(start_syn);
    free(stop_syn);

    /*
     *  that's all
     */
    return;
}
