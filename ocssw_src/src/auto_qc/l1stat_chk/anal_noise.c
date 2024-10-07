#include <string.h>
#include <libgen.h>
#include "l1stat.h"
#include "l1stat_proto.h"
extern int32 stat_status;

void anal_noise(int32 rec, int32 nrec, int32 nscans, int32 nsamp,
        int32 nbands, int16 *i16buf, int *spike_cnt, float *line_sd,
        int32_t *cnt_coin_jmp, int32_t *jmp_hist)
/*******************************************************************

   anal_noise

   purpose: check a buffer of level-1a data line by line and band by band
            to see if there is evidence of spike noise or an encrypted line
            only information is generated.  later, this will get combined
            with checks of other SDS data to evaluate the total noise

   Returns type: void - none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32             rec              I      line # of this current 
                                                batch of level-1a lines
      int32             nrec             I      # of lines in this array
      int32             nscans           I      total # scan lines
      int32             nsamp            I      # pixels per line
      int32             nbands           I      # of bands
      int16 *           i16buf           I      line by pixel by band
                                                array of raw counts
      int *             spike_cnt        O      line by band array of
                                                spike count 
      float *           line_sd          O      line by band array of
                                                std deviation for line
      int32_t *            cnt_coin_jmp     O      size 8 count of co-incidences
                                                (# times 1-8 bands for a 
                                                pixel are on at once)
      int32_t *            jmp_hist         O      size 1024 histogram
                                                of jump size

      access spike_cnt and line_sd to get band b, line l with
      index = b * nrec + l

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       31-Jul-1998     Original development

 *******************************************************************/ {
    int rec_num, irec, ibnd, ipix;
    short *spike_jmp, pixm1, pix, pixp1, isum;
    double apx, sum, sumsq, pix_jmp, sd;
    float jfact = 1.3; /* 1.3 for HRPT, 5 for GAC?  It seems that GAC
                          and HRPT still get 'spikes' for good, cloud pix,
                          so usefullness may be limited at this time  */
    int32_t neighbor_diff;

    if ((spike_jmp = (short *) malloc(nbands * nsamp * sizeof ( short)))
            == NULL) {
        printf("\n*****anal_noise: program error, unable to allocate space\n");
        stat_status = stat_status | 1;
    } else {
        /*
         *  loop through the lines in this batch
         */
        rec_num = rec;

        for (irec = 0; irec < nrec; irec++, rec_num++) {
            for (ibnd = 0; ibnd < nbands; ibnd++) {
                /*
                 *  Locate the spikes.  The general idea is to see if a pixel 
                 *  is wildly different from it's next store neighbors
                 *  a band-dependent minimum excursion will be imposed (to avoid
                 *  classifying noise jumps as a spike) and a maximun
                 *  neighbor change will also be imposed (to avois electronics
                 *  overshoot points)
                 */
                sum = 0.;
                sumsq = 0.;
                for (ipix = 1; ipix < (nsamp - 1); ipix++) {
                    pixm1 = *(i16buf + ibnd + nbands * ((ipix - 1) + nsamp * irec));
                    pix = *(i16buf + ibnd + nbands * (ipix + nsamp * irec));
                    pixp1 = *(i16buf + ibnd + nbands * ((ipix + 1) + nsamp * irec));
                    neighbor_diff = abs(pixm1 - pixp1);
                    pix_jmp = fabs(pix - (pixm1 + pixp1) / 2.);
                    /*
                     * do not do points that have a large neighbor jump or a small
                     * pixel jump
                     */

                    *(spike_jmp + ibnd + nbands * ipix) = 0;
                    if (pix_jmp > 15 && neighbor_diff < 200) {
                        if (pix_jmp > jfact * neighbor_diff) {
                            (*(spike_cnt + ibnd + nbands * rec_num))++;
                            *(spike_jmp + ibnd + nbands * ipix) = pix_jmp;
                            /*
                             printf( "rec: %d, pix: %d, bnd: %d, jump: %f\n",
                                rec_num, ipix, ibnd, pix_jmp );
                             */
                            if (pix_jmp < 1024)
                                (*(jmp_hist + (int) pix_jmp))++;
                        }
                    }
                    /*
                     *  accumulate the sum and sum square for the statistics
                     */
                    apx = (double) pix;
                    sum += apx;
                    sumsq += apx * apx;
                } /* end pixel loop */
                /*
                 *  compute the standard deviation for the line and band
                 */
                sd = sumsq / nsamp - sum * sum / (nsamp * nsamp);
                sd = (sd > 0) ? sqrt(sd) : 0.;
                *(line_sd + ibnd + nbands * (rec_num)) = (float) sd;
                /*
                 printf( "rec: %d, bnd: %d, sd: %f\n", irec, ibnd, sd );
                 */
            } /* end band loop */
            /*
             *  find frequency of pops occuring at same pixel
             */
            for (ipix = 1; ipix < nsamp - 1; ipix++) {
                isum = 0;
                for (ibnd = 0; ibnd < nbands; ibnd++) {
                    if (*(spike_jmp + ibnd + nbands * ipix) > 0)
                        isum++;
                }
                if (isum > 0) {
                    /*
                     printf( "rec: %d, pix: %d, co-incident jumps: %d\n",
                           irec, ipix, isum );
                     */
                    (* (cnt_coin_jmp + isum - 1))++;
                }
            }
        } /* end record loop */
        /*
         *  remove space for spike jumps
         */
        free(spike_jmp);
    }
}
