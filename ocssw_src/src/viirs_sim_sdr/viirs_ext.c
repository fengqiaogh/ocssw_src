#include "viirs_sim_sdr.h"
#include <string.h>
#include <stdlib.h>

int viirs_ext(ctl_struc *ctl, in_rec_struc *in_rec)
/*-----------------------------------------------------------------------------
    Program:   viirs_ext.c

    Description:  Add the VIIRS electrical crosstalk artifact to the simulated 
      dn field (derived from the radiance field)

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        ctl_struc *  ctl        I    input controls
        in_rec_struc * in_rec    I/O  controls for input record reading

    Modification history:

    W. Robinson, SAIC  18 Jan 2011  Original development

----------------------------------------------------------------------------*/ {
    static int ext_init = 0;
    static float *dn_tmp[MAX_BND];
    static double *coeff;
    float dn_mod;
    int nsgain, nsbnd, nsdet, nrgain, nrbnd, nrdet, nrpix, nrlin, nspix;
    int irbnd, irlin, irpix, ispix, isdet, isbnd, isgain, irgain, irdet;
    int reg_off, nval;
    int16_t nsfgain, nsfbnd, nsfdet, nrfgain, nrfbnd, nrfdet;
    h5io_str ext_fid;
    char tstr[500];
    /* registration offset for the re-ordered sender bands:
       (where is band M? relative to M1?)
                             M1  M2  M3  M4   M5   M6   M7   I1   I2   */
    static int snd_off[] = {0, -3, -9, -6, -18, -21, -15, -12, -14};
    /*
     *  The storage of the info for the gain state in the ext corff array is 
     *  in the order ( low, high ) while in the in_rec->gain_bit it is
     *  ( high, low ), so this will translate that assignment
     */
    static int gain_lut[] = {0, 1};

    nrlin = in_rec->ndet_scan;
    nrpix = in_rec->npix;
    nspix = nrpix;
    nrbnd = N_VNIR_BND;
    nrdet = NDET;
    nrgain = 2;
    nsgain = 2;
    nsbnd = 7;
    nsdet = NDET;
    /*
     *  we could put in code to de-allocate the dn_tmp and coeff
     *  if in_rec->lat = NULL as scan_cvt.c does.  You'd call this, when
     *  appropriate in fin_sdr.c
     */
    if (in_rec->lat == NULL) {
        free(coeff);
        for (irbnd = 0; irbnd < nrbnd; irbnd++)
            free(dn_tmp[irbnd]);
        return 0;
    }
    /*
     *  do initialization
     */
    if (ext_init == 0) {
        ext_init = 1;

        printf("%s, %d: Entering routine viirs_ext\n", __FILE__, __LINE__);
        printf("ext_coef = %s\n", ctl->ext_coeff_file);
        /*
         *  allocate temporary, un-modified dn storage for the dn scan
         */
        for (irbnd = 0; irbnd < nrbnd; irbnd++)
            if ((dn_tmp[irbnd] = (float *) malloc(nrlin * nrpix * sizeof ( float)))
                    == NULL) {
                printf("%s, %d: Unable to allocate dn tmp storage\n", __FILE__,
                        __LINE__);
                return 1;
            }
        /*
         *  Get the coefficients for the crostalk
         */
        nval = nsgain * nrgain * nsbnd * nsdet * nrbnd * nrdet;
        if ((coeff = (double *) calloc(nval, sizeof ( double))) == NULL) {
            printf("%s, %d: Unable to allocate electronic crosstalk storage\n",
                    __FILE__, __LINE__);
            return 1;
        }
        if (strcmp(ctl->ext_coeff_file, "Unspecified") == 0) {
            /*  default coeffs of 0 = no crosstalk, are used if ext coeff
             *  file is unspecified
             *  perhaps some code will be needed here, so have space here
             */
        } else {
            /*  open file  */
            if (h5io_openr(ctl->ext_coeff_file, 0, &ext_fid) != 0) {
                printf("%s, %d - could not open HDF 5 Electrical crosstalk file: %s\n",
                        __FILE__, __LINE__, ctl->ext_coeff_file);
                return 1;
            }
            /*  read info on storage size, note that most sizes are a fixed value */
            strcpy(tstr, "number of sender gains");
            if (h5io_rd_attr(&ext_fid, tstr, &nsfgain) != 0) {
                printf("%s, %d - failed to read Electrical Xtalk attrib %s\n",
                        __FILE__, __LINE__, tstr);
                return 1;
            }
            strcpy(tstr, "number of receiver gains");
            if (h5io_rd_attr(&ext_fid, tstr, &nrfgain) != 0) {
                printf("%s, %d - failed to read Electrical Xtalk attrib %s\n",
                        __FILE__, __LINE__, tstr);
                return 1;
            }
            strcpy(tstr, "number of sender bands");
            if (h5io_rd_attr(&ext_fid, tstr, &nsfbnd) != 0) {
                printf("%s, %d - failed to read Electrical Xtalk attrib %s\n",
                        __FILE__, __LINE__, tstr);
                return 1;
            }
            strcpy(tstr, "number of sender detectors");
            if (h5io_rd_attr(&ext_fid, tstr, &nsfdet) != 0) {
                printf("%s, %d - failed to read Electrical Xtalk attrib %s\n",
                        __FILE__, __LINE__, tstr);
                return 1;
            }
            strcpy(tstr, "number of receiver bands");
            if (h5io_rd_attr(&ext_fid, tstr, &nrfbnd) != 0) {
                printf("%s, %d - failed to read Electrical Xtalk attrib %s\n",
                        __FILE__, __LINE__, tstr);
                return 1;
            }
            strcpy(tstr, "number of receiver detectors");
            if (h5io_rd_attr(&ext_fid, tstr, &nrfdet) != 0) {
                printf("%s, %d - failed to read Electrical Xtalk attrib %s\n",
                        __FILE__, __LINE__, tstr);
                return 1;
            }
            /*
             *  check dimensions out to match
             */
            if (nsfgain != nsgain) {
                printf(
                        "%s, %d - Electrical crosstalk # sender gains: %d != expected #: %d\n",
                        __FILE__, __LINE__, nsfgain, nsgain);
                return 1;
            }
            if (nrfgain != nrgain) {
                printf(
                        "%s, %d - Electrical crosstalk # receiver gains: %d != expected #: %d\n",
                        __FILE__, __LINE__, nrfgain, nrgain);
                return 1;
            }
            if (nsfbnd != nsbnd) {
                printf(
                        "%s, %d - Electrical crosstalk # sender bands: %d != expected #: %d\n",
                        __FILE__, __LINE__, nsfbnd, nsbnd);
                return 1;
            }
            if (nrfbnd != nrbnd) {
                printf(
                        "%s, %d - Electrical crosstalk # receiver bands: %d != expected #: %d\n",
                        __FILE__, __LINE__, nrfbnd, nrbnd);
                return 1;
            }
            if (nsfdet != nsdet) {
                printf(
                        "%s, %d - Electrical crosstalk # sender detectors: %d != expected #: %d\n",
                        __FILE__, __LINE__, nsfdet, nsdet);
                return 1;
            }
            if (nrfdet != nrdet) {
                printf(
                        "%s, %d - Electrical crosstalk # receiver detectors: %d != expected #: %d\n",
                        __FILE__, __LINE__, nrfdet, nrdet);
                return 1;
            }
            /*
             *  read data to the coefficient storage
             */
            if (h5io_grab_ds(&ext_fid, "electronic crosstalk coefficients",
                    (void *) coeff) != 0) {
                printf("%s, %d: Unable to read electronic crosstalk coefficients\n",
                        __FILE__, __LINE__);
                return 1;
            }

            /*  close the coeff file */
            if (h5io_close(&ext_fid) != 0) {
                printf(
                        "%s, %d: Unable to close electronic crosstalk coefficient file\n",
                        __FILE__, __LINE__);
                return 1;
            }
        }
    }
    /*
     *  fill temporary storage with the scan dn
     *  Note that all the dn data is being used - no track margin cutoff
     */
    nval = nrpix * nrlin * sizeof ( float);
    for (irbnd = 0; irbnd < nrbnd; irbnd++)
        memcpy(dn_tmp[irbnd], in_rec->dn[irbnd], nval);
    /*
     *  loop through the lines, pixels, and bands in the scan
     */
    for (irlin = 0; irlin < nrlin; irlin++) {
        irdet = irlin - in_rec->margin[0];
        if (irdet < 0) irdet = 0;
        if (irdet >= NDET) irdet = NDET - 1;
        for (irpix = 0; irpix < in_rec->npix; irpix++) {
            for (irbnd = 0; irbnd < nrbnd; irbnd++) {
                /*
                 *  zero out the dn modification
                 *  and accumulate all the dn modifications from the senders
                 */
                dn_mod = 0.;
                for (isbnd = 0; isbnd < nsbnd; isbnd++) {
                    reg_off = snd_off[isbnd] - snd_off[irbnd];
                    ispix = irpix + reg_off;
                    if ((ispix >= 0) && (ispix < in_rec->npix)) {
                        for (isdet = 0; isdet < nsdet; isdet++) {
                            isgain = gain_lut[ *(in_rec-> gain_bit[isbnd] +
                                    ispix + nspix * isdet) ];
                            irgain = gain_lut[ *(in_rec-> gain_bit[irbnd] +
                                    irpix + nspix * irlin) ];

                            dn_mod += *(dn_tmp[isbnd] + ispix + isdet * nspix) *
                                    *(coeff + isgain + nsgain * (irgain + nrgain *
                                    (isbnd + nsbnd * (isdet + nsdet *
                                    (irbnd + nrbnd * irdet)))));
                        }
                    }
                }
                /*
                 *  add the modified dn to the dn for the receiver pix, lin, band
                 */
                *(in_rec->dn[irbnd] + irpix + nrpix * irlin) =
                        *(dn_tmp[irbnd] + irpix + nrpix * irlin) + dn_mod;
            }
        }
    }
    /*
     *  scan is modified, return
     */
    return 0;
}
