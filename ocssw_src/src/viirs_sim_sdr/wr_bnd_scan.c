#include "viirs_sim_sdr.h"
#include "l12_parms.h"
#include <stdlib.h>

int wr_bnd_scan(int iscn, out_rec_struc *out_rec)
/*-----------------------------------------------------------------------------
    Program:   wr_bnd_scan

    Description:  write a scan of band data to the VIIRS SDRs

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       iscn          I    scan to write
        out_rec_struc * out_rec   I/O  output record controls

    Modification history:

    W. Robinson, SAIC  15 Dec 2008  Original development
    W. Robinson, SAIC  15 Mar 2010  re-orient for scan based writing
    W. Robinson, SAIC  19 Nov 2010  collect, report the rad, refl/bbt 
                min, max, note q1 is primary quantity carried in out_rec->lt,
                while q2 is derived from it, currently q1 = radiance for 
                reflective and bbt for emissive bands, q2 is reflectance 
                for reflective and radiance for emissive

----------------------------------------------------------------------------*/ {
    static unsigned short *bnd_short;
    static float *bnd_float;
    float tmpval, tmpval2, scale, offset, mu0, lam_um;
    static int entry = 0;
    int ilin, nlin, loc, ipix, ibnd, start[2], count[2];
    h5io_str lcl_id;
    static float min_q1[MAX_BND], max_q1[MAX_BND],
            min_q2[MAX_BND], max_q2[MAX_BND];
    static int min_sq1[MAX_BND], max_sq1[MAX_BND],
            min_sq2[MAX_BND], max_sq2[MAX_BND];
    /*
     *  initialize the record of min, max values
     */
    if (iscn == 0) {
        for (ibnd = 0; ibnd < out_rec->nbnd; ibnd++) {
            min_q1[ibnd] = 999.;
            max_q1[ibnd] = -999.;
            min_q2[ibnd] = 999.;
            max_q2[ibnd] = -999.;
            min_sq1[ibnd] = 99999;
            max_sq1[ibnd] = -99999;
            min_sq2[ibnd] = 99999;
            max_sq2[ibnd] = -99999;
        }
    }
    /*
     *  prepare the output short and float buffers for final output scaled
     *  and unscaled data
     */
    nlin = out_rec->ndet_scan;
    if (entry == 0) {
        entry = 1;
        if ((bnd_short = (unsigned short *)
                malloc(nlin * out_rec->npix * sizeof ( unsigned short))) == NULL) {
            printf("%s, %d: Unable to allocate bnd_short storage\n",
                    __FILE__, __LINE__);
            return 1;
        }
        if ((bnd_float = (float *)
                malloc(nlin * out_rec->npix * sizeof ( float))) == NULL) {
            printf("%s, %d: Unable to allocate bnd_float storage\n",
                    __FILE__, __LINE__);
            return 1;
        }
    }

    ilin = out_rec->ndet_scan * iscn;
    start[0] = ilin;
    start[1] = 0;
    count[0] = out_rec->ndet_scan;
    count[1] = out_rec->npix;
    /*
     *  output all the band data
     */
    for (ibnd = 0; ibnd < out_rec->nbnd; ibnd++) {
        lam_um = (float) out_rec->lam_band[ibnd] / 1000.;
        /*
         *  either convert the primary radiance quantity: Lt / BBT to unsigned 
         *  short and output or output directly
         */
        lcl_id = out_rec->bnd_dat_id[1][ibnd]; /* set for BBT out */
        scale = out_rec->refl_scale[ibnd];
        offset = out_rec->refl_offset[ibnd];

        if (out_rec->meas_typ[ibnd] == 0) {
            lcl_id = out_rec->bnd_dat_id[0][ibnd]; /* or set for rad out */
            scale = out_rec->scale[ibnd];
            offset = out_rec->offset[ibnd];
        }
        if (out_rec->out_bnd_typ[ibnd] == 0) /* for scaled output of rad / BBT */ {
            for (ilin = 0; ilin < nlin; ilin++) {
                for (ipix = 0; ipix < out_rec->npix; ipix++) {
                    loc = ilin * out_rec->npix + ipix;
                    if (*(out_rec->bnd_q[ibnd] + loc) == 1)
                        *(bnd_short + loc) = ONBOARD_PT_UINT16_FILL;
                    else if (*(out_rec->bnd_q[ibnd] + loc) == 2)
                        *(bnd_short + loc) = MISS_UINT16_FILL;
                    else {
                        tmpval2 = *(out_rec->bnd_lt[ibnd] + loc);
                        tmpval = (tmpval2 - offset) / scale;
                        if ((tmpval < 0.) || (tmpval >= SOUB_UINT16_FILL)) {
                            *(bnd_short + loc) = SOUB_UINT16_FILL;
                        } else {
                            *(bnd_short + loc) = (unsigned short) tmpval;
                            if (tmpval < min_sq1[ibnd]) min_sq1[ibnd] = (int) tmpval;
                            if (tmpval > max_sq1[ibnd]) max_sq1[ibnd] = (int) tmpval;

                            if (tmpval2 < min_q1[ibnd]) min_q1[ibnd] = tmpval2;
                            if (tmpval2 > max_q1[ibnd]) max_q1[ibnd] = tmpval2;
                        }
                    }
                }
            }
            if (h5io_wr_ds_slice(&lcl_id, start, count, bnd_short) != 0) {
                printf(
                        "%s, %d: radiance (short) scan write failure, scan %d, band %d\n",
                        __FILE__, __LINE__, iscn, ibnd);
            }
        } else
            /*
             * for unscaled output of rad / BBT
             */ {
            for (ilin = 0, loc = 0; ilin < nlin; ilin++)
                for (ipix = 0; ipix < out_rec->npix; ipix++, loc++)
                    if (*(out_rec->bnd_q[ibnd] + loc) == 1)
                        *(bnd_float + loc) = ONBOARD_PT_FLOAT32_FILL;
                    else if (*(out_rec->bnd_q[ibnd] + loc) == 2)
                        *(bnd_float + loc) = MISS_FLOAT32_FILL;
                    else {
                        tmpval = *(out_rec->bnd_lt[ibnd] + loc);
                        *(bnd_float + loc) = tmpval;
                        if (tmpval < min_q1[ibnd]) min_q1[ibnd] = tmpval;
                        if (tmpval > max_q1[ibnd]) max_q1[ibnd] = tmpval;
                    }

            if (h5io_wr_ds_slice(&lcl_id, start, count, bnd_float) != 0) {
                printf(
                        "%s, %d: radiance (float) scan write failure, scan %d, band %d\n",
                        __FILE__, __LINE__, iscn, ibnd);
            }
        }
        /*
         *  end primary quantity output, start on the secondary
         *
         *  convert the radiance to reflectance or bbt to radiance
         */
        if (out_rec->meas_typ[ibnd] == 0) {
            /*
             *  for reflective bands, convert radiance to reflectance, scale and output
             *  Note assumption: all reflectances are scaled
             */
            lcl_id = out_rec->bnd_dat_id[1][ibnd];
            for (ilin = 0; ilin < nlin; ilin++) {
                for (ipix = 0; ipix < out_rec->npix; ipix++) {
                    loc = ilin * out_rec->npix + ipix;
                    if (*(out_rec->bnd_q[ibnd] + loc) == 1)
                        *(bnd_short + loc) = ONBOARD_PT_UINT16_FILL;
                    else if (*(out_rec->bnd_q[ibnd] + loc) == 2)
                        *(bnd_short + loc) = MISS_UINT16_FILL;
                    else {
                        mu0 = cos(*(out_rec->solz + loc) / RADEG);
                        tmpval2 = *(out_rec->bnd_lt[ibnd] + loc) * PI /
                                (out_rec->f0[ibnd] * mu0);
                        tmpval = (tmpval2 - out_rec->refl_offset[ibnd])
                                / out_rec->refl_scale[ibnd];
                        if ((tmpval < 0.) || (tmpval >= SOUB_UINT16_FILL))
                            *(bnd_short + loc) = SOUB_UINT16_FILL;
                        else {
                            *(bnd_short + loc) = (unsigned short) tmpval;
                            if (tmpval < min_sq2[ibnd]) min_sq2[ibnd] = (int) tmpval;
                            if (tmpval > max_sq2[ibnd]) max_sq2[ibnd] = (int) tmpval;

                            if (tmpval2 < min_q2[ibnd]) min_q2[ibnd] = tmpval2;
                            if (tmpval2 > max_q2[ibnd]) max_q2[ibnd] = tmpval2;
                        }
                    }
                }
            }
            if (h5io_wr_ds_slice(&lcl_id, start, count, bnd_short) != 0) {
                printf("%s, %d: reflectance scan write failure, scan %d, band %d\n",
                        __FILE__, __LINE__, iscn, ibnd);
            }
        } else {
            /*
             *  for emissive bands, convert bbt to radiance and output
             *  Note assumption: if primary quantity is scaled, then this
             *  will also be scaled
             */
            lcl_id = out_rec->bnd_dat_id[0][ibnd];
            for (ilin = 0; ilin < nlin; ilin++) {
                for (ipix = 0; ipix < out_rec->npix; ipix++) {
                    loc = ilin * out_rec->npix + ipix;
                    /*
                     *  make the radiance from brightness temp
                     */
                    if (*(out_rec->bnd_q[ibnd] + loc) == 1)
                        *(bnd_float + loc) = ONBOARD_PT_FLOAT32_FILL;
                    else if (*(out_rec->bnd_q[ibnd] + loc) == 2)
                        *(bnd_float + loc) = MISS_FLOAT32_FILL;
                    else {
                        tmpval = bbt_2_rad(*(out_rec->bnd_lt[ibnd] + loc), lam_um);
                        *(bnd_float + loc) = tmpval;
                        if (tmpval < min_q2[ibnd]) min_q2[ibnd] = tmpval;
                        if (tmpval > max_q2[ibnd]) max_q2[ibnd] = tmpval;
                    }
                    /*
                     *  convert if scaled output
                     */
                    if (out_rec->out_bnd_typ[ibnd] == 0) {
                        if (*(out_rec->bnd_q[ibnd] + loc) == 1)
                            *(bnd_short + loc) = ONBOARD_PT_UINT16_FILL;
                        else if (*(out_rec->bnd_q[ibnd] + loc) == 2)
                            *(bnd_short + loc) = MISS_UINT16_FILL;
                        else {
                            tmpval = (*(bnd_float + loc) - out_rec->offset[ibnd])
                                    / out_rec->scale[ibnd];
                            *(bnd_short + loc) = (unsigned short) tmpval;
                            if (tmpval < min_sq2[ibnd]) min_sq2[ibnd] = (int) tmpval;
                            if (tmpval > max_sq2[ibnd]) max_sq2[ibnd] = (int) tmpval;
                        }
                    }
                }
            }
            /*
             *  output the radiance
             */
            if (out_rec->out_bnd_typ[ibnd] == 0) {
                if (h5io_wr_ds_slice(&lcl_id, start, count, bnd_short) != 0) {
                    printf(
                            "%s, %d: emissive radiance scan write failure, scan %d, band %d\n",
                            __FILE__, __LINE__, iscn, ibnd);
                }
            } else {
                if (h5io_wr_ds_slice(&lcl_id, start, count, bnd_float) != 0) {
                    printf(
                            "%s, %d: emissive radiance scan write failure, scan %d, band %d\n",
                            __FILE__, __LINE__, iscn, ibnd);
                }
            }
        }
        /*
         *  write the scan of pixel level quality: QF1_VIIRSMBNDSDR
         */
        if (h5io_wr_ds_slice(&(out_rec->qual1_m_id[ibnd]), start, count,
                out_rec->qual1_m[ibnd]) != 0) {
            printf(
                    "%s, %d: scan write to QF1_VIIRSMBNDSDR failed. scan %d, band %d\n",
                    __FILE__, __LINE__, iscn, ibnd);
        }
    }
    /*
     *  remove any allocations for the transfer and report the min, max encountered
     */
    if (iscn == out_rec->nscan - 1) {
        free(bnd_short);
        free(bnd_float);

        printf("%s, %d: info - Data minimum, maximum\n\n", __FILE__, __LINE__);
        printf("       Radiance (W / m^2 sr um) BBT (deg K for 12-16)\n");
        printf("       unscaled             scaled\n");
        printf("Band   Minimum   Maximum    Minimum   Maximum\n");
        for (ibnd = 0; ibnd < out_rec->nbnd; ibnd++)
            printf("%4d  %8.4f  %8.4f    %6d    %6d\n", ibnd + 1, min_q1[ibnd],
                max_q1[ibnd], min_sq1[ibnd], max_sq1[ibnd]);

        printf(
                "\n       Reflectance (no units) Rad (W / m^2 sr um for 12-16)\n\n");
        printf("       unscaled             scaled\n");
        printf("Band   Minimum   Maximum    Minimum   Maximum\n");
        for (ibnd = 0; ibnd < out_rec->nbnd; ibnd++)
            printf("%4d  %8.6f  %8.6f    %6d    %6d\n", ibnd + 1, min_q2[ibnd],
                max_q2[ibnd], min_sq2[ibnd], max_sq2[ibnd]);
    }
    return 0;
}
