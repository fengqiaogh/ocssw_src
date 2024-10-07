#include "viirs_sim_sdr.h"
#include <math.h>

int viirs_decal(ctl_struc *ctl, in_rec_struc *in_rec)
/*-----------------------------------------------------------------------------
    Program:   viirs_decal.c

    Description:  convert VIIRS radiances for a scan into dn values
     The dn can be either dark subtracted (dn) or not (DN) based on the
     cal_dark_sub value (if 1 [default] make it dark sutract, 0 no subtract)

    Returns type: int - 0 if good

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        ctl_struc *  ctl        I    input controls
        in_rec_struc * in_rec  I/O   controls for input record reading

    Modification history:

    W. Robinson, SAIC  16 Sep 2010  Original development

----------------------------------------------------------------------------*/ {
    static int first_time = 0, rvs_note = 0;
    static float dn_max = 4095.; /* max representable dn in 12 bits */
    static float dn_trans_pct = 0.83; /* general % of max dn when gain 
                                        transition occurs, NICST memo
                                        NICST_MEMO_10_012_*.doc  */
    int ilin, idet, ismp, ipix, ibnd, igain, iham, dn_sat;
    int npix, nbnd;
    float c0, c1, c2, a0, a1, a2, dark, lt, aoi, rvs_val, dn, big_dn;
    float min_rvs[MAX_BND], max_rvs[MAX_BND];
    static vir_gain_struc gain;
    static vir_rvs_struc rvs;
    /*
     *  set up the gain and rvs strucs from the source files
     *  only do the vis NIR (7 bands)
     */
    iham = in_rec->ham_side;
    nbnd = in_rec->nbnd;
    if (first_time == 0) {
        if (vset_cal_gain(ctl->count_decal_gain_file, &gain) != 0)
            return 1;
        /*
         *  check gain struct so it satisfys the needs
         */
        if (gain.nham != 2) {
            printf("%s, %d: # of HAM sides in gain struct (%d) != 2\n", __FILE__,
                    __LINE__, gain.nham);
            printf("\tgain file: %s\n", ctl->count_decal_gain_file);
            return 1;
        }
        if (gain.ndet != NDET) {
            printf(
                    "%s, %d: # of detectors in gain struct (%d) != the VIIRS # (%d)\n",
                    __FILE__, __LINE__, gain.ndet, NDET);
            printf("\tgain file: %s\n", ctl->count_decal_gain_file);
            return 1;
        }
        if (gain.ngain != 2) {
            printf("%s, %d: # of gain states in gain struct (%d) != 2\n", __FILE__,
                    __LINE__, gain.nham);
            printf("\tgain file: %s\n", ctl->count_decal_gain_file);
            return 1;
        }
        if (gain.nbnd < N_VNIR_BND) {
            printf(
                    "%s, %d: # of bands in gain struct (%d) not >= # VIS / NIR bands (%d)\n",
                    __FILE__, __LINE__, gain.nbnd, N_VNIR_BND);
            printf("\tgain file: %s\n", ctl->count_decal_gain_file);
            return 1;
        }

        if (vset_cal_rvs(ctl->count_decal_rvs_file, &rvs) != 0)
            return 1;
        /*
         *  check rvs struct so it satisfys the needs
         */
        if (rvs.nham != 2) {
            printf("%s, %d: # of HAM sides in RVS struct (%d) != 2\n", __FILE__,
                    __LINE__, rvs.nham);
            printf("\trvs file: %s\n", ctl->count_decal_rvs_file);
            return 1;
        }
        if (rvs.ndet != NDET) {
            printf(
                    "%s, %d: # of detectors in RVS struct (%d) != the VIIRS # (%d)\n",
                    __FILE__, __LINE__, rvs.ndet, NDET);
            printf("\trvs file: %s\n", ctl->count_decal_rvs_file);
            return 1;
        }
        if (rvs.nbnd < N_VNIR_BND) {
            printf("%s, %d: # of bands in RVS struct (%d) not >= # needed # (%d)\n",
                    __FILE__, __LINE__, rvs.nbnd, nbnd);
            printf("\trvs file: %s\n", ctl->count_decal_rvs_file);
            return 1;
        }

        for (ibnd = 0; ibnd < N_VNIR_BND; ibnd++) {
            min_rvs[ibnd] = 999.;
            max_rvs[ibnd] = -999.;
        }

        first_time = 1;
    }
    /*
     *  loop through the pixels, lines, and bands to get the lt
     */
    npix = in_rec->npix;
    for (ipix = 0; ipix < npix; ipix++) {
        /* get sample */
        ismp = ipix - in_rec->margin[1];

        /*  get the aoi from the sample */
        vir_xf_scan((float) ismp, VIR_SCAN_UASMP, VIR_SCAN_AOI, &aoi);
        if ((rvs_note == 0) && ((ipix % 200) == 0)) {
            float scn_ang;
            vir_xf_scan((float) ismp, VIR_SCAN_UASMP, VIR_SCAN_ANG, &scn_ang);
            printf("%s, %d: info, sample: %d, scan angle: %f, aoi: %f\n",
                    __FILE__, __LINE__, ipix, scn_ang, aoi);
        }

        for (ilin = 0; ilin < in_rec->ndet_scan; ilin++) {
            /*  get the detector from the line */
            idet = ilin - in_rec->margin[0];
            if (idet < 0) idet = 0;
            if (idet >= NDET) idet = NDET - 1;

            for (ibnd = 0; ibnd < N_VNIR_BND; ibnd++) {
                /*
                 *  if the lt is not flagged, get it and de-cal
                 */
                if (*(in_rec->bnd_q[ibnd] + ipix + npix * ilin) == 0) {
                    /*
                     *  Assume a high gain range and check later
                     */
                    lt = *(in_rec->bnd_lt[ibnd] + ipix + npix * ilin);
                    igain = 1;
                    dn_sat = 0;
                    dark = *(gain.dark + ibnd + gain.nbnd * (igain + gain.ngain *
                            (idet + gain.ndet * iham)));
                    /*
                     *  get the gain, rvs coeffs
                     */
                    c0 = *(gain.c0 + ibnd + gain.nbnd * (igain + gain.ngain *
                            (idet + gain.ndet * iham)));
                    c1 = *(gain.c1 + ibnd + gain.nbnd * (igain + gain.ngain *
                            (idet + gain.ndet * iham)));
                    c2 = *(gain.c2 + ibnd + gain.nbnd * (igain + gain.ngain *
                            (idet + gain.ndet * iham)));

                    a0 = *(rvs.a0 + ibnd + rvs.nbnd * (idet + rvs.ndet * iham));
                    a1 = *(rvs.a1 + ibnd + rvs.nbnd * (idet + rvs.ndet * iham));
                    a2 = *(rvs.a2 + ibnd + rvs.nbnd * (idet + rvs.ndet * iham));

                    rvs_val = a0 + a1 * aoi + a2 * aoi * aoi;
                    if (rvs_note == 0) {
                        /*
                         *  record the min and max rvs for each band
                         */
                        if (rvs_val < min_rvs[ibnd])
                            min_rvs[ibnd] = rvs_val;
                        if (rvs_val > max_rvs[ibnd])
                            max_rvs[ibnd] = rvs_val;
                    }
                    /*
                     *  apply the gain, rvs to de cal
                     */
                    if (c2 == 0) /* just linear fit */ {
                        dn = (lt * rvs_val - c0) / c1;
                    } else /* a quadratic fit */ {
                        dn = (-c1 + sqrt(c1 * c1 - 4. * c2 * (c0 - lt * rvs_val)))
                                / (2. * c2);
                    }
                    /*
                     *  if a dual gain band, see if dn is above the transition
                     */
                    if (gain.gain_ct[ibnd] == 2) {
                        if (dn > (dn_max * dn_trans_pct)) {
                            igain = 0;
                            dark = *(gain.dark + ibnd + gain.nbnd * (igain + gain.ngain *
                                    (idet + gain.ndet * iham)));
                            /*
                             *  get the low gain coefficients
                             */
                            c0 = *(gain.c0 + ibnd + gain.nbnd * (igain + gain.ngain *
                                    (idet + gain.ndet * iham)));
                            c1 = *(gain.c1 + ibnd + gain.nbnd * (igain + gain.ngain *
                                    (idet + gain.ndet * iham)));
                            c2 = *(gain.c2 + ibnd + gain.nbnd * (igain + gain.ngain *
                                    (idet + gain.ndet * iham)));
                            /*
                             *  apply the gain, rvs to de cal
                             */
                            if (c2 == 0) /* just linear fit */ {
                                dn = (lt * rvs_val - c0) / c1;
                            } else /* a quadratic fit */ {
                                dn = (-c1 + sqrt(c1 * c1 - 4. * c2 * (c0 - lt * rvs_val)))
                                        / (2. * c2);
                            }
                            big_dn = dn + dark;
                            /*
                             *  if it saturates at low gain too, note it
                             */
                            if (big_dn > dn_max) {
                                dn_sat = 1;
                                big_dn = dn_max;
                            }
                        } else /* below, dual gain transition, set the raw dn */ {
                            big_dn = dn + dark;
                        }
                        dn = big_dn - dark;
                    } else /* single gain treatment  */ {
                        big_dn = dn + dark;
                        if (big_dn > dn_max) {
                            dn_sat = 1;
                            big_dn = dn_max;
                        }
                        dn = big_dn - dark;
                    }
                    /*
                     *  if DN is wanted use big_dn, else, remove the dark from big_dn
                     */
                    if (ctl->count_dark_opt == 0)
                        dn = big_dn - dark;
                    else
                        dn = big_dn;
                } else /* bad lt val, leave dn 0  */ {
                    dn = 0.;
                    dn_sat = -1;
                }
                /*
                 *  place the count... in the structure
                 */
                *(in_rec->dn[ibnd] + ipix + npix * ilin) = dn;
                *(in_rec->dn_sat[ibnd] + ipix + npix * ilin) = dn_sat;
                *(in_rec->gain_bit[ibnd] + ipix + npix * ilin) = igain;
            }
        }
    }
    /*
     *  report the rvs min, max
     */
    if (rvs_note == 0) {
        rvs_note = 1;
        printf("\n\n%s, %d:info\n", __FILE__, __LINE__);
        printf("\nDe-calibration RVS min and max / band:\n");
        printf("(file: %s\n", ctl->count_decal_rvs_file);
        printf("\nM band    Min RVS    Max RVS\n");
        for (ibnd = 0; ibnd < N_VNIR_BND; ibnd++)
            printf("   %3d  %9.6f  %9.6f\n", ibnd + 1,
                min_rvs[ibnd], max_rvs[ibnd]);
        printf("----------------------------------------\n\n");
    }
    /* 
     *  end loops and done
     */
    return 0;
}
