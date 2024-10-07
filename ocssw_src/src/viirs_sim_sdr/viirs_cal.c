#include "viirs_sim_sdr.h"
#include <math.h>

int viirs_cal(ctl_struc *ctl, in_rec_struc *in_rec)
/*-----------------------------------------------------------------------------
    Program:   viirs_cal.c

    Description:  convert dn values (either dark subtracted (dn) or 
      not (DN) based on the cal_dark_sub value) and gain indicator 
      into VIIRS radiances

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
    int ilin, idet, ipix, ibnd, igain, iham;
    int npix, nbnd;
    float c0, c1, c2, a0, a1, a2, dark, dn, lt, aoi, rvs_val;
    float min_rvs[MAX_BND], max_rvs[MAX_BND];
    static vir_gain_struc gain;
    static vir_rvs_struc rvs;
    /*
     *  set up the gain and rvs strucs from the source files
     *  for now, only do vis, NIR (7 bands)
     */
    iham = in_rec->ham_side;
    nbnd = in_rec->nbnd;
    if (first_time == 0) {
        if (vset_cal_gain(ctl->count_cal_gain_file, &gain) != 0)
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
            printf("%s, %d: # of bands in gain struct (%d) not >= # needed # (%d)\n",
                    __FILE__, __LINE__, gain.nbnd, nbnd);
            printf("\tgain file: %s\n", ctl->count_decal_gain_file);
            return 1;
        }

        if (vset_cal_rvs(ctl->count_decal_rvs_file, &rvs) != 0)
            return 1;
        /*
         *  check rvs struct so it satisfys the needs
         */
        if (rvs.nham != 2) {
            printf("%s, %d: # of HAM sides in gain struct (%d) != 2\n", __FILE__,
                    __LINE__, rvs.nham);
            printf("\trvs file: %s\n", ctl->count_decal_rvs_file);
            return 1;
        }
        if (rvs.ndet != NDET) {
            printf(
                    "%s, %d: # of detectors in gain struct (%d) != the VIIRS # (%d)\n",
                    __FILE__, __LINE__, rvs.ndet, NDET);
            printf("\trvs file: %s\n", ctl->count_decal_rvs_file);
            return 1;
        }
        if (rvs.nbnd < N_VNIR_BND) {
            printf("%s, %d: # of bands in gain struct (%d) not >= # needed # (%d)\n",
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
        /*  get the aoi from the sample */
        vir_xf_scan(ipix, VIR_SCAN_UASMP, VIR_SCAN_AOI, &aoi);

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
                     *  get the gain and dn from the input structure
                     */
                    igain = *(in_rec->gain_bit[ibnd] + ipix + npix * ilin);
                    dn = *(in_rec->dn[ibnd] + ipix + npix * ilin);
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
                    if ((rvs_val > 2.) || (rvs_val < .5)) {
                        printf("%s, %d: rvs val %f, odd at pix: %d, lin: %d, bnd: %d\n",
                                __FILE__, __LINE__, rvs_val, ipix, ilin, ibnd);
                    }
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
                     *  remove dark count if required, apply the gain and rvs to calibrate
                     */
                    if (ctl->count_dark_opt == 1) {
                        dark = *(gain.dark + ibnd + gain.nbnd * (igain + gain.ngain *
                                (idet + gain.ndet * iham)));
                        dn = dn - dark;
                    }
                    lt = (c0 + c1 * dn + c2 * dn * dn) / rvs_val;
                    /*
                     *  place the radiance back in the structure
                     */
                    *(in_rec->bnd_lt[ibnd] + ipix + npix * ilin) = lt;
                }
                /* else bad lt val, leave lt alone  */
            }
        }
    }
    /*
     *  report the rvs min, max
     */
    if (rvs_note == 0) {
        rvs_note = 1;
        printf("\n\n%s, %d: info:\n", __FILE__, __LINE__);
        printf("\nCalibration RVS min and max / band:\n");
        printf("(file: %s\n", ctl->count_cal_rvs_file);
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
