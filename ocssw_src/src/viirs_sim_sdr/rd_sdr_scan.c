#include "viirs_sim_sdr.h"
#include <stdlib.h>

int rd_sdr_scan(int iscn, ctl_struc *ctl, sdr_info_struc *sdr_info,
        in_rec_struc *in_rec)
/*-----------------------------------------------------------------------------
    Program:   rd_sdr_scan.c

    Description:  get all scan oriented data from simulated data sources

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       iscn          I    scan to read
        ctl_struc *  ctl        I    processing controls
        sdr_info_struc * sdr_info I  general sdr information
        in_rec_struc * in_rec   I/O  input record controls

    Modification history:

    W. Robinson, SAIC  21 Oct 2008  Original development
    W. Robinson, SAIC  17 Feb 2010  update readL2 call seq for added last
       (and undefined) parameter
    W. Robinson, SAIC  12 Mar 2010  adapt to process whole scans

----------------------------------------------------------------------------*/ {
    int ipix, ibnd, ilin, idet;
    static float *rhos = NULL;
    float rad;
    /*
     *  place the ham side in the in_rec
     */
    in_rec->ham_side = *(sdr_info->ham_side + iscn);
    /*
     *  read the geoloc SDR
     */
    if (rd_geo_scan(iscn, sdr_info, in_rec) != 0) return 1;
    /*
     *  Proceed to reading up the band files
     */
    switch (ctl->l2_use) {
    case 0:
        /*
         printf( "%s-temp: no Lt data file option, making constant line\n",
           __FILE__ );
         */
        if (gen_const_rad_scn(sdr_info, 0, in_rec) != 0) {
            printf("%s, %d: call to gen_const_rad_scn failed on line %d\n",
                    __FILE__, __LINE__, iscn);
            return 1;
        }
        break;
    case 1:
        /*
         *  get a scan of reflectance if needed
         */
        if (ctl->rhos_use == 1)
            if (rd_rhos_scan(ctl->rhos_file, in_rec->npix, iscn,
                    in_rec->ndet_scan, &rhos) != 0) return 1;
        /*
         *  get the L2 line and move the radiances into the input record
         *  convert from milliwatt cm^-2... to w m^-2...
         *  (some other l2 info may also go there in the future)
         */
        for (idet = 0; idet < in_rec->ndet_scan; idet++) {
            ilin = idet + iscn * in_rec->ndet_scan;
            if (readL2(&(in_rec->l2_str), 0, ilin, -1, NULL) != 0) {
                printf("%s, %d: call to readL2 failed on line %d\n",
                        __FILE__, __LINE__, ilin);
                return 1;
            }
            for (ibnd = 0; ibnd < N_VNIR_BND; ibnd++) {
                for (ipix = 0; ipix < in_rec->npix; ipix++) {
                    /*  remember to check for a bad value and note it */
                    rad = in_rec->l2_str.l2_data[ibnd][ipix];
                    if ((rad == BAD_FLT) ||
                            (rad < -0.5)) {
                        *(in_rec->bnd_q[ibnd] + idet * in_rec->npix + ipix) = 2;
                        *(in_rec->bnd_lt[ibnd] + idet * in_rec->npix + ipix) = 0.;
                    } else {
                        *(in_rec->bnd_lt[ibnd] + idet * in_rec->npix + ipix) =
                                rad * RAD_CGS_2_MKS;
                        *(in_rec->bnd_q[ibnd] + idet * in_rec->npix + ipix) = 0;
                    }
                }
            }
            /*
             *  If requested, convert reflectance to Lt using L2 fields for 
             *  that line and fill in the SDR radiances
             */
            if (ctl->rhos_use == 1) {
                if (rhos_to_lt(ctl->rhos_opt, &rhos, in_rec, idet, sdr_info) != 0)
                    return 1;
            }
            /*
             *  finally, if the rest of the bands are done, put in typical values
             */
            if (in_rec->nbnd > N_VNIR_BND) {
                if (gen_const_rad_scn(sdr_info, N_VNIR_BND, in_rec) != 0) {
                    printf("%s, %d: call to gen_const_rad_scn failed on line %d\n",
                            __FILE__, __LINE__, iscn);
                    return 1;
                }
            }
        }
        break;
    default:
        printf("%s: Bad Lt input source, Exiting\n",
                __FILE__);
        return 1;
        break;
    }
    return 0;
}
