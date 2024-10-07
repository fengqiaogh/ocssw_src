#include "viirs_sim_sdr.h"

int fin_sdr(ctl_struc *ctl, in_rec_struc *in_rec, out_rec_struc *out_rec)
/*-----------------------------------------------------------------------------
    Program:   fin_sdr

    Description:  finish up on the simulation processing: close input
      and output dataset IDs

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        ctl_struc * ctl         I    processing controls
        in_rec_struc * in_rec   I    input record controls
        out_rec_struc * out_rec   I/O  output record controls

    Modification history:

    W. Robinson, SAIC  21 Oct 2008  Original development

----------------------------------------------------------------------------*/ {
    int iid, jid, isdr, ibnd;
    /*
     *  Input dataset controls
     *  geolocation
     */
    for (iid = 0; iid < 6; iid++) {
        if (h5io_close(&(in_rec->geo_dat_id[ iid ])) != 0) {
            printf("%s: Unable to close input geoloc dataset # %d\n",
                    __FILE__, iid);
            return 1;
        }
    }

    if (h5io_close(&(in_rec->geo_fid)) != 0) {
        printf("%s: Unable to close input geoloc file\n", __FILE__);
        return 1;
    }
    /*
     *  compute the granule latitude, longitude limits
     *  this could be placed in the dataset, if needed, in the future
     */
    printf("latitude limits, S: %f, N: %f\n", in_rec->ll_lims[0],
            in_rec->ll_lims[1]);
    printf("(-) longitude limits, W: %f, E: %f\n", in_rec->ll_lims[2],
            in_rec->ll_lims[3]);
    printf("(+) longitude limits, W: %f, E: %f\n", in_rec->ll_lims[4],
            in_rec->ll_lims[5]);
    /* w_lim = -180.;
    e_lim = 180.;
     we'll have to make some kind of actual W, E limits */
    /*
     *  Output dataset controls
     *  geolocation and then band data ids
     */
    for (iid = 0; iid < 6; iid++) {
        if (h5io_close(&(out_rec->geo_dat_id[ iid ])) != 0) {
            printf("%s: Unable to close output geoloc dataset id # %d\n",
                    __FILE__, iid);
            return 1;
        }
    }
    for (isdr = 0; isdr < out_rec->nbnd; isdr++) {
        for (iid = 0; iid < 2; iid++) {
            if (h5io_close(&(out_rec->bnd_dat_id[iid][isdr])) != 0) {
                printf("%s: Unable to close output band file # %d, dataset # %d\n",
                        __FILE__, isdr, iid);
                return 1;
            }
        }
        /*
         *  QF1_VIIRSMBANDSDR dataset id close
         */
        if (h5io_close(&(out_rec->qual1_m_id[isdr])) != 0) {
            printf(
                    "%s, %d: Failed closing output QF1_VIIRSMBANDSDR dataset file # %d\n",
                    __FILE__, __LINE__, isdr);
            return 1;
        }
    }
    /*  open groups  */
    for (isdr = 0; isdr < out_rec->nbnd + 1; isdr++) {
        for (iid = 0; iid < 2; iid++) {
            jid = 1 - iid;
            if (h5io_close(&(out_rec->sdr_dat_gid[ jid ][isdr])) != 0) {
                printf("%s: Unable to close output sdr group id # %d, sdr # %d\n",
                        __FILE__, jid, isdr);
                return 1;
            }
        }
    }
    /* open files  */
    for (isdr = 0; isdr < out_rec->nbnd + 1; isdr++) {
        if (h5io_close(&(out_rec->sdr_fid[isdr])) != 0) {
            printf("%s: Unable to close output sdr file for sdr# %d\n",
                    __FILE__, isdr);
            return 1;
        }
    }
    /* close the L2 file if used */
    if (ctl->l2_use == 1)
        if (closeL2(&in_rec->l2_str, 0) != 0) {
            printf("%s: Unable to close input L2 file: %s\n",
                    __FILE__, in_rec->l2_str.filename);
            return 1;
        }
    /*
     *  free data storage, note that the in_rec->lat of NULL signals frees 
     *  in scan_cvt
     */
    free(in_rec->lat);
    in_rec->lat = NULL;
    if (scan_cvt(in_rec, out_rec) != 0)
        return 1;
    free(in_rec->lon);
    free(in_rec->sena);
    free(in_rec->senz);
    free(in_rec->sola);
    free(in_rec->solz);
    for (ibnd = 0; ibnd < out_rec->nbnd; ibnd++) {
        free(in_rec->bnd_lt[ibnd]);
        free(in_rec->bnd_q[ibnd]);
        free(out_rec->qual1_m[ibnd]);
        free(in_rec->gain_bit[ibnd]);
        free(in_rec->dn_sat[ibnd]);
        if (ctl->count_cal_opt != 0) {
            free(in_rec->dn[ibnd]);
        }
    }
    /*
     *  clean up for artifact use
     *  free the electronic crosstalk internal storage
     */
    if (ctl->ext_opt != 0)
        if (viirs_ext(ctl, in_rec) != 0)
            return 1;
    /*  for the noise and stray light generation  */
    if (ctl->noise_mode == 1)
        viirs_noise(ctl, in_rec, 1);
    if (ctl->stray_opt == 1)
        viirs_straylt(ctl, in_rec, 1);

    return 0;
}
