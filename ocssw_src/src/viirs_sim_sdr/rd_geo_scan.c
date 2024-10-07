#include "viirs_sim_sdr.h"

int rd_geo_scan(int iscn, sdr_info_struc *sdr_info, in_rec_struc *in_rec)
/*-----------------------------------------------------------------------------
    Program:   rd_geo_scan.c

    Description:  get all scan oriented data from simulated geolocation dataset

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       iscn          I    scan to read
        sdr_info_struc * sdr_info I  general sdr information
        in_rec_struc * in_rec   I/O  input record controls

    Modification history:

    W. Robinson, SAIC  21 Oct 2008  Original development
    W. Robinson, SAIC  17 Mar 2010  adapt for scan oriented use

----------------------------------------------------------------------------*/ {
    int start[2], count[2], ipix, idet, st_line, nlin, loc;
    float lon;
    /*
     *  Start with the latitude scan
     */
    st_line = iscn * in_rec->ndet_scan;
    nlin = in_rec->ndet_scan;
    start[0] = st_line;
    start[1] = 0;
    count[0] = nlin;
    count[1] = in_rec->npix;
    if (h5io_rd_ds_slice(&(in_rec->geo_dat_id[0]), start, count,
            (void *) (in_rec->lat)) != 0) {
        printf("%s, %d: lat input failed on scan: %d\n", __FILE__, __LINE__,
                iscn);
        return 1;
    }
    /*  longitude */
    if (h5io_rd_ds_slice(&(in_rec->geo_dat_id[1]), start, count,
            (void *) (in_rec->lon)) != 0) {
        printf("%s, %d: lon input failed on scan: %d\n", __FILE__, __LINE__,
                iscn);
        return 1;
    }
    /*
     *  update the latitude and longitude limits
     */
    for (idet = 0; idet < nlin; idet++) {
        for (ipix = 0; ipix < in_rec->npix; ipix++) {
            loc = idet * in_rec->npix + ipix;
            if (*(in_rec->lat + loc) < in_rec->ll_lims[0])
                in_rec->ll_lims[0] = *(in_rec->lat + loc);
            if (*(in_rec->lat + loc) > in_rec->ll_lims[1])
                in_rec->ll_lims[1] = *(in_rec->lat + loc);
            if ((lon = *(in_rec->lon + loc)) >= 0) {
                if (lon < in_rec->ll_lims[2]) in_rec->ll_lims[2] = lon;
                if (lon > in_rec->ll_lims[3]) in_rec->ll_lims[3] = lon;
            } else {
                if (lon < in_rec->ll_lims[4]) in_rec->ll_lims[4] = lon;
                if (lon > in_rec->ll_lims[5]) in_rec->ll_lims[5] = lon;
            }
        }
    }
    /*  sensor azimuth */
    if (h5io_rd_ds_slice(&(in_rec->geo_dat_id[2]), start, count,
            (void *) (in_rec->sena)) != 0) {
        printf("%s, %d: sena input failed on scan: %d\n", __FILE__, __LINE__
                , iscn);
        return 1;
    }
    /*  sensor zenith */
    if (h5io_rd_ds_slice(&(in_rec->geo_dat_id[3]), start, count,
            (void *) (in_rec->senz)) != 0) {
        printf("%s, %d: senz input failed on scan: %d\n", __FILE__, __LINE__
                , iscn);
        return 1;
    }
    /*  solar azimuth */
    if (h5io_rd_ds_slice(&(in_rec->geo_dat_id[4]), start, count,
            (void *) (in_rec->sola)) != 0) {
        printf("%s, %d: sola input failed on scan: %d\n", __FILE__, __LINE__,
                iscn);
        return 1;
    }
    /*  solar zenith */
    if (h5io_rd_ds_slice(&(in_rec->geo_dat_id[5]), start, count,
            (void *) (in_rec->solz)) != 0) {
        printf("%s, %d: solz input failed on scan: %d\n", __FILE__, __LINE__,
                iscn);
        return 1;
    }
    return 0;
}
