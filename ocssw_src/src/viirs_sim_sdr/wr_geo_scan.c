#include "viirs_sim_sdr.h"

int wr_geo_scan(int iscn, out_rec_struc *out_rec)
/*-----------------------------------------------------------------------------
    Program:   wr_sdr_scan.c

    Description:  write scan of data to the VIIRS geolocation SDR

    Arguments:
        Type      Name         I/O   Description
        ----      ----         ---   -----------
        int       iscn          I    scan to write
        out_rec_struc * out_rec   I/O  output record controls

    Modification history:

    W. Robinson, SAIC  21 Oct 2008  Original development
    W. Robinson, SAIC  15 Mar 2010  reconfigure to scan oriented output

----------------------------------------------------------------------------*/ {
    int ilin, start[2], count[2];
    /*
     *  latitude
     */
    ilin = iscn * out_rec->ndet_scan;
    start[0] = ilin;
    start[1] = 0;
    count[0] = out_rec->ndet_scan;
    count[1] = out_rec->npix;
    if (h5io_wr_ds_slice(&(out_rec->geo_dat_id[0]), start, count,
            out_rec->lat) != 0) {
        printf("%s: lat output failed on scan: %d\n", __FILE__, iscn);
        return 1;
    }
    /* longitude */
    if (h5io_wr_ds_slice(&(out_rec->geo_dat_id[1]), start, count,
            out_rec->lon) != 0) {
        printf("%s: lon output failed on scan: %d\n", __FILE__, iscn);
        return 1;
    }
    /* sensor azimuth */
    if (h5io_wr_ds_slice(&(out_rec->geo_dat_id[2]), start, count,
            out_rec->sena) != 0) {
        printf("%s: sena output failed on scan: %d\n", __FILE__, iscn);
        return 1;
    }
    /* sensor zenith */
    if (h5io_wr_ds_slice(&(out_rec->geo_dat_id[3]), start, count,
            out_rec->senz) != 0) {
        printf("%s: senz output failed on scan: %d\n", __FILE__, iscn);
        return 1;
    }
    /* solar azimuth */
    if (h5io_wr_ds_slice(&(out_rec->geo_dat_id[4]), start, count,
            out_rec->sola) != 0) {
        printf("%s: sola output failed on scan: %d\n", __FILE__, iscn);
        return 1;
    }
    /* solar zenith */
    if (h5io_wr_ds_slice(&(out_rec->geo_dat_id[5]), start, count,
            out_rec->solz) != 0) {
        printf("%s: solz output failed on scan: %d\n", __FILE__, iscn);
        return 1;
    }
    return 0;
}
