#include <stdio.h>
#include <fcntl.h>
#include "proto.h"
#include "mfhdf.h"
#include "calib_cal_l1a.h"

void free_all(void *dark_rest, void *gain, void *tdi, void *scan_temp, void *side, void *msec);

int compute_l1b(int32 sd_id, int iline, float *l1b_data[], char *calibtable, char newfile) {
    static int16 syear, sday, eday;
    static int32 smsec, Nlines, Npixels;
    static int16 *dark_rest = NULL, *gain = NULL, *tdi = NULL, *scan_temp = NULL, *side = NULL;
    static int32 *msec = NULL;
    static uint8 *sflags = NULL;
    static float32 dark_mean[NBANDS];
    static int32 sds_id, start[6], edge[6];
    static char dtype[4];
    int32 hdf_ind, rank, dimsize[6], numtype, nattrs;
    int16 *l1a_data_bip;
    float32 *l1b_data_bip;
    char name[128];
    int kline, ipix, ib;
    /*
    cal_mod_struc cal_mod;
    cal_mod.flag = 0;
     */

    /* Get attributes */
    if (newfile) {
        if (read_attr(sd_id, "Number of Scan Lines", &Nlines) == -1 ||
                read_attr(sd_id, "Pixels per Scan Line", &Npixels) == -1 ||
                read_attr(sd_id, "Data Type", &dtype) == -1)
            return -1;

        /* Read Level 1A data */
        if ((hdf_ind = SDnametoindex(sd_id, "l1a_data")) < 0) {
            fprintf(stderr, "can't locate SDS l1a_data\n");
            return -1;
        }
        if ((sds_id = SDselect(sd_id, hdf_ind)) < 0) {
            fprintf(stderr, "can't select SDS l1a_data\n");
            return -1;
        }
        if (SDgetinfo(sds_id, name, &rank, dimsize, &numtype, &nattrs) < 0) {
            fprintf(stderr, "can't get info for SDS l1a_data\n");
            return -1;
        }

        dark_rest = (int16 *) realloc(dark_rest, Nlines * NBANDS * sizeof (int16));
        gain = (int16 *) realloc(gain, Nlines * NBANDS * sizeof (int16));
        tdi = (int16 *) realloc(tdi, Nlines * NBANDS * sizeof (int16));
        scan_temp = (int16 *) realloc(scan_temp, Nlines * NBANDS * sizeof (int16));
        side = (int16 *) realloc(side, Nlines * sizeof (int16));
        msec = (int32 *) realloc(msec, Nlines * sizeof (int32));
        sflags = (uint8 *) realloc(sflags, Nlines * 4 * sizeof (uint8));

        if (get_calib_sds(sd_id, dark_rest, gain, tdi, scan_temp, side, msec) == -1) {
            fprintf(stderr, "can't compute L1B data\n");
            return -1;
        }
        if (read_sds(sd_id, "s_flags", sflags) == -1) return -1;

        /* Compute dark_mean */
        for (ib = 0; ib < NBANDS; ib++) {
            dark_mean[ib] = 0.;
            for (kline = 0; kline < Nlines; kline++) dark_mean[ib] += (float) dark_rest[kline * NBANDS + ib];
            dark_mean[ib] = dark_mean[ib] / Nlines;
        }

        edge[0] = 1;
        edge[1] = dimsize[1];
        edge[2] = dimsize[2];
        start[1] = 0;
        start[2] = 0;

    }

    if (sflags[4 * iline] > 127) {
        printf("Scan %d is corrupted. Skipped.\n", iline);
        return +1;
    }

    l1a_data_bip = (int16 *) malloc(Npixels * NBANDS * sizeof (int16));
    l1b_data_bip = (float32 *) malloc(Npixels * NBANDS * sizeof (float32));
    start[0] = iline;
    if (SDreaddata(sds_id, start, NULL, edge, l1a_data_bip) < 0) {
        fprintf(stderr, "can't read Level 1A data\n");
        return -1;
    }

    calibrate_l1a(calibtable, syear, sday, smsec, eday, msec[iline], dtype, 1, Npixels, &dark_rest[iline * NBANDS], dark_mean, &gain[iline * NBANDS], &tdi[iline * NBANDS], &scan_temp[iline * NBANDS], side[iline], l1a_data_bip, l1b_data_bip, NULL);

    for (ib = 0; ib < NBANDS; ib++)
        if (l1b_data[ib])
            for (ipix = 0; ipix < Npixels; ipix++)
                l1b_data[ib][ipix] = l1b_data_bip[ipix * NBANDS + ib];

    free(l1a_data_bip);
    free(l1b_data_bip);

    return 0;
}
