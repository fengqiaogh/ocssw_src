#include <limits.h>
#include "l1czcs.h"

int cz_mov_scn(int init, int ds_num, char *file, mstr_struc *mstr,
        int mstr_last, int nscan_out, gattr_struc *gattr, l1_data_struc *l1_data)
/*******************************************************************

   cz_mov_scn

   purpose: move the proper data scan information from the input L1
      czcs file to the output storage

   Returns type: int - return status: 0 is good

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               init             I      set up global attributes
                                                flag: 1, set gattr from 1st
                                                file.  0, only modify with
                                                incoming attrs
      int               ds_num           I      dataset id #
      char *            file             I      L1 czcs file to transfer 
                                                scans from
      mstr_struc *      mstr             I      master control structure
      int               mstr_last        I      last entry in the master list
      int               nscan_out        I      # lines in merged file
      gattr_struc *     gattr           I/O     global attributes for output
      l1_data_struc *   l1_data         I/O     data structure for output

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 10 Aug 2004     Original development
      W. Robinson, SAIC 20 Dec 2005     add the moving of the pos_err or
                                        position error values

 *******************************************************************/
 {
    l1_data_struc in_dat;
    gattr_struc in_attr;
    int i, j, n_ctl, nscan, time_reset;
    short in_scn;
    int32_t out_scn, sum;
    /*
     *  read in the L1 data and attributes from the file
     *  fill output attr struct if first file
     */
    if (init == 1) {
        if (cz_l1_read(file, 0, gattr, &in_dat) != 0) return -1;
        n_ctl = gattr->n_ctl_pt;
        nscan = gattr->scan_lines;
        /*
         *  copy ctl_pt_cols entirely the first time
         */
        memcpy(l1_data->ctl_pt_cols, in_dat.ctl_pt_cols,
                (n_ctl * sizeof ( int)));
    } else {
        if (cz_l1_read(file, 0, &in_attr, &in_dat) != 0) return -1;
        n_ctl = in_attr.n_ctl_pt;
        nscan = in_attr.scan_lines;
        /*
         *  Update some of the global attributes, Note that some items are 
         *  updated in cz_meta_adj and also in czl1merge as well as in cz_l1_write
         *
         *  end time as the latest end time seen
         */
        time_reset = 0;
        if (in_attr.end_year > gattr->end_year)
            time_reset = 1;
        else if (in_attr.end_day > gattr->end_day)
            time_reset = 1;
        else if (in_attr.end_msec > gattr->end_msec)
            time_reset = 1;

        if (time_reset == 1) {
            gattr->end_year = in_attr.end_year;
            gattr->end_day = in_attr.end_day;
            gattr->end_msec = in_attr.end_msec;
            strcpy(gattr->end_time, in_attr.end_time);
        }

        /*
         *  file metrics will best represent the cumulative scan info
         *
         *  the ILT flags and parameter presence will be the bitwise or
         *  of the current and new values
         */
        gattr->ilt_flags = gattr->ilt_flags | in_attr.ilt_flags;
        gattr->parm_presence = gattr->parm_presence | in_attr.parm_presence;
        /*
         *  for the other file metrics, sum together till short limit is reached
         */
        sum = (int32_t) gattr->n_miss_scans + (int32_t) in_attr.n_miss_scans;
        gattr->n_miss_scans = (sum > SHRT_MAX) ? SHRT_MAX : (short) sum;

        for (i = 0; i < 6; i++) {
            sum = (int32_t) gattr->n_scan_mis_chan[i] + (int32_t) in_attr.n_scan_mis_chan[i];
            gattr->n_scan_mis_chan[i] = (sum > SHRT_MAX) ? SHRT_MAX : (short) sum;
        }

        sum = (int32_t) gattr->n_hdt_sync_loss + (int32_t) in_attr.n_hdt_sync_loss;
        gattr->n_hdt_sync_loss = (sum > SHRT_MAX) ? SHRT_MAX : (short) sum;

        sum = (int32_t) gattr->n_hdt_parity_err + (int32_t) in_attr.n_hdt_parity_err;
        gattr->n_hdt_parity_err = (sum > SHRT_MAX) ? SHRT_MAX : (short) sum;

        sum = (int32_t) gattr->n_wbvt_sync_loss + (int32_t) in_attr.n_wbvt_sync_loss;
        gattr->n_wbvt_sync_loss = (sum > SHRT_MAX) ? SHRT_MAX : (short) sum;

        sum = (int32_t) gattr->n_wbvt_slips + (int32_t) in_attr.n_wbvt_slips;
        gattr->n_wbvt_slips = (sum > SHRT_MAX) ? SHRT_MAX : (short) sum;

    }
    printf("cz_mov_scn:  dataset # (ds_num) = %d,  # scans: %d\n",
            ds_num, nscan);
    printf("cz_mov_scn:  file: %s\n", file);
    /*
     *  transfer each portion of the scan info as indicated by the master list
     *  note that mstr_last is pointer to last, hence <=
     */
    for (i = 0; i <= mstr_last; i++) {
        if ((mstr->ds_num[i] == ds_num) && (mstr->exist[i] == 1)) {
            in_scn = mstr->scan[i];
            out_scn = mstr->out_scan[i];
            /*  for detail of output
             printf( 
         "cz_mov_scn:  moving scan line %d  to output line %ld, mstr loop indx: %d\n", 
               in_scn, out_scn, i );
             */
            for (j = 0; j < 6; j++)
                memcpy(l1_data->counts[j] + (NCZCS_PIX * out_scn),
                    in_dat.counts[j] + (NCZCS_PIX * in_scn),
                    (NCZCS_PIX * sizeof (unsigned char)));

            l1_data->msec[out_scn] = in_dat.msec[in_scn];
            l1_data->tilt[out_scn] = in_dat.tilt[in_scn];
            l1_data->slat[out_scn] = in_dat.slat[in_scn];
            l1_data->slon[out_scn] = in_dat.slon[in_scn];
            l1_data->clat[out_scn] = in_dat.clat[in_scn];
            l1_data->clon[out_scn] = in_dat.clon[in_scn];
            l1_data->elat[out_scn] = in_dat.elat[in_scn];
            l1_data->elon[out_scn] = in_dat.elon[in_scn];
            memcpy(l1_data->cal_sum + (5 * out_scn),
                    in_dat.cal_sum + (5 * in_scn), (5 * sizeof (unsigned char)));
            memcpy(l1_data->cal_scan + (6 * out_scn),
                    in_dat.cal_scan + (6 * in_scn), (6 * sizeof (unsigned char)));
            memcpy(l1_data->orb_vec + (3 * out_scn),
                    in_dat.orb_vec + (3 * in_scn), (3 * sizeof ( float)));
            memcpy(l1_data->att_ang + (3 * out_scn),
                    in_dat.att_ang + (3 * in_scn), (3 * sizeof ( float)));
            l1_data->pos_err[out_scn] = in_dat.pos_err[in_scn];
            memcpy(l1_data->slope + (6 * out_scn), in_dat.slope + (6 * in_scn),
                    (6 * sizeof ( float)));
            memcpy(l1_data->intercept + (6 * out_scn),
                    in_dat.intercept + (6 * in_scn), (6 * sizeof ( float)));
            l1_data->gain[out_scn] = in_dat.gain[in_scn];
            memcpy(l1_data->ctl_pt_lat + (out_scn * n_ctl),
                    in_dat.ctl_pt_lat + (in_scn * n_ctl), (n_ctl * sizeof (float)));
            memcpy(l1_data->ctl_pt_lon + (out_scn * n_ctl),
                    in_dat.ctl_pt_lon + (in_scn * n_ctl), (n_ctl * sizeof (float)));
            l1_data->ctl_pt_rows[out_scn] = out_scn + 1;
        }
    }
    /*
     *  remove allocated arrays in data struct
     */
    cz_dat_free(&in_dat, 0);
    return 0;
}
