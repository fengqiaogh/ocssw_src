#include <stdio.h>
#include "l1czcs.h"

void hdr_2_gattr(HEADER2_TYPE header2, gattr_struc *gattr)
/*******************************************************************

   hdr_2_gattr

   purpose: extract and convert the CRTT header data into values needed 
            for the output file global attributes

   Returns type: void

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      HEADER2TYPE       header2          I      structure with CRTT header
      gattr_struc *     gattr           I/O     structure with final 
                                                attributes

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson        9-Mar-2004     Original development

 *******************************************************************/ {
    int i;
    unsigned short us;
    /*
     *  Just fill the needed values, converting when necessary
     */
    gattr->start_year = header2.starting_yr;
    gattr->start_day = header2.starting_day;
    gattr->start_msec = header2.start_msec;

    time_str(header2.starting_yr, header2.starting_day, header2.start_msec,
            gattr->start_time);
    /*
     *  get the scene center time
     */
    time_str(header2.scene_cntr_yr, header2.scene_cntr_doy,
            header2.scene_cntr_msec, gattr->center_time);

    gattr->orbit = (unsigned short) header2.orbit;
    gattr->gain = header2.gain;
    gattr->thresh = header2.threshold;
    gattr->tilt = header2.tilt / 1000.;

    gattr->center_lat = header2.lat_cntr / 100. - 90.;
    us = (unsigned short) header2.long_cntr;
    gattr->center_lon = (us > 18000) ? us / 100. - 360. : us / 100.;

    gattr->pix_per_scan = 1968;
    gattr->scan_lines = header2.no_of_scans;

    /* North limit, choose among top corners
  printf( "long corners: fitr %d, litr %d, fitl %d, litl %d\n", header2.long_crnr_fitr, header2.long_crnr_litr, header2.long_crnr_fitl, header2.long_crnr_litl );
  printf( "lat corners: fitr %d, litr %d, fitl %d, litl %d\n", header2.lat_crnr_fitr, header2.lat_crnr_litr, header2.lat_crnr_fitl, header2.lat_crnr_litl );
     */

    gattr->limits[0] = (header2.lat_crnr_litl > header2.lat_crnr_litr) ?
            header2.lat_crnr_litl / 100. : header2.lat_crnr_litr / 100.;
    /* South */
    gattr->limits[1] = (header2.lat_crnr_fitl < header2.lat_crnr_fitr) ?
            header2.lat_crnr_fitl / 100. : header2.lat_crnr_fitr / 100.;
    /* West */
    gattr->limits[2] = (header2.long_crnr_fitl < header2.long_crnr_litl) ?
            header2.long_crnr_fitl / 100. : header2.long_crnr_litl / 100.;
    /* East */
    gattr->limits[3] = (header2.long_crnr_fitr > header2.long_crnr_litr) ?
            header2.long_crnr_fitr / 100. : header2.long_crnr_litr / 100.;

    /*
     * slope and intercept
     */
    for (i = 0; i < 6; i++) {
        gattr->slope[i] = fixed_pt_2_floating_pt(header2.chan_line[i].slope, 24);
        gattr->intercept[i] =
                fixed_pt_2_floating_pt(header2.chan_line[i].intercept, 24);
    }
    /*
     *  pitch, roll and yaw
     */
    gattr->roll = header2.cntr_roll / 1000.;
    gattr->pitch = header2.cntr_pitch / 1000.;
    gattr->yaw = header2.cntr_yaw / 1000.;
    /*
     *  file metric info (is any of this really used, we'll find out)
     */
    gattr->ilt_flags = header2.ilt_flags;
    gattr->parm_presence = header2.parameter_presence;
    gattr->n_miss_scans = header2.no_of_missing_scans_all;
    for (i = 0; i < 6; i++)
        gattr->n_scan_mis_chan[i] = header2.no_of_missing_scans[i];
    gattr->n_hdt_sync_loss = header2.no_of_hdt_sync_losses;
    gattr->n_hdt_parity_err = header2.no_of_hdt_parity_errs;
    gattr->n_wbvt_sync_loss = header2.no_of_wbvt_sync_losses;
    gattr->n_wbvt_slips = header2.no_of_wbvt_bit_slips;
    return;
}
