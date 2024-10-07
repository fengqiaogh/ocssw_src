#include "l1czcs.h"

int cz_clean(gattr_struc *gattr, l1_data_struc *l1_data)
/*******************************************************************

   cz_clean

   purpose: clean up the czcs data for the following
     - bad time in the lines (solution - remove lines)

   Returns type: int - 0 no problems, else no good lines left in file

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      gattr_struc *     gattr           I/O     structure with final 
                                                attributes
      l1_data_struc *   l1_data         I/O     arrays data counts, lat, lons

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       23-Sep-2004     Original development, 
      W. Robinson, SAIC 19 Dec 2005     add the pos_err array to this

 *******************************************************************/ {
    int i, linout = 0, nlin, npix, nctl_pix, msec, icompr, j, ncompr;

    /*
     *  loop through the scan lines, any bad lines are compressed out
     *  and the # scan lines is adjusted accordingly
     */
    nlin = gattr->scan_lines;
    npix = gattr->pix_per_scan;
    nctl_pix = gattr->n_ctl_pt;
    icompr = 0;

    for (i = 0; i < nlin; i++) {
        /*
         *  currently, the only check is on the msec in each line
         */
        msec = *(l1_data->msec + i);
        /*
         *  if a bad msec found, don't incriment the output line and note that 
         *  compression must be done
         */
        if ((msec < 0) || (msec > 86399999)) {
            icompr = 1;
        } else {
            /*
             *  once the compress switch is on, copy lines to proper place
             */
            if (icompr == 1) {
                for (j = 0; j < 6; j++)
                    memcpy(l1_data->counts[j] + npix * linout,
                        l1_data->counts[j] + npix * i, npix);
                *(l1_data->msec + linout) = *(l1_data->msec + i);
                memcpy(l1_data->ctl_pt_lat + nctl_pix * linout,
                        l1_data->ctl_pt_lat + nctl_pix * i,
                        nctl_pix * sizeof ( float));
                memcpy(l1_data->ctl_pt_lon + nctl_pix * linout,
                        l1_data->ctl_pt_lon + nctl_pix * i,
                        nctl_pix * sizeof ( float));
                /*  note that ctl_pt_cols, ctl_pt_rows don't need change  */
                *(l1_data->tilt + linout) = *(l1_data->tilt + i);
                *(l1_data->slat + linout) = *(l1_data->slat + i);
                *(l1_data->slon + linout) = *(l1_data->slon + i);
                *(l1_data->clat + linout) = *(l1_data->clat + i);
                *(l1_data->clon + linout) = *(l1_data->clon + i);
                *(l1_data->elat + linout) = *(l1_data->elat + i);
                *(l1_data->elon + linout) = *(l1_data->elon + i);
                memcpy(l1_data->cal_sum + 5 * linout,
                        l1_data->cal_sum + 5 * i, 5 * sizeof ( unsigned char));
                memcpy(l1_data->cal_scan + 6 * linout,
                        l1_data->cal_scan + 6 * i, 6 * sizeof ( unsigned char));
                memcpy(l1_data->orb_vec + 3 * linout,
                        l1_data->orb_vec + 3 * i, 3 * sizeof ( float));
                memcpy(l1_data->att_ang + 3 * linout,
                        l1_data->att_ang + 3 * i, 3 * sizeof ( float));
                *(l1_data->pos_err + linout) = *(l1_data->pos_err + i);
                memcpy(l1_data->slope + 6 * linout,
                        l1_data->slope + 6 * i, 6 * sizeof ( float));
                memcpy(l1_data->intercept + 6 * linout,
                        l1_data->intercept + 6 * i, 6 * sizeof ( float));
                *(l1_data->gain + linout) = *(l1_data->gain + i);
#ifdef GEOM_CAL
                memcpy(l1_data->sen_zen + npix * linout,
                        l1_data->sen_zen + npix * i, npix * sizeof ( float));
                memcpy(l1_data->sen_az + npix * linout,
                        l1_data->sen_az + npix * i, npix * sizeof ( float));
                memcpy(l1_data->sol_zen + npix * linout,
                        l1_data->sol_zen + npix * i, npix * sizeof ( float));
                memcpy(l1_data->sol_az + npix * linout,
                        l1_data->sol_az + npix * i, npix * sizeof ( float));
                memcpy(l1_data->all_lat + npix * linout,
                        l1_data->all_lat + npix * i, npix * sizeof ( float));
                memcpy(l1_data->all_lon + npix * linout,
                        l1_data->all_lon + npix * i, npix * sizeof ( float));
                memcpy(l1_data->Lt_443 + npix * linout,
                        l1_data->Lt_443 + npix * i, npix * sizeof ( float));
                memcpy(l1_data->Lt_520 + npix * linout,
                        l1_data->Lt_520 + npix * i, npix * sizeof ( float));
                memcpy(l1_data->Lt_550 + npix * linout,
                        l1_data->Lt_550 + npix * i, npix * sizeof ( float));
                memcpy(l1_data->Lt_670 + npix * linout,
                        l1_data->Lt_670 + npix * i, npix * sizeof ( float));
                memcpy(l1_data->Lt_750 + npix * linout,
                        l1_data->Lt_750 + npix * i, npix * sizeof ( float));
                memcpy(l1_data->Lt_11500 + npix * linout,
                        l1_data->Lt_11500 + npix * i, npix * sizeof ( float));
#endif
            }
            linout++;
        }
    }
    /*
     *  re-set the # lines and # control lines
     */
    gattr->scan_lines = linout;
    gattr->n_ctl_lin = linout;
    /*
     *  and end
     */
    if (icompr != 0) {
        ncompr = nlin - linout;
        printf("cz_clean: %d bad lines (of %d) were compressed out of the L1 file\n",
                ncompr, nlin);
    }
    if (linout == 0) {
        printf("cz_clean: No good lines were found in the file\n");
        return -1;
    } else
        return 0;
}
