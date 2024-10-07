#include "l1czcs.h"

void cz_sd_set(l1_data_struc *l1_data, gattr_struc *gattr)
/*******************************************************************

   cz_sd_set

   purpose: Propagate lat, lon extremes and some global
            attribute values into the line-by-line values

   Returns type: void

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      l1_data_struc *   l1_data          I      arrays data counts, lat, lons
      gattr_struc *     gattr           I/O     structure with final 
                                                attributes

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       7-Sep-2004      Original development

 *******************************************************************/ {
    int i, j, nctlpix, nlin, ctr_ctl_pt;
    /*
     *  set up some sizes
     */
    nctlpix = gattr->n_ctl_pt;
    nlin = gattr->scan_lines;
    ctr_ctl_pt = (nctlpix - 1) / 2;
    /*
     *  set up the lat, lon, tilt, attitude, slope, intercept and gain per line
     */
    for (i = 0; i < nlin; i++) {
        l1_data->slat[i] = l1_data->ctl_pt_lat[ i * nctlpix ];
        l1_data->slon[i] = l1_data->ctl_pt_lon[ i * nctlpix ];
        l1_data->clat[i] = l1_data->ctl_pt_lat[ i * nctlpix + ctr_ctl_pt ];
        l1_data->clon[i] = l1_data->ctl_pt_lon[ i * nctlpix + ctr_ctl_pt ];
        l1_data->elat[i] = l1_data->ctl_pt_lat[ i * nctlpix + nctlpix - 1 ];
        l1_data->elon[i] = l1_data->ctl_pt_lon[ i * nctlpix + nctlpix - 1 ];

        l1_data->tilt[i] = gattr->tilt;
        l1_data->att_ang[ i * 3 ] = gattr->yaw;
        l1_data->att_ang[ i * 3 + 1 ] = gattr->roll;
        l1_data->att_ang[ i * 3 + 2 ] = gattr->pitch;
        for (j = 0; j < 6; j++) {
            l1_data->slope[ i * 6 + j ] = gattr->slope[j];
            l1_data->intercept[ i * 6 + j ] = gattr->intercept[j];
        }
        l1_data->gain[i] = (short) gattr->gain;
    }
    return;
}
