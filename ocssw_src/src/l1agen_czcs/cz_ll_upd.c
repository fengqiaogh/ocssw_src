#include "l1czcs.h"

void cz_ll_upd(l1_data_struc *l1_data, gattr_struc *gattr)
/*******************************************************************

   cz_ll_upd

   purpose: Update the latitude, longitude information using the 
            improved navigation.  

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
      W. Robinson        5-Apr-2004     Original development

 *******************************************************************/ {
    int i, j, nctlpix, nlin, ctr_ctl_pt;
    int32_t loc;
    float maxlat, minlat, low_neg, low_pos, hi_neg, hi_pos;
    /*
     *  the western longitude is the longitude of the first pixel of the last line
     *  and the eastern is the last pixel of the first line
     */
    nctlpix = gattr->n_ctl_pt;
    nlin = gattr->scan_lines;
    /*
     *  We need to use the most general method to get the dataset bounds.
     *  get the min, max lat and for longitude, get min, max in (+) and
     *  (-) long ranges.  The W and E get determined from these
     */
    low_neg = 999.;
    low_pos = 999.;
    hi_neg = -999.;
    hi_pos = -999.;
    maxlat = -999.;
    minlat = 999.;
    loc = 0;
    for (i = 0; i < nlin; i++) {
        /*
         *  Avoid getting information from lines where the data is missing in
         *  the vis and NIR bands and where scan quality(4, 5) are bad
         */
        if ((l1_data->cal_sum[ i * 5 + 3 ] == 0) &&
                (l1_data->cal_sum[ i * 5 + 4 ] == 0) &&
                (l1_data->cal_scan[ i * 6 ] == 0) &&
                (l1_data->cal_scan[ i * 6 + 1 ] == 0) &&
                (l1_data->cal_scan[ i * 6 + 2 ] == 0) &&
                (l1_data->cal_scan[ i * 6 + 3 ] == 0) &&
                (l1_data->cal_scan[ i * 6 + 4 ] == 0)) {
            loc = i * nctlpix;
            for (j = 0; j < nctlpix; j++) {
                if (l1_data->ctl_pt_lat[ loc ] > maxlat)
                    maxlat = l1_data->ctl_pt_lat[ loc ];
                if (l1_data->ctl_pt_lat[ loc ] < minlat)
                    minlat = l1_data->ctl_pt_lat[ loc ];
                if (l1_data->ctl_pt_lon[ loc ] >= 0.) {
                    if (l1_data->ctl_pt_lon[ loc ] > hi_pos)
                        hi_pos = l1_data->ctl_pt_lon[ loc ];
                    if (l1_data->ctl_pt_lon[ loc ] < low_pos)
                        low_pos = l1_data->ctl_pt_lon[ loc ];
                } else {
                    if (l1_data->ctl_pt_lon[ loc ] > hi_neg)
                        hi_neg = l1_data->ctl_pt_lon[ loc ];
                    if (l1_data->ctl_pt_lon[ loc ] < low_neg)
                        low_neg = l1_data->ctl_pt_lon[ loc ];
                }
                loc++;
            }
        }
    }
    /*
     *  latitude is done, so set it.  for longitude, there are 5 conditions
     *  limits are 0 - N, 1 - S, 2 - W, 3 - E
     */
    gattr->limits[0] = maxlat;
    gattr->limits[1] = minlat;
    if ((low_neg == 999.) && (low_pos != 999.)) {
        gattr->limits[2] = low_pos; /* all + long */
        gattr->limits[3] = hi_pos;
    } else if ((low_neg != 999.) && (low_pos == 999.)) {
        gattr->limits[2] = low_neg; /* all - long */
        gattr->limits[3] = hi_neg;
    } else if ((hi_neg > -90.) && (low_pos < 90)) {
        gattr->limits[2] = low_neg; /* long straddles 0 degrees */
        gattr->limits[3] = hi_pos;
    } else if ((low_neg <= -90.) && (hi_pos >= 90)) {
        gattr->limits[2] = low_pos; /* long straddles date line */
        gattr->limits[3] = hi_neg;
    } else {
        gattr->limits[2] = -180.; /* indeterminent */
        gattr->limits[3] = 180.;
    }
    /*
     *  corner and center start, end line lat and lon info 
     *  note that using aqua convention, this is dataset array relative,
     *  ie, top left is first line, pixel in array
     */
    gattr->up_lft_lat = l1_data->ctl_pt_lat[ 0 ];
    gattr->up_lft_lon = l1_data->ctl_pt_lon[ 0 ];
    gattr->lo_lft_lat = l1_data->ctl_pt_lat[ (nlin - 1) * nctlpix ];
    gattr->lo_lft_lon = l1_data->ctl_pt_lon[ (nlin - 1) * nctlpix ];
    gattr->up_rgt_lat = l1_data->ctl_pt_lat[ nctlpix - 1 ];
    gattr->up_rgt_lon = l1_data->ctl_pt_lon[ nctlpix - 1 ];
    gattr->lo_rgt_lat =
            l1_data->ctl_pt_lat[ (nlin - 1) * nctlpix + (nctlpix - 1) ];
    gattr->lo_rgt_lon =
            l1_data->ctl_pt_lon[ (nlin - 1) * nctlpix + (nctlpix - 1) ];

    ctr_ctl_pt = (nctlpix - 1) / 2;
    gattr->start_cntr_lat = l1_data->ctl_pt_lat[ ctr_ctl_pt ];
    gattr->start_cntr_lon = l1_data->ctl_pt_lon[ ctr_ctl_pt ];
    gattr->end_cntr_lat =
            l1_data->ctl_pt_lat[ (nlin - 1) * nctlpix + ctr_ctl_pt ];
    gattr->end_cntr_lon =
            l1_data->ctl_pt_lon[ (nlin - 1) * nctlpix + ctr_ctl_pt ];
    gattr->center_lat =
            l1_data->ctl_pt_lat[ (nlin - 1) / 2 * nctlpix + ctr_ctl_pt ];
    gattr->center_lon =
            l1_data->ctl_pt_lon[ (nlin - 1) / 2 * nctlpix + ctr_ctl_pt ];

    return;
}
