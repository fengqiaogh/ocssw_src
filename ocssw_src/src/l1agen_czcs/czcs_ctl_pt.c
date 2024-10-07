#include <stdlib.h>
#include "l1czcs.h"
#define N_ANCHOR 77
#define PI 3.14159265359
#define RADEG 180. / PI

void czcs_ctl_pt(DATA_REC_TYPE data, gattr_struc *gattr, int line,
        l1_data_struc *l1_data)
/*******************************************************************

   czcs_ctl_pt

   purpose: take the unevenly spaced anchor points from the CZCS CRTT file,
            correct the navigation and make an evenly spaced set of
            control points

   Returns type: void

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      DATA_REC_TYPE     data             I      CRTT data record
      gattr_struc *     gattr            I      structure with final 
                                                attributes
      int               line             I      line being processed
      l1_data_struc *   l1_data         I/O     data structure containing 
                                                the control point grid

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       19-Mar-2004     Original development

 *******************************************************************/ {
    int j, cpix, start_year, start_day;
    float tmplon, alon[N_ANCHOR], alat[N_ANCHOR], xlon[NCZCS_PIX],
            ylat[NCZCS_PIX], senz[NCZCS_PIX], sena[NCZCS_PIX];
    float suna, sunz, gmt;
    double pi, radeg;
    /*
     * First, convert the scaled values to floating point
     * and convert lon from 0 - 360 system to -180 - +180 system
     */
    for (j = 0; j < N_ANCHOR; j++) {
        alat[j] = fixed_pt_2_floating_pt(data.lat_anchr_pts[j], 22);
        tmplon = fixed_pt_2_floating_pt(data.long_anchr_pts[j], 22);
        if (tmplon > 180.) tmplon = tmplon - 360.;
        alon[j] = tmplon;
    }
    /*
     * expand anchor points to every pixel and initially correct the 
     * latitudes, longitudes
     */
    lonlat_(alon, alat, xlon, ylat);
    /*
     * get the sensor zenith angle (needed for final nav correction)
     */
    pi = PI;
    radeg = RADEG;
    satang_(&pi, &radeg, &(gattr->tilt), &(gattr->roll), &(gattr->pitch),
            &(gattr->yaw), xlon, ylat, senz, sena);
    /*
     * do the final navigation correction
     */
    lladjust_(&(gattr->tilt), &(gattr->roll), &(gattr->pitch),
            &(gattr->yaw), senz, xlon, ylat);
    /*
     *  as test for Fred, code to call satang using updated lat, lon
     */
    /*
     satang_( &pi, &radeg, &(gattr->tilt), &(gattr->roll), &(gattr->pitch),
       &(gattr->yaw), xlon, ylat, senz, sena );
     */
    /*
     *  get the solar zenith for the center and for testing, get 
     *  the solar data for the whole line
     */
    if (line == gattr->scan_lines / 2) {
        cpix = gattr->pix_per_scan / 2;
        gmt = (float) l1_data->msec[line] / (1000. * 3600.);
        start_year = gattr->start_year; /* start year, day are shorts */
        start_day = gattr->start_day; /* and sunangs need ints      */
        sunangs_(&start_year, &start_day, &gmt,
                &xlon[ cpix ], &ylat[ cpix ], &sunz, &suna);
        gattr->cntr_sol_zen = sunz;
    }
    /*
     *  for extra outputs of geometry and calibrated radiance data
     */
#ifdef GEOM_CAL
    {
        int off;
        /*
         *  derive the solar geometry for all points and fill struct with
         *  it and sensor info
         */
        gmt = (float) l1_data->msec[line] / (1000. * 3600.);
        start_year = gattr->start_year;
        start_day = gattr->start_day;
        off = (line * gattr->pix_per_scan);
        for (j = 0; j < gattr->pix_per_scan; j++) {
            sunangs_(&start_year, &start_day, &gmt,
                    &xlon[ j ], &ylat[ j ], &sunz, &suna);
            *(l1_data->sol_zen + off) = sunz;
            *(l1_data->sol_az + off) = suna;
            *(l1_data->sen_zen + off) = *(senz + j);
            *(l1_data->sen_az + off) = *(sena + j);
            *(l1_data->all_lat + off) = *(ylat + j);
            *(l1_data->all_lon + off) = *(xlon + j);
            off++;
        }
    }
#endif
    /*
     *  set up the lat, lon control points
     */
    for (j = 0; j < gattr->n_ctl_pt; j++) {
        *(l1_data->ctl_pt_lat + (line * gattr->n_ctl_pt) + j) =
                ylat[ l1_data->ctl_pt_cols[j] - 1 ];
        *(l1_data->ctl_pt_lon + (line * gattr->n_ctl_pt) + j) =
                xlon[ l1_data->ctl_pt_cols[j] - 1 ];
    }
    /*
     *  and finish
     */
    return;
}
