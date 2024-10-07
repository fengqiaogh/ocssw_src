#include "l1czcs.h"
#include <navigation.h>

void cz_meta_adj(l1_data_struc *l1_data, gattr_struc *gattr)
/*******************************************************************

   cz_meta_adj

   purpose: adjust the metadata for the newly merged CZCS file
    (note that it starts with the metadata from the first L1)

   Returns type: void

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      l1_data_struc *   l1_data          I      level-1 data struct
      gattr_struc *     gattr           I/O     global attribute struct

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 2 Sep 2004      Original development

 *******************************************************************/
 {
    int32_t cmsec;
    int clin, year, day, hr, min, msec;
    float gmt, suna, sunz;
    /*
     * Note that some items are updated in cz_mov_scn
     */
    /*
     *  center time need doing and
     *  # scan control points
     */
    gattr->n_ctl_lin = gattr->scan_lines;
    /*
     * center time is gotten from msec of center line
     */
    clin = (gattr->scan_lines + 1) / 2;
    cmsec = l1_data->msec[clin];
    if (cmsec < l1_data->msec[0]) {
        year = gattr->end_year;
        day = gattr->end_day;
    } else {
        year = gattr->start_year;
        day = gattr->start_day;
    }
    hr = cmsec / 3600000;
    min = (cmsec / 60000) % 60;
    msec = cmsec % 60000;
    sprintf(gattr->center_time, "%4d%3.3d%2.2d%2.2d%5.5d", year, day,
            hr, min, msec);

    /*
     *  get the center solar zenith from the center lat, lon, and time
     */
    cdata_();
    gmt = (float) gattr->start_msec / (1000. * 3600.);
    year = (int) gattr->start_year;
    day = (int) gattr->start_day;
    sunangs_(&year, &day, &gmt, &gattr->center_lon, &gattr->center_lat,
            &sunz, &suna);
    gattr->cntr_sol_zen = sunz;

    return;
}
