#include <math.h>

int b128_box_num(float lat, float lon, float *lat_off, float *lon_off)
/*******************************************************************

   b128_box_num

   purpose: for the 128 / degree land mask file, find the 
      index of the degree box from the lat, lon

   Returns type: int - index (0 rel) of the 1 degree box

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      float             lat              I      latitude of point from
                                                -90. to 90.
      float             lon              I      longitudeof point from
                                                -180 to 180.
      float *           lat_off          O      all positive latitude used
                                                during call to 
                                                b128_wd_bit to get mask
                                                at a 128th of a degree point
      float *           lon_off          O      all positive longitude

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       8-Oct-1997      Original development

 *******************************************************************/ {
    int box_num, lat_index, lon_index;
    /*
     *  make positive lat and lon and check them
     */
    *lat_off = 90. + lat;
    *lon_off = 180. + lon;

    if (*lat_off < 0.) *lat_off = 0.;
    if (*lat_off > 180.) *lat_off = 180.;
    if (*lon_off < 0.) *lon_off = 0.;
    if (*lon_off > 360.) *lon_off = 360.;
    /*
     *  Take care of case of 90 lat and 180 lon properly
     */
    lat_index = (*lat_off == 180.) ? 179 : (int) (*lat_off);
    lon_index = (*lon_off == 360.) ? 359 : (int) (*lon_off);
    /*
     *  compute the box #
     */
    box_num = lon_index + 360 * lat_index;

    return box_num;
}
