#include <math.h>

int b128_wd_bit(float lat_off, float lon_off, int *box_wd, int *box_bit)
/*******************************************************************

   b128_box_num

   purpose: for the 128 / degree land mask file, find the 
      index of the degree box from the lat, lon

   Returns type: int - nothing now, later, maybe an error condition

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      float             lat_off          I      all positive latitude 
                                                generated in call to 
                                                b128_box_num 
      float             lon_off          I      all positive longitude
                                                as above
      int *             box_wd           O      word # in array of 
                                                128 sq box to find
                                                desired mask value
      int *             box_bit          O      bit # in word of
                                                desired mask value

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       8-Oct-1997      Original development

 *******************************************************************/ {
    int lat_sub, lon_sub;
    int lon_wd;
    double dumb;

    /*
     *  find distance in to the bit from edge of 128 by 128 array
     */
    lat_sub = (int) (modf(lat_off, &dumb) * 128);
    if (lat_off == 180.) lat_sub = 127;

    lon_sub = (int) (modf(lon_off, &dumb) * 128);
    if (lon_off == 360.) lon_sub = 127;
    /*
     *  get the word in longitude direction and bit and then get linear word
     */
    lon_wd = lon_sub / 16;
    *box_bit = lon_sub - lon_wd * 16;
    *box_wd = lon_wd + lat_sub * 8;

    return 0;
}
