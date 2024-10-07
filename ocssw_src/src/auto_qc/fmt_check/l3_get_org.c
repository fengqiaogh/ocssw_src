#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fmt_check.h"
#define PI      3.141592653589793
#define RADCONV PI/180.L

void l3_get_org(int resolve, l3_org_str *l3_org)
/*******************************************************************

   fmt_rd_dim

   purpose: set up the organization information for a L3 based on the 
      resolve value

   Returns type: none

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               resolve          I      approximate bin size in km
                                                currently valid values:
                                                1, 2, 4, 9, 36
      l3_org_str *      l3_org          I/O     organization information

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 5 Apr 2005      Original development

 *******************************************************************/
 {
    int res_val[5] = {1, 2, 4, 9, 36};
    int numrows[5] = {17280, 8640, 4320, 2160, 540};
    float64 latbin;
    int i;

    /*
     *  first, find the # rows for the resolve value
     */
    for (i = 0; i < 5; i++) {
        if (resolve == res_val[i]) {
            l3_org->numrows = numrows[i];
            break;
        }
    }

    /*
     *  get the bins at the equator and the vertical size
     */
    l3_org->bins_eq = (int32) l3_org->numrows * 2;
    l3_org->vsize = (float64) (180.L / l3_org->numrows);

    /*
     *  For each row, set up the rest
     */
    l3_org->hsize = (float64 *) malloc(l3_org->numrows * sizeof ( float64));
    l3_org->start_bin = (int32 *) malloc(l3_org->numrows * sizeof ( int32));
    l3_org->max_bin = (int32 *) malloc(l3_org->numrows * sizeof ( int32));

    l3_org->start_bin[0] = 1;
    latbin = (90.0L / l3_org->numrows) - 90.0L;
    l3_org->max_bin[0] = (int32) (cos(latbin * RADCONV) *
            l3_org->bins_eq + 0.5);
    l3_org->hsize[0] = 360. / l3_org->max_bin[0];

    /*
     printf( "WDR debug %s: line 0, lat, start, # = %f  %d  %d\n", __FILE__,
     latbin, l3_org->start_bin[0], l3_org->max_bin[0] );
     */

    for (i = 1; i < l3_org->numrows; i++) {
        latbin = ((float64) i + 0.5) * (180.0L / l3_org->numrows) - 90.0L;
        l3_org->max_bin[i] = (int32) (cos(latbin * RADCONV) *
                l3_org->bins_eq + 0.5);
        l3_org->start_bin[i] = l3_org->start_bin[i - 1] +
                l3_org->max_bin[i - 1];
        l3_org->hsize[i] = 360. / l3_org->max_bin[i];
        /*
         printf( "WDR debug %s: line %d, lat, start, # = %f  %d  %d\n", __FILE__,
         i, latbin, l3_org->start_bin[i], l3_org->max_bin[i] );
         */
    }
    return;
}
