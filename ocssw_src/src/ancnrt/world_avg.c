#include <mfhdf.h>
#include "ancil.h"

void world_avg(int16 *inarr, int32 nlin, int32 npix, int avg_pix,
        int avg_lin, int16 *outarr)
/*******************************************************************

   world_avg

   purpose: average all non 0 array values inside a box, cross the
       longitude seam to do this

   Returns type: void no value

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int16 *           inarr           I       input array
      int32             nlin            I       # lines in the input array
      int32             npix            I       # pixels in the input array
      int               avg_pix         I       pixel size of box to average 
                                                with, will be forced to 1 
                                                if <1 and to npix - 1 
                                                if > npix - 1
      int               avg_lin         I       line size of box to average 
                                                with, will be forced to 1 
                                                if <1 and to nlin - 1 
                                                if > nlin - 1
      int16 *           outarr          O       output array

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       3-Jul-1997      Original development

 *******************************************************************/
 {
    int pix_minus, pix_plus, lin_minus, lin_plus;
    int num, avg, ipx, iln, jpx, jln, jpxm;
    /*
     *  check / adjust the box size
     */
    if (npix < 2 || nlin < 2) {
        printf("world_avg: either #pix (%d) or # lin (%d) < 2\n", npix, nlin);
        printf("           no averaging done (probably an error)\n");
        return;
    }
    if (avg_pix < 1) {
        printf("world_avg: avg_pix (%d) < 1, resetting to 1\n", avg_pix);
        avg_pix = 1;
    }
    if (avg_lin < 1) {
        printf("world_avg: avg_lin (%d) < 1, resetting to 1\n", avg_lin);
        avg_lin = 1;
    }
    if (avg_pix > npix - 1) {
        printf(
                "world_avg: avg_pix (%d) > # pixels - 1 (# pixels = %d), resetting to %d\n",
                avg_pix, npix, (npix - 1));
        avg_pix = npix - 1;
    }
    if (avg_lin > nlin - 1) {
        printf(
                "world_avg: avg_lin (%d) > # lines - 1 (# lines = %d), resetting to %d\n",
                avg_lin, nlin, (nlin - 1));
        avg_pix = npix - 1;
    }
    /*
     *  get the box extent in +, - direction from box size
     */
    pix_minus = (avg_pix - 1) / 2;
    pix_plus = avg_pix / 2;

    lin_minus = (avg_lin - 1) / 2;
    lin_plus = avg_lin / 2;
    /*
     *  loop over the lines and pixels and only smooth the non-zero (data) points
     */
    for (iln = 0; iln < nlin; iln++) {
        for (ipx = 0; ipx < npix; ipx++) {
            if (*(inarr + ipx + iln * npix) > 0) {
                /*
                 *  for this data pixel, average the surrounding box of values
                 */
                num = 0;
                avg = 0;
                for (jln = iln - lin_minus; jln <= iln + lin_plus; jln++) {
                    if (jln >= 0 && jln < nlin) {
                        for (jpx = ipx - pix_minus; jpx <= ipx + pix_plus; jpx++) {
                            /*
                             *  adjust the pixel to wrap around the seam
                             */
                            if (jpx < 0)
                                jpxm = npix + jpx;
                            else if (jpx >= npix)
                                jpxm = jpx - npix;
                            else
                                jpxm = jpx;

                            /*
                             *  sum the pixel if it is > 0
                             */
                            if (*(inarr + jpxm + jln * npix) > 0) {
                                num++;
                                avg += *(inarr + jpxm + jln * npix);
                            }
                        }
                    }
                }
                *(outarr + ipx + iln * npix) = avg / num;
            } else {
                *(outarr + ipx + iln * npix) = 0;
            }
        }
    }
    /*
     *  return the output
     */
    return;
}
