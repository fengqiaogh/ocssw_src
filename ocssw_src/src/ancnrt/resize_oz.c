void resize_oz(short *input, int in_pix, int in_lin,
        int out_pix, int out_lin, short *output)
/*********************************************************************

   Function: resize_oz

   Returns type: void

   Description:  convert a short 2D array of values (ozone but it could
       be anything) into a larger / smaller short array.  expand the
       pixels / lines by whatever factors are needed

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      short *           input            I      input array
      int               in_pix            I      input # pixels
      int               in_lin            I      input # lines
      int               out_pix           I      output # pixels
      int               out_lin           I      output # lines
      short *           output           O      finished array

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       9 Jul 98        Original development

 *********************************************************************/ {
    float fact_lin, fact_pix, wt_p, wt_l, v1, v2, pix_pos, lin_pos;
    int ipx, iln, lin_index, pix_index;

    /*
     *  find the ratio between output and input pixels and lines
     */

    fact_pix = (float) (out_pix - 1) / (float) (in_pix - 1);
    fact_lin = (float) (out_lin - 1) / (float) (in_lin - 1);

    /*
     *  expand the lixels and lines using interpolation
     *  start with a line loop
     */
    for (iln = 0; iln < out_lin; iln++) {
        /*
         * get the line # of the first interpolation point and the weight 
         * to apply
         */
        lin_pos = (float) iln / fact_lin;
        lin_index = (short) lin_pos;
        wt_l = 1. - ((lin_pos - (float) lin_index) * fact_lin);

        /*
         *  loop over the pixels also
         */
        for (ipx = 0; ipx < out_pix; ipx++) {
            /*
             * get the pixel # of the first interpolation point and the weight 
             * to apply
             */
            pix_pos = (float) ipx / fact_pix;
            pix_index = (short) pix_pos;
            wt_p = 1. - ((pix_pos - (float) pix_index) * fact_pix);
            /*
             *  interpolate first in line direction.  Note that at end of rows 
             *  and columns, the next interpolate point is beyond the array
             *  bounds, hence the ifs below
             */
            v1 = (float) *(input + in_pix * lin_index + pix_index) * wt_l;
            if (wt_l < 1.) {
                v1 += (float) *(input + in_pix * (lin_index + 1) + pix_index) *
                        (1. - wt_l);
            }
            if (wt_p < 1.) {
                v2 = (float) *(input + in_pix * lin_index + (pix_index + 1))
                        * wt_l;
                if (wt_l < 1.) {
                    v2 += (float) *(
                            input + in_pix * (lin_index + 1) + (pix_index + 1)
                            ) * (1. - wt_l);
                }
            } else
                v2 = 0.;
            /*
             *  interpolate in pixel now
             */
            *(output + out_pix * iln + ipx) = wt_p * v1 + (1. - wt_p) * v2;
        } /* end pixel loop */
    } /* end line loop */

}
