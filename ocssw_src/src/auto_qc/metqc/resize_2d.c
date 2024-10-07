int resize_2d(float *in_arr, int inpix, int inlin,
        int outpix, int outlin, float *out_arr)
/*******************************************************************

   resize_2d

   purpose: resize a 2d float array to a new size.  all arrays must
            be pre-allocated

   Returns type: int - 0 if all went well, 
                       -1 otherwise

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      float *           in_arr           I      input array size [npix, nlin]
      int               inpix            I      # input pixels (columns)
      int               inlin            I      # input lines (rows)
      int               outnpix          I      # output pixels (columns)
      int               outnlin          I      # output lines (rows)
      float *           out_arr          O      final output array

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       6-Feb-1997      Original development

 *******************************************************************/
 {
    float ratio_l, ratio_p;
    int op, oln, iln, ipx, ip, opx;
    /*
     *  do some quick checks
     */
    if (inpix <= 0 || inlin <= 0 || outpix <= 0 || outlin <= 0) {
        return -1;
    }
    /*
     *  find the ratio between the input and output size in lines and pixels
     *  (the -1 will convert the 0 origin indicies correctly)
     */
    ratio_l = (float) (inlin - 1) / (float) (outlin - 1);
    ratio_p = (float) (inpix - 1) / (float) (outpix - 1);

    /*
     *  loop thru the file and assign the closest input pixel to each output
     *  the input line, pixel are computed as the input * ratio.  .5 will 
     *  make the truncation go to the closest value
     */
    op = 0;
    for (oln = 0; oln < outlin; oln++) {
        iln = (int) ((float) oln * ratio_l + .5);

        for (opx = 0; opx < outpix; opx++) {
            ipx = (int) ((float) opx * ratio_p + .5);
            ip = ipx + iln * inpix;

            /*  out_arr[op++] = in_arr[ip];  */

            *(out_arr + op++) = *(in_arr + ip);
        }
    }
    return 0;
}
