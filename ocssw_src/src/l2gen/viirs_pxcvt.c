#include "l12_proto.h"
/*
 *  viirs_pxcvt.c contains routines for conversion of unaggregated to 
 *  aggregated VIIRS pixels and visa versa.  Also, it finds an 
 *  aggregated pixel that is x unaggregated samples away.
 */
static int ag_px_st[] = {0, 640, 1008, 2192, 2560, 3200};
static int ag_fact[] = {1, 2, 3, 2, 1, -1};
static int uag_px_st[] = {0, 640, 1376, 4928, 5664, 6304};

void viirs_pxcvt_2uag(int in_pix, int *out_pix, int *nag)
/*******************************************************************

   viirs_pxcvt_2uag

   purpose: convert pixel number from aggregated to unaggregated

   Returns type: void

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               in_pix           I      aggregated pixel number
      int *             out_pix          O      unaggregated pixel number
      int *             nag              O      number of samples aggregated
                                                at this agregated pixel

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       15-Aug-2011     Original development

  Note that the mapping from aggregated to unaggregated makes a range of 
  pixel numbers from out_pix -> out_pix + nag - 1

 *******************************************************************/ {
    int found, nzone = 6, iz;
    /*
     *  Find the zone for the pixel and convert it to the unaggregated value
     */
    found = 0;
    for (iz = 0; iz < (nzone - 1); iz++) {
        if ((in_pix >= ag_px_st[iz]) && (in_pix < ag_px_st[ iz + 1 ])) {
            *out_pix = uag_px_st[iz] + (in_pix - ag_px_st[iz]) * ag_fact[iz];
            *nag = ag_fact[iz];
            found = 1;
            break;
        }
    }
    /*
     *  extrapolate for values outside actual range
     */
    if (found != 1) {
        *nag = 1;
        if (in_pix < ag_px_st[0])
            *out_pix = in_pix;
        else
            *out_pix = uag_px_st[ nzone - 1 ] + in_pix - ag_px_st[ nzone - 1 ];
    }
}

void viirs_pxcvt_2ag(int in_pix, int *out_pix)
/*******************************************************************

   viirs_pxcvt_2ag

   purpose: convert pixel number from unaggregated to aggregated

   Returns type: void

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               in_pix           I      unaggregated pixel number
      int *             out_pix          O      aggregated pixel number

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       15-Aug-2011     Original development

 *******************************************************************/ {
    int found, nzone = 6, iz;
    /*
     *  Find the zone for the pixel and convert it to the aggregated value
     */
    found = 0;
    for (iz = 0; iz < (nzone - 1); iz++) {
        if ((in_pix >= uag_px_st[iz]) && (in_pix < uag_px_st[ iz + 1 ])) {
            *out_pix = ag_px_st[iz] + (in_pix - uag_px_st[iz]) / ag_fact[iz];
            found = 1;
            break;
        }
    }
    /*
     *  extrapolate for values outside actual range
     */
    if (found != 1) {
        if (in_pix < uag_px_st[0])
            *out_pix = in_pix;
        else
            *out_pix = ag_px_st[ nzone - 1 ] + in_pix - uag_px_st[ nzone - 1 ];
    }
}

void viirs_pxcvt_agdel(int in_pix, int del, int *out_pix)
/*******************************************************************
  
   viirs_pxcvt_agdel

   purpose: add unaggregated samples to aggregated pixel numbers

   Returns type: void

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int               in_pix           I      aggregated pixel number
      int               del              I      unaggregated sample offset
      int *             out_pix          O      final aggregated pixel number

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       15-Aug-2011     Original development

 *******************************************************************/ {
    int uag_px, ag;
    /*
     *  get the unaggregated pixel
     */
    viirs_pxcvt_2uag(in_pix, &uag_px, &ag);
    /*
     *  add the offset, if in (+) direction, make sure it is from highest 
     *  unaggregated pixel # possible
     */
    uag_px = (del < 0) ? uag_px + del : uag_px + (ag - 1) + del;
    /*
     *  convert back to aggregated and done
     */
    viirs_pxcvt_2ag(uag_px, out_pix);
}
