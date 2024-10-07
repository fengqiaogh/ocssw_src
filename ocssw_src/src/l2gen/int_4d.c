#include <stdlib.h>
#include <stdio.h>
#include "l12_proto.h"
#include "l2prod.h"
#include "l1.h"
#include "get_ctht.h"

/*
  file int_4d.c, get index and weights to do a 4-d interpolation
  Also for now, have the 3-d interp hard coded as int_3d
*/
int32_t int_4d( float *sxax[4], int32_t snax[4], int64_t n_block, float x[4], 
  int64_t ar_off[16], float ar_wt[16], float orig_wt[4] )

 /*
    int_4d - get index and weights to do a 4-d interpolation

    Returns, int32_t status 0 is good 1 - outside range

   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     float *[4]        sxax             I      array of pointers to each 
                                               dimension axis
     int32_t [4]       snax             I      size of each axis 
     int64_t           n_block          I      size of all dimensions before
                                               the ones here being used
     float [4]         x                I      value of point on each axis
     int64_t [16]      ar_off           O      linear offset to each corner pt
     float [16]        ar_wt            O      weight for each corner pt
     float [4]         orig_wt          O      the weights from the low point

   the multi-dim array, sxax has dims [ snax[0],... snax[3] ] with a the 
   combination of the fastest dimensions Based on x's position in each 
   sxax axis, there will be a index in each axis of the start point for 
   the interpolation and also a fraction of the distance to the next 
   point.  That is converted into a linear position and multiplied by a 
   to properly index the data for the offset and to a combined weight 
   to apply to all 16 corner points

  W. Robinson, SAIC, 3 Aug 2023

*/
  {
  int32_t ndim = 4, base_ix[4];
  int32_t idim, idim0, idim1, idim2, idim3;
  int32_t int_off0, int_off1, int_off2, int_off3;
  int64_t base_off, ixoff;
  float int_wt0, int_wt1, int_wt2, int_wt3;
  float base, frac_ix[4], delta[4], bounds[2];

  for( idim = 0; idim < ndim; idim++ )
    {
    //  check for bounds violation
    if( sxax[idim][0] < sxax[idim][snax[idim]-1] )
      {
      bounds[0] = sxax[idim][0];
      bounds[1] = sxax[idim][snax[idim]-1];
      }
    else
      {
      bounds[0] = sxax[idim][snax[idim]-1];
      bounds[1] = sxax[idim][0];
      }
    if( ( x[idim] < bounds[0] ) || ( x[idim] > bounds[1] ) )
      {
      //printf( "%s - %d: I: a point is outside table bounds\n", __FILE__, 
      // __LINE__ );
      return 1;
      }
    delta[idim] = sxax[idim][1] - sxax[idim][0];
    base = ( x[idim] - sxax[idim][0] ) / delta[idim];
    // *** for Andy's sol zen, it starts at 0.001 so I'll make it 0
    //  to match what he does FOR NOW
    if( idim == 1 ) {
      delta[idim] = sxax[idim][2] - sxax[idim][1];
      base = x[idim] / delta[idim];
    }
    base_ix[idim] = (int32_t) base;
    frac_ix[idim] = fmod( (double)base, 1. );
    orig_wt[idim] = frac_ix[idim];
    }
  //  get offset from 0 to base indicies
  base_off = base_ix[0] + snax[0] * 
    ( base_ix[1] + snax[1] * ( base_ix[2] + snax[2] * base_ix[3] ) );
 /*
  *  compile the locations of corners and weights here 
  */
  ixoff = 0;
  for( idim0 = 0; idim0 < 2; idim0++ )
    {
    int_wt0 = ( idim0 == 0 ) ? 1 - frac_ix[0] : frac_ix[0];
    int_off0 = ( idim0 == 0 ) ? 0 : 1;
    for( idim1 = 0; idim1 < 2; idim1++ )
      {
      int_wt1 = ( idim1 == 0 ) ? 1 - frac_ix[1] : frac_ix[1];
      int_off1 = int_off0 + ( ( idim1 == 0 ) ? 0 : snax[0] );
      for( idim2 = 0; idim2 < 2; idim2++ )
        {
        int_wt2 = ( idim2 == 0 ) ? 1 - frac_ix[2] : frac_ix[2];
        int_off2 = int_off1 + ( ( idim2 == 0 ) ? 0 : snax[0] * snax[1] );
        for( idim3 = 0; idim3 < 2; idim3++ )
          {
          int_wt3 = ( idim3 == 0 ) ? 1 - frac_ix[3] : frac_ix[3];
          int_off3 = int_off2 + ( ( idim3 == 0 ) ? 
            0 : snax[0] * snax[1] * snax[2] );
 
          ar_off[ixoff] = n_block * ( base_off + int_off3 );
          ar_wt[ixoff++] = int_wt3 * int_wt2 * int_wt1 * int_wt0;
          }
        }
      }
    }
  return 0;
  }

int32_t iint_3d( int32_t snax[3], double x[3], int64_t ar_off[9], 
  float ar_wt[9], float orig_wt[3] )

 /*
    iint_3d - get index and weights to do a 3-d interpolation.  this is done 
      assuming that the input x values are in units of each axis with the
      fractional part as the fraction of the span from 1 grid point to the next
      (I'm not sure of the universality of this but I think it will reproduce
      what Andy is getting)

    Returns, int32_t status 0 is good 1 - outside range

   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     int32_t [3]       snax             I      size of each axis
     double [3]        x                I      value of point on each axis
     int64_t [9]       ar_off           O      linear offset to each corner pt
     float [9]         ar_wt            O      weight for each corner pt
     float [3]         orig_wt          O      the weights from the low point

  W. Robinson, SAIC, 22 Aug 2023

*/
  {
  int32_t ndim = 3, base_ix[3];
  int32_t idim, idim0, idim1, idim2;
  int32_t int_off0, int_off1, int_off2;
  int64_t base_off, ixoff;
  float int_wt0, int_wt1, int_wt2;
  float base, frac_ix[3], bounds[2];

  for( idim = 0; idim < ndim; idim++ )
    {
    bounds[0] = 0;
    bounds[1] = snax[idim] - 1;
    //  check for bounds violation
    if( ( x[idim] < bounds[0] ) || ( x[idim] > bounds[1] ) )
      {
      // printf( "%s - %d: E: a point is outside table bounds\n", __FILE__,
      // __LINE__ );
      // exit(1);
      // quick fix is to move point in-bounds
      if( x[idim] < bounds[0] ) x[idim] = bounds[0];
      if( x[idim] > bounds[1] ) x[idim] = bounds[1];
      }
    //  steps are 1, so delta = 1
    base = x[idim];
    base_ix[idim] = (int32_t) base;
    frac_ix[idim] = fmod( (double)base, 1. );
    orig_wt[idim] = frac_ix[idim];
    }
  //  get offset from 0 to base indicies
  base_off = base_ix[0] + snax[0] * ( base_ix[1] + snax[1] * base_ix[2] );
 /*
  *  compile the locations of corners and weights here
  */
  ixoff = 0;
  for( idim0 = 0; idim0 < 2; idim0++ )
    {
    int_wt0 = ( idim0 == 0 ) ? 1 - frac_ix[0] : frac_ix[0];
    int_off0 = ( ( idim0 == 0 ) || ( base_ix[0] == snax[0] - 1 ) ) ? 0 : 1;
    for( idim1 = 0; idim1 < 2; idim1++ )
      {
      int_wt1 = ( idim1 == 0 ) ? 1 - frac_ix[1] : frac_ix[1];
      int_off1 = int_off0 + 
        ( ( ( idim1 == 0 ) || ( base_ix[1] == snax[1] - 1 ) ) ? 0 : snax[0] );
      for( idim2 = 0; idim2 < 2; idim2++ )
        {
        int_wt2 = ( idim2 == 0 ) ? 1 - frac_ix[2] : frac_ix[2];
        int_off2 = int_off1 + 
          ( ( ( idim2 == 0 ) || ( base_ix[2] == snax[2] - 1 ) ) ? 
          0 : snax[0] * snax[1] );
        ar_off[ixoff] = ( base_off + int_off2 );
        ar_wt[ixoff++] = int_wt2 * int_wt1 * int_wt0;
        }
      }
    }
  return 0;
  }
