#include "get_ctht.h"

/*
axis_interp - emulate the idl interpolate function with all its drawbacks

   Returns - float - the interpolate on the axis values

   Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     float *           axis_arr_vals    I      the axis values
     int32_t           n_ax             I      # of values in axis
     float             ind_frac         I      axis index and fraction to the
                                                 next index

   W. Robinson, SAIC, 8 Sep 2023
 */

float axis_interp( float *axis_arr_vals, int32_t n_ax, float ind_frac )
  {
  float retval, fr;
  int32_t indx;

  // first, just return min or max if we're outside
  if( ind_frac <= 0 )
    {
    retval = axis_arr_vals[0];
    }
  else if( ind_frac >= n_ax - 1 )
    {
    retval = (float) n_ax - 1;
    retval = axis_arr_vals[ n_ax - 1 ];
    }
  else
    {
    indx = (int32_t) ind_frac;
    fr = ind_frac - indx;
    retval = fr * axis_arr_vals[ indx + 1 ] + ( 1. - fr ) * axis_arr_vals[indx];
    }
  return retval;
  }

