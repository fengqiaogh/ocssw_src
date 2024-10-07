#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "fmt_check.h"

int var_decode(char *str, int32 hdf_type, void *array, int32 index, int32 adj_flg)
/*******************************************************************

   var_decode

   purpose: decode a value in a string to the proper hdf type specified

   Returns type: int - 0 if all went well, -1 if unable to read from the
          string, -2 if value read does not fit in storage type, -3 if
          hdf_type is not a valid one (or one handled currently)

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            str              I      variable in the string to decode
      int32             hdf_type         I      type of parameter to put data 
                                                into (DFNT_...)
      void *            array           I/O     array to place decoded value
      int32             index            I      position in array to place value
                                                (not used for char output)
      int32             adj_flg          I      limit adjust flag.  If 0, do not adjust
                                                the variable decoded.  If 1, consider 
                                                index = 0 to mean low
                                                limit and index = 1 to mean high limit
                                                and adjust limits for rounding error

   Note on adj_flg. var_decode is used to get the actual values of attributes
   that should come in.  It is also used to get the valid range for sds elements.
   In this case, it was found that the sds values could exceed the limits by small
   amounts.  The adjustment will be to expand the range limits by the uncertainty
   in the float values of the limits.  In other words 

    low limit = (input low limit) - (input low limit) * 10 ** ( - # digits accuracy )
    hi  limit = (input hi  limit) + (input hi  limit) * 10 ** ( - # digits accuracy )

   This is only done for float types

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       1-Oct-1996      Original development
      W. Robinson       31-Oct-1996     upgrade to handle float64 data type
      W. Robinson       18-Feb-1997     modify range limits in float32 and 64
                                        to include excess due to rounding error

 *******************************************************************/ {
    int32 i32val, *ptri32;
    float32 f32val, *ptrf32;
    float64 f64val, *ptrf64;
    int8 i8val, *ptri8;
    uint8 ui8val, *ptrui8;
    int16 i16val, *ptri16;
    char *ptrc;

    /*
     *  be able to point to the array in all these ways
     */
    ptri32 = (int32 *) array;
    ptrf32 = (float32 *) array;
    ptrf64 = (float64 *) array;
    ptri8 = (int8 *) array;
    ptrc = (char *) array;
    ptrui8 = (uint8 *) array;
    ptri16 = (int16 *) array;

    /*
     *  do 32 bit integer
     */
    if (hdf_type == DFNT_INT32) {
        if (sscanf(str, "%d", &i32val) == EOF) {
            return -1;
        }
        *(ptri32 + index) = i32val;
    }
        /*
         *  do 32 bit float
         */
    else if (hdf_type == DFNT_FLOAT32) {
        if (sscanf(str, "%f", &f32val) == EOF) {
            return -1;
        }
        /*
         *  for low limit, expand it lower by the value * 10 ** -( # bits of significance),
         *  for high limit, expand it higher in same way
         */
        if (adj_flg != 0)
            f32val = (index == 0) ?
            f32val - fabs(f32val) * pow(10., -(FLT_DIG)) :
            f32val + fabs(f32val) * pow(10., -(FLT_DIG));
        *(ptrf32 + index) = f32val;
    }
        /*
         *  do 64 bit float
         */
    else if (hdf_type == DFNT_FLOAT64) {
        if (sscanf(str, "%lf", &f64val) == EOF) {
            return -1;
        }
        /*
         *  for low limit, expand it lower by the value * 10 ** -( # bits of significance),
         *  for high limit, expand it higher in same way
         */
        if (adj_flg != 0)
            f64val = (index == 0) ?
            f64val - fabs(f64val) * pow(10., -(DBL_DIG)) :
            f64val + fabs(f64val) * pow(10., -(DBL_DIG));
        *(ptrf64 + index) = f64val;
    }        /*
  *  do 8 bit integer
  */
    else if (hdf_type == DFNT_INT8) {
        if (sscanf(str, "%d", &i32val) == EOF) {
            return -1;
        }
        if (i32val > 127 || i32val < -128) {
            return -2;
        }
        i8val = i32val;
        *(ptri8 + index) = i8val;
    }
        /*
         *  do character
         */
    else if (hdf_type == DFNT_CHAR) {
        strcpy(ptrc, str);
    }
        /*
         *  do 8 bit unsigned integer
         */
    else if (hdf_type == DFNT_UINT8) {
        if (sscanf(str, "%d", &i32val) == EOF) {
            return -1;
        }
        if (i32val > 255 || i32val < 0) {
            return -2;
        }
        ui8val = i32val;
        *(ptrui8 + index) = ui8val;
    }
        /*
         *  do 16 bit integer
         */
    else if (hdf_type == DFNT_INT16) {
        if (sscanf(str, "%hd", &i16val) == EOF) {
            return -1;
        }
        *(ptri16 + index) = i16val;
    }
        /*
         *  for an unplanned hdf_type
         */
    else {
        return -3;
    }

    return 0;
}
