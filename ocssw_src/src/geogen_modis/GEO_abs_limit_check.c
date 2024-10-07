/* file: GEO_abs_limit_check.c */

#include "GEO_geo.h"
#include "GEO_validation.h"

int GEO_abs_limit_check(
	double data_samples[],
	int const number_of_samples,
	double data_limits[2],
	int sample_flags[]
	)
{
/*
!C*****************************************************************************
!Description:   
		Routine in Input group of the Level-1A geolocation
                software to validate a set of data samples against 
		absolute limits.  Flags set to BAD_DATA for all samples
		not within limits.

!Input Parameters:
                double data_samples[] - array of input data samples
		int number_of_samples - number of samples to check
		double data_limits[2] - lower and upper limits

!Output Parameters:
                int sample_flags[] - array of validation flags

Return parameter:
                int err_stat - error status

Global variables:
		None

Call functions:
		None

!Revision History:
$Log: GEO_abs_limit_check.c,v $
Revision 1.6  1997/07/21 16:24:34  kuyper
Baselined Version 1

 * Revision 1.6  1997/03/26  17:33:13  fhliang
 * Initial revision of SDST delivery of GEO_abs_limit_check.c.
 *
Revision 1.5  1996/07/29 14:45:51  kuyper
Used BAD_DATA for flags, rather than FAIL (same value).

 * Revision 1.4  1996/07/24  20:53:10  kuyper
 * Declared number_of_samples const.
 *
Revision 1.3  1996/07/23 22:21:52  kuyper
Inserted required '!' in comments.

Revision 1.2  1996/07/18 15:28:46  kuyper
Added self-checking include file.

                6/14/95
                Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
                Finished coding.

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/

  /* Local variable declarations */

  int num_valid = 0; /* number of valid data samples found */
  int i = 0; /* loop index */

  /* Begin program logic */
  /* Loop through samples */
  for (i = 0; i < number_of_samples; i++) {

    /* Check samples against limits */
    if ((data_samples[i] < data_limits[0]) || 
	(data_samples[i] > data_limits[1])) {

      sample_flags[i] = BAD_DATA; 
    }
    else {
      num_valid += 1;
    }
  }

  /* If no valid samples, return FAIL */
  if (num_valid == 0) 
    return FAIL;

  return SUCCESS;
}
