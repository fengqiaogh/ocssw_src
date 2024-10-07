/* file: GEO_vec_length3.c */

#include "GEO_basic.h"
#include "GEO_util.h"

void GEO_vec_length3(
	double vec[3],
	double * const length
        )
{
/*
!C*****************************************************************************
!Description:   subroutine in the utility functions of the Level-1A geolocation
                software to calculate the length of a 3 element vector.

!Input Parameters:
		double vec - a three element vector

!Output Parameters:
		double *length - the length of the three element vector

!Revision History:
		$Log: GEO_vec_length3.c,v $
		Revision 1.4  1997/07/21 16:24:34  kuyper
		Baselined Version 1

 * Revision 1.4  1997/03/26  18:16:39  fhliang
 * Initial revision of SDST delivery of GEO_vec_length3.c.
 *
		Revision 1.3  1996/07/24 21:47:10  kuyper
		Changed to return void.
		Declared arguments const.
		Inserted required '!'s in comments.
		Converted constants to double, to skip implied conversion.

		Revision 1.2  1996/07/18 21:25:35  kuyper
		Included self-checking header file.


		4/5/95
		Ruiming Chen
		Finished coding.

		6/12/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Revised prolog

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/

  int i = 0;  /* iteration parameter */

  *length = 0.0;
  for (i = 0; i < 3; ++i)
    *length = *length + vec[i] * vec[i];

  *length = sqrt(*length);


  return;
}
