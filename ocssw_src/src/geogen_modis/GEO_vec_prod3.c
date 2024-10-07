/* GEO_vec_prod3.c */

#include "GEO_basic.h"
#include "GEO_util.h"

void GEO_vec_prod3(
	double  vec1[3],
	double  vec2[3],
	double  * const prod
        )
{
/*
!C*****************************************************************************
!Description:   subroutine in the utility functions of the Level-1A geolocation
                software to calculate the products of two 3 element vectors.

!Input Parameters:
		double vec1[3] - the first 3 element vector
		doubel vec2[3] - the second 3 element vector

!Output Parameters:
		double *prod - the product of the two input 3 element vectors

!Revision History:
		$Log: GEO_vec_prod3.c,v $
		Revision 1.4  1997/07/21 16:24:34  kuyper
		Baselined Version 1

 * Revision 1.4  1997/03/26  18:17:00  fhliang
 * Initial revision of SDST delivery of GEO_vec_prod3.c.
 *
		Revision 1.3  1996/07/24 21:52:58  kuyper
		Changed to return void.
		Declared arguments const.
		Inserted required '!'s in comments.
		Converted constants to double, to skip implied conversion.
		Removed ret_val.

		Revision 1.2  1996/07/18 21:27:33  kuyper
		Included self-checking header file.


		4/5/95
		Ruiming Chen (rchen@ltpmail.gsfc.nasa.gov)
		Finished coding.

		6/12/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Revised prolog

Call function:
		None

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/

  int i = 0;  /* iteration parameter */

  /* calculate the product of vec1 and vec2 */	
  *prod = 0.0;
  for (i = 0; i < 3; ++i)
    *prod = *prod + vec1[i] * vec2[i];

  return;
}
