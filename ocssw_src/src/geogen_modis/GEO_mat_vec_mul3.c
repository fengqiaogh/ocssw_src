/* file: GEO_mat_vec_mul3.c */

#include "GEO_basic.h"
#include "GEO_util.h"

void GEO_mat_vec_mul3(
	double mat[3][3],
	double vec1[3],
	double vec[3]
	)
{
/*
!C*****************************************************************************
!Description:   subroutine in the utility functions of the Level-1A geolocation
               software to perform multiplication of a matrix and a vector.

!Input Parameters:
                double mat[3][3] - the 3x3 matrix 
                double vec1[3] - the vector (with 3 elements)

!Output Parameters:
                double vec[3] - output vector with 3 elements

!Revision History:
		$Log: GEO_mat_vec_mul3.c,v $
		Revision 1.4  1997/07/21 16:24:34  kuyper
		Baselined Version 1

 * Revision 1.4  1997/03/26  18:10:14  fhliang
 * Initial revision of SDST delivery of GEO_mat_vec_mul3.c.
 *
		Revision 1.3  1996/07/23 23:28:50  kuyper
		Inserted required '!'s in comments.
		Changed to return void, not int.
		Changed constants to double, to avoid conversions.

		Revision 1.2  1996/07/18 19:49:36  kuyper
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
  int j = 0;  /* iteration parameter */

  /* multiplication of a matrix mat[3][3] and a vector vec1[3] */	
  for (i = 0; i < 3; ++i) {
    vec[i] = 0.0;
    for (j = 0; j < 3; ++j)
      vec[i] += mat[i][j] * vec1[j];	
  }

  return;
}
