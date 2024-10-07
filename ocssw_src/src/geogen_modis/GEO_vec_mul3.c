/* file: GEO_vec_mul3.c */

#include "GEO_basic.h"
#include "GEO_util.h"

int  GEO_vec_mul3(
	double vec1[3],
	double vec2[3],
	double vec[3]
        )
{
/*
!C*****************************************************************************
!Description:   subroutine in the utility functions of the Level-1A geolocation
                software to calculate the cross product of two 3 element
		vectors.

!Input Parameters:
		double vec1[3]	the first vector of 3 elements
		double vec2[3]	the second vector of 3 elements

!Output Parameters:
		double vec[3]	the output vector of 3 elements

!Revision History:
		$Log: GEO_vec_mul3.c,v $
		Revision 1.4  1997/07/21 16:24:34  kuyper
		Baselined Version 1

 * Revision 1.4  1997/03/26  18:24:32  fhliang
 * Initial revision of SDST delivery of GEO_vec_mul3.c.
 *
		Revision 1.3  1996/07/24 21:51:10  kuyper
		Inserted required '!'s in comments.

		Revision 1.2  1996/07/18 21:26:27  kuyper
		Included self-checking header file.


		4/5/95
		Ruiming Chen (rchen@ltpmail.gsfc.nasa.gov)
		Finished coding.

		6/12/95
		Frederick S> Patt (patt@modis-xl.gsfc.nasa.gov)
		Revised prolog

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/

	vec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
        vec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
        vec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

	return SUCCESS;
}
