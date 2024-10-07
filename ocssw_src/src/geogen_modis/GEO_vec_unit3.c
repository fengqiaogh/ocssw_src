/* file: GEO_vec_unit3.c */

#include <float.h>
#include "PGS_MODIS_35251.h"
#include "smfio.h"
#include "GEO_basic.h"
#include "GEO_util.h"

int  GEO_vec_unit3(
	double vec[3],
	double unit_vec[3]
        )
/*
!C*****************************************************************************
!Description:   subroutine in the utility functions of the Level-1A geolocation
                software to calculate the unit vector of a 3 element vector.

!Input Parameters:
		double vec[3] - the input 3-vector 

!Output Parameters:
	        double unit_vec[3] - the unit vector from the input 3-vector

Call functions: modsmf(MODIS_X_MNEMONIC_STRING, "user message string", "function,
                GEO_vec_unit3.c") - writes error status messages to log

!Revision History:
		$Log: GEO_vec_unit3.c,v $
		Revision 1.6  1997/07/21 16:24:34  kuyper
		Baselined Version 1

 * Revision 1.6  1997/03/26  18:17:18  fhliang
 * Initial revision of SDST delivery of GEO_vec_unit3.c.
 *
		Revision 1.5  1997/03/10 17:30:41  kuyper
		Changed to inequality test on vec_length.

		Revision 1.4  1997/02/13 19:59:52  kuyper
		Merged seed files
		Added RCS Log keyword.

		James Kuyper Jr. (kuyper@ltpmail.gsfc.nasa.gov)

		6/30/95
		Tracey W. Holmes (holmes@modis-xl.gfsc.nasa.gov)
		Added SDP error messages.

		6/12/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Revised prolog.

		4/5/95
		Ruiming Chen (rchen@ltpmail.gsfc.nasa.gov)
		Finished coding.

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/
{
  double vec_length = 0.0;  /* vector length */
  int i = 0;  /* iteration parameter */

  GEO_vec_length3(vec, &vec_length);
  if (vec_length < DBL_MIN) {
    /* write to SDP event message */
    modsmf(MODIS_E_UTIL_VEC, "", "GEO_vec_unit3, GEO_vec_unit3"); 

    return FAIL;
  }

  for (i = 0; i < 3; ++i)
    unit_vec[i] = vec[i] / vec_length;

  return SUCCESS;
}
