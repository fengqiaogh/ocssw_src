/* file: poly_fit.c */ 

#include "PGS_MODIS_35251.h"
#include "smfio.h"
#include "GEO_geo.h"
#include "GEO_util.h"

int GEO_poly_fit(
	double const	*coef,	
	double		const in_x,	
	int		const degree,
	double		* const out_fit
)
/*
!C*****************************************************************************
!Description:   subroutine in the utility functions of the Level-1A geolocation
                software to evaluate a polynomial.

!Input Parameters:
		double *coef - the array of coefficients for the polynomial
		double in_x - the input independent variable x
		int degree - the degree of the polynomial

!Output Parameters:
		double *out_fit	- the output value

Call functions: modsmf(MODIS_X_MNEMONIC_STRING, "user message string", "function,
                GEO_poly_fit.c") - writes error status messages to log

!Revision History:
		$Log: GEO_poly_fit.c,v $
		Revision 6.1  2011/02/14 21:34:00  kuyper
		Corrected const qualification of parameter which points at input data.

		Revision 1.6  1997/07/21 16:24:34  kuyper
		Baselined Version 1

 * Revision 1.6  1997/03/26  18:10:52  fhliang
 * Initial revision of SDST delivery of GEO_poly_fit.c.
 *
		Revision 1.5  1997/02/13 19:52:19  kuyper
		Merged seed files.

		Revision 1.4  1996/07/24 21:09:07  kuyper
		Standardized order of #include files.
		Declared arguments const.

		Revision 1.3  1996/07/23 23:30:28  kuyper
		Inserted required '!'s in comments.
		Changed constants to double, to avoid conversions.

		Revision 1.2  1996/07/18 19:54:18  kuyper
		Included self-checking header file.
		Added other required header file.
		James Kuyper Jr (kuyper@ltpmail.gsfc.nasa.gov)

		10/10/95
		Tracey W. Holmes
		Added debug option. 

		6/30/95
		Tracey W. Holmes (holmes@modis-xl.gsfc.nasa.gov)
		Added SDP error messages.

		6/12/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Revised prolog and initialized internal variables
                
		4/5/95
		Ruiming Chen
		Finished coding

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/
{

  double x_power_degree = 1.0;  /* degree power of x */
  int i = 0;  /* iteration parameter */

  if (degree <= 0)  {
    /* write to SDP event message */
    modsmf(MODIS_E_UTIL_POLY, "", "GEO_poly_fit, GEO_poly_fit.c");

    return FAIL;
  }

  if (degree > MAX_POLY_DEGREE)  {
    /* write to SDP event message */
    modsmf(MODIS_E_UTIL_MAX_POLY, "", "GEO_poly_fit, GEO_poly_fit.c");

    return FAIL;
  }

  /* initialize out_fit */
  *out_fit = 0.0;

  /* calculate the fit number */	
  for (i = 0; i <= degree; i++)  {
    *out_fit += coef[i]*x_power_degree;	
    x_power_degree = x_power_degree * in_x;
  }

  return SUCCESS;
}

