#include "GEO_basic.h"
#include "GEO_util.h"

int GEO_poly_coef1(
	double const var_val1,
	double const var_val2,
	double const derivative1,
	double const derivative2,
	double const t1,
	double const t2,
	double fit_coef[CUBIC])

/*
!C*****************************************************************************
!Description:   
		subroutine in utility group of the Level-1A geolocation
                software to generate cubic polynomial coefficients using
		values and slopes at two times.

!Input Parameters:
                double var_val1 - value of first point
	        double var_val2 - value of second point
	        double derivative1 - slope of first point
	        double derivative2 - slope of second point
	        double t1 - time of first point
	        double t2 - time of second point
                int scan_number - the scan number

!Output Parameters:
                double fit_coef[CUBIC] - cubic polynomial coefficients

Return parameter:
                int err_stat - error status

Global variables:
	        None

Call functions:
		None

!Revision History:
 * $Log: GEO_poly_coef1.c,v $
 * Revision 2.1  1997/10/21 18:16:22  kuyper
 * Returned from ClearCase
 *
 * Revision /main/GEO_V2_DEV/1 1997/10/15 kuyper
 * Changed to use CUBIC macro.
 *
 * Revision 1.5  1997/07/21  16:24:34  kuyper
 * Baselined Version 1
 *
 * Revision 1.5  1997/03/26  18:10:35  fhliang
 * Initial revision of SDST delivery of GEO_poly_coef1.c.
 *
		Revision 1.4  1996/07/24 21:08:47  kuyper
		Declared arguments const.

		Revision 1.3  1996/07/23 23:29:41  kuyper
		Inserted required '!'s in comments.
		Changed constants to double, to avoid conversions.

		Revision 1.2  1996/07/18 19:52:42  kuyper
		Included self-checking header file.


		6/2/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Finished coding.

		6/12/95
		Frederick S. Patt
		Revised prolog

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/

{
/* Local variable declarations */

  double dt = 1.0;

  /* Begin program logic */

  dt = t2 - t1;

  /* calculate the polynominal coefficients */
  fit_coef[0] = var_val1;
  fit_coef[1] = derivative1*dt;
  fit_coef[2] = 3.0*var_val2 - 3.0*var_val1 - 2.0*derivative1*dt -
	  derivative2*dt;
  fit_coef[3] = 2.0*var_val1 - 2.0*var_val2 + derivative1*dt + derivative2*dt;

  return SUCCESS;
} 
