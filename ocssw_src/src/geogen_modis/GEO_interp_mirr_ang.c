/* file: GEO_interp_mirr_ang.c */

#include "PGS_MODIS_35251.h"
#include "smfio.h"
#include "GEO_inst.h"
#include "GEO_parameters.h"
#include "GEO_util.h"

int GEO_interp_mirr_ang(
	double const sample_enc,
	double * const sample_mirr_ang	
	)
/*
!C*****************************************************************************
!Description:   Routine in Instrument Model of the Level-1A geolocation
                software to calculate the sample mirror angle as a polynomial
		function of the input encoder value.

!Input Parameters:
		int sample_number the sample number in the ideal band
		double sample_enc - sample encoder value
		double *sample_mirr_ang - sample mirror_angle

!Output Parameters:
		double *sample_mirr_ang - sample mirror angle

Return parameter:
                int err_stat - error status

Global variables:
		int poly_degree - degree of the mirror polynomial
		double poly_coef[MAX_POLY_DEGREE+1] - coefficient of the 
		  polynomial

Call functions:
		int GEO_poly_fit(double, double, int, double) - 
		  evaluate a polynomial
                modsmf(MODIS_X_MNEMONIC_STRING, "user message string", "function,
                GEO_interp_mirr_ang.c") - writes error status messages to log  

!Revision History:
		$Log: GEO_interp_mirr_ang.c,v $
		Revision 1.6  1997/07/21 16:24:34  kuyper
		Baselined Version 1

 * Revision 1.6  1997/03/26  18:07:54  fhliang
 * Initial revision of SDST delivery of GEO_interp_mirr_ang.c.
 *
		Revision 1.5  1997/02/13 20:10:48  kuyper
		Merged seed files.

		Revision 1.4  1996/07/24 21:04:56  kuyper
		Standardized order of #include files.
		Declared arguments const.

		Revision 1.3  1996/07/23 23:17:46  kuyper
		Inserted required '!' in comments.
		Removed ret_val.

		Revision 1.2  1996/07/18 19:20:39  kuyper
		Included self-checking header file.
		Replaced extern declarations with corresponding header file.
		Added needed header file.
		James Kuyper Jr. (kuyper@ltpmail.gsfc.nasa.gov)

		10/10/95
		Tracey W. Holmes
		Added debug option.

		6/30/95
		Tracey W. Holmes (holmes@modis-xl.gsfc.nasa.gov)
		Added SDP error messages. 

		6/12/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Fixed wrap-around lines and added.

		4/25/95
		Ruiming Chen
		Finished coding.

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/
{

  /* fit the polynomial to get the mirr angle for the sample pixel */
  if (GEO_poly_fit(poly_coef, sample_enc, poly_degree, sample_mirr_ang) == FAIL)
  {
    /* call SDP function to report error */
    modsmf(MODIS_E_UTIL_POLY, "", "GEO_poly_fit, GEO_interp_mirr_ang.c");

    return FAIL;
  }

  return SUCCESS;
}	     		

