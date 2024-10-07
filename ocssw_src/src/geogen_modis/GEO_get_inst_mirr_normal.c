/* file: GEO_get_inst_mirr_normal.c */

#include "PGS_MODIS_35251.h"
#include "smfio.h"
#include "GEO_inst.h"
#include "GEO_parameters.h"
#include "GEO_util.h"

int GEO_get_inst_mirr_normal(
	const int scan_number,
	const int sample_number,
	double n_inst_normal[3]
        )
/*
!C*****************************************************************************
!Description:	
		Routine in Instrument Model of the Level-1A geolocation
                software to calculate the mirror normal vector for a sample.
		It calls functions to interpolate the mirror encoder data 
		to the sample and convert the encoder value to an angle.  
		It then computes the mirror normal vector using the mirror
		wedge angles (alpha, beta, gammaa) and axis error, and 
                transforms the vector to the instrument coordinate frame.
		Reference Geolocation ATBD pp. 34-49.
		(Note that spelling of gammaa avoids confusion with gamma()
		function of math library.)

!Input Parameters:
                int scan_number - the scan number
                int sample_number - the sample number in the ideal band

!Output Parameters:
		double n_inst_normal[3] - the unit mirr normal vector in 
		                         instrument coordinates

Return parameter:
                int  err_stat - error status

Global variables:
		double alpha - mirr wedge error alpha
		double beta - mirr wedge error beta
		double gammaa - mirr wedge error gamma
		double mirr_side1_range[2] - range of side 1 mirror angles
		double T_mirr2inst[3][3] - transformation matrix from mirror 
		 to instrument coord.

Call functions:
		int GEO_interp_mirr_enc(int, int, double) - interpolate mirror 
		  encoder to sample time
		int GEO_interp_mirr_ang(double, double)	- compute mirror angle
		  from encoder value
		int GEO_vec_unit3(double vec[3], vec_unit[3]) - get unit 
		  vector of a vec
		int GEO_mat_vec_mul3(double mat[3][3], double vec1[3],
		  double vec[3]) - Multiply a vector by a matrix.
                modsmf(MODIS_X_MNEMONIC_STRING, "user message string", "function,
                GEO_get_inst_mirr_normal.c") - writes error status messages to log

!Revision History:
		$Log: GEO_get_inst_mirr_normal.c,v $
		Revision 1.7  1997/07/21 16:24:34  kuyper
		Baselined Version 1

 * Revision 1.7  1997/03/26  18:05:45  fhliang
 * Initial revision of SDST delivery of GEO_get_inst_mirr_normal.c.
 *
		Revision 1.6  1997/02/13 19:36:08  kuyper
		Merged seed files.

		Revision 1.5  1996/09/23 19:34:56  kuyper
		Corrected sign of gammaa in n_normal_side[0] calculation.

		Revision 1.4  1996/07/24 21:01:57  kuyper
		Standardized order of #include files.
		Declared arguments const.

		Revision 1.3  1996/07/23 23:06:06  kuyper
		Inserted required '!' in comments.
		Removed ret_val, i, j.

		Revision 1.2  1996/07/18 16:35:12  kuyper
		Included self-checking header file.
		Replaced extern declarations with corresponding header files.
		Added required header files.
		Converted int constants to double, to avoid conversion lint warnings.
		James Kuyper kuyper@ltpmail.gsfc.nasa.gov


		4/25/95
		Ruiming Chen
		Finished coding.

		6/12/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Expanded description.
		 
		6/21/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Added mirror side angle range as a parameter.

		6/30/95
		Tracey W. Holmes (holmes@modis-xl.gsfc.nasa.gov)
		Added SDP error messages.

		10/10/95
		Tracey W. Holmes
		Added debug option. 
	

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/
{


/***************************************************************************
declare local variables
***************************************************************************/

  double sample_mirr_ang = 0.0;  /* mirr angle for the sample */
  double sample_enc = 0.0;  /* encoder number for the sample */
  double T_rot[3][3] = {0.0};  /* The mirror rotation matrix */ 
  double n_normal_side[3] = {0.0};  /* mirr normal side vector */
  double n_mirr_normal[3] = {0.0};  /* mirr normal */
  double n_inst[3] = {0.0};  /* mirr normal in the instrument coord. */ 


/***************************************************************************
calculation from here
***************************************************************************/

  /* calculate encoder position */
  if (GEO_interp_mirr_enc(scan_number, sample_number, &sample_enc) == FAIL) {
    /* call SDP function to report error */
    modsmf(MODIS_E_GEO_INT_MIRR_ENC, "", 
           "GEO_interp_mirr_enc, GEO_get_inst_mirr_normal.c");

    return FAIL;
    }

    /* compute the scan mirror angle */
    if (GEO_interp_mirr_ang(sample_enc, &sample_mirr_ang) == FAIL)
    {
      /* call SDP function to report error */
      modsmf(MODIS_E_GEO_INT_MIRR_ANG, "", 
             "GEO_interp_mirr_ang, GEO_get_inst_mirr_normal.c");

      return FAIL;
    }

    /* construct T_rot */
    T_rot[0][0] = 1.0;
    T_rot[1][1] = cos(sample_mirr_ang);
    T_rot[1][2] = - sin(sample_mirr_ang);
    T_rot[2][1] = sin(sample_mirr_ang);
    T_rot[2][2] = cos(sample_mirr_ang);

    /* construct n_normal_side */
    if (sample_mirr_ang > mirr_side1_range[0] && 
	sample_mirr_ang <= mirr_side1_range[1]) {
      n_normal_side[0] = - sin(beta/2.0 + gammaa);
      n_normal_side[1] = sin(alpha/2.0) * cos(beta/2.0 + gammaa);
      n_normal_side[2] = cos(alpha/2.0) * cos(beta/2.0 + gammaa);
    }

    else {
      n_normal_side[0] = - sin(beta/2.0 - gammaa);
      n_normal_side[1] = sin(alpha/2.0) * cos(beta/2.0 - gammaa);
      n_normal_side[2] = - cos(alpha/2.0) * cos(beta/2.0 - gammaa);
    }

    /* get n_mirr_normal */
    GEO_mat_vec_mul3(T_rot, n_normal_side, n_mirr_normal);

    /* get n_inst */
    GEO_mat_vec_mul3(T_mirr2inst, n_mirr_normal, n_inst);

    /* calculate n_inst_normal */
    if (GEO_vec_unit3(n_inst, n_inst_normal) == FAIL) {
      /* call SDP function to report error */
      modsmf(MODIS_E_GEO_VEC_UNIT3, "",  
             "GEO_vec_unit3, GEO_get_inst_mirr_normal.c");

      return FAIL;
    }

    return SUCCESS;
}
