#include "PGS_MODIS_35251.h"
#include "GEO_inst.h"
#include "GEO_parameters.h"
#include "GEO_product.h"
#include "GEO_util.h"
#include "smfio.h"

int GEO_get_view_vec(
	const int scan_number,
	const int sample_number,
	const int num_detectors,
	double u_inst[MAX_DETECTORS][3]
	)
/*
!C*****************************************************************************
!Description: 
       Routine in Instrument Model of the Level-1A geolocation software to
       calculate the view vectors for a sample.  This routine calls functions
       to compute the mirror normal vector.  It then reflects the telescope
       view vectors off of the mirror to determine the view vectors in the
       instrument frame. 

!Input Parameters:
       scan_number                     The scan number
       sample_number                   The sample number
       num_detectors                   The number of detectors in the sample

!Output Parameters:
       u_inst                          View vectors in instrument coordinates

Return Values:
       FAIL                            If any GEO_get_inst_mirr_normal() fails
       SUCCESS                         Otherwise

Externally Defined:
        FAIL                            "GEO_basic.h"
        MODIS_E_GEO                     "PGS_MODIS_35251.h"
        SUCCESS                         "GEO_basic.h"
        T_tel2inst                      "GEO_parameters.h"
        u_tel                           "GEO_parameters.h"
 
Called by:
        GEO_earth_location
 
Routines Called:
        GEO_get_inst_mirr_normal        "GEO_inst.h"
        GEO_mat_vec_mul3                "GEO_util.h"
        GEO_vec_prod3                   "GEO_util.h" 
        modsmf                          "smfio.h"
                
!Revision History:
		$Log: GEO_get_view_vec.c,v $
		Revision 4.2  2003/11/21 21:38:29  vlin
		Global variables added

		Revision 4.1  2003/08/08 16:25:42  vlin
		Moved call to GEO_get_sample_time() to a higher level.
		Changed to use generic message mnemonic.

		Revision 3.1  2002/06/11 13:13:34  kuyper
		Corrected off-by-one loop error.

		Revision 2.1  1997/10/21 18:16:22  kuyper
		Returned from ClearCase

		4/25/95
		Ruiming Chen (rchen@ltpmail.gsfc.nasa.gov)
		Finished coding.

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END
*****************************************************************************/

{
  double u_img[3] = {0.0}; /* viewing vector in instrument coordinate system */
  double n_inst_normal[3] = {0.0}; /* normalize n_inst */
  double un_prod = 0.0;  /* product of u_img and n_inst_normal */ 
  int i = 0, j = 0;
  char filefunc[] = __FILE__ ", GEO_get_view_vec";


  if (GEO_get_inst_mirr_normal(scan_number, sample_number, n_inst_normal)
	!= SUCCESS) {
    modsmf(MODIS_E_GEO, "GEO_get_inst_mirr_normal()", filefunc);
    return FAIL;
  }

  for (i = 0; i < num_detectors; i++) {
  /* transfer global variable u_tel to instrument coordinate system u_img */
    GEO_mat_vec_mul3(T_tel2inst, &u_tel[i][0], u_img);

    /* calculate the viewing vector u_inst */
    GEO_vec_prod3(u_img, n_inst_normal, &un_prod);
    for (j = 0; j < 3; ++j)
      u_inst[i][j] = u_img[j] - 2.0*un_prod*n_inst_normal[j];
  }
    return SUCCESS;
}
