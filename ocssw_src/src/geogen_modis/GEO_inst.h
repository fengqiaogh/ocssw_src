/*
!C-INC*************************************************************************
!Description:   the prototype of the modules which implement the instrument 
                model for the Level-1A geolocation software

!Input Parameters: N/A

!Output Parameters: N/A

!Revision History:
 * $Log: GEO_inst.h,v $
 * Revision 6.1  2006/10/10 18:41:10  kuyper
 * Changed GEO_get_sample_time to take a focal_plane_geometry parameter, and to
 *   return a double.
 *
 * Revision 5.1  2005/03/16 21:35:28  kuyper
 * Changed header guard macro name to avoid reserved name space.
 *
 * Revision 4.1  2002/12/05 19:25:09  kuyper
 * Removed obsolete constant.
 *
 * Revision 2.2  1999/03/12 17:48:37  kuyper
 * Capitalized Prolog Sections
 *
 * Revision 2.1  1997/10/21  18:15:47  kuyper
 * Returned from ClearCase
 *
 * Revision /main/GEO_V2_DEV/2 1997/08/14 kuyper
 * Added ORDER_CHEBY
 *
 * Revision 1.4.1.1  1997/07/18  22:40:58  kuyper
 * Merged in out-of-sequence changes.
 *
 * Revision 1.4  1997/07/18  21:58:00  kuyper
 * Baselined Version 1
 *
 *Parallel development:
 * Revision 1.5  1997/04/22  18:26:52  fhliang
 * commented un used argument (LL.42-44).
 *
 * Revision 1.4  1997/03/26  19:10:17  fhliang
 * Initial revision of SDST delivery of GEO_inst.h.
 *
		Revision 1.3  1996/07/24 22:16:04  kuyper
		Inserted required '!'s in comments.
		Declared arguments const.

		Revision 1.2  1996/07/18 22:16:58  kuyper
		Included GEO_geo.h for MAX_DETECTORS definition.


		4/27/95
		Ruiming Chen (rchen@ltpmail.gsfc.nasa.gov)
		Finished coding

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END***************************************************************************
*/

#ifndef GEO_INST_H
#define GEO_INST_H
#include "GEO_parameters.h"

/*********************************************************
function prototypes
*********************************************************/

double GEO_get_sample_time(
	focal_plane_geometry_struct const * geometry_params,
	const int sample_number
	);

int GEO_interp_mirr_enc(
	int const scan_number,
	int const sample_number,
	double * const sample_enc
	);

int GEO_interp_mirr_ang(
	double const sample_enc,
	double * const sample_mirr_ang
	);

int GEO_get_inst_mirr_normal(
	int const scan_number,
	int const sample_number,
	double n_inst_normal[3]
	);

int GEO_get_view_vec(
	int const scan_number,
	int const sample_number,
	int const num_detectors,
	double u_inst[][3]
	);

#endif

