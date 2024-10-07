/*
!C-INC*************************************************************************
!Description:   define functions for the geolocation software main module

!Input Parameters: N/A

!Output Parameters: N/A

!Revision History:
 * $Log: GEO_main_func.h,v $
 * Revision 6.2  2011/02/14 21:01:06  kuyper
 * Corrected const qualification of geo_params parameter of
 *   GEO_locate_one_granule().
 *
 * Revision 6.1  2010/05/26 18:41:03  kuyper
 * Dropped obsolete macro, structure type, and function parameters of that type.
 *
 * Revision 5.5  2005/03/16 21:35:44  kuyper
 * Changed header guard macro name to avoid reserved name space.
 *
 * Revision 5.4  2004/10/25 19:27:24  vlin
 * typo "GEO_in_maneuver" fixed.   vlin@saicmodis.com
 *
 * Revision 5.3  2004/10/15 21:34:40  kuyper
 * Corrected definition of maneuver_list_struct.
 *
 * Revision 5.2  2004/10/14 20:47:37  kuyper
 * Corrected declaration of spacecraft_ID argument of GEO_read_maneuver_file().
 *
 * Revision 5.1  2004/08/24 15:27:23  vlin
 * 1. maneuver_list_struct, GEO_read_maneuver_file(), and
 *    GEO_in_maneuver() added.
 * 2. GEO_locate_one_granule() updated.

!Team-unique Header:
		This software is developed by the MODIS Science Data Support
		Team for the National Aeronautics and Space Administration,
		Goddard Space Flight Center, under contract NAS5-32373.

!END***************************************************************************
*/

#ifndef GEO_MAIN_FUNC_H
#define GEO_MAIN_FUNC_H

#include <PGS_PC.h>
#include "GEO_parameters.h"
#include "smfio.h"

/* function prototypes */

PGSt_SMF_status GEO_in_maneuver(
	PGSt_double                   tai1993
	);

PGSt_SMF_status GEO_read_param_file(GEO_param_struct * const param);

PGSt_SMF_status GEO_locate_one_granule(
	char			/* l1a_file_name */ [/* PGSd_PC_FILE_PATH_MAX */],
        char			/* geo_file_name */ [/* PGSd_PC_FILE_PATH_MAX */],
        GEO_param_struct const	* /* geo_params */,
        int			const /* version */
);

/*  End of function prototypes */
 
#endif
/*  End of include file */
