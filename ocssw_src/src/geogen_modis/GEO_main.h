/* file: GEO_main.h */

/*
!C-INC*************************************************************************
!Description:   define parameters for the geolocation software main module

!Input Parameters: N/A

!Output Parameters: N/A

!Revision History:
		$Log: GEO_main.h,v $
		Revision 5.1  2005/03/16 21:35:36  kuyper
		Changed header guard macro name to avoid reserved name space.

		Revision 2.2  1999/03/12 17:48:37  kuyper
		Capitalized Prolog Sections

 * Revision 2.1  1997/10/21  18:15:47  kuyper
 * Returned from ClearCase
 *
 * Revision /main/GEO_V2_DEV/1 29-Sept 1997  
 * modifying to correct defects found in W/T 9/25. 
 *
 * Revision 1.6.1.1  1997/07/18  22:40:58  kuyper
 * Merged in out-of-sequence changes.
 *
 * Revision 1.6  1997/07/18  21:58:00  kuyper
 * Baselined Version 1
 *
 * Revision 1.6  1997/06/20  20:16:27  fhliang
 * changed '#define l1aID 500000' --> '#define l1aID           500100';
 * changed '#define l1a_MCFID 600110' --> '#define l1a_MCFID       500500'.
 *
 * Revision 1.5  1997/03/26  19:10:35  fhliang
 * Initial revision of SDST delivery of GEO_main.h.
 *
		Revision 1.4  1996/12/23 18:41:02  kuyper
		Added geo_MCFID.

		Revision 1.3  1996/10/10 20:05:56  kuyper
		Corrected, rearranged, and added to LUNs.

		Revision 1.2  1996/07/24 22:16:35  kuyper
		Inserted required '!'s in comments.

                6/30/95
                Tracey W. Holmes (holmes@modis-xl.gsfc.nasa.gov)
                Finished coding.

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END***************************************************************************
*/

#ifndef GEO_MAIN_H
#define GEO_MAIN_H

#include "PGS_MODIS_35251.h"

/*********************************************************
#defines
*********************************************************/
#define l1aID		500100
#define geoID		600000
#define PARAM		600020
#define l1a_MCFID	500500
#define geo_MCFID	600111

/*  End of definitions */

#endif

/*  End of include file */
