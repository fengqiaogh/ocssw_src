#ifndef GEO_DEM_H
#define GEO_DEM_H
/*!C-INC************************************************************************

!Description:   
	Shared information about the set of DEM file resolutions to use.

!Input Parameters:
	N/A

!Output Parameters:
	N/A

Return Values:
	N/A

Externally Defined:
	N/A

Called by:
	N/A

Routines Called:
	N/A

!Revision History:
$Log: GEO_DEM.h,v $
Revision 6.2  2011/02/08 18:29:25  kuyper
Added external header guards for PGS_DEM.h, to compensate for lack of
  internal ones.

Revision 6.1  2010/08/12 22:46:17  kuyper
Initial Revision.

James Kuyper Jr.		James.R.Kuyper@NASA.gov

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits

Design Notes
!END**************************************************************************
*/
/* PGS_DEM.h does not have internal header guards, must provide externally. */
#ifndef H_PGS_DEM
    #define H_PGS_DEM
    #include "PGS_DEM.h"
#endif
#define RESOLUTIONS 1
extern PGSt_DEM_Tag DEM_resolutions[RESOLUTIONS];
extern double min_lat, max_lon;

#endif
