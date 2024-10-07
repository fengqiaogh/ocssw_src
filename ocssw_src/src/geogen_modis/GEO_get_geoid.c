/* DEBUG */
#include <math.h>

#include "GEO_earth.h"
#include "GEO_geo.h"
#include "PGS_MODIS_35251.h"
#include "smfio.h"

int GEO_get_geoid
(
	double latitude,
	double longitude
)
/*
!C*****************************************************************************
!Description:   
		Routine in the DEM group of the Level-1A geolocation
                software to retrieve and store the geoid height information.

!Input Parameters:
                latitude	latitude of the requested location (radians)
                longitude	longitude of the requested location (radians)

!Output Parameters:
                None

Return value:
		The height of the geoid in meters, or BAD_GEOID if it could not
		be retrieved, bilinearly interpolated to the specified latitude
		and longitude.

Global variables:
		None

Called by:
		GEO_terrain_correct()

Call functions:
		modsmf()
		PGS_DEM_GetSize()
		PGS_DEM_GetQualityData()

!Revision History:
# $Log: GEO_get_geoid.c,v $
# Revision 6.1  2011/02/15 18:44:31  kuyper
# Changed to bilinearly interpolated geoid heights, to avoid artifacts at the
#   arbitrary boundaries used by older nearest-neighbor algorithm.
#
# Revision 4.1  2003/10/27 01:14:30  vlin
# Message buffer size expanded.
#
# Revision 2.3  1999/03/12 17:37:49  fhliang
# Capitalized Prolog section names.
#
 * Revision 2.2  1999/02/02  18:04:58  kuyper
 * Make some implicit conversions explicit.
 *
 * Revision 2.1  1999/02/02  17:14:59  lma
 * *** empty log message ***
 *
#
# Liqun Ma <lma@ltpmail.gsfc.nasa.gov>

!Requirements:
		PR03-F-3.3.1-1

!Team-unique Header:
		This software is developed by the MODIS Science Data Support
		Team for the National Aeronautics and Space Administration,
		Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/

{
    #define ROWS 180
    #define COLS 360

    /* Arbitrary limits imposed by poor design of SDP Toolkit, for calls to
     * PGS_DEM_GetQualityData().
     */
    #define MIN_LAT (-90.0*(1.0-1.0/(double)ROWS))
    #define MAX_LON (180.0*(1.0-1.0/(double)COLS))

    typedef int16 geoid_height_t[COLS];
    static geoid_height_t *geoid_height = NULL;
    if(!geoid_height) {
        geoid_height = calloc(ROWS, sizeof(geoid_height_t));
        geoid_height[0][0] = BAD_GEOID;
    }
    int lat_lo, lat_hi, lon_lo, lon_hi;
    double dlat, dlon;
    char msgbuf[128] = "";
    static char filefunc[] = __FILE__ ", GEO_get_geoid"; 

    latitude  *= RAD2DEG;
    longitude *= RAD2DEG;

    /* Accept arbitrary longitudes, adjust into range. */
    /* Algorithm is inefficient if fabs(longitude) is large; it's only intended
     * to cover possibilities such as use of the range 0<=longitude<360 */
    while( longitude < -180.0 )
	longitude += 360.0;
    while( longitude >= 180.0 )
	longitude -= 360.0;

    lon_lo = (int)floor(longitude+179.5);
    dlon = longitude + 179.5 - lon_lo;
    lon_hi = lon_lo + 1;

    /* wrap around at 180E/180W line. */
    if(lon_lo < 0)
	lon_lo = COLS-1;
    else if(lon_hi > COLS-1)
	lon_hi = 0;
    
    lat_lo = (int)floor(89.5-latitude);
    dlat = 89.5 - latitude - lat_lo;
    lat_hi = lat_lo + 1;

    /* Flatten out at poles. */
    if(lat_lo < 0)
	dlat = lat_lo = lat_hi = 0;
    else if(lat_hi > ROWS-1)
    {
	lat_lo = lat_hi = ROWS-1;
	dlon = 0.0;
    }

    if( geoid_height[0][0] == BAD_GEOID )
    {	/* First call to this routine. */
	PGSt_integer numPixVertical;
	PGSt_integer numPixHorizontal;
	PGSt_integer sizeDataType;
	PGSt_double  longitude_arry[2] = { -180.0,MAX_LON};
	PGSt_double  latitude_arry[2] = {90.0, MIN_LAT};

	if( PGS_DEM_GetSize(PGSd_DEM_30ARC, PGSd_DEM_GEOID, PGSd_DEM_DEGREE,
	    latitude_arry, longitude_arry, &numPixVertical, &numPixHorizontal,
	    &sizeDataType ) != PGS_S_SUCCESS )
	{
	    modsmf(MODIS_E_GEO, "PGS_DEM_GetSize", filefunc);

	    return BAD_GEOID;
	}


	if( (numPixVertical != ROWS) || (numPixHorizontal != COLS) ||
	    ((size_t)sizeDataType != sizeof(geoid_height[0][0])) )
	{
	   sprintf(msgbuf,
	       "\nnumPixVertical= %ld numPixHorizontal= %ld sizeDataType= %ld", 
	       (long)numPixVertical,(long)numPixHorizontal,(long)sizeDataType);
	   modsmf(MODIS_E_GEO_BAD_SIZE, msgbuf, filefunc);

	   return BAD_GEOID;
	}


	if( PGS_DEM_GetQualityData(PGSd_DEM_30ARC, PGSd_DEM_GEOID,
	    PGSd_DEM_DEGREE, latitude_arry, longitude_arry,
	    (void *)geoid_height) != PGS_S_SUCCESS )
	{
	   modsmf(MODIS_E_GEO, "PGS_DEM_GetQualityData", filefunc);
	   return BAD_GEOID;
	}
    }

    return (int)floor( (1-dlon)*((1.0-dlat)*geoid_height[lat_lo][lon_lo] +
	dlat*geoid_height[lat_hi][lon_lo]) +
	dlon*((1.0-dlat)*geoid_height[lat_lo][lon_hi] +
	dlat*geoid_height[lat_hi][lon_hi]) + 0.5); 
}

