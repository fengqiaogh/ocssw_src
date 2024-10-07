/*
 * Contains: GEO_close_DEM(), GEO_DEMalloc(), GEO_initialize_DEM(),
 *	GEO_latlon2height(), and GEO_read_DEM().
 *
 * Revision History:
 * $Log: GEO_DEM.c,v $
 * Revision 6.4  2011/12/08 19:58:03  kuyper
 * Corrected calculation of minimum latitude considered valid by the SDP
 *   Toolkit's DEM routines.
 *
 * Revision 6.3  2011/02/15 18:44:17  kuyper
 * Changed to call GEO_get_geoid() whenever needed, rather than calling it only
 * once and storing the value. This needs to be done, now that the geiod height
 * is being bilinearly interpolated, and will therefore not be exactly
 * constant.
 *
 * Revision 6.2  2010/12/14 22:01:54  kuyper
 * Changed to return heights above ellipsoid, rather than geoid.
 *
 * Revision 6.1  2010/08/12 22:09:49  kuyper
 * Changed to use 15 arc second DEM files, falling back to 30 arc second if
 *   necessary.
 * Improved error messaging.
 * Changed to use GEO_DEM.h
 * Changed GEO_initialize_DEM() to return status codes.
 *
 * Revision 4.1  2003/08/11 15:40:07  vlin
 * GEO_read_DEM():  Removed conversion of fill values, specified initialization.
 *
 * Revision 3.3  2002/06/17 20:07:20  kuyper
 * Added messages to accompany all FAIL returns from GEO_read_DEM().
 *
 * Revision 3.2  2002/06/13 22:46:42  kuyper
 * Removed unneeded NCSA acknowledgement
 *
 * Revision 3.1  2002/06/11 13:16:28  kuyper
 * Increased the border size.
 *
 * Revision 2.0  1997/07/28  16:29:41  kuyper
 * Merged together GEO_DEM.h, GEO_initialize_DEM.c, GEO_read_DEM.c, and
 *   GEO_latlon2height.c
 * Added GEO_close_DEM().
 * Total rewrite of GEO_intialize_DEM() and GEO_read_DEM() to make use of
 *   new Toolkit DEM utilities.
 */

#include <float.h>
#include <ctype.h>
#include <stddef.h>
#include "PGS_MODIS_35251.h"
#include "GEO_DEM.h"
#include "GEO_earth.h"
#include "smfio.h"

/* Values from Toolkit DEM metadata used to identify spatial units. */
#define DEGREE_STRING "Decimal Degree (DD)"
#define METER_STRING "METERS"
#define DEMCACHE	50     /* number of entries for tile_cache_struct */

/* Toolkit DEM interface parameters. */
/* PGSt_DEM_Tag DEM_resolutions[RESOLUTIONS] = {PGSd_DEM_15ARC, PGSd_DEM_30ARC}; */
PGSt_DEM_Tag DEM_resolutions[RESOLUTIONS] = {PGSd_DEM_30ARC};
double min_lat = LAT_FVALUE;
double max_lon = LONG_FVALUE;

/* We may increase the number of layers, with time. */
#define LAYERS 1
static PGSt_integer DEM_layers[LAYERS]={PGSd_DEM_ELEV};
static PGSt_double fillValue[LAYERS]={-9999.0};
static PGSt_double offset[LAYERS]={0.0};
static PGSt_double scaling[LAYERS]={1.0};

/* DEM tile parameters. */
static int num_rows;	/* Number of tile rows. */
static int hor_size; 	/* Pixels per pixel row. */
static int ver_size;	/* Number of pixel rows per tile. */
static int hor_start;	/* Pixel column corresponding to corner_lon. */
static int ver_start;	/* Pixel row corresponding to corner_lat. */
static int num_tiles;	/* Total number of tiles in all rows. */
static double tile_ver_size=8;	/* Latitude spacing of tiles (radians) */
/* Initialization of tile_ver_size and row_num_tiles avoids division by zero in
 * GEO_read_DEM(), in case it is accidentally called without calling
 * GEO_initialize_DEM().
 */

/* DEM file row indicies */
/* Number of first tile in each tile row. */
static short int row_start_tile[MAX_DEM_ROWS];
/* Number of tiles in each tile row. */
static short int row_num_tiles[MAX_DEM_ROWS]={1};

/* DEM tile record parameters */      
static short int *DEM_record; /* Indicates which index to use on
		      height_min and height_max for each tile. Obsolescent,
		      but used by GEO_latlon2height() as an error indicator. */
static short int *height_min; /* Minimum height within a tile. */
static short int *height_max; /* Maximum height within a tile. */

/* Current tile parameters */
static short int *DEM_tile_data=NULL;	/* Elevation (meters). */
static int current_tile=-1;		/* Number of current tile. */
static double corner_lat; /* Latitude (radians) corresponding to ver_start */
static double corner_lon; /* Longitude (radians) corresponding to hor_start */
static double lat_inc=DEG2RAD/120.0; /* Radians of latitude between pixels */
static double lon_inc=DEG2RAD/120.0; /* Radians of Longitude between pixels */


typedef struct {
	int	tile;
	int	hor_size;
	int	ver_size;
	int	hor_start;
	int	ver_start;
	double  corner_lat;
	double	corner_lon;
	short int *DEM_tile_data;
}tile_cache_struct;

tile_cache_struct	*tile_cache;
unsigned short int	tiles_cached=0;


void geoDemInitGlobals() {
    static int firstRun = 1;
    if(firstRun) {
        firstRun = 0;
        DEM_record = (short int*) calloc(MAX_DEM_TILES, sizeof(short int));
        height_min = (short int*) calloc(MAX_DEM_TILES, sizeof(short int));
        height_max = (short int*) calloc(MAX_DEM_TILES, sizeof(short int));
        tile_cache = (tile_cache_struct*) calloc(DEMCACHE, sizeof(tile_cache_struct));
    }
}

PGSt_SMF_status GEO_initialize_DEM(
	void
)
/*
!C*****************************************************************************
!Description:   
	Routine in earth location group of the Level-1A geolocation software to
	open the DEM file and read the file parameters into global variables.

!Input Parameters:
	None

!Output Parameters:
	None

Return parameter:
	MODIS_E_GEO		If any subroutine returns an error code
	MODIS_E_DEM_METADATA	If a DEM metadata field has an unexpected value
	PGS_S_SUCCESS		otherwise

Externally defined:
    File scope objects defined in GEO_DEM.c:
	Inputs:
	    DEM_layers
	    DEM_resolutions	"GEO_DEM.h"
		
	Outputs:
	    fillValue
	    height_max
	    height_min
	    lat_inc
	    lon_inc
	    num_rows
	    num_tiles
	    max_lon		"GEO_DEM.h"
	    min_lat		"GEO_DEM.h"
	    offset
	    row_num_tiles
	    row_start_tile
	    scaling
	    tile_ver_size

    Macros:
	DEG2RAD		"GEO_geo.h"
	EQUATOR_TILES	"GEO_DEM.c"
	LAYERS		"GEO_DEM.c"
	MAX_DEM_ROWS	"GEO_geo.h"
	MAX_DEM_TILES	"GEO_geo.h"
	NO_DEM_DATA	"GEO_geo.h
	PGS_S_SUCCESS	"PGS_SMF.h"
	PGSd_DEM_30ARC	"PGS_DEM.h"
	PGSd_DEM_ELEV	"PGS_DEM.h"
	RESOLUTIONS	"GEO_DEM.h"

Functions called:
	PGS_DEM_Open() - Opens Toolkit DEM files.
	PGS_DEM_GetMetadata() - Gets metadata from DEM files.
	modsmf() - writes status messages to log

!Revision History:
 See top of file for revision history.


!Requirements:
                PR03-F-3.3.3-1


!Team-unique Header:
		This software is developed by the MODIS Science Data Support
		Team for the National Aeronautics and Space Administration,
		Goddard Space Flight Center, under contract NAS5-32373.


!END*************************************************************************
*/
/* Adjustable parameters chosen to minimize average I/O load. */
#define TILES_BASE 90.0
#define COSLAT_OFFSET 0.1059

{
    geoDemInitGlobals();
    
  /* Local variable declarations */
  /* The first element of each pix*Info pair is the size of the pixels in
    degrees. The second element is the offset of the pixel center (the point at
    which the retrieved data is most relevant) from the NW corner of the pixel.
   */
  PGSt_double	pixLatInfo[2]={1.0/120.0, 1.0/240.0};
  PGSt_double	pixLonInfo[2]={1.0/120.0, 1.0/240.0};
  double	top_lat; /* Latitude of tile edge farthest from equator */
  double	coslat;	/* cos(top_lat) */
  int	row, tile;	/* Loop indices. 	*/
  char	positionUnits[30]="", dataUnits[30]="";
  char	msgbuf[sizeof(positionUnits)+sizeof(dataUnits)+64];
  char	*p1, *p2; 	/* First non-blank characters of Units strings*/
  char	filefunc[] = __FILE__ ", GEO_initialize_DEM";

  /* Begin program logic */
  if(PGS_DEM_Open(DEM_resolutions, RESOLUTIONS, DEM_layers, LAYERS)
    != PGS_S_SUCCESS)
  {
    sprintf(msgbuf, "PGS_DEM_Open(%d,%d)", RESOLUTIONS, LAYERS);
    modsmf(MODIS_E_GEO, msgbuf, filefunc);

    return MODIS_E_GEO;
  }

 
  if(PGS_DEM_GetMetadata(DEM_resolutions[0], DEM_layers[0], pixLatInfo,
      pixLonInfo, positionUnits, &scaling[0], &offset[0], &fillValue[0],
      dataUnits, NULL, NULL) !=PGS_S_SUCCESS)
  {
    sprintf(msgbuf, "PGS_DEM_GetMetadata(%ld,%ld)",
      (long)DEM_resolutions[0], (long)DEM_layers[0]);
    modsmf(MODIS_E_GEO, msgbuf, filefunc);

    return FAIL;
  }
  
  for(p1=positionUnits; p1 < positionUnits+sizeof(positionUnits) &&
    isspace((int)*p1); p1++);
  for(p2=dataUnits; p2 < dataUnits+sizeof(dataUnits) &&
    isspace((int)*p2); p2++);

  if(strncmp(p1, DEGREE_STRING,
      sizeof(positionUnits)-(size_t)(p1-positionUnits)) ||
    strncmp(p2, METER_STRING, sizeof(dataUnits)-(size_t)(p2-dataUnits)) ||
    fabs(scaling[0]-1.0) > DBL_EPSILON || fabs(offset[0]) > DBL_MIN )
  {
    sprintf(msgbuf, "\"%.*s\", \"%.*s\" scale:%g offset:%g",
   	 (int)sizeof(positionUnits), positionUnits,
	 (int)sizeof(dataUnits), dataUnits, scaling[0], offset[0]);
    modsmf(MODIS_E_DEM_METADATA, msgbuf, "GEO_DEM.c, GEO_initialize_DEM");

    return FAIL;
  }

  /* Setting these variables here allows GEO_latlon2height() to remain
    unchanged from version 1. */
  lon_inc = pixLonInfo[0]*DEG2RAD;
  lat_inc = -pixLatInfo[0]*DEG2RAD;
  min_lat = pixLatInfo[0]*0.5 - 90.0;
  max_lon = 180.0 - 0.5 * pixLonInfo[0];
  num_rows = MAX_DEM_ROWS;
  tile_ver_size = PGS_PI / (double)num_rows;

  row_start_tile[0] = 0;
  for(row=0; row<num_rows; row++)
  {	/* Calculate number of tiles in row. */

    /* Calculate latitude of the tile edge farthest from equator. */
    top_lat = -0.5*PGS_PI+(double)row*tile_ver_size;
    if(top_lat>0.0)
      top_lat += tile_ver_size;
    
    coslat = cos(top_lat);
    /* Calculate the optimum number of tiles in this row. This formula is an
      empirical fit to exact values calculated using numerical methods and
      timing data from an SGI Power Challenger, running version 5.2 of the
      Toolkit. For details, check "Derivation of Optimum Tile Layout" in the
      Software Development Folder.*/
    row_num_tiles[row] =
      (short)floor(TILES_BASE*coslat/sqrt(coslat+COSLAT_OFFSET)+0.5);
    /* GEO_latlon2height() will fail if there is not at least one tile in each
      tile row, and may fail if there is more than one tile in a polar tile
      row. */
    if(row_num_tiles[row] < 1 || row==0 || row==num_rows-1)
      row_num_tiles[row] = 1;
    
    if(row==num_rows-1)
      num_tiles = row_start_tile[row]+row_num_tiles[row];
    else
      row_start_tile[row+1] = (short)(row_start_tile[row]+row_num_tiles[row]);
  }

  for(tile=0; tile<num_tiles; tile++)
  {	/* Set values indicating no data. */
    height_min[tile] = SHRT_MAX;
    height_max[tile] = SHRT_MIN;
    DEM_record[tile] = NO_DEM_DATA;
  }

  return SUCCESS;
}

/*============================================================================*/

void	*GEO_DEMalloc(
	   size_t	size
	)

/*
!C*****************************************************************************
!Description:   
	 	Rourine in the DEM group of the Level-1A geolocation
		software to allocate memory for new tiles, freeing older tiles
		if necessary.

!Input Parameters:
		size		Number of bytes to allocate

!Output Parameters:
		None

Return parameter:
		NULL			If the requested memory could not be
					allocated.
		A pointer to the	Otherwise.
		allocated memory.

Global variables:
		unsigned		tiles_cached
		tile_cache_struct	tile_cache

Call functions:
		malloc()
		free()
Called by:
                GEO_read_DEM

!Revision History:
    See top of file for revision history.

!Requirements:
		PR03-F-3.3.3-1

!Team-unique Header:
		This software is developed by the MODIS Science Data Support
		Team for the National Aeronautics and Space Administration,
		Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/
{
    geoDemInitGlobals();

 	void	*p;

	if( tiles_cached > DEMCACHE )
	  return NULL;	/* This should never happen.	*/

	p = malloc(size);

	while( p == NULL && tiles_cached > 0 )
	 {
	
	  tiles_cached--;

	  free(tile_cache[tiles_cached].DEM_tile_data);
	  tile_cache[tiles_cached].DEM_tile_data = (short *)NULL;

	  p = malloc(size);

	 }

  return p;

}

/*============================================================================*/

int GEO_read_DEM(
	double lat,
	double lon,
	int * const hgtmin,
	int * const hgtmax
)
/******************************************************************************
!C
!Description:   
	Routine in Earth location group of the Level-1A geolocation software to
	determine the DEM tile number for the input lat/lon, load the min, max
	and flag and read the DEM record into a global variable.

!Input Parameters:
	lat		the input latitude (radians)
	lon		the input longitude (radians)

!Output Parameters:
	hgtmin		the minimum DEM height above the ellipsoid for the tile
			(meters)
	hgtmax		the maximum DEM height above the ellipsoid for the tile
			(meters)

Return Values:
	FAIL		if NULL pointers are passed, or Toolkit DEM functions
			    or GEO_DEMalloc() fails, or if the Toolkit DEM
			    dataSize has changed.
	SUCCESS		otherwise

Externally Defined:
    File globals defined in GEO_initialize_DEM:
	Input:
	    DEM_resolutions	"GEO_DEM.h"
	    fillValue
	    lat_inc
	    lon_inc
	    num_rows
	    offset
	    row_start_tile
	    row_num_tiles
	    scaling
	    tile_ver_size
	Input/Output:
	    DEM_flags
	    DEM_record
	    height_max
	    height_min
    File globals defined here:
	Input/Output:
	    current_tile	index of current tile in memory,
				    initialize to -1
	    DEM_tile_data	DEM data for current tile, initialize to NULL
	Output:
	    corner_lat		latitude of first element of tile (radians)
	    corner_lon		longitude of first element of tile (radians)
	    hor_size		Horizontal dimension of tile array, including
				    any overlap with adjacent tiles (pixels)
	    ver_size		Vertical dimension of tile array, including any
				    overlap with adjacent tiles (pixels)
	    hor_start		Horizontal array index aligned with tile corner
	    ver_start		Vertical array index aligned with tile corner
	    tile_cache		details about cached tiles.
	    tiles_cached	Count of valid entries in tile_cache

    Macros:
	BAD_GEOID			"GEO_earth.h"
	DEG2RAD				"GEO_geo.h"
	DEMCACHE			"GEO_DEM.c"
	FAIL				"GEO_basic.h"
	LAT_FVALUE			"GEO_geo.h"
	LAYERS				"GEO_DEM.c"
	LONG_FVALUE			"GEO_DEM.c"
	MAX_DEM_ROWS			"GEO_geo.h"
	MAX_DEM_TILES 			"GEO_geo.h"
	max_lon				"GEO_DEM.h"
	min_lat				"GEO_DEM.h"
	MODIS_E_BAD_DEM			"PGS_MODIS_35251.h"
	MODIS_E_DEM_DATA_SIZE		"PGS_MODIS_35251.h
	MODIS_E_GEO			"PGS_MODIS_35251.h"
	PGS_FALSE			"PGS_SMF.h"
	PGS_S_SUCCESS			"PGS_SMF.h"
	PGSd_DEM_DEGREE			"PGS_DEM.h"
	PGSd_DEM_ELEV			"PGS_DEM.h"
	PGSd_DEM_NEAREST_NEIGHBOR	"PGS_DEM.h"
	RAD2DEG				"GEO_geo.h"
	RESOLUTIONS			"GEO_DEM.h"
	SUCCESS				"GEO_basic.h"

Called by:
	GEO_terrain_correct
		
Routines called:
	GEO_get_geoid()			"GEO_earth.h"
	GEO_DEMalloc()			"GEO_earth.h"
	PGS_DEM_GetSize()		"PGS_DEM.h"
	PGS_DEM_GetRegion()		"PGS_DEM.h"
	PGS_SMF_TestErrorLevel()	"PGS_SMF.h"
 	modsmf()			"smfio.h"

Requirements:
	PR03-F-3.3.3-1
	PR03-F-3.3.3-2	(Not implemented yet - waived)

!Revision History:
      See top of file for revision history. 

		6/26/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Finished coding.
	
!Team-unique Header:
		This software is developed by the MODIS Science Data Support
		Team for the National Aeronautics and Space Administration,
		Goddard Space Flight Center, under contract NAS5-32373.

Limitations: The scale and offset coefficients provided for scaling the
		DEM data are not applied to them.

!END
*****************************************************************************/

/* The number of radians of latitude corresponding to 26km along the Earth's
   surface. Oblateness is too small to matter, but if it did, this should be
   the value for locations near the poles
*/
#define BORDER	0.004084

{
  geoDemInitGlobals();
  /* Parameters describing forbidden zone for DEM files with version 5.2 of
  Toolkit. See the README file in the DEM file directory. Hopefully later
  versions will correct this. */

  /* Labels for the parts of tiles which cross longitude +/-180 */
  enum {WEST, EAST, PARTS};

  double		tile_hor_size=1.0;    /* Hor. dimension of tile */
  /* tile_ver_size is constant, calculated by GEO_initialize_DEM(). */
  PGSt_double		latitude[2], longitude[PARTS][2];
  PGSt_double		firstElement[2]={0.0};
  double		border=0.0;
                        /* Radians of longitude corresponding to 25km */
  long			pole = 0L;
  PGSt_integer		sizeDataType=0;
  PGSt_integer		numPixVertical=0, numPixHorizontal[PARTS];
  int16			*tempdata[PARTS]={NULL,NULL};
  int			row, col;	/* Row and column tile indices.*/
  int			i, j;		/* Pixel indices within tiles.*/
  int			tile;		/* tile index. */
  int16			*inrow, *outrow;
  char			msgbuf[64];
  unsigned short	int ci;
  tile_cache_struct	temp_cache;
  char filefunc[] = __FILE__ ", GEO_read_DEM";


  if( hgtmin == NULL || hgtmax == NULL ||
    (fabs(lat) > PGS_PI*0.5*(1.0+DBL_EPSILON)) )
  {
    sprintf(msgbuf, "lat:%g hgtmin:%p hgtmax:%p", lat, (void *)hgtmin, 
           (void *)hgtmax);
    modsmf(MODIS_E_BAD_DEM, msgbuf, filefunc);

    return FAIL;
  }

  /* Compute tile index for input lat/lon */
  /* Compute row number */
  /* Force lat into range -PI/2 <= lat <= PI/2 */
  if( lat >= PGS_PI*0.5 )
      lat = PGS_PI*0.5-DBL_EPSILON*10.0;
  else if( lat < -PGS_PI*0.5 )
	   lat = -PGS_PI*0.5;

  row = (int)floor((lat + PGS_PI*0.5)/tile_ver_size);

  /* Compute column number */
  tile_hor_size = PGS_PI*2.0/(double)row_num_tiles[row];

  /* Force lon to be in the range -PI <= lon < PI */
  lon -= 2.0*PGS_PI*floor(lon/(PGS_PI*2.0)+0.5);
  col = (int)floor((lon + PGS_PI)/tile_hor_size);

  /* Compute tile number */
  tile = row_start_tile[row] + col;

  /* search for tile_cache[ci].tile == tile */
  for( ci = 0; ci < tiles_cached && tile_cache[ci].tile!=tile; ci++ );

  if(ci<tiles_cached)
  {	/* Matching entry was found.	*/
    temp_cache = tile_cache[ci];

    if( ci > 0 )
    {
      hor_size		= temp_cache.hor_size;
      ver_size		= temp_cache.ver_size;
      hor_start		= temp_cache.hor_start;
      ver_start		= temp_cache.ver_start;
      corner_lat	= temp_cache.corner_lat;
      corner_lon	= temp_cache.corner_lon;
      DEM_tile_data	= temp_cache.DEM_tile_data;
    }
  }
  else
  {	/* Need to load in a new tile of data */
    if( ci == DEMCACHE )
    { /* Drop least recently used cache entry, to make room for new one.*/
      free(tile_cache[DEMCACHE-1].DEM_tile_data);

      /* Initialize pointer if calls to Free or malloc are made later */
      tiles_cached--;
      ci = tiles_cached;
    }

    current_tile = tile;

    /* set output to error values, in case of early returns or errors*/
    *hgtmin	= SHRT_MAX;
    *hgtmax	= SHRT_MIN;
    hor_size = ver_size = 0;
    /* Guarantees early exit from GEO_latlon2height(). */

    DEM_tile_data = (short int *)NULL;


    /*From here on in returning with DEM_tile_data==NULL indicates an error.*/
    /* Calculate latitude boundaries. Nominal tiles are defined by the
     * calculations above. Full tile boundaries must include every point
     * within a great-circle distance of BORDER radians from any point in the
     * nominal tile area.
     */

    if( row == 0 )  
    {	/* south polar tile. */
      latitude[0] = -90.0 + (tile_ver_size+BORDER)*RAD2DEG;
      latitude[1] = min_lat;
    }
    else if( row == num_rows-1 )
    {	/* North polar tile. */
      latitude[0] = 90.0;
      latitude[1] = 90.0 - (tile_ver_size+BORDER)*RAD2DEG;
    }
    else
    {	/* Non-polar tile. */
      corner_lat  = (double)row*tile_ver_size - PGS_PI/2.0;
      latitude[0] = (corner_lat+tile_ver_size+BORDER)*RAD2DEG; 
      latitude[1] = (corner_lat-BORDER)*RAD2DEG;

      /* Use the edge farthest from the equator to calculate the border size. */
      if( corner_lat < 0.0 )
	border = BORDER/cos(corner_lat);
      else
	border = BORDER/cos(corner_lat+tile_ver_size);
    }

    /* Calculate the longitude boundaries. For tiles whose border cross +/-180
     * degrees longitude, two sets of logitude boundaries are needed, one for
     * each side of that barrier. If only one set is needed, it is stored in
     * the WEST part.
     * Note that the WEST and EAST parts of the tile are named from the point
     * of view of the tile, not from the prime meridian; the WEST part is in
     * the eastern hemisphere, and just east of it is the EAST part, which is
     * in the western hemisphere.
     */

    if( row_num_tiles[row] == 1 )
    {	/*Tile covers entire range of longitude, doesn't need a border.*/
      longitude[WEST][0] = -180.0;
      longitude[WEST][1] = max_lon;
    }
    else
    {	/* Calculate East/West border sizes. */
      longitude[WEST][0] = ((double)col*tile_hor_size-border)*RAD2DEG-180.0;
      longitude[WEST][1] = longitude[WEST][0] + 
			    (tile_hor_size+2.0*border)*RAD2DEG;
      if( longitude[WEST][0] < -180.0 )
      {	/* Western border overlaps -180.0, must be split */
	longitude[EAST][0] = -180.0;
	longitude[EAST][1] = longitude[WEST][1];
	longitude[WEST][0] += 360.0;
	longitude[WEST][1] = max_lon;
      }
      else if(longitude[WEST][1] > 180.0 )
      {	/* Eastern border overlaps 180.0, must be split */
	longitude[EAST][0] = -180.0;
	longitude[EAST][1] = longitude[WEST][1]-360.0;
	longitude[WEST][1] = max_lon;
      }
    }

    /* Retrieve data for first part of tile */
    if( PGS_DEM_GetSize(DEM_resolutions[0], PGSd_DEM_ELEV, PGSd_DEM_DEGREE,
      latitude, longitude[WEST], &numPixVertical, &numPixHorizontal[WEST],
      &sizeDataType)
      != PGS_S_SUCCESS || numPixVertical<=0 || numPixHorizontal[WEST]<=0)
    {
	sprintf(msgbuf, "PGS_DEM_GetSize({%g,%g}, {%g,%g})", latitude[0],
	    latitude[1], longitude[WEST][0], longitude[WEST][1]);
	modsmf(MODIS_E_GEO, msgbuf, "GEO_DEM.c, GEO_read_DEM");

	return FAIL;
    }

    if( (size_t)sizeDataType != sizeof(tempdata[0][0]) )
    {	/* Unexpected data size.        */
      sprintf(msgbuf, "%ld", (long)sizeDataType);
      modsmf(MODIS_E_DEM_DATA_SIZE, msgbuf, "GEO_DEM.c, GEO_read_DEM");
      return FAIL;
    }

    ver_size  = numPixVertical;
    hor_size  = numPixHorizontal[WEST];

    tempdata[WEST] = (int16 *)GEO_DEMalloc(
      (size_t)(numPixVertical*numPixHorizontal[WEST]*sizeDataType) );

    if( tempdata[WEST] != NULL )
    {	/* Get data for WEST part. */
      if( PGS_SMF_TestErrorLevel(PGS_DEM_GetRegion(DEM_resolutions, RESOLUTIONS,
	PGSd_DEM_ELEV, PGSd_DEM_DEGREE, PGSd_DEM_NEAREST_NEIGHBOR,
	latitude, longitude[WEST], tempdata[WEST], NULL, firstElement,
	NULL))==PGS_FALSE )
      {	/* Got elevation data for WEST part */
	corner_lat = firstElement[0]*DEG2RAD;
	corner_lon = firstElement[1]*DEG2RAD;
	hor_start  = ver_start = 0;

	if( row_num_tiles[row]==1 )
	{
	  /* Need extra columns on East and West edges, to allow
	   * interpolation
	   */
	  hor_size += 2;

	  /* Change column which corresponds to firstElement[1] */
	  hor_start++;

	  /* Polar tile rows must never contain more than 1 tile, or
	   * GEO_latlon2height() could fail. Therefore, we only need
	   * to test the following cases if row_num_tiles[row]==1
	   */
	  if( row == 0 )
	  {
	    /* Need an extra row on South edge, to allow interpolation
	     *  all the way to South pole
	     */
	    ver_size++;

	    pole = 0L;
	    inrow = tempdata[WEST] + (numPixVertical-2)*numPixHorizontal[WEST];
	    for( i = 0; i < numPixHorizontal[WEST]; i++ )
	      pole += (long)inrow[i];
	    pole /= (long)numPixHorizontal[WEST];
	  }
	  else if( row == num_rows-1 )
	  {
	    /* Need an extra row at the North edge, to allow
	     * interpolation all the way to north pole
	     */
	    ver_size++;
	    /* Change row which corresponds to firstElement[0] */
	    ver_start++;

	    pole = 0L;
	    inrow = tempdata[WEST];
	    for( i = 0; i < numPixHorizontal[WEST]; i++ )
	    pole += (long)inrow[i];
	    pole /= (long)numPixHorizontal[WEST];
	  }

	  DEM_tile_data = (int16 *)GEO_DEMalloc((size_t)
	    (ver_size*hor_size*sizeDataType));

	  if( DEM_tile_data != NULL )
	  {	/* Copy tempdata[WEST] to DEM_tile_data. */
	    for( i = 0; i < numPixVertical; i++ )
	    {
	      inrow = tempdata[WEST]+i*numPixHorizontal[WEST];
	      
	      /* Note that copy is offset by ver_start, for North polar row. */
	      outrow = DEM_tile_data + (i+ver_start)*hor_size;

	      /* Fill in extra buffer columns needed for interpolation. */
	      outrow[0] = inrow[numPixHorizontal[WEST]-1];
	      outrow[hor_size-1] = inrow[0];

	      for( j = 0; j < numPixHorizontal[WEST]; j++ )
		outrow[j+1] = inrow[j];
	    }

	    if( row == 0 )
	      outrow = DEM_tile_data + (ver_size-1)*hor_size;
	    else if(row == num_rows-1)
	      outrow = DEM_tile_data;
	    else
	      outrow = NULL;

	    if( outrow != NULL )	/* Fill in average near-polar value. */
	      for( j = 0; j < hor_size; j++ )
		outrow[j] = (int16)pole;
	  }
	}
	else if( col == 0 || col == row_num_tiles[row]-1 )
	{
	  /* This tile straddles longitude +/-180 degrees. */
	  if( col == 0 )
	    corner_lon -= 2.0*PGS_PI;

	  if( PGS_DEM_GetSize(DEM_resolutions[0], PGSd_DEM_ELEV,
	    PGSd_DEM_DEGREE, latitude, longitude[EAST], &numPixVertical,
	    &numPixHorizontal[EAST], NULL)
	    == PGS_S_SUCCESS && numPixVertical>0 && numPixHorizontal[EAST]>0)
	  {
	    tempdata[EAST] = (int16 *)GEO_DEMalloc((size_t)
	      (numPixVertical*numPixHorizontal[EAST]*sizeDataType));
	  }

	  if( tempdata[EAST] != NULL )
	  {	/* Get elevation for EAST part. */
	    if( PGS_SMF_TestErrorLevel(PGS_DEM_GetRegion( DEM_resolutions,
		RESOLUTIONS, PGSd_DEM_ELEV, PGSd_DEM_DEGREE,
		PGSd_DEM_NEAREST_NEIGHBOR, latitude, longitude[EAST],
		tempdata[EAST], NULL, NULL, NULL) ) == PGS_FALSE )
	    {
	      /* Got elevation for EAST part. */
	      hor_size += numPixHorizontal[EAST];
	      DEM_tile_data = (int16 *)GEO_DEMalloc((size_t)
		(ver_size*hor_size*sizeDataType));
	    }
          }

   	  if( DEM_tile_data != NULL )
	  {	/* Merge data from the two parts. */
	    for( i = 0; i < ver_size; i++ )
	    {
	      /* copy tempdata[WEST] into first numPixHorizontal[WEST]
	       * columns of DEM_tile_data.
	       */
	      outrow = DEM_tile_data + i*hor_size;
	      inrow  = tempdata[WEST] + i*numPixHorizontal[WEST];
	      for( j = 0; j < numPixHorizontal[WEST]; j++ )
	        outrow[j] = inrow[j];

	      /* copy tempdata[EAST] into last numPixHorizontal[EAST]
	       * columns of DEM_tile_data.
	       */
	      outrow += numPixHorizontal[WEST];
	      inrow   = tempdata[EAST] + i*numPixHorizontal[EAST];
	      for( j = 0; j < numPixHorizontal[EAST]; j++ )
	        outrow[j] = inrow[j];
	    }
          }
	}
	else
	{	/* No special case for this tile. */
	  DEM_tile_data  = tempdata[WEST];
	  tempdata[WEST] = NULL;
	}

	if( DEM_tile_data != NULL )
	{
	  int ver, hor;

	  for(i = ver = 0; ver < ver_size; ver++)
	      for(hor = 0; hor < hor_size; hor++,i++)
	  {	
      /* In principle, the following lines should be inserted here:
       *
       *   SET DEM_tile_data[i] = DEM_tile_data[i]*scaling[0]+offset[0]
       *  (rounded to nearest integer)
       *
       * However, this code assumes scaling[0]==1 and offset[0]==0, which
       * is currently the case. Those assumptions are checked in
       * GEO_initialize_DEM(), in case they change in the future.
       */ 
	    if( DEM_record[tile] == NO_DEM_DATA )
	    {	/* First time this tile has been loaded. */
		/* Cumulate height_min, height_max.     */
	      double lat = (ver - ver_start)*lat_inc + corner_lat;
	      double lon = (hor - hor_start)*lon_inc + corner_lon;
	      int geoid_height = GEO_get_geoid(lat, lon);

	      if(geoid_height == BAD_GEOID)
	      {
		  sprintf(msgbuf, "GEO_get_geoid(%g,%g) = %d",
		      lat, lon, geoid_height);
		  modsmf(MODIS_E_GEO, msgbuf, filefunc);
		  return FAIL;
	      }

	      DEM_tile_data[i] += geoid_height;
	      if( DEM_tile_data[i] < height_min[tile] )
		height_min[tile] = DEM_tile_data[i];

	      if( DEM_tile_data[i] > height_max[tile] )
		height_max[tile] = DEM_tile_data[i];
	    }
	  }
	  /* DEM_record is set solely to avoid changing GEO_latlon2height();
	   * the actual value doesn't matter, as long as it's not NO_DEM_DATA.
	   */
	  DEM_record[tile] = (short)tile;

	  temp_cache.tile    		= tile;
	  temp_cache.hor_size		= hor_size;
	  temp_cache.ver_size		= ver_size;
	  temp_cache.hor_start		= hor_start;
	  temp_cache.ver_start		= ver_start;
	  temp_cache.corner_lat		= corner_lat;
	  temp_cache.corner_lon		= corner_lon;
	  temp_cache.DEM_tile_data	= DEM_tile_data;

	  tiles_cached++;

	  if( ci == 0 )
	    tile_cache[0] = temp_cache;

	}	/* END if DEM_tile_data has been successfully allocated */
      }	/* END if first call to PGS_DEM_Get_Region() succeeded  */

      free(tempdata[WEST]);
      free(tempdata[EAST]);

    }	/* end if tempdata[WEST] has been successfully allocated */
  }	/* end if DEM tile not in memory */

  if( DEM_tile_data == NULL )
  {
    modsmf(MODIS_E_GEO, "GEO_DEMalloc()", "GEO_DEM.c, GEO_read_DEM");

    return FAIL;
  }

  /* Can return valid values. */
  if( ci > 0 )
  {	/* overlapping copy; must use memmove(), not memcpy().	*/
    memmove(&tile_cache[1], tile_cache, (size_t)(ci)*sizeof(tile_cache[0]));
    tile_cache[0] = temp_cache;
  }

  *hgtmin   = (int)height_min[tile];
  *hgtmax   = (int)height_max[tile];

  return SUCCESS;
}

/*===========================================================================*/

int GEO_latlon2height(
	double const lat,
	double const lon,
	double * const h
	)

/*
!C******************************************************************************
!Description:   
		Routine in earth location group of the Level-1A geolocation
		software to determine terrain height given lat/lon.  It 
		uses an array of terrain heights for a tile previously read
		into memory by GEO_DEM_read.  The offset into the tile is
		computed for the lat/lon and the height is computed by 
		bilinear interpolation from the four adjacent points

!Input Parameters:
		double lat - latitude 
		double lon - longitude

!Output Parameters:
		double *h - height 

Return parameter:
		int err_stat - error status

Global variables:
		int hor_size - Horizontal dimension of tile array (includes any 
		  overlap with adjacent tiles)
		int ver_size - Vertical dimension of tile array (includes any 
		  overlap with adjacent tiles)
		int hor_start - Horizontal array index aligned with tile corner
		int ver_start - Vertical array index aligned with tile corner
		short int DEM_record[MAX_DEM_TILES] - Pointers to DEM records
		  (-1 = no DEM record for tile) 
		short int height_min[MAX_DEM_TILES] - Minimum terrain height 
		  for tile (meters)
		short int height_max[MAX_DEM_TILES] - Maximum terrain height 
		  for tile (meters)
		short int DEM_TILE_DATA[MAX_DEM_HORIZONTAL*MAX_DEM_VERTICAL] -
		  DEM data for current tile
		int current_tile - index of current tile in memory
		double corner_lat - latitude of tile lower-left corner
		double corner_lon - longitude of tile lower-left corner
		double lat_inc - latitude increment of data in tile
		double lon_inc - longitude increment of data in tile


Call functions:
		modsmf(MODIS_X_MNEMONIC_STRING, "user message string", 
		  "function", GEO_latlon2height.c") - writes error 
		  status messages to log

!Revision History:
   see top of file for revision history. 

		6/27/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Finished coding.
		 
!Team-unique Header:
		This software is developed by the MODIS Science Data Support
		Team for the National Aeronautics and Space Administration,
		Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/
{
  /* Local variable definitions */

  double hor = 0.0; /* relative horizontal location in tile */
  double ver = 0.0; /* relative vertical location in tile */
  double dihor = 0.0; /* integral part of horizontal position */ 
  double diver = 0.0; /* integral part of vertical position */ 
  double dhor = 0.0; /* fractional part of horizonal position */
  double dver = 0.0; /* fractional part of vertical position */
  double dlon = 0.0; /* difference in longitude */
  double dlat = 0.0; /* difference in latitude */
  int ihor = 0; /* horizontal index into tile */
  int iver = 0; /* vertical index into tile */
  int near_DEM[2][2] = {0}; /* DEM points surrounding lat/lon */  

  /* Begin program logic */


  /* Compute offset relative to current tile */

  dlon = lon - corner_lon;
  dlat = lat - corner_lat;

  /* Check for longitude rollover */

  if(dlon > 2.0*PGS_PI) 
    dlon = dlon - 2.0*PGS_PI;

  if(dlon < 0.0) 
    dlon = dlon + 2.0*PGS_PI;

  /* Compute indices into tile data */
  
  hor = dlon / lon_inc + (double)hor_start;
  ver = dlat / lat_inc + (double)ver_start;

  /* Check if we have required data */

  if( (hor < 0.0) || (hor >= (double)hor_size) ||
 		(ver < 0.0) || (ver >= (double)ver_size)) {
    /* call SDP function to report error */
    modsmf(MODIS_E_DEM_IN_MEM, "",
    "GEO_DEM.c, GEO_latlon2height");

    return FAIL;
  }

  /* Check flag to see if there is data for this tile */
  
  if(DEM_record[current_tile] == NO_DEM_DATA) {
    *h = 0.0;
    
  }
  else {

    /* Compute integral and fractional parts of index */

    dhor = modf(hor, &dihor);
    dver = modf(ver, &diver);
    ihor = (int)dihor;
    iver = (int)diver;

    /* Get four nearest DEM points */

    near_DEM[0][0] = (int)DEM_tile_data[iver * hor_size + ihor];
    near_DEM[0][1] = (int)DEM_tile_data[iver * hor_size + ihor + 1];
    near_DEM[1][0] = (int)DEM_tile_data[(iver + 1) * hor_size + ihor];
    near_DEM[1][1] = (int)DEM_tile_data[(iver + 1) * hor_size + ihor + 1];

    /* Compute height as weighted sum */

    *h = (double)near_DEM[0][0] * (1.0 - dhor) * (1.0 - dver)
       + (double)near_DEM[0][1] * dhor * (1.0 - dver)
       + (double)near_DEM[1][0] * (1.0 - dhor) * dver
       + (double)near_DEM[1][1] * dhor * dver;

  } /* End if there is data for this tile */
  return SUCCESS;
}

/*===========================================================================*/

int GEO_close_DEM(void)

/*
!C*****************************************************************************
!Description:   
		Routine in the DEM group of the Level-1A geolocation
		software to close access to the DEM files and free space
		allocated to store it

!Input Parameters:
		None

!Output Parameters:
		None

Return value:
		FAIL if PGS_DEM_Close() call fails
		SUCCESS otherwise


Externally defined:
                DEM_resolution
                DEM_layers
                tile_cache	 details about cached tiles.
                tiles_cached	Count of valid entries in tile_cache

Called by:
                main()


Call functions:
                PGS_DEM_Close()		"PGS_DEM.h"

!Revision History:
 * Revision /main/GEO_V2_DEV/1 1997/08/26 kuyper
 * Corrected SUCCEED to SUCCESS.

!Requirements:
		PR03-F-3.3.3-1

!Team-unique Header:
		This software is developed by the MODIS Science Data Support
		Team for the National Aeronautics and Space Administration,
		Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/
{
    char msgbuf[64];
    char filefunc[] = __FILE__ ", GEO_close_DEM";

    if( tiles_cached > DEMCACHE )
    {
	/* Something is seriously wrong. Safest approach is to not clean up. */
	sprintf(msgbuf, "tiles_cached = %u", tiles_cached);
	modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
	return MODIS_E_BAD_INPUT_ARG;
    }

    while( tiles_cached > 0 )
	free(tile_cache[--tiles_cached].DEM_tile_data);

    if(PGS_DEM_Close(DEM_resolutions, RESOLUTIONS, DEM_layers, 1)
	!= PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "PGS_DEM_Close(%d, 1)", RESOLUTIONS);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);

	return MODIS_E_GEO;
    }

    return PGS_S_SUCCESS;
}

