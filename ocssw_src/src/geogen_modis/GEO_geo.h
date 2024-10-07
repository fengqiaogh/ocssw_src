/*
!C-INC*************************************************************************
!Description:   define parameters for the geolocation software

!Input Parameters: N/A

!Output Parameters: N/A

!Revision History:
	$Log: GEO_geo.h,v $
	Revision 6.7  2014/08/08 22:21:21  jkuyper
	Removed macro no longer in use.

	Revision 6.6  2011/02/09 19:06:42  kuyper
	Added external header guard for "PGS_DEM.h" to compesate for lack of
	  internal one.
	Added enumeration of values for land/sea mask.
	Added explanation of best guess at reason for odd value of landsea mask
	  fill value.

	Revision 6.5  2010/05/19 20:28:27  kuyper
	Helped resolve Bug 1969 by adding THERMCORR_FVALUE macro.

	James Kuyper Jr.	James.R.Kuyper@NASA.gov

	Revision 6.4  2009/09/15 20:15:46  kuyper
	Removed obsolete macros.

	Revision 6.3  2009/03/27 20:39:00  kuyper
	Changed macro name to MAX_PADDED.

	Revision 6.2  2009/03/18 20:53:15  kuyper
	Corrected new definition of MAX_DETECTORS.
	Added HIRES_FVALUE.

	Revision 6.1  2009/02/20 22:22:09  kuyper
	Added new macros needed to support position averaging and high-resolution
	offsets.

	Revision 5.4  2008/12/15 15:41:18  kuyper
	Added LUN to control orbit number validation.

	Revision 5.3  2005/03/16 21:35:11  kuyper
	Changed header guard macro name to avoid reserved name space.

	Revision 5.2  2004/11/05 01:46:54  kuyper
	Corrected definition of L_SMASK_FVALUE.

	Revision 5.1  2004/08/24 18:54:45  vlin
	L_SMASK_FVALUE changed to match the fill value used in DEM file.

	Revision 4.1  2002/11/19 23:11:55  kuyper
	Added TEMP_FVALUE

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END***************************************************************************
*/

#ifndef GEO_GEO_H
#define GEO_GEO_H

#include "GEO_basic.h"
/* PGS_DEM.h does not have internal header guards, must provide externally. */
#ifndef H_PGS_DEM
    #define H_PGS_DEM
    #include "PGS_DEM.h"
#endif

/*********************************************************
#defines
*********************************************************/

/* General geolocation definitions */
#define MAX_SCAN_NUMBER		1000
#define MAX_FRAMES		1354
#define SAMPLE_RATIO_HKM	2
#define SAMPLE_RATIO_QKM	4
#define SAMPLES_HKM		(MAX_FRAMES*SAMPLE_RATIO_HKM)
#define SAMPLES_QKM		(MAX_FRAMES*SAMPLE_RATIO_QKM)
#define MAX_PADDED		(SAMPLES_QKM+3)
#define DETECTORS_1KM		10
#define DETECTORS_HKM		(DETECTORS_1KM*SAMPLE_RATIO_HKM)
#define DETECTORS_QKM		(DETECTORS_1KM*SAMPLE_RATIO_QKM)
#define MAX_SCAN_SAMPLE		MAX_FRAMES	/* Obsolescent */
#define MAX_DETECTORS		DETECTORS_1KM	/* Obsolescent */

#define MAX_IMPULSE_NUMBER 25
#define MAX_POLY_DEGREE 4
#define EA_SOURCE_SELECT_LUN 600280
#define VALIDATE_ORBIT_NO_LUN 600281
#define TERRAIN_CORRECT_LUN 600310

/* Instrument model */
#define MAX_BAND_NUMBER 36
#define FIRST_BAND 30
#define LATCH_TO_CENTER 0.5
#define BASE_SAMPLES	DETECTORS_1KM		/* Obsolescent */
#define MAX_SAMPLES	SAMPLE_RATIO_QKM	/* Obsolescent */
#define ROLL 0
#define PITCH 1
#define YAW 2
#define POSITION 0
#define VELOCITY 3
#define ELEC_SIDES 2

/* Input L1A data definitions */
#define ENCODER_LENGTH 78
#define SECTOR_LENGTH 40
#define ANCIL_LENGTH 64
#define NUM_L1A_QUALITY_FLAGS 4
#define NUM_L1A_SECTOR_VIEWS 6

/* Flag values */
#define INVALID_INPUT_DATA 128
#define NO_ELLIPSE_INTERSECT 64
#define BAD_TERRAIN 32
#define NEAR_LIMB 4
#define INVALID_SENSOR_ANGLES 8
#define INVALID_SOLAR_ANGLES 8
#define TAI_FLAG 1.e10
#define TERRAIN_CORRECT "TRUE"

/* Land/Water mask values. */
enum { SHALLOW_OCEAN, DRYLAND, COAST, SHALLOW_INLAND, EPHEMERAL, DEEP_INLAND,
    CONTINENTAL, DEEP_OCEAN, NUM_LWMASK,
    L_SMASK_FVALUE = (0x10000 + PGSd_DEM_NO_FILLVALUE)>>8};
/* That peculiar expression for L_SMASK_FVALUE has a value of 221. It is the
 * actual value we find for the _FillValue attribute in the SDP Toolkit's DEM
 * files. The expression given above reflects my best guess as to where that
 * peculiar value came from.
 * PGSd_DEM_NO_FILLVALUE (-8888) is the value that the SDP Toolkit DEM routines
 * use for the Land/Water mask, apparently indicating that it was not supposed
 * to have fill values. If someone called the appropriate HDF function to give
 * it a fill value anyway, and used a 16 bit 2's complement big-endian integer
 * with a value of -8888 to do so, that would explain why the actual value came
 * out as 221. That's because an HDF _FillValue attribute is always of the same
 * type as the SDS, and this SDS is stored as uint8, so that HDF routine would
 * have looked only at the first byte of that integer.
 */

/* Fill values. */
#define LONG_FVALUE (-999.0)
#define LAT_FVALUE (-999.0)
#define HGHT_FVALUE (-32767)
#define SENSORZEN_FVALUE (-32767)
#define SENSORAZIM_FVALUE (-32767)
#define RANGE_FVALUE 0
#define SOLARZEN_FVALUE (-32767)
#define SOLARAZIM_FVALUE (-32767)
#define TEMP_FVALUE (-999.0)
#define THERMCORR_FVALUE (-999.0)
#define GFLAGS_FVALUE 255
#define HIRES_FVALUE (-128)

/* DEM definitions */
#define MAX_DEM_ROWS 108
#define MAX_DEM_TILES 26410
#define MAX_DEM_HORIZONTAL 200
#define MAX_DEM_VERTICAL 200
#define NO_DEM_DATA (-1)

/* DEM validation definitions */
/* Nominal values with about 10% margin, in m */
#define MIN_TERRAIN_HEIGHT (-450.0)
#define MAX_TERRAIN_HEIGHT 9600.0

/*	Mathematical constants (that should have been in ANSI math.h)	*/
#define PGS_PI 3.14159265358979323846
#define RAD2DEG (180.0/PGS_PI)
#define DEG2RAD (PGS_PI/180.0)
#define MAX_UINT16_VAL (0xFFFF)

/*  End of definitions */

/* Structures	*/
/* Structure for bounding coordinate ECS core metadata (radians)	*/
typedef struct {
  double  northcoord;	/* Northern-most granule bounding latitude	*/
  double  eastcoord;	/* Eastern-most granule bounding longitude	*/
  double  southcoord;	/* Southern-most granule bounding latitude	*/
  double  westcoord;	/* Western-most granule bounding longitude	*/
  double  easthemi_ebc;	/* Eastern-most granule bounding longitude
			    in the Eastern hemisphere.	*/
  double  easthemi_wbc; /* Western-most granule bounding longitude
			    in the Eastern hemisphere.	*/
  double  westhemi_ebc;	/* Eastern-most granule bounding longitude
			    in the Western hemisphere.	*/
  double  westhemi_wbc; /* Western-most granule bounding longitude
			    in the Western hemisphere.	*/
} GEO_bcoord_struct;

#endif

