#include "PGS_MODIS_35251.h"
#include "GEO_global_arrays.h"  /* Just for LAT_IDX, LON_IDX, HT_IDX */
#include "GEO_validation.h"
#include "mapi.h"
#include "smfio.h"
 
PGSt_SMF_status GEO_validate_earth_location (
	int 		     const mirr_impulse_flag,
        double 	             terrain_sample_position[][MAX_PADDED][3],
	int                  const sample_number,
	int 		     const num_detectors,
	unsigned char        sample_flags[][MAX_PADDED])

/*
!C******************************************************************************
!Description:   
		subprogram in validation group of the Level-1A geolocation
                software to validate the pixel earth locations.  The
		pixel_flags are set for any earth location where validity 
		checks fail.  

!Input Parameters:
		mirr_impulse_flag - mirror impulse encoder flag for this scan
		terrain_sample_position
	        	- computed geodetic coordinates for pixels
		sample_number - the sample number within the scan
		num_detectors - number of detectors in the sample

!Output Parameters:
		sample_flags - geolocation flags for pixels

Return Values:
        	MODIS_E_BAD_INPUT_ARG   If any input argument fails validity tests.
	        PGS_S_SUCCESS           Otherwise

Externally Defined:
	        BAD_TERRAIN             "GEO_geo.h"
   	    	DETECTORS_QKM           "GEO_geo.h"
 	        HT_IDX                  "GEO_global_arrays.h"
	        INVALID_INPUT_DATA      "GEO_geo.h"
	        LAT_IDX                 "GEO_global_arrays.h"
	        LON_IDX                 "GEO_global_arrays.h"
	        MAX_TERRAIN_HEIGHT      "GEO_geo.h"
	        MIN_TERRAIN_HEIGHT      "GEO_geo.h"
  	        MODIS_E_BAD_INPUT_ARG   "PGS_MODIS_35251.h"
	        NO_ELLIPSE_INTERSECT    "GEO_geo.h"
 	        MAX_PADDED              "GEO_geo.h"
 	        PGS_S_SUCCESS           "PGS_SMF.h"

Called by:
     		GEO_locate_one_scan     "GEO_output.h"

Routines Called:
     		modsmf                  "smfio.h"

!Revision History:
 * $Log: GEO_validate_earth_location.c,v $
 * Revision 6.3  2010/12/15 23:26:52  kuyper
 * Removed setting of terrain_sample_position fill values which was redundant
 *   in C5, and pointless in C6.
 *
 * Revision 6.2  2009/05/29 14:22:05  xgeng
 * Removed 3 locally defined constants and use global constants instead.
 *
 * Revision 6.1  2009/05/19 22:09:22  xgeng
 * Replaced MAX_SCAN_SAMPLE with PADDED_SAMPLES
 * Replaced MAX_DETECTORS with DETECTORS_QKM.
 * Renamed frame_number to sample_number, and pixel_flags to sample_flags.
 * Changed to return a status code.
 * Changed MAX_SCAN_SAMPLE to MAX_PADDED.
 * The above changes match with PDL revision 6.1 and 6.2.
 *
 *
 * Xu Geng (xu.geng@saic.com)
 *	
 * Revision 4.2  2003/10/24 21:16:25  vlin
 * expanded the size of message buffer.
 *
 * Revision 4.1  2003/02/21 22:51:30  kuyper
 * Corrected to use void* pointers with %p format code.

		5/22/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Finished coding (stub version).

Requirements:
		PR03-F-3.2.3-2
		PR03-F-4.1-1

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

Design Notes:
                The design specifies the same fill values as
                GEO_initialize_product for terrain_sample_position latitude,
                longitude, and height. The structures that are used for that
                initialization are not availible to this routine. If
                GEO_initialize_product is changed then LAT_FILL, LON_FILL, and
                HT_FILL will have to be changed to correspond.

!END*************************************************************************
*/

{
  /* Local variable declarations */
  int det; /* loop index */
  char msgbuf[128];
  char   filefunc[] = __FILE__ ", GEO_validate_earth_location";

  /* Begin program logic */
  
  if((sample_flags == NULL) || (terrain_sample_position == NULL) || 
    (num_detectors < 0) || (num_detectors > DETECTORS_QKM) || 
    (sample_number < 0) || (sample_number >= MAX_PADDED))
  {
      sprintf(msgbuf, "sample_flags: %p,terrain_sample_position: %p, "
	  "num_detectors: %d, sample_number = %d", (void*)sample_flags,
	  (void*)terrain_sample_position, num_detectors, sample_number);
      modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
      return MODIS_E_BAD_INPUT_ARG;
  } 
  
  /* Check for invalid mirror data */
  if (mirr_impulse_flag == BAD_DATA) {
    for (det = 0; det < num_detectors; det++) { 
      sample_flags[det][sample_number] |= INVALID_INPUT_DATA;
    }
  }
   
  for (det = 0; det < num_detectors; det++) {  
    if( (sample_flags[det][sample_number] &
	(INVALID_INPUT_DATA | NO_ELLIPSE_INTERSECT)) == 0 &&
	( (terrain_sample_position[det][sample_number][HT_IDX] <
	  MIN_TERRAIN_HEIGHT) ||
	  (terrain_sample_position[det][sample_number][HT_IDX] >
	  MAX_TERRAIN_HEIGHT) ) )
    {
      sample_flags[det][sample_number] |= BAD_TERRAIN;
    }
  }
  
  return PGS_S_SUCCESS;
  
} /* END GEO_validate_earth_location */
