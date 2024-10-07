#include "PGS_MODIS_35251.h"
#include "GEO_global_arrays.h"
#include "GEO_output.h"


#define MAX_ANGLE (PGS_PI/2.0)

int GEO_get_bounding_coords(
	 double      terrain_sample_position[MAX_DETECTORS][MAX_SCAN_SAMPLE][3],
	 unsigned char	pixel_flags[MAX_DETECTORS][MAX_SCAN_SAMPLE],
	 int    const num_detectors,
         int    const num_frames,
	 GEO_bcoord_struct * const bounding_coords )
/*
!C***************************************************************************** 
!Description:   
		Routine in Output group of the Level-1A geolocation
                software to calculate the granule bounding coordinates
		for the ECS Core Metadata.  
		
		Every pixel in the granule is inspected for bounding
		coordinates, scan by scan.  Track is kept of the most
		eastern and western coordinates in both the Eastern
		and Western Hemisphere.  These four longitudes are
		then compared to identify the granule's most eastern
		and western bounding coordinates.  If the granule is
		found to lie in both Eastern and Western hemispheres
		and near both the 0 and 180 meridians, it is identified
		as a 'polar' granule with the east bounding and west 
		bounding coordinates set to +PI and -PI, respectively,
		but the north and south bounding coordinates are still
		set from the coordinates of the northern-most and southern-
		most pixels found in the granule.

		The bounding_coords structure must be initialized to
		the following state prior to beginning the granule's
		bounding coordinate search:

		northcoord   = -PI/2
		southcoord   =  PI/2
		eastcoord    =  PI
		westcoord    = -PI
		easthemi_ebc =   0
		easthemi_wbc =  PI
		westhemi_ebc = -PI
		westhemi_wbc =   0
		
!Input Parameters:
		terrain_sample_position[MAX_DETECTORS][MAX_SCAN_SAMPLE][3]
                        - computed geodetic coordinates for pixels
			  (radians, meters)
		pixel_flags[MAX_DETECTORS][MAX_SCAN_SAMPLE] - geolocation
			pixel quality flags.
		num_detectors - the number of detectors per sample
                                        (usually 10).  Each scan has
                                        num_detectors scan lines.
		num_frames - the number of (sample) frames in the scan.


!Output Parameters:
		bounding_coords - structure of the granules bounding coordinates

Return parameter:
                SUCCESS if all bounding coordinates are identified,
		FAIL otherwise.

Global variables: None

Call functions:
                modsmf(MODIS_X_MNEMONIC_STRING, "user message string", 
		  "function, GEO_write_parameters.c");

!Revision History:
 * $Log: GEO_get_bounding_coords.c,v $
 * Revision 2.2  1997/11/06 23:18:17  kuyper
 * Changed pixel_flags to unsigned char, to match gflags.
 *
 * Revision 2.1  1997/10/21  18:16:22  kuyper
 * Returned from ClearCase
 *
 * Revision 1.5.1.1  1997/07/21  22:20:17  kuyper
 * Merged in out-of-sequence changes.
 *
 * Revision 1.5  1997/07/21  16:24:34  kuyper
 * Baselined Version 1
 *
 *Parallel Development:
 * Revision 1.7  1997/04/15  19:38:46  fhliang
 * added file name (L.1);
 * fixed prolog: '!Input Parameter' --> '!Input Parameters', and
 *               '!Output Parameter' --> '!Output Parameters'.
 *
 * Revision 1.6  1997/04/10  19:35:30  fhliang
 * Initial revision of SDST re-delivery of GEO_get_bounding_coords.c.
 *
 * Revision 1.5  1997/04/08  22:16:23  kuyper
 * Fixed program prologue.
 *
 * Revision 1.4  1997/03/26  18:05:17  fhliang
 * Initial revision of SDST delivery of GEO_get_bounding_coords.c.
 *
		Revision 1.3  1997/01/21 21:02:37  kuyper
		Merged seed files.

		Revision 1.3  1997/01/21 21:02:37  kuyper
		Merged seed files.

		Revision 1.2  1997/01/15 13:01:17  fshaw
		updated walkthrough corrections


Requirements:
                PR03-F-4.1-1
                PR03-F-4.1-2

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END***************************************************************************
*/
{

  int detector, frame;

  if ((terrain_sample_position == NULL) || (pixel_flags  == NULL)
       || (bounding_coords == NULL) || (num_detectors < 0)
       || (num_detectors > MAX_DETECTORS) || (num_frames < 0)
       || (num_frames > MAX_SCAN_SAMPLE)){
      modsmf(MODIS_E_BAD_INPUT_ARG, "",
         "GEO_get_bounding_coords, GEO_get_bounding_coords.c():");
      return FAIL;
  }

  for(detector = 0; detector < num_detectors; detector++){
     for( frame = 0; frame < num_frames; frame++){

        if ((  (pixel_flags[detector][frame] & INVALID_INPUT_DATA) == 0) && 
               ( (pixel_flags[detector][frame] & NO_ELLIPSE_INTERSECT) == 0)){ 

	    if ( terrain_sample_position[detector][frame][LAT_IDX] >
                 bounding_coords->northcoord)
	      bounding_coords->northcoord =
                 terrain_sample_position[detector][frame][LAT_IDX];

	    if( terrain_sample_position[detector][frame][LAT_IDX] <
                 bounding_coords->southcoord)
	      bounding_coords->southcoord =
                 terrain_sample_position[detector][frame][LAT_IDX];

            if( terrain_sample_position[detector][frame][LON_IDX] > 0.0 ){
                     /* ( it is in the eastern hemisphere)   */
               if( terrain_sample_position[detector][frame][LON_IDX] >
                   bounding_coords->easthemi_ebc)
		 bounding_coords->easthemi_ebc =
                    terrain_sample_position[detector][frame][LON_IDX]; 

               if( terrain_sample_position[detector][frame][LON_IDX] <
                     bounding_coords->easthemi_wbc)
		 bounding_coords->easthemi_wbc =
                     terrain_sample_position[detector][frame][LON_IDX]; 

            }else { /* (it's in the western hemisphere)*/

	        if( terrain_sample_position[detector][frame][LON_IDX] >
                      bounding_coords->westhemi_ebc)
		  bounding_coords->westhemi_ebc =
                      terrain_sample_position[detector][frame][LON_IDX]; 

	        if( terrain_sample_position[detector][frame][LON_IDX] <
                       bounding_coords->westhemi_wbc)
		  bounding_coords->westhemi_wbc =
                       terrain_sample_position[detector][frame][LON_IDX]; 
           }
        }	
      }
    }
    /* update eastern and western bounding coordinates */
    if (bounding_coords->westhemi_wbc > bounding_coords->westhemi_ebc){
     /* (the inspected portion of the granule is not in the Western hemisphere) */

      bounding_coords->eastcoord = bounding_coords->easthemi_ebc;
      bounding_coords->westcoord = bounding_coords->easthemi_wbc;

    } else if( bounding_coords->easthemi_wbc > bounding_coords->easthemi_ebc){
     /* (the inspected portion of the granule is not in the Eastern hemisphere) */

      bounding_coords->eastcoord = bounding_coords->westhemi_ebc;
      bounding_coords->westcoord = bounding_coords->westhemi_wbc;

    } else if( (bounding_coords->easthemi_wbc - bounding_coords->westhemi_ebc) <
          MAX_ANGLE){
     /* (the granule lies in both hemispheres or neither) */

      if( 2.0 * PGS_PI - (bounding_coords->easthemi_ebc -
          bounding_coords->westhemi_wbc) < MAX_ANGLE){
        /*  (the granule straddles both 0 and 180 meridians; it is 'polar')   */

          bounding_coords->eastcoord =  PGS_PI;
      	  bounding_coords->westcoord = -PGS_PI;
          bounding_coords->easthemi_ebc =  PGS_PI;
      	  bounding_coords->westhemi_wbc = -PGS_PI;
          bounding_coords->easthemi_wbc =  0.0;
      	  bounding_coords->westhemi_ebc =  0.0;

      } else{ /* (the granule straddles the 0 meridian)  */
          bounding_coords->eastcoord = bounding_coords->easthemi_ebc;
          bounding_coords->westcoord = bounding_coords->westhemi_wbc;
          bounding_coords->easthemi_wbc =  0.0;
      	  bounding_coords->westhemi_ebc =  0.0;
      }
    } else if( 2.0 * PGS_PI - (bounding_coords->easthemi_ebc - 
         bounding_coords->westhemi_wbc) < MAX_ANGLE){
      /* (the granule straddles the 180 meridian)   */

      bounding_coords->eastcoord = bounding_coords->westhemi_ebc;
      bounding_coords->westcoord = bounding_coords->easthemi_wbc;
      bounding_coords->easthemi_ebc =  PGS_PI;
      bounding_coords->westhemi_wbc = -PGS_PI;
   }
  return SUCCESS;
}
