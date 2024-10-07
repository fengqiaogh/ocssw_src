#include "GEO_DEM.h"
#include "GEO_output.h"
#include "GEO_earth.h"
#include "PGS_MODIS_35251.h"

PGSt_SMF_status GEO_landsea_mask(
	int	num_samples,
	int	num_detectors,
        double	terrain_sample_position[][MAX_PADDED][3],
        uint8	sample_flags[][MAX_PADDED],
	uint8	*land_seamask_qaflag,
	uint8	sample_landsea[][MAX_PADDED]
)

/*
!C***************************************************************************
!Description:   
        subroutine in output group of the Level-1A geolocation software to
        retrieve the EOS format land/sea mask data for the scan from the SDP
        Toolkit.

!Input Parameters:
	num_samples		The number of samples to be processed.
        terrain_sample_position	High resolution ground locations in geodetic
                                coordinates.
        sample_flags		High resolution pixel flags.

!Output Parameters:
	landseamask_qaflag	QA flag for land/water mask issues.
	sample_landsea		High resolution Land/Water mask values
 
Return Values:
        MODIS_E_BAD_INPUT_ARG   If any argument is invalid
        MODIS_E_GEO             If any subroutine fails
        PGS_S_SUCCESS           Otherwise
 
Externally Defined:
	DEM_resolutions		"GEO_DEM.h"
        DETECTORS_QKM           "GEO_geo.h"
        L_SMASK_FVALUE          "GEO_geo.h"
        LAT_IDX			"GEO_earth.h"
        LON_IDX		        "GEO_earth.h"
	max_lon			"GEO_DEM.h"
        MAX_PADDED              "GEO_geo.h"
	min_lat			"GEO_DEM.h"
        MODIS_E_BAD_INPUT_ARG   "PGS_MODIS_35251.h"
        MODIS_E_GEO             "PGS_MODIS_35251.h"
        NO_ELLIPSE_INTERSECT    "GEO_geo.h"
        PGS_FALSE               "PGS_SMF.h"
        PGS_S_SUCCESS           "PGS_SMF.h"
 
Called by:
        GEO_locate_one_scan	"GEO_earth.h"
 
Routines Called:
        PGS_DEM_GetPoint        "PGS_DEM.h"
        PGS_SMF_TestErrorLevel  "PGS_SMF.h" 

!Revision History:
   $Log: GEO_landsea_mask.c,v $
   Revision 6.6  2011/02/09 19:42:14  kuyper
   Changed to work at high resolution, rather than 1km resolution.
   Change to pass landsea information directly, rather than through scan_data.

   Revision 6.5  2010/12/16 22:08:37  kuyper
   Switched from using *_index to *_IDX.

   Revision 6.4  2010/08/03 20:18:18  kuyper
   Changed to use 15 arc second DEM files, falling back to 30 arc second if
     necessary.
   Changed to use GEO_DEM.h

   Revision 6.3  2009/09/11 18:28:54  kuyper
   Corrected location of calculation of the pixel index.

   Revision 6.2  2009/05/31 01:30:11  ltan
   Minor corrections.

   Revision 6.1  2009/05/29 15:03:25  ltan
   Moved inputs to proper section of prolog.
   Changed MAX_SCAN_SAMPLE to MAX_FRAMES, MAX_DETECTORS to DETECTORS_1KM.
   Changed num_samples to num_frames.
   Changed to use terrain_frame_position and frame_flags, rather than the corresponding members of scan_data.
   Changed to return status codes.

   Revision 5.1  2004/11/03 17:53:01  vlin
   unsigned integer conversion for scan_data->land_seamask.
   vlin@saicmodis.com

   Revision 4.2  2003/08/15 14:23:53  vlin
   Changed to process by detector within a frame, rather than by frame within
   a detector, in hopes of improving tile caching by PGS_DEM_GetPoint().

   Revision 4.1  2003/02/21 22:07:48  kuyper
   Corrected to use void* pointers with %p format code.
  
Requirements:	PR03-F-3.8-1

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END***************************************************************************
*/

{
  PGSt_double	*latitude = (PGSt_double*) malloc(DETECTORS_QKM*MAX_PADDED*sizeof(PGSt_double));
  PGSt_double	*longitude = (PGSt_double*) malloc(DETECTORS_QKM*MAX_PADDED*sizeof(PGSt_double));
  PGSt_integer	numPoints=0;	/* The number of point to retrieve.*/
  int		samp, det;
  int8  	*land_seamask = (int8*) malloc(DETECTORS_QKM*MAX_PADDED*sizeof(int8));
  char  	msgbuf[256];
  char          filefunc[] = __FILE__ ", GEO_landsea_mask";

  if (terrain_sample_position == NULL || sample_landsea == NULL || 
      sample_flags == NULL)
  {
     sprintf(msgbuf, "terrain_sample_position:%p sample_landsea:%p\n"
	"sample_flags:%p\n", (void*)terrain_sample_position,
	(void*)sample_landsea, (void*)sample_flags);
     modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);

     return MODIS_E_BAD_INPUT_ARG;
  }
  else if (num_samples <0 || num_samples>MAX_PADDED ||
      num_detectors <=0 || num_detectors>DETECTORS_QKM)
  {
      sprintf(msgbuf, "samples:%d detectors:%d", num_samples, num_detectors);
      modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);

      return MODIS_E_BAD_INPUT_ARG;
  }

  numPoints = 0;

/* The LandWater SDS in the DEM files is chunked with a chunk size of 600x800,
 * corresponding to 20 degrees by 26 degrees 40 arcmin. Therefore, since MODIS
 * scans are long and skinny, processing by detector within a sample, rather
 * than by sample within a detector, improves chunk caching by
 * PGS_DEM_GetPoint().
 */
  for (samp = 0; samp < num_samples; samp++)
  {
    for (det = 0; det < num_detectors; det++)
    {
      if (sample_flags[det][samp] < NO_ELLIPSE_INTERSECT)
      {	/* Convert lat/lon to degrees and pack earth-viewing locations
         for input to PGS_DEM_GetPoint: */
         latitude[numPoints] =
	     terrain_sample_position[det][samp][LAT_IDX] * RAD2DEG;
         longitude[numPoints] =
	     terrain_sample_position[det][samp][LON_IDX] * RAD2DEG;

         /* Adjust to avoid region that the SDP Toolkit handles incorrectly*/
         if (latitude[numPoints] < min_lat)
	     latitude[numPoints] = min_lat;
         if (longitude[numPoints] > max_lon)
	     longitude[numPoints] = max_lon;
         numPoints++;
      }
    }
  }

  if (numPoints == 0)
  {
      *land_seamask_qaflag = 1;     /* indicating no valid data */

      /* Normal processing may include 0-sample scans, or scans with bad data*/
      return PGS_S_SUCCESS;
  }

  if (PGS_SMF_TestErrorLevel( PGS_DEM_GetPoint(DEM_resolutions, RESOLUTIONS,
      PGSd_DEM_WATER_LAND, PGSd_DEM_DEGREE, latitude, longitude, numPoints,
      PGSd_DEM_NEAREST_NEIGHBOR, land_seamask)) != PGS_FALSE)
  {
      *land_seamask_qaflag = 1;     /* indicating invalid data. */
      return MODIS_E_GEO;
  }
  else
  {	/* indicate (at least some) valid land/sea mask data */
      *land_seamask_qaflag = 0;
      numPoints = 0;
      for (samp = 0; samp < num_samples; samp++)
        for (det = 0; det < num_detectors; det++)
	{	/* unpack the land_seamask data to match the other scan data
		 * for output.
		 */
          if (sample_flags[det][samp] >= NO_ELLIPSE_INTERSECT){
	      sample_landsea[det][samp] = (uint8)L_SMASK_FVALUE;
          }
          else {
	      if (land_seamask[numPoints] >= 0)	
		 sample_landsea[det][samp] = (uint8)land_seamask[numPoints];
	      numPoints++;
          }
        }
  }

  free(latitude);
  free(longitude);
  free(land_seamask);
          
  return PGS_S_SUCCESS;
}
