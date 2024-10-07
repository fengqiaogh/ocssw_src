#include "PGS_MODIS_35251.h"
#include "GEO_global_arrays.h"
#include "GEO_product.h"
#include "GEO_output.h"
#include "GEO_geo.h"

PGSt_SMF_status GEO_write_scan_data(
	MODFILE			* const geo_file,
	int                     const scan_number,
        swath_elem_struct const	* swath_elem,
        GEO_param_struct const	* GEO_param,
        double                  terrain_frame_position[][MAX_FRAMES][3],
        double                  frame_to_sensor[][MAX_FRAMES][3],
        double                  frame_solar_angles[][MAX_FRAMES][2],
        uint8                   frame_flags[][MAX_FRAMES],
	uint8			frame_landsea[][MAX_FRAMES],
	uint8			frame_waterpresent[][MAX_FRAMES],
        int8                    hires_offsets[][DETECTORS_QKM][SAMPLES_QKM]
	)
/*
!C******************************************************************************
!Description:   
		subroutine in output group of the Level-1A geolocation
                software to write one scan of spatial element geolocation
		data to the geolocation product file.  

!Input Parameters:
       geo_file                MAPI structure for geolocation product to be written to
       scan_number             The scan number to be written.
       swath_elem              The swath values to be written.
       GEO_param               Contains the angle and range scale values.
       terrain_frame_position  Geodetic coordinates of ground position
       frame_to_sensor         Polar coordinates of spacecraft from ground.
       frame_solar_angles      Solar azimuth and zenith angles from ground.
       frame_flags             Pixel level flags
       hires_offsets           High resolution offsets from bilinear
                                    interpolation of low resolution positions.
       frame_landsea	       Land/Water mask.
       frame_waterpresent      Count of high-resolution pixels with water.


!Output Parameters:    None

Return Values:
       MODIS_E_BAD_INPUT_ARG   If any argument is invalid.
       MODIS_E_GEO             If all scan data could not be written.
       PGS_S_SUCCESS           Otherwise
 
Externally Defined:
       azimuth_index           "GEO_global_arrays.h"
       DETECTORS_QKM           "GEO_geo.h"
       DETECTORS_1KM           "GEO_geo.h"
       GFLAGS                  "GEO_product.h"
       HEIGHT_OFFSET           "GEO_product.h"
       HT_IDX		       "GEO_global_arrays.h"
       INVALID_INPUT_DATA      "GEO_geo.h"
       LAND_SEAMASK            "GEO_product.h"
       LAT_FVALUE              "GEO_geo.h"
       LAT_IDX	               "GEO_global_arrays.h"
       LATITUDE                "GEO_product.h"
       LON_IDX	               "GEO_global_arrays.h"
       LONG_FVALUE             "GEO_geo.h"
       LONGITUDE               "GEO_product.h"
       MAPIOK                  "mapi.h"
       MAX_FRAMES              "GEO_geo.h"
       MAX_SCAN_NUMBER         "GEO_geo.h"
       MAX_UINT16_VAL          "GEO_geo.h"
       NO_ELLIPSE_INTERSECT    "GEO_geo.h"
       PGS_S_SUCCESS           "PGS_SMF.h"
       MODIS_E_BAD_INPUT_ARG   "PGS_MODIS_35251.h"
       MODIS_E_GEO             "PGS_MODIS_35251.h"
       RAD2DEG                 "GEO_geo.h"
       RANGE_FVALUE            "GEO_geo.h"
       range_index             "GEO_global_arrays.h"
       SCAN_GRP                "GEO_product.h"
       SCAN_OFFSET             "GEO_product.h"
       SEN_AZIMUTH             "GEO_product.h"
       SEN_ZENITH              "GEO_product.h"
       SENSORAZIM_FVALUE       "GEO_geo.h"
       SENSORZEN_FVALUE        "GEO_geo.h"
       SOL_AZIMUTH             "GEO_product.h"
       SOL_ZENITH              "GEO_product.h"
       SOLARAZIM_FVALUE        "GEO_geo.h"
       SOLARZEN_FVALUE         "GEO_geo.h"
       TRACK_OFFSET            "GEO_product.h"
       z_angle_index           "GEO_global_arrays.h"
 
Called by: 
       GEO_write_one_scan
 
Routines Called:
       putMODISarray           "mapi.h"
       modsmf                  "smfio.h"

!Revision History:
	$Log: GEO_write_scan_data.c,v $
	Revision 6.7  2014/08/01 21:38:36  jkuyper
	Corrected landsea and waterpresent dimensions to actual maximum value.

	Revision 6.6  2011/03/08 20:21:28  kuyper
	Corrected bug introduced in 6.1, by writing sensor and solar angles in
	  degrees, rather than radians.

	Revision 6.5  2011/02/18 22:03:04  kuyper
	In order to Resolve feature request Bug 3446, changed to get landsea mask
	  directly, rather than through swath_elem, and to also get waterpresent.
	Corrected declaration of parameters that point at input data.

	Revision 6.4  2010/12/15 23:19:20  kuyper
	Changed to use *_IDX instead of *_index.
	Corrected to set fillvalue for Height along with Latitude and Longitude.

	Revision 6.3  2010/06/18 20:49:36  kuyper
	Corrected message about invalid input arguments.

	Revision 6.2  2009/06/12 16:21:48  kuyper
	Corrected to write the entire scan of each high-resolution offset SDS.

	Revision 6.1  2009/05/28 21:29:20  xgeng
	Changed MAX_SCAN_SAMPLE and MAX_DETECTORS to MAX_FRAMES and DETECTORS_1KM.
	Changed num_samples to num_frames.
	Changed to write high-resolution offsets, controlled by num_detectors;
	  all other SDS now use DETECTORS_1KM, instead.
	Changed to reorganize and convert data at one time, avoiding redundant copy
	  that used to be stored in swath_elem.
	Changed to return status codes.
	The above code changes match with PDL revision 6.1 and 6.2.

        Xu Geng (xu.geng@saic.com)
 
	Revision 4.3  2004/03/11 00:08:29  kuyper
	Corrected macro name.

	Revision 4.2  2003/10/30 23:52:40  kuyper
	Corrected scaling of range.

	Revision 4.1  2003/08/15 20:01:50  vlin
	Resolve DDTs MODxl01751: Don't scale fill values.

	Revision 3.1  2002/06/13 22:57:47  kuyper
	Removed unnecessary NCSA acknowledgement.

	Revision 2.4  1998/02/08 22:31:59  jjb
	Merged from V2.0 DAAC Delivery.

		5/25/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Finished coding.

Requirements:
        PR03-F-4.3-1
        PR03-F-4.3-1
        PR03-I-1
        PR03-I-3
        PR03-S-1
        PR03-S-2 

!Team-unique Header:

        This software is developed by the MODIS Science Data Support
        Team for the National Aeronautics and Space Administration,
        Goddard Space Flight Center, under contract NAS5-32373.

References and Credits:    None
 
Design Notes:              None

!END
*******************************************************************************/
{
  int   obj, pixel, scanline, frame;
  int   objs = 0, ier, N_samp, EV_frames;
  PGSt_SMF_status ret_val = PGS_S_SUCCESS;
  uint8 gflags[DETECTORS_1KM*MAX_FRAMES]; 
  int16 height[DETECTORS_1KM*MAX_FRAMES];
  int16 sensorazim[DETECTORS_1KM*MAX_FRAMES]; 
  int16 sensorzen[DETECTORS_1KM*MAX_FRAMES];
  int16 solarazim[DETECTORS_1KM*MAX_FRAMES];
  int16 solarzen[DETECTORS_1KM*MAX_FRAMES];
  uint16 range[DETECTORS_1KM*MAX_FRAMES];
  static long int dims[2] = {DETECTORS_1KM, MAX_FRAMES};
  static long int start[2] = {0L};
  static long int hires_start[2] ; 
  static long int hires_dims[2]; 
  int8 scan_offset[DETECTORS_QKM*SAMPLES_QKM]; 
  int8 track_offset[DETECTORS_QKM*SAMPLES_QKM];
  int8 height_offset[DETECTORS_QKM*SAMPLES_QKM];
  uint8 landsea[DETECTORS_1KM*MAX_FRAMES];
  uint8 waterpresent[DETECTORS_1KM*MAX_FRAMES];
  float32 latit[DETECTORS_1KM*MAX_FRAMES];
  float32 longit[DETECTORS_1KM*MAX_FRAMES];
  char msgbuf[PGS_SMF_MAX_MSGBUF_SIZE];
  char filefunc[] = __FILE__ ", GEO_write_scan_data";

  struct {
     char *name;  /*text string of the product file SDS object to write to.  */
     void *data;  /* data to write to the SDS. */ 
     long *start;  /*Pointer to array of start indices */
     long *dims;   /*Pointer to array of dimensions */
  } Scan_metadata[] = {
     {LATITUDE, NULL, start, dims},
     {LONGITUDE, NULL, start, dims},
     {HEIGHT, NULL, start, dims},
     {SEN_AZIMUTH, NULL, start, dims}, 
     {SEN_ZENITH, NULL, start, dims}, 
     {SOL_AZIMUTH, NULL, start, dims},
     {SOL_ZENITH, NULL, start, dims}, 
     {RANGE, NULL, start, dims}, 
     {LAND_SEAMASK, NULL, start, dims}, 
     {WATER_PRESENT, NULL, start, dims},
     {GFLAGS, NULL, start, dims}, 
     {SCAN_OFFSET, NULL, hires_start, hires_dims}, 
     {TRACK_OFFSET, NULL, hires_start, hires_dims}, 
     {HEIGHT_OFFSET, NULL, hires_start, hires_dims}, 
  }; 

  if (geo_file == NULL || swath_elem == NULL ||
      GEO_param == NULL || terrain_frame_position == NULL || 
      frame_to_sensor == NULL || frame_solar_angles == NULL ||
      frame_flags == NULL || frame_landsea == NULL ||
      frame_waterpresent == NULL || hires_offsets == NULL || 
      scan_number < 0 || scan_number >= MAX_SCAN_NUMBER)
  {
      sprintf(msgbuf, "swath_elem = %p, geo_file = %p, scan_number = %d\n"
	  "\tterrain_frame_position = %p, frame_to_sensor = %p\n"
	  "\tframe_solar_angles = %p, frame_flags = %p frame_landsea=%p\n"
	  "\frame_waterpresent = %p, thires_offsets = %p",
	  (void*)swath_elem, (void *)geo_file, scan_number,
	  (void*)terrain_frame_position, (void*)frame_to_sensor,
	  (void*)frame_solar_angles, (void*)frame_flags, (void*)frame_landsea,
	  (void*)frame_waterpresent, (void*)hires_offsets);
      modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
      return MODIS_E_BAD_INPUT_ARG;
  }

  N_samp = GEO_param->geometry_params.N_samp[GEO_param->geometry_params.band_number];
  if (N_samp < 1 || N_samp > 4 || swath_elem->num_frames < 0 
     || swath_elem->num_frames > MAX_FRAMES 
     || swath_elem->num_detectors <= 0 
     || swath_elem->num_detectors > DETECTORS_1KM) { 
     sprintf(msgbuf, "N_samp=%d, num_frames = %d, num_detectors = %d\n",
	     N_samp, swath_elem->num_frames, swath_elem->num_detectors);
     modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
     return MODIS_E_BAD_INPUT_ARG;
  }    

  if (swath_elem->num_frames == 0)
     return PGS_S_SUCCESS;

  start[0] = (long)(scan_number * DETECTORS_1KM);
  dims[1] = (long)swath_elem->num_frames;
  EV_frames = swath_elem->num_frames; 

  for(scanline = 0; scanline < DETECTORS_1KM; scanline++)
  {
    memcpy(&gflags[EV_frames*scanline], frame_flags[scanline], EV_frames);
    for(frame = 0; frame < EV_frames; frame++)
    {
      pixel = EV_frames*scanline+frame;  
      latit[pixel] = terrain_frame_position[scanline][frame][LAT_IDX]*RAD2DEG; 
      longit[pixel] = terrain_frame_position[scanline][frame][LON_IDX]*RAD2DEG; 
      height[pixel] = (int16)floor(terrain_frame_position[scanline][frame][HT_IDX]+0.5); 
      landsea[pixel] = frame_landsea[scanline][frame];
      waterpresent[pixel] = frame_waterpresent[scanline][frame];

      if(frame_to_sensor[scanline][frame][azimuth_index] == SENSORAZIM_FVALUE)
        sensorazim[pixel] = SENSORAZIM_FVALUE;
      else sensorazim[pixel] = (int16)floor(
        frame_to_sensor[scanline][frame][azimuth_index]*RAD2DEG/GEO_param->angle_scale + 0.5);

      if(frame_to_sensor[scanline][frame][z_angle_index] == SENSORZEN_FVALUE)
        sensorzen[pixel] = SENSORZEN_FVALUE ;  
      else sensorzen[pixel] = (int16)floor( 
        frame_to_sensor[scanline][frame][z_angle_index]*RAD2DEG/GEO_param->angle_scale+0.5); 

      if(frame_to_sensor[scanline][frame][range_index] /
            GEO_param->range_scale > MAX_UINT16_VAL)
        range[pixel] = RANGE_FVALUE;
      else range[pixel] = (uint16)floor(  
        frame_to_sensor[scanline][frame][range_index]/GEO_param->range_scale + 0.5);
      
      if(frame_solar_angles[scanline][frame][azimuth_index] == SOLARAZIM_FVALUE)
	solarazim[pixel] = SOLARAZIM_FVALUE;
      else solarazim[pixel] = (int16)floor(
        frame_solar_angles[scanline][frame][azimuth_index]*RAD2DEG/GEO_param->angle_scale + 0.5);  

      if(frame_solar_angles[scanline][frame][z_angle_index] == SOLARZEN_FVALUE)
        solarzen[pixel] = SOLARZEN_FVALUE ;
      else solarzen[pixel] = (int16)floor(
        frame_solar_angles[scanline][frame][z_angle_index]*RAD2DEG/GEO_param->angle_scale + 0.5);

      if (gflags[pixel] & (INVALID_INPUT_DATA | NO_ELLIPSE_INTERSECT))
      {
         latit[pixel] = (float32)LAT_FVALUE;
         longit[pixel] = (float32)LONG_FVALUE;
	 height[pixel] = (int8)HGHT_FVALUE;
      } 

    }
  }

  if (!swath_elem->lat_qaflag){
      Scan_metadata[objs].data = latit;
  } 
  objs++;

  if (!swath_elem->lon_qaflag){
     Scan_metadata[objs].data = longit;
  }
  objs++;

  if (!swath_elem->height_qaflag){
     Scan_metadata[objs].data = height;
  }
  objs++;

  if (!swath_elem->sensorazimuth_qaflag){
     Scan_metadata[objs].data = sensorazim;
  }
  objs++;

  if (!swath_elem->sensorzenith_qaflag){
    Scan_metadata[objs].data = sensorzen;
  }
  objs++;

  if (!swath_elem->solarazimuth_qaflag){
     Scan_metadata[objs].data = solarazim;
  }
  objs++;

  if (!swath_elem->solarzenith_qaflag){
     Scan_metadata[objs].data = solarzen;
  }
  objs++;

  if (!swath_elem->range_qaflag){
     Scan_metadata[objs].data = range;
  }
  objs++;

  if (!swath_elem->land_seamask_qaflag){
     Scan_metadata[objs].data = landsea;
     Scan_metadata[objs+1].data = waterpresent;
  }
  objs+=2;

  Scan_metadata[objs].data = gflags;
  objs++;

  if(N_samp > 1) {
    hires_start[0] = start[0] * N_samp; 
    hires_dims[0] = DETECTORS_1KM * N_samp; 
    hires_dims[1] = dims[1] * N_samp; 
    for(scanline = 0; scanline < hires_dims[0]; scanline++){ 
      memcpy(&scan_offset[scanline*hires_dims[1]], 
             &hires_offsets[0][scanline], hires_dims[1]); 
      memcpy(&track_offset[scanline*hires_dims[1]], 
             &hires_offsets[1][scanline], hires_dims[1]); 
      memcpy(&height_offset[scanline*hires_dims[1]], 
             &hires_offsets[2][scanline], hires_dims[1]); 
    }
    Scan_metadata[objs].data = scan_offset;
    objs++;
    Scan_metadata[objs].data = track_offset;
    objs++;
    Scan_metadata[objs].data = height_offset;
    objs++;
  }

  for (obj = 0; obj < objs; obj++){
     if (Scan_metadata[obj].data != NULL) {
        ier = putMODISarray(geo_file, Scan_metadata[obj].name, SCAN_GRP, 
              Scan_metadata[obj].start, Scan_metadata[obj].dims, 
              Scan_metadata[obj].data);
        if (ier != MAPIOK){
           sprintf(msgbuf, "putMODISarray(%s, %s)", geo_file->filename, 
                   Scan_metadata[obj].name);
           modsmf(MODIS_E_GEO, msgbuf, filefunc);
           ret_val = MODIS_E_GEO;
        }
     }
  }

  return ret_val;
}
