/*  GEO_cumulate_GRing() and GEO_get_GRing_points() */

#include <float.h>
#include "GEO_input.h"
#include "GEO_earth.h"
#include "imsl_wrap.h"
#include "PGS_SMF.h"
#include "PGS_MET.h"
#include "PGS_CSC.h"
#include "GEO_geo.h"
#include "PGS_MODIS_35251.h"
#include "smfio.h"

static double tmin=DBL_MAX;
static double smin=DBL_MAX;	/* Upper and lower bounds on*/
static double tmax=-DBL_MAX;
static double smax=-DBL_MAX;	/* projected points. */
static PGSt_double project[3][3];		/* Projection axes. */

PGSt_SMF_status GEO_cumulate_GRing(
	GEO_param_struct const		* const params,
	int32 const			num_frames,
	frame_state_struct const	sc_ev_frame_state[],
	unsigned char			frame_flags[DETECTORS_1KM][MAX_FRAMES],
	double			ecr_frame_position[DETECTORS_1KM][MAX_FRAMES][3]
)

/**************************************************************************
!C  

!Description:  Routine in the output group of the Level-1A geolocation 
               software to cumulate the information needed to determine 
               the GRing for a granule.
  
!Input Parameters:
    params	        Geolocation parameters
    num_frames	        number of frames in this scan of data
    sc_ev_frame_state	An array containing, among other things, the
  			spacecraft position at the time of each frame.
    frame_flags		An array indicating any problems that may have
  			occurred while geolocating a pixel.
    ecr_frame_position	The ECR coordinates for each pixel.
  
!Output Parameters:
      None
  
Return Parameters:
	MODIS_E_BAD_INPUT_ARG	If there's any problem with the inputs
	MODIS_E_GEO			If the matrix inversion fails
	PGS_S_SUCCESS		Otherwise
  
Externally Defined:
	DBL_EPSILON		<float.h>
  	DBL_MAX			<float.h>
  	DETECTORS_1KM		"GEO_geo.h"
	IMSL_INVERSE_USER	"imsl.h"
	IMSL_INVERSE_ONLY	"imsl.h"
  	MAX_FRAMES		"GEO_geo.h"
  	MODIS_E_BAD_INPUT_ARG	"PGS_MODIS_35251.h"
	MODIS_E_GEO		"PGS_MODIS_35251.h"
  	NO_ELLIPSE_INTERSECT	"GEO_geo.h"
        PGS_S_SUCCESS           "PGS_SMF.h"
  
Called by:
  	GEO_locate_one_scan
  
Routines Called:
	imsl_d_lin_sol_gen
	imsl_error_code
  	modsmf
  	PGS_CSC_crossProduct
	PGS_CSC_dotProduct
  	PGS_CSC_Norm
	PGS_CSC_quatRotate
  
!Revision History:
 * $Log: GEO_get_GRing_points.c,v $
 * Revision 6.2  2010/06/18 21:10:37  kuyper
 * Removed leading 'const' qualifiers from parameters that are pointers to
 *   arrays.
 * Corrected error message for invalid input arguments.
 *
 * Revision 6.1  2009/05/31 20:12:27  kuyper
 * Changes to GEO_cumulate_GRing():
 * Use DETECTORS_1KM and MAX_FRAMES where needed.
 * Changed names with "sample" to use "frame" instead.
 * Removed validation of params->num_detectors; no longer needed.
 * Expanded message buffer.
 *
 * Revision 4.3  2003/10/27 01:07:12  vlin
 * expanded message buffer size.
 *
 * Revision 4.2  2003/08/01 18:22:42  kuyper
 * Corrected error message format.
 *
 * Revision 4.1  2003/07/28 15:19:34  kuyper
 * Changed to use imsl_wrap.h.
 * Relaxed radius test for division by zero.
 * Corrected test for an error return from imsl_d_lin_sol_gen().
 * Corrected to cast pointer arguments for "%p" to (void*).
 *
 * Revision 3.8  2001/12/06 16:11:35  kuyper
 * Relaxed division by zero guards.
 *
  
Requirements:
  	PR03-F-4.1-1
  	PR03-F-4.4-1
  	PR03-I-1
  	PR03-I-2
  	PR03-S-1
  
!Team-unique Header:

  	This software is developed by the MODIS Science Data Support
  	Team for the National Aeronautics and Space Administration,
  	Goddard Space Flight Center, under contract NAS5-32373.
  
References and Credits
  
!END**************************************************************************/

#define T_HALF	150.0	/* Half the length of a normal granule (seconds) */
{
  static char	filefunc[] = __FILE__ ", GEO_cumulate_GRing";
  static double	inverse[3][3];	/* Inverse of project[][]. */

  PGSt_SMF_status	retval=PGS_S_SUCCESS;	/* Value to be returned.*/
  int		det, frame;	/* detector number, frame number	*/
  int		endframe;	/* Frame number of last good pixel.	*/
  PGSt_double	positionECR[3];	/* Spacecraft position.			*/
  PGSt_double	velocityECR[3];	/* Spacecraft velocity.			*/
  PGSt_double	quat[4];	/* Rotation quaternion.			*/
  PGSt_double	temp[3];	/* Temporary storage for a vector.	*/
  double	radius;		/* magnitude of the positionECR vector	*/
  double	pcrossv;	/* magnitude of the cross-product of the*/
				/* positionECR and velocityECR vectors	*/ 
  double	theta;		/* Rotation angel for T_HALF.		*/
  double	length=0.0;	/* Length of a vector.			*/
  double	matrix[3][3];	/* Projection axes.			*/
  double	rts[3];		/* projected components of ecr_frame_position*/
  int		row, col;
  char		msgbuf[PGS_SMF_MAX_MSGBUF_SIZE]="";
  
  if (params == NULL || sc_ev_frame_state == NULL || frame_flags == NULL
      || ecr_frame_position == NULL || num_frames > MAX_FRAMES)
  {
      sprintf(msgbuf, "params:%p, num_frames:%ld, sc_ev_frame_state:%p, "
	  "frame_flags:%p, ecr_frame_position:%p",
	  (void*)params, (long)num_frames, (void*)sc_ev_frame_state,
	  (void*)frame_flags, (void*)ecr_frame_position);
      modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);

      return MODIS_E_BAD_INPUT_ARG;
  }

  for (frame = 0; frame < num_frames; frame++)
  {
    for (det = 0; det < DETECTORS_1KM; det++)
    {
      if (frame_flags[det][frame] < (unsigned char) NO_ELLIPSE_INTERSECT) 
	/* the pixel was successfully geolocated */
      {
	if (tmin >= DBL_MAX)
	{	/*  this is first pixel processed for this granule. */
	  /*  Determine the projection matrix  */

	  for (col=0; col<3; col++) {
	     positionECR[col] =
	       (PGSt_double)sc_ev_frame_state[frame].positionECR[col]; 
	     velocityECR[col] =
	       (PGSt_double)sc_ev_frame_state[frame].velocityECR[col];
	  }
	  PGS_CSC_crossProduct(positionECR, velocityECR, quat+1);
	  radius = PGS_CSC_Norm(positionECR);
	  pcrossv = PGS_CSC_Norm(quat+1);

	  if (radius < 0.5 * params->orbit_valid_params.position_mag_limit[0] ||
	  pcrossv < params->orbit_valid_params.ang_mom_limit[0]) 
	  {
	     sprintf(msgbuf, "sc_ev_frame_state[%d].positionECR too small",
		frame);
	     modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);

	     retval = MODIS_E_BAD_INPUT_ARG;
	     break;	/* break out of detector loop. This problem will affect
			 * every pixel in same frame. */
	  }

	  memcpy(project[1], velocityECR, sizeof(project[1]));

	  /* Search for last good frame in same scanline */
	  for(endframe=num_frames-1;
	    endframe>frame && frame_flags[det][endframe]>=NO_ELLIPSE_INTERSECT;
	    endframe--);
	  if(endframe > frame)
	  {	/* Normal case: found. */
	    for (col=0; col<3; col++)
	      project[2][col] = ecr_frame_position[det][endframe][col]
			      - ecr_frame_position[det][frame][col];
	    /* Calculate normal to projection plane from scan and track */
	    PGS_CSC_crossProduct(project[1], project[2], project[0]);
	    length = PGS_CSC_Norm(project[0]);
	  }
	  if(length < DBL_EPSILON)
	  {	/* length was either not calculated, or is too short. */
	    /* Either there was exactly one good pixel in the scanline, or the
	     * yaw was exactly 90 degrees. Either case should be pretty rare.
	     * Substitute the rotation axis for the scan direction. */
	    memcpy(project[0], positionECR, sizeof(project[0]));
	    memcpy(project[2], quat+1, sizeof(project[2]));	
	  }
	  else if(PGS_CSC_dotProduct(project[0], positionECR, 3) < 0.0)
	  {	/* spacecraft flying backwards; unlikely but possible.	*/
	    for (col=0; col<3; col++)
	    {	/* Reverse two axes, to maintain intended orientation.	*/
	      project[0][col] *= -1.0;	
	      project[2][col] *= -1.0;	
	    }
	  }

	  /* Rotate the vectors forward corresponding to T_HALF seconds of
	   * orbital rotation. */
	  theta = T_HALF*pcrossv/(radius*radius);
	  quat[0] = cos(0.5*theta);
	  for (col=1; col<4; col++)
	      quat[col] *= sin(0.5*theta)/pcrossv;
	  for(row=0; row<3; row++)
	  {
	    if(PGS_CSC_quatRotate(quat, project[row], temp)!=PGS_S_SUCCESS)
	    {
	      modsmf(MODIS_E_GEO, "PGS_CSC_quatRotate()", filefunc);
	      break;
	    }
	    length = PGS_CSC_Norm(temp);
	    if(length < DBL_EPSILON)
	    {
	      modsmf(MODIS_E_GEO, "PGS_CSC_Norm()", filefunc);
	      break;
	    }
	    for(col=0; col<3; col++)
	      project[row][col] = temp[col]/length;
	  }
	  if(row<3)
	    continue; 
	  
	  for (row=0; row<3; row++)
	    for (col=0; col<3; col++)
	      matrix[row][col] = (double)project[row][col];
	  
	  /* Invert matrix. */
	  imsl_d_lin_sol_gen( 3, matrix[0], NULL,
	    IMSL_INVERSE_USER, inverse[0],
	    IMSL_INVERSE_ONLY,
	    0);
	  if(imsl_error_code())
	  {	/* Matrix singular or not positive definite */
	    modsmf(MODIS_E_GEO, "imsl_d_lin_sol_gen()", filefunc);
	    continue;
	  }
	}  /* if tmin >= DBL_MAX ended */

	for (row=0; row<3; row++)
	{
	  rts[row] = 0.0;
	  for (col=0; col<3; col++)
	    rts[row] += inverse[col][row] *
		ecr_frame_position[det][frame][col];
	}

	if (rts[0] > (0.01 * params->orbit_valid_params.position_mag_limit[0]))
	{
	  rts[1] /= rts[0];
	  rts[2] /= rts[0];
	  tmin = (tmin < rts[1]) ? tmin : rts[1];
	  tmax = (tmax > rts[1]) ? tmax : rts[1];
	  smin = (smin < rts[2]) ? smin : rts[2];
	  smax = (smax > rts[2]) ? smax : rts[2];
	}
	else
	{
	  sprintf(msgbuf,"ecr_frame_position[%d][%d]: rts[0]=%g",
	      det, frame, rts[0]);
	  modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
	  retval = MODIS_E_BAD_INPUT_ARG;
	}
      }	/* end of if the pixel was successfully geolocated */
    }	/* for det loop ended */
  }	/* for frame loop ended */
  return retval;
}

/******************************************************************************/

PGSt_SMF_status GEO_get_GRing_points(
	GEO_GRing_struct	* const GRing_points)

/************************************************************************ 
!C

!Description:   
   	Routine in Output group of the Level-1A geolocation software to
  	calculate the G-Ring polygon vertices for the ECS Core Metadata.
   	At this point, GEO_cumulate_GRing() has already projected each ground
   	location from the center of the earth onto a plane that's roughly
  	oriented with the granule. It has determined the bounding box in that
  	plane containing all of those projected points. This routine projects
  	the corners of that bounding box toward the center of the earth,
  	stopping at the surface of the earth. As result, the great circle arcs
  	that connect those points correspond precisely to the edges of the
  	bounding box, and will therefore enclose all of the ground locations.
   
!Input Parameters:
   	None
   
!Output Parameters:
   	GRing_points         G-Ring polygon latitudes and longitudes
   
Return Parameters:
  	MODIS_E_BAD_INPUT_ARG   If GRing_points is NULL
  	MODIS_E_GEO		If projection failed for any vertex
  	MODIS_W_NO_GEO		If no geolocatable pixels were found
  	PGS_S_SUCCESS		If projection was successful
   
Global variables: None
   
Externally Defined:
  	DBL_MAX			<float.h>
  	MODIS_E_GEO		"PGS_MODIS_35251.h"
  	MODIS_E_BAD_INPUT__ARG	"PGS_MODIS_35251.h"
  	MODIS_W_NO_GEO		"PGS_MODIS_35251.h"
  	EARTH_MODEL		"GEO_earth.h"
  	PGS_S_SUCCESS		"PGS_SMF.h"
  	PGSCSC_W_HIT_EARTH	"PGS_CSC.h"
  	RAD2DEG			"GEO_geo.h"
  
Called by:
   	GEO_write_granule_metadata
   
Routines Called:
  	modsmf			Logs a status message
  	PGS_CSC_GrazingRay	Projects from a point along a line of 
                                sight to the surface of the Earth.
  	PGS_CSC_ECRtoGEO	Converts from ECR to Geodetic coordinates.
   
!Revision History:
        Please see prologue in "GEO_cumulate_GRing.c"
  
Requirements:
  	PR03-F-4.1-1
  	PR03-F-4.4-1
  	PR03-I-1
  	PR03-I-2
  	PR03-S-1
  
!Team-unique Header:
        This software is developed by the MODIS Science Data Support
        Team for the National Aeronautics and Space Administration,
        Goddard Space Flight Center, under contract NAS5-32373.                 

!END
***************************************************************************/

{
#define CORNERS	4
#define BORDER	1e-5     /* roughly 1/2 the angular width in radians 
                           of a pixel at nadir */
  PGSt_double corner_point[CORNERS][3], posECR[3];
  PGSt_double length, ray[3];
  PGSt_double missAlt, slantRg;
  PGSt_double posNEAR[3], posSURF[3];
  int cor, j;
  GEO_GRing_struct new_points;
  PGSt_SMF_status rtnfun;

  if (GRing_points == NULL) 
  {
     modsmf(MODIS_E_BAD_INPUT_ARG,".",
     "GEO_get_GRing_points.c, GEO_get_GRing_points");
     return MODIS_E_BAD_INPUT_ARG;
  }

  new_points = *GRing_points;

  if (tmin >= DBL_MAX)
  {
    /* There were no geolocatable pixels in the granule. This is a normal
     * occurrance in the face of bad data or an empty granule; not a defect in
     * this program.*/
    modsmf(MODIS_W_NO_GEO,"", "GEO_get_GRing_points.c, GEO_get_GRing_points");
    return MODIS_W_NO_GEO;
  }

/* The following adjustments guarantee that even if there's only one
   geolocatable pixel in the granule, the GRing will cover a finite area - 
   GRings with zero area are not permitted by the archiving system. */

   tmax += BORDER;
   smax += BORDER;
   tmin -= BORDER;
   smin -= BORDER;

   for (j=0; j<3; j++) {
   corner_point[0][j] = project[0][j] + tmin*project[1][j] + smin*project[2][j];
   corner_point[1][j] = project[0][j] + tmin*project[1][j] + smax*project[2][j];
   corner_point[2][j] = project[0][j] + tmax*project[1][j] + smax*project[2][j];
   corner_point[3][j] = project[0][j] + tmax*project[1][j] + smin*project[2][j];
   }

   for (cor = 0; cor < CORNERS; cor++)
   {
     length = 1.0;
     for (j=0; j<3; j++) {
	 ray[j] = -corner_point[cor][j];
	 posECR[j] = 7.0e6 * corner_point[cor][j];
	 if (fabs(ray[j]) >= length) 
	     length = fabs(ray[j]);
     }

     if (length > 1.0)
	 for (j=0; j<3; j++) 
	     ray[j] /=length; 

      rtnfun = PGS_CSC_GrazingRay(EARTH_MODEL, posECR, ray, 
	  &new_points.latitude[cor], &new_points.longitude[cor],
	  &missAlt, &slantRg, posNEAR, posSURF);  
			      /* project along a ray parrallel to the corner's 
				 position vector, up to the earth's surface, 
				 starting at the corner's position */

      if (rtnfun != PGSCSC_W_HIT_EARTH)
      {
         modsmf(MODIS_E_GEO,"PGS_CSC_GrazingRay()",
         "GEO_get_GRing_points.c, GEO_get_GRing_points");
         return MODIS_E_GEO;
      }

      rtnfun = PGS_CSC_ECRtoGEO(posSURF, EARTH_MODEL,
	  &new_points.longitude[cor], &new_points.latitude[cor], &missAlt);

      if (rtnfun != PGS_S_SUCCESS)
      {
         modsmf(MODIS_E_GEO,"PGS_CSC_ECRtoGEO()",
         "GEO_get_GRing_points.c, GEO_get_GRing_points");
         return MODIS_E_GEO;
      }               
      new_points.latitude[cor] *= RAD2DEG;
      new_points.longitude[cor] *= RAD2DEG;

   }    /* for cor loop ended */

  tmin=DBL_MAX;
  smin=DBL_MAX;
  tmax=-DBL_MAX;
  smax=-DBL_MAX; 

  *GRing_points = new_points;

  return PGS_S_SUCCESS;
}
