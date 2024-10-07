#include "PGS_CSC.h"
#include "PGS_MODIS_35251.h"
#include "GEO_global_arrays.h"
#include "GEO_earth.h"
#include "GEO_util.h"
#include "smfio.h"

PGSt_SMF_status GEO_terrain_correct(
	int const sample_number,
	int const num_detectors,
        double  sample_view_vec[][MAX_PADDED][3],
        double  ecr_sample_position[][MAX_PADDED][3],
        double  ellip_sample_position[][MAX_PADDED][3],
        unsigned char sample_flags[][MAX_PADDED],
        double  terrain_sample_position[][MAX_PADDED][3]
	)
/*
!C******************************************************************************
!Description:   
		subroutine in Earth Location Group of the Level-1A geolocation
                software to perform terrain correction.  The algorithm
		establishes a regular set of points along the line 
		of sight from the spacecraft to the ellipsoid position
		of each pixel.  The geographic heights are interpolated
		from the DEM data to identify fhe line segment in which
		terrain pierce occurs.  A further (linear) interpolation
		is then performed using the DEM data points bounding pierce
		to estimate the position and height of the pixel.

!Input Parameters:
                int sample_number     -   the sample number
		int num_detectors     -   the number of detectors
                sample_view_vec       -   ECR sample view vec
                ecr_sample_position   -   ECR sample position to the ellipsoid
                ellip_sample_position -   ellip sample position
                sample_flags          -   geolocation pixel data status bits

!Output Parameters:
                terrain_sample_position  -  terrain corrected sample position 

Return values:
                MODIS_E_BAD_INPUT_ARG   If any argument is invalid
                MODIS_E_GEO             If any subroutine failed.
                PGS_S_SUCCESS           Otherwise

Externally Defined:
                BAD_TERRAIN             "GEO_geo.h"
                DETECTORS_QKM           "GEO_geo.h"
                EARTH_MODEL             "GEO_earth.h"
                HT_IDX                  "GEO_global_arrays.h"
                LAT_IDX                 "GEO_global_arrays.h"
                LON_IDX                 "GEO_global_arrays.h"
                MODIS_E_BAD_INPUT_ARG   "PGS_MODIS_35251.h"
                MODIS_E_GEO             "PGS_MODIS_35251.h"
                MODIS_E_UTIL_VEC        "PGS_MODIS_35251.h"
                NEAR_LIMB               "GEO_geo.h"
                NO_ELLIPSE_INTERSECT    "GEO_geo.h"
                MAX_PADDED              "GEO_geo.h"
                PGS_S_SUCCESS           "PGS_SMF.h"
                PGSd_PC_LINE_LENGTH_MAX "PGS_PC.h"
                SUCCESS                 "GEO_basic.h"
                TERRAIN_CORRECT         "GEO_geo.h"
                TERRAIN_CORRECT_LUN     "GEO_geo.h"

Call functions:
		GEO_read_DEM - read in tile containing specified location
		GEO_latlon2height - read DEM get height from lat lon
                PGS_CSC_dotProduct - compute dot product of two vectors 
                PGS_CSC_Norm - perform vector normalization   
		PGS_PC_GetConfigData - Get parameter from PCF
		PGS_CSC_ECRtoGEO - Convert ECR coordinates to geodetic ones.
		modsmf - writes error status messages to log
			
Called by:	GEO_earth_location

!Requirements:
                PR03-F-3.3.1-1


!Revision History:
 * $Log: GEO_terrain_correct.c,v $
 * Revision 6.6  2010/12/14 22:02:06  kuyper
 * Changed to perform all calculations using the height above the ellipsoid.
 *
 * Revision 6.5  2010/06/24 21:36:49  kuyper
 * Corrected geoid height adjustment to match PDL.
 *
 * Revision 6.4  2010/03/31 19:41:44  kuyper
 * Resolved Bug 2937 by always subtracting the geoid_height, regardless of
 *   how close the terrain-corrected position comes to the ellipsoid.
 *
 * James Kuyper		james.kuyper@sigmaspace.com
 *
 * Revision 6.3  2009/09/15 16:34:34  kuyper
 * Corrected calculation of cos_view.
 * Corrected handling of geoid_height for locations near the ellipsoid.
 *
 * James Kuyper		James.R.Kuyper@nasa.gov
 *
 * Revision 6.2  2009/05/29 13:34:00  xgeng
 * Corrected a typo and an error message.
 *
 * Revision 6.1  2009/05/18 17:21:52  xgeng
 * Changed return type to PGSt_SMF_status.
 * Renamed pixel_flags to sample_flags.
 * Changed sample_view_vec, ecr_sample_position, ellip_sample_position, and
 * pixel_flags->sample_flags into input pointer parameters.
 * Changed terrain_sample_position to an output pointer parameter.
 * Changed MAX_SCAN_SAMPLE to MAX_PADDED, and MAX_DETECTORS to
 * DETECTORS_QKM.
 * Added validation for new pointer parameters.
 * Replaced calls to GEO_vec* functions with calls to PGS_CSC_dotProduct() and PGS_CSC_Norm().
 * Cleaned up error messaging.
 * The above changes match with PDL revision 6.1 to 6.3.
 *
 *
 * Xu Geng (xu.geng@saic.com)
 *
 * Revision 4.1  2003/02/21 22:13:15  kuyper
 * Corrected to ensure 'D' is initialized before use.
 *
 * Revision 3.1  2001/11/30 17:31:34  kuyper
 * Corrected comments to match code.
 *
 * Revision 2.10  1999/03/01 20:28:14  kuyper
 * Changed to protect against overly large step sizes, including division by 0.
 *
 * Revision 2.9  1999/02/16  17:57:56  lma
 *  Additional modification during unit testing.
 * Changed calculation of num_steps.
 * removed "if(num_steps>2)" branches
 * corrected initialization of h_prime and
 * corrected calculation of D during interpolation
 *
 * Revision 2.8  1999/02/08  19:12:26  kuyper
 * Changed implicit type conversions to explicit ones.
 *
 * Revision 2.7  1999/02/02  16:17:50  lma
 * Corrected descriptions of hgtmin, hgtmax.
 * Inserted correction for cos_view>1.0, from code.
 * Corrected calculation of delta_D in the 2-step case.
 * Changed design to avoid excess work in the case where hgtmin is close to
 *  hgtmax.
 * Corrected application of geoid_height to final output.
 *
 * Revision 2.6  1998/03/04  00:14:11  jjb
 * CCR-391: Replaced IMSL matrix routine with GEO_location routine.
 * Made search_grid a static constant.
 * Added 'called by' line.
 * Added setting of BAD_TERRAIN flag
 * 	when no terrain correction selected.
 * Modified awkward phrase in 'Description'.
 * Removed reference to obsolete 'scan_number' parameter.
 *
 * Revision 2.5  1997/11/03  17:09:47  jjb@ltpmail.gsfc.nasa.gov
 * Removed possibility of delta_D = 0 causing an infinite loop.
 *
 * Revision 2.4  1997/10/30  20:14:28  jjb@ltpmail.gsfc.nasa.gov
 * Corrected (infinite) loop condition to handle all-ocean DEM tiles.
 *
 * Revision 2.3  1997/10/24  13:21:57  kuyper@ltpmail.gsfc.nasa.gov
 * Included header file which declares terrain_sample_position[][].
 * Initialized DEMflag to NO_DEM_DATA.
 *
 * Revision 2.2  1997/10/23  19:22:35  ding@ltpmail.gsfc.nasa.gov
 * Moved pixel_flags check up to right below the for loop.
 *
 * Revision 2.1  1997/10/21  18:16:22  kuyper
 * Returned from ClearCase
 *
 * Revision /main/GEO_V2_DEV/2 1997/10/03 kuyper
 * Updated call to GEO_read_DEM().
 *
 * Revision 1.9.1.1  1997/07/21  22:20:17  kuyper
 * Merged in out-of-sequence changes.
 *
 * Revision 1.9  1997/07/21  16:24:34  kuyper
 * Baselined Version 1
 *
 *Parallel development:
 * Revision 1.12  1997/04/21  22:37:32  fhliang
 * commented unused argument (LL.19-21).
 *
 * Revision 1.11  1997/04/15  19:44:51  fhliang
 * fixed prolog: 'Description' --> '!Description';
 * commented one 'printf' statement (LL.369-372).
 *
 * Revision 1.10  1997/04/10  19:36:13  fhliang
 * Initial revision of SDST re-delivery of GEO_terrain_correct.c.
 *
 * Revision 1.9  1997/04/09  19:37:11  kuyper
 * Corrected ==FAIL to !=SUCCESS, three places.
 * Simplified detection of single pass through loop.
 * #defined HGT_TOL, used consistently.
 *
 * Revision 1.8  1997/03/26  18:14:47  fhliang
 * Initial revision of SDST delivery of GEO_terrain_correct.c.
 *
		Revision 1.7  1997/03/11 18:02:14  kuyper
		Allowed for some roundoff error in loop
		  termination conditions.

		Revision 1.6  1997/02/04 20:57:35  kuyper
		Corrected, modified loop termination condition.
		
		Revision 1.5  1996/11/15 23:47:16  kuyper
		Corrected return value when TERRAIN_CORRECT_LUN not found.
		Changed earthEllipsTag to EARTH_MODEL
		Corrected handling of pixel_flags.
		Corrected cos_view error message.
		Enforced cos_view <= 1.0
		Removed redundant copy through x[i]

		Revision 1.4  1996/11/06 16:02:24  kuyper
		Use PGS_CSC_ECRtoGEO() rather than GEO_ecr2latlon().
		Use imsl_d_mat_mul_rect() rather than GEO_vec_prod3()
		Get terrain_correct from LUN parameter, not parameter file.
		Always initialize terrain_sample_position with copy of ellip_sample_position
		Validate input arguments.
		Implement minimum_height, MIN_COS_VIEW checks
		Always iterate main loop at least once.
		Rearrange delta_D check to avoid division by 0 problems.
		Added detector number to log status messages.
		Set BAD_TERRAIN if GEO_vec_unit3() fails.
		Moved SMF messages to 35256 seed file.
		Updated function call list.
		Name Changes: idet=>detector_index, nu=>cos_view, *_final=>previous_*
		Added LAT_IDX, LON_IDX, HT_IDX macros.
		Fixed typos in prologue.

		Revision 1.3  1996/07/24 21:37:10  kuyper
		Standardized order of #include files.
		Declared arguments const.
		Inserted required '!'s in comments.
		Converted constants to double, to skip implied conversion.
		Removed ret_val.
		Made implicit casts explicit.

		Revision 1.2  1996/07/18 21:08:08  kuyper
		Included self-checking header file.
		Converted terrain_sample_position from extern declaration to definition.
		Corrected call to GEO_vec_unit3, and comment about GEO_ecr2latlon.
		Added needed header files.

		10/10/95
		Tracey W. Holmes
		Added debug option.

		9/29/95
		Frederick S. Patt
		Corrected calculation of unit_n; saved hgtmin, hgtmax and
		  DEMflag as static variables; corrected bug in use of
		  h_prime and h_final_prime

		9/28/95
		Frederick S. Patt
		Correction to unit_n index.

		8/2/95
		Tracey W. Holmes (holmes@modis-xl.gsfc.nasa.gov)
		Changed extern correct_terrain from type double to int.

		7/3/95
		Tracey W. Holmes (holmes@modis-xl.gsfc.nasa.gov)
		Added SDP error messages. 

		4/25/95
		Ruiming Chen (rchen@ltpmail.gsfc.nasa.gov)
		Finished coding.

!Team-unique Header:
		This software is developed by the MODIS Science Data Support
		Team for the National Aeronautics and Space Administration,
		Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/
{

#define MIN_COS_VIEW 0.087155743
/* Minimum cosine of the LOS wrt ellipsoid normal. = cos(85.0*PGS_PI/180.0)*/

/*******************************************************
declare local variables
*******************************************************/

  double search_grid = 500.0;   /* terrain search grid*/
  double step_factor = 2.0;     /*Maximum ratio of step size to tile thickness.*/
  double unit_n[3];	/* local ellipsoid normal unit vec */
  PGSt_double u[3];	/* ECR vector from the round point to the inst. */
  double range;         /* magnitude of 'u' */ 
  PGSt_double unit_u[3];	/* unit u */ 
  double cos_view;	/* Local vertical component of satellite view vector */
  double sin_view;	/* sin of the viewing angle */
  double D;  		/* distance along satellite vec to achieve a height h */
  double D0;  		/* the starting value of D in the intersection loop */
  double D_max;		/* maximum D */
  double D_min;		/* minimum D */
  double delta_D;	/* increment of D */
  int step;		/* the index in the intersection loop */
  int num_steps;	/* the number of steps to be executed in the intersection loop */
  PGSt_double x[3];		/* ECR position of one iteration */
  double height; 	/* DEM height */
  double previous_height;	/* Previous value of height.	*/
  PGSt_double h_prime;	/* height on the search line */
  double previous_h_prime = 0.0;/* Previous value of h_prime	*/
  PGSt_double latitude;	/* Geodetic location of sample point on LOS */
  PGSt_double longitude;
  double interpolation_weight = 0.0;  /* weight */
  double minimum_height = 1.0;	/* minimum along view vector terrain deviation*/
  static int hgtmin = 0;	/* minimum height for tile */
  static int hgtmax = 0;	/* maximum height for tile */
  int det = 0;	                /* detector index */
  int i = 0;			/* iteration parameter */
  static char correct_terrain[PGSd_PC_LINE_LENGTH_MAX]="";
  char msgbuf[PGS_SMF_MAX_MSGBUF_SIZE];		/* To store 2nd argument to modsmf calls. */
  char   filefunc[] = __FILE__ ", GEO_terrain_correct";


/********************************************************
calculation from here
********************************************************/

  /* Validate arguments	*/
  if( sample_number<0 || sample_number>=MAX_PADDED || 
      num_detectors<0 || num_detectors>DETECTORS_QKM ||
      sample_view_vec == NULL || ecr_sample_position == NULL ||
      ellip_sample_position == NULL || sample_flags == NULL ||
      terrain_sample_position == NULL
    )
  {
     sprintf(msgbuf,"\nsample_number: %d, num_detectors: %d"
                    "sample_view_vec: %p, ecr_sample_position: %p"
                    "ellip_sample_position: %p, sample_flags: %p" 
                    "terrain_sample_position: %p", sample_number,num_detectors,
                    (void *)sample_view_vec, (void *)ecr_sample_position, 
	            (void *)ellip_sample_position, (void *)sample_flags, 
                    (void *)terrain_sample_position ); 		   
     modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf,filefunc);
     return MODIS_E_BAD_INPUT_ARG;
  }

  /* Initialize output to position of ellipsoid height input	*/
  for (det = 0; det < num_detectors; det++) {
     for (i = 0; i < 3; i++)
        terrain_sample_position[det][sample_number][i] = 
	ellip_sample_position[det][sample_number][i];
  }

  /* If no terrain correction required, copy data and return */
  if ((correct_terrain[0] == '\0')&&(PGS_PC_GetConfigData(
      TERRAIN_CORRECT_LUN, correct_terrain) != PGS_S_SUCCESS))
  {
     sprintf(msgbuf, "PGS_PC_GetConfigData(%ld)", (long)TERRAIN_CORRECT_LUN);
     modsmf(MODIS_E_GEO, msgbuf, filefunc);
     return MODIS_E_GEO;
  } 

  if(strncmp(correct_terrain, TERRAIN_CORRECT, sizeof(correct_terrain))!=0)
  {                	/* No terrain correction required.	*/
     for (det = 0; det < num_detectors; det++)
 	sample_flags[det][sample_number] |= BAD_TERRAIN;
     return PGS_S_SUCCESS;
  }

  /* Loop through detectors */

  for (det = 0; det < num_detectors; det++)
  {
     if(sample_flags[det][sample_number]<NO_ELLIPSE_INTERSECT)
     {
        if(GEO_read_DEM(
	ellip_sample_position[det][sample_number][LAT_IDX],
	ellip_sample_position[det][sample_number][LON_IDX],
	&hgtmin, &hgtmax) != SUCCESS)
        {
           sprintf(msgbuf, "GEO_read_DEM(%f, %f) for detector %d",
                   ellip_sample_position[det][sample_number][LAT_IDX],
                   ellip_sample_position[det][sample_number][LON_IDX], det);
	   modsmf(MODIS_E_GEO, msgbuf, filefunc);
	   sample_flags[det][sample_number] |= BAD_TERRAIN;
        }
        else 
        {
           /* calculate local ellipsoid normal unit vector */

           unit_n[0] =
	     cos(ellip_sample_position[det][sample_number][LAT_IDX])* 
	     cos(ellip_sample_position[det][sample_number][LON_IDX]);
	   unit_n[1] =
	     cos(ellip_sample_position[det][sample_number][LAT_IDX])* 
	     sin(ellip_sample_position[det][sample_number][LON_IDX]);
	   unit_n[2] =
	     sin(ellip_sample_position[det][sample_number][LAT_IDX]);

	   /* compute the ECR unit vector from the ground point to the satellite */

	   for (i = 0; i < 3; ++i)
	      u[i] = - sample_view_vec[det][sample_number][i];
     
           range = PGS_CSC_Norm(u); 
           if (fabs(range) < DBL_MIN) {
             sprintf(msgbuf, "u for detector %d", det);
             modsmf(MODIS_E_UTIL_VEC, msgbuf, filefunc);
             sample_flags[det][sample_number] |= BAD_TERRAIN;
             continue;
           }

	   /* compute the component of the satellite vector in the 
	      local vertical direction (cos_view)*/
	   for (i = 0; i < 3; ++i)
              unit_u[i] = u[i]/range; 

           cos_view = PGS_CSC_dotProduct(unit_n, unit_u, 3);
	   if(cos_view < MIN_COS_VIEW)
	   {
	      sample_flags[det][sample_number] |= (BAD_TERRAIN | NEAR_LIMB);
	      continue;
	   }
           if(cos_view > 1.0)  cos_view =1.0;
            
           /* compute D_max and D_min (distance along line-of-sight vector to 
	      hgtmax and hgtmin) */

           D_max = (double)hgtmax / cos_view;
           D_min = (double)hgtmin / cos_view;

	   /* Compute search step size along view vector */
	   sin_view = sqrt(1.0 - cos_view*cos_view);

           if( (double)(hgtmax - hgtmin) < minimum_height)
           {
              /* The terrain is so flat that interpolation is not necessary. */
              D = 0.5*(D_max + D_min);
           }
           else
           {
	      if(sin_view*(D_max-D_min) < search_grid/step_factor)
		/* protect against division by 0. */
		delta_D = (D_max-D_min)*step_factor;
	      else	/* horizontal component == search_grid */
		delta_D = search_grid / sin_view;

	      /* Must cover full range, but don't want any step to correspond
	       * exactly to one of the boundary values.
	       */
              num_steps = 1+(int)ceil( (D_max-D_min)/delta_D + 0.1);

                                        /* Initialization before
                                           iteration */
              /* D_max+delta_D/2  >  D0                        >= D_max,
               * D_min            >= D0-(num_steps-1)*delta_D  > D_min-delta_D/2
              */
              D0 = 0.5*( D_max + D_min + (double)(num_steps-1)*delta_D );
	      D = D0;

                                        /* Begin stepping down line of sight,
                                           searching for intersection with
                                           terrain. */
              height = 0.0;
              h_prime = 1.0;

              for(step = 0; (step < num_steps) && (height<h_prime) &&
                 !(sample_flags[det][sample_number] & BAD_TERRAIN); step++)
              {
                 previous_height  = height;  /* Save current values for heights */
                 previous_h_prime = h_prime;

                                        /* Decrement distance along
                                           line-of-sight */
                 D = D0 - (double)step*delta_D;
                                        /* calculate x, the ECR position
                                           vector for this initial
                                           location */
                 for (i = 0; i < 3; i++)
                    x[i] = ecr_sample_position[det][sample_number][i]
                           + D*unit_u[i];

                 /* Compute lat, lon, height */
                 if(PGS_CSC_ECRtoGEO(x, EARTH_MODEL, &longitude, &latitude, &h_prime)
                    != PGS_S_SUCCESS)
                 {
                    /* call SDP toolkit function to report error */
                    sprintf(msgbuf, "PGS_CSC_ECRtoGEO() on detector %d", det);
                    modsmf(MODIS_E_GEO, msgbuf, filefunc);
                    sample_flags[det][sample_number] |= BAD_TERRAIN;
                 }

                 /* Compute terrain height */
                 else if (GEO_latlon2height(latitude, longitude, &height) != SUCCESS)
                 {
                    /* call SDP toolkit function to report error */
                    sprintf(msgbuf, "GEO_latlon2height(%f, %f) on detector %d", 
                            latitude, longitude, det);
                    modsmf(MODIS_E_GEO, msgbuf, filefunc);
                    sample_flags[det][sample_number] |= BAD_TERRAIN;
                 }
              }

              if(sample_flags[det][sample_number] & BAD_TERRAIN)
                 continue;

              if(step > 0) /*the search loop was performed more than once */
              {
                                        /* Now (should) have line of sight
                                           and terrain heights for two
                                           locations.  Interpolate between
                                           these two positions to the
                                           estimated line of sight intersection
                                           with the terrain. */
                                        /* calculate the interpolation_weight */
                 interpolation_weight = (height-h_prime) /
                       (height - previous_height + previous_h_prime - h_prime);

                  /* calculate the final ECR coordinates */
                  D = D0 + ( interpolation_weight-(double)(step-1) ) * delta_D;
              }
           }  /* end else hgtmax-hgtmin is less than minimum height */
     
           /* Check for correction less than 1 meter */
           if (fabs(D) > minimum_height) 
           {
              for (i = 0; i < 3; ++i)
              ecr_sample_position[det][sample_number][i]
                  += D*unit_u[i];

              /* calculate final geodetic coordinates */
              if(PGS_CSC_ECRtoGEO(ecr_sample_position[det][sample_number],
                   EARTH_MODEL, &longitude, &latitude, &height) != PGS_S_SUCCESS)
              {
                 /* call SDP toolkit function to report error */
                 sprintf(msgbuf, "PGS_CSC_ECRtoGEO() on detector %d", det);
                 modsmf(MODIS_E_GEO, msgbuf, filefunc);
                 sample_flags[det][sample_number] = BAD_TERRAIN;
                 continue;
              }

              /* assign the final position to the global variable */
              terrain_sample_position[det][sample_number][LAT_IDX] =
	          latitude;
              terrain_sample_position[det][sample_number][LON_IDX] =
	          longitude;
              terrain_sample_position[det][sample_number][HT_IDX] =
	          height;
           }

        } /* End of if DEM data */

     } /* End of if pixel_flags<No_ELLIPS_INTERSECT */
    
  } /* End of detector loop */
  
  return SUCCESS;
}
