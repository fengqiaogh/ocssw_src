#include "PGS_CSC.h" 
#include "PGS_TD.h"
#include "PGS_MODIS_35251.h"
#include "GEO_earth.h"
#include "GEO_output.h"
#include "GEO_inst.h"
#include "GEO_util.h"
#include "GEO_validation.h"
#include "smfio.h"

PGSt_SMF_status GEO_solar_and_lunar_vectors(
        PGSt_double const		frame_time[],
	frame_data_struct const		* const frame_data,
	GEO_param_struct const		* const geo_params,
	fill_values_struct const 	* const fill_values,
	celestial_bodies_struct		* const cb_vectors)
/******************************************************************************
!C

!Description:   
    Subroutine in output group of the Level-1A geolocation software to compute 
    the celestial body vectors for one scan.  Outputs produced are:
          spacecraft state in ECR frame at Earth View center
          sun unit vector in ECR frame at Earth View center
          sun azimuth and zenith angle to the Solar Diffuser at Solar Diffuser 
                                          scan center
          lunar unit vector in MODIS instrument frame at Space Scan view center.
    (The Earth View center time is defined as the midpoint of the
     integration of the middle frame of the Earth View scan. If the Earth
     View has an even number of frames, the integration midpoint of the lower
     number frame adjacent to the scan center is used.)
 
!Input Parameters:  
        frame_data              scan start times and frame sizes
        geo_params              geolocation parameters read from input file
	fill_values		a list of the fill values from the L1A file
        frame_time              The time from the start of the scan to the
                                collection of data for each frame.
 
!Output Parameters:
	cb_vectors		spacecraft, solar and lunar position vectors
 
Return Values:
	MODIS_E_BAD_INPUT_ARG	If mandatory pointer arguments are NULL
	MODIS_E_GEO		If any subroutine failed.
	PGS_S_SUCCESS		Otherwise.
 
Externally Defined:
        BAD_DATA                        "GEO_validation.h"
        GEO_DOUBLE_FILLVALUE            "GEO_output.h"
        GOOD_DATA                       "GEO_validation.h"
        MAX_FRAMES                      "GEO_geo.h"
        MODIS_E_BAD_INPUT_ARG           "PGS_MODIS_35251.h"
        MODIS_E_GEO                     "PGS_MODIS_35251.h"
        PGS_FALSE                       "PGS_SMF.h"
        PGS_S_SUCCESS                   "PGS_SMF.h"
        PGSd_GEO_ERROR_VALUE            "PGS_TD.h"
        PGSd_MOON                       "PGS_CBP.h" 
        PGSd_SUN                        "PGS_CBP.h"
        SUCCESS                         "GEO_basic.h"
        TIMECODEASIZE                   "smfio.h"
 
Called by:
        GEO_write_one_scan
 
Routines Called:
        GEO_interp_ECR                  "GEO_earth.h"
        GEO_vec_unit3                   "GEO_utils.h"
        modsmf                          "smfio.h"
        PGS_CBP_Earth_CB_Vector         "PGS_CBP.h"
        PGS_CSC_ECItoECR                "PGS_CSC.h"
        PGS_SMF_TestWarningLevel        "PGS_SMF.h"

!Revision History:
   $Log: GEO_solar_and_lunar_vectors.c,v $
   Revision 6.5  2017/03/27 18:46:03  jkuyper
   Corrected to treat a return of PGSCSC_W_BAD_TRANSFORM_VALUE from
     PGS_CSC_ECItoECR() as an error.

   Revision 6.4  2010/06/29 18:48:30  kuyper
   Converted floating point equality tests into inequalities.

   Revision 6.3  2010/05/27 15:18:27  kuyper
   Changed pointers to input-only data to be pointers to const.

   Revision 6.2  2010/05/19 22:10:00  kuyper
   Helped resolve Bug 1969 by creating rpy and passing it to GEO_interp_ECR().
   Helped resolve Bug 2472 by passing params rather a couple of it's members
     to GEO_interp_ECR, as well as by creating sample_quality, and passing it
     as well.

   James Kuyper Jr.		James.R.Kuyper@NASA.gov

   Revision 6.1  2009/05/22 16:54:28  xgeng
   Change MAX_SCAN_SAMPLE to MAX_FRAMES.
   Change global array sample_time to a pointer parameter named frame_time.

   Xu Geng (xu.geng@saic.com) 

   Revision 5.4  2006/07/21 17:36:42  kuyper
   Changed order of enumeration, so that GEO_interp_ECR() will use the
     center time of the Earth view for calculating solar elevation correction,
     matching the value used for the frame-level calculations.

   Revision 5.3  2004/10/20 19:14:15  kuyper
   Corrected handling of return value from GEO_interp_ECR().

   Revision 5.2  2004/10/19 15:46:07  kuyper
   Added fill_values to list of pointer parameters needing to be checked.

   Revision 5.1  2004/09/23 15:40:26  kuyper
   Removed obsolete parameter.
   Added validity check of EV start time versus L1A fill value.
   Changed to return status codes, rather than a pointer.

   Revision 4.7  2004/04/09 18:16:13  kuyper
   Corrected to make MODIS_E_GEO returns from GEO_interp_ECR() non-fatal.

   Revision 4.6  2003/10/07 21:32:17  kuyper
   Corrected error handling of solar angles.

   Revision 4.5  2003/08/28 16:00:09  kuyper
   Corrected prolog.

   Revision 4.4  2003/08/26 20:07:47  kuyper
   Corrected return value for a bad position vector.

   Revision 4.3  2003/08/25 21:41:14  kuyper
   Changed to calculate sunSD using T_sc2SD and T_sc2ecr, rather than T_inst2SD
     and T_inst2ecr.
   Corrected to actually calculate Solar Diffuser angles.

   Revision 4.2  2003/08/19 21:27:06  kuyper
   Corrected some type-changing copies to not use memcpy.
   Corrected handling of MOON_AT_SV case.

   Revision 4.1  2003/08/14 15:50:18  vlin
   Simplified by moving common code with GEO_locate_one_scan() out
   into GEO_interp_ECR().  Cleaned up messages.

   Revision 2.5  1999/05/12 22:23:59  kuyper
   Corrected to produce spacecraft position and velocity in ECR coordinates.
  
   Revision 2.4  1999/01/29  21:18:48  kuyper
   Changed to log error messages if GEO_interp_ephemeris returns a WARNING.
   Changed to call GEO_get_sample_time() to load sample_time[].
   Removed max_extrap argument of GEO_interp_ephemeris_attitude().
  
	Revision 1.1  1996/12/10 20:08:22  fshaw
	Initial revision

Requirements:
		PR03-F-3.5.2-1
		PR03-F-3.5.3-1
		PR03-F-3.5.3-2
		PR03-F-3.5.3-3
		PR03-F-3.6.2-1
		PR03-F-3.6.3-1
		PR03-F-3.6.3-2
		PR03-I-1
		PR03-I-2
		
!Team-unique Header:

           This software is developed by the MODIS Science Data Support
           Team for the National Aeronautics and Space Administration,
           Goddard Space Flight Center, under contract NAS5-32373.

References and Credits:   None

!END
*******************************************************************************/

{
  enum {SUN_AT_SD, SUN_AT_EV, MOON_AT_SV, NUM_TIMES};
  static celestial_bodies_struct const empty_cb={
    {	/* frame state structure */
      GEO_DOUBLE_FILLVALUE,	/* time, position, velocity, euler angles */
      {GEO_DOUBLE_FILLVALUE, GEO_DOUBLE_FILLVALUE, GEO_DOUBLE_FILLVALUE},
      {GEO_DOUBLE_FILLVALUE, GEO_DOUBLE_FILLVALUE, GEO_DOUBLE_FILLVALUE},
      {GEO_DOUBLE_FILLVALUE, GEO_DOUBLE_FILLVALUE, GEO_DOUBLE_FILLVALUE},
      {	/* T_inst2ecr */
	{GEO_DOUBLE_FILLVALUE, GEO_DOUBLE_FILLVALUE, GEO_DOUBLE_FILLVALUE},
	{GEO_DOUBLE_FILLVALUE, GEO_DOUBLE_FILLVALUE, GEO_DOUBLE_FILLVALUE},
	{GEO_DOUBLE_FILLVALUE, GEO_DOUBLE_FILLVALUE, GEO_DOUBLE_FILLVALUE}
      },
      BAD_DATA  /* qa_flag */
    },
    /* Sun unit vector. */
    {GEO_DOUBLE_FILLVALUE, GEO_DOUBLE_FILLVALUE, GEO_DOUBLE_FILLVALUE},
    GEO_DOUBLE_FILLVALUE, GEO_DOUBLE_FILLVALUE, /* Solar Zenith, Azimuth */
    /* Moon unit vectory. */
    {GEO_DOUBLE_FILLVALUE, GEO_DOUBLE_FILLVALUE, GEO_DOUBLE_FILLVALUE},
  };
  sc_state_struct sc_state[NUM_TIMES]; /* the SC kinematic state record */
  double sunSD[3] = {0.0};         /* Solar Diffuser frame solar vector */
  double sunsc[3];             /* MODIS instrument frame solar vector */
  double moonecr[3] = {0.0};              /* the ECR frame lunar vector */
  double mooninst[3] = {0.0};    /* MODIS instrument frame lunar vector */
  double T_sc2ecr[NUM_TIMES][3][3] = {0.0}; 
  double T_inst2ecr[NUM_TIMES][3][3] = {0.0}; 
                             /* instrument frame to ECR frame transform */
  PGSt_double sunecr[3] = {0.0};              /* ECR frame solar vector */
  PGSt_double moon_vector[1][3] = {0.0};      /* ECI frame lunar vector */
  PGSt_double sun_vector[2][3] = {0.0};   /* the ECI frame solar vector */
  PGSt_double toff[NUM_TIMES] = {0.0}; /* time offsets from CCSDS_EV_start */
  PGSt_double tveci[NUM_TIMES][6] = {0.0};
  /* ECI position/velocity vector pair*/
  PGSt_double tvecr[NUM_TIMES][6] = {0.0};  
  /* ECR position/velocity vector pair.*/
  PGSt_double positionECR[NUM_TIMES][3], velocityECR[NUM_TIMES][3];
  PGSt_SMF_status PGS_error_code = PGS_S_SUCCESS;
  int ecr, eci, inst, sc, sd;
  PGSt_SMF_status 	retval=PGS_S_SUCCESS;
  char CCSDS_EV_start[TIMECODEASIZE] = "";
  PGSt_double	rpy[3];
  uint32	sample_quality[NUM_TIMES][2];
  char msgbuf[PGS_SMF_MAX_MSGBUF_SIZE] = "";
  static char filefunc[] = __FILE__ ", GEO_solar_and_lunar_vectors";

  if (frame_data == NULL || geo_params == NULL || fill_values==NULL ||
      cb_vectors == NULL || frame_time == NULL)
  {
      sprintf(msgbuf,
	  "frame_data = %p, geo_params = %p, fill_values=%p, cb_vectors = %p"
          "frame_time = %p", (void*)frame_data, (void*)geo_params, 
          (void*)fill_values, (void*)cb_vectors, (void*)frame_time); 
      modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
      return MODIS_E_BAD_INPUT_ARG;
  }

  *cb_vectors = empty_cb;

  /* Spacecraft state at Earth View center */
  /* The Earth View center time is defined as the midpoint
    of the integration of the middle frame of the Earth View 
    scan.  If the Earth View has an even number of frames, 
    the integration midpoint of the lower numbered frame 
    adjacent to the scan center is used. */

  if(frame_data->EV_start == fill_values->EV_start_time)
  {
      modsmf(MODIS_E_BAD_INPUT_ARG, "EV start time is a Fill value.", filefunc);
      return MODIS_E_BAD_INPUT_ARG;
  }

  if (frame_data->EV_frames < 0 || frame_data->EV_frames > MAX_FRAMES)
      toff[SUN_AT_EV] = frame_time[MAX_FRAMES/2];
  else
      toff[SUN_AT_EV] = frame_time[frame_data->EV_frames/2];

  /* Note: we could check whether SD_start or SV_start are fill values.
   * However, that would require splitting up the next function call into two
   * calls. Since the current fill values should produce error condition which
   * should set the corresponding outputs to PGSd_GEO_ERROR_VALUE, it should
   * still have the right effect.
   */
  toff[SUN_AT_SD] = (frame_data->SD_start - frame_data->EV_start) +
      geo_params->geometry_params.t_frame * (double)frame_data->SD_frames / 2.0;
  toff[MOON_AT_SV] = (frame_data->SV_start - frame_data->EV_start) +
      geo_params->geometry_params.t_frame * (double)frame_data->SV_frames / 2.0;

  PGS_error_code = GEO_interp_ECR(frame_data->EV_start, NUM_TIMES, toff,
      geo_params, CCSDS_EV_start, sc_state, T_sc2ecr, T_inst2ecr, positionECR,
      velocityECR, rpy, sample_quality);
  if (PGS_error_code != PGS_S_SUCCESS) {
      sprintf(msgbuf, "GEO_interp_ECR(%.3f)", frame_data->EV_start);
      modsmf(MODIS_E_GEO, msgbuf, filefunc);
      return PGS_error_code == MODIS_E_GEO ? PGS_S_SUCCESS : MODIS_E_GEO;
  }

  if (sc_state[SUN_AT_EV].position[0] < PGSd_GEO_ERROR_VALUE)
  {
      cb_vectors->sc.time = (float64)sc_state[SUN_AT_EV].time;
      for (ecr=0; ecr<3; ecr++)
      {
	  cb_vectors->sc.positionECR[ecr] =
	      (float64)positionECR[SUN_AT_EV][ecr];
	  cb_vectors->sc.velocityECR[ecr] =
	      (float64)velocityECR[SUN_AT_EV][ecr];
	  cb_vectors->sc.eulerAngles[ecr] = 
	      (float64)sc_state[SUN_AT_EV].eulerAngles[ecr];
          for (inst=0; inst<3; inst++)
          cb_vectors->sc.T_inst2ecr[ecr][inst] =
	      (float64)T_inst2ecr[SUN_AT_EV][ecr][inst];
      }
      cb_vectors->sc.qa_flag = GOOD_DATA;
  }

  if (positionECR[SUN_AT_EV][0] < PGSd_GEO_ERROR_VALUE ||
      positionECR[SUN_AT_SD][0] < PGSd_GEO_ERROR_VALUE)
  {
      PGS_error_code = PGS_CBP_Earth_CB_Vector(2, CCSDS_EV_start, toff, 
                       PGSd_SUN, sun_vector); 
      if (PGS_error_code != PGS_S_SUCCESS) {
          sprintf(msgbuf, "PGS_CBP_Earth_CB_Vector(%s, SUN)", CCSDS_EV_start);
          modsmf(MODIS_E_GEO, msgbuf, filefunc);
          retval = MODIS_E_GEO;
      }
      else
      {
	  /* Sun vectors retrieved OK.  Sun Vector from spacecraft at SD center
	   *  time tveci[SUN_AT_SD] =
	   *  sun_vector[SUN_AT_SD] - sc_state[SUN_AT_SD].position:
	   * Adjustment of sun vector from wrt/Earth center to wrt/ spacecraft
	   * sufficiently small to be ignored         
	   */

	  memcpy(tveci[SUN_AT_SD], sun_vector[SUN_AT_SD],
	      sizeof(sun_vector[SUN_AT_SD]));
	  memcpy(tveci[SUN_AT_EV], sun_vector[SUN_AT_EV],
	      sizeof(sun_vector[SUN_AT_EV]));
						
	  PGS_error_code =
	      PGS_CSC_ECItoECR(2, CCSDS_EV_start, toff, tveci, tvecr);
	  if (PGS_error_code != PGS_S_SUCCESS &&
	      PGS_error_code != PGSCSC_W_PREDICTED_UT1)
	  {
	      sprintf(msgbuf, "PGS_CSC_ECItoECR(%s, SUN)", CCSDS_EV_start);
	      modsmf(MODIS_E_GEO, msgbuf, filefunc);
	      retval = MODIS_E_GEO;
	  }
	  else	/* ECR-frame sun vectors computed OK */
	  {
	      if(positionECR[SUN_AT_EV][0] < PGSd_GEO_ERROR_VALUE &&
		  tvecr[SUN_AT_EV][0] < PGSd_GEO_ERROR_VALUE &&
		  cb_vectors->sc.qa_flag == GOOD_DATA)
	      {
		  for(ecr=0; ecr<3; ecr++)
		      sunecr[ecr] = (double)tvecr[SUN_AT_EV][ecr];

		  if (GEO_vec_unit3(sunecr, cb_vectors->sun_unit_vector)
		      != SUCCESS)
		  {
		      modsmf(MODIS_E_GEO, "GEO_vec_unit3(sun_unit_vector)",
			  filefunc);
		      retval = MODIS_E_GEO;
		  }
	      }

	      /* complete computation of sun on solar diffuser */

	      if(positionECR[SUN_AT_SD][0] < PGSd_GEO_ERROR_VALUE &&
		  tvecr[SUN_AT_SD][0] < PGSd_GEO_ERROR_VALUE)
	      {
		  for(ecr=0; ecr<3; ecr++)
		      sunecr[ecr] = (double)tvecr[SUN_AT_SD][ecr];

		  /* transform the sunecr vector in ECR coordinates at the
		   * SD_center_time to sunsc, the vector in the MODIS instrument
		   * frame using the transpose of the T_inst2ecr matrix.
		   */
		  for (sc=0; sc<3; sc++)
		  {
		      sunsc[sc] = 0.0;
		      for (ecr=0; ecr<3; ecr++)
			  sunsc[sc] +=
			      T_sc2ecr[SUN_AT_SD][ecr][sc] * sunecr[ecr];
		  }
		      
		  /* transform the suninst vector in instrument coordinates at
		   * the SD center time to sunSD, the vector in the Solar
		   * Diffuser frame instrument frame using the
		   * geo_params->coord_trans.T_inst2SD matrix:
		   */
		  for (sd=0; sd<3; sd++)
		  {
		      sunSD[sd] = 0.0;
		      for (sc=0; sc<3; sc++)
			  sunSD[sd] += geo_params->coord_trans.T_sc2SD[sd][sc]
			      * sunsc[sc];
		  }
		      
		  /* Note: the sign of the arctangent arguments is used to
		   * determine the quadrant of the angle returned
		   */
		  cb_vectors->SD_sun_zenith = PGS_PI/2.0 - atan2(sunSD[2],
		      sqrt(sunSD[0] * sunSD[0] + sunSD[1] * sunSD[1]));
		  cb_vectors->SD_sun_azimuth = atan2(sunSD[1], sunSD[0]);
	      }
	  }	/* ECR-frame sun vectors computed OK */
      }	/* sun vectors retrieved OK */
  }

  /* Lunar vectors */

  PGS_error_code = PGS_CBP_Earth_CB_Vector(1, CCSDS_EV_start, &toff[MOON_AT_SV],
                   PGSd_MOON, moon_vector);
  if (PGS_error_code != PGS_S_SUCCESS) {
      sprintf(msgbuf, "PGS_CBP_Earth_CB_Vector(%s, MOON)", CCSDS_EV_start);
      modsmf(MODIS_E_GEO, msgbuf, filefunc);
      retval = MODIS_E_GEO;
  }
  else
  {  /*  sc_state retrieved OK  */
    for (eci = 0; eci < 3; eci++)
      tveci[MOON_AT_SV][eci] =
	  moon_vector[0][eci] - sc_state[MOON_AT_SV].position[eci];

    PGS_error_code = PGS_CSC_ECItoECR(1, CCSDS_EV_start, toff+MOON_AT_SV,
	tveci+MOON_AT_SV, tvecr+MOON_AT_SV);
    if (PGS_error_code != PGS_S_SUCCESS &&
	PGS_error_code != PGSCSC_W_PREDICTED_UT1 && (
	tvecr[MOON_AT_SV][0] >= PGSd_GEO_ERROR_VALUE ||
        tvecr[MOON_AT_SV][1] >= PGSd_GEO_ERROR_VALUE ||
        tvecr[MOON_AT_SV][2] >= PGSd_GEO_ERROR_VALUE))
    {
       sprintf(msgbuf, "PGS_CSC_ECItoECR(%s, MOON)", CCSDS_EV_start);
       modsmf(MODIS_E_GEO, msgbuf, filefunc);
       retval = MODIS_E_GEO;
    }
    else
    {
       for(ecr=0; ecr<3; ecr++)
	  moonecr[ecr] = (double)tvecr[MOON_AT_SV][ecr];
       /* transform the moonecr vector in ECR coordinates at the SV center
         * time to mooninst, the vector in the MODIS instrument frame using the
         * transpose of the T_inst2ecr matrix. */

       for(inst=0; inst<3; inst++) {
         mooninst[inst] = 0.0;
         for (ecr=0; ecr<3; ecr++)
           mooninst[inst] += T_inst2ecr[MOON_AT_SV][ecr][inst] * moonecr[ecr];
       }

       if (GEO_vec_unit3(mooninst, cb_vectors->moon_unit_vector) != SUCCESS) {
           modsmf(MODIS_E_GEO, "GEO_vec_unit3(moon_unit_vector)", filefunc);
           retval = MODIS_E_GEO;
       }
    } 
  }	/* sc_state retrieved OK */

  return retval;
}
