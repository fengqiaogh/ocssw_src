#include "GEO_earth.h"
#include "PGS_MODIS_35251.h"

PGSt_SMF_status GEO_interp_ECR(
	PGSt_double const	base_time,
	PGSt_integer const	numValues,
	const PGSt_double	offsets[],
	GEO_param_struct const	*params,
	char              	asciiUTC[],
	sc_state_struct   	sc_state[],
	double            	T_sc2ecr[][3][3],
	double            	T_inst2ecr[][3][3],
	double            	positionECR[][3],
	double            	velocityECR[][3],
	PGSt_double		rpy[],
	uint32			sample_quality[][2]
){
/******************************************************************************
!C

!Description:   
     Routine in the Earth Location group of the Level-1A geolocation
     software to interpolate the ephemeris, attitude, and orientation of the
     satellite in ECR coordinates at a specified list of offsets from a
     specified base time.

!Input Parameters:
     base_time    TAI reference start time
     numValues    The number of values requested
     offsets      An array of numValues offsets, in seconds, from asciiUTC
     params	  Geolocation parameters struct

!Output Parameters:
	asciiUTC	The UTC time string corresponding to base_time
	sc_state	An array of numValues spacecraft kinematic state
			records to be filled in for the specified times
	T_sc2ecr	An array of numValues transformation matrices from the
			spacecraft to ECR coordinates for the spacecraft state.
	T_inst2ecr	An array of numValues transformation matrices from the
			instrument to ECR coordinates for the spacecraft state.
	positionECR	An array of numValues spacecraft positions from
			sc_state converted to ECR coordinates. Optional output
			(may be set to NULL).
	velocityECR	An array of numValues spacecraft velocities from
			sc_state converted to ECR coordinates. Optional output
			(may be set to NULL).
	rpy		An array of roll, pitch, and yaw thermal corrections to
			be filled in for the specified times.
	sample_quality	An array of numValues ephemeris and attitude quality
			flags to be filled in for the specified times.
  
Return Values:
     MODIS_E_GEO                 If any subroutine fails.
     MODIS_E_GEO_BAD_INPUT_ARG   If any pointer argument is null.
     PGS_S_SUCCESS               Otherwise

Externally Defined:
     MODIS_E_GEO                 "PGS_MODIS_35251.h"
     MODIS_E_GEO_BAD_INPUT_ARG   "PGS_MODIS_35251.h"
     PGS_S_SUCCESS               "PGS_SMF.h"
     PGS_TRUE                    "PGS_SMF.h"
  
Called by:
     GEO_locate_one_scan
     GEO_solar_and_lunar_vectors
  
Routines Called:
     GEO_get_T_inst2ecr              "GEO_earth.h"
     GEO_interp_ephemeris_attitude   "GEO_earth.h"
     PGS_SMF_TestWarningLevel        "PGS_SMF.h"
     PGS_TD_TAItoUTC                 "PGS_TD.h"

!Revision History:
  $Log: GEO_interp_ECR.c,v $
  Revision 6.3  2013/06/17 16:49:23  jkuyper
  Dropped attitQuat, which is no longer needed.

  Revision 6.2  2010/05/19 22:06:36  kuyper
  Helped resolve Bug 1969 by: Adding rpy as an output parameter filled in by call
    to GEO_get_T_inst2ecr().
  Helped resolve Bug 2742: by replacing max_extrap and sol_elev_cor with params,
    adding sample_quality as an output parameters, and passing params and
    sample_quality to GEO_interp_ephemeris_attitude().
  Helped resolve Bug 2473 by defining scTagInfo, whose eulerAngleOrder member is
    filled in by calling GEO_interp_ephemeris_attitude(), and then passed to
    GEO_get_T_inst2ecr.

  James Kuyper Jr.		James.R.Kuyper@NASA.gov

  Revision 6.1  2009/05/06 19:10:12  xgeng
  Code updated to match the PDL change made in revision 6.1 and 6.2

  Revision 5.2  2004/11/02 01:15:35  kuyper
  Added T_sc2ecr to the list of outputs.

  Revision 5.1  2004/09/23 15:08:09  kuyper
  Added attitQuat and sol_elev_cor

  Revision 4.3  2003/11/17 23:18:32  vlin
  data type for base_time is changed to float

  Revision 4.2  2003/08/22 22:04:27  kuyper
  Added T_sc2ecr as a parameter to be passed to GEO_get_T_inst2ecr().

  Revision 4.1  2003/08/01 20:54:11  vlin
  vlin@saicmodis.com


Requirements:                        None

!Team-unique Header:

   This software is developed by the MODIS Science Data Support
   Team for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

Design Note:                       None

!END
******************************************************************************/

  PGSt_SMF_status status = PGS_S_SUCCESS;
  PGSt_scTagInfo	scTagInfo;
  char msgbuf[128] = "";
  char filefunc[] = __FILE__ ", GEO_interp_ECR";

  if (offsets == NULL || asciiUTC == NULL ||
      sc_state == NULL || T_inst2ecr == NULL) {
      sprintf(msgbuf,"offsets: %p, asciiUTC: %s, sc_state: %p, T_inst2ecr: %p",
             (void *)offsets,asciiUTC,(void *)sc_state,(void *)T_inst2ecr);
      modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
      return MODIS_E_BAD_INPUT_ARG;
  } 

  scTagInfo.eulerAngleOrder[0] = -1;	/* Mark as uninitialized */
  status = PGS_TD_TAItoUTC(base_time, asciiUTC);
  if (status != PGS_S_SUCCESS && PGS_SMF_TestWarningLevel(status) != PGS_TRUE){
      sprintf(msgbuf, "PGS_TD_TAItoUTC(%.2f)", base_time);
      modsmf(MODIS_E_GEO, msgbuf, filefunc);
      return MODIS_E_GEO; 
  }

  if (GEO_interp_ephemeris_attitude(numValues, asciiUTC, base_time, offsets,
      params, sc_state, sample_quality, &scTagInfo) != PGS_S_SUCCESS)
  {
      sprintf(msgbuf, "GEO_interp_ephemeris_attitude(%s)", asciiUTC);
      modsmf(MODIS_E_GEO, msgbuf, filefunc); 
      return MODIS_E_GEO; 
  }

  if (GEO_get_T_inst2ecr(numValues, asciiUTC, offsets, sc_state,
      (double (*)[3])params->sol_elev_cor, scTagInfo.eulerAngleOrder, T_sc2ecr,
      T_inst2ecr, positionECR, velocityECR, rpy) != PGS_S_SUCCESS)
  {  
      sprintf(msgbuf, "GEO_get_T_inst2ecr(%s)", asciiUTC);
      modsmf(MODIS_E_GEO, msgbuf, filefunc);
      return MODIS_E_GEO;
  }
  return PGS_S_SUCCESS;
}

