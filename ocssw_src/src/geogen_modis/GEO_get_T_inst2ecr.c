/* 
 * $Log: GEO_get_T_inst2ecr.c,v $
 * Revision 6.7  2017/03/27 17:01:11  jkuyper
 * Corrected to treat a return of PGSCSC_W_BAD_TRANSFORM_VALUE from
 *   PGS_CSC_ECItoECR() as an error.
 *
 * Revision 6.6  2013/07/24 22:20:58  jkuyper
 * Corrected initialization of sol_sc.
 *
 * Revision 6.5  2013/06/18 20:18:47  jkuyper
 * Resolved Bug 2473 by using T_sc2ecr rather than attitQuat to convert sol_eci
 *   to sol_sc. This rendered attitQuat unnecessary, and required substantial
 *   rearrangement of the code.
 *
 * Revision 6.4  2010/06/29 18:47:03  kuyper
 * Changed floating point equality tests into inequalities.
 *
 * Revision 6.3  2010/03/31 20:20:37  kuyper
 * Helped resolve Bug 1969 by making rpy an output argument.
 * Helped resolve Bug 2473 by making eulerAngleOrder an input argument.
 *
 * James Kuyper	james.kuyper@sigmaspace.com
 *
 * Revision 6.2  2009/05/15 16:52:15  xgeng
 * Minor change.
 * 
 * Xu Geng ( xu.geng@saic.com ) 
 * Revision 6.1  2009/05/11 15:28:30  xgeng
 * Removed meaningless leading dimensions of "array" parameters.
 * Changed macro name to MAX_PADDED.
 * The code updates match the PDL change made in revision 6.1 and 6.2.
 *
 * Revision 5.5  2005/01/06 21:58:56  mash
 * Corrected errors in algorithm for applying solar elevation correction.
 *
 * Revision 5.4  2004/11/02 17:38:34  kuyper
 * Added T_sc2ecr to list of outputs, and initialized with fill values.
 * Corrected error in handling of H_eci.
 *
 * Revision 5.3  2004/10/25 14:02:46  kuyper
 * Corrected some typos.
 * Corrected calculation and use of sol_idx.
 *
 * Revision 5.2  2004/10/18 16:57:31  kuyper
 * Added sol_elev_cor and attitQuat to list of input pointer parameters to be
 *   checked.
 *
 * Revision 5.1  2004/10/13 21:36:29  kuyper
 * Changed to apply an adjustment to T_inst2sc based upon the solar "elevation"
 *    angle.
 *
 * Revision 4.7  2003/10/07 21:01:06  kuyper
 * Corrected to set and check for PGSd_GEO_ERROR_VALUE.
 *
 * Revision 4.6  2003/08/28 15:59:18  kuyper
 * Corrected prolog
 *
 * Revision 4.5  2003/08/26 20:31:45  kuyper
 * Corrected offsets to be a pointer to const.
 *
 * Revision 4.4  2003/08/22 22:22:15  kuyper
 * Changed T_sc2ecr into an externally supplied array, to be filled in.
 *
 * Revision 4.3  2003/08/22 15:56:46  kuyper
 * Corrected calculation of sum, to match PDL.
 *
 * Revision 4.2  2003/07/30 21:07:08  vlin
 * updated after code walkthrough.
 *
 * Revision 4.1  2003/07/01 18:37:22  vlin
 * Changed to work from an array of offsets from an ASCII UTC time,
 * rather than from a single TAI time extracted from the sc_state
 * Moved creation of T_sc2orb after call to PGS_CSC_ECItoECR(),
 * allowing that array to be only 3x3.  Re-wrote comments and messages.
 *
 * Revision 3.2  2002/05/29 21:07:37  kuyper
 * Removed test for warning status from UTCtoTAI call.
 * Corrected rpy_times subscripts.
 * Changed to handle even "impossible" cases safely.
 * Corrected range of loop indices.
 *
 * Revision 1.2  1996/07/18 16:28:23  kuyper
 * Included self-checking header file.
 * Changed extern declarations of T_inst2ecr, ecr_sc_position and
 * ecr_sc_velocity to definitions. Replaced all other extern declarations
 * with corresponding header files.
 * Added required header files.
 * Fixed pointer level problem with GEO_mat_mul call.
 */

#include "GEO_earth.h"
#include "GEO_geo.h"
#include "GEO_global_arrays.h"
#include "PGS_CBP.h"
#include "PGS_CSC.h"
#include "PGS_MODIS_35251.h"
#include "smfio.h"

/* Transformation matrix from instrument to spacecraft 
   frame.                                             */
static double T_inst2sc[3][3]; 


PGSt_SMF_status GEO_set_T_inst2sc(
	const internal_coord_trans_struct	* const coord_trans,
	const ECS_metadata_struct		* const ECS_metadata
){
/*******************************************************************************
!C
!Description:   Interpolates T_inst2sc to the time specified in the ECS
	metadata for the beginning of a granule.

!Input Parameters:
	coord_trans		Contains the base T_inst2sc, and time series
				data for the interpolation of roll/pitch/yaw
				corrections to that matrix.
	ECS_metadata		Contains the RangeBeginningDate and Time
				metadata.

!Output Parameters:
	None.

Return Values:
	MODIS_E_BAD_INPUT_ARG	If either argument is NULL
	MODIS_E_GEO		If any subroutine failed.
	PGS_S_SUCCESS		Otherwise

Externally Defined:
	MODIS_E_BAD_INPUT_ARG	PGS_MODIS_35251.h
	MODIS_E_GEO		PGS_MODIS_35251.h
	PGS_S_SUCCESS		PGS_SMF.h
	TIMECODEASIZE		smfio.h

Called by:
	GEO_locate_one_granule()

Routines Called:
	modsmf				"smfio.h"
	PGS_TD_UTCtoTAI			"PGS_TD.h"

!Revision History:
See top of file.
  
Requirements:
See "Implementation of Time dependent Matrix T_inst2sc by Geolocation Parameter
File" by Masahiro Nishihama, 2002-05-10. A copy is is in the SDF for this
routine. The requirements process has died through lack of customer support,
and this will therefore not be being converted into a formal requirement.
  
!Team-unique Header:
  	This software is developed by the MODIS Science Data Support
  	Team for the National Aeronautics and Space Administration,
  	Goddard Space Flight Center, under contract NAS5-32373.
  
Design Notes:
	T_inst2sc is a 3x3 double precision global with file scope, shared with
	GEO_get_T_inst2ecr().

!END
*******************************************************************************/

    double		adjust[3][3];
    double		rpy[3]={0.0};
    double		cosine[3], sine[3];
    double		fraction;
    PGSt_double		TAI_time;
    char		asciiUTC[TIMECODEASIZE];
    char		msgbuf[256];
    char		filefunc[] = __FILE__ ", GEO_set_T_inst2sc";
    int			i, j, k, axis;
    int			first, last;

    if(coord_trans == NULL || ECS_metadata == NULL)
    {
	sprintf(msgbuf, "coord_trans:%p ECS_metadata:%p",
	    (void *)coord_trans, (void *)ECS_metadata);
	modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);

	return MODIS_E_BAD_INPUT_ARG;
    }

    if(coord_trans->rpy_count<1)
    {
	memcpy(T_inst2sc, coord_trans->T_inst2sc, sizeof(T_inst2sc));

	return PGS_S_SUCCESS;
    }

    sprintf(asciiUTC, "%sT%s", ECS_metadata->rangebeginningdate,
	ECS_metadata->rangebeginningtime);

    if(PGS_TD_UTCtoTAI(asciiUTC, &TAI_time)!= PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "PGS_TD_UTCtoTAI(\"%s\")", asciiUTC);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);

	return MODIS_E_GEO;
    }

    for(first=-1; first < coord_trans->rpy_count-1 &&
	coord_trans->rpy_times[first+1] < TAI_time; first++);
    for(last = coord_trans->rpy_count;
	last > 0 && coord_trans->rpy_times[last-1] >= TAI_time; last--);

    if(last<=0)
	memcpy(rpy, coord_trans->rpy_inst2sc[0], sizeof(rpy));
    else if(first >= coord_trans->rpy_count-1)
	memcpy(rpy, coord_trans->rpy_inst2sc[ coord_trans->rpy_count-1],
	    sizeof(rpy));
    else
    {
	/* NOT protected against division by 0: we trust the parameter file.*/
	fraction = (TAI_time-coord_trans->rpy_times[first])/
	    (coord_trans->rpy_times[last] - coord_trans->rpy_times[first]);
	for(axis = 0; axis<3; axis++)
	    rpy[axis] = fraction*coord_trans->rpy_inst2sc[last][axis] +
		(1.0-fraction)*coord_trans->rpy_inst2sc[first][axis];
    }

    for(axis = 0; axis<3; axis++)
    {
	cosine[axis] = cos(rpy[axis]);
	sine[axis] = sin(rpy[axis]);
    }

    adjust[0][0] = cosine[1]*cosine[2]-sine[0]*sine[1]*sine[2];
    adjust[0][1] = -sine[2]*cosine[0];
    adjust[0][2] = cosine[2]*sine[1] + sine[2]*sine[0]*cosine[1];
    adjust[1][0] = sine[2]*cosine[1] + cosine[2]*sine[0]*sine[1];
    adjust[1][1] = cosine[2]*cosine[0];
    adjust[1][2] = sine[2]*sine[1] - cosine[2]*sine[0]*cosine[1];
    adjust[2][0] = -cosine[0]*sine[1];
    adjust[2][1] = sine[0];
    adjust[2][2] = cosine[0]*cosine[1];

    for(i=0; i<3; i++)
	for(j=0; j<3; j++)
	{
	    T_inst2sc[i][j] = 0.0;
	    for(k=0; k<3; k++)
		T_inst2sc[i][j] += adjust[i][k]*coord_trans->T_inst2sc[k][j];
	}
  
    return PGS_S_SUCCESS;
}

/*===========================================================================*/

PGSt_SMF_status GEO_get_T_inst2ecr(
        PGSt_integer            numValues,
        char                    asciiutc[],
        const PGSt_double	offsets[],
        sc_state_struct const   sc_state[],
	double			sol_elev_cor[][3],
	PGSt_integer const	eulerAngleOrder[],
        double                  T_sc2ecr[][3][3],
        double                  T_inst2ecr[][3][3],
        double                  ecr_position[][3],
        double                  ecr_velocity[][3],
	PGSt_double		rpy[]
)

/******************************************************************************
!C

!Description:  Routine in the instrument model of the Level-1A geolocation
               software to calculate transformation matrix from instrument
               to ECR coordinates.  The algorithm is based on the
               Version 3.0 Geolocation ATBD (section 3.1.4.2).

!Input Parameters:
        numValues       The number of time values for which the matrix is to be
                        converted.
        asciiUTC        The base time for the conversion
        offsets         An array of numValues offsets in seconds from the base
                        time.
        sc_state        An array of numValues spacecraft kinematic state
                        records for the specified times. 
	sol_elev_cor	Solar elevation-based lookup table for corrections to
			T_inst2sc.
	eulerAngleOrder	Specification of the order in which euler angles are
			applied.

!Output Parameters:
	T_sc2ecr	An array of numValues tranformation matrices from the
			spacecraft to ECR coordinates for the spacecraft state.
        T_inst2ecr      An array of numValues transformation matrices from the
                        instrument to ECR coordinates for the spacecraft state.
        ecr_position    An array of numValues spacecraft positions from
                        sc_state converted to ECR coordinates.
                        Optional output (may be null).
        ecr_velocity    An array of numValues spacecraft velocities from
                        sc_state converted to ECR coordinates.
                        Optional output (may be null).
	rpy		Array containing the roll, pitch, and yaw corrections
			based upon the current solar "elevation" angle.

Return value:
        MODIS_E_BAD_INPUT_ARG   If any mandatory input pointer is null
        MODIS_E_GEO             If any subroutine fails.
        PGS_S_SUCCESS           Otherwise

Externally defined:
        MODIS_E_BAD_INPUT_ARG   "PGS_MODIS_35251.h"
        MODIS_E_GEO             "PGS_MODIS_35251.h"
	NUM_SOL_ELEV		"GEO_parameters.h"
        MAX_PADDED              "GEO_geo.h"
	PGS_PI			"GEO_geo.h"
        PGS_S_SUCCESS           "PGS_SMF.h"
        PGS_TRUE                "PGS_SMF.h"
	PGSd_SUN		"PGS_CBP.h"
	PGSd_GEO_ERROR_VALUE	"PGS_CSC.h"
        PITCH                   "GEO_geo.h"
	POSITION		"GEO_geo.h"
        ROLL                    "GEO_geo.h"
        T_inst2sc               "GEO_parameters.h"
        VELOCITY                "GEO_geo.h"
        YAW                     "GEO_geo.h"

Called by:
        GEO_interp_ECR
 
Routines called:
        PGS_CSC_crossProduct            "PGS_CSC.h"
        PGS_CSC_ECItoECR                "PGS_CSC.h"
	PGS_CSC_EulerToQuat		"PGS_CSC.h"
	PGS_CSC_Norm			"PGS_CSC.h"
	PGS_CSC_quatRotate		"PGS_CSC.h"
        PGS_SMF_TestWarningLevel        "PGS_SMF.h"
        PGS_TD_UTCtoTAI                 "PGS_TD.h"
        modsmf                          "smfio.h"

!Revision History:
See top of file

!Team-unique Header:
   This software is developed by the MODIS Science Data Support
   Team for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

Requirements:
		PR03-F-3.2.2-1

!END
****************************************************************************/

{
  PGSt_double     tveci[MAX_PADDED*2][6] = {0.0}; 
  PGSt_double     tvecr[MAX_PADDED*2][6] = {0.0}; 
                              /* ECR vector array for transformation routine */
  PGSt_double     toff[MAX_PADDED*2] = {0.0};    
                                      /* time offset array for toolkit calls */
  double T_inst2sc_adj[3][3] = {0.0}; /* adjusted T_inst2sc matrix. */
  int	i, val;  /* iteration parameters */ 
  PGSt_double sol_eci[6]={0.0};
  int sc, ecr, inst, orb;	/* Loop counters. */

  char	msgbuf[PGS_SMF_MAX_MSGBUF_SIZE] = "";
  PGSt_SMF_status      status;
  char  filefunc[] = __FILE__ ", GEO_get_T_inst2ecr";

  /*************************************************************
  calculation from here
  *************************************************************/

  if (numValues < 0 || numValues > MAX_PADDED || asciiutc == NULL ||
      offsets == NULL || sc_state == NULL || sol_elev_cor == NULL ||
      eulerAngleOrder == NULL || T_sc2ecr == NULL || T_inst2ecr == NULL)
  { 
      sprintf(msgbuf, "numValues = %d, asciiutc = %p, offsets = %p, \n"
	  "sc_state = %p, sol_elev_cor=%p, eulerAngleOrder=%p, T_sc2ecr=%p "
	  "T_inst2ecr = %p", numValues, (void*)asciiutc, (void*)offsets,
	  (void*)sc_state, (void*)sol_elev_cor, (void*)eulerAngleOrder,
	  (void*)T_sc2ecr, (void*)T_inst2ecr);
      modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
      return MODIS_E_BAD_INPUT_ARG;
  }
 
  if (numValues == 0)
      return PGS_S_SUCCESS;
  
  if (ecr_position!=NULL)
     for (val=0; val<numValues; val++)
	 for (ecr=0; ecr<3; ecr++)
	      ecr_position[val][ecr] = PGSd_GEO_ERROR_VALUE;

  if (ecr_velocity!=NULL)
     for (val=0; val<numValues; val++)
	 for (ecr=0; ecr<3; ecr++)
	      ecr_velocity[val][ecr] = PGSd_GEO_ERROR_VALUE;

  for (val=0; val<numValues; val++)
      for (ecr=0; ecr<3; ecr++)
      {
	  for (inst=0; inst<3; inst++) 
	  {
	      T_sc2ecr[val][ecr][inst] = PGSd_GEO_ERROR_VALUE;
	      T_inst2ecr[val][ecr][inst] = PGSd_GEO_ERROR_VALUE;
	  }
      }

  for (val=0; val<numValues; val++)
  {	/* construct the spacecraft to orbit coord. transformation matrix
	   T_sc2orb for the target sample. Get trigonometry terms from euler
	   angles */ 
      PGSt_double     H_eci[3]   = {0.0};
      /* ECI vector parallel to SC orbital angular momentum in the ECI frame */

      PGS_CSC_crossProduct((PGSt_double*)sc_state[val].position,
	  (PGSt_double*)sc_state[val].velocity, H_eci);
      memcpy(&tveci[2*val][POSITION], sc_state[val].position, 
             sizeof(sc_state[val].position));
      memcpy(&tveci[2*val][VELOCITY], sc_state[val].velocity, 
             sizeof(sc_state[val].velocity));
      memcpy(&tveci[2*val+1][POSITION], H_eci, sizeof(H_eci)); 
      toff[2*val] = offsets[val];
      toff[2*val+1] = offsets[val];

      if(val==numValues/2 && PGS_CBP_Earth_CB_Vector(1, asciiutc,
	  (PGSt_double*)offsets+numValues/2, PGSd_SUN,
	  (PGSt_double(*)[3])&sol_eci) !=PGS_S_SUCCESS)
      {
	  sprintf(msgbuf, "PGS_CBP_Earth_CB_Vector(%s)", asciiutc);
	  modsmf(MODIS_E_GEO, msgbuf, filefunc);
	  return MODIS_E_GEO;
      }
  }

  status = PGS_CSC_ECItoECR(2*numValues, asciiutc, toff, tveci, tvecr);
  if (status != PGS_S_SUCCESS && status != PGSCSC_W_PREDICTED_UT1){
      sprintf(msgbuf, "PGS_CSC_ECItoECR(%s)", asciiutc);
      modsmf(MODIS_E_GEO, msgbuf, filefunc);
      return MODIS_E_GEO;
  }

  for (val=0; val<numValues; val++)
  {
      double T_orb2ecr[3][3] = {0.0};  /* orbit to ECR transformation matrix */
      double T_sc2orb[3][3];	/* Spacecraft to orbit transformation matrix */
      PGSt_double b1[3]; /* Cross product of the following pair: */
      PGSt_double b2[3]; /* negative of angular momentum vector direction */
      PGSt_double b3[3]; /* nadir vector direction */
      double b2_length, b3_length;
      double sin_roll,  cos_roll;   /* sin, cos of the roll  */
      double sin_pitch, cos_pitch;  /* sin, cos of the pitch */
      double sin_yaw,   cos_yaw;    /* sin, cos of the yaw   */

      if(tveci[2*val][0] >= PGSd_GEO_ERROR_VALUE || 
	  tvecr[2*val][0] >= PGSd_GEO_ERROR_VALUE)
	  continue;

      if (ecr_position != NULL)
          for (ecr = 0; ecr < 3; ecr++)
              ecr_position[val][ecr] = tvecr[2*val][POSITION+ecr];
      if (ecr_velocity != NULL)
          for (ecr = 0; ecr < 3; ecr++)
              ecr_velocity[val][ecr] = tvecr[2*val][VELOCITY+ecr];
      b3_length = PGS_CSC_Norm(tvecr[2*val]+POSITION);
      b2_length = PGS_CSC_Norm(tvecr[2*val+1]+POSITION);
      for (ecr =0; ecr < 3; ecr++) {
          b3[ecr] = -tvecr[2*val][POSITION+ecr]/b3_length;
          b2[ecr] = -tvecr[2*val+1][POSITION+ecr]/b2_length;
      }
      PGS_CSC_crossProduct(b2, b3, b1);

      /* construct the orbit to ECR transformation matrix T_orb2ecr with the b1,
         b2, b3 orbital coordinate unit vectors in the ecr coordinate frame.  */
      for (ecr = 0; ecr < 3; ecr++) {
           T_orb2ecr[ecr][0] = (double)b1[ecr];
           T_orb2ecr[ecr][1] = (double)b2[ecr];
           T_orb2ecr[ecr][2] = (double)b3[ecr];
      }

      /* construct the SC to orbit coord. transformation matrix T_sc2ecr */
      sin_roll  = sin(sc_state[val].eulerAngles[ROLL]);
      cos_roll  = cos(sc_state[val].eulerAngles[ROLL]);
      sin_pitch = sin(sc_state[val].eulerAngles[PITCH]);
      cos_pitch = cos(sc_state[val].eulerAngles[PITCH]);
      sin_yaw   = sin(sc_state[val].eulerAngles[YAW]);
      cos_yaw   = cos(sc_state[val].eulerAngles[YAW]);

      /* get trigonometry terms from euler angles */
      T_sc2orb[0][0] = cos_yaw * cos_pitch - sin_yaw * sin_roll * sin_pitch;
      T_sc2orb[0][1] = -cos_roll * sin_yaw;
      T_sc2orb[0][2] = cos_yaw * sin_pitch + sin_yaw * sin_roll * cos_pitch;
      T_sc2orb[1][0] = sin_yaw * cos_pitch + cos_yaw * sin_roll * sin_pitch;
      T_sc2orb[1][1] = cos_yaw * cos_roll;
      T_sc2orb[1][2] = sin_yaw * sin_pitch - cos_yaw * sin_roll * cos_pitch;
      T_sc2orb[2][0] = -cos_roll * sin_pitch;
      T_sc2orb[2][1] = sin_roll;
      T_sc2orb[2][2] = cos_roll * cos_pitch;

      /* construct the composite transformation matrix T_inst2ecr. Construct
       * the SC frame to ECR coordinate frame transform matrix, T_sc2ecr, from
       * the matrix multiplication of the T_sc2orb and T_orb2ecr transform
       * matrices.
       */
      for (ecr=0; ecr<3; ecr++)
	  for (sc=0; sc<3; sc++)
	  {
	      T_sc2ecr[val][ecr][sc] = 0.0;
	      for(orb=0; orb<3; orb++)
		  T_sc2ecr[val][ecr][sc] +=
		      T_orb2ecr[ecr][orb] * T_sc2orb[orb][sc];
	  }

      if(val==numValues/2)
      {
	  PGSt_double sol_ecr[6];
	  PGSt_double sol_sc[3]={0.0};
	  int sc, ecr;
	  double sol_elev, extra;
	  int sol_idx;
	  PGSt_double  adjust_quat[4];

	  status = PGS_CSC_ECItoECR(1, asciiutc, (PGSt_double*)offsets+val,
	      &sol_eci, &sol_ecr);
	  if (status != PGS_S_SUCCESS && status != PGSCSC_W_PREDICTED_UT1)
	  {
	      sprintf(msgbuf, "PGS_CSC_ECItoECR(%s) new", asciiutc);
	      modsmf(MODIS_E_GEO, msgbuf, filefunc);
	      return MODIS_E_GEO;
	  }
	  else for(sc=0; sc<3; sc++)
	      for(ecr=0; ecr<3; ecr++)
		  sol_sc[sc] += T_sc2ecr[val][ecr][sc]*sol_ecr[ecr];

	  sol_elev = atan2(-sol_sc[2], sol_sc[0]);
	  sol_idx =
	      (int)floor((sol_elev+PGS_PI)*(double)NUM_SOL_ELEV*0.5/PGS_PI);
	  extra = (sol_elev + PGS_PI)/(2.0*PGS_PI/(double)NUM_SOL_ELEV)
	      - (double)sol_idx; 
	  for(i=0; i<3; i++)	/* interpolate to sol_elev. */
	      rpy[(i+1)%3] = sol_elev_cor[sol_idx%NUM_SOL_ELEV][i]*(1.0-extra)
		  + extra*sol_elev_cor[(sol_idx+1)%NUM_SOL_ELEV][i];
	  if(PGS_CSC_EulerToQuat(rpy, (PGSt_integer*)eulerAngleOrder,
	      adjust_quat) !=PGS_S_SUCCESS)
	  {
	      modsmf(MODIS_E_GEO,"PGS_CSC_EulerToQuat(adjustment)" , filefunc);
	      return MODIS_E_GEO;
	  }

	  for(sc=0; sc<3; sc++)
	  {
	      for(i = 0; i < 3; i++) /*  Pulling column vectors out */
		  b1[i] = T_inst2sc[i][sc];
	      if(PGS_CSC_quatRotate(adjust_quat, b1, b2) != PGS_S_SUCCESS)
	      {
		  modsmf(MODIS_E_GEO,"PGS_CSC_quatRotate(adjustment)" ,
		      filefunc);
		  return MODIS_E_GEO;
	      }
	      for(i =0; i < 3; i++)
		  T_inst2sc_adj[i][sc] = b2[i];
	  }
      }
  }

  /* construct the instrument frame to ECR coordinate frame transform matrix
     for the target spacecraft state, T_inst2ecr, from the matrix
     multiplication of the T_sc2ecr and T_inst2sc_adj transform matrices. */
  for (val=0; val<numValues; val++)
      for (ecr=0; ecr<3; ecr++)
	  for (inst=0; inst<3; inst++) 
	  {
	      T_inst2ecr[val][ecr][inst] = 0.0;
	      for (sc=0; sc<3; sc++)
		   T_inst2ecr[val][ecr][inst] +=
		       T_sc2ecr[val][ecr][sc] * T_inst2sc_adj[sc][inst];
	  }

  return PGS_S_SUCCESS;
}
