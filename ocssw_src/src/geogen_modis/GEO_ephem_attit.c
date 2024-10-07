/* GEO_ephem_attit.c includes geolocation functions GEO_get_ancil_packet_time(), 
 * GEO_interp_ephemeris_attitude() and GEO_prepare_ancil_data().
 *
 * $Log: GEO_ephem_attit.c,v $
 * Revision 6.11  2015/05/20 21:39:44  jkuyper
 * Corrected resolution of Bug 5334.
 *
 * Revision 6.10  2014/08/18 22:36:00  jkuyper
 * Resolved Bug 5535 by changing to delta <= 0.
 *
 * Revision 6.9  2013/07/29 19:39:17  jkuyper
 * Corrected defective subscripting of eulerAngles introduced in Revision 6.3.
 *
 * Revision 6.8  2013/06/24 22:39:54  jkuyper
 * Corrected to zero-initialize sample_quality.
 * Corrected subscript.
 * Corrected to search qualityFlags for good attitude records.
 *
 * Revision 6.7  2013/06/18 20:12:14  jkuyper
 * Reverted attitQuat to a local variable.
 *
 * Revision 6.6  2010/06/29 19:27:02  kuyper
 * GEO_prepare_ancil_data(): Corrected typo, and off-by-one error.
 * GEO_interp_ephemeris_attitude(): made many minor fixes to the code for
 * handling entrained ephemeris and attitude data, without resolving known bugs.
 * I'll have to deliver it as is, for a science test, which won't make use of
 * entrained data.
 *
 * Revision 6.5  2010/06/18 20:37:10  kuyper
 * Corrected to initialize last_qcount.
 * Corrected syntax error in bit-field test.
 *
 * Revision 6.4  2010/06/01 17:45:18  kuyper
 * Changed GEO_prepapre_ancil_data() to return a status code.
 *
 * Revision 6.3  2010/05/14 22:19:16  kuyper
 * Simplified interface of GEO_prepare_ancil_data().
 * Changed to pass ephemeris and attitude quality flags returned by
 *   PGS_EPH_EphemAttit back to caller.
 * Changed to set the quality flags when using entrained ephemeris and attitude
 *   information, and to use those quality flags the same way as when they're
 *   retrieved from PGS_EPH_EphemAttit. As a result, a major rearrangement of
 *   orbit and attitude validation procedures was required.
 * Removed global validation parameters.
 * Changed to use "Safe Mode" rather than "Science Abnormal state" to reject
 *   ancillary data packets.
 * Corrected to set attitQuat by calling PGS_CSC_EulerToQuat when using entrained
 *   attitude data.
 *
 * Revision 6.2  2009/05/31 01:29:07  ltan
 * Minor corrections.
 *
 * Revision 6.1  2009/05/18 22:18:11  ltan
 * Changed macro name to MAX_PADDED. Removed leading dimensions of "array" parameters. Changed MAX_SCAN_SAMPLE to PADDED_SAMPLES.
 *
 * Revision 5.6  2009/05/15 19:00:41  kuyper
 * Corrected to initialize sc_state[].position[0], to resolve Bug 2403.
 *
 * James Kuyper Jr.	James.R.Kuyper@nasa.gov
 *
 * Revision 5.5  2006/12/06 23:02:13  kuyper
 * Corrected macro names.
 *
 * Revision 5.4  2006/11/08 22:24:05  kuyper
 *  Added PM-1 specific quality flags.
 * Changed to use PGSd_PLATFORM_FATAL to filter data.
 * Net effect will be to drop data points collected with missing status words
 * for the attitude files.
 *
 * Revision 5.3  2004/10/27 16:20:27  kuyper
 * Changed to produce a valid, but (negligibly) incorrect quaternion if run
 *   in "MODIS packet" mode.
 *
 * Revision 5.2  2004/10/18 14:51:51  kuyper
 * Added attitQuat to the pointer validity test check list.
 *
 * Revision 5.1  2004/09/20 19:14:35  kuyper
 * Made attitQuat an output parameter.
 *
 * Revision 4.3  2003/10/30 22:46:07  kuyper
 * Corrected re-ordering of attitude angles.
 *
 * Revision 4.2  2003/08/26 20:27:54  kuyper
 * Corrected offsets to be a constant.
 *
 * Revision 4.1  2003/07/17 18:56:32  vlin
 * updated after code walkthrough.
 *
 * Revision 4.0  2003/05/23 18:47:03  vlin
 * Changes to function GEO_interp_ephemeris_attitude:
 *   Take an array of time offsets for input, and to write to an output array.
 *   Treat missing eph/att files as an error condition, while treating missing data
 *   within those files as a normal condition.  Changed to return a status value,
 *   use automatically allocated arrays for storing PGS_EPH_EphemAttit() output
 *   handle numValues==0 as a valid value, allow excess ephemeris/attitude files.
 * Added new function GEO_get_ancil_packet_time.
 *
 * Revision 3.1  2002/06/13 22:48:25  kuyper
 * Removed unnecessary NCSA acknowledgement.
 *
 * Revision 2.6  2000/08/15 18:56:26  kuyper
 * Corrected strcasecmp() to strcmp().
 *
 * Revision 2.5  2000/08/07  14:56:55  fhliang
 * Changes to GEO_prepare_ancil_data():
 *   Changed to use the spacecraft_ID to determine which function should be used
 *   to convert packet header timestamps to TAI times.
 *   Changed to set sc_evolution.spacecraftTag, based upon the value of
 *   spacecraft_ID.
 *   Corrected to not bother loading or ancillary data if it comes from Aqua.
 * Changes to GEO_interp_ephemeris_attitude():
 *   Changed to use sc_evolution.spacecraftTag in call to PGS_EPH_EphemAttit.
 *   Made the message printed if that call fails somewhat more useful.
 *
 * Revision 2.4  1999/02/11  22:27:45  kuyper
 * Changes to GEO_prepare_ancil_data():
 *   Moved in setting of orbit/attit validation globals from
 *   GEO_prepare_l1a_data().
 *   Added max_extrap member to sc_evolution.
 *   Changed to take GEO_param_struct parameter, rather than ancil_scale_struct.
 *   Expanded index names to clarify code.
 *   Moved initialization of num_samples to just before corresponding changes.
 *   Corrected to not use SCI_STATE to validate packets.
 *   Changed checks of ==FAIL to !=SUCCESS.
 * Changes to GEO_interp_ephemeris_attitude():
 *   Removed max_extrap parameter.
 *   Added check of qFlag values.
 *
 * Check the RCS directory at modular:/L1A/INHOUSE/PGE01/MOD_PR03 
 * for earlier revision history 
 */

#include <float.h>
#include "PGS_CSC.h"
#include "smfio.h"
#include "GEO_earth.h"
#include "GEO_input.h"
#include "GEO_product.h"
#include "GEO_util.h"
#include "PGS_MODIS_35251.h"
#ifndef STATIC
#define STATIC static
#endif

double position_abs_limit[2];	/* orbit position absolute limits */
double position_mag_limit[2];	/* orbit position magnitude limits */
double velocity_abs_limit[2];	/* orbit velocity absolute limits */
double velocity_mag_limit[2];	/* orbit velocity magnitude limits */
double ang_mom_limit[2];	/* angular momentum magnitude limits */
double ang_mom_z_limit[2]; /* angular momentum Z component absolute limits */
double orbit_consistency; /* orbit position/velocity consistency limit */
double angvel_abs_limit[2];      /* angular velocity absolute limits */
double angvel_del_limit;		    /* angular velocity delta limits */
double attit_consistency;  /* angle/angular velocity consistency limit */
double attitude_abs_limit[2];    /* attitude angle absolute limits */
double attitude_del_limit;       /* attitude angle delta limits */

typedef struct {
	/* number of validated ephemeris and attitude samples */
	int		num_samples;
	double		max_extrap;
	sc_state_struct sc_state[MAX_SCAN_NUMBER * 2];
        PGSt_tag	spacecraftTag;
} sc_evolution_t;

static sc_evolution_t *sc_evolution = NULL;

PGSt_double GEO_get_ancil_packet_time( PGSt_double      in_time)

/*****************************************************************************
!C

!Description:  
      Determines the time of the last ancillary data packet in sc_evolution, 
      with a timestamp less than or equal to the specified time.

!Input Parameters:
      in_time             Time value that controls the search.

!Output Parameters:       None

Return Value:
      -DBL_MAX            If there are no ancillary data packets, or if the 
                          requested time comes before the earliest one.
      the requested time  If there are some ancillary data packets.

Externally Defined:
      DBL_MAX             "float.h"
      sc_evolution        GEO_ephem_attit.c

Called by:
      GEO_prepare_mirr_data

Routines Called:          None

!Revision History:
      Revision history is located at the beginning of this file.

Requirements:

!Team-unique Header:

      This software is developed by the MODIS Science Data Support
      Team for the National Aeronautics and Space Administration,
      Goddard Space Flight Center, under contract NAS5-32373.

References and Credits:   None

Design Notes:
      Assumes that sc_evolution contains records in ascending order by time.

!END
*****************************************************************************/

{
  int  samp;

  for (samp = sc_evolution->num_samples - 1 ; samp > -1 ; samp--) {
       if (sc_evolution->sc_state[samp].time <= in_time)
           break;
  }

  if (samp < 0)
      return -DBL_MAX;
  else
      return sc_evolution->sc_state[samp].time;
}

/*============================================================================*/

/* Value for SS_CP_MODE field in L1A file that indicates instrument is in
 * safe mode. */
#define SAFE_MODE	2

PGSt_SMF_status GEO_prepare_ancil_data(
	int const			number_of_scans,
	const GEO_param_struct		* params,
	const sc_ancil_struct		sc_ancillary_data[2],
	const uint16			ss_cp_mode[]
	)
/*
!C*****************************************************************************
!Description:   
	Routine in Input group of the Level-1A geolocation software to unpack,
	convert and validate the spacecraft orbit and attitude data from the
	ancillary packets.

!Input Parameters:
	number_of_scans		the number of scans in the granule
	params			Geolocation paramters structure.
	sc_ancillary_data	Spacecraft ancillary data
	ss_cp_mode          	Command processor mode

!Output Parameters:
	None

Return Values:
	MODIS_E_BAD_INPUT_ARG	If any pointer argument is invalid.
	MODIS_E_GEO_ANCIL_DATA	If there are too few valid samples to permit
				    interpolation.
	PGS_S_SUCCESS		Otherwise.

Externally Defined:
	ATT_IDX			"GEO_global_arrays.h"
	EPH_IDX			"GEO_global_arrays.h"
	EPH_LONG_FOLLOW		"GEO_global_arrays.h"
	EPH_LONG_PRECEED	"GEO_global_arrays.h"
	EPH_SAFE_MODE		"GEO_global_arrays.h"
	EPH_SHORT_FOLLOW	"GEO_global_arrays.h"
	EPH_SHORT_PRECEED	"GEO_global_arrays.h"
	EPH_YELLOWHI		"GEO_global_arrays.h"
	EPH_YELLOWLO		"GEO_global_arrays.h"
	MAX_SCAN_NUMBER		"GEO_geo.h"
	MODIS_E_BAD_INPUT_ARG	"PGS_MODIS_35251.h"
	MODIS_E_GEO		"PGS_MODIS_35251.h"
	MIN_TIME_OFFSET		"GEO_geo.h"
	PGS_PI			"GEO_geo.h"
	PGS_S_SUCCESS		"PGS_SMF.h"
	PGSd_EOS_AM1		"PGS_TD.h"
	PGSd_EOS_PM1		"PGS_TD.h"
	PGSd_NO_DATA		"PGS_EPH.h"
	PGSd_PLATFORM_FATAL	"PGS_EPH.h"
	QFL_IDXS		"GEO_global_arrays.h"
	SC_CURRENT		"GEO_input.h"
	SC_PRIOR		"GEO_input.h"
	SS_CP_MODE		"L1a_data.h"
	TAI_FLAG		"GEO_geo.h"

Called by:
	GEO_prepare_l1a_data()

Routines called:
	PGS_TD_EOSAMtoTAI	Converts AM timestamps to TAI
	modsmf			writes error status messages to log


!Revision History:
   See top of file for revision history after merge into GEO_ephem_attit.c

Requirements:
		PR03-F-2.2-1
                PR03-F-2.2-2
                PR03-F-3.2.1-3
		PR03-I-1
		PR03-I-2

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/
{

  /* Factor for conversion of arcseconds to radians */
  static double const arcsec_to_rad = PGS_PI/(180.0*3600.0);
  /* spacecraft attitude rate */
  int scan = 0; /* scan loop index */
  int pack = 0; /* prior/current index */
  int cmb = 0; /* combined scan& prior/current index.	*/
  int samp = 0; /* index in selected packets */
  int eci = 0; /* dimension loop index */

  /* ancillary packet quality*/
  int last_qcount = INT_MAX;		/* Last value for qcount[cmb] */
  double delta;				/* Time difference between packets. */
  int last_att;
  static char filefunc[] = __FILE__ ", GEO_prepare_ancil_data";
  char msgbuf[128];

  typedef double double3_t[3];
  static double3_t *sc_position; /* spacecraft position */
  static double3_t *sc_velocity; /* spacecraft velocity */
  static double3_t *sc_attit; /* spacecraft attitude */
  static double3_t *sc_xyzRotRates;
  
  typedef uint32 ancil_flags_t[MAX_SCAN_NUMBER*2];
  static ancil_flags_t *ancil_flags;
  static int *qcount;	/* Overall packet quality */

  static int firstRun = 1;
  if(firstRun) {
    firstRun = 0;
      
    sc_position = (double3_t*) calloc(MAX_SCAN_NUMBER*2, sizeof(double3_t)); /* spacecraft position */
    sc_velocity = (double3_t*) calloc(MAX_SCAN_NUMBER*2, sizeof(double3_t)); /* spacecraft velocity */
    sc_attit = (double3_t*) calloc(MAX_SCAN_NUMBER*2, sizeof(double3_t)); /* spacecraft attitude */
    sc_xyzRotRates = (double3_t*) calloc(MAX_SCAN_NUMBER*2, sizeof(double3_t));

    ancil_flags = (ancil_flags_t*) calloc(QFL_IDXS, sizeof(ancil_flags_t));
    qcount = (int*) calloc(MAX_SCAN_NUMBER*2, sizeof(int));	/* Overall packet quality */

  }
  bzero(ancil_flags, QFL_IDXS * sizeof(ancil_flags_t));
  bzero(qcount, MAX_SCAN_NUMBER*2 * sizeof(int));	/* Overall packet quality */
  
  if(number_of_scans<0 || number_of_scans>MAX_SCAN_NUMBER || params==NULL ||
    sc_ancillary_data==NULL || ss_cp_mode==NULL)
  {
    sprintf(msgbuf, "number_of_scans:%d, params:%p, ancil_scale_factors:%p "
	   "ss_cp_mode:%p", number_of_scans, (void*)params, 
           (void*)sc_ancillary_data, (void*)ss_cp_mode);
    modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
    return MODIS_E_BAD_INPUT_ARG;
  }

  if(number_of_scans==0)
    return PGS_S_SUCCESS;

  if(!sc_evolution) {
      sc_evolution = calloc(1, sizeof(sc_evolution_t));
  } else {
    memset(sc_evolution, 0, sizeof(sc_evolution_t));
  }
  
  if (strcmp(params->spacecraft_ID, "Aqua") == 0) {
    sc_evolution->spacecraftTag = PGSd_EOS_PM;
    return PGS_S_SUCCESS; /* There is no valid ancillary data available.*/
  }
  else
  {  /* Assume that the spacecraft_ID is "Terra"; it's been validated *
      * elsewhere, there's no need for a redundant validity check.    */
      sc_evolution->spacecraftTag = PGSd_EOS_AM;
  }

  /* First extract data from ancillary data message */
  /* Loop through scans */
  for (scan = 0; scan < number_of_scans; scan++)
  {
    /* Loop for previous and current data */
    for (pack=SC_PRIOR; pack<=SC_CURRENT; pack++)
    {
      double ang_mom[3];	/* Angular momentum pseudo-vector. */
      double pos2_sum=0.0, vel2_sum=0.0, ang_mom2_sum=0.0; /* Cumulators */
      double pos_mag, vel_mag, ang_mag;	/* vector magnitudes */

      cmb = 2*scan+pack;

      for(eci =  0; eci < 3; eci++)
      {
	  /* Validity tests on pre-scaled data */
	  if(sc_ancillary_data[pack].attit_angvel[scan][eci] <
	      params->attit_valid_params.att_valid_range[0] ||
	      params->attit_valid_params.att_valid_range[1] <
	      sc_ancillary_data[pack].attit_angvel[scan][eci])
	  {
	      ancil_flags[ATT_IDX][cmb] |= PGSd_NO_DATA;
	      qcount[cmb] += 3;
	  }

	  /* Apply scale factors. */
	  sc_position[cmb][eci] =
	      (double)sc_ancillary_data[pack].posvel[scan][eci]
	      * params->ancil_params.ancil_scale_factors.S_position;
	  sc_velocity[cmb][eci] =
	      (double)sc_ancillary_data[pack].posvel[scan][eci+3]
	      * params->ancil_params.ancil_scale_factors.S_velocity;
	  sc_attit[cmb][eci] =
	      (double)sc_ancillary_data[pack].attit_angvel[scan][eci]
	      * params->ancil_params.ancil_scale_factors.S_attitude;
	  sc_xyzRotRates[cmb][eci] =
	      (double)sc_ancillary_data[pack].attit_angvel[scan][eci+3]
	      * params->ancil_params.ancil_scale_factors.S_angvel;
	  
	  pos2_sum += sc_position[cmb][eci]* sc_position[cmb][eci];
	  vel2_sum += sc_velocity[cmb][eci]* sc_velocity[cmb][eci];

	  /* Per-component single-packet validity tests. */
	  if(sc_position[cmb][eci] <
	      params->orbit_valid_params.position_abs_limit[0])
	  {
	      ancil_flags[EPH_IDX][cmb] |= EPH_YELLOWLO;
	      qcount[cmb] += 3;
	  }
	  else if( params->orbit_valid_params.position_abs_limit[1] <
	      sc_position[cmb][eci])
	  {
	      ancil_flags[EPH_IDX][cmb] |= EPH_YELLOWHI;
	      qcount[cmb] += 3;
	  }

	  if(sc_velocity[cmb][eci] <
	      params->orbit_valid_params.velocity_abs_limit[0])
	  {
	      ancil_flags[EPH_IDX][cmb] |= EPH_YELLOWLO;
	      qcount[cmb] += 3;
	  }
	  else if( params->orbit_valid_params.velocity_abs_limit[1] <
	      sc_velocity[cmb][eci])
	  {
	      ancil_flags[EPH_IDX][cmb] |= EPH_YELLOWHI;
	      qcount[cmb] += 3;
	  }

	  if(sc_attit[cmb][eci] <
	      params->attit_valid_params.attitude_abs_limit[0])
	  {
	      ancil_flags[ATT_IDX][cmb] |= EPH_YELLOWLO;
	      qcount[cmb] += 3;
	  }
	  else if( params->attit_valid_params.attitude_abs_limit[1] <
	      sc_attit[cmb][eci])
	  {
	      ancil_flags[ATT_IDX][cmb] |= EPH_YELLOWHI;
	      qcount[cmb] += 3;
	  }

	  if(sc_xyzRotRates[cmb][eci] <
	      params->attit_valid_params.angvel_abs_limit[0])
	  {
	      ancil_flags[ATT_IDX][cmb] |= EPH_YELLOWLO;
	      qcount[cmb] += 3;
	  }
	  else if( params->attit_valid_params.angvel_abs_limit[1] <
	      sc_xyzRotRates[cmb][eci])
	  {
	      ancil_flags[ATT_IDX][cmb] |= EPH_YELLOWHI;
	      qcount[cmb] += 3;
	  }

      }
      /* Calculate ang_mom as a cross product. */
      GEO_vec_mul3(sc_position[cmb], sc_velocity[cmb], ang_mom);

      for(eci = 0; eci < 3; eci++)
	  ang_mom2_sum += ang_mom[eci]*ang_mom[eci];
      
      pos_mag = sqrt(pos2_sum);
      vel_mag = sqrt(vel2_sum);
      ang_mag = sqrt(ang_mom2_sum);

      /* Other single packet validity checks. */
      if(pos_mag < params->orbit_valid_params.position_mag_limit[0])
      {
	  ancil_flags[EPH_IDX][cmb] |= EPH_YELLOWLO;
	  qcount[cmb]++;
      }
      else if( params->orbit_valid_params.position_mag_limit[1] < pos_mag)
      {
	  ancil_flags[EPH_IDX][cmb] |= EPH_YELLOWHI;
	  qcount[cmb]++;
      }

      if(vel_mag < params->orbit_valid_params.velocity_mag_limit[0])
      {
	  ancil_flags[EPH_IDX][cmb] |= EPH_YELLOWLO;
	  qcount[cmb]++;
      }
      else if( params->orbit_valid_params.velocity_mag_limit[1] < vel_mag)
      {
	  ancil_flags[EPH_IDX][cmb] |= EPH_YELLOWHI;
	  qcount[cmb]++;
      }

      if(ang_mag < params->orbit_valid_params.ang_mom_limit[0])
      {
	  ancil_flags[EPH_IDX][cmb] |= EPH_YELLOWLO;
	  qcount[cmb]++;
      }
      else if( params->orbit_valid_params.ang_mom_limit[1] < ang_mag)
      {
	  ancil_flags[EPH_IDX][cmb] |= EPH_YELLOWHI;
	  qcount[cmb]++;
      }

      if(ang_mom[2] < params->orbit_valid_params.ang_mom_z_limit[0])
      {
	  ancil_flags[EPH_IDX][cmb] |= EPH_YELLOWLO;
	  qcount[cmb]++;
      }
      else if( params->orbit_valid_params.ang_mom_z_limit[1] < ang_mom[2])
      {
	  ancil_flags[EPH_IDX][cmb] |= EPH_YELLOWHI;
	  qcount[cmb]++;
      }
    }
  }

  sc_evolution->max_extrap = params->max_extrap;

  for (scan = 0; scan < number_of_scans; scan++)
  {
    if(ss_cp_mode[scan] == SAFE_MODE)
	continue;	/* Safe mode packets should not be used. */

    /* Loop for previous and current data */
    for (pack=SC_PRIOR; pack<=SC_CURRENT; pack++)
    {

      PGSt_double time_stamp[MAX_SCAN_NUMBER*2]; /* spacecraft time stamp */
      /* Secondary packet header for Terra */
      PGSt_IO_L0_SecPktHdrEOS_AM AM_packet;
      cmb = 2*scan+pack;

      /* Unpack time stamp */
      memcpy(AM_packet.scTime, sc_ancillary_data[pack].second_header[scan],
             sizeof(AM_packet.scTime));
      /* This memcpy() is done solely to cover the remote possiblity that:
       * the type of scTime is changed to be different from the type of
       * second_header. memcpy() is used rather than a pointer cast, in case
       * scTime is changed to a type with stricter alignment requirements than
       * second_header; an even more unlikely possibilty, and one that is
       * irrelevant so long as both fields are the first ones in their
       * respective structures. With the current typedefs and structure
       * definitions, both fields are arrays of unsigned char.
       */

      /* Convert from spacecraft time to TAI */

      if (PGS_TD_EOSAMtoTAI(AM_packet.scTime, &time_stamp[cmb])
	  != PGS_S_SUCCESS)
      {
	/* call SDP function to report error */
	sprintf(msgbuf, "PGS_TD_EOSAMtoTAI() for scan:%d pack:%d", scan, pack);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);

	continue;	/* Skip packets with invalid time tags. */
      }

      if(sc_evolution->num_samples > 0)
      {
	delta = time_stamp[cmb] -
	  sc_evolution->sc_state[sc_evolution->num_samples-1].time;
	if(delta <= 0)
	{	/* Duplicate packet */
	  if(qcount[cmb] < last_qcount)
	  {
	      sc_evolution->num_samples--;	/* Replace existing packet */
	      if(sc_evolution->num_samples > 0)
		  delta = time_stamp[cmb] -
		    sc_evolution->sc_state[sc_evolution->num_samples-1].time;
	      sc_evolution->sc_state[sc_evolution->num_samples].
		  qualityFlags[ATT_IDX] = 0;
	      sc_evolution->sc_state[sc_evolution->num_samples].
		  qualityFlags[EPH_IDX] = 0;
	  }
	  else
	      continue;	/* Skip packet */
	}

	if(params->max_non_gap < delta)
	{
	    if(params->orbit_valid_params.eph_max_short_gap < delta)
	    {
		sc_evolution->sc_state[sc_evolution->num_samples-1].
		    qualityFlags[EPH_IDX] |= EPH_LONG_FOLLOW;
		sc_evolution->sc_state[sc_evolution->num_samples].
		    qualityFlags[EPH_IDX] |= EPH_LONG_PRECEED;
	    }
	    else{
		sc_evolution->sc_state[sc_evolution->num_samples-1].
		    qualityFlags[EPH_IDX] |= EPH_SHORT_FOLLOW;
		sc_evolution->sc_state[sc_evolution->num_samples].
		    qualityFlags[EPH_IDX] |= EPH_SHORT_PRECEED;
	    }

	    if(params->attit_valid_params.att_max_short_gap < delta)
	    {
		sc_evolution->sc_state[sc_evolution->num_samples-1].
		    qualityFlags[ATT_IDX] |= EPH_LONG_FOLLOW;
		sc_evolution->sc_state[sc_evolution->num_samples].
		    qualityFlags[ATT_IDX] |= EPH_LONG_PRECEED;
	    }
	    else{
		sc_evolution->sc_state[sc_evolution->num_samples-1].
		    qualityFlags[ATT_IDX] |= EPH_SHORT_FOLLOW;
		sc_evolution->sc_state[sc_evolution->num_samples].
		    qualityFlags[ATT_IDX] |= EPH_SHORT_PRECEED;
	    }
	}
      }

      last_qcount = qcount[cmb];
      sc_evolution->sc_state[sc_evolution->num_samples].qualityFlags[ATT_IDX]
	  |= ancil_flags[ATT_IDX][cmb];
      sc_evolution->sc_state[sc_evolution->num_samples].qualityFlags[EPH_IDX]
	  |= ancil_flags[EPH_IDX][cmb];
      
      sc_evolution->sc_state[sc_evolution->num_samples].time = time_stamp[cmb];

      for (eci = 0; eci < 3; eci++)
      {

	  sc_evolution->sc_state[sc_evolution->num_samples].position[eci] =
	      sc_position[cmb][eci];
	  sc_evolution->sc_state[sc_evolution->num_samples].velocity[eci] =
	      sc_velocity[cmb][eci];

	  /* Need to convert attitude from arcseconds to radians */
	  sc_evolution->sc_state[sc_evolution->num_samples].eulerAngles[eci] = 
	      sc_attit[cmb][eci];
	  sc_evolution->sc_state[sc_evolution->num_samples].xyzRotRates[eci] = 
	      sc_xyzRotRates[cmb][eci];
      } /* End of vector loop */

      sc_evolution->num_samples++;
    }
  }

  /* Set unused time samples to flag values to be safe */
  for (samp = sc_evolution->num_samples; samp < MAX_SCAN_NUMBER*2; samp++)
    sc_evolution->sc_state[samp].time = TAI_FLAG;

  if (sc_evolution->num_samples < 2) {

    /* call SDP function to report error */
    sprintf(msgbuf, "%d", sc_evolution->num_samples);
    modsmf(MODIS_E_GEO_ANCIL_DATA, msgbuf, filefunc);
    sc_evolution->num_samples = 0;
    return MODIS_E_GEO_ANCIL_DATA;
  }

  last_att = sc_evolution->sc_state[samp].qualityFlags[ATT_IDX] & PGSd_NO_DATA
      ? -1 : 0;
  for(samp=0; samp < sc_evolution->num_samples-1; samp++)
  {
      double del;	/* Difference between interpolated and actual value.*/

      delta = sc_evolution->sc_state[samp+1].time -
	      sc_evolution->sc_state[samp].time;

      for(eci = 0; eci < 3; eci++)
      {
	  del = fabs(sc_evolution->sc_state[samp+1].position[eci] -
		     sc_evolution->sc_state[samp  ].position[eci] - 
		    (sc_evolution->sc_state[samp+1].velocity[eci] +
	             sc_evolution->sc_state[samp  ].velocity[eci]) * delta/2.0);
	  if(params->orbit_valid_params.orbit_consistency < del)
	  {
	      sc_evolution->sc_state[samp  ].qualityFlags[EPH_IDX]
		  |= EPH_YELLOWHI;
	      sc_evolution->sc_state[samp+1].qualityFlags[EPH_IDX]
		  |= EPH_YELLOWHI;
	  }

      }

      if((sc_evolution->sc_state[samp+1].qualityFlags[ATT_IDX] & PGSd_NO_DATA)
	  == 0)
      {
	  if(last_att > -1)
	  {
	      delta = sc_evolution->sc_state[samp+1].time -
		      sc_evolution->sc_state[last_att].time;
	      for(eci = 0; eci < 3; eci++)
	      {
		  double del_att =
		      fabs(sc_evolution->sc_state[samp+1  ].eulerAngles[eci] -
		           sc_evolution->sc_state[last_att].eulerAngles[eci]);
		  double del_angv =
		      fabs(sc_evolution->sc_state[samp+1  ].xyzRotRates[eci] -
			   sc_evolution->sc_state[last_att].xyzRotRates[eci]);
		  del = fabs(sc_evolution->sc_state[samp + 1].position[eci] -
			     sc_evolution->sc_state[last_att].position[eci] -
			    (sc_evolution->sc_state[samp + 1].velocity[eci] +
			     sc_evolution->sc_state[last_att].velocity[eci])
			    * delta/2.0);
		  if(del_att > 
		      params->attit_valid_params.attitude_del_limit * delta ||
		      del_angv  >
		      params->attit_valid_params.angvel_del_limit * delta ||
		      del >
		      params->attit_valid_params.attit_consistency * delta)
		  {
		      sc_evolution->sc_state[last_att].qualityFlags[ATT_IDX]
			  |= EPH_YELLOWHI;
		      sc_evolution->sc_state[samp + 1].qualityFlags[ATT_IDX]
			  |= EPH_YELLOWHI;
		  }
	      }
	  }
      }
  }

  /* The validity criteria are in arcseconds, so this conversion had to be
   * delayed */
  for(samp=0; samp < sc_evolution->num_samples; samp++)
      for(eci = 0; eci < 3; eci++)
	  sc_evolution->sc_state[samp].eulerAngles[eci] *= arcsec_to_rad;

  return PGS_S_SUCCESS;
}

/*============================================================================*/

PGSt_SMF_status  GEO_interp_ephemeris_attitude(
	PGSt_integer    	numValues,
	char            	asciiUTC[],
	PGSt_double     	target_time,
	const PGSt_double	offsets[],
	GEO_param_struct const	*params,
	sc_state_struct		sc_state[],
	uint32			sample_quality[][QFL_IDXS],
	PGSt_scTagInfo		*scTagInfo
) 
/*******************************************************************************
!C

!Description:   
       Routine in Earth Location group of the Level-1A geolocation software 
       to interpolate the ephemeris and attitude to the sample time.  Spacecraft 
       kinematic state data is retrieved from either the spacecraft ancillary 
       data message in the L1A or through the SDP toolkit Spacecraft and 
       Attitude Data Access Tools.  Data from the spacecraft (SC) ancillary data 
       message are interpolated using a cubic polynomial based on successive
       values and their first derivatives described on pp. 3-1 - 3-3 in the 
       "Theoretical Basis of the SDP Toolkit Geolocation Package for the ECS 
       Project" (445-TP-002-002).  It first checks whether the requested time 
       is within the range of coefficients previously computed.  If not, it 
       finds the bracketing ephemeris and attitude samples and recomputes the 
       interpolation coefficients.  Finally it interpolates the orbit position, 
       orbit velocity attitude angles and body rates to the sample time.

!Input Parameters:
       numValues          The number of values requested
       asciiUTC           UTC reference start time
       target_time        TAI reference start time, should match asciiUTC
       offsets            An array of offsets, in seconds, from asciiUTC
       params		  Geolocation parameters structure.

!Output Parameters:
       sc_state           SC kinematic state records to be filled in for 
                          the specified times
       sample_quality	  Ephemeris/attitude quality flags for each sample

Return Values:
       MODIS_E_GEO                     If any subroutine fails.
       MODIS_E_BAD_INPUT_ARG           If sc_state is a null pointer or
                                       numValues is too large.
       MODIS_E_INSUFFICIENT_PKTS       If there's too few ancillary pacekts.
       MODIS_E_UNKNOWN_PARAMETER       If EA source has an invalid value.
       PGS_S_SUCCESS                   Otherwise

Externally Defined:
       CUBIC                           "GEO_util.h
       EA_SOURCE                       "GEO_product.h"
       EA_SOURCE_SELECT_LUN            "GEO_geo.h"
       L1A_PKT_EA                      "GEO_product.h"
       MODIS_E_GEO                     "PGS_MODIS_35251.h"
       MODIS_E_BAD_INPUT_ARG           "PGS_MODIS_35251.h"
       MODIS_E_INSUFFICIENT_PKTS       "PGS_MODIS_35251.h"
       MODIS_E_SAMP_TIME_OOR           "PGS_MODIS_35251.h"
       MODIS_E_UNKNOWN_PARAMETER       "PGS_MODIS_35251.h"
       MODIS_N_GEO_EXCESS_FILES        "PGS_MODIS_35251.h"
       MAX_PADDED                      "GEO_geo.h"
       PGS_S_SUCCESS                   "PGS_SMF.h"
       PGSd_NO_DATA                    "PGS_EPH.h"
       PGSd_PC_LINE_LENGTH_MAX         "PGS_PC.h"
       PGSd_GEO_ERROR_VALUE            "PGS_TD.h"
       PGSEPH_W_BAD_EPHEM_VALUE        "PGS_EPH_5.h"
       PGSEPH_E_NO_SC_EPHEM_FILE       "PGS_EPH_5.h"
       sc_evolution                    GEO_ephem_attit.c
       SUCCESS                         "GEO_basic.h"
       TOOLKIT_EA                      "GEO_product.h"
 
Routines Called:
       GEO_check_ea_headers	"GEO_earth.h"
       GEO_poly_coef1		"GEO_util.h"
       GEO_poly_fit		"GEO_util.h"
       modsmf			"smfio.h"
       PGS_EPH_EphemAttit	"PGS_EPH.h"
       PGS_PC_GetConfigData	"PGS_PC.h"

Called by:
      GEO_interp_ECR()

!Revision History:
   See top of file for revision history after merge into GEO_ephem_attit.c

!Team-unique Header:

           This software is developed by the MODIS Science Data Support
           Team for the National Aeronautics and Space Administration,
           Goddard Space Flight Center, under contract NAS5-32373.

Requirements:
		PR03-F-3.2.1-1
		PR03-F-3.2.1-2
		PR03-F-3.2.1-3
		PR03-I-1
		PR03-I-2
		PR03-S-1

Design Notes:
       Assumes sc_evolution global structure contains temporally monitonically 
       increasing data.
       Assumes the L1A file does not temporally overlap.

!END
******************************************************************************/
{
    static sc_state_struct const empty_sc={0.0};
    double t_norm = 0.0;	/* normalized sample time */
    static double position_coef[3][CUBIC]={0.0};
    static int eph = -1;
    static int attlo, atthi;
    static char ea_source[PGSd_PC_LINE_LENGTH_MAX]="";
    PGSt_SMF_status PGS_error_code = PGS_S_SUCCESS;
    PGSt_SMF_status retval = PGS_S_SUCCESS;
    PGSt_integer	qualityFlags[MAX_PADDED][2] = {0};
    PGSt_double positionECI[MAX_PADDED][3], eulerAngles[MAX_PADDED][3],
	       velocityECI[MAX_PADDED][3], xyzRotRates[MAX_PADDED][3];
    PGSt_double		attitQuat[MAX_PADDED][4];
    int samp = 0;  /* iteration parameter */
    char msgbuf[PGS_SMF_MAX_MSGBUF_SIZE] = "";
    static char filefunc[] = __FILE__ ", GEO_interp_ephemeris_attitude";

    static const unsigned EPH_CRITERIA=(PGSd_NO_DATA|PGSd_PLATFORM_FATAL);

    /* Check output structure available */
    if (sc_state == NULL || offsets == NULL || params==NULL ||
	sample_quality == NULL || scTagInfo == NULL)
    {
	sprintf(msgbuf, "sc_state = %p, offsets = %p, params=%p, \n"
	    "sample_quality=%p, scTagInfo = %p",
	    (void *)sc_state, (void *)offsets, (void*)params,
	    (void*)sample_quality, (void*)scTagInfo);
      modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
      return MODIS_E_BAD_INPUT_ARG;
    }

    *sc_state = empty_sc;

    /* Check if EA_SOURCE has already been read in to static memory */
    if (ea_source[0]=='\0') {	/* Select source of SC kinematic data */
      PGS_error_code=PGS_PC_GetConfigData(EA_SOURCE_SELECT_LUN, ea_source);
      if (PGS_error_code!=PGS_S_SUCCESS) {
	  sprintf(msgbuf,"PGS_PC_GetConfigData(%ld)",(long)EA_SOURCE_SELECT_LUN);
	  modsmf(MODIS_E_GEO, msgbuf, filefunc);
	  return MODIS_E_GEO;
      }
    }

    if(scTagInfo->eulerAngleOrder[0] < 0)
    {	/* scTagInfo has not yet been filled in. */
	if (PGS_EPH_GetSpacecraftData(sc_evolution->spacecraftTag, "",
	    PGSe_TAG_SEARCH, scTagInfo) != PGS_S_SUCCESS)
	{
	    modsmf(MODIS_E_GEO, "PGS_EPH_GetSpacecraftData()", filefunc);
	    return MODIS_E_GEO;
	}
    }

    if (strcmp(ea_source, TOOLKIT_EA)==0)
    {  /* Use SDP toolkit to retrieve SC state */
	PGS_error_code = GEO_check_ea_headers(target_time, scTagInfo);
	if (PGS_error_code != PGS_S_SUCCESS && 
	    PGS_error_code != MODIS_N_GEO_EXCESS_FILES) {
	    sprintf(msgbuf, "GEO_check_ea_headers(%s)", asciiUTC);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);
	    retval = MODIS_E_GEO;
	}

	if (numValues > MAX_PADDED)
	{
	    sprintf(msgbuf, "numValues: %ld", (long)numValues);
	    modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc); 
	    return MODIS_E_BAD_INPUT_ARG;
	}

	if (numValues == 0)
	    return PGS_S_SUCCESS;

	PGS_error_code = PGS_EPH_EphemAttit(sc_evolution->spacecraftTag,
	    (PGSt_integer)numValues, asciiUTC, (PGSt_double*)offsets, PGS_TRUE,
	    PGS_TRUE, qualityFlags, positionECI, velocityECI, eulerAngles,
	    xyzRotRates, attitQuat);

	if (PGS_error_code != PGS_S_SUCCESS  &&
	    PGS_error_code != PGSEPH_W_BAD_EPHEM_VALUE)
	{
	    if (PGS_error_code != PGSEPH_E_NO_SC_EPHEM_FILE)
	    {	/* Because GEO_check_ea_headers() succeeded, a value of 
		 * PGSEPH_E_NO_SC_EPHEM_FILE woud merely indicate a data gap in 
		 * the files, not an actual operator error of failing to stage the
		 * right files. However, all other error codes are serious:
		 */
		sprintf(msgbuf, "PGS_EPH_EphemAttit(%s)", asciiUTC);
		modsmf(MODIS_E_GEO, msgbuf, filefunc);
		retval = MODIS_E_GEO;
	    }
	    /*  Since the output values from PGS_EPH_EphemAttit are meaningless
	     *  there's no point in copying them over - except for the error code
	     *  that is used later on to determine whether or not the data is
	     *  useful. Resolves bug 5334.
	     */
	    for (samp = 0; samp < numValues; samp++)
		sc_state[samp].position[0] = PGSd_GEO_ERROR_VALUE;
	}
	else
	{	/* Successful retrieval of SC state.  Copy to output records */
	    for (samp = 0; samp < numValues; samp++)
	    {
		int i = 0;

		sc_state[samp].time = target_time + offsets[samp];
		for (i = 0; i < 3; i++)
		{
		    sc_state[samp].position[i] = positionECI[samp][i];
		    sc_state[samp].velocity[i] = velocityECI[samp][i];
		}
		sc_state[samp].eulerAngles[0] = eulerAngles[samp][1];
		sc_state[samp].eulerAngles[1] = eulerAngles[samp][2];
		sc_state[samp].eulerAngles[2] = eulerAngles[samp][0];
		if ( ((unsigned)qualityFlags[samp][ATT_IDX] |
		    (unsigned)qualityFlags[samp][EPH_IDX]) & EPH_CRITERIA )
		    sc_state[samp].position[0] = PGSd_GEO_ERROR_VALUE;
		sample_quality[samp][ATT_IDX] =
		    (uint32)qualityFlags[samp][ATT_IDX];
		sample_quality[samp][EPH_IDX] =
		    (uint32)qualityFlags[samp][EPH_IDX];
	    }
	}
    }
    else if (strcmp(ea_source, L1A_PKT_EA) == 0)
    {
	/* Interpolate spacecraft state from data in L1A ancillary data packets.
	   This is essentially the algorithm from the beta version  */
	if (sc_evolution->num_samples < 2)
	{
	    modsmf(MODIS_E_INSUFFICIENT_PKTS, "", filefunc);
	    retval = MODIS_E_INSUFFICIENT_PKTS;
	}
	else
	{
	    memset(sample_quality, 0, numValues * sizeof *sample_quality);
	    for (samp = 0; samp < numValues; samp++)
	    {
		int eci;	/* ECI coordinates loop counter */
		sc_state[samp].time = target_time + offsets[samp];
		sc_state[samp].position[0] = 0.0;

		/* Check if sample time is outside the period bracketed by the
		 * spacecraft state records used for the previous sample.  */

		if (eph < 0 || eph >= sc_evolution->num_samples-2
		   || sc_state[samp].time < sc_evolution->sc_state[eph].time
		   || sc_state[samp].time >= sc_evolution->sc_state[eph+1].time) 
		{	/* New polynomial parameters for the interpolation
		    must be calculated if in range of available data. */
		    if (sc_state[samp].time >= sc_evolution->sc_state[0].time &&
		       sc_state[samp].time < 
		       sc_evolution->sc_state[sc_evolution->num_samples-1].time)
		    { 
		       /* Find spacecraft records bracketing the sample time */
		       for(eph = 0; sc_state[samp].time >= 
			      sc_evolution->sc_state[eph+1].time; eph++);
		    }
		    else if (fabs(sc_evolution->sc_state[0].time -
			sc_state[samp].time) <= params->max_extrap)
			eph = 0;
		    else if (fabs(sc_state[samp].time - sc_evolution->sc_state
			[sc_evolution->num_samples - 1].time) <
			params->max_extrap)
			eph = sc_evolution->num_samples - 2;
		    else
		    {
			sprintf(msgbuf, "%.6f for samp = %d is more than %.3f "
			    "seconds outside the range from %.6f to %.6f\n",
			    sc_state[samp].time, samp, params->max_extrap,
			    sc_evolution->sc_state[0].time,
			    sc_evolution->sc_state
			    [sc_evolution->num_samples-1].time);
			modsmf(MODIS_W_SAMP_TIME_OOR, msgbuf, filefunc);
			sc_state[samp].position[0] = PGSd_GEO_ERROR_VALUE;
			sample_quality[samp][EPH_IDX] |= PGSd_NO_DATA; /* ?? */
			eph = -1;
			continue;
		    }
		    /* get coefficients for the cubic polynomial to interpolate
		     * the spacecraft's position */
		    for (eci = 0; retval==SUCCESS && eci < 3; ++eci) 
		       retval = GEO_poly_coef1(
			   sc_evolution->sc_state[eph  ].position[eci],
			   sc_evolution->sc_state[eph+1].position[eci],
			   sc_evolution->sc_state[eph  ].velocity[eci],
			   sc_evolution->sc_state[eph+1].velocity[eci],
			   sc_evolution->sc_state[eph  ].time,
			   sc_evolution->sc_state[eph+1].time,
			   position_coef[eci]);
		}   /* end if sample time is outside the period */

		if (sc_state[samp].position[0] < PGSd_GEO_ERROR_VALUE)
		{	/* normalize target sample time to a fraction of the
			 * period between SC state records */
		    int efat;	/* Ephemeris/attitude loop counter. */

		    t_norm = (sc_state[samp].time -
			sc_evolution->sc_state[eph].time) /
			(sc_evolution->sc_state[eph+1].time -
			 sc_evolution->sc_state[eph  ].time);
		    for (eci = 0; retval==SUCCESS && eci < 3; eci++)
		    {
			retval = GEO_poly_fit(position_coef[eci], t_norm, CUBIC-1,
			    &sc_state[samp].position[eci]);
			sc_state[samp].velocity[eci] = (position_coef[eci][1] +
			    2*position_coef[eci][2]*t_norm +
			    3*position_coef[eci][3]*t_norm*t_norm) /
			    (sc_evolution->sc_state[eph+1].time -
			     sc_evolution->sc_state[eph  ].time) ;
		    } 

		    if (eci != 3)
		    {
			 modsmf(MODIS_E_GEO, "GEO_poly_fit()", filefunc);
			 retval = MODIS_E_GEO;
			 sc_state[samp].position[0] = PGSd_GEO_ERROR_VALUE;
		    }

		    if(sc_state[samp].time <= sc_evolution->sc_state[eph].time)
			sample_quality[samp][EPH_IDX] =
			    sc_evolution->sc_state[eph  ].qualityFlags[EPH_IDX];
		    else if(sc_evolution->sc_state[eph+1].time
			<= sc_state[samp].time)
			sample_quality[samp][EPH_IDX] =
			    sc_evolution->sc_state[eph+1].qualityFlags[EPH_IDX];
		    else
		    {
			sample_quality[samp][EPH_IDX] =
			    sc_evolution->sc_state[eph  ].qualityFlags[EPH_IDX]|
			    sc_evolution->sc_state[eph+1].qualityFlags[EPH_IDX];
			if(sc_evolution->sc_state[eph].time +
			    params->max_non_gap <
			    sc_evolution->sc_state[eph+1].time)
			    sample_quality[samp][EPH_IDX] |= EPH_REPAIRED;
		    }

		    /* Interpolate attitude data */
		    /* Find previous bracketing good attitude samples. */
		    for(attlo = eph; attlo > -1 &&
			sc_evolution->sc_state[attlo].qualityFlags[ATT_IDX] &
			PGSd_NO_DATA; attlo--);
		    for(atthi = eph + 1; atthi < sc_evolution->num_samples &&
			sc_evolution->sc_state[atthi].qualityFlags[ATT_IDX] &
			PGSd_NO_DATA; atthi++);
		    if(attlo < 0 || attlo >= sc_evolution->num_samples)
		    {
			sample_quality[samp][ATT_IDX] |= PGSd_NO_DATA;
			sc_state[samp].position[0] = PGSd_GEO_ERROR_VALUE;
		    }
		    else
		    {
			if(attlo + 1 < atthi)
			{	/* t_norm needs to be recalculated. */
			    t_norm = (sc_state[samp].time -
				sc_evolution->sc_state[attlo].time) /
				(sc_evolution->sc_state[atthi].time -
				 sc_evolution->sc_state[attlo].time);
			}

			/* Set quality flags. */
			if(sc_state[samp].time <=
			    sc_evolution->sc_state[attlo].time)
			    sample_quality[samp][ATT_IDX] = sc_evolution->
				sc_state[attlo].qualityFlags[ATT_IDX];
			else if(sc_state[samp].time >=
			    sc_evolution->sc_state[atthi].time)
			    sample_quality[samp][ATT_IDX] = sc_evolution->
				sc_state[atthi].qualityFlags[ATT_IDX];
			else
			{
			    sample_quality[samp][ATT_IDX] = sc_evolution->
				sc_state[attlo].qualityFlags[ATT_IDX] |
				sc_evolution->sc_state[atthi].
				qualityFlags[ATT_IDX] | PGSd_INTERPOLATED_POINT;
			    if(sc_evolution->sc_state[attlo].time +
				params->max_non_gap <
				sc_evolution->sc_state[atthi].time)
				sample_quality[samp][ATT_IDX] |= EPH_REPAIRED;
			}

			for (eci = 0; retval==SUCCESS && eci < 3; eci++)
			{
			    sc_state[samp].eulerAngles[eci] =
				sc_evolution->sc_state[attlo].eulerAngles[eci] *
				(1-t_norm) +
				sc_evolution->sc_state[atthi].eulerAngles[eci] *
				t_norm;
			}
		    }

		    /* Adjustments to make flag values internally consistent. */
		    for(efat = 0; efat < QFL_IDXS; efat++)
		    {
			if(sample_quality[samp][efat] &&
			    sample_quality[samp][efat] != PGSd_NO_DATA)
			{
			    if(sample_quality[samp][efat] & PGSd_NO_DATA-1)
				sample_quality[samp][efat] |= EPH_DATA;
			    sample_quality[samp][efat] |= EPH_OVERALL;
			}
		    }
		}
	    }     /*  end for samp = 0 to numValues - 1  */
	}         /*  end else */
    }
    else
    {      /* ea_source is neither TOOLKIT_EA nor L1A_PKT_EA  */
      modsmf(MODIS_E_UNKNOWN_PARAMETER, ea_source, filefunc);
      retval = MODIS_E_UNKNOWN_PARAMETER;
    }

    return retval;
}
