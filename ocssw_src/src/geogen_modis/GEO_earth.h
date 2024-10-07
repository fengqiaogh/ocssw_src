/*
!C-INC*************************************************************************
!Description:   the prototype of the modules which performs the earth location
	        for the Level-1A geolocation software

!Input Parameters: N/A

!Output Parameters: N/A

!Revision History:
$Log: GEO_earth.h,v $
Revision 6.17  2013/06/18 19:56:01  jkuyper
Removed attitQuat from function interfaces.

Revision 6.16  2013/05/23 18:40:08  jkuyper
Change l1a_data argument of GEO_earth_location() to be a pointer to
  non-const, to allow initialization of mirr_impulse_enc.

Revision 6.15  2011/02/14 21:15:14  kuyper
Changed GEO_aggregate to work on landsea mask, and to fill in waterpresent
  array.

Revision 6.14  2010/08/12 22:13:55  kuyper
changed GEO_initialize_DEM to return a status code.

Revision 6.13  2010/05/27 20:38:34  kuyper
Corrected declaration for GEO_get_ephatt_inputs().

Revision 6.12  2010/05/27 17:45:20  kuyper
Corrected dimensions of frame_quality parameter for GEO_aggregate().

Revision 6.11  2010/05/27 15:20:07  kuyper
Changed input pointers for GEO_solar_and_lunar_vectors() into pointers to
  const.

Revision 6.10  2010/05/19 22:05:54  kuyper
Simplified interfaces by passing a GEO_parameters_struct, rather than it's
  members.

Revision 6.9  2010/04/09 20:31:24  kuyper
Removed obsolete argument from GEO_check_ea_headers().

Revision 6.8  2010/04/02 21:30:23  kuyper
Dropped unused argument of GEO_aggregate().
Corrected several typos.

Revision 6.7  2010/03/31 15:28:36  kuyper
Helped resolve Bug 1969 by adding rpy as an argument to two functions.
Helped resolve Bug 2471 by adding declaration for GEO_get_ephatt_inputs().
Helped resolve Bug 2472 by adding sample_quality, frame_quality, and
  qualityFlags arguments to several functions.
Helped resolve Bug 2473 by adding scTagInfo and eulerAngelOrder arguments.
James Kuyper	james.kuyper@sigmaspace.com

Revision 6.6  2009/05/30 17:23:45  kuyper
Corrected parameter list for GEO_earth_location().

Revision 6.5  2009/05/22 16:17:07  xgeng
Corrected prototype for GEO_solar_and_lunar_vectors.

Revision 6.4  2009/04/28 19:59:51  kuyper
Corrected order of parameters for GEO_aggregate().

Revision 6.3  2009/03/27 21:38:50  kuyper
Removed pointless leading dimensions of "array" parameters.

Revision 6.2  2009/03/20 21:00:15  kuyper
Corrected prototype for GEO_solar_and_lunar_vectors.

Revision 6.1  2009/03/13 14:23:45  kuyper
Replaced MAX_SCAN_SAMPLE with PADDED_SAMPLES.
Replaced MAX_DETECTORS with DETECTORS_QKM.
Changed parameter lists and return types for GEO_ellip_position(),
  GEO_earth_location(), and GEO_terrain_correct().
Added GEO_aggregate() and GEO_hires().

Revision 5.2  2005/03/16 21:37:03  kuyper
Changed header guard macro name to avoid reserved name space.

Revision 5.1  2004/10/14 17:47:55  vlin
Function prototypes updated for collection 5.

Revision 4.5  2003/08/26 21:26:50  kuyper
Changed offsets parameters to be pointers to const.

Revision 4.4  2003/08/22 22:03:25  kuyper
Added T_sc2ecr arguments to two functions.

Revision 4.3  2003/08/11 20:53:54  vlin
GEO_interp_ECR() added.

Revision 4.2  2003/07/23 20:16:03  vlin
several functions updated.

Revision 4.1  2003/03/13 21:40:34  kuyper
Updated two functions to take multiple time values.
Added GEO_check_ea_headers.
James Kuyper Jr. kuyper@saicmodis.com


!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END***************************************************************************
*/

#ifndef GEO_EARTH_H
#define GEO_EARTH_H
#include "PGS_EPH.h"
#include "smfio.h"
#include "GEO_parameters.h"
#include "mapi.h"
#include "GEO_geo.h"
#include "GEO_global_arrays.h"

#define EARTH_MODEL "WGS84"
#define BAD_GEOID -9999
/*********************************************************
function prototypes
*********************************************************/

PGSt_SMF_status GEO_aggregate(
        int32   	/* EV_frames */,
        uint16  	/* N_samp */,
	double		/* hires_scale */,
        unsigned char 	/* sample_flags */[][MAX_PADDED],
        double  	/* ecr_sample_position */[][MAX_PADDED][3],
	double		/* ecr_sc_sample_position */[][3],
        double  	/* terrain_sample_position */[][MAX_PADDED][3],
        uint32		/* sample_quality */[][2],
	uint8		/* sample_landsea */[][MAX_PADDED],
        double  	/* ecr_frame_position */[][MAX_FRAMES][3],
        double  	/* terrain_frame_position */[][MAX_FRAMES][3],
        uint8		/* frame_flags */[][MAX_FRAMES],
	double		/* ecr_sc_frame_position */[][3],
	int8		/* hires_offsets */[][DETECTORS_QKM][SAMPLES_QKM],
        uint32	   	/* frame_quality */[][MAX_FRAMES],
	uint8           /* frame_landsea */[][MAX_FRAMES],
	uint8           /* frame_waterpresent */[][MAX_FRAMES]
);

PGSt_SMF_status GEO_interp_ECR(
	PGSt_double		const base_time,
	PGSt_integer		const numValues,
	const PGSt_double	offsets[],
	GEO_param_struct const	*params,
	char			asciiUTC[],
	sc_state_struct		sc_state[],
	double			T_sc2ecr[][3][3],
	double			T_inst2ecr[][3][3],
	double			positionECR[][3],
	double			velocityECR[][3],
	PGSt_double		rpy[],
        uint32			sample_quality[][2]
);

PGSt_SMF_status GEO_interp_ephemeris_attitude(
	PGSt_integer		numValues,
	char			asciiUTC[],
	PGSt_double		target_time,
	const PGSt_double	offsets[],
	GEO_param_struct const	*params,
	sc_state_struct		sc_state[],
	uint32			qualityFlags[][2],
	PGSt_scTagInfo		*scTagInfo
);

void *GEO_DEMalloc(size_t);

int GEO_close_DEM(void);

PGSt_SMF_status GEO_ellip_position(
	int const scan_number,
	int const sample_number,
	int const num_detectors,
	double u_inst[][3],
	double ecr_sc_sample_position[][3],
	double ecr_sc_velocity[][3],
	double T_inst2ecr[][3][3],
	double ecr_sample_position[][MAX_PADDED][3],
	double ellip_sample_position[][MAX_PADDED][3],
	unsigned char sample_flags[][MAX_PADDED],
	double sample_view_vec[][MAX_PADDED][3]
);

PGSt_SMF_status GEO_terrain_correct(
	int const sample_number,
	int const num_detectors,
	double sample_view_vec[][MAX_PADDED][3],
	double ecr_sample_position[][MAX_PADDED][3],
	double ellip_sample_position[][MAX_PADDED][3],
	unsigned char sample_flags[][MAX_PADDED],
	double terrain_sample_position[][MAX_PADDED][3]
);

PGSt_SMF_status GEO_hires(
        uint16  N_samp,
	int	padded_samples,
	double	hires_scale,
        double  terrain_sample_position[][MAX_PADDED][3],
        double  ecr_sample_position[][MAX_PADDED][3],
        double  ecr_frame_position[][MAX_FRAMES][3],
	uint8	sample_flags[][MAX_PADDED],
	uint8	frame_flags[][MAX_FRAMES],
	int8	hires_offsets[][DETECTORS_QKM][SAMPLES_QKM]
);

PGSt_SMF_status GEO_initialize_DEM(
	void
);

int GEO_latlon2height(
	double const lat,
	double const lon,
	double * const h);

int GEO_read_DEM(
	double lat,
	double lon,
	int * const hgtmin,
	int * const hgtmax);

PGSt_SMF_status GEO_earth_location(
	int const		/* scan_number */,
	int const		/* sample_number */,
        GEO_param_struct const	* /* geo_params */,
	double			/* sample_time */,
	l1a_data_struct		* /* l1a_data */,
	double			/* ecr_sc_sample_position */ [][3],
	double			/* ecr_sc_velocity */ [][3],
	double			/* T_inst2ecr */ [][3][3],
	unsigned char		/* sample_flags */ [][MAX_PADDED],
	double			/* ecr_sample_position */ [][MAX_PADDED][3],
	double			/* terrain_sample_position */ [][MAX_PADDED][3]
);

PGSt_SMF_status GEO_get_ephatt_inputs(
	PGSt_PC_Logical		file_logical,
	PGSt_integer		file_version,
	char			universal_references[][PGSd_UR_FIELD_SIZE]
);

PGSt_SMF_status GEO_get_T_inst2ecr(
        PGSt_integer            numValues,
        char                    asciiutc[],
        const PGSt_double	offsets[],
        sc_state_struct const   sc_state[],
	double                  sol_elev_cor[][3],
	PGSt_integer const	eulerAngleOrder[],
        double                  T_sc2ecr[][3][3],
        double                  T_inst2ecr[][3][3],
        double                  ecr_position[][3],
        double                  ecr_velocity[][3],
	PGSt_double	  	rpy[]
);

PGSt_SMF_status GEO_set_T_inst2sc(
	const internal_coord_trans_struct	* const,
	const ECS_metadata_struct		* const);

PGSt_SMF_status GEO_solar_and_lunar_vectors(
	PGSt_double const		frame_time[],
	frame_data_struct const		* const frame_data,
	GEO_param_struct const		* const geo_params,
	fill_values_struct const	* const fill_values,
	celestial_bodies_struct		* const cb_vectors
);

PGSt_SMF_status GEO_get_GRing_points(
	GEO_GRing_struct	* const GRing_points);

int GEO_get_geoid(
        double latitude,
        double longitude);

PGSt_double	GEO_get_ancil_packet_time (
	PGSt_double in_time);

PGSt_SMF_status GEO_check_ea_headers (
	PGSt_double	in_time,
	PGSt_scTagInfo	*scTagInfo
);
#endif
