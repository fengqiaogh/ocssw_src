#ifndef GEO_OUTPUT_H
#define GEO_OUTPUT_H
#include "GEO_parameters.h"
#include "GEO_main_func.h"
#include "smfio.h"
#include "hdfi.h"
#include "PGS_PC.h"
#include "mapi.h"

/*
!C-INC*************************************************************************

!Description:	the prototypes for the output functions in the Level-1A 
	geolocation software

!Input Parameters: N/A

!Output Parameters: N/A

!Revision History:
 * $Log: GEO_output.h,v $
 * Revision 6.6  2011/02/18 21:58:48  kuyper
 * In order to resolve feature-request Bug 3446, changed interfaces to pass
 *   landsea mask around directly, rather than through swath_elem, and added
 *   waterpresent array.
 * Corrected declarations of parameters that point at input data.
 *
 * Revision 6.5  2010/06/18 20:09:13  kuyper
 * Changed GEO_locate_one_scan(), GEO_write_scan_metadata(), GEO_write_granule_metadata(), GEO_write_geospecific_metadata(), to return a
 *   status code.
 * Corrected parameters that were pointers to array by removing their leading
 *   'const' qualifiers.
 * Helped resolve Bug 2470 by removing the relevant structure from all places
 *   where it was used.
 * Helped resolve Bug 2472 by adding rpy and frame_quality parameters to GEO_write_one_scan() and GEO_write_scan_metadata().
 *
 * Revision 6.4  2009/05/27 14:30:29  xgeng
 * Removed sc_ev_frame_state and changed frame_time's data type
 * for GEO_write_one_scan.
 *
 * Revision 6.3  2009/04/28 19:56:37  kuyper
 * Corrected two minor typos.
 *
 * Revision 6.2  2009/04/17 17:15:02  kuyper
 * Removed arrays in polar coordinates from swath_elem struct, as well as gflags
 *   and hires_offsets.
 * Changed interfaces of GEO_landsea_mask() and GEO_write_scan_data() to allow
 *   removal of those arrays.
 *
 * Revision 6.1  2009/03/23 19:44:15  kuyper
 * Changed MAX_DETECTORS and MAX_SCAN_SAMPLE to DETECTORS_1KM and MAX_FRAMES.
 * Changed "sample" to "frame".
 * Changed pixel_flags to frame_flags.
 * Added hires_offsets member to swath_elem_struct.
 * Changed return type and argument list for GEO_write_one_scan().
 * CHanged return type for GEO_create_swath() and GEO_derived_products().
 *
 * James Kuyper Jr. 	James.R.Kuyper@nasa.gov
 *
 * Revision 5.4  2008/09/21 18:51:33  kuyper
 * Added sc_tag parameter to GEO_write_input_metadata().
 *
 * Revision 5.3  2004/10/08 19:33:13  vlin
 * GEO_write_granule_metadata and GEO_update_L1A_metadata updated
 *
 * Revision 5.2  2004/08/31 19:48:33  vlin
 * Functions GEO_write_scan_metadata, GEO_write_geospecific_metadata,
 *           GEO_update_L1A_metadata, GEO_write_one_scan, and
 *           GEO_locate_one_scan updated.
 *
 * Revision 5.1  2004/08/17 16:45:46  vlin
 * GEO_initialize_product() updated.
 *
 * Revision 4.6  2004/08/16 21:13:39  vlin
 * GEO_DOUBLE_FILLVALUE updated
 *
 * Revision 4.5  2004/04/09 22:10:26  kuyper
 * Added FILL_INT8.
 *
 * Revision 4.4  2003/08/20 14:49:51  kuyper
 * Reinstated redundant num_detectors parameter for GEO_initialize_product();
 *   the interface is frozen for now.
 *
 * Revision 4.3  2003/08/12 21:46:08  kuyper
 * Added entry for GEO_write_input_metadata().
 *
 * Revision 4.2  2003/08/12 14:58:07  vlin
 * function GEO_derived_products() updated.
 *
 * Revision 4.1  2002/12/14 18:51:11  kuyper
 * Changed GEO_initialize_product, GEO_derived_products,
 *   GEO_update_l1a_metadata, and GEO_write_ECS_metadata to return status codes.
 * Changed argument lists for GEO_initialize_product and GEO_derived_products.

!Team-unique Header:

	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.


!END***************************************************************************
*/

/* Macro definitions	*/
#define		GEO_DOUBLE_FILLVALUE	9.9692099683868690e+36	
/* If char is unsigned, FILL_BYTE has the wrong value. */
#define FILL_INT8 ((int8)-127)

/* Structure definitions */
typedef struct
{
 	int	num_frames;
	int	num_detectors;
	uint8	lat_qaflag;
	uint8	lon_qaflag;
	uint8	height_qaflag;
	uint8	sensorazimuth_qaflag;
	uint8	sensorzenith_qaflag;
	uint8	solarazimuth_qaflag;
	uint8	solarzenith_qaflag;
	uint8	range_qaflag;
	uint8	land_seamask[DETECTORS_1KM*MAX_FRAMES];
		/* EOS Land/Sea mask format */
	uint8	land_seamask_qaflag;
} swath_elem_struct;

typedef struct {
	char	header[80];
	float64	polar_motion[3];
} utcpole_metadata_struct;

/* function prototype */

PGSt_SMF_status GEO_initialize_product(	/* Initialize file for HDF output */
	int				const /* number_of_scans */,
	fill_values_struct const	* /* fill_values */,
	MODFILE				* const /* geo_file */,
	int32				const /* swfid */,
	GEO_param_struct const		* /* GEO_param */
	);

PGSt_SMF_status GEO_derived_products(
	/* Compute derived geolocation parameters */
	int const,
	int const,
	double[][MAX_FRAMES][3],
	double[][3],
	unsigned char[][MAX_FRAMES],
	double[][MAX_FRAMES][3],
	char[28],
	PGSt_double[],
	double[][MAX_FRAMES][2],
	double[][MAX_FRAMES][3]
);

PGSt_SMF_status GEO_get_utcpole_metadata( /* Fill in utcpole_metadata */
	const ECS_metadata_struct  *ECS_metadata,
 	utcpole_metadata_struct    *utcpole_metadata
	);

PGSt_SMF_status GEO_get_version_metadata( /* Loading most of version_metadata*/
	version_metadata_struct *version_metadata
	);

PGSt_SMF_status GEO_landsea_mask(
	int			/* num_samples */,
	int			/* num_detectors */,
	double			/* terrain_sample_position*/[][MAX_PADDED][3],
	uint8			/* sample_flags */[][MAX_PADDED],
	uint8			* /*land_seamask_qaflag*/,
	uint8			/* sample_landsea */[][MAX_PADDED]
);

PGSt_SMF_status GEO_locate_one_scan(
	/* Perform geolocation for one scan */
	GEO_param_struct const	* const geo_params,
	l1a_data_struct		* const l1a_data,
	int			const scan_number,
	qa_metadata_struct	* const qa_metadata,
	GEO_bcoord_struct	* const bounding_coords,
	MODFILE			* const geo_file
);

PGSt_SMF_status GEO_update_L1A_metadata(
	MODFILE             * const l1a_file,
	GEO_GRing_struct    * const GRing_points,
	EPH_metadata_struct * const EPH_metadata,
	GEO_bcoord_struct   * const bounding_coords
	);

PGSt_SMF_status GEO_write_ECS_metadata(
	MODFILE                 * const geo_file,
	ECS_metadata_struct     * const ECS_metadata,
	EPH_metadata_struct     * const EPH_metadata,
	GEO_bcoord_struct       * const bounding_coords,
	GEO_GRing_struct        * const GRing_points,
 	pointer_metadata_struct * const pointer_metadata,
	qa_metadata_struct      * const qa_metadata,
	char                    * const sci_state,  
	char                    * const sci_abnorm  
	);

PGSt_SMF_status GEO_write_one_scan(
	/* Write data to product for one scan */
	PGSt_double const	/* frame_time */[],
	int const               /* scan_number */,
	l1a_data_struct const	* const /* l1a_data */,
	GEO_param_struct const	* const /* geo_params */,
	double                  /* frame_to_sensor */[][MAX_FRAMES][3],
	double                  /* frame_solar_angles */[][MAX_FRAMES][2],
	double			/* terrain_frame_position */[][MAX_FRAMES][3],
	uint8			/* frame_flags */[][MAX_FRAMES],
	int8			/* hires_offsets */
					[][DETECTORS_QKM][SAMPLES_QKM],
	PGSt_double const	/* rpy */[],
	uint32			/* frame_quality */[][MAX_FRAMES],
	uint8                   /* frame_landsea */[][MAX_FRAMES],
	uint8                   /* frame_waterpresent */[][MAX_FRAMES],
	uint8                   /* land_seamask_qaflag */,
	MODFILE          	* const /* geo_file */
);

PGSt_SMF_status GEO_write_scan_data(           /* Write geolocation scan data */
	MODFILE			* const /* geo_file */,
	int			const /* scan_number */,
	swath_elem_struct const	* /* swath_elem */,
	GEO_param_struct const	* /* GEO_param */,
	double			/* terrain_frame_position */[][MAX_FRAMES][3],
	double			/* frame_to_sensor */[][MAX_FRAMES][3],
	double			/* frame_solar_angles */[][MAX_FRAMES][2],
	uint8			/* frame_flags */[][MAX_FRAMES],
	uint8			/* frame_landsea */[][MAX_FRAMES],
	uint8			/* frame_waterpresent */[][MAX_FRAMES],
	int8			/* hires_offsets */
					[3][DETECTORS_QKM][SAMPLES_QKM]
);

PGSt_SMF_status GEO_write_scan_metadata( /* Write geolocation scan metadata */
	int				const scan_number,
	l1a_data_struct const		* const l1a_data,
	celestial_bodies_struct const	* const cb_vectors,
	MODFILE				* const geo_file,
	PGSt_double const		rpy[],
	uint32				frame_quality[][MAX_FRAMES]
);

PGSt_SMF_status GEO_write_granule_metadata(
	/* Write geolocation granule metadata */
	MODFILE			* const geo_file,
	MODFILE			* const l1a_file,
	GEO_param_struct const	* const geo_parameter,
	GEO_bcoord_struct const	* const bounding_coords,
	int			const version,
	l1a_data_struct		* const l1a_data,
	qa_metadata_struct	* const qa_metadata
);

int GEO_write_parameters(          /* Write geolocation parameters */
	MODFILE * const		/* geo_file */,
	GEO_param_struct const	* /* parameter */
	);
 
PGSt_SMF_status GEO_write_geospecific_metadata(
	MODFILE				* const geo_file,
	l1a_metadata_struct const	* const granule_metadata,
	int				const number_of_scans,
	GEO_param_struct const		* const geo_parameter,
	qa_metadata_struct const	* const qa_metadata,
	utcpole_metadata_struct const	* const utcpole_metadata
);


int GEO_write_input_metadata(
	MODFILE			* const geo_file,
	int			const l1a_version,
	PGSt_integer		const param_version,
	pointer_metadata_struct	* const pointer_metadata,
	PGSt_tag		sc_tag
	);

int GEO_get_bounding_coords(
	double	terrain_frame_position[DETECTORS_1KM][MAX_FRAMES][3],
	uint8	frame_flags[DETECTORS_1KM][MAX_FRAMES],
	int	const num_detectors,
	int	const num_frames,
	GEO_bcoord_struct * const bounding_coords
);

PGSt_SMF_status GEO_create_swath(
	int	const number_of_scans,
	int	const num_detectors,
	int32	const swfid
);

#endif

