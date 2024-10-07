/* file: GEO_input.h */

/*
!C-INC*************************************************************************

!Description:	the header file for the input functions in the Level-1A
	        geolocation software

!Input Parameters: N/A

!Output Parameters: N/A

!Revision History:
$Log: GEO_input.h,v $
Revision 6.4  2010/06/18 20:41:41  kuyper
Changed GEO_prepare_ancil_data(), GEO_prepare_l1a_data(), and
  GEO_read_L1Apacket_data() to return a status code.
Corrected parameters that were pointers to array by removing
  their leading 'const' qualifiers.

Revision 6.3  2010/04/23 18:32:29  kuyper
Removed parameter from GEO_prepare_ancil_data().

James Kuyper	James.R.Kuyper@NASA.gov

Revision 6.2  2010/04/08 18:42:06  kuyper
Helped address Bug 2472 by adding argument to GEO_prepare_ancil_data.
Simplified interface to GEO_prepare_mirr_data.

Revision 6.1  2010/03/15 20:53:31  gbritz
Add ephemeris/attitude quality information to MOD03 file.

Georgios Britzolakis georgios.britzolakis-1@nasa.gov

Revision 5.1  2005/03/16 21:35:19  kuyper
Changed header guard macro name to avoid reserved name space.

Revision 4.2  2003/08/07 18:48:39  kuyper
Removed never-implemented extra argument.

Revision 4.1  2003/03/13 22:11:08  kuyper
Added frame_data parameter to GEO_prepare_mirr_data().
Added average temperatures parameter to GEO_read_L1Apacket_data().
Added GEO_read_L1Atemp_data().

Revision 3.7  2001/04/13 15:52:48  kuyper
Reverted earth_encoder_times back to uint16.

 * Revision 3.6  2001/04/03  00:24:15  kuyper
 * Changed mirr_params parameter of GEO_prepare_mirr_data to 'const *'.
 *
 * Revision 3.5  2001/04/02  20:36:15  kuyper
 * Removed 'const' from 'l1a_file' parameter; it's not const.
 *
 * Revision 3.4  2001/03/31  16:06:25  seaton
 * Changed type from uint16 to int16 for earth_encoder_times in GEO_read_L1Apacket_data.c
 *
 * Revision 3.3  2001/03/31  15:29:11  seaton
 * Changed declaration for GEO_prepare_l1a_data for geo_params.
 *
 * Revision 3.2  2001/03/28  18:20:10  vlin
 * data type of view_sector_start is changed to int16
 *
 * Revision 3.1  2001/03/20  22:33:35  kuyper
 * Changed *params to const for GEO_cumulate_GRing.
 *
 * Revision 2.9  2001/03/14  18:03:42  seaton
 * Added new arguments to function GEO_read_L1Apacket_data.c
 *
 * Revision 2.8  2001/03/06  20:15:45  vlin
 * 		GEO_prepare_mirr_data() updated & GEO_cumulate_GRing() added
 *

		4/5/95
		Ruiming Chen
		Finished coding

!Team-unique Header:
	        This software is developed by the MODIS Science Data Support
	        Team for the National Aeronautics and Space Administration,
	        Goddard Space Flight Center, under contract NAS5-32373.

!END***************************************************************************
*/

#ifndef GEO_INPUT_H
#define GEO_INPUT_H
#include "hdfi.h"
#include "PGS_PC.h"
#include "PGS_IO_L0.h"
#include "mapi.h"
#include "GEO_parameters.h"

enum {SC_PRIOR, SC_CURRENT};

typedef struct{
	uint8   second_header[MAX_SCAN_NUMBER]
		[sizeof( ((PGSt_IO_L0_SecPktHdrEOS_AM *)NULL)->scTime)];
	int32   posvel[MAX_SCAN_NUMBER][6];
	int16   attit_angvel[MAX_SCAN_NUMBER][6];
} sc_ancil_struct;

/* function prototype */

PGSt_SMF_status GEO_prepare_ancil_data( 
	/* Prepare ancillary data for one granule */
	int const			number_of_scans,
	const GEO_param_struct		* params,
	const sc_ancil_struct		sc_ancillary_data[2],
	const uint16			ss_cp_mode[]
);

PGSt_SMF_status GEO_prepare_l1a_data(
	/* Read and prepare L1a data for one granule */
	MODFILE * const l1a_file,
	GEO_param_struct const * geo_params,
	l1a_data_struct * const l1a_data
);

PGSt_SMF_status GEO_prepare_mirr_data( /* Prepare mirror data for one granule */
	uint16 [MAX_SCAN_NUMBER][ENCODER_LENGTH],
	int16 [MAX_SCAN_NUMBER][SECTOR_LENGTH],
	const mirror_preparation_struct * const,
	double const,
	l1a_data_struct *
);

int GEO_read_L1AECS_metadata(
       MODFILE * const l1a_file,
       ECS_metadata_struct * const ECS_metadata
       );

PGSt_SMF_status GEO_read_L1Apacket_data(
	MODFILE * l1a_file,
	int const number_of_scans,
	uint16 earth_encoder_times[MAX_SCAN_NUMBER][ENCODER_LENGTH],
	sc_ancil_struct sc_ancillary_data[2],
	int16 view_sector_start[MAX_SCAN_NUMBER][SECTOR_LENGTH],
	uint16 FRside[MAX_SCAN_NUMBER][ELEC_SIDES],
	uint16 SAside[MAX_SCAN_NUMBER][ELEC_SIDES],
	uint16 ss_cp_mode[]
);

int GEO_read_L1Ascan_metadata(
	MODFILE  * const l1a_file,
	int const number_of_scans,
	l1a_metadata_struct const * const granule_metadata,
	frame_data_struct frame_data[MAX_SCAN_NUMBER],
	int16 mirr_side[MAX_SCAN_NUMBER]
	);

int GEO_read_L1Aspecific_metadata(
	MODFILE             * const l1a_file,
	l1a_metadata_struct * const granule_metadata,
	int 		    * const number_of_scans
	);

PGSt_SMF_status GEO_read_L1Atemp_data(
	const GEO_param_struct	* const,
	MODFILE			* const,
	int			const,
	float32			[]
	);

PGSt_SMF_status GEO_cumulate_GRing(
    GEO_param_struct const	* const params,
    int32 const			num_frames,
    frame_state_struct const	sc_ev_frame_state[],
    unsigned char		pixel_flags[MAX_DETECTORS][MAX_SCAN_SAMPLE],
    double		ecr_sample_position[MAX_DETECTORS][MAX_SCAN_SAMPLE][3]
);

#endif
