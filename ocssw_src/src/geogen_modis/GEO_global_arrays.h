/*
!C-INC**************************************************************************
!Description:	the declaration of the global array variables for the
		Level-1A geolocation software

!Input Parameters: N/A

!Output Parameters: N/A

!Revision History:
 * $Log: GEO_global_arrays.h,v $
 * Revision 6.4  2010/12/15 22:56:47  kuyper
 * Merged *_IDX macros with *_index enumeration constants into *_IDX
 *   enumeration constants.
 *
 * Revision 6.3  2010/05/03 18:15:38  kuyper
 * Reinstated xyzRotRates member of sc_state_struct.
 *
 * James Kuyper		James.R.Kuyper@NASA.gov
 *
 * Revision 6.2  2010/04/01 19:53:32  kuyper
 * Helped resolve Bug 2472 by adding qualityFlags to sc_state_struct.
 *
 * Revision 6.1  2009/04/17 16:24:57  kuyper
 * Removed most global arrays.
 * Changed MAX_SCAN_SAMPLE to MAX_PADDED.
 *
 * Revision 4.1  2003/08/22 14:10:25  kuyper
 * Removed sc_attitude.
 *
 * Revision 3.1  2001/04/03 00:25:26  kuyper
 * Removed obsolete global mirr_side.
 *
 * Revision 2.6  2001/03/01  19:38:34  rtian
 * change mirr_impulse_enc and mirr_impulse_time to be array of pointer to double.
 * rtian@gscmail.gsfc.nasa.gov
 *
 * Revision 2.5  1999/03/12  17:48:37  kuyper
 * Capitalized Prolog Sections
 *
 * Revision 2.4  1999/01/22  17:00:19  kuyper
 * Removed sc_evolution and it's typedef.
 * Removed xyzRotRates from sc_state_struct.
 *
 * Revision 2.3  1997/11/07  16:11:21  kuyper
 * Changed pixel_flags to unsigned char, to match gflags.
 * Made include guard macro names ANSI compliant.
 *
 * Revision 2.2  1997/10/24  13:25:57  kuyper
 * Removed ecs_sc_* arrays.
 *
 * Revision 2.1  1997/10/21  18:15:47  kuyper
 * Returned from ClearCase
 *
  Revision /main/GEO_V2_DEV/2 on 30-Sep-97.00:43:14
	by Jeffrey Blanchette (jjb@modis-xl.gsfc.nasa.gov)
  "Added sample_to_sensor enumerated indices.
   Removed range_scale global."

 * Revision 1.8  1997/07/18  21:58:00  kuyper
 * Baselined Version 1
 *
 * Revision 1.8  1997/03/26  19:09:37  fhliang
 * Initial revision of SDST delivery of GEO_global_arrays.h.
 *
		Revision 1.7  1997/01/02 16:29:05  kuyper
		Added sc_angvel.

		Revision 1.6  1996/11/07 16:19:22  kuyper
		Added array index macros for lat,lon,height vectors.

		Revision 1.5  1996/09/17 19:33:22  kuyper
		Moved various global arrays into sc_evolution struct.

		Revision 1.4  1996/07/24 22:14:31  kuyper
		Inserted required '!'s in comments.

		Revision 1.3  1996/07/18 22:14:41  kuyper
		Included GEO_geo.h for MAX_DETECTORS definition.

		Revision 1.2  1996/07/18 21:56:13  kuyper
		Converted definitions to extern declarations.
		Removed samples_per_scan, sector_start_enc.


		5/22/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Created this include file from GEO_global.h

		6/21/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Clean up and verify consistency with code

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END****************************************************************************
*/


#ifndef GEO_GLOBAL_ARRAYS_H
#define GEO_GLOBAL_ARRAYS_H
#include "hdfi.h"
#include "PGS_TYPES.h"
#include "GEO_geo.h"


/**********************************************************
variables for one granule of the input data based on scan
**********************************************************/

/* Earth sector start time for the scans in a granule */
extern double scan_start_time[MAX_SCAN_NUMBER];

enum { EPH_IDX, ATT_IDX, QFL_IDXS};
typedef struct {
        PGSt_double	time;		/* record time, TAI seconds */
        PGSt_double	position[3];	/* ECI postion, meters */
        PGSt_double	velocity[3];	/* ECI velocity, meters/sec */
        PGSt_double	eulerAngles[3];	/* attitude euler angles (roll, pitch,
	    				  yaw) wrt the orbital reference frame
					  (?), radians. */
	PGSt_double	xyzRotRates[3];	/* attitude rates (roll, pitch, yaw)
					  in the spacecraft reference frame.
					  radians/sec. */
	uint32		qualityFlags[QFL_IDXS]; /* Ephemeris and atittude
						  quality flags. */
} sc_state_struct;

/* total number of good impulse data */
extern int num_impulse[MAX_SCAN_NUMBER];

/* validation flags for mirror data (each scan) */
extern int mirr_impulse_flag[MAX_SCAN_NUMBER];

/* mirror impulse encoder for the scans in a granule */
extern double * mirr_impulse_enc[MAX_SCAN_NUMBER];

/* mirror impulse time for scans in a granule */
extern double * mirr_impulse_time[MAX_SCAN_NUMBER];

/* sample time for all the pixels in a scan (relative to sector start time */
extern double *sample_time;

/* offsets into terrain_*_position */
enum {LAT_IDX, LON_IDX, HT_IDX};

/* offsets into sensor, solar arrays. */
enum {azimuth_index, z_angle_index, range_index};

#endif

