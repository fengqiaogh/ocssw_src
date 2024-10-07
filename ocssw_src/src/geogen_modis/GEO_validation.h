#ifndef GEO_VALIDATION_H
#define GEO_VALIDATION_H

/*
!C-INC*************************************************************************

!Description:	the prototypes for the validation functions in the Level-1A 
                geolocation software

!Input Parameters: N/A

!Output Parameters: N/A

!Revision History:
 * $Log: GEO_validation.h,v $
 * Revision 6.3  2010/04/23 18:23:00  kuyper
 * Removed obsolete functions.
 *
 * Revision 6.2  2009/03/27 20:40:31  kuyper
 * Removed pointless leading dimensions of "array" parameters.
 *
 * Revision 6.1  2009/03/18 21:02:20  kuyper
 * Removed MAX_DETECTORS, changed MAX_SCAN_SAMPLE to PADDED_SAMPLES.
 * Changed return type for GEO_validate_derived_products() and
 *   GEO_validate_earth_location().
 * Changed parameter list for GEO_validate_derived_products().
 *
 * James Kuyper Jr.	James.R.Kuyper@nasa.gov
 *
 * Revision 2.6  2000/03/22 16:42:17  lma
 * delete qa_metadata arg in GEO_validate_earth_location definition
 *
 * Revision 2.5  2000/01/28  16:14:35  lma
 * modified BAD_DATA from -1 to 1
 *
 * Revision 2.4  1999/04/09  15:18:20  kuyper
 * Removed duplicate Log keyword, and the duplicate revision history created by
 * it.
 *
 * Revision 2.3  1999/03/12  17:48:37  kuyper
 * Capitalized Prolog Sections
 *
 * Revision 2.2  1997/11/07  16:10:22  kuyper
 * Changed pixel_flags to unsigned char, to match gflags.
 *
 * Revision 2.1  1997/10/21  18:15:47  kuyper
 * Returned from ClearCase
 *
  Revised /main/GEO_V2_DEV/2 on 30-Sep-97.00:35:12 
	by Jeffrey Blanchette (jjb@modis-xl.gsfc.nasa.gov)
  "Revised GEO_validate_derived_products's prototype."

 * Revision 1.5.1.1  1997/07/18  22:40:58  kuyper
 * Merged in out-of-sequence changes.
 *
 * Revision 1.5  1997/07/18  21:58:00  kuyper
 * Baselined Version 1
 *
 * Revision 1.3  1997/04/22  18:29:22  fhliang
 * commented unused argument (LL.67-69).
 *
 * Revision 1.2  1997/04/11  19:37:34  fhliang
 * added log header
 *
 * Revision 1.1  1997/03/26  19:12:15  fhliang
 *Initial revision of SDST delivery of GEO_validation.h
 *
		6/20/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Finished coding


!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END***************************************************************************
*/

#define	BAD_DATA	1
#define GOOD_DATA	0
#include "GEO_parameters.h"

/* function prototype */

int GEO_abs_limit_check (           /* Check samples against absolute limits */
        double data_samples[],
        int const number_of_samples,
        double data_limits[],
        int sample_flags[]
        );

int GEO_del_limit_check (           /* Check samples against delta limits */
        double data_samples[],
        int const number_of_samples,
        double const delta_limit,
        int sample_flags[]
	);

int GEO_find_next_flag (            /* Find next unflagged sample */
        int sample_flags[],
        int const number_of_samples,
        int const start_sample
        );

PGSt_SMF_status GEO_validate_derived_products (
	/* Validate derived geolocation values */
        int const	sample_number,
        int const	num_detectors,
	double const	range_scale,
	double		frame_to_sensor[][MAX_FRAMES][3],
	uint8		frame_flags[][MAX_FRAMES]
);

PGSt_SMF_status GEO_validate_earth_location ( /* Validate geolocation values */
	int		const mirr_impulse_flag,
        double		terrain_sample_position[][MAX_PADDED][3],
	int		const frame_number,
	int		const num_detectors,
	unsigned char	sample_flags[][MAX_PADDED]
);

#endif

