
#include "PGS_MODIS_35251.h"
#include "smfio.h"
#include "GEO_global_arrays.h"
#include "GEO_validation.h"

PGSt_SMF_status GEO_validate_derived_products(
	int const	frame,
	int const	num_detectors,
	double const	range_scale,
	double		frame_to_sensor[][MAX_FRAMES][3],
	uint8		frame_flags[][MAX_FRAMES]
)
/*
!C******************************************************************************
!Description:
	subroutine in validation group of the Level-1A geolocation software to
	validate the derived products.  It checks only that the slant-range to
	the observed frame is within the range capable of being stored in a
	2-byte unsigned integer (for output) and sets the pixel flags
	appropriately. This routine will only be called if the zenith angle and
	azimuth angles for the sun and spacecraft from the frame location have
	been successfully calculated using a standard arctangent function. The
	checks for proper angle range would never fail, therefore, and have
	been removed (as have the flags for these errors in the output product).

!Input Parameters:
	frame		the frame number in the scan
	num_detectors	the number of detectors in the frame
	range_scale	factor (meters/count) for scaling pixel slant-range
			    measurements for output.
	frame_to_sensor	sensor aximuth, zenith and range for each pixel

!Output Parameters:
	frame_flags	geolocation flag for each pixel

Returns:
	MODIS_E_BAD_INPUT_ARG   If there are invalid input arguments,
	PGS_S_SUCCESS           Otherwise.

Externally Defined:
	DETECTORS_1KM           "GEO_geo.h"
	INVALID_SENSOR_ANGLES   "GEO_geo.h"
	MAX_FRAMES              "GEO_geo.h"
	MAX_UINT16_VAL		"GEO_geo.h"
	MODIS_E_BAD_INPUT_ARG   "PGS_MODIS_35251.h"
	PGS_S_SUCCESS           "PGS_SMF.h"
	range_index             "GEO_global_arrays.h"

Called by:
	GEO_locate_one_scan	"GEO_output.h"

Routines Called:
	modsmf			"smfio.h"

!Revision History:
  $Log: GEO_validate_derived_products.c,v $
  Revision 6.1  2009/05/28 22:48:14  kuyper
  Changed macro names to MAX_FRAMES, DETECTORS_1KM.
  Changed sample_number to frame, pixel_flags to frame_flags, and replaced
    "sample" with "frame" in all other variable names.
  Changed frame_to_sensor into an input pointer parameter, and frame_flags
    into an output pointer parameter. Validate both pointers as non-null.
  Changed to return status code.
  Changed name of loop index for improved searchability.
  Cleaned up prolog.

  James Kuyper Jr.	James.R.Kuyper@nasa.gov

  Revision 2.2  1997/11/06 21:48:50  kuyper
  Made implicit conversions explicit.

 * Revision 2.1  1997/10/21  18:16:22  kuyper
 * Returned from ClearCase
 *
  Revision /main/GEO_V2_DEV/3
  checked out 29-Sep-97.23:29:00 by Jeffrey Blanchette (jjb.user@modis-xl)
  by view: jjb ("modispc:/cc/cc_view/jjb.vws")
  "Corrected defects found in code walkthrough:
        Added range_scale parameter.
        Added check of input parameters."

 * Revision /main/GEO_V2_DEV/2  1997/09/15  Ding
 * Removed extraneous code lines.
 * ding@ltpmail.gsfc.nasa.gov

 * Revision 1.6.1.1  1997/07/21  22:20:17  kuyper
 * Merged in out-of-sequence changes.
 *
 * Revision 1.6  1997/07/21  16:24:34  kuyper
 * Baselined Version 1
 *
 *Parallel development:
 * Revision 1.6  1997/04/21  22:38:48  fhliang
 * commented unused argument (LL.10-12).
 *
 * Revision 1.5  1997/03/26  18:15:45  fhliang
 * Initial revision of SDST delivery of GEO_validate_derived_products.c.
 *
 * Revision 1.5  1997/06/02  17:21:26  kuyper
 * Merged seed files.
 *
		Revision 1.4  1996/08/06 19:20:57  kuyper
		Removed #include "PGS_CONSTS.h"

		Revision 1.3  1996/07/24 21:42:03  kuyper
		Standardized order of #include files.
		Declared arguments const.
		Inserted required '!'s in comments.
		Removed ret_val.
		Converted constants to double, to skip implied conversion.

		Revision 1.2  1996/07/18 21:17:12  kuyper
		Included self-checking header file.
		Replaced extern declarations with corresponding header files.
		Added needed header file.
		James Kuyper Jr. <kuyper@lptmail.gsfc.nasa.gov>


		5/22/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Finished coding.

		6/20/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Clean up code and documentation.

		7/3/95
		Tracey W. Holmes (holmes@modis-xl.gsfc.nasa.gov)
		Added SDP error messages.

		9/11/95
		Tracey W. Holmes (holmes@modis-xl.gsfc.nasa.gov)
		Increased limit value.

		10/10/95
		Tracey W. Holmes
		Added DEBUG option.


!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/
{
    int det = 0;   /* loop index */
    char msgbuf[PGS_SMF_MAX_MSGBUF_SIZE];

    if( num_detectors<0 || num_detectors> DETECTORS_1KM || frame<0 ||
	frame>=MAX_FRAMES || frame_to_sensor==NULL || frame_flags==NULL)
    {
	sprintf(msgbuf, "num_detectors = %d, frame = %d, "
	    "frame_to_sensor = %p, frame_flags = %p",
	    num_detectors, frame, (void*)frame_to_sensor, (void*)frame_flags);
	modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf,
	      "GEO_validate_derived_products, GEO_validate_derived_products.c");
	return MODIS_E_BAD_INPUT_ARG;
    }

    for (det = 0; det < num_detectors; det++) /* Range checking */
    {
	if (frame_to_sensor[det][frame][range_index] < 1.0 ||
	    frame_to_sensor[det][frame][range_index] >
	    (double)MAX_UINT16_VAL*range_scale)
	    frame_flags[det][frame] |= INVALID_SENSOR_ANGLES;
    }
    return PGS_S_SUCCESS;
}

