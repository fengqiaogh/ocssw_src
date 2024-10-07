/* file: GEO_get_sample_time.c */

#include "GEO_geo.h"
#include "GEO_global_arrays.h"
#include "GEO_inst.h"
#include "GEO_parameters.h"
#include "GEO_output.h"
#include "PGS_MODIS_35251.h"

PGSt_double GEO_get_sample_time( 
        focal_plane_geometry_struct const *geometry_params,
	const int sample_number
	)
{
/*
!C*****************************************************************************
!Description:   subroutine in Instrument Model of the Level-1A geolocation 
                software to calculate the sample time offset of a pixel from
		the sector start time.  It incorporates the time offset from
		the first band collected (Band 30) to the selected band and 
		also corrects for the latch-up time (trailing edge) to the 
		sample center time.  The equations are derived in the 
		Geolocation ATBD (pp. 32-33).  

!Input Parameters: 
                geometry_params - Numbers describing the MODIS instrument
		sample_number - the sample number 
		
!Output Parameters:
                None

Return values:
	        The sample time offset in seconds from the EV sector start time.	

Externally Defined:
                FIRST_BAND              "GEO_geo.h"
                GEO_DOUBLE_FILLVALUE    "GEO_output.h"
                LATCH_TO_CENTER         "GEO_geo.h"
                MODIS_E_BAD_INPUT_ARG   "PGS_MODIS_35251.h"

Call functions: None

Called by:      GEO_locate_one_scan     "GEO_output.h"

!Revision History:
		$Log: GEO_get_sample_time.c,v $
		Revision 6.2  2009/05/15 16:42:35  xgeng
		Corrected the formula to calculate the sample time

                Xu Geng (xu.geng@saic.com)
		Revision 6.1  2009/05/06 19:10:18  xgeng
		Changed all global variables into pointer parameters.
		Considered sample time for high resolution bands.
            
		Revision 2.2  1999/02/08 22:31:12  kuyper
		Corrected to reflect fact that EV start time corresponds to the start of
		  integration for band 30, not the end.

 * Revision 2.1  1997/10/21  18:16:22  kuyper
 * Returned from ClearCase
 *
 * Revision 1.5.1.1  1997/07/21  22:20:17  kuyper
 * Merged in out-of-sequence changes.
 *
 * Revision 1.5  1997/07/21  16:24:34  kuyper
 * Baselined Version 1
 *
 *Parallel development:
 * Revision 1.6  1997/04/21  22:33:23  fhliang
 * commented unused argument (LL.12-14).
 *
 * Revision 1.5  1997/03/26  18:06:20  fhliang
 * Initial revision of SDST delivery of GEO_get_sample_time.c.
 *
		Revision 1.4  1996/07/24 21:02:24  kuyper
		Declared arguments const.

		Revision 1.3  1996/07/23 23:06:58  kuyper
		Inserted required '!' in comments.
		Changed constants to double to avoid conversions.
		Made implicit casts explicit.

		Revision 1.2  1996/07/18 16:40:53  kuyper
		Included self-checking header file.
		Replaced extern declaration of sample_time with definition.
		Replaced other extern declarations with corresponding header files.


		4/25/95
		Ruiming Chen
		Finished coding.

		6/12/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Revised to change hardcoded values to defined constants 
		and improve documentation.

		6/22/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Added band_number as global variable.

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/			

  double t_sample = 1.0/3000.0; 
  char filefunc[] = __FILE__ ", GEO_get_sample_time"; 

  if (geometry_params == NULL) { 
    modsmf(MODIS_E_BAD_INPUT_ARG, "geometry_params is null", filefunc);
    return GEO_DOUBLE_FILLVALUE;
  }

  /* calculate the sample time */
  t_sample = geometry_params->t_frame / 
     (double)geometry_params->N_samp[geometry_params->band_number];

  return t_sample*(double)(sample_number + 1 -  
         geometry_params->N_samp[geometry_params->band_number]) +  
         geometry_params->t_frame*(
               geometry_params->band_position[geometry_params->band_number] - 
               geometry_params->F_offset[geometry_params->band_number] -  
	       geometry_params->band_position[FIRST_BAND]) + 
         LATCH_TO_CENTER*(t_sample + geometry_params->t_reset);

} 

