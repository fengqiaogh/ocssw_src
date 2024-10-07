/* file: GEO_find_next_flag.c */

#include "GEO_geo.h"
#include "GEO_validation.h"

int GEO_find_next_flag(
	int sample_flags[],
	int const number_of_samples, 
	int const start_sample
	)
{
/*
!C*****************************************************************************
!Description:   
	        Routine in Validation group of the Level-1A geolocation
                software to find the next unflagged sample in a flag array.

!Input Parameters:
                int sample_flags [] - array of flags for data samples
		int number_of_samples - array size
		int start_sample - start index for search

!Output Parameters:

Return parameter:
                int - index of next unflagged value if found, or
		  FAIL if end of array reached.

Global variables:
		None

Call functions:
		None

!Revision History:
		$Log: GEO_find_next_flag.c,v $
		Revision 1.6  1997/07/21 16:24:34  kuyper
		Baselined Version 1

 * Revision 1.6  1997/03/26  18:04:07  fhliang
 * Initial revision of SDST delivery of GEO_find_next_flag.c.
 *
		Revision 1.5  1996/07/29 14:54:25  kuyper
		Used GOOD_DATA for flags, rather than SUCCESS (value unchanged).

		Revision 1.4  1996/07/24  21:00:34  kuyper
		Declared arguments const.

		Revision 1.3  1996/07/23 23:02:03  kuyper
		Inserted required '!' in comments.

		Revision 1.2  1996/07/18 16:24:14  kuyper
		Included self-checking header file.


        6/15/95
        Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
        Finished coding.

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/

  /* Local variable declarations */

  int ind = 0; /* index */


  /* Begin program logic */

  ind = start_sample;
  for (ind = start_sample; ind < number_of_samples; ind++) {
    if (sample_flags[ind] == GOOD_DATA)
      return ind;
  }

  return FAIL;
}



