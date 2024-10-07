/* file: GEO_basic.h */

/*
!C-INC************************************************************************
!Description:   the basic .h file for the  Level-1A geolocation
                software

!Input Parameters: N/A

!Output Parameters: N/A

!Revision History:
		$Log: GEO_basic.h,v $
		Revision 5.1  2005/03/16 21:34:52  kuyper
		Changed header guard macro name to avoid reserved name space.

		Revision 2.3  1999/03/12 17:48:37  kuyper
		Capitalized Prolog Sections

 * Revision 2.2  1999/01/22  16:58:05  kuyper
 * Added definition of WARNING.
 *
 * Revision 2.1  1997/10/21  18:15:47  kuyper
 * Returned from ClearCase
 *
 * Revision 1.7.1.1  1997/07/18  22:40:58  kuyper
 * Merged in out-of-sequence changes.
 *
 * Revision 1.7  1997/07/18  21:58:00  kuyper
 * Baselined Version 1
 *
 *Parallel development:
 * Revision 1.6  1997/05/12  18:18:56  kuyper
 * Changed definition of TRUE for compatibility
 * with n32 and 64 bit versions of Toolkit.
 *
 * Revision 1.6  1997/03/26  19:08:41  fhliang
 * Initial revision of SDST delivery of GEO_basic.h.
 *
		Revision 1.5  1996/10/16 22:09:34  kuyper
		Modified definition of FAIL to be compatible with hdf.h.

		Revision 1.4  1996/08/06 15:43:24  kuyper
		Changed to more conventional definitions of TRUE/FALSE; only external effect
		is the value of correct_terrain in parameters file.

		Revision 1.3  1996/07/30 19:09:45  kuyper
		Added definition of EXIT_FAIL, to comply with standards.


		4/5/95
		Ruiming Chen (rchen@ltpmail.gsfc.nasa.gov)
		Finished coding

!Team-unique Header:
		This software is developed by the MODIS Science Data Support
		Team for the National Aeronautics and Space Administration,
		Goddard Space Flight Center, under contract NAS5-32373.

!END**************************************************************************
*/

#ifndef GEO_BASIC_H
#define GEO_BASIC_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define EXIT_FAIL	1	/* Per MODIS software standard 5.1.1 .	*/
#define EXIT_SUCCESS	0

#define FAIL	(-1)
#define SUCCESS 0
#define WARNING 1

#ifndef FALSE
#define FALSE	0
#endif
#ifndef TRUE
#define TRUE	(!FALSE)
#endif

#endif


