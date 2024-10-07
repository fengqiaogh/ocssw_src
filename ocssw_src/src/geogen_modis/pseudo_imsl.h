#ifndef PSEUDO_IMSL_H
#define PSEUDO_IMSL_H
/*
!C-INC*************************************************************************
!Description:
	This header file contains declarations that are shared between the
	different modules that make up the pseudo-IMSL library, but which are
	not meant to be accessible from outside of it.

!Input Parameters:
	N/A

!Output Parameters:
	N/A

Return Value:
	N/A

Externally Defined:
	None

Called by:
	N/A

Routines Called:
	N/A

!Revision History:
$Log: pseudo_imsl.h,v $
Revision 4.2  2003/08/28 16:07:33  kuyper
Corrected prolog.

Revision 4.1  2003/05/19 18:46:27  kuyper
Initial Revision.

James Kuyper Jr. kuyper@saicmodis.com

Requirements:
	None

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits:

!END***************************************************************************
*/

/* The maximum spline order that can be handled by the spline routines. */
#define MAX_SPLINE_ORDER 4

/* The global variable used to store the exit status for the most recently
 * called pseudo-IMSL function.
 */
extern Imsl_code IMSL_error_code;

#endif
