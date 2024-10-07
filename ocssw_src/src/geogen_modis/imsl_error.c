/*
$Log: imsl_error.c,v $
Revision 4.3  2003/08/28 16:06:26  kuyper
Corrected prolog.

Revision 4.2  2003/06/05 15:56:52  kuyper
Corrected lists of externals, and location of revision history.

Revision 4.1  2003/05/19 19:11:28  kuyper
Reverse engineered (with immense simplification) from IMSL library
  documentation.

James Kuyper Jr. kuyper@saicmodis.com
*/

#include <stdarg.h>
#include "imsl_wrap.h"
#include "pseudo_imsl.h"

Imsl_code IMSL_error_code;

#if defined(__alpha) && defined(__osf__)
    int
#else
    long
#endif
    IMSL_DECL imsl_error_code(
	void
)
/*******************************************************************************
!C
!Description:   
	Gets the code corresponding to the error message from the last
	pseudo-IMSL function called.

!Input Parameters:
	None.

!Output Parameters:
	None

Return Value:
	This function returns the error message code from the last pseudo-IMSL
	function called. The include file imsl_wrap.h defines a name for each
	error code.

Externally Defined:
	__alpha			Set by certain compilers
	__osf__			Set by certain compilers
	IMSL_DECL		"imsl_wrap.h"
	IMSL_error_code		imsl_error.c

Called by:
	GEO_cumulate_Gring
	GEO_interp_mirr_enc
	GEO_prepare_mirr_data

Routines Called:
	None

!Revision History:
See top

Requirements:

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits
	Reverse engineered from the "IMSL C/Math/Library", copyright Visual
	Numerics, Inc. 1999-2001

!END**************************************************************************
*/
{
    return IMSL_error_code;
}

void IMSL_DECL imsl_error_options(
	int arg,
	...
)
/*******************************************************************************
!C
!Description:   
	Ignores attempts to set various error handling options. This function
	exists solely to maintain compatibility with code written to use the
	real IMSL function. The only valid arguments for this function are the
	ones that turn off message display, and those arguments are valid
	because message display will be turned off, whether or not this
	function is ever called.
	The only method the real IMSL function has to indicate invalid
	arguments is by printing an error message and aborting the program;
	this is unacceptable behavior for delivered code, so there's absolutely
	nothing useful for this version to do.
	A function was used, rather than a do-nothing function-like macro,
	solely to force error messages if this function is called with no
	arguments, or with a first argument that can't be implicitly converted
	to an 'int'.

!Input Parameters:
	arg	Identifies the type of the first set of optional arguments
	...	Matches an arbitrary list of optional arguments, which will also
		be ignored.

!Output Parameters:
	None

Return Values:
	None

Externally Defined:
	IMSL_DECL		"imsl_wrap.h"

Called by:
	main

Routines Called:
	None

!Revision History:
See top

James Kuyper Jr. kuyper@saicmodis.com

Requirements:
	None

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits
	"IMSL C/Math/Library", copyright 1990-2001 Visual Numerics, Inc.

Design Notes
!END**************************************************************************
*/
{
    /* The only reason there's any body to this function, is to turn off the
     * warning messages about a function that doesn't actually do anything with
     * any of its parameters.
     */
    va_list ap;
    va_start(ap, arg);
    arg = va_arg(ap, int);
    va_end(ap);
}

