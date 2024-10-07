#ifdef NOIMSL
#ifndef IMSL_H
/*
!C-INC*************************************************************************
!Description:
	This header allows us to use reverse-engineered substitutes for the
	IMSL library functions, when running on platforms where that library is
	unavailable, while using the actual functions when they are available.

!Input Parameters:
	None

!Output Parameters:
	None

Return Value:
	None

Externally Defined:
	_WIN32		Defined by the compiler, on 32-bit Windows platforms.

Called by:
	None

Routines Called:
	None

!Revision History:
$Log: imsl_wrap.h,v $
Revision 4.5  2003/08/28 16:06:39  kuyper
Corrected prolog.

Revision 4.4  2003/07/30 20:07:52  kuyper
Corrected value for IMSL_NEGATIVE_ORDER.

Revision 4.3  2003/07/25 19:24:31  kuyper
Changed error code used for bad spline derivatives.
Added code for xvec not increasing.

Revision 4.2  2003/05/21 19:15:41  kuyper
Added code values and option values needed by imsl_d_spline_interp() and
imsl_d_lin_sol_gen().

Revision 4.1  2002/11/06 23:52:28  kuyper
Copied select portions from the IMSL header files.

James Kuyper Jr. kuyper@saicmodis.com

Requirements:
	None

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits:
	"IMSL C/Math/Library", copyright 1990-2001 Visual Numerics, Inc.
	Copied from the imsl.h and imslerr.h header files, with modifications.

!END***************************************************************************
*/
    #define IMSL_H

    #ifdef _WIN32
	#define IMSL_DECL	__cdecl
    #else
	#define IMSL_DECL
    #endif

    typedef enum {
	    IMSL_UNKNOWN_OPTION			=    103,
	    IMSL_OPTIONAL_ARG_NULL_1		=    104,
	    IMSL_INTEGER_OUT_OF_RANGE		=    132,
	    IMSL_UNEXPECTED_NULL_POINTER	=    150,
	    IMSL_OUT_OF_MEMORY			=    200,
	    IMSL_OUT_OF_MEMORY_2		=    202,
	    IMSL_ILL_CONDITIONED		=   1003,
	    IMSL_SINGULAR_MATRIX		=   1004,
	    IMSL_NEGATIVE_ORDER			=   1010,
	    IMSL_XDATA_NOT_INCREASING		=   3024,
	    IMSL_KNOT_MULTIPLICITY		=   3028,
	    IMSL_KNOT_NOT_INCREASING		=   3029,
	    IMSL_SPLINE_BAD_ORDER		=   3031,
	    IMSL_SPLINE_BAD_COEFFS		=   3032,
	    IMSL_SPLINE_ORDER_DERIV		=   3033,
	    IMSL_DUPLICATE_XDATA_VALUES		=   3034,
	    IMSL_SPLINE_NEED_DATA_PTS		=   3035,
	    IMSL_XDATA_TOO_LARGE		=   3052,
	    IMSL_XDATA_TOO_SMALL		=   3053,
	    IMSL_KNOT_DATA_INTERLACING		=   3101,
	    IMSL_XVEC_NOT_INCREASING		=   3122,
	    IMSL_XVEC_LENGTH			=   3123
    } Imsl_code;

    typedef enum {
	    IMSL_WARNING	= 3,
	    IMSL_FATAL		= 4,
	    IMSL_TERMINAL	= 5
    } Imsl_error;

    typedef struct {
	int	domain_dim;
	int	target_dim;
	int	*order;
	int	*num_coef;
	int     *num_knots;
	double	**knots;
	double	**coef;
    } Imsl_d_spline;

    double * IMSL_DECL	imsl_d_lin_sol_gen(int, double*, double*, ...);
    Imsl_d_spline * IMSL_DECL
			imsl_d_spline_interp(int, double[], double[], ...);
    double IMSL_DECL	imsl_d_spline_value(double, Imsl_d_spline*, ...);

#if defined(__alpha) && defined(__osf__)
    int
#else
    long
#endif
    IMSL_DECL	imsl_error_code(void);

    void IMSL_DECL	imsl_error_options (int, ...);
    double IMSL_DECL	imsl_d_machine(int);

    enum Imsl_keyword {
	    IMSL_DERIV		= 10028,
	    IMSL_KNOTS         	= 10035,
	    IMSL_ORDER         	= 10036,
	    IMSL_INVERSE	= 10152,
	    IMSL_INVERSE_USER	= 10153,
	    IMSL_INVERSE_ONLY	= 10155,
	    IMSL_SET_PRINT	= 10188,
	    IMSL_SET_STOP	= 10190,
	    IMSL_GRID_USER	= 11051
    };
    #endif /* IMSL_H */

#else

    #include "imsl.h"

#endif /* NOIMSL */

