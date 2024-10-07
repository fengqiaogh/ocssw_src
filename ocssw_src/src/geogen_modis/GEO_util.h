/*
!C-INC*************************************************************************

!Description:	the .h file for the utility functions in the Level-1A 
                geolocation software

!Input Parameters: N/A

!Output Parameters: N/A

!Revision History:
 * $Log: GEO_util.h,v $
 * Revision 6.1  2011/02/14 21:36:35  kuyper
 * Corrected const-qualification of arguments to GEO_poly_fit().
 *
 * Revision 5.1  2005/03/16 21:36:13  kuyper
 * Changed header guard macro name to avoid reserved name space.
 *
 * Revision 4.1  2002/12/05 19:26:59  kuyper
 * Dropped obsolete chebyshev polynomial functions.
 *
 * Revision 2.3  1999/03/12 17:48:37  kuyper
 * Capitalized Prolog Sections
 *
 * Revision 2.2  1997/10/24  13:25:20  kuyper
 * Removed GEO_lagrange() declaration.
 *
 * Revision 2.1  1997/10/21  18:15:47  kuyper
 * Returned from ClearCase
 *
 * Revision /main/GEO_V2_DEV/2 1997/09/19 kuyper
 * Added CUBIC macro.
 *
 * Revision /main/GEO_V2_DEV/1 1997/08/14 kuyper
 * Added declarations for Chebyshev functions and globals.
 *
 * Revision 1.4  1997/07/18  21:58:00  kuyper
 * Baselined Version 1
 *
 * Revision 1.4  1996/07/24  22:19:42  fhliang
 * Initial revision of SDST delivery of GEO_util.h.
 *
		Revision 1.3  1996/07/24 22:19:42  kuyper
		Inserted required '!'s in comments.
		Declared GEO_mat_vec_mul3, GEO_vec_prod3, GEO_vec_length3 as void.
		Declared arguments const.

		Revision 1.2  1996/07/18 22:07:16  kuyper
		Removed GEO_factorial() and GEO_combination() function declarations.


		4/5/95
		Ruiming Chen
		Finished coding

		6/20/95
		Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
		Removed unused routines

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END***************************************************************************
*/

#ifndef GEO_UTIL_H
#define GEO_UTIL_H
/* Macros.	*/
#define CUBIC 4

/* function prototypes */
int GEO_vec_mul3(              /* cal. the cross product of two 3-dim vecs */
	double vec1[3],              /* input vec 1 */
	double vec2[3],              /* input vec 2 */
	double vec[3]		     /* output vec */
        );

void GEO_mat_vec_mul3(         /* cal. a vec of the mul. of a matrix and  vec */
	double  mat[3][3],           /* input 2-dim matrix */
	double  vec1[3],               /* input vec */
	double  vec[3]		     /* output vec */
        );

int GEO_mat_mul3(              /* cal. matrix of the mul. of 2 matrices */
	double  mat1[3][3],          /* input 2-dim matrix 1 */
	double  mat2[3][3],           /* input 2-dim matrix 2 */
	double  mat[3][3]	     /* output 2-dim matrix */
        );

void GEO_vec_prod3(         /* cal. the dot product of two vecs */
	double  vec1[3],        /* input vec 1 */
	double  vec2[3],        /* input vec 2 */
	double  * const prod	/* output vec */
	);

void GEO_vec_length3(		/* length of vec with size 3*/
	double vec[3],			/* input vec */
	double * const length	/* output length */
	);

int GEO_vec_unit3(          /* normalizes a vector */
	double vec[3],			/* input vec */
	double unit_vec[3]		/* output unit_vec */
	);

int GEO_poly_coef1(         /* determine cubic polynomial coefficients */
	double const var_val1,  /* displacement at t1 */
	double const var_val2,  /* displacement at t2 */
	double const val_derv1, /* derivative at t1 */
	double const val_derv2, /* derivative at t1 */
	double const t1,        /* first time */ 
	double const t2,        /* second time */
	double coef[CUBIC]	/* cubic polynomial coefficients */
	);

int GEO_poly_fit(              		/* evaluate polynomial */
	double const	*,		/* array of polynomial coefficients */
	double		const,		/* independent variable */
	int		const,		/* order of the polynomial */
	double		* const y	/* the pointer to the fitted number */
);

#endif

