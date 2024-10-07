/*
$Log: imsl_d_lin_sol_gen.c,v $
Revision 4.3  2003/08/28 16:04:04  kuyper
Corrected prolog.

Revision 4.2  2003/07/31 20:56:20  kuyper
Corrected some loop boundary errors.
Replaced incorrect test for ill-condition matrix in ludcmp3() with a more
  correct estimate in imsl_d_lin_sol_gen().
Changed to scale singular rows by 0.0.

Revision 4.1  2003/07/27 20:59:53  kuyper
Initial revision.

James Kuyper Jr. (kuyper@saicmodis.com)
 */
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <float.h>
#include "imsl_wrap.h"
#include "pseudo_imsl.h"
#define THREED 3
#define TINY 1.0E-99

static void ludcmp3(
	double	a[THREED][THREED],
	int	indx[THREED],
	double	*d
){
/******************************************************************************
!C

!Description:
	Replaces the THREED by THREED array 'a' with an L/U decomposition of
	itself. This function may be used in combination with lubksb3() to
	solve linear equations.

!Input Parameters:
	None

!Output Parameters:
	indx	The permutation of the rows that was performed as a result of
		partial pivoting
	d	+/- 1, depending upon whether the number of row interchanges
		in indx is even or odd.

!Input/Output Parameters:
	a	On input, the two-dimensional array to be decomposed. On output,
		the elements of that array below the diagonal are the non-zero
		off-diagonal elements of the array L, and the elements on and
		above the diagonal are the non-zero elements of the array U,
		with the diagonal elements of L being 1.0, such that A = L*U.

Return Parameters:
	None.

Externally Defined:
	IMSL_error_code		"pseudo_imsl.h"
	IMSL_SINGULAR_MATRIX	"imsl_wrap.h"
	THREED			"imsl_d_lin_sol_gen.c"

Called by:
	imsl_d_lin_sol_gen

Routines Called:
	None

!Revision History:
See top of file.

Requirements:
	None

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits
	Based upon the ludcmp() function described in section 2.3 of "Numerical
	Recipes in C", by Preus et. al., 1988.

Design Notes
	Changes from NR version by:
	Use an IMSL return code to indicate errors, rather than calling nrerror.
	Specialize for an array of size THREEDxTHREED, avoiding the need to
	dynamically allocate vv.
	Use 0-based indexing.
	Use double precision variables instead of 'float', with a corresponding
	    change to the value for TINY.
	Automatically adjusts singular arrays to be non-singular.

	Input arguments are not checked, because this is a static function that
	will be only be called by routines in the same source code file, which
	will be responsible for checking them before calling this function.
!END**************************************************************************
*/
    int row, col, k;	/* loop counters. */
    int rowmax;	/* the index of the best row for a partial pivot */
    double big, dum, sum,  temp;	/* intermediate results */
    double vv[THREED]={0.0};	/* the implicit scaling of each row. */

    *d = 1.0;	/* No interchanges yet */

    for(row=0; row<THREED; row++)
    {	/* To get implicit scaling information*/

	big = fabs(a[row][0]);
	for(col=0; col<THREED; col++)
	{ /* find the largest of the absolute values of elements in a[row] */
	    temp = fabs(a[row][col]);
	    if(temp > big)
		big = temp;
	}

	if(big < TINY)
	    IMSL_error_code = IMSL_SINGULAR_MATRIX;
	else
	    vv[row] = 1.0/big;		/* Save the scaling */
    }

    for(col=0; col<THREED; col++)	/* Crout's method*/
    {
	for(row=0; row<col; row++)	/* equation 2.1.13 except for i==j */
	{
	    sum = a[row][col];
	    for(k=0; k<row; k++)
		sum -= a[row][k]*a[k][col];
	    a[row][col] = sum;
	}

	big = 0.0;	/* Initialize for the search for the largest pivot. */
	for(row=col; row<THREED; row++)
	{ /* This is i==j of equation 2.3.12 and i=j+1...N of equation 2.3.13 */
	    sum = a[row][col];
	    for(k=0; k<col; k++)
		sum -= a[row][k]*a[k][col];
	    a[row][col] = sum;
	    dum = vv[row]*fabs(sum);
	    if(dum >= big)
	    { /* Figure of merit for the pivot is better than best so far. */
		big = dum;
		rowmax = row;
	    }
	}

	if(col != rowmax)
	{				/* Need to interchange rows. */
	    for(k=0; k<THREED; k++)
	    {
		temp = a[rowmax][k];
		a[rowmax][k] = a[col][k];
		a[col][k] = temp;
	    }
	    *d *= -1.0;		/* And change the parity of d */
	    vv[rowmax] = vv[col];	/* Also interchange the scale factor. */
	    /* Not swapped?? */
	}

	indx[col] = rowmax;
	if(fabs(a[col][col]) < TINY)
	{
	    /* Remove singularities and near-singularities */
	    a[col][col] = TINY;
	    IMSL_error_code = IMSL_SINGULAR_MATRIX;
	}

	if(col != THREED-1) /* Now, finally, divide by the pivot element */
	    for(row=col+1; row<THREED; row++)
		a[row][col] /= a[col][col];
    }

}

static void lubksb3(
	double		a[THREED][THREED],
	const int	indx[THREED],
	double		b[THREED]
){
/*******************************************************************************
!C
!Description:   
	Solves the set of 3 linear equations A*X = B. The 'a' given as input
	contains the L/U decomposition of A, as performed by ludcmp().

!Input Parameters:
	a	The elements of 'a' below the diagonal are the non-zero
		off-diagonal elements of L; the elements of 'a' on and above
		the diagonal are the non-zero elements of U, with the diagonal
		elements of L being 0.0, such that A = L*U.
	idnx	The permutation vector returned by ludcmp()

!Output Parameters:

!Input/Output Parameters:	(Remove this section if none exist)
	b	On input, contains the values of B. On output, contains the
		values of the corresponding X vector.

Return Parameters:
	None

Externally Defined:
	THREED			"imsl_d_lin_sol_gen.c"

Called by:
	imsl_d_lin_sol_gen	"imsl_wrap.h"

Routines Called:
	None

!Revision History:
See top of file.

Requirements:

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits
	Based upon the lubksb() function described in section 2.3 of "Numerical
	Recipes in C", by Preuss et. al., 1988.

Design Notes
	Changes from NR version:
	Uses 0-based indexing. Specialized for arrays of size THREED.

	Input arguments are not checked, because this is a static function that
	will only be called by routines in the same source code file, which
	will be responsible for checking them before calling this function.
!END**************************************************************************
*/
    double sum;
    int row, col;	/* loop indices */
    int ip, ii=-1;

    for(row=0; row<THREED; row++)
    {
	ip = indx[row];
	sum = b[ip];
	b[ip] = b[row];
	if(ii >= 0)
	    for(col=ii; col<row; col++)
		sum -= a[row][col]*b[col];
	else if(fabs(sum) > TINY)
	    ii = row;

	b[row] = sum;
    }

    for(row = THREED-1; row>=0; row--)
    {
	sum = b[row];
	for(col=row+1; col<THREED; col++)
	    sum -= a[row][col]*b[col];
	b[row] = sum/a[row][row];
    }
}

double * IMSL_DECL imsl_d_lin_sol_gen(
	int	const dim,
	double	matrix[],
	double	not_used[],
	...
){
/*******************************************************************************
!C
!Description:
	Computes the inverse of the matrix provided.

!Input Parameters:
	dim		The number of rows and columns in the matrix. The only
			value permitted by this version is THREED.
	matrix		Array of size dim * dim containing the matrix.
	not_used	Not used, since IMSL_INVERSE_ONLY is mandatory.

!Output Parameters:
	...			Accepts a variable list of options using the C
				<stdarg.h> interface. They come in sets, with
				the first argument of each set being an int that
				identifies the type of set. Multiple occurances
				of each set type are permitted; the later
				values replace the earlier ones. The following
				sets of "optional" arguments are the only ones
				supported by this version, and are both
				mandatory:

	int IMLS_INVERSE_USER	Indicates that the following argument is:
	double inverse[]	A user-allocated array of size dim*dim
				for storing the inverse of matrix.

	int IMSL_INVERSE_ONLY	Compute the inverse of 'matrix'. The argument
				'not used' is ignored, and the return value is
				NULL


!Input/Output Parameters:	(Remove this section if none exist)

Return Values:
	NULL

Externally Defined:
	DBL_EPSILON			<float.h>
	IMSL_DECL			"imsl_wrap.h"
	IMSL_error_code			"pseudo_imsl.h"
	IMSL_OPTIONAL_ARG_NULL_1	"imsl_wrap.h"
	IMSL_ILL_CONDITIONED		"imsl_wrap.h"
	IMSL_INVERSE			"imsl_wrap.h"
	IMSL_INVERSE_ONLY		"imsl_wrap.h"
	IMSL_INVERSE_USER		"imsl_wrap.h"
	IMSL_NEGATIVE_ORDER		"imsl_wrap.h"
	IMSL_UNEXPECTED_NULL_POINTER	"imsl_wrap.h"
	IMSL_UNKNOWN_OPTION		"imsl_wrap.h"

Called by:
	GEO_cumulate_GRing

Routines Called:
	lubksb3		"imsl_d_lin_sol_gen.c"
	ludcmp3		"imsl_d_lin_sol_gen.c"

!Revision History:
See top of file.

Requirements:
	None

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits
	Based upon the documented interface for the correspondingly named real
	IMSL library function.

Design Notes
	Differs from the real IMSL function in that it supports only the
	IMSL_INVERSE_USER and IMSL_INVERSE_ONLY options, and makes both of them
	mandatory, and it only works for dim==THREED.

!END**************************************************************************
*/
    va_list ap;
    double col[THREED];
    int indx[THREED];
    double factor[THREED][THREED];
    double *inverse=NULL;
    int row, c, arg, inverse_only=0;
    double d;
    double biggest=0.0, biginv=0.0, temp;

    IMSL_error_code = (Imsl_code)0;

    /* Collect arguments. */
    va_start(ap, not_used);
    while(( arg=va_arg(ap, int) )) /* = rather than == is deliberate */
    {
	switch(arg)
	{
	    case IMSL_INVERSE_USER:
		inverse = va_arg(ap,double *);	break;

	    case IMSL_INVERSE_ONLY:
		inverse_only = 1;	break;

	    default:
		IMSL_error_code = IMSL_UNKNOWN_OPTION;	break;
	}
    }
    va_end(ap);

    if(IMSL_error_code)
	return NULL;

    if(inverse == NULL)
    {
	IMSL_error_code = IMSL_OPTIONAL_ARG_NULL_1;
	return NULL;
    }

    if(inverse_only == 0)
    {
	IMSL_error_code = IMSL_UNKNOWN_OPTION;
	return NULL;
    }

    if(matrix == NULL)
    {
	IMSL_error_code = IMSL_UNEXPECTED_NULL_POINTER;
	return NULL;
    }

    if(dim <= 0)
    {
	IMSL_error_code = IMSL_NEGATIVE_ORDER;	/* Misnamed message mnemonic. */
	return NULL;
    }

    if(dim != THREED)
    {
	IMSL_error_code = IMSL_INTEGER_OUT_OF_RANGE;
	return NULL;
    }

    memcpy(factor, matrix, sizeof(factor));

    ludcmp3(factor, indx, &d);

    for(row=0; row<THREED; row++)
    {
	for(c=0; c<THREED; c++) 	/* Create identity matrix row */
	    col[c] = (row==c) ? 1.0 : 0.0;

	lubksb3(factor, indx, col);

	for(c=0; c<THREED; c++)
	{
	    inverse[c*THREED+row] = col[c];
	    temp = fabs(col[c]);
	    if(temp > biginv)
		biginv = temp;
	    temp = fabs(matrix[row*THREED+c]);
	    if(temp > biggest)
		biggest = temp;
	}
    }

    /* biggest*biginv is a rough estimate of the condition number of the
     * matrix.
     */
    if(IMSL_error_code == 0 && biggest*DBL_EPSILON*biginv > 1.0)
	IMSL_error_code = IMSL_ILL_CONDITIONED;

    return NULL;
}

