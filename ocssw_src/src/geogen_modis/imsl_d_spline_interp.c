#include <stddef.h>
#include <math.h>
#include <float.h>
#include <stdarg.h>
#include "imsl_wrap.h"
#include "pseudo_imsl.h"
#include "GEO_geo.h"

/*
Revision History:
$Log: imsl_d_spline_interp.c,v $
Revision 6.1  2010/06/18 20:25:23  kuyper
Corrected declaration of a banslv() parameter that is a pointer to an array by
  removing 'const'.

Revision 4.3  2003/08/28 16:21:03  kuyper
Corrected prolog.

Revision 4.2  2003/07/01 20:11:25  kuyper
Corrected off-by-one errors in bsplvb(), banfac(), and splint().
Corrected descriptive text.
Changed imsl_d_spline_interp() to allocate memory as soon as the numbers
  controlling the amount of memory needed have been validate, and to return
  a pointer to that memory even when error conditions occur.
Corrected location pointed at by pspline->coef.

Revision 4.1  2003/05/22 19:57:00  kuyper
Reverse-engineered from documentation and testing.

James Kuyper Jr. (kuyper@saicmodis.com)
 */

static void bsplvb(
	const double	t[],
	int		jhigh,
	double		x,
	int		left,
	double		biatx[]
)

/*******************************************************************************
!C
!Description:   
	Calculates the value of all possibly nonzero b-splines at x of order
	jhigh with knot sequence t.

!Input Parameters:
	t	Knot sequence of length left+jhigh, assumed to be nondecreasing.
		ASSUMPTION: t[left] < t[left+1]
		DIVISION BY ZERO will result if t[left] == t[left+1]
	jhigh	The order of the b-splines whose values at x are to be returned.
		WARNING: the restriction jhigh<=MAX_SPLINE_ORDER is imposed
		arbitrarily by the declaration of deltal and deltar below, but
		is NOWHERE CHECKED for.
	x	The point at which the b-splines are to be evaluated.
	left	An integer chosen (usually) so that
			t[left] <= x < t[left+1]
		
!Output Parameters:
	biatx	An array of length jhigh, with biatx[i] containing the value at
		x of the polynomial of order jhigh which agrees with the
		b-spline b[left-jhigh+i][jout](t) on the interval (t[left],
		t[left+1]).

Return Value:
	None

Externally Defined:
	MAX_SPLINE_ORDER	"pseudo_imsl.h"

Called by:
	splint()

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
	This routine is a port to C of the Fortran program of the same name as
	described by Carl de Boor in "A Practical Guide to Splines", 2001. The
	original routine can be found in bsplvb.f, which is part of the pppack
	package which can be downloaded from
	<ftp://netlib.bell-labs.com/netlib/pppack.tar>

Design Notes
	The recurrence relation
			 x - t[i]		  t[i+j+1] - x
	b[i][j+1](x) = -----------b[i][j+1](x) + ---------------b[i+1][j](x)
		       t[i+j]-t[i]		 t[i+j+1]-t[i+1]

	is used (repeatedly) to generate the (j+1)-vector b[left-j][j+1](x),
	...,b[left][j+1](x) from the j-vector b[left-j+1][j](x),...,
	b[left][j](x), storing the new values in biatx over the old. The facts
	that
		b[i][0] == 1 if t[i] <= x < t[i+1]
	and that
		b[i][j](x) == 0  unless t[i] <= x < t[i+j]
	are used. The particular organization of the calculations follows
	algorithm (8) in chapter X of de Boor (2001).

	Unlike the original, this version does not support incremental
	calculation of the spline values.
!END**************************************************************************
*/
{
    int j;
    int i, jp1;
    double saved, term;
    double deltal[MAX_SPLINE_ORDER], deltar[MAX_SPLINE_ORDER]; 

    biatx[0] = 1.0;

    for(j=0; j < jhigh-1; j=jp1)
    {
	jp1 = j+1;
	deltar[j] = t[left+jp1] - x;
	deltal[j] = x - t[left-j];
	saved = 0.0;

	for(i=0; i<=j; i++)
	{
	    term = biatx[i]/(deltar[i] + deltal[j-i]);
	    biatx[i] = saved + deltar[i]*term;
	    saved = deltal[j-i]*term;
	}

	biatx[jp1] = saved;
    }
}


static void banfac(
	double	*w[/* nrow */],
	int	nrow,
	int	nband
)

/*******************************************************************************
!C
!Description:   
	From  * A Practical Guide to Splines *  by C. de Boor
	Returns in w the LU-factorization (without pivoting) of the banded
	matrix 'a'  of order nrow with (2*nband + 1) bands or diagonals stored
	in the work array w.

!Input Parameters:
	nrow	Number of rows of data, assumed to be positive.
	nband	Number of bands of a above and below the main diagonal

!Output Parameters:
	None

!Input/Output Parameters:
	w	On input, contains nrow pointers to arrays containing the
		interesting part of a banded matrix 'a', with the diagonals or
		bands of 'a' stored in the columns of w, while rows of 'a'
		correspond to rows of w. This is the storage mode used in
		linpak and results in efficient innermost loops.
		Explicitly, 'a' has nband bands below the diagonal
			       +      1   (main) diagonal
			       +    nband bands above the diagonal
		and thus, a[j][i+j] is in w[j][i+nband]
		for i=-nband,...,nband-1 and j=0,...,nrow-1.
		On output, contains the LU-factorization of  'a'  into a unit
		lower triangular matrix L and an upper triangular matrix U
		(both banded) and stored in customary fashion over the
		corresponding entries of  a. This makes it possible to solve
		any particular linear system  a*x = b  for  x  by calling

			banslv ( w, nrow, nband, b )

		with the solution x contained in b on return .

Return Values:
	N/A

Externally Defined:
	None

Called by:
	splint()

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

References and Credits:
	This routine is a port to C of the Fortran program of the same name as
	described by Carl de Boor in "A Practical Guide to Splines", 2001. The
	original routine can be found in banfac.f, which is part of the pppack
	package which can be downloaded from
	<ftp://netlib.bell-labs.com/netlib/pppack.tar>

Design Notes:
	Gauss eliminitation WITHOUT pivoting is used. The routine is intended
	for use with matrices 'a' which do not require row interchanges during
	factorization, especially for the TOTALLY POSITIVE matrices which occur
	in spline calculations. The routine should not be used for an arbitrary
	banded matrix.

	This routine differs from the original in that it's specialized for a
	banded matrix with the same number of bands above and below the
	diagonal, and that nrow is positive.
	The nroww parameter was dropped, because the C version of the code does
	not need it for the reasons that the fortran version did.
	imsl_d_spline_interp() performs validity tests on the inputs which are
	sufficient to guarantee that the matrix will not be singular.
	Therefore, singularity checks have been removed, and so has the status
	argument, which was used only to report the results of those checks.
	Finally, since we can't use C99 yet, the 'w' argument can't be handled
	as a variable-length array, which would have been an exact match to the
	way it's handled in the Fortran version.

	This routine assumes that 'w' and all the pointers in the array it
	points at are valid, and that the integral arguments all correctly
	describe the memory that that those pointers point at. Ensuring these
	validity conditions is the responsibility of the calling function. This
	routine is declared static, so as to ensure that it can only be called
	by functions that ensure the validity of those assumptions.
!END**************************************************************************
*/
{
    int i, j, k, ipk, jmax, midmk;
    int nrowm1=nrow-1;
    double factor, pivot;

    if(nrowm1<1 || nband <1)
	return;

    /* 'a' in not just a triangular matrix. Construct LU factorization. */
    for(i=0; i<nrowm1; i++)
    {
	pivot = w[i][nband];
	jmax = (nband < nrow-i-1) ? nband : nrow-i-1;
	for(j=0; j<jmax; j++)
	    w[i][nband+j+1] /= pivot;

	for(k=0; k<jmax; k++)
	{
	    ipk = i+k+1;
	    midmk = nband-k;
	    factor = w[ipk][midmk-1];
	    for(j=0; j<jmax; j++)
		w[ipk][midmk+j] -= w[i][nband+j+1]*factor;
	}
    }

    return;
}


static void banslv(
	double	*w[/* nrow */],
	int	nrow,
	int	nband,
	double	b[/* nrow */]
)

/*******************************************************************************
!C
!Description:   
	Companion function to banfac(). It returns the solution x of the
	linear system a*x = b in place of b, given the LU-factorization for 'a'
	in the workarray w.

!Input Parameters:
	w	An array of nrow pointers to double precision arrays containing
		the LU-factorization of a banded matrix 'a' of order nrow, as
		constructed in banfac(). For details, see banfac.c
	nrow	The number of rows in w.
	nband	The number of bands above and below the main diagonal in 'a'.

!Output Parameters:
	None

!Input/Output Parameters:
	b	Input: Right side of the system to be solved.
		Output: Solution vector

Return Parameters:
	N/A

Externally Defined:
	None

Called by:
	splint()

Routines Called:

!Revision History:
See top of file

Requirements:
	None

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits
	This routine is a port to C of the Fortran program of the same name as
	described by Carl de Boor in "A Practical Guide to Splines", 2001. The
	original routine can be found in banslv.f, which is part of the pppack
	package which can be downloaded from
	<ftp://netlib.bell-labs.com/netlib/pppack.tar>

Design Notes
	(With 'a'=L*U, as stored in w,) the unit lower triangular system
	L(U*x) = b is solved for y = u*x, and y stored in b. Then the upper
	triangular system u*x = y is solved for x. The calculations are
	arranged so that the innermost loops stay within rows.

	This routine differs from the original in that it's designed for arrays
	where the number of bands above and below the diagonal are equal. This
	simplifies our testing, since that's the only kind that
	imsl_d_spline_interp() will be passing to this routine.
	The nroww parameter was dropped, because the C version of the code does
	not need it for the purpose that the fortran version did.

	This function assumes that it's pointer arguments are valid pointers,
	that it's numeric arguments correctly describe the arrays that the
	pointers point at, and that w[0][nband] is not too close to 0.0. For
	nband<=0, it assumes that w[i][0] is not too close to 0.0. Those
	validity checks are the responsibility of the calling routine. This
	function is declared 'static', to ensure that it can only be called 
	directly by routines which do ensure those requirements.
!END**************************************************************************
*/
{
    int i, j, jmax, nrowm1;

    if(nrow != 1)
    {
	nrowm1 = nrow-1;

	if(nband>0)
	{
	    /* Forward pass. */
	    for(i=0; i<nrowm1; i++)
	    {
		jmax = (nband<nrow-i-1)? nband : nrow-i-1;
		for(j=1; j<=jmax; j++)
		    b[i+j] -= b[i]*w[i][nband+j];
	    }
	}
	else
	{
	    for(i=0; i<nrow; i++)
		b[i] /= w[i][0];

	    return;
	}

	/* Backward pass. */
	for(i=nrow-1; i>0; i--)
	{
	    b[i] /= w[i][nband];
	    jmax = (nband<i)? nband : i;
	    for(j=1; j<=jmax; j++)
		b[i-j] -= b[i]*w[i][nband-j];
	}
    }

    b[0] /= w[0][nband];

}


static Imsl_code splint (
	const double	tau[],
	const double	gtau[],
	const double	t[],
	int		n,
	int		k,
	double		*q[],
	double		bcoef[]
)
/*******************************************************************************
!C
!Description:   Computes a spline interpolant.

!Input Parameters:
	tau		Array of length n, containing data point abscissae
			ASSUMPTION: tau is strictly increasing
	gtau		Corresponding array of length n, containing data point
			ordinates.
	t		Knot sequence, of length n+k
	n		Number of data points, and dimension of the spline
			space s(k,t)
	k		Order of the spline.

!Output Parameters:
	q		Array of n pointers to double arrays of length 2*k-1,
			containing the triangular factorization of the
			coefficient matrix of the linear system for the
			b-coefficients of the spline interpolant.
			The b-coeffs for the interpolant of an additional data
			set (tau[i], htau[i]), i=1,...,n with the same data
			abscissae can be obtained without going through all of
			the calculations in this routine, simply by loading
			htau into bcoef and then executing the call
			banslv(q, n, k-1, bcoef).
	bcoef		The b-coefficients of the interpolant, of length n.


Return Value:
	0				Success
	IMSL_KNOT_DATA_INTERLACING	The data and the knots are not properly
					interlaced

Externally Defined:
	IMSL_KNOT_DATA_INTERLACING	"imsl_wrap.h"

Called by:
	imsl_d_spline_interp()

Routines Called:
	bsplvb				imsl_d_spline_interp.c
	banfac				imsl_d_spline_interp.c
	banslv				imsl_d_spline_interp.c

!Revision History:
See top of file

Requirements:

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits
	"A Practical Guide to Splines - Revised edition", Carl de Boor,
	copyright 2001
	<http://netlib.bell-labs.com/netlib/master/readme.htm>

Design Notes:
	The is routine is a port to C from the FORTRAN routine named splint()
	described by de Boor (2001). That routine is provide in the splint.f
	file which is part of the pppack package which can be downloaded from
	the netlib ftp site listed above.

	This function's implementation differs from the netlib splint() by:
	Using the return value to indicate failure.

	This routine performs almost no validity checks on its arguments.
	Ensuring their validity is the responsibility of the calling routine.
	Therefore, this routine is static, so it can only be called from within
	the same file, making it easier to verify that calling routine has
	lived up to its responsibilities.
!END**************************************************************************
*/

{
    int i, ilp1mx, j, kpkm2, left;
    double taui;

    kpkm2 = 2*(k-1);
    left = k-1;

    for(i=0; i<n; i++)
	for(j=0; j<=kpkm2; j++)
	    q[i][j] = 0.0;

    /* Fill in the banded matrix with B-spline values at tau. */
    for(i=0; i<n; i++)
    {
	taui = tau[i];
	ilp1mx = (i+k<n) ? i+k-1 : n-1;
	left = (left>i) ? left : i;

	/* Find 'left' in the closed interval (i, ilp1mx-1) such that
	 *		t[left] <= tau[i] < t[left+1]
	 * matrix is singular if this is not possible.
	 */
	if(taui < t[left])
	    return IMSL_KNOT_DATA_INTERLACING;

	while(taui >= t[left+1] && left < ilp1mx)
	    left++;

	if(taui > t[left+1])
	    return IMSL_KNOT_DATA_INTERLACING;

	bsplvb(t, k, taui, left, bcoef);

	for(j=0; j<k; j++)
	    q[left-k+j+1][2*k-2+i-left-j] = bcoef[j];
    }

    banfac(q, n, k-1);
    memcpy(bcoef, gtau, n*sizeof(double));
    banslv(q, n, k-1, bcoef);

    return (Imsl_code)0;
}


Imsl_d_spline * IMSL_DECL imsl_d_spline_interp (
	int		ndata,
	double		xdata[],
	double		fdata[], 
	...
)
/*******************************************************************************
!C
!Description:   Computes a spline interpolant.

!Input Parameters:
	ndata		The number of data points.
	xdata		Array with ndata components containing the abscissas of
			the interpolation problem. Unlike the real IMSL
			function, this version requires that they be in
			ascending order.
	fdata		Array with ndata components containing the ordinates of
			the interpolation problem.
	...		Accepts a variable list of options, using the C
			<stdarg.h> interface. They come in sets, with the first
			argument of each set being an int that identifies the
			type of set. Multiple occurences of each set type are
			permitted; the later values replace the earlier ones.
			The sets of optional arguments permitted are:

	int IMSL_ORDER	Indicates that the next variable argument is:
	int order	The order of the spline subspace. Defaults to 4.

	int IMSL_KNOTS	Indicates that the next variable argument is:
	const double
	    knots[]	A list of ndata+order knots, provided by the user. 

	int 0		Terminates the variable arguments list.

!Output Parameters:
	None

Return Value:
	A pointer to a dynamically allocated structure that represents the
	spline interpolant. If an interpolant cannot be computed, then NULL is
	returned. The structure should be deallocated by calling free() when
	the caller is done with it.

Externally Defined:
	__STDC_VERSION__		set by compiler
	DBL_MAX				<float.h>
	IMSL_DECL			"imsl_wrap.h"
	IMSL_DUPLICATE_XDATA_VALUES	"imsl_wrap.h"
	IMSL_error_code			"pseudo_imsl.h"
	IMSL_KNOT_MULTIPLICITY		"imsl_wrap.h"
	IMSL_KNOT_NOT_INCREASING	"imsl_wrap.h"
	IMSL_KNOTS			"imsl_wrap.h"
	IMSL_ORDER			"imsl_wrap.h"
	IMSL_OUT_OF_MEMORY		"imsl_wrap.h"
	IMSL_OUT_OF_MEMORY_2		"imsl_wrap.h"
	IMSL_UNEXPECTED_NULL_POINTER	"imsl_wrap.h"
	IMSL_SPLINE_NEED_DATA_PTS	"imsl_wrap.h"
	IMSL_UNKNOWN_OPTION		"imsl_wrap.h"
	IMSL_XDATA_NOT_INCREASING	"imsl_wrap.h"
	IMSL_XDATA_TOO_LARGE		"imsl_wrap.h"
	IMSL_XDATA_TOO_SMALL		"imsl_wrap.h"
	MAX_IMPULSE_NUMBER		"GEO_geo.h"
	MAX_SPLINE_ORDER		"pseudo_imsl.h"

Called by:
	GEO_prepare_mirr_data		"GEO_input.h"

Routines Called:
	splint				imsl_d_spline_interp.c

!Revision History:
See top of file.

Requirements:

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits
	"IMSL C/Math/Library", copyright 1990-2001, Visual Numerics, Inc.

Design Notes:
	The interface is intended to mimic a subset of the capabilities of the 
	correspondingly named IMSL library function. That will allow us to
	make the use of the actual IMSL library optional, without requiring a
	multiple versions of the code that uses the function. The actual work
	is done by the splint() routine

	This function's interface differs from that of the real
	imsl_d_spline_interp in that:
	1. The xdata are required to be already sorted.
	2. No default value is provided for the knots. They must be provided by
	the user.
	3. ndata can't exceed MAX_IMPULSE_NUMBER, and order can't exceed
	MAX_SPLINE_ORDER.
	4. It doesn't try to de-reference null pointer arguments.
	5. Knots that are in decreasing order for some reason cause memory
	segment violations for the IMSL routine. That shouldn't happen here.
	(it shouldn't happen there, either).

!END**************************************************************************
*/

#if __STDC_VERSION__ < 199901L
    /* This definition makes the struct hack work on most C90 compilers, but
     * produces undefined behavior in C99.
     */
    #define HACK_LENGTH 1
#else
    /* This definition makes the struct hack work on all conforming C99
     * compilers, but will produce a syntax error in C90.
     */
    #define HACK_LENGTH
#endif
{
    typedef struct{
	Imsl_d_spline spline;
	int order, num_coef, num_knots;
	double *knots;
	double *coef;
	double hack[HACK_LENGTH];
    } spline_wrap;
    spline_wrap *pspline;

    va_list ap;
    int order=4;
    const double *knots=NULL;
    double w[MAX_IMPULSE_NUMBER][2*MAX_SPLINE_ORDER-1];
    double *q[MAX_IMPULSE_NUMBER];
    double epsilon;
    int i, arg;

    IMSL_error_code = (Imsl_code)0;

    va_start(ap, fdata);
    while((arg=va_arg(ap, int))) /* = is deliberate, not == */
    {
	if(arg==IMSL_ORDER)
	    order = va_arg(ap, int);
	else if(arg==IMSL_KNOTS)
	    knots = va_arg(ap, double*);
	else
	{
	    IMSL_error_code = IMSL_UNKNOWN_OPTION;	/* undocumented */
	    break;
	}
    }
    va_end(ap);

    /* Validity checks */

    if(IMSL_error_code)
	return NULL;

    if(order<1 || order>MAX_SPLINE_ORDER)
    {
	/* The upper limit is specific to pseudo-IMSL . */
	IMSL_error_code = IMSL_SPLINE_BAD_ORDER;	/* undocumented. */
	return NULL;
    }

    if(ndata < order)
    {
	IMSL_error_code = IMSL_SPLINE_NEED_DATA_PTS;	/* undocumented. */
	return NULL;
    }

    if(ndata > MAX_IMPULSE_NUMBER)
    {
	IMSL_error_code = IMSL_OUT_OF_MEMORY;	/* specific to pseudo-IMSL */
	return NULL;
    }

    /* Note: the real IMSL function computes default values for the knots, if
     * IMSL_KNOTS isn't specified by the caller. The fact that this version
     * doesn't do so makes IMSL_KNOTS mandatory.
     */

    if(xdata==NULL || fdata==NULL || knots==NULL)
    {
	/* The IMSL routine doesn't bother checking this. */
	IMSL_error_code = IMSL_UNEXPECTED_NULL_POINTER;
	return NULL;
    }

    /* Set up the spline structure. */
    /* The following technique is well known as the "struct hack". IMSL appears
     * to be using it in imsl_d_spline_interp().
     */
    pspline = (spline_wrap*) malloc(
	offsetof(spline_wrap,hack)+(2*ndata+order)*sizeof(double));
    if(!pspline)
    {
	IMSL_error_code = IMSL_OUT_OF_MEMORY_2;	/* undocumented. */
	return NULL;
    }

    pspline->spline.domain_dim = 1;
    pspline->spline.target_dim = 1;
    pspline->spline.order = &pspline->order;
    pspline->spline.num_coef = &pspline->num_coef;
    pspline->spline.num_knots = &pspline->num_knots;
    pspline->spline.knots = &pspline->knots;
    pspline->spline.coef = &pspline->coef;

    pspline->order = order;
    pspline->num_coef = ndata;
    pspline->num_knots = order+ndata;
    pspline->knots = pspline->hack;

    /* The next line had debatable legality in C90, but it works correctly on
     * every real compiler. It is fully legal in C99, so long as hack is
     * declared as "double hack[];", a declaration that would have been illegal
     * in C90.
     */
    pspline->coef = pspline->hack+ndata+order;
    memcpy(pspline->hack, knots, (ndata+order)*sizeof(double));

    epsilon = (knots[ndata+order-1]-knots[0]) / pow(DBL_MAX,1.0/(double)order);

    /* Calculations farther down will evaluate expressions whose value might
     * overflow, but which should always be smaller in magnitude than
     * pow((knots[ndata+order-1]-knots[0])/denominator,(double)order) for some
     * particular value of 'denominator'. The validity checks below that use
     * 'epsilon' are intended to ensure that 'denominator' is never small
     * enough to cause those expressions to overflow. The real IMSL routines
     * both do the equivalent of setting epsilon=0.0, which is a little riskier.
     */

    for(i=0; i<ndata-1; i++)
    {
	if(xdata[i+1] < xdata[i])
	{
	    /* This doesn't come up for the IMSL function; it sorts them. */
	    IMSL_error_code = IMSL_XDATA_NOT_INCREASING;
	    return &pspline->spline;
	}

	if(xdata[i+1] - xdata[i] < epsilon)
	{
	    IMSL_error_code = IMSL_DUPLICATE_XDATA_VALUES;
	    return &pspline->spline;
	}
    }

    for(i=0; i<ndata+order-1; i++)
	if(knots[i+1] < knots[i])
	{
	    /* This error code mnemonic is misnamed, but it is the documented
	     * one.
	     */
	    IMSL_error_code = IMSL_KNOT_NOT_INCREASING;
	    return &pspline->spline;
	}

    for(i=0; i<ndata; i++)
    {
	if(knots[i+order]-knots[i] <= epsilon)
	{
	    IMSL_error_code = IMSL_KNOT_MULTIPLICITY;
	    return &pspline->spline;
	}

	/* The INTERLACE checks in the main code, combined with the
	 * NOT_INCREASING checks above, are sufficient to cover the next two
	 * cases, so in principle they don't need to be done separately.
	 * However, the IMSL routine does the equivalent of these redundant
	 * checks, and therefore, so must this version.
	 */
	if(xdata[i] > knots[ndata])
	{
	    IMSL_error_code = IMSL_XDATA_TOO_LARGE;
	    return &pspline->spline;
	}

	if(xdata[i] < knots[order-1])
	{
	    IMSL_error_code = IMSL_XDATA_TOO_SMALL;
	    return &pspline->spline;
	}
    }

    for(i=0; i<ndata; i++)
	q[i] = w[i];

    IMSL_error_code = splint(xdata, fdata, knots, ndata, order, q,
	pspline->coef);

    return &pspline->spline;
}

