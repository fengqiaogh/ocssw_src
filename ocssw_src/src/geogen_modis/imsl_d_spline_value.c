/* Revision History:
$Log: imsl_d_spline_value.c,v $
Revision 6.1  2014/08/18 22:23:32  jkuyper
Corrected to zero-initial vector length.

Revision 4.4  2003/08/28 16:05:59  kuyper
Corrected prolog.

Revision 4.3  2003/08/26 21:23:33  kuyper
Corrected handling of optional arguments, to prevent use of uninitialized
  values.

Revision 4.2  2003/07/26 16:33:36  kuyper
Changed interv() to assume that calls come in non-decreasing order.
Changed imsl_d_spline_value() to validate that xvec is increasing.
Corrected sign error in bvalue().
Changed to return IMSL_SPLINE_ORDER_DERIV for bad values of deriv.
Corrected imsl_d_spline_value() to set IMSL_error_code after call to
  imsl_d_machine(), which re-sets the code.

Revision 4.1  2003/07/16 18:11:10  kuyper
Created by reverse engineering from documentation and testing of the real
IMSL library routine.

*/
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include "imsl_wrap.h"
#include "pseudo_imsl.h"

static int interv (
	const double    xt[],
	int             lxt,
	double          x,
	int             *left
){

/*******************************************************************************
!C
!Description:   
	Sets left to the maximum value for i such that
		xt[i] < xt[lxt-1] && xt[i] <= x

!Input Parameters:
	xt	A double array of length lxt, assumed to be non-decreasing.
	lxt	Number of elements in the array xt.
	x	The point whose location with respect to the sequence xt is to
		be determined.

!Output Parameters:
	None

!Input/Output Parameters:
	left	On input, points to an initial guess, assumed to be
		non-negative, for the final value to be written over that guess.
		This code assumes 0<= *left < lxt.
		On output, the int that it points at is set to:
		0 if x < xt[0]
		i if xt[i] <= x < xt[i+1]
		i if xt[i] < x < xt[i+1]  == xt[lxt1-1]
		i if xt[i] < xt[i+1] == xt[lxt-1] < x

	The asymmetric treatment of the intervals is due to the decision to
	make all pp functions continuous from the right.
		

Return Value:
	-1	if x < xt[0]
	1	if xt[lxt-1] < x	
	0	otherwise	(the usual case)

	By returning 0 even if x==xt[lxt-1], there is the option of having the
	computed pp function continuous from the left at xt[lxt].

Externally Defined:
	None

Called by:
	bvalue

Routines Called:
	None

!Revision History:
See top of file

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits
	This routine is based loosley upon the Fortran routine of the same name
	as described by Carl de Boor in "A Practical Guide to Splines - Revised
	Edition", 2001. The original routine can be found in interv.f, which is
	part of the pppack package which can be downloaded from
	<ftp://netlib.bell-labs.com/netlib/pppack.tar>

Design Notes
	The function is designed to be efficient in the common situation that
	it is called repeatedly, with x taken from an increasing sequence. This
	will happen, e.g., when a pp function is to be graphed.  The first guess
	pointed at by 'left' should  therefore normally be the value returned
	through that same parameter from the previous call. Then, if
	xt[ilo] <= x < xt[ilo+1], we set left = ilo and are done after just
	three comparisons.
	Otherwise, we repeatedly double the difference istep = ihi - ilo while
	also moving ihi in the direction of x, until
		xt[ilo] <= x < xt[ihi]
	after which we use bisection to get, in addition, ilo+1 = ihi.
	*left = ilo is then returned.

	This version differs from the original in that *left is used to
	initialize ilo, rather than making ilo static. This allows multiple
	different searches through different tables (or even the same table) to
	alternate with each other without interference. More importantly, it
	allows *left to be externally re-initialized when the table changes.
	Since it is externally initialized, *left is assumet to always be in a
	valid range for the current table. Finally, it is assumed that each
	table is accessed in increasing order.
!END**************************************************************************
*/
    int ilo = *left;
    int ihi=ilo+1, middle;

    while(ihi < lxt-1 && x >= xt[ihi])
	/* Increase ihi to capture x. */
	ihi = 2*ihi-ilo;

    if(ihi > lxt-1)
	ihi = lxt-1;

    while(ihi > ilo+1)
    {
	/* Narrow the interval. */
	middle = (ihi+ilo)/2;
	if(x > xt[middle])
	    ilo = middle;
	else
	    ihi = middle;
    }

    *left = ilo;
    /* Set return flag. */
    if(x < xt[0])
	return -1;

    if(x > xt[lxt-1])
	return 1;

    return 0;

}

static double bvalue(
	const double	t[/*n+k*/],
	const double	bcoef[/*n*/],
	int		n,
	int		k,
	double		x,
	int		jderiv,
	int		*pi
){
/*******************************************************************************
!C
!Description:   
	Calculates the value at x of the jderiv-th derivative of the spline
	described by t, bcoef, n, and k. The spline is taken to be continuous
	from the right, EXCEPT at the rightmost knot, where it is taken to be
	continuous from the left.

!Input Parameters:
	t	knot sequence, of length n+k, assumed nondecreasing
	bcoef	b-coefficient sequence, of length n.
	n	length of bcoef and dimension of spline(k,t), ASSUMED positive.
	k	order of the spline.
		WARNING: the restriction k < MAX_SPLINE_ORDER is imposed
		arbitrarily  by the declarations of aj, dl, dr below, but is
		NOWHERE CHECKED for.
	x	The point at which to evaluate.
	jderiv	Integer giving the order of the derivative to be evaluated.
		ASSUMED to be zero or positive.

!Output Parameters:
	None

!Input/Output Parameters:
	pi	On input: points at an estimate of the location for x in t[],
		usually the value returned by a previous call to this function.
		On output: The left side of the bracketing location for x, as
		returned by interv().

Return Value:
	The value of the (jderiv)-th derivative of f at x.

Externally Defined:
	MAX_SPLINE_ORDER		"pseudo_imsl.h"

Called by:
	imsl_d_spline_value

Routines Called:
	interv				imsl_d_spline_value.c

!Revision History:
See top of file.

Requirements:
	None

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits
	"A Practical Guide to Splines - Revised edition", Carl de Boor,
	copyright 2001
	<http://netlib.bell-labs.com/netlib/master/readme.htm>

Design Notes
	This is a port to C of the FORTRAN routine named bvalue() described in
	de Boor (2001). That routine is provided in the bvalue.f file which is
	part of the pppack package which can be downloaded from the netlib web
	site listed above.

	The nontrivial knot interval t[i], t[i+1] containing x is located with
	the aid of interv(). The k b-coeffs of f relevant for this interval are
	then obtained from bcoef (or taken to be zero if not explicitly
	available) and are then differenced jderiv times to obtain the b-coeffs
	of (d**jderiv)f relevant for that interval. Precisely, with j=jderiv,
	we have from X.12 of de Boors (2001) that

		(d**j)f = sum ( bcoef[.,j]*b[.,k-j,t] )

	where

			    | bcoef[.]			     j==0
			    |
		bcoef[.,j] =| bcoef[.,j-1] - bcoef[.-1,j-1]
			    | -----------------------------, j > 0
			    |    (t[.+k-j] - t[.])/(k-j)

	Then, we use repeatedly the fact that

	sum( a[.]*b[.,m,t](x) )  = sum(a[.,x]*b[.,m-1,t](x) )

	with

		 (x - t[.])*a[.] + (t[.+m-1] - x)*a[.-1]
	a[.,x] = ---------------------------------------
		 (x - t[.])      + (t[.+m-1] - x)

	to write (d**j)f(x) eventually as a linear combination of b-splines of
	order 1, and the coefficient for b[i,1,t](x) must then be the desired
	number (d**j)f(x). (See X.17-X.19 of de Boors (2001)).

	This function makes many unchecked assumptions, including both the ones
	explicitly mentioned above, and the assumption that its pointer
	arguments are valid. The validity of those assumptions is the
	responsibility of the calling function. Therefore, this function has
	been made static, so it can only be called by other functions in the
	same file. This will make it easier to verify that the the calling
	function fulfills its responsibility.
!END**************************************************************************
*/
    double aj[MAX_SPLINE_ORDER], dl[MAX_SPLINE_ORDER], dr[MAX_SPLINE_ORDER];
    double fkmj;
    int i, ilo, imk, j, jc, jcmin, jcmax, jj, kmj, km1, nmi;

    if(jderiv >= k)
	return 0.0;

    /* Find i s.t. 0 <= i < n+k-1 AND t[i] < t[i+1] and t[i] <= x < t[i+1]. If
     * no such i can be found, x lies outside the support of the spline f,
     * hence bvalue = 0.
     * The asymmetry in this choice of i makes f rightcontinuous, except at
     * t[n+k-1] where it is leftcontinuous.)
     */
    if(interv(t, n+k, x, pi)!=0)
	return 0.0;

    i = *pi;
    km1 = k - 1;
    if(km1 <= 0)
	return bcoef[i];

    /* Store the k b-spline coefficients relevant for the knot interval
     * (t[i], t[i+1]) in aj[0],...,aj[k-1] and compute dl[j] = x - t[i+1-j],
     * from input to zero. Set any t's not obtainable equal to t[0] or to t[n+k]
     * appropriately.
     */
    jcmin = 0;
    imk = i - k + 1;
    if(imk < 0)
    {
	jcmin = -imk;
	for(j=0; j<=i; j++)
	    dl[j] = x - t[i-j];
	for(j=i; j<km1; j++)
	{
	    aj[k-j-2] = 0.0;
	    dl[j] = dl[i];
	}
    }
    else
	for(j=0; j<km1; j++)
	    dl[j] = x - t[i-j];

    jcmax = k;
    nmi = n - i;
    if(nmi < 1)
    {
	jcmax = k + nmi - 1;
	for(j=0; j<jcmax; j++)
	    dr[j] = t[i+j+1] - x;
	for(j=jcmax-1; j< km1; j++)
	{
	    aj[j+1] = 0.0;
	    dr[j] = dr[jcmax-1];
	}
    }
    else
	for(j=0; j<km1; j++)
	    dr[j] = t[i+j+1] - x;

    for(jc=jcmin; jc<jcmax; jc++)
	aj[jc] = bcoef[imk+jc];

    /* Difference the coefficients jderiv times */
    for(j=1; j<=jderiv; j++)
    {
	kmj = k - j;
	fkmj = (double)kmj;
	ilo = kmj;
	for(jj=1; jj<=kmj; jj++)
	    aj[jj-1] = ((aj[jj] - aj[jj-1])/(dl[--ilo] + dr[jj-1]))*fkmj;
    }

    /* Compute value at x in (t[i],t[i+1]) of jderiv-th derivative,
     * give its relevant b-spline coeffs in aj[0],...,aj[k-jderiv].
     */
    for(j=jderiv; j<km1; j++)	/* Double check these limits. */
    {
	kmj = k - j - 1;
	ilo = kmj;
	for(jj=0; jj<kmj; jj++) /* Double check limits. */
	{
	    --ilo;
	    aj[jj] = (aj[jj+1]*dl[ilo] + aj[jj]*dr[jj])/(dl[ilo] + dr[jj]);
	}

    }

    return aj[0];
}

double IMSL_DECL imsl_d_machine(
	int n
){
/*******************************************************************************
!C
!Description:   
	Returns information describing the computer's double precision
	floating-point arithmetic.

!Input Parameters:
	n	Index indicating which value is to be returned. The index must
		be from 1 to 8, inclusive.

!Output Parameters:
	None

Return Value:
	n	Value returned
	----	-----------------
	1	The smallest positive number.
	2	The largest number
	3	The smallest relative spacing.
	4	The largest relative spacing.
	5	log10(floating point radix)
	6	NaN (not a number)
	7	positive machine infinity
	8	negative machine infinity
	other	NaN

	On machines that do not support NaN, imsl_d_machine(6) returns a value
	larger than imsl_d_machine(2). On machines which do not have a special
	representation for infinity, imsl_d_machine(2)==imsl_d_machine(7).

Externally Defined:
	__STDC_VERSION__		Defined by compiler (C94 or later)
	DBL_EPSILON			<float.h>
	DBL_MAX				<float.h>
	DBL_MIN				<float.h>
	FLT_RADIX			<float.h>
	HUGE_VAL			<math.h>
	Imsl_code			"imsl_wrap.h"
	IMSL_DECL			"imsl_wrap.h"
	IMSL_error_code			"pseudo_imsl.h"
	IMSL_INTEGER_OUT_OF_RANGE	"imsl_wrap.h"

Called by:
	GEO_prepare_mirr_data
	GEO_interp_mirr_enc
	imsl_d_spline_value

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
	Reverse engineered from:
	"IMSL C/Math/Library", copyright Visual Numerics, Inc. 1999-2001

Design Notes
!END**************************************************************************
*/

/* The nan() function was introduced in C99. It's defined to return a quiet NaN
 * if the implementation supports quiet NaNs, and to otherwise return 0.0.
 */
#if __STDC_VERSION__ < 199901L
#define IMSL_NAN HUGE_VAL
#else
#define IMSL_NAN (nan("")==0.0 ? HUGE_VAL : nan(""))
#endif

    IMSL_error_code = (Imsl_code)0;

    switch(n)
    {
	case 1: return DBL_MIN;

	case 2:
	    if(IMSL_NAN!=HUGE_VAL || DBL_MAX < HUGE_VAL)
		return DBL_MAX;
#if __STDC_VERSION__ < 199901L
	    /* Note: for C90, __STDC_VERSION__ is not defined, which counts as
	     * a 0 in #if statements, which achieves exactly the desired result.
	     */
	    return DBL_MAX*(1.0-DBL_EPSILON);
#else
	    return nextafter(DBL_MAX, 0.0);
#endif

	case 3: return DBL_EPSILON/FLT_RADIX;
	case 4: return DBL_EPSILON;
	case 5: return log10(FLT_RADIX);

	case 7:
	    if(IMSL_NAN!=HUGE_VAL || DBL_MAX < HUGE_VAL)
		return HUGE_VAL;
#if __STDC_VERSION__ < 199901L
	    return HUGE_VAL*(1.0-DBL_EPSILON);
#else
	    return nextafter(HUGE_VAL, 0.0);
#endif

	case 8: return -HUGE_VAL;

	default:	
	    IMSL_error_code = IMSL_INTEGER_OUT_OF_RANGE; /* undocumented */
	    /* FALL through to next case */
	case 6:
	    return IMSL_NAN;
    }
}

double IMSL_DECL imsl_d_spline_value(
	double		x,
	Imsl_d_spline	*sp,
	...
){
/*******************************************************************************
!C
!Description:   
	Calculates the value of a spline or the value of one of its derivatives.

!Input Parameters:
	x	Evaluation point for the spline. Ignored if the IMSL_GRID_USER
		option is selected (undocumented). 
	sp	Pointer to the structure that represents the spline.
	...	One or more sets of optional arguments are accepted, using the
		C <stdarg.h> interface, each starting with an int argument that
		identifies the type of the set. The permitted sets of optional
		arguments are:

	int	IMSL_DERIV	Indicates that the following argument is:
	int	deriv		The derivative of the spline to be computed.
				default = 0.

	int	0		Indicates the end of the list of arguments

!Output Parameters:
	None.

!Input/Output Parameters: 
	int	IMSL_GRID_USER	Indicates that the next three arguments are:
	int	n		The length of the xvec array
	const double xvec[]	The evaluation points for the spline
	double	value_user[]	The location where the values are to be written

Return Value:
	The value of the spline or one of its derivatives at the point x.
	If no value can be computed, imsl_d_machine(6) is returned instead.
	That same value is also returned if the the IMSL_GRID_USER option is
	chosen (undocumented).
	imsl_d_machine(6) is a NaN, if possible. If not, it's a normally valid
	value that the IMSL functions treat as if it were a true NaN.

Externally Defined:
	DBL_MAX				<float.h>
	Imsl_code			"imsl_wrap.h"
	IMSL_DECL			"imsl_wrap.h"
	IMSL_DERIV			"imsl_wrap.h"
	IMSL_error_code			"pseudo_imsl.h"
	IMSL_GRID_USER			"imsl_wrap.h"
	IMSL_KNOT_MULTIPLICITY		"imsl_wrap.h"
	IMSL_KNOT_NOT_INCREASING	"imsl_wrap.h"
	IMSL_SPLINE_BAD_COEFFS		"imsl_wrap.h"
	IMSL_SPLINE_BAD_ORDER		"imsl_wrap.h"
	IMSL_SPLINE_ORDER_DERIV		"imsl_wrap.h"
	IMSL_UNEXPECTED_NULL_POINTER	"imsl_wrap.h"
	IMSL_UNKNOWN_OPTION		"imsl_wrap.h"
	IMSL_XVEC_LENGTH		"imsl_wrap.h"
	IMSL_XVEC_NOT_INCREASING	"imsl_wrap.h"
	MAX_SPLINE_ORDER		"pseudo_imsl.h"

Called by:
	GEO_interp_mirr_enc
	GEO_prepare_mirr_data

Routines Called:
	bvalue				imsl_d_spline_value.c
	imsl_d_machine			"imsl_wrap.h"

!Revision History:
See top of file

Requirements:
	None

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits
	"IMSL C/Math/Library", copyright 1990-2001.
	"A Practical Guide to Splines - Revised edition", Carl de Boor,
	copyright 2001

Design Notes
	This function is a wrapper for a modified C version of the Fortran
	bvalue() function described in de Boor (2001).

	The interface is intended to mimic a subset of the features of the
	correspondingly named IMSL library function. That will allow us to make
	the use of the actual IMSL library optional, without requiring a
	re-write of the code that uses the function. 

	This function's interface differs from that of the real
	imsl_d_spline_value() in that the IMSL_GRID option is not supported.
!END**************************************************************************
*/
    va_list ap;	/* see <stdarg.h> */
    const double *xvec=NULL;
    double *value_user=NULL, epsilon, imsl_nan=imsl_d_machine(6);
    int i=0, grid=0, deriv=0;
    int arg, knot, n=0, num_knots, val;

    IMSL_error_code = (Imsl_code)0;
    va_start(ap, sp);

    while(( arg=va_arg(ap,int) )) /* = rather than == is deliberate */
    {
	if(arg == IMSL_DERIV)
	    deriv = va_arg(ap,int);
	else if(arg == IMSL_GRID_USER)
	{
	    n = va_arg(ap,int);
	    xvec = (const double *)va_arg(ap,double *);
	    value_user = va_arg(ap, double*);
	    grid = 1;
	}
	else
	{
	    IMSL_error_code = IMSL_UNKNOWN_OPTION;	/* undocumented. */
	    return imsl_nan;
	}
    }
    va_end(ap);

    if(sp==NULL || sp->order==NULL || sp->num_coef==NULL || sp->knots==NULL ||
	sp->coef==NULL || sp->knots[0]==NULL || sp->coef[0] == NULL)
    {
	/* IMSL doesn't actually test this. */
	IMSL_error_code = IMSL_UNEXPECTED_NULL_POINTER;
	return imsl_nan;
    }

    /* The real IMSL function does not verify that sp->domain_dim==1 or
     * sp->target_dim==1, therefore we don't have to.
     */

    if(sp->order[0] < 1 || sp->order[0] > MAX_SPLINE_ORDER)
    {
	/* The upper limit is specific to pseudo-IMSL. */
	IMSL_error_code = IMSL_SPLINE_BAD_ORDER;	/* undocumented. */
	return imsl_nan;
    }

    if(sp->num_coef[0] < sp->order[0])
    {
	IMSL_error_code = IMSL_SPLINE_BAD_COEFFS;	/* undocumented. */
	return imsl_nan;
    }

    if(deriv < 0)
    {
	IMSL_error_code = IMSL_SPLINE_ORDER_DERIV;	/* undocumented. */
	return imsl_nan;
    }

    num_knots = sp->num_coef[0] + sp->order[0];
    /* The above calculation removes the need to validate that
     *	sp->num_knots[0] == sp->num_coef[0]+sp->order[0]
     * a validity test which the real IMSL function doesn't make.
     */

    /* The real IMSL function uses the equivalent of epsilon=0.0, in order to
     * avoid division by zero. Setting epsilon to this particular value allows
     * us to also avoid division by numbers that are so small that the division
     * could overflow.
     */
    epsilon = pow((sp->knots[0][num_knots-1] - sp->knots[0][0])/DBL_MAX,
	1.0/(double)sp->order[0]);

    /* Validate that knots are non-decreasing. */
    for(knot=0; knot<num_knots-1; knot++)
    {
	if(sp->knots[0][knot] > sp->knots[0][knot+1])
	{
	    /* This mnemonic is mis-named. */
	    IMSL_error_code = IMSL_KNOT_NOT_INCREASING;
	    return imsl_nan;
	}
    }

    /* Validate that there are never more than sp->order[0] knots at (almost)
     * the same location.
     */
    for(knot=0; knot<sp->num_coef[0]; knot++)
    {
	if(sp->knots[0][knot+sp->order[0]] - sp->knots[0][knot] <= epsilon)
	{
	    IMSL_error_code = IMSL_KNOT_MULTIPLICITY;
	    return imsl_nan;
	}
    }

    if(grid==0)
    {
    	/* Only need to evaluate for x. */
	return bvalue(sp->knots[0], sp->coef[0], sp->num_coef[0], sp->order[0],
	    x, deriv, &i);
    }

    if(xvec==NULL || value_user==NULL)
    {
	/* IMSL doesn't actually make this test. */
	IMSL_error_code = IMSL_UNEXPECTED_NULL_POINTER;
	return imsl_nan;
    }

    if(n < 1)
    {
	IMSL_error_code = IMSL_XVEC_LENGTH;	/* Undocumented. */
	return imsl_nan;
    }

    for(val=0; val<n; val++)
    {
	if(val && xvec[val] <= xvec[val-1])
	{
	    IMSL_error_code = IMSL_XVEC_NOT_INCREASING;	/* Undocumented */
	    break;
	}

	value_user[val] = bvalue(sp->knots[0], sp->coef[0], sp->num_coef[0],
	    sp->order[0], xvec[val], deriv, &i);
    }

    return imsl_nan;
}

