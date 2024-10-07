/*
 *----------------------------------------------------------------------
 * @(#) gmha.c		1.1	20 Mar 93	<shc>
 * Copyright (c) 1993, CSIRO Division of Oceanography
 *----------------------------------------------------------------------
 *
 * gmha --
 *
 *	Greenwich Mean Hour angle.
 *
 * Results:
 *
 *	gmha() returns the Greenwich mean hour angle at a given time,
 *	which must be in modified Julian date format (number of days since
 *	noon on the 31st Dec 4714 minus 2,400,000).  It does not exhibit
 *	the slight day-to-day discontinuity of algorithms based on day and
 *	fraction.  Lifted and adapted from Chris Rizos's code which was
 *	derived from a U.S. Naval Observatory circular.  dgmha() returns
 *	the derivative of mean hour angle in radians/sec.
 *
 * Side effects:
 *	None.
 *
 * History:
 *	1.0	28 Nov 90	<shc>
 *	   Initial FORTRAN version.
 *
 *	1.1	14 May 92	<shc>
 *	   Changed gha() to gmha().
 *
 *	1.1	20 Mar 93	<shc>
 *	   Converted from FORTRAN to C.
 *
 *	1.2	27 Jun 95	<shc>
 *	   Added dgmha().
 *
 *----------------------------------------------------------------------
 */

#include <math.h>
#include "orbit.h"


/*
 * Polynomial coefficients for GMST from USNO circular #163, P.A3,
 * converted to radians and divided by 36525.0**N (N=0..3)
 */

static double gmstc[] = {
    0.4894961212823058751375704430e+01,
    0.6300388098984893552276513720e+01,
    0.5075209994113591478053805523e-14,
    -0.9253097568194335640067190688e-23
};

/*
 * Greenwich mean hour angle (radians)
 */

double
gmha(tjd)
double tjd; /* Time in Julian days mod 2400000 */
{
    double t, sid, gmha;

    t = tjd - J2000;
    sid = gmstc[0] + t * (gmstc[1] + t * (gmstc[2] + t * gmstc[3]));
    gmha = fmod(fmod(sid, D2PI) + D2PI, D2PI);

    return (gmha);
}

/*
 * Greenwich mean hour angle derivative (radians/day)
 */

double
gmhadot(tjd)
double tjd; /* Time in Julian days mod 2400000 */
{
    double t, dsid;

    t = tjd - J2000;
    dsid = gmstc[1] + t * (2.0 * gmstc[2] + t * 3.0 * gmstc[3]);

    return (dsid);
}
