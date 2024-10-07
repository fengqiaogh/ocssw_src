/*
 *----------------------------------------------------------------------
 * @(#) tconv.c		1.0	20 Oct 93	<shc>
 * Copyright (c) 1993, CSIRO Division of Oceanography
 *----------------------------------------------------------------------
 *
 * tconv --
 *
 *	Time conversions between UT, UT1 and TDT, which is required for
 *	calculation of solar/lunar positions and precession/nutation
 *	conversions.  Rigorously, Barycentric Dynamical Time (TDB) should be
 *	used, but it differs from TDT by at most 2msec.
 *
 *	The difference between UTC and TDT can be specified accurately using
 *	the initialisation routine TCSET.  If no difference is specified, a
 *	polynomial approximation is used, but any resulting errors in
 *	miscalculation of solar/lunar/precession/nutation effects are
 *	probably in the noise.  The approximation is derived from a cubic fit
 *	to the UT/TDT difference from 1950 to 1992 plus some speculative
 *	points for 1993, 2000 and 2010.  A cubic, rather than a higher order
 *	polynomial, is used so the divergence isn't too wild past 2010.
 *
 * Results:
 *	Given a time in modified Julian days and a conversion code string
 *	of the form "<in_type>:<out_type>", return the appropriate time.
 *	On error, ERRSTR is set to point to an error message and FP_ERRVAL
 *	is returned.
 *
 * Side effects:
 *	None.
 *
 * References:
 *	"Astronomical Algorithms", Jean Meeus, Willmann-Bell Inc, 1991.
 *
 * History:
 *	1.0	20 Oct 93	<shc>
 *	   Converted from FORTRAN to C.
 *
 *	1.0	12 Dec 94	<shc>
 *	   Added errmsg stuff for integration with Python
 *
 *----------------------------------------------------------------------
 */

#include <stdio.h>
#include <string.h>
#include "orbit.h"


/* User specified corrections (if any) */
static double dtmutc = 0, /* TDT-UTC */
        ut1mtc = 0; /* UT1-UTC */

double
tconv(tjd, conv, errmsg)
double tjd;
char *conv;
char **errmsg;
{
    double c, u, tdt;

    if (strlen(conv) != 7 || conv[3] != ':') {
        if (errmsg)
            *errmsg = "conversion spec should be <inp_type>:<out_type>";
        return FP_ERRVAL;
    }

    /* Get or calculate a rough value for TDT-UTC in days */
    if (dtmutc != 0)
        u = dtmutc;
    else {
        c = (tjd - J2000) / 36525.0;
        u = ((-44.177 * c + 26.763) * c + 104.011) * c + 68.215;
        u /= 86400.0;
    }

    /* Convert input time to TDT */
    if (conv[0] == 'u' && conv[1] == 't' && conv[2] == 'c')
        tdt = tjd + u;
    else if (conv[0] == 'u' && conv[1] == 't' && conv[2] == '1')
        tdt = tjd - ut1mtc + u;
    else if (conv[0] == 't' && conv[1] == 'd' && conv[2] == 't')
        tdt = tjd;
    else {
        if (errmsg)
            *errmsg = "input type should be utc | ut1 | tdt";
        return FP_ERRVAL;
    }


    /* Convert TDT to output time */
    if (conv[4] == 'u' && conv[5] == 't' && conv[6] == 'c')
        return tdt - u;
    else if (conv[4] == 'u' && conv[5] == 't' && conv[6] == '1')
        return tdt - u + ut1mtc;
    else if (conv[4] == 't' && conv[5] == 'd' && conv[6] == 't')
        return tdt;
    else {
        if (errmsg)
            *errmsg = "output type should be utc | ut1 | tdt";
        return FP_ERRVAL;
    }
}

/*
 * This routine allows specification of the value of TDT-UT in seconds.
 */
int
tcset(x, what, errmsg)
double x;
char *what;
char **errmsg;
{
    if (strcmp(what, "tdt-utc") == 0)
        dtmutc = x / 86400.0;
    if (strcmp(what, "ut1-utc") == 0)
        ut1mtc = x / 86400.0;
    else {
        if (errmsg)
            *errmsg = "identifier should be tdt-utc | ut1-utc";
        return -1;
    }

    return 0;
}
