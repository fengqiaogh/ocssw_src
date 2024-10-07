/*
 *----------------------------------------------------------------------
 * @(#) obliq.c		1.0	20 Oct 93	<shc>
 * Copyright (c) 1993, CSIRO Division of Oceanography
 *----------------------------------------------------------------------
 *
 * obliq --
 *
 *	Calculate mean obliquity of date.
 *
 * Results:
 *
 *	Given a modified Julian dynamical time TJD, this function
 *	returns	the mean obliquity of date of the ecliptic and its
 *	time derivative in radians per second.  Each value pointer
 *	can be NULL if the particular value is not required.
 *
 * Side effects:
 *	None.
 *
 * References:
 *	"Astronomical Algorithms", Jean Meeus, Willmann-Bell Inc, 1991,
 *	Chapter 45.
 *
 * History:
 *	1.0	20 Oct 93	<shc>
 *	   Initial version.
 *
 *----------------------------------------------------------------------
 */

#include "orbit.h"

void
obliq(tjd, mood, dmood)
double tjd;
double *mood, *dmood;
{
    double t, a0, a1;

    /* Julian centuries from 2000 Jan 1.5 TDB */
    t = (tjd - J2000) / 36525.0;

    if (mood != (double *) NULL) {
        a0 = ((1.813e-3 * t - 5.9e-4) * t - 46.8150) * t + 84381.4119;
        *mood = DTOR(a0 / 3600.0);
    }

    if (dmood != (double *) NULL) {
        a1 = (3 * 1.813e-3 * t - 2 * 5.9e-4) * t - 46.8150;
        *dmood = DTOR(a1 / 3600.0) / 36525.0;
    }

    return;
}
