/*
 *----------------------------------------------------------------------
 *  @(#) ctogd.c	1.0	30 Jan 95	<shc>
 *  Copyright (c) 1995, CSIRO Division of Oceanography.
 *----------------------------------------------------------------------
 *
 * ctogd --
 *
 *	Rectangular (Cartesian) to geodetic coordinate conversion.
 *
 * Results:
 *
 *	Given the rectangular (Cartesian) coordinates of an object,
 *	calculate its geodetic latitude, longitude, height and, optionally,
 *	the derivatives of these wrt rectangular.  The GHA argument allows
 *	the Greenwich Hour Angle to be taken into account if desired.
 *	Algorithm and formulations are from GTDS manual, section 3.3.6.3.
 *	Yet again, GTDS is incorrect in its formulations of the partial
 *	derivatives.  The ones here have been derived by analytically
 *	inverting the Jacobian for the geodetic to Cartesian conversion.
 *
 *	Units are radians, metres and seconds.
 *
 *	Lat, long and height are returned in GEOD.  If "DGEOD" is
 *	non-NULL, the 3x3 Jacobian of geodetic coordinates wrt
 *	rectangular coordinates is returned in DGEOD.
 *
 * Side effects:
 *	None.
 *
 * History:
 *	30 Jan 95  <shc>
 *	   Converted from FORTRAN to C
 *
 *----------------------------------------------------------------------
 */

#include <math.h>

#include "orbit.h"
#include "earth.h"

#define EPS 1.0e-14  /* Convergence tolerance */

void
ctogd(r, gha, geod, dgeod)
double r[3]; /* Rectangular coordinates of object (metres) */
double gha; /* Greenwich Hour Angle correction (radians) */
struct GEODETIC *geod; /* Geodetic coordinates of object (radians, metres) */
double dgeod[3][3]; /* Jacobian of (lat, lon, hgt) wrt (x, y, z) */
{

    double esq, eqrsq, zi, ziold, zib, enph, sinlat, en, h, eqr,
            omesq, a, b;

    esq = FLT * (2.0 - FLT);
    eqrsq = r[0] * r[0] + r[1] * r[1];
    zi = -esq * r[2];

    do {
        ziold = zi;
        zib = r[2] - ziold;
        enph = sqrt(eqrsq + zib * zib);
        sinlat = zib / enph;
        en = EQRAD / sqrt(1.0 - esq * sinlat * sinlat);
        zi = -en * esq * sinlat;
    } while (fabs(zi - ziold) > EPS * fabs(zi));

    /*
     * Spacecraft latitude, longitude and height.
     */

    geod->lat = asin(sinlat);
    geod->lon = ANPI(atan2(r[1], r[0]) - gha);
    geod->height = h = enph - en;

    /*
     * Derivatives of latitude, longitude and height wrt Cartesian
     * coordinates.  First row is derivatives of latitude wrt
     * x (dgeod[0][0])), y (dgeod[0][1]) and z (dgeod[0][2]).  Second
     * row is derivatives of longitude, third is derivatives of height.
     */

    if (dgeod != NULL) {

        eqr = sqrt(eqrsq);
        omesq = (1.0 - esq);
        a = h + en * omesq / (1.0 - esq * sinlat * sinlat);
        b = sinlat / (a * eqr);
        dgeod[0][0] = -b * r[0];
        dgeod[0][1] = -b * r[1];
        dgeod[0][2] = eqr / (a * enph);

        dgeod[1][0] = -r[1] / eqrsq;
        dgeod[1][1] = r[0] / eqrsq;
        dgeod[1][2] = 0.0;

        dgeod[2][0] = r[0] / enph;
        dgeod[2][1] = r[1] / enph;
        dgeod[2][2] = r[2] / enph;

    }
}
