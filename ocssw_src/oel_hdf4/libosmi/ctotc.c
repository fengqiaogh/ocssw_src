/*
 *----------------------------------------------------------------------
 *  @(#) ctotc.c	1.0	20 Mar 93	<shc>
 *  Copyright (c) 1993, CSIRO Division of Oceanography.
 *----------------------------------------------------------------------
 *
 * ctotc --
 *
 *	Rectangular (Cartesian) to topocentric coordinate conversion.
 *
 * Results:
 *
 *	Given the rectangular (Cartesian) earth-centred coordinates of an
 *	object, the latitude, longitude and height above the geoid of
 *	an observer and the Greenwich Hour angle, calculate the topocentric
 *	coordinates of the object (azimuth, elevation and range) and,
 *	optionally, derivatives of these wrt Cartesian.
 *	Units are radians, metres and seconds.
 *
 *	Azimuth, elevation and range are returned in TOP.  If "DAER"
 *	is non-NULL, the 3x3 Jacobian of topocentric coordinates wrt
 *	rectangular coordinates is returned in DAER.
 *
 * Side effects:
 *	None.
 *
 * History:
 *	17 Dec 90  <shc>
 *	   Initial version
 *
 *	23 Jan 91  <shc>
 *	   Reduced number of trig routine calls
 *
 *	20 Mar 93  <shc>
 *	   Converted from FORTRAN to C
 *
 *----------------------------------------------------------------------
 */

#include "orbit.h"

void
ctotc(r, obs, gha, top, daer)
double r[3]; /* Rectangular coordinates of object (metres) */
struct GEODETIC *obs; /* Observer's lat, lon (radians) and height (m) */
double gha; /* Greenwich hour angle (radians) */
struct TOPOCENTRIC *top; /* Resulting az, el and range */
double daer[3][3]; /* Jacobian of (az,el,range) wrt (x,y,z) */
{
    int i, j;
    double olat, olon, ohgt, clat, clon, slat, slon;
    double ens, xo, yo, zo, xos, yos, zos, xlt, ylt, zlt;
    double rngsq, range;
    double emlt[3][3], dlt[3][3];
    /*
     * Observer geodetic coordinates (latitude, longitude and height) in
     * radians and metres.
     */
    olat = obs->lat;
    olon = AN2PI(obs->lon + gha);
    ohgt = obs->height;

    clat = cos(olat);
    clon = cos(olon);
    slat = sin(olat);
    slon = sin(olon);
    /*
     * Rotation matrix to bring geocentric reference frame
     * parallel to local tangent reference frame.
     */
    emlt[0][0] = -slon;
    emlt[0][1] = clon;
    emlt[0][2] = 0.0;
    emlt[1][0] = -slat*clon;
    emlt[1][1] = -slat*slon;
    emlt[1][2] = clat;
    emlt[2][0] = clat*clon;
    emlt[2][1] = clat*slon;
    emlt[2][2] = slat;
    /*
     * Geocentric vector to observer
     */
    ens = EQRAD / sqrt(1.0 - FLT * (2.0 - FLT) * slat * slat);
    xo = (ens + ohgt) * clat * clon;
    yo = (ens + ohgt) * clat * slon;
    zo = (ens * (1.0 - FLT)*(1.0 - FLT) + ohgt) * slat;
    /*
     * Vector from observer to object
     */
    xos = r[0] - xo;
    yos = r[1] - yo;
    zos = r[2] - zo;
    /*
     * Position of object in local tangent reference frame
     */
    xlt = emlt[0][0] * xos + emlt[0][1] * yos + emlt[0][2] * zos;
    ylt = emlt[1][0] * xos + emlt[1][1] * yos + emlt[1][2] * zos;
    zlt = emlt[2][0] * xos + emlt[2][1] * yos + emlt[2][2] * zos;
    /*
     * Spacecraft azimuth, elevation and range
     */
    rngsq = xlt * xlt + ylt * ylt + zlt*zlt;
    range = sqrt(rngsq);
    top->az = AN2PI(atan2(xlt, ylt));
    top->el = asin(zlt / range);
    top->range = range;
    /*
     * Optionally calculate derivatives of azimuth, elevation and range
     * with respect to the object rectangular coordinates.
     */
    if (daer != NULL) {
        double a1, a2;
        /*
         * Derivatives of azimuth (DLT[0][*]), elevation (DLT[1][*]) and
         * range (DLT[2][*]) wrt object position in local tangent frame.
         */
        a1 = xlt * xlt + ylt*ylt;
        a2 = sqrt(a1);
        dlt[0][0] = ylt / a1;
        dlt[0][1] = -xlt / a1;
        dlt[0][2] = 0.0;
        dlt[1][0] = -zlt * xlt / (rngsq * a2);
        dlt[1][1] = -zlt * ylt / (rngsq * a2);
        dlt[1][2] = a2 / rngsq;
        dlt[2][0] = xlt / range;
        dlt[2][1] = ylt / range;
        dlt[2][2] = zlt / range;
        /*
         * Convert derivatives wrt local tangent frame to derivatives wrt
         * geocentric reference frame: Product[dlt, emlt].
         */
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                daer[i][j] = dlt[i][0] * emlt[0][j] +
                        dlt[i][1] * emlt[1][j] +
                        dlt[i][2] * emlt[2][j];
            }
        }
    }

    return;
}
