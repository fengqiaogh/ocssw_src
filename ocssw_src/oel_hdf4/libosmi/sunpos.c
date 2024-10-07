/*
 *----------------------------------------------------------------------
 * @(#) sunpos.c	1.0	22 Oct 93	<shc>
 * Copyright (c) 1993, CSIRO Division of Oceanography
 *----------------------------------------------------------------------
 *
 * sunpos --
 *
 *	Position of the Sun.
 *
 * Results:
 *
 *	Given a modified Julian dynamical time TJD, this function
 *	set the elements of R to the mean-of-date rectangular coordinates
 *	(in metres) of the sun and returns 0.  Otherwise, ERRSTR is set to
 *	point to an error message and -1 is returned.
 *
 * Side effects:
 *	None.
 *
 * References:
 *	"Astronomical Formulae for Calculators", Jean Meeus.
 *
 * History:
 *	1.0	22 Oct 93	<shc>
 *	   Initial version.
 *
 *	1.0	22 Oct 93	<shc>
 *	   Added errmsg stuff for integration with Python
 *
 *----------------------------------------------------------------------
 */

#include "orbit.h"

/* Prototypes of internal functions */
static void oearth(double tjd, struct KEPLER *kep);
static void cearth(double tjd, double *el, double *er);

int
sunpos(tjd, r, errmsg)
double tjd;
double r[3];
char **errmsg;
{
    double ea, ta, helon, rad, x, y, z, obl, cobl, sobl;
    struct KEPLER kep;

    /* Get equinox-of-date orbital elements for the earth */
    (void) oearth(tjd, &kep);

    /* Calculate eccentric anomaly and true anomaly */
    ea = eanom(kep.man, kep.ecc, errmsg);
    if (ea == FP_ERRVAL)
        return -1;
    ta = 2 * atan(sqrt((1 + kep.ecc) / (1 - kep.ecc)) * tan(ea / 2));

    /* Calculate heliocentric ecliptic longitude and radius */
    helon = AN2PI(ta + kep.arp);
    rad = kep.sma * (1 - kep.ecc * kep.ecc) / (1 + kep.ecc * cos(ta));

    /* Add perturbations */
    (void) cearth(tjd, &helon, &rad);

    /* Calculate heliocentric ecliptic cartesian coordinates */
    x = rad * cos(helon);
    y = rad * sin(helon);
    z = 0;

    /*
     * Rotate around X-axis through obliquity of date to convert to
     * heliocentric equatorial coordinates.  Take negatives to get
     * coordinates of sun relative to earth and convert from AUs
     * to metres.
     */
    (void) obliq(tjd, &obl, (double *) NULL);
    cobl = cos(obl);
    sobl = sin(obl);

    r[0] = -AU * x;
    r[1] = -AU * (cobl * y - sobl * z);
    r[2] = -AU * (sobl * y + cobl * z);

    return 0;
}

/*
 *  This routine returns the equinox-of-date orbital elements for
 *  the earth.  Time is barycentric dynamical.
 */
static void
oearth(tjd, kep)
double tjd;
struct KEPLER *kep;
{
    double t, a0, a1;

    /* Centuries from 1900 January 0, 12h UT */
    t = (tjd - 15020e0) / 36525e0;

    /* Mean anomaly of the earth (and sun) and argument of perihelion */
    a0 = ((-3.3e-6 * t - 1.5e-4) * t + 35999.04975e0) * t + 358.47583e0;
    a1 = (3.025e-4 * t + 36000.76892e0) * t + 99.69668e0 - a0;

    kep->sma = 1.0000002e0;
    kep->ecc = (-1.26e-7 * t - 4.18e-5) * t + 1.675104e-2;
    kep->inc = 0.0;
    kep->arp = DTOR(AN360(a1));
    kep->ran = 0.0;
    kep->man = DTOR(AN360(a0));

    return;
}

/*
 *  Perturbations of the earth's orbit added in after solving Kepler's
 *  equation.  These adjust the earth's ecliptic longitude and radius.
 */
static void
cearth(tjd, el, er)
double tjd;
double *el, *er;
{
    double t, a, b, c, d, e, h, a0;

    /* Centuries from 1900 January 0, 12h UT */
    t = (tjd - 15020e0) / 36525e0;

    /* Perturbations due to Venus: */
    a = DTOR(22518.7541e0 * t + 153.23e0);
    b = DTOR(45037.5082e0 * t + 216.57e0);

    /* Due to Jupiter: */
    c = DTOR(32964.3577e0 * t + 312.69e0);
    h = DTOR(65928.7155e0 * t + 353.40e0);

    /* Due to the moon: */
    d = DTOR((-1.44e-3 * t + 445267.1142e0) * t + 350.74e0);

    /* Long period perturbation: */
    e = DTOR(20.20e0 * t + 231.19e0);

    /* Corrections to earth's longitude */
    a0 = 1.34e-3 * cos(a) +
            1.54e-3 * cos(b) +
            2.00e-3 * cos(c) +
            1.79e-3 * sin(d) +
            1.78e-3 * sin(e);

    *el = *el + DTOR(a0);

    /* Corrections to earth's radius vector (in AU) */
    a0 = 5.430e-6 * sin(a) +
            1.575e-5 * sin(b) +
            1.627e-5 * sin(c) +
            3.076e-5 * cos(d) +
            9.270e-6 * sin(h);

    *er = *er + a0;
}
