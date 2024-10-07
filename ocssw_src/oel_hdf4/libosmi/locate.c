/*
 *----------------------------------------------------------------------
 * @(#) locate.c	1.0	11 Aug 95	<shc>
 *  Copyright (c) 1995, CSIRO Division of Oceanography.
 *----------------------------------------------------------------------
 *
 * Locate the rectangular coordinates of the point of intersection of
 * a satellite scanner ray and the surface of a spheriod.  A simple
 * port to C of the Fortran routine in Puccinelli's paper.
 *
 * Inputs:
 *	svec  - Satellite state vector (X, Y, Z, Xdot, Ydot, Zdot)
 *	sat   - Satellite orientation.  Rotations are made in the
 *		order pitch, roll, yaw, and each rotation is made
 *		about the axis system formed by the previous rotation.
 *	scan  - Orientation of the scanner look vector relative to the
 *		spacecraft axes.  Equivalent to rotating the nadir
 *		vector about the pitch axis through the pitch angle
 *		and then about the roll axis by the roll angle.  Yaw
 *		is ignored.
 *
 * The reference frame axes are as follows (origin at the satellite):
 *
 *	Yaw:		normal to spheroid
 *	Velocity:	coincident with velocity vector
 *	Pitch:		yaw x velocity
 *
 * Rotations have the normal sense about the axes (e.g. positive pitch
 * causes the scanner to look ahead of the spacecraft.
 *
 * Outputs:
 *	Function returns 1 if the scanner misses the earth, 2 if it
 *	looks directly away from the earth, and 0 if the scanner sees
 *	the earth, in which case the X, Y, Z coordinates of the intercept
 *	of the scanner look vector with the surface of the earth are
 *	returned in "loc".
 *
 * Reference:
 *
 *	'Ground Location of Satellite Scanner Data',
 *	E. Puccinelli, Photogrammetric Engineering and
 *	Remote Sensing, Vol 42, April 1976, pp 537-543.
 */

#include <math.h>

#include "orbit.h"

#define TINY (1.0e-10)

int
locate(svec, sat, scan, loc)
double svec[6]; /* Spacecraft state vector */
struct ATTITUDE *sat; /* Satellite attitude */
struct ATTITUDE *scan; /* Scanner vector angles */
double loc[3]; /* Intersection coordinates */
{
    double a, adr, anbla1, anbla2, b, c, cp, cr, csp, csr, cy, d11, d12,
            d13, d21, d22, d23, d31, d32, d33, delta, ersq, f, gx, gy, gz,
            m1, m2, m3, phi, prsq, r1, r2, r3, rad, rxy, sdelta, sp, sr,
            ssp, ssr, sy, u3, vl, ynorm[3], ynorml;
    ersq = EQRAD * EQRAD;
    prsq = POLRAD * POLRAD;

    /* Normal to the spheroid (Desai's algorithm) */

    ynorm[0] = -svec[0];
    ynorm[1] = -svec[1];
    ynorm[2] = -svec[2];

    if (fabs(svec[2]) > TINY) {
        rxy = sqrt(svec[0] * svec[0] + svec[1] * svec[1]);
        if (fabs(rxy) > TINY) {
            rad = sqrt(svec[0] * svec[0] +
                    svec[1] * svec[1] + svec[2] * svec[2]);
            delta = atan2(svec[2], rxy);
            sdelta = sin(delta);
            adr = EQRAD / rad;
            anbla1 = adr * sin(2.0 * delta) / 2.0;
            anbla2 = adr + (0.5 - 2.0 * adr) * sdelta*sdelta;
            phi = delta + anbla1 * ECCSQ * (1.0 + anbla2 * ECCSQ);
            ynorm[2] = -rxy * tan(phi);
        }
    }

    /* Satellite attitude trig values */

    cy = cos(sat->yaw);
    cp = cos(sat->pitch);
    cr = cos(sat->roll);
    sy = sin(sat->yaw);
    sp = sin(sat->pitch);
    sr = sin(sat->roll);

    /* Scanner orientation trig values */

    csp = cos(scan->pitch);
    csr = cos(scan->roll);
    ssp = sin(scan->pitch);
    ssr = sin(scan->roll);

    /* Vector m */

    m1 = csr * ssp;
    m2 = -ssr;
    m3 = csp * csr;

    /* R = M . m */

    r1 = m1 * (cp * cy + sp * sr * sy) + m2 * (sp * sr * cy - cp * sy) + m3 * cr*sp;
    r2 = m1 * sy * cr + m2 * cr * cy - m3*sr;
    r3 = m1 * (cp * sr * sy - sp * cy) + m2 * (sp * sy + cp * sr * cy) + m3 * cp*cr;

    /* D matrix */

    ynorml = sqrt(ynorm[0] * ynorm[0] +
            ynorm[1] * ynorm[1] + ynorm[2] * ynorm[2]);

    d13 = ynorm[0] / ynorml;
    d23 = ynorm[1] / ynorml;
    d33 = ynorm[2] / ynorml;

    d12 = d23 * svec[5] - d33 * svec[4];
    d22 = d33 * svec[3] - d13 * svec[5];
    d32 = d13 * svec[4] - d23 * svec[3];

    vl = sqrt(d12 * d12 + d22 * d22 + d32 * d32);
    d12 = d12 / vl;
    d22 = d22 / vl;
    d32 = d32 / vl;

    d11 = d22 * d33 - d23*d32;
    d21 = d32 * d13 - d33*d12;
    d31 = d12 * d23 - d13*d22;

    /* Unit vector g = D * M * m */

    gx = r1 * d11 + r2 * d12 + r3*d13;
    gy = r1 * d21 + r2 * d22 + r3*d23;
    gz = r1 * d31 + r2 * d32 + r3*d33;

    /* Values for calculating distance from scanner to intercept point */

    a = prsq * (gx * gx + gy * gy) + ersq * gz*gz;
    b = prsq * (svec[0] * gx + svec[1] * gy) + ersq * svec[2] * gz;
    c = prsq * (svec[0] * svec[0] + svec[1] * svec[1]) +
            ersq * (svec[2] * svec[2] - prsq);
    f = b * b - a*c;

    /* If f is negative, the scanner misses the earth */

    if (f < 0)
        return (1);

    /*
     * Calculate the distance from the scanner to the intercept point.
     * If u3 is negative, the scanner is looking away from the body.
     */

    u3 = -(b + sqrt(f)) / a;
    if (u3 < 0)
        return (2);

    /* Calculate the intercept point */
    loc[0] = svec[0] + u3*gx;
    loc[1] = svec[1] + u3*gy;
    loc[2] = svec[2] + u3*gz;
    return (0);
}


#ifdef __TEST__

main() {
    int miss;
    double svec[6], loc[3];
    struct ATTITUDE sat, scanner;

    scanf("%lf %lf %lf %lf %lf %lf", &svec[0], &svec[1], &svec[2],
            &svec[3], &svec[4], &svec[5]);
    scanf("%lf %lf %lf", &sat.pitch, &sat.roll, &sat.yaw);
    scanf("%lf %lf %lf", &scanner.yaw, &scanner.pitch, &scanner.roll);

    miss = locate(svec, &sat, &scanner, loc);

    printf("%lf %lf %lf %d\n", loc[0], loc[1], loc[2], miss);
}
#endif
