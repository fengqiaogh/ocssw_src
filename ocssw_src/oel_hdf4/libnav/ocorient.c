#include "libnav.h"
#include "nav.h"
#include "ocorient.h"
#include "math_utils.h"
void matmpy(const float xm1[3][3], const float xm2[3][3], float xm3[3][3]) {
    int i, j, m, n;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            xm3[i][j] = 0;
            for (m = 0; m < 3; m++)
                xm3[i][j] += xm1[i][m] * xm2[m][j];
        }
    }
}

void matmpy_(const float xm1[3][3], const float xm2[3][3], float xm3[3][3]) {
    float xsave3[3][3];
    float xsave1[3][3];
    float xsave2[3][3];
    transpose3d(xm1, xsave1);
    transpose3d(xm2, xsave2);
    matmpy(xsave1, xsave2, xsave3);
    transpose3d(xsave3, xm3);
}
void oceuler(float a[3], float xm[3][3]) {
    float xm1[3][3], xm2[3][3], xm3[3][3], xmm[3][3];

    //  Initialize all matrix elements to zero.
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            xm1[i][j] = xm2[i][j] = xm3[i][j] = 0.0f;
        }
    }
    //    c  Compute sines and cosines
    float c1 = cos(a[0] / radeg);
    float s1 = sin(a[0] / radeg);
    float c2 = cos(a[1] / radeg);
    float s2 = sin(a[1] / radeg);
    float c3 = cos(a[2] / radeg);
    float s3 = sin(a[2] / radeg);
    // c Convert individual rotations to matrices xm1(1, 1) = 1.d0;
    xm1[1 - 1][1 - 1] = 1.0f;
    xm1[2 - 1][2 - 1] = c1;
    xm1[3 - 1][3 - 1] = c1;
    xm1[2 - 1][3 - 1] = s1;
    xm1[3 - 1][2 - 1] = -s1;
    xm2[2 - 1][2 - 1] = 1.0f;
    xm2[1 - 1][1 - 1] = c2;
    xm2[3 - 1][3 - 1] = c2;
    xm2[3 - 1][1 - 1] = s2;
    xm2[1 - 1][3 - 1] = -s2;
    xm3[3 - 1][3 - 1] = 1.0f;
    xm3[2 - 1][2 - 1] = c3;
    xm3[1 - 1][1 - 1] = c3;
    xm3[1 - 1][2 - 1] = s3;
    xm3[2 - 1][1 - 1] = -s3;
    matmpy(xm2, xm3, xmm);
    matmpy(xm1, xmm, xm);
}

/**
 * @brief
 *      c  This subroutine performs a simple calculation of the sensor
        c  orientation from the orbit position vector and input values of the
        c  attitude offset angles.  The calculations assume that the angles
        c  represent the roll, pitch and yaw offsets between the local vertical
        c  reference frame (at the spacecraft position) and the sensor frame.
        c  Sensor tilt angles are assumed to be included in the pitch angle.
        c  The outputs are the matrix which represents the transformation from
        c  the geocentric rotating to sensor frame, and the coefficients which
        c  represent the Earth scan track in the sensor frame.
        c  The reference ellipsoid uses an equatorial radius of 6378.137 km and
        c  a flattening factor of 1/298.257 (WGS 1984).
        c       Subprograms Called (attached):
        c       CROSSP          Compute cross product of two vectors
        c       EULER           Compute matrix from Euler angles
        c       MATMPY          Multiply two 3x3 matrices
        c
        c       Program written by:     Frederick S. Patt
        c                               General Sciences Corporation
        c                               July 22, 1992
        c
        c       Modification History:
        c       Added improved calculation of local vertical reference frame and
        c        modified calling argument names.
        c        F. S. Patt, September 30, 1992
        c       Expanded vector normalization in-line to eliminate subroutine
        c        call.  F. S. Patt, October 19, 1992
        c       Eliminated redundant calculations, changed to correspond to
        c       paper, "Exact closed-form geolocation algorithm for earth
        c       survey sensors", International Journal of Remote Sensing, Patt
        c       and Gregg, 1993.  W. Gregg, 4/5/93.
        c       Modified to support three-dimensional view vectors by computing
        c       coefficients array with all 10 ellipsoid terms.
        c       F. S. Patt, November 25, 1996
 * @param pos I      Orbit Position Vector ECEF (km)
 * @param vel        Orbit Velocity Vector ECEF (km/sec)
 * @param att        Attitude Offsets (Euler angles in  degrees); referenced to local vertical   coordinates;
 order is roll, pitch, yaw (X, Y, Z)
 * @param rm        Sensor orientation matrix
 * @param coef       Scan path coefficients
 */
void ocorient_(float pos[3], float vel[3], float att[3], float rm[3][3], float coef[10]) {
    float vc[3], xpri[3], ypri[3], zpri[3], yrp[3][3], rn[3][3];
    /**
     * c  Compute correction to orbit velocity vector in Earth-centered
      c   Earth-fixed (ECEF) frame; this involves subtracting effect of Earth
      c   rotation rate on velocity to get correct scan plane orientation in
      c   ECEF frame.
     */

    vc[1 - 1] = vel[1 - 1] - omegae * pos[2 - 1];
    vc[2 - 1] = vel[2 - 1] + omegae * pos[1 - 1];
    vc[3 - 1] = vel[3 - 1];
    /**
     * c  Determine nadir frame reference axes
       c  Uses method of local ellipsoid approximation good to 0.3 arcsecond
       c   Compute Z axis as local nadir vector
     */

    double pm = sqrt((double)(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]));
    double omf2p = (omf2 * rem + pm - rem) / pm;
    double pxy = pos[0] * pos[0] + pos[1] * pos[1];
    double temp = sqrt((double)(pos[2] * pos[2] + omf2p * omf2p * pxy));
    zpri[0] = -omf2p * pos[0] / temp;
    zpri[1] = -omf2p * pos[1] / temp;
    zpri[2] = -pos[2] / temp;
    // c   Compute Y axis along negative orbit normal
    crossp_(vc, zpri, ypri);
    double yprim = sqrt((double)(ypri[0] * ypri[0] + ypri[1] * ypri[1] + ypri[2] * ypri[2]));
    // normalization to length 1
    for (int i = 0; i < 3; i++)
        ypri[i] /= -yprim;
    // c   Compute X axis to complete orthonormal triad
    crossp_(ypri, zpri, xpri);
    // c  Store in matrix
    for (int i = 0; i < 3; i++) {
        rn[0][i] = xpri[i];
        rn[1][i] = ypri[i];
        rn[2][i] = zpri[i];
    }
    // c  Convert attitude (Euler) angles to YRP matrix
    oceuler(att, yrp);
    // c  Compute sensor orientation matrix
    float rm_temp[3][3];
    matmpy(yrp, rn, rm);
    memcpy(rm_temp,rm, sizeof rm_temp);
//    transpose3d(rm_temp,rm);
    // c  Compute coefficients of intersection ellipse in scan plane
    double rd = 1.e0 / omf2;
    *(coef + 1 - 1) = 1.e0 + (rd - 1.e0) * rm[1 - 1][3 - 1] * rm[1 - 1][3 - 1];
    *(coef + 2 - 1) = 1.e0 + (rd - 1.e0) * rm[2 - 1][3 - 1] * rm[2 - 1][3 - 1];
    *(coef + 3 - 1) = 1.e0 + (rd - 1.e0) * rm[3 - 1][3 - 1] * rm[3 - 1][3 - 1];
    *(coef + 4 - 1) = (rd - 1.e0) * rm[1 - 1][3 - 1] * rm[2 - 1][3 - 1] * 2.e0;
    *(coef + 5 - 1) = (rd - 1.e0) * rm[1 - 1][3 - 1] * rm[3 - 1][3 - 1] * 2.e0;
    *(coef + 6 - 1) = (rd - 1.e0) * rm[2 - 1][3 - 1] * rm[3 - 1][3 - 1] * 2.e0;
    *(coef + 7 - 1) =
        (rm[1 - 1][1 - 1] * pos[1 - 1] + rm[1 - 1][2 - 1] * pos[2 - 1] + rm[1 - 1][3 - 1] * pos[3 - 1] * rd) *
        2.e0;
    *(coef + 8 - 1) =
        (rm[2 - 1][1 - 1] * pos[1 - 1] + rm[2 - 1][2 - 1] * pos[2 - 1] + rm[2 - 1][3 - 1] * pos[3 - 1] * rd) *
        2.e0;
    *(coef + 9 - 1) =
        (rm[3 - 1][1 - 1] * pos[1 - 1] + rm[3 - 1][2 - 1] * pos[2 - 1] + rm[3 - 1][3 - 1] * pos[3 - 1] * rd) *
        2.e0;
    *(coef + 10 - 1) = pos[1 - 1] * pos[1 - 1] + pos[2 - 1] * pos[2 - 1] + pos[3 - 1] * pos[3 - 1] * rd -
                       re_const * re_const;
    transpose3d(rm_temp,rm);
}