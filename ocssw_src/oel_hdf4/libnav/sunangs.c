//
// Created by alex on 11/25/23.
//
#include "libnav.h"
#include "nav.h"
#include "sun2000.h"
/**
 * $brief c  Given year, day of year, time in hours (GMT) and latitude and
c  longitude, returns an accurate solar zenith and azimuth angle.
c  Based on IAU 1976 Earth ellipsoid.  Method for computing solar
c  vector and local vertical from Patt and Gregg, 1993, Int. J.
c  Remote Sensing.c
c  Subroutines required: sun2000
c                        gha2000
c                        jd
c
 * @param year
 * @param day
 * @param gmt
 * @param lon
 * @param lat
 * @param sunz
 * @param suna
 */
void sunangs_(int *iyr, int *iday, float *gmt, float *xlon, float *ylat, float *sunz, float *suna) {
    //  Compute sun vector
    //   Compute unit sun vector in geocentric inertial coordinates
    double sec = *gmt * 3600.0e0;
    float rs, suni[3], sung[3], up[3], ea[3], no[3];
    sun2000(*iyr, *iday, sec, suni, &rs);
    double day = *iday + sec / 3600.00 / 24.0;
    double gha;
    //   Get Greenwich mean sidereal angle
    gha2000(*iyr, day, &gha);
    const double ghar = gha / radeg;
    //   Transform Sun vector into geocentric rotating frame
    sung[0] = suni[0] * cos(ghar) + suni[1] * sin(ghar);
    sung[1] = suni[1] * cos(ghar) - suni[0] * sin(ghar);
    sung[2] = suni[2];
    // c Convert geodetic lat/lon to Earth-centered, earth-fixed (ECEF)
    // c vector (geodetic unit vector)
    const double rlon = *xlon / radeg;
    const double rlat = *ylat / radeg;
    const double cosy = cos(rlat);
    const double siny = sin(rlat);
    const double cosx = cos(rlon);
    const double sinx = sin(rlon);
    up[0] = cosy * cosx;
    up[1] = cosy * sinx;
    up[2] = siny;
    // c
    // c  Compute the local East and North unit vectors
    const double upxy = sqrt(up[0] * up[0] + up[1] * up[1]);
    ea[0] = -up[1] / upxy;
    ea[1] = up[0] / upxy; 
    ea[2] = 0.0e0;
    no[0] = up[1] * ea[2] - up[2] * ea[1];  // !cross product
    no[1] = up[2] * ea[0] - up[0] * ea[2];
    no[2] = up[0] * ea[1] - up[1] * ea[0];
    //    c  Compute components of spacecraft and sun vector in the
    // c  vertical (up), North (no), and East (ea) vectors frame
    double sunv = 0.0;
    double sunn = 0.0;
    double sune = 0.0;
    for (int j = 0; j < 3; j++) {
        sunv = sunv + sung[j] * up[j];
        sunn = sunn + sung[j] * no[j];
        sune = sune + sung[j] * ea[j];
    }
    // c
    // c  Compute the solar zenith and azimuth
    *sunz = radeg * atan2(sqrt(sunn * sunn + sune * sune), sunv);
    // c  Check for zenith close to zero
    if (*sunz > 0.05e0)
        *suna = radeg * atan2(sune, sunn);
    else
        *suna = 0.0e0;
    if (*suna < 0.0e0)
        *suna = *suna + 360.0e0;
}