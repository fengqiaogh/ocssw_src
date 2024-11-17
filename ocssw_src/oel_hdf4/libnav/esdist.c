#include "libnav.h"
#include "nav.h"
/**
 * @brief c  This subroutine computes the earth-sun distance in AU. It uses the model
            c  referenced in The Astronomical Almanac for 1984, Section S (Supplement).

            c       Subprograms referenced:
            c
            c       JD              Computes Julian day from calendar date
            c
            c       Coded by:  Frederick S. Patt, GSC, November 2, 1992
            c       Adapted from sun2000 to esdist by B. Franz, Oct 2003.

 * @param year
 * @param day
 * @param msec
 * @return
 */
 // something wrong with fsol? Need to check
 // check jd
double esdist_(int *iyr, int *iday, int *msec) {
    /**
     *  c   Compute floating point days since Jan 1.5, 2000
        c   Note that the Julian day starts at noon on the specified date
     */
    double t = jd(*iyr, 1, *iday) - 2451545.0e0 + (*msec / 1000.e0 - 43200.e0) / 86400.e0;
    // c  Compute mean anomaly
    double gs = 357.52772e0 + 0.9856002831e0 * t;
    //   Compute solar distance (AU)
    const static double radeg_to_match_fortran = 57.295780181884766e0;
    return 1.00014e0 - 0.01671e0 * cos(gs / radeg_to_match_fortran) - 0.00014e0 * cos(2.0e0 * gs / radeg_to_match_fortran);
}
