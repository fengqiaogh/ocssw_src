#include "libnav.h"
#include "sun2000.h"
#include "nav.h"
/**
 * @brief c  Computes unit Sun vector in geocentric rotating coodinates, using
          c  subprograms to compute inertial Sun vector and Greenwich hour angle
c       Subprograms referenced:
c
c       SUN2000         Computes inertial Sun vector
c       GHA2000         Computes Greenwich sidereal angle
c
c       Coded by:  Frederick S. Patt, GSC, September 29, 1992
c
c       Modification History:
c
c       Modifified to use new Sun and hour angle routines
c       Frederick S. Patt, November 3, 1992
c
c       Removed internal jd() function, since it is available as an
c       independent module.  B. A. Franz, GSC, November 14, 1997.
 * @param iyr input      Year, four digits (i.e, 1993)
 * @param idoy input  Day of year (1-366)
 * @param sec input Seconds of day
 * @param sunr output Unit Sun vector in geocentric rotating coordinates
 * @param rs Earth-to-Sun distance (AU)
 */

void l_sun_(int *iyr, int *iday, double *sec, float sunr[3], float *rs) {
    float su[3];
    // Get unit Sun vector in geocentric inertial coordinates
    sun2000(*iyr, *iday, *sec, su, rs);
    //  Get Greenwich mean sideral angle
    double day = *iday + *sec / 864.e2;
    double gha;
    gha2000(*iyr, day, &gha);
    double ghar = gha / radeg;
    // Transform Sun vector into geocentric rotating frame
    sunr[0] = su[0]*cos(ghar) + su[1]*sin(ghar);
    sunr[1] = su[1]*cos(ghar) - su[0]*sin(ghar);
    sunr[2] = su[2];
}
