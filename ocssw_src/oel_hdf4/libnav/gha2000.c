#include "sun2000.h"
#include "nav.h"
extern double dpsi, eps;
extern int nutime;
void gha2000(int iyr, double day, double* gha) {
    const static int imon = 1;
    //  Compute days since J2000
    int iday = (int)day;
    double fday = day - iday;
    int jday = jd(iyr, imon, iday);
    double t = jday - 2451545.5e0 + fday;
    //  Compute Greenwich Mean Sidereal Time (degrees)
    double gmst = 100.4606184e0 + 0.9856473663e0 * t + 2.908e-13 * t * t;
    // Check if need to compute nutation correction for this day
    int nt = (int)t;
    if (nt != nutime) {
        nutime = nt;
        double xls, gs, xlm, omega;
        ephparms(t, &xls, &gs, &xlm, &omega);
        nutate(t, xls, gs, xlm, omega, &dpsi, &eps);
    }
    //  Include apparent time correction and time-of-day
    *gha = gmst + dpsi * cos(eps / radeg) + fday * 360.e0;
    *gha = dmod(*gha, 360.e0);
    if (*gha < 0.e0)
        *gha = *gha + 360.e0;
}