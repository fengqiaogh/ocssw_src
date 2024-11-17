//
// Created by alex on 11/25/23.
//
#include "sun2000.h"
#include "nav.h"
double dpsi, eps;
int nutime = -9999;
void sun2000(int iyr, int iday, double sec, float sun[3], float* rs) {
    const static int imon = 1;
    const double static xk = 0.0056932;  //              !Constant of aberration
    //  Compute floating point days since Jan 1.5, 2000
    //    Note that the Julian day starts at noon on the specified date
    double t = jd(iyr, imon, iday) - 2451545.0e0 + (sec - 43200.e0) / 86400.e0;
    double xls, gs, xlm, omega;
    //    extern double dpsi, eps;
    //    extern int nutime;
    //  Compute solar ephemeris parameters
    ephparms(t, &xls, &gs, &xlm, &omega);
    //  Check if need to compute nutation corrections for this day
    int nt = (int)t;
    if (nt != nutime) {
        nutime = nt;
        nutate(t, xls, gs, xlm, omega, &dpsi, &eps);
    }

    //    c  Compute planet mean anomalies
    //        c   Venus Mean Anomaly
    double g2 = 50.40828e0 + 1.60213022e0 * t;
    g2 = dmod(g2, 360.e0);

    // c Mars Mean Anomaly
    double g4 = 19.38816e0 + 0.52402078e0 * t;
    g4 = dmod(g4, 360.e0);

    //      c Jupiter Mean Anomaly
    double g5 = 20.35116e0 + 0.08309121e0 * t;
    g5 = dmod(g5, 360.e0);

    // Compute solar distance(AU)
    *rs = 1.00014e0 - 0.01671e0 * cos(gs / radeg) - 0.00014e0 * cos(2.0e0 * gs / radeg);

    //  c  Compute Geometric Solar Longitude
    double dls = (6893.0e0 - 4.6543463e-4 * t) * sin(gs / radeg) + 72.0e0 * sin(2.0e0 * gs / radeg) -
                 7.0e0 * cos((gs - g5) / radeg) + 6.0e0 * sin((xlm - xls) / radeg) +
                 5.0e0 * sin((4.0e0 * gs - 8.0e0 * g4 + 3.0e0 * g5) / radeg) -
                 5.0e0 * cos((2.0e0 * gs - 2.0e0 * g2) / radeg) - 4.0e0 * sin((gs - g2) / radeg) +
                 4.0e0 * cos((4.0e0 * gs - 8.0e0 * g4 + 3.0e0 * g5) / radeg) +
                 3.0e0 * sin((2.0e0 * gs - 2.0e0 * g2) / radeg) - 3.0e0 * sin(g5 / radeg) -
                 3.0e0 * sin((2.0e0 * gs - 2.0e0 * g5) / radeg);  //! arcseconds
    double xlsg = xls + dls / 3600.e0;
    // c  Compute Apparent Solar Longitude; includes corrections for nutation
    //   in longitude and velocity aberration
    double xlsa = xlsg + dpsi - xk / (*rs);
    //   Compute unit Sun vector
    sun[0] = cos(xlsa / radeg);
    sun[1] = sin(xlsa / radeg) * cos(eps / radeg);
    sun[2] = sin(xlsa / radeg) * sin(eps / radeg);
}