#include "sun2000.h"
#include "nav.h"
void ephparms(double t, double* xls, double* gs, double* xlm, double* omega) {
    // Sun Mean Longitude
    *xls = 280.46592e0 + 0.9856473516e0 * t;
    *xls = dmod(*xls, 360e0); //fmod(xls, (double)360);

    // Sun Mean Anomaly
    *gs = 357.52772e0 + 0.9856002831e0 * t;
    *gs = dmod(*gs, 360.e0); //  fmod(xls, (double)360);

    //  Moon Mean Longitude
    *xlm = 218.31643e0 + 13.17639648e0 * t;
    *xlm = dmod(*xlm, 360.e0); //  fmod(xls, (double)360);

    // Ascending Node of Moon's Mean Orbit
    *omega = 125.04452e0 - 0.0529537648e0 * t;
    *omega = dmod(*omega, 360.e0); //  fmod(xls, (double)360);
}