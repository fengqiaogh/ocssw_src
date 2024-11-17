//
// Created by alex on 11/25/23.
//
#include "sun2000.h"
#include "nav.h"

void nutate(double t, double xls, double gs, double xlm, double omega, double *dpsi, double *eps) {
    //    c  Nutation in Longitude
    *dpsi = -17.1996e0 * sin(omega / radeg) + 0.2062e0 * sin(2.0e0 * omega / radeg) -
            1.3187e0 * sin(2.0e0 * xls / radeg) + 0.1426e0 * sin(gs / radeg) -
            0.2274e0 * sin(2.0e0 * xlm / radeg);
    //    c  Mean Obliquity of the Ecliptic
    double epsm = 23.439291e0 - 3.560e-7 * t;
    //  Nutation in Obliquity
    double deps = 9.2025e0 * cos(omega / radeg) + 0.5736e0 * cos(2.0e0 * xls / radeg);
    // True Obliquity of the Ecliptic
    *eps = epsm + deps / 3600.e0;
    *dpsi = *dpsi / 3600.e0;
}