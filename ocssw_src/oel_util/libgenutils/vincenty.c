#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#ifndef VINCENTY_TOLERANCE
/** @brief Used to control the accuracy requirement. 1e−12 corresponds to approximately 0.06mm. */
#define VINCENTY_TOLERANCE 1e-12
#endif

#ifndef VINCENTY_MAX_LOOP
/** @brief Used to prevent infinite loops caused by non-convergence. */
#define VINCENTY_MAX_LOOP 1e3
#endif

#ifndef USE_HELMERTS
/** @brief Use Helmert's expansion something-or-other.

 https://en.wikipedia.org/wiki/Vincenty's_formulae#Vincenty.27s_modification
 */
#define USE_HELMERTS 1
#endif

#define pi 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899L
#define flattening (1 / 298.257223563) // flattening of the ellipsoid
#define semimajor 6378137.0 // abbreviated a, radius at the equator
#define semiminor ((1 - flattening) * semimajor) // abbreviated b, radius at the poles

inline static double deg2rad(const double deg) {
    return (pi * deg / 180.0);
}
//inline static double rad2deg(const double rad) {
//    return (rad * 180.0 / pi);
//}

inline static void vincenty_AB_equations(const double u2, double *A, double *B) {
#ifdef USE_HELMERTS
    const double k1 = (sqrt(1 + u2) - 1) / (sqrt(1 + u2) + 1);
    *A = (1 + 1 / 4 * pow(k1, 2));
    *B = k1 * (1 - 3 / 8 * pow(k1, 2));
# else
    *A = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)));
    *B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)));
#endif
}

// https://en.wikipedia.org/wiki/Vincenty's_formulae#Inverse_problem
double vincenty_distance(double lat1, double lon1, double lat2, double lon2) {
    if (lat1 == lat2 && lon1 == lon2) { // returns -nan in this case
        return 0;
    }
    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    lat2 = deg2rad(lat2);
    lon2 = deg2rad(lon2);

    // reduced latitudes (latitudes on the auxiliary sphere)
    const double U1 = atan((1 - flattening) * tan(lat1));
    const double U2 = atan((1 - flattening) * tan(lat2));
    const double sin_U1 = sin(U1), cos_U1 = cos(U1);
    const double sin_U2 = sin(U2), cos_U2 = cos(U2);
    const double L = lon2 - lon1;

    double lambda = L;

    double old_lambda, sin_sigma, cos_sigma, sigma, sin_alpha, cos2_alpha, cos_2sigmam = 0, C;
    int loop_count = 0;
    do {
        old_lambda = lambda;
        const double cos_lambda = cos(lambda), sin_lambda = sin(lambda);
        sin_sigma = sqrt(pow(cos_U2 * sin_lambda, 2) + pow(cos_U1 * sin_U2 - sin_U1 * cos_U2 * cos_lambda, 2));
        cos_sigma = sin_U1 * sin_U2 + cos_U1 * cos_U2 * cos_lambda;
        sigma = atan2(sin_sigma, cos_sigma);
        sin_alpha = (cos_U1 * cos_U2 * sin_lambda) / sin_sigma;
        cos2_alpha = 1 - pow(sin_alpha, 2);
        if (cos2_alpha == 0) {
            lambda = L + flattening * sin_alpha * sigma;
        } else {
            cos_2sigmam = cos_sigma - (2 * sin_U1 * sin_U2) / cos2_alpha;
            C = (flattening / 16) * cos2_alpha * (4 + flattening * (4 - 3 * cos2_alpha));
            lambda = L + (1 - C) * flattening * sin_alpha * (sigma + C * sin_sigma * (cos_2sigmam + C * cos_sigma * (-1 + 2 * cos_2sigmam)));
        }
    } while (fabs(lambda - old_lambda) > VINCENTY_TOLERANCE && ++loop_count < VINCENTY_MAX_LOOP);

    const double u2 = cos2_alpha * (pow(semimajor, 2) - pow(semiminor, 2)) / pow(semiminor, 2);
    double A, B;
    vincenty_AB_equations(u2, &A, &B);
    const double cos2_2sigmam = pow(cos_2sigmam, 2);
    const double delta_sigma = B * sin_sigma * (cos_2sigmam + 1 / 4 * B * (-1 + 2 * cos2_2sigmam) - 1 / 6 * B * cos_2sigmam * (-3 + 4 * pow(sin_sigma, 2)) * (-3 + 4 * cos2_2sigmam));
    const double s = semiminor * A * (sigma - delta_sigma);

    // I think this stuff calculates the bearings, not required for this function but useful for a geolib.
//	const double cos_lambda = cos(lambda);
//	const double sin_lambda = sin(lambda);
//	const double cosU1_sin_U2 = cos_U1 * sin_U2;
//	const double alpha1 = atan((cos_U2 * sin_lambda) / (cosU1_sin_U2 - (sin_U1 * cos_U2 * cos_lambda)));
//	const double alpha2 = atan((cos_U1 * sin_lambda) / ((-sin_U1 * cos_U2) - (cosU1_sin_U2 * cos_lambda)));

//	printf("%f, %f, %f, %f = %.64f\n", rad2deg(lat1), rad2deg(lon1), rad2deg(lat2), rad2deg(lon2), s);
    return s;
}
