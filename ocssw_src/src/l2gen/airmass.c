#include "l12_proto.h"

/* airmass for a spherical atmosphere, Kasten & Young, 1989 */
float ky_airmass(float theta) {
    float mu = MAX(cos(theta / RADEG), 0.01);
    return (1. / (mu + pow(0.50572 * (96.07995 - theta), (-1.6364))));
}

float ky_airmass_(float *theta) {
    return (ky_airmass(*theta));
}

/* airmass for a plane parallel atmosphere */
float pp_airmass(float theta) {
    float mu = MAX(cos(theta / RADEG), 0.01);
    return (1. / mu);
}

float pp_airmass_(float *theta) {
    return (pp_airmass(*theta));
}
