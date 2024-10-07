#include <math.h>

double fintexp3(float xtau);
double fintexp1(float xtau);

float csalbr(float xtau) {
    return (3 * xtau - fintexp3(xtau)*(4 + 2 * xtau) + 2 * expf(-xtau)) / (4 + 3 * xtau);
}

double fintexp3(float xtau) {
    return (expf(-xtau)*(1. - xtau) + xtau * xtau * fintexp1(xtau)) / 2.;
}

double fintexp1(float xtau) {
    double xx, xftau;
    int i;
    const float a[6] = {-.57721566, 0.99999193, -0.24991055,
        0.05519968, -0.00976004, 0.00107857};
    xx = a[0];
    xftau = 1.;
    for (i = 1; i < 6; i++) {
        xftau *= xtau;
        xx += a[i] * xftau;
    }
    return xx - logf(xtau);
}
