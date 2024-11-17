#include "nav.h"

double norm(const double* vec, size_t n) {
    return sqrt(square(vec, n));
}

double square(const double* vec, size_t n) {
    double norm_ = 0.0;
    for (size_t i = 0; i < n; i++) {
        norm_ += vec[i] * vec[i];
    }
    return norm_;
}

double dot(const double* vec1, const double* vec2, size_t n) {
    double ans = 0;
    for (size_t i = 0; i < n; i++) {
        ans += vec1[i] * vec2[i];
    }
    return ans;
}

double dmod(double a, double p) {
    return a - ((int32_t)(a / p) * p);
}




double julian(double tin[2]) {
    const static double d0 = 2415020.5e0, d1 = 1.e-4, d2 = 1.e-2, d3 = 24.e0, d4 = 1.44e3, d5 = 8.64e4;
    double tout = d0;
    int32_t iy = (int)(tin[0] * d1 - 1.9e3);
    int32_t im = (int)(tin[0] * d2 - 1.9e5) - iy * 100;
    int32_t id = (int)(tin[0] - 1.9e7) - iy * 10000 - im * 100;
    int32_t ihour = (int)(tin[1] * d1);
    int32_t imin = (int)(tin[1] * d2) - ihour * 100;
    double sec = tin[1] - ihour * 10000 - imin * 100;
    int32_t jd = iy * 365 + (iy - 1) / 4;
    int32_t im1 = im - 1;
    switch (im1) {
        case 1:
            jd += 31;
            break;
        case 2:
            jd += 59;
            break;
        case 3:
            jd += 90;
            break;
        case 4:
            jd += 120;
            break;
        case 5:
            jd += 151;
            break;
        case 6:
            jd += 181;
            break;
        case 7:
            jd += 212;
            break;
        case 8:
            jd += 243;
            break;
        case 9:
            jd += 273;
            break;
        case 10:
            jd += 304;
            break;
        case 11:
            jd += 334;
            break;
        default:
            break;
    }
    if (iy / 4 * 4 - iy == 0 && im > 2)
        jd += 1;
    jd += id - 1;
    tout += (double)jd + (double)ihour / d3 + (double)imin / d4 + sec / d5;
    return tout;
}


void jdate(int jd, int* i, int* k) {
    // c       Compute days since January 0, 1900
    int l = jd - 2415020;

    // c       Compute years since 1900
    *i = 4 * l / 1461;

    // c       Compute day-of-year
    *k = l - 1461 * (*i - 1) / 4 - 365;

    // c       Add first two digits of year
    *i = *i + 1900;
}

void jddate(int jd, int* i, int* j, int* k) {
    int l = jd + 68569;
    int n = 4 * l / 146097;
    l = l - (146097 * n + 3) / 4;
    *i = 4000 * (l + 1) / 1461001;
    l = l - 1461 * *i / 4 + 31;
    *j = 80 * l / 2447;
    *k = l - 2447 * *j / 80;
    l = *j / 11;
    *j = *j + 2 - 12 * l;
    *i = 100 * (n - 49) + *i + l;
}

void ymdhms2jul(int32_t year, int32_t month, int32_t day, int32_t hour, int32_t minute, double sec,
                double* jul) {
    *jul = jd(year, month, day) + (hour * 3600.e0 + minute * 60.e0 + sec) / 86400.e0;
}


int32_t jd(int32_t i, int32_t j, int32_t k) {
    return 367 * i - 7 * (i + (j + 9) / 12) / 4 + 275 * j / 9 + k + 1721014;
    // c  This additional calculation is needed only for dates outside of the
    // c   period March 1, 1900 to February 28, 2100
    // c       jd = jd + 15 - 3*((i+(j-9)/7)/100+1)/4
};