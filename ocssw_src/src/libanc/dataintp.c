#include <math.h>
#include "anc.h"

#define DSIGN(A, B) (B >= 0 ? fabs(A) : -fabs(A))
#define ERROR_RET(x, err) {err = x; if(err) return err;}
#define MAX_BAND 15
int timeint(float f1[][MAX_BAND], float f2[][MAX_BAND], double dt1, double dt2, int32_t ipt, int32_t *nband, float *mimx, float *def,
            float *fout, float *func, int32_t ir) {
    double t1 = fabs(dt1);
    double t2 = fabs(dt2);
    if (t1 == 0.0e0 || t2 == 0.0e0) t2 = 1.e0;
    double w1 = t2 / (t1 + t2);
    double w2 = t1 / (t1 + t2);
    if (w1 > 1.0e0) w2 = 0.0e0;
    if (w2 > 1.0e0) w1 = 0.0e0;
    for (int n = 0; n < *nband; n++) {
        int ij = n * ir;
        for (int i = 0; i < ipt; i++) {
            int k = ij + i;
            int flag1 = f1[n][i] >= mimx[0] && f1[n][i] <= mimx[1];
            int flag2 = f1[n][i] >= mimx[0] && f1[n][i] <= mimx[1];
            int sw = flag1 * 2 + flag2;
            switch (sw) {
                case 0:
                    fout[k] = def[n];//
                    func[k] = 0;
                    break;
                case 1:
                    fout[k] = f2[n][i];//
                    func[k] = 0;
                    break;
                case 2:
                    fout[k] = f1[n][i];//
                    func[k] = 0;
                    break;
                case 3:
                    fout[k] = w1 * f1[n][i] + w2 * f2[n][i];
                    func[k] = fabsf(f1[n][i] - f2[n][i]);
                    break;
                default:
                    return 1;
            }

        }
    }
    return 0;
}

int spaceint(float *ll, float *lat, float *lon, float **f, int32_t *ipt, int32_t *nband, float *mimx, float *def,
             float *fout,
             int32_t *int_bad) {
    const static int nc[4] = {1, 10, 100, 1000};
    const static int n1[4][3] = {{1,    2, 4},
                                 {10,   1, 3},
                                 {100,  2, 4},
                                 {1000, 1, 3}};
    const static int n2[6][5] = {{11,   1, 2, 4, 3},
                                 {101,  1, 3, 4, 2},
                                 {1001, 1, 4, 2, 3},
                                 {110,  2, 3, 1, 4},
                                 {1010, 2, 4, 3, 1},
                                 {1100,
                                        3, 4, 2, 1}};
    const static int n3[4][2] = {{111,  4},
                                 {1110, 1},
                                 {1101, 2},
                                 {1011, 3}};
    float llft[2], urht[2];
    const int npt = *ipt;
    const static double error_tolerance = 1e-9;
    llft[1 - 1] = lat[1 - 1];
    llft[2 - 1] = lon[1 - 1];
    urht[1 - 1] = lat[3 - 1];
    urht[2 - 1] = lon[3 - 1];
    // C: Rectangular bi-linear interpolation.
    float xp = ll[1] - llft[1];
    if (fabsf(xp) > 180.0f) { xp = DSIGN(360.0f - fabsf(xp), xp); }
    float yp = ll[0] - llft[0];
    float dx = urht[1] - llft[1];
    if (fabsf(dx) > 180.f) dx = DSIGN(360.0f - fabsf(dx), dx);
    float dy = urht[0] - llft[0];
    float dd = dx * dy;
    for (int n = 0; n < *nband; n++) {
        int ncc = 0;
        int nng = 0;
        *int_bad = 0;
        float g[npt];
        for (int i = 0; i < npt; i++) {
            g[i] = f[n][i];
            if (g[i] < mimx[0] || g[i] > mimx[1]) {
                ncc += nc[i];
                nng++;
            }
        }
        switch (nng) {
            case 0:
                break;
            case 1:
                for (int i = 0; i < npt; i++) {
                    if (ncc == n1[i][0]) {
                        g[i] = (g[n1[i][1]] + g[n1[i][2]]) / 2.f;
                    }
                }
                break;
            case 2:
                for (int i = 0; i < 6; i++) { // ! WDR replace IPT with 6 to correct
                    if (ncc == n2[i][0]) {
                        g[n2[i][1]] = g[n2[i][3]];
                        g[n2[i][2]] = g[n2[i][4]];
                    }
                }
                break;
            case 3:
                for (int i = 0; i < npt; i++) {
                    if (ncc == n3[i][0]) {
                        fout[n] = g[n3[i][1]];
                        continue;
                    }
                }
                break;
            case 4:
                fout[n] = def[n];
                *int_bad = 1;
                continue;
            default:
                return 1;
        }
        float a[4];
        a[0] = g[0];
        if (fabsf(dx) < error_tolerance)
            a[1] = 0.0f;
        else
            a[1] = (g[3] - a[0]) / dx;
        if (fabsf(dy) < error_tolerance)
            a[2] = 0.0f;
        else
            a[2] = (g[1] - a[0]) / dy;
        if (fabsf(dd) < error_tolerance)
            a[3] = 0.0f;
        else
            a[3] = (a[0] - g[1] + g[2] - g[3]) / dd;
        fout[n] = a[0] + a[1] * xp + yp * a[2] + yp * xp * a[3];
    }
    return 0;
}


int dataintp(float in_latlon[2], float *lat, float *lon,
             float *data_list1, double *dt1, float *data_list2, double *dt2,
             int32_t *ipt, int32_t *nband, float rng[2], float *def,
             int32_t *intporder, float *dummy, float *dataout, float *unc,
             int32_t *int_bad, int32_t *row, int32_t *col) {
    float mimx[2];
    float fout1[1][MAX_BAND];
    float fout2[1][MAX_BAND];
    const int ntot = *nband * *row;
    int int_bad1, int_bad2;
    if (rng[1] > rng[0]) {
        mimx[1] = rng[1];
        mimx[0] = rng[0];
    } else {
        mimx[1] = rng[0];
        mimx[0] = rng[1];
    }
    *int_bad = 0;
    for (int i = 0; i < ntot; i++) {
        unc[i] = 0.0f;
    }
    int ret_ = 0;
    switch (*intporder) {
        case 112: ERROR_RET(
                spaceint(in_latlon, lat, lon, &data_list1, ipt, nband, mimx, def, (float *) fout1, &int_bad1), ret_);
            ERROR_RET(spaceint(in_latlon, lat, lon, &data_list2, ipt, nband, mimx, def, (float *) fout2, &int_bad2),
                      ret_);
            ERROR_RET(
                    timeint( fout1, fout2, *dt1, *dt2, 1, nband, mimx, def, dataout, unc, *row),
                    ret_);
            *int_bad = int_bad1 || int_bad2;
            break;
        case 110: ERROR_RET(spaceint(in_latlon, lat, lon, &data_list1, ipt, nband, mimx, def, dataout, int_bad),
                            ret_);
            break;
        case 101: ERROR_RET(
                timeint( fout1, fout2, *dt1, *dt2, *ipt, nband, mimx, def, dataout, unc, *row),
                ret_);
            break;
        case 121: ERROR_RET(
                timeint(fout1,  fout2, *dt1, *dt2, *ipt, nband, mimx, def, dummy, unc, *ipt),
                ret_);
            ERROR_RET(spaceint(in_latlon, lat, lon, &dummy, ipt, nband, mimx, def, dataout, int_bad),
                      ret_);
        default:
            return 1;
    }
}
