#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "amoeba.h"

#define GET_PSUM \
for (j=0;j<ndim;j++) \
{ \
  for (sum=0.0,i=0;i<ndim+1;i++) \
  sum += *(pnts + i*ndim +j); \
  psum[j] = sum;\
}

short amoeba(double *pnts, FITSTRUCT * auxdata,
        double (*func)(FITSTRUCT *, double []), float tol) {
    short i, j, k, isml = 0, inxt, ibig, ndim;
    double y[16], ptry[16], psum[16], ytry, ysave, sum, rtol;
    double reflect(double *, double [], double [], short, FITSTRUCT *,
            double (*func)(FITSTRUCT *, double []), short, float);
    int32_t itermax = auxdata->niter;


    ndim = auxdata->nfunc;


    GET_PSUM

            /* Determine initial functional values */
    for (i = 0; i < ndim + 1; i++) {
        for (j = 0; j < ndim; j++) ptry[j] = pnts[i * ndim + j];
        y[i] = (*func)(auxdata, ptry);
    }

    for (k = 0; k < itermax; k++) {

        auxdata->niter = k + 1;

        /* Order by smallest, second-largest, largest */
        if (y[0] > y[1]) ibig = 0;
        else ibig = 1;
        isml = 1 - ibig;
        inxt = isml;

        for (i = 2; i < ndim + 1; i++) {
            if (y[i] > y[ibig]) {
                inxt = ibig;
                ibig = i;
            } else if (y[i] < y[isml]) {
                isml = i;
            } else if (y[i] > y[inxt] && y[i] < y[ibig]) inxt = i;
        }

        rtol = 2.0 * fabs(y[ibig] - y[isml]) /
                (fabs(y[ibig]) + fabs(y[isml]));

        /*    printf("%f %f\n", rtol,tol);*/

        if (rtol < tol) break;

        ytry = reflect(pnts, y, psum, ndim, auxdata, func, ibig, -1.0);

        if (ytry <= y[isml]) {
            ytry = reflect(pnts, y, psum, ndim, auxdata, func, ibig, 2.0);
        }
        else if (ytry >= y[inxt]) {
            ysave = y[ibig];
            ytry = reflect(pnts, y, psum, ndim, auxdata, func, ibig, 0.5);
            if (ytry >= ysave) {
                for (i = 0; i < ndim + 1; i++) {
                    if (i != isml) {
                        for (j = 0; j < ndim; j++)
                            pnts[i * ndim + j] = psum[j] =
                                0.5 * (pnts[i * ndim + j] + pnts[isml * ndim + j]);
                    }
                }
                GET_PSUM
            }
        }

    }

    auxdata->niter = k + 1;

    return (isml);
}

double reflect(double *pnts, double y[], double psum[], short ndim,
        FITSTRUCT *auxdata, double (*func)(FITSTRUCT *, double []),
        short ibig, float fac) {

    short j;
    float fac1, fac2;
    double ytry, ptry[16];

    fac1 = (1.0 - fac) / ndim;
    fac2 = fac1 - fac;


    for (j = 0; j < ndim; j++) ptry[j] = psum[j] * fac1 - pnts[ibig * ndim + j] * fac2;

    ytry = (*func)(auxdata, ptry);

    if (ytry < y[ibig]) {
        y[ibig] = ytry;
        for (j = 0; j < ndim; j++) {
            psum[j] += ptry[j] - pnts[ibig * ndim + j];
            pnts[ibig * ndim + j] = ptry[j];
        }
    }

    return (ytry);

}









