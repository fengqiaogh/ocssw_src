/* -----------------------------------------------------------
   simple linear interpolation (BAF, 6/99)

 Inputs:
   float  xin (nin )    - input independent variable
   float  yin (nin )    - input dependent variable
   int32_t   nin           - number of input points
   float  xout(nout)    - output independent variable
   int32_t   nout          - number of output points

 Outputs:
   float  yin (nout)    - output dependent variable

 Warning: xin, xout assumed to be ascending order

 ------------------------------------------------------------- */

#include <stdint.h>

void lspline(float xin [], float yin [], int32_t nin,
        float xout[], float yout[], int32_t nout) {
    float a, b, x;
    int32_t i, j;

    /* Can't interpolate if only 1 input point */

    if (nin < 2) {
        for (i = 0; i < nout; i++)
            yout[i] = yin[0];
        return;
    }

    j = 0;

    for (i = 0; i < nout; i++) {

        while (xin[j] < xout[i]) {
            j = j + 1;
            if (j >= nin) {
                j = nin - 1;
                break;
            }
        }

        if (xout[i] == xin[j]) {
            yout[i] = yin[j];
        } else {
            if (j < 1) j = 1;
            a = (yin[j] - yin[j - 1]) / (xin[j] - xin[j - 1]);
            b = yin[j - 1];
            x = xout[i] - xin[j - 1];
            yout[i] = a * x + b;
        }
    }

    return;
}

float linterp(float xin [], float yin [], int32_t nin, float xout) {
    float yout;
    lspline(xin, yin, nin, &xout, &yout, 1);
    return (yout);
}

#ifdef LSPINE_TEST

int main(int argc, char *argv[]) {
    int i;
    int32_t order = 2;
    int32_t n = 5;

    float x1[5] = {1, 2, 3, 4, 5};
    float y1[5] = {2.5, 4.5, 6.5, 8.5, 10.5};
    float x2[9] = {-0.5, 0.5, 1.5, 1.9, 2.1, 3.5, 4.5, 5.0, 6.5};
    float y2[9];

    float ylgint_(float x1[], float y1[], int32_t *n, float *x2, int32_t * order);

    lspline(x1, y1, 5, x2, y2, 8);

    for (i = 0; i < 7; i++)
        printf("%d %f %f %f %f\n", i, x2[i], y2[i], linterp(x1, y1, n, x2[i]), ylgint_(x1, y1, &n, &x2[i], &order));

    exit(0);
}
#endif
