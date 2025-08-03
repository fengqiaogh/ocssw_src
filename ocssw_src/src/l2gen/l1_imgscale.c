#include "l12_proto.h"
#include "input_struc.h"
#define INT32   int32_t 
#define FLOAT32 float
#define BYTE    unsigned char

/* -------------------------------------------------------------- */
float toa_reflect(l1str *l1rec, int32_t ip, int32_t ib) {
    int32_t ipb;
    float mu0;
    static float pi = OEL_PI;
    float rho;

    mu0 = l1rec->csolz[ip];
    ipb = ip * l1rec->l1file->nbands + ib;
    rho = l1rec->Lt[ipb] * pi / l1rec->Fo[ib] / mu0;
    return (rho);
}


/* -------------------------------------------------------------- */
/* linscale() - linear reflectance scaling                        */

/* -------------------------------------------------------------- */
BYTE linscale(FLOAT32 x, FLOAT32 minrad, FLOAT32 maxrad) {

    if (x >= maxrad || x <= 0)
        return (255);

    if (x <= minrad)
        return (0);

    return (BYTE) ((x - minrad)*255.0 / maxrad);
}


/* -------------------------------------------------------------- */
/* linscale() - linear scaling                                    */

/* -------------------------------------------------------------- */
BYTE linscale_old(FLOAT32 x, int band) {
    float maxrad = 15.0;
    float minrad = 0.0;

    switch (band) {
    case 0: minrad = 0;
        maxrad = 15;
        break;
    case 4: minrad = 0;
        maxrad = 12;
        break;
    case 5: minrad = 0;
        maxrad = 10;
        break;
        /*
            case 2: minrad=0; maxrad=30; break;
            case 4: minrad=0; maxrad=24; break;
            case 5: minrad=0; maxrad=20; break;
         */
    default: break;
    }

    if (x >= maxrad || x <= 0)
        return (250);

    if (x <= minrad)
        return (0);

    return (BYTE) ((x - minrad)*250.0 / maxrad);
}

/* -------------------------------------------------------------- */
/* scale() - log scaling of reflectance                           */

/* -------------------------------------------------------------- */
BYTE logscale(FLOAT32 val, FLOAT32 min, FLOAT32 max) {
    int32_t scaleval;
    /* Save a few calls to the log function. */
    static double logmaxmin = 1, lastmin = 0, lastmax = 0;

    if (min != lastmin || max != lastmax) {

        logmaxmin = log(max / min);
        lastmin = min;
        lastmax = max;
    }

    if (val > min)
        /*scaleval = ((log10(val) + 2.0)/(2.0/255));*/
        scaleval = rint(255 * log(val / min) / logmaxmin);
    else
        scaleval = 0;
    if (scaleval < 0) return ( 0);
    if (scaleval > 255) return (255);
    return ((BYTE) scaleval);
}
