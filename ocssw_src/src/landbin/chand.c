#include <math.h>

#define fac 0.0174532925199  /* PI/180 */

void chand(float xphi, float xmuv, float xmus, float *xtau, float *rhoray, double *trup, double *trdown, int nbands, char *process) {
    /*
    input parameters: xphi,xmus,xmuv,xtau
    xphi: azimuthal difference between sun and observation (xphi=0,
          in backscattering) and expressed in degree (0.:360.)
    xmus: cosine of the sun zenith angle
    xmuv: cosine of the observation zenith angle
    xtau: molecular optical depth
    output parameter: xrray : molecular reflectance (0.:1.)
    constant : xdep: depolarization factor (0.0279)
               xfd = (1-xdep/(2-xdep)) / (1 + 2*xdep/(2-xdep)) = 2 * (1 - xdep) / (2 + xdep) = 0.958725775
     */
    const double xfd = 0.958725775;
    const float xbeta2 = 0.5;
    float pl[5];
    double fs01, fs02, fs0, fs1, fs2;
    const float as0[10] = {0.33243832, 0.16285370, -0.30924818, -0.10324388, 0.11493334,
        -6.777104e-02, 1.577425e-03, -1.240906e-02, 3.241678e-02, -3.503695e-02};
    const float as1[2] = {.19666292, -5.439061e-02};
    const float as2[2] = {.14545937, -2.910845e-02};
    float phios, xcosf1, xcosf2, xcosf3;
    float xph1, xph2, xph3, xitm1, xitm2;
    float xlntau, xitot1, xitot2, xitot3;
    int i, ib;

    phios = xphi + 180;
    xcosf1 = 1.;
    xcosf2 = cosf(phios * fac);
    xcosf3 = cosf(2 * phios * fac);
    xph1 = 1 + (3 * xmus * xmus - 1) * (3 * xmuv * xmuv - 1) * xfd / 8.;
    xph2 = -xfd * xbeta2 * 1.5 * xmus * xmuv * sqrtf(1 - xmus * xmus) * sqrtf(1 - xmuv * xmuv);
    xph3 = xfd * xbeta2 * 0.375 * (1 - xmus * xmus) * (1 - xmuv * xmuv);
    pl[0] = 1.;
    pl[1] = xmus + xmuv;
    pl[2] = xmus * xmuv;
    pl[3] = xmus * xmus + xmuv * xmuv;
    pl[4] = xmus * xmus * xmuv * xmuv;
    fs01 = fs02 = 0;
    for (i = 0; i < 5; i++) fs01 += pl[i] * as0[i];
    for (i = 0; i < 5; i++) fs02 += pl[i] * as0[5 + i];
    for (ib = 0; ib < nbands; ib++) {
        if (process[ib]) {
            xlntau = logf(xtau[ib]);
            fs0 = fs01 + fs02 * xlntau;
            fs1 = as1[0] + xlntau * as1[1];
            fs2 = as2[0] + xlntau * as2[1];
            trdown[ib] = expf(-xtau[ib] / xmus);
            trup[ib] = expf(-xtau[ib] / xmuv);
            xitm1 = (1 - trdown[ib] * trup[ib]) / 4. / (xmus + xmuv);
            xitm2 = (1 - trdown[ib]) * (1 - trup[ib]);
            xitot1 = xph1 * (xitm1 + xitm2 * fs0);
            xitot2 = xph2 * (xitm1 + xitm2 * fs1);
            xitot3 = xph3 * (xitm1 + xitm2 * fs2);
            rhoray[ib] = xitot1 * xcosf1 + xitot2 * xcosf2 * 2 + xitot3 * xcosf3 * 2;
        }
    }

}
