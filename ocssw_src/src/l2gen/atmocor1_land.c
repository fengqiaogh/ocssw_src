/* ================================================================ */
/* atmocor1_land - computes atmospheric components over land.       */
/*                                                                  */
/* Written By: B. A. Franz, SAIC GSC, NASA/SIMBIOS, April 2000      */
/*     Based on code provided by Jacques Descloitres                */
/*                                                                   */
/* Revised by: M.Zhang, December 2023                                */
/* add in a trace gas absorption correction for all species corrected */
/* for over ocean.                                                   */
/*                                                                  */
/* ================================================================ */

#include <math.h>
#include "l12_proto.h"


void chand(float xphi, float xmuv, float xmus, float *xtau,
        float *rhoray, double *trup, double *trdown, int nbands);

int atmocor1_land(l1str *l1rec, int32_t ip) {
    static int32_t firstRun = 1;
    static float p0 = 1013.25;

    float delphi = l1rec->delphi[ip] + 180.0;
    float mu0 = l1rec->csolz[ip];
    float mu = l1rec->csenz[ip];
    float airmass = 1.0 / mu0 + 1.0 / mu;
    int nbands = l1rec->l1file->nbands;

    static float *Taur;
    static float *rhor;
    static double *trup;
    static double *trdown;

    int32_t ib, ipb;
    if (firstRun) {
        firstRun = 0;
        if ((Taur = (float *) calloc(nbands, sizeof (float))) == NULL) {
            printf("-E- : Error allocating memory to Taur\n");
            exit(FATAL_ERROR);
        }
        if ((rhor = (float *) calloc(nbands, sizeof (float))) == NULL) {
            printf("-E- : Error allocating memory to rhor\n");
            exit(FATAL_ERROR);
        }
        if ((trup = (double *) calloc(nbands, sizeof (double))) == NULL) {
            printf("-E- : Error allocating memory to trup\n");
            exit(FATAL_ERROR);
        }
        if ((trdown = (double *) calloc(nbands, sizeof (double))) == NULL) {
            printf("-E- : Error allocating memory to trdown\n");
            exit(FATAL_ERROR);
        }
    }

    /*                                                              */
    /* First, correct Rayleigh optical thickness for pressure.      */
    /*                                                              */
    for (ib = 0; ib < nbands; ib++) {
        Taur[ib] = l1rec->l1file->Tau_r[ib] * l1rec->pr[ip] / p0;
    }

    /*                                                              */
    /* Compute Rayleigh path reflectances, also total diff. trans.  */
    /*                                                              */
    chand(delphi, mu, mu0, Taur, rhor, trup, trdown, nbands);

    /* Compute gaseous transmittance functions */
    gaseous_transmittance(l1rec, ip);

    /*                                                              */
    /* Loop through bands and compute transmittance and reflectance */
    /*                                                              */
    for (ib = 0; ib < l1rec->l1file->nbands; ib++) {

        ipb = ip * l1rec->l1file->nbands + ib;

        /* Compute and store Rayleigh radiance and diff trans */
        l1rec->Lr[ipb] = rhor[ib] * l1rec->Fo[ib] * mu0 / OEL_PI;
        l1rec->t_sol[ipb] =
                ((2.0 / 3.0 + mu0) + (2.0 / 3.0 - mu0) * trdown[ib])
                / (4.0 / 3.0 + Taur[ib]);
        l1rec->t_sen[ipb] =
                ((2.0 / 3.0 + mu) + (2.0 / 3.0 - mu) * trup [ib])
                / (4.0 / 3.0 + Taur[ib]);

        /* no way to comput polarization correction over land */
        l1rec->polcor[ipb] = 1.0;

        if (l1rec->l1file->sensorID == SEAWIFS && ((input->gas_opt & GAS_TRANS_TBL_BIT) == 0)) {
                l1rec->t_h2o[ipb] = water_vapor(ib, l1rec->wv[ip], airmass);
        }

    }

    return (0);
}


/* ------------------------------------------------------------------------ */
/* Code below was provided by Jacques Descloitres, March 2000               */
/* ------------------------------------------------------------------------ */

#define fac OEL_DEGRAD  /* PI/180 */

void chand(float xphi, float xmuv, float xmus, float *xtau, float *rhoray, double *trup, double *trdown, int nbands) {
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

    chAnged all instances of expf, logf, and cosf to exp, log, cos, to fix problems
    with 64-bit solaris systems. BAF, 8/2002.
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

    phios = 180 - xphi;
    xcosf1 = 1.;
    xcosf2 = cos(phios * fac);
    xcosf3 = cos(2 * phios * fac);
    xph1 = 1 + (3 * xmus * xmus - 1) * (3 * xmuv * xmuv - 1) * xfd / 8.;
    xph2 = -xfd * xbeta2 * 1.5 * xmus * xmuv * sqrt(1 - xmus * xmus) * sqrt(1 - xmuv * xmuv);
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
        xlntau = log(xtau[ib]);
        fs0 = fs01 + fs02 * xlntau;
        fs1 = as1[0] + xlntau * as1[1];
        fs2 = as2[0] + xlntau * as2[1];
        trdown[ib] = exp(-xtau[ib] / xmus);
        trup[ib] = exp(-xtau[ib] / xmuv);
        xitm1 = (1 - trdown[ib] * trup[ib]) / 4. / (xmus + xmuv);
        xitm2 = (1 - trdown[ib]) * (1 - trup[ib]);
        xitot1 = xph1 * (xitm1 + xitm2 * fs0);
        xitot2 = xph2 * (xitm1 + xitm2 * fs1);
        xitot3 = xph3 * (xitm1 + xitm2 * fs2);
        rhoray[ib] = xitot1 * xcosf1 + xitot2 * xcosf2 * 2 + xitot3 * xcosf3 * 2;
    }

}
