/* ================================================================ */
/* get_rhos() - computes surface reflectance.                       */
/*                                                                  */
/* Note to users:    This surface reflectance is used               */
/*  by some derived product algorithms which may not require        */
/*  full atmospheric correction.  The same calculations are         */
/*  repeated for standard atmospheric correction over ocean.        */
/*                                                                  */
/* Written By: B. A. Franz, SAIC GSC, NASA/SIMBIOS, April 2000      */
/*     Based on code provided by Jacques Descloitres                */
/* W. Robinson, SAIC, 7 Nov 2012, for radiance < 0, return          */
/*     reflectance of BAD_FLT                                       */
/*                                                                  */
/* ================================================================ */

#include <math.h>
#include "l12_proto.h"

double fintexp3(float taur);
double fintexp1(float taur);
float csalbr(float taur);

int get_rhos(l1str *l1rec, int32_t ip) {
    static float pi = 3.141592654;
    static float radeg = 180. / 3.141592654;
    static float p0 = 1013.25;
    static int firstCall = TRUE;
    static float sphalb[5000];

    float mu0 = cos(l1rec->solz[ip] / radeg);

    float Ka = 0.8;
    float Taur;
    int32_t i, ib, ipb;
    int32_t cirrus_opt = input->cirrus_opt;

    /*                                                              */
    /* Derive spherical albedo correction table                     */
    /*                                                              */
    if (firstCall) {
        firstCall = FALSE;
        sphalb[0] = 0.0;
        for (i = 1; i < 5000; i++)
            sphalb[i] = csalbr(i / 10000.);
    }
    /*                                                              */
    /* Loop through bands and compute transmittance and reflectance */
    /*                                                              */
    for (ib = 0; ib < l1rec->l1file->nbands; ib++) {

        ipb = ip * l1rec->l1file->nbands + ib;

        /* Compute surface reflectance */
        if (l1rec->Lt[ipb] > 0.) {
            if (cirrus_opt)
                l1rec->rhos[ipb] = pi / l1rec->Fo[ib] / mu0
                    * ((l1rec->Lt[ipb] / l1rec->tg_sol[ipb] / l1rec->tg_sen[ipb]
                    - l1rec->Lr[ipb]) - l1rec->rho_cirrus[ip] / Ka)
                / l1rec->t_sol[ipb] / l1rec->t_sen[ipb] / l1rec->t_o2[ipb] / l1rec->t_h2o[ipb];
            else
                l1rec->rhos[ipb] = pi / l1rec->Fo[ib] / mu0
                    * (l1rec->Lt[ipb] / l1rec->tg_sol[ipb] / l1rec->tg_sen[ipb]
                    - l1rec->Lr[ipb])
                / l1rec->t_sol[ipb] / l1rec->t_sen[ipb] / l1rec->t_o2[ipb] / l1rec->t_h2o[ipb];

            /* Correct for spherical albedo (atmosphere-surface interaction).*/
            /* Note that for Ocean processing, these effects are already     */
            /* included in the Rayleigh radiances, based on a priori         */
            /* knowledge of the ocean surface reflectance as fn of windspeed.*/
            if (l1rec->land[ip]) {
                Taur = l1rec->l1file->Tau_r[ib] * l1rec->pr[ip] / p0;
                l1rec->rhos[ipb] /=
                        (1 + sphalb[(int) (Taur * 10000 + 0.5)] * l1rec->rhos[ipb]);
            }
        } else
            l1rec->rhos[ipb] = BAD_FLT;
    }
    return (0);
}


/* ------------------------------------------------------------------------ */
/* Code below was provided by Jacques Descloitres, March 2000               */
/*                                                                          */
/* changed expf, logf to exp, log, to fix problems with 64-bit solaris      */
/* BAF, 8/2002.                                                             */

/* ------------------------------------------------------------------------ */

float csalbr(float xtau) {
    return (3 * xtau - fintexp3(xtau)*(4 + 2 * xtau) + 2 * exp(-xtau)) / (4 + 3 * xtau);
}

double fintexp3(float xtau) {
    return (exp(-xtau)*(1. - xtau) + xtau * xtau * fintexp1(xtau)) / 2.;
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
    return xx - log(xtau);
}

