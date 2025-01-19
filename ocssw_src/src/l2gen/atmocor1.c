#include "l12_proto.h"
#include "glint.h"
/* ========================================================================== */
/* module atmocor1() - computes pre-aerosol atmospheric components            */
/*                                                                            */
/* Written By: B. A. Franz, SAIC, SeaWiFS Project, February 1998              */
/* Conversion to C, May 2006                                                  */
/* ========================================================================== */


/* ------------------------------------------------------------------- */
/* correction factor to remove oxygen absorption from La(765)          */

/* ------------------------------------------------------------------- */
float oxygen_aer(float airmass) {
    /* base case: m80_t50_strato  visibility (550nm ,0-2km):25km */
    static float a[] = {-1.0796, 9.0481e-2, -6.8452e-3};
    return (1.0 + pow(10.0, a[0] + airmass * a[1] + airmass * airmass * a[2]));
}


/* ------------------------------------------------------------------- */
/* correction factor to replace oxygen absorption to Lr(765)           */

/* ------------------------------------------------------------------- */
float oxygen_ray(float airmass) {
    /* base case is the 1976 Standard atmosphere without aerosols */
    static float a[] = {-1.3491, 0.1155, -7.0218e-3};
    return (1.0 / (1.0 + pow(10.0, a[0] + airmass * a[1] + airmass * airmass * a[2])));
}


/* ------------------------------------------------------------------- */
/* main function, loads various quanitities into l1 record for 1 pixel */

/* ------------------------------------------------------------------- */

void atmocor1(l1str *l1rec, int32_t ip) {
    static float p0 = STDPR;

    int32_t sensorID = l1rec->l1file->sensorID;
    int32_t nwave = l1rec->l1file->nbands;

    float solz = l1rec->solz[ip];
    float senz = l1rec->senz[ip];
    float raz = l1rec->delphi[ip];
    float mu0 = l1rec->csolz[ip];
    float mu = l1rec->csenz[ip];

    float ws = l1rec->ws[ip];
    float pr = l1rec->pr[ip];
    float wv = l1rec->wv[ip];

    float *Fo = l1rec->Fo;
    float *Tau_r = l1rec->l1file->Tau_r;

    float zero = 0.0;
    int32_t ib765 = -1;
    float airmass;
    float a_o2;
    float glint_coef_q;
    float glint_coef_u;
    int32_t ib, ipb;
    float scaleRayleigh;

    airmass = 1.0 / mu0 + 1.0 / mu;
    ipb = ip*nwave;

    /* Initialize output values */
    if (sensorID == SEAWIFS || sensorID == OCTS) {
        for (ib = 0; ib < nwave; ib++) {
            // Explicitly only doing the Ding and Gordon correction if the sensor has a band AT 765nm
            // So, you know, SeaWiFS, OCTS, OSMI, and any others to be named later...
            if (input->oxaband_opt == 1 && l1rec->l1file->iwave[ib] == 765) {
                ib765 = ib;
                l1rec->t_o2[ipb + ib765] = 1.0 / oxygen_aer(airmass);
            }

            if (sensorID == SEAWIFS && ((input->gas_opt & GAS_TRANS_TBL_BIT) == 0)) {
                /* this is for rhos only, but effects seawifs land products and cloud   */
                /* will modify get_rhos() to use gaseous transmittance in future update */
                l1rec->t_h2o[ipb + ib] = water_vapor(ib, wv, airmass);
            }
        }
    }

    /* apply gaseous transmittance */
    gaseous_transmittance(l1rec, ip);

    /* white-cap radiances at TOA */
    whitecaps(sensorID, l1_input->evalmask, nwave, ws, input->wsmax, &l1rec->rhof[ipb]);

    for (ib = 0; ib < nwave; ib++) {
        l1rec->t_sol [ipb + ib] = exp(-0.5 * pr / p0 * Tau_r[ib] / mu0);
        l1rec->t_sen [ipb + ib] = exp(-0.5 * pr / p0 * Tau_r[ib] / mu);
        l1rec->tLf[ipb + ib] = l1rec->rhof[ipb + ib] * l1rec->t_sen[ipb + ib] * l1rec->t_sol[ipb + ib] * Fo[ib] * mu0 / M_PI;
    }

    /* Rayleigh scattering */
    if (sensorID != AVHRR && sensorID != OCRVC) {
        rayleigh(l1rec, ip);
    }

    // If we are running with Ding and Gordon, need to do the rayleigh part as well...
    if (input->oxaband_opt == 1 && ib765 > -1) {
        a_o2 = oxygen_ray(airmass);
        l1rec->Lr [ipb + ib765] *= a_o2;
        l1rec->L_q[ipb + ib765] *= a_o2;
        l1rec->L_u[ipb + ib765] *= a_o2;
    }

    //Scale by the altitude of the sensor, assuming height of atmosphere=100km
    if (l1rec->alt > 0) {

        scaleRayleigh = 1.0 - exp(-l1rec->alt / 10); // Assume 10km is e-folding height

        for (ib = 0; ib < nwave; ib++) {
            l1rec->Lr [ipb + ib] *= scaleRayleigh;
            l1rec->L_q[ipb + ib] *= scaleRayleigh;
            l1rec->L_u[ipb + ib] *= scaleRayleigh;
        }
    }

    /* glint coefficients and approximate glint radiances */
    /* also add glint to polarization components          */

    /* for avhrr, l1_aci_hdf.c calls avhrrsub5h.f which calls getglint */
    if (sensorID != AVHRR) {

        getglint_iqu(senz, solz, raz, ws, zero,
                &l1rec->glint_coef[ip], &glint_coef_q, &glint_coef_u);

        // getglint_iqu_(&senz, &solz, &raz, &ws, &zero,
        //         &l1rec->glint_coef[ip], &glint_coef_q, &glint_coef_u);

        for (ib = 0; ib < nwave; ib++) {
            l1rec->TLg[ipb + ib] = l1rec->glint_coef[ip] * exp(-(Tau_r[ib] + 0.1) * airmass) * Fo[ib];
            l1rec->L_q[ipb + ib] += (glint_coef_q * l1rec->TLg[ipb + ib]);
            l1rec->L_u[ipb + ib] += (glint_coef_u * l1rec->TLg[ipb + ib]);
        }
    }

}