#include "l12_proto.h"
#include <gsl/gsl_fit.h>

// ============================================================================
// Compute refractive index of seawater given temperature and salinity
// ============================================================================

float seawater_nsw(float wave, float sst, float sss, float *dnswds) {
    float n0 = 1.31405;
    float n1 = 1.779e-4;
    float n2 = -1.05e-6;
    float n3 = 1.6e-8;
    float n4 = -2.02e-6;
    float n5 = 15.868;
    float n6 = 0.01155;
    float n7 = -0.00423;
    float n8 = -4382.0;
    float n9 = 1.1455e6;

    float n_air;
    float n_sw;
    float S;
    float T;
    float T2;
    float wave2;
    float wave3;

    S = sss;
    T = sst;
    T2 = T*T;
    wave2 = wave*wave;
    wave3 = wave2*wave;

    // refractive index of air is from Ciddor (1996, Applied Optics)
    n_air = 1.0 + (5792105.0 / (238.0185 - 1 / (wave2 / 1e6)) + 167917.0 / (57.362 - 1 / (wave2 / 1e6))) / 1e8;

    // refractive index of seawater is from Quan and Fry (1994, Applied Optics)
    n_sw = (n0 + (n1 + n2 * T + n3 * T2) * S + n4 * T2 + (n5 + n6 * S + n7 * T) / wave + n8 / wave2 + n9 / wave3) * n_air;

    if (dnswds != NULL) {
        *dnswds = (n1 + n2 * T + n3 * T2 + n6 / wave) * n_air;
    }

    return (n_sw);
}

// ============================================================================
// Compute seawater isothermal compressibility
// ============================================================================

float seawater_betat(float sst, float sss) {
    float S = sss;
    float T = sst;
    float T2, T3, T4;
    float kw;
    float a0, b0, Ks;
    float betat;

    T2 = T*T;
    T3 = T2*T;
    T4 = T2*T2;

    // pure water secant bulk Millero (1980, Deep-sea Research)

    kw = 19652.21 + 148.4206 * T - 2.327105 * T2 + 1.360477e-2 * T3 - 5.155288e-5 * T4;

    // seawater secant bulk

    a0 = 54.6746 - 0.603459 * T + 1.09987e-2 * T2 - 6.167e-5 * T3;
    b0 = 7.944e-2 + 1.6483e-2 * T - 5.3009e-4 * T2;

    Ks = kw + a0 * S + b0 * S*sqrt(S);

    // calculate seawater isothermal compressibility from the secant bulk

    betat = 1 / Ks * 1e-5; // unit is Pa

    return (betat);
}


// ============================================================================
// Compute density of seawater, unit is kg/m3, from UNESCO,38,1981
// ============================================================================

float seawater_density(float sst, float sss) {
    float S = sss;
    float T = sst;
    float T2, T3, T4, T5;
    float density_w, density_sw;

    float a0 = 8.24493e-1;
    float a1 = -4.0899e-3;
    float a2 = 7.6438e-5;
    float a3 = -8.2467e-7;
    float a4 = 5.3875e-9;
    float a5 = -5.72466e-3;
    float a6 = 1.0227e-4;
    float a7 = -1.6546e-6;
    float a8 = 4.8314e-4;
    float b0 = 999.842594;
    float b1 = 6.793952e-2;
    float b2 = -9.09529e-3;
    float b3 = 1.001685e-4;
    float b4 = -1.120083e-6;
    float b5 = 6.536332e-9;

    T2 = T*T;
    T3 = T2*T;
    T4 = T2*T2;
    T5 = T3*T2;

    // density for pure water 

    density_w = b0 + b1 * T + b2 * T2 + b3 * T3 + b4 * T4 + b5*T5;

    // density for pure seawater

    density_sw = density_w + ((a0 + a1 * T + a2 * T2 + a3 * T3 + a4 * T4) * S + (a5 + a6 * T + a7 * T2) * S*sqrt(S) + a8 * S * S);

    return (density_sw);
}


// ============================================================================
// Compute activity wrt salinity
// ============================================================================

float seawater_dlnaswds(float sst, float sss) {
    float S = sss;
    float T = sst;
    float T2, T3;
    float dlnaswds;

    T2 = T*T;
    T3 = T2*T;

    // water activity data of seawater is from Millero and Leung (1976,American
    // Journal of Science,276,1035-1077). Table 19 was reproduced using
    // Eqs.(14,22,23,88,107) then were fitted to polynominal equation.
    // dlnawds is partial derivative of natural logarithm of water activity
    // w.r.t.salinity

    dlnaswds = (-5.58651e-4 + 2.40452e-7 * T - 3.12165e-9 * T2 + 2.40808e-11 * T3) +
            1.5 * (1.79613e-5 - 9.9422e-8 * T + 2.08919e-9 * T2 - 1.39872e-11 * T3) * sqrt(S) +
            2 * (-2.31065e-6 - 1.37674e-9 * T - 1.93316e-11 * T2) * S;

    return (dlnaswds);
}


// ============================================================================
// Compute derivative of nsw wrt density via PMH model
// ============================================================================

float seawater_dnswdrho(float n_sw) {
    float n_sw2;
    float dnswdrho;

    n_sw2 = n_sw*n_sw;

    dnswdrho = (n_sw2 - 1)*(1 + 2. / 3. * (n_sw2 + 2) * pow((n_sw / 3 - 1. / 3. / n_sw), 2));

    return (dnswdrho);
}



// ============================================================================
// Compute backscattering coeff for seawater at temperature and salinity
//
// Ref:  X. Zhang, L. Hu, and M-X. He, "Scattering by pure seawater: Effect of salinity," 
//       Optics Express 17(7), 5698-5710 (2009).
// ============================================================================

float seawater_bb(float wave, float sst, float sss, double delta) {
    double S = sss;
    double T = sst;

    // constants

    double Na = 6.0221417930e23; // Avagadro
    double Kbz = 1.3806503e-23; // Boltzmann
    double M0 = 18e-3; // molecular weight of water in kg/mol
    // double delta = 0.039; // depolarization ratio from Farinato and Roswell 1975

    double Tk;
    float n_sw, dnswds;
    float isocomp;
    float density_sw;
    float dlnawds;
    float dfri;
    float beta_df;
    float flu_con;
    float beta_cf;
    float beta90sw;
    float bbsw;

    Tk = T + 273.15;

    // n_sw: absolute refractive index of seawater
    // dnswds: partial derivative of seawater refractive index w.r.t. salinity

    n_sw = seawater_nsw(wave, T, S, &dnswds);

    // isothermal compressibility is from Lepple & Millero 1971
    // the error ~ 0.004e-6 bar-1

    isocomp = seawater_betat(T, S);

    // density of water and seawater, unit is kg/m3, from UNESCO,38,1981

    density_sw = seawater_density(T, S);
    // derivative of natural logarithm of water activity w.r.t.salinity

    dlnawds = seawater_dlnaswds(T, S);

    // density derivative of refractive index from PMH model

    dfri = seawater_dnswdrho(n_sw);

    // volume scattering at 90 degree due to the density fluctuation
    float power_wave = 1/(wave * 1e-9);
    power_wave*=power_wave;
    power_wave*=power_wave;
    beta_df = OEL_PI * OEL_PI / 2 * power_wave * Kbz * Tk * isocomp * dfri * dfri * (6 + 6 * delta) / (6 - 7 * delta);

    // volume scattering at 90 degree due to the concentration fluctuation

    flu_con = S * M0 * dnswds * dnswds / density_sw / (-dlnawds) / Na;
    beta_cf = 2 * OEL_PI * OEL_PI * power_wave * n_sw * n_sw * (flu_con)*(6 + 6 * delta) / (6 - 7 * delta);

    // total volume scattering at 90 degrees

    beta90sw = beta_df + beta_cf;

    // total backscattering coefficient

    bbsw = 8 * OEL_PI / 3 * beta90sw * (2 + delta) / (1 + delta) * 0.5;

    return (bbsw);
}

/**
 * function to return spectral slope parameter for bbw
 *
 * @param l2rec l2rec structure
 * @param p product sturcture
 * @param prod output bbw slope parameter
 *
 * The slope is computed on the log transformed data
 * for the visible wavelengths only.  This works,
 * because bbw can be approximated by a power function.
 *
 */
void get_bbws(l2str *l2rec, l2prodstr *p, float prod[]) {
    int i, ip;
    int32_t nbands = l2rec->l1rec->l1file->nbands;
    float *fwave = l2rec->l1rec->l1file->fwave;

    int nvisbands = 0;

    double x[nbands];
    double y[nbands];

    for (i = 0; i < nbands; i++)
        if (fwave[i] > 400 && fwave[i] < 700.) {
            x[i] = log(fwave[i]);
            nvisbands++;
        }

    for (ip = 0; ip < l2rec->l1rec->npix; ip++) {
        prod[ip] = p->badData;

        double c0, c1, cov00, cov01, cov11, sumsq;
        for (i = 0; i < nvisbands; i++) {
            int ibp = ip * nbands + i;
            y[i] = log(l2rec->l1rec->sw_bb[ibp]);
        }
        //        gsl_fit_mul(x,1,y,1,nbands, &c1, &cov11,&sumsq);
        gsl_fit_linear(x, 1, y, 1, nvisbands,
                &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

        prod[ip] = (float) c1;
    }
}
