#include "l12_proto.h"
#include "pml_iop.h"
#include "pml_iop_config.h" 
#include "pml_iop_tables.h"
#include "pml_iop_calculate.h"

/* Module: pml_iop_calculate */
/* Authors: Tim Smyth and Gerald Moore */
/* Date: 13/2/2006 */
/* Version: 1.0 */
/* Description */

/* Function: iop_model */
int iop_model(double rho_w[], float sun_theta, float sen_theta, float dphi, double a[],
        double bbp[], double ady[], double ap[], int MODIS, int CASEII) {

    int i, result, lower = 0, upper = 1;
    int iop_model_error = 0; /* successful iop calculation */
    float eps_a, ady_0;
    float adyu, adyl, apl, au, al, ysbpa_sc;
    float x1, x2;

    eps_a = eps_a_init;

    /* If the iw531 exists (set to 3) then use MODIS scalings */
    if (MODIS > 0)
        eps_a = eps_a_init_modis;

    /* calculate the total absorption and backscatter */
    /* only do calculation if rho_w values are sensible */

    if (rho_w[1] > 0. && rho_w[2] > 0. && rho_w[3] > 0. &&
            rho_w[4] > 0.) {
        result = mod_iter(rho_w, sun_theta, sen_theta, dphi, eps_a, a, bbp, MODIS, CASEII);
        if (result != 0) {
            iop_model_error = 1; /* unsuccessful primary IOP calculation */
            for (i = 0; i < NB; i++) {
                a[i] = 0.;
                bbp[i] = 0.;
            }
            return iop_model_error;
        }

        /* using the spectral slopes between 412 and 443 work out ady */
        /* only calculate if a_412 and a_443 are greater than zero */
        if (a[0] > 0. && a[1] > 0. && a[2] > 0.) {

            /* new formulation to account for non-linear behavior */
            /* of spectral slope in ady */
            if (ysbpa_l == 443.0) {
                lower = 1; /* 443.0 */
                upper = 2; /* 490.0 */
            }

            al = a[lower];
            au = a[upper];

            result = biogeochem_iter(al, au, &adyu, &adyl, &apl);
            if (result != 0) {
                iop_model_error = 2;
                for (i = 0; i < NB; i++) {
                    ady[i] = 0.;
                    ap[i] = 0.;
                }
                return iop_model_error;
            }

            ady[lower] = adyl;
            ady[upper] = adyu;
            ap[lower] = apl;

            if (ysbpa_l == 443.0)
                ady_0 = ady[lower];
            else
                ady_0 = ady[upper];

            /* calculate the spectral slope of CDOM from info at two shorter wavelengths */
            ysbpa_sc = log(ady[lower] / ady[upper]) / (lambda[lower] - lambda[upper]);

            for (i = 0; i < NB; i++) {
                /* Apply gelbstoff slope and intercept (at 443) to give entire spectra */
                ady[i] = ady_0 * exp(ysbpa_sc * (lambda[i] - ysbpa_0));
                /* Calculate spectral ap (absorption due to pigments) using the remainder */
                ap[i] = a[i] - ady[i];
            }

            /* check aph443 range (this sanity check is an empirical fix from the Lee model */
            x1 = ap[1] / a[1];
            if (x1 < 0.15 || x1 > 0.6) {
                /* ratio is between 443 and 412 */
                x2 = -0.8 + 1.4 * (a[1] / a[0]);
                if (x2 < 0.15)
                    x2 = 0.15;
                if (x2 > 0.6)
                    x2 = 0.6;
                ap[1] = a[1] * x2;
                ady_0 = a[1] - ap[1];

                for (i = 0; i < NB; i++) {
                    /* Use the preset slope for gelbstoff */
                    ady[i] = ady_0 * exp(ysbpa_s * (lambda[i] - ysbpa_0));
                    ap[i] = a[i] - ady[i];
                }
            }
            /* End of the empirical fix */
        } else {
            for (i = 0; i < NB; i++) {
                ady[i] = 0.;
                ap[i] = 0.;
            }
            iop_model_error = 2; /* unsuccessful bio-geochemical IOP calculation */
        }
    }

    return iop_model_error;
}

/* Function: mod_iter */
int mod_iter(double rho_w[], float sun_theta, float sen_theta, float dphi, float eps_a, double a[], double bb[], int MODIS, int CASEII) {
    int i, iter, fail;
    int mod_iter_init = 0;
    double b[NB], FC[NB], F[NB];
    double rho_wT[2], awT[2], bbwT[2], FT[2], df;
    double ab[2];
    double temp, eps_b;
    float bbw[NB]; /* backscatter due to pure water */
    float init_a[NB], init_b[NB], bn[NB]; /* Initial guess at IOP values */

    /* Modification for MODIS - reference scattering wavelength */
    if (MODIS > 0)
        scat_l = scat_l_modis;

    /* Initial guess for scattering and absorption */
    /* (based on init_chl) */
    if (mod_iter_init == 0) {
        temp = scat_a + scat_b * lambda[bp[1]] + scat_c * pow(lambda[bp[1]] / scat_l, scat_n);
        for (i = 0; i < NB; i++) {
            init_a[i] = geo2iop(ch_lev, ac, i, init_chl, ch_n);
            init_b[i] = geo2iop(ch_lev, bc, i, init_chl, ch_n);
            bbw[i] = b_w[i] * b_tilde_w;
            bn[i] = scat_a + scat_b * lambda[bp[1]] + scat_c * pow(lambda[i] / scat_l, scat_n);
            bn[i] = bn[i] / temp;
        }
        /* This comes out with spectrally flat eps_b */
        eps_b = bn[bp[0]];

        /* This is a flag raised once the iterations have been initialised */
        mod_iter_init = 1;
    }

    if (NFLAG == 1) {
        /* Re-evaluate scattering slope */
        /* Set here to scat_n - but originally a user defined input 'n' */
        /* Re-evaluate if neccessary */
        temp = pow(lambda[bp[1]] / scat_l, scat_n);
        for (i = 0; i < NB; i++) {
            bn[i] = pow(lambda[i] / scat_l, scat_n);
            bn[i] = bn[i] / temp;
        }
        /* The result of this is eps_b of 1.0202 for SeaWiFS (490:510) */
        eps_b = bn[bp[0]];
    }

    /* Set the relevant geometry */
    setgeom(sun_theta, sen_theta, dphi);
    /* Prepare the iterations */
    for (i = 0; i < NB; i++) {
        /* set the reflectance to zero if it is negative */
        if (rho_w[i] < 0.0) rho_w[i] = 0.0;
        a[i] = init_a[i];
        b[i] = init_b[i];
        bb[i] = b[i] * b_tilde_p;

        if (CASEII)
            FC[i] = f_ab(init_a[i], init_b[i], i) / (1. + bb[i] / a[i]);
        else
            FC[i] = f_ab(init_a[i], init_b[i], i);

    }

    iter = fail = 0;

    /* Iteration loop */
    do {
        /* set the F(=pi*R*(f/Q)) for each band - */
        /* this is reset during the iterations */
        for (i = 0; i < NB; i++)
            F[i] = FC[i];

        /* temporary arrays*/
        rho_wT[0] = rho_w[bp[0]];
        rho_wT[1] = rho_w[bp[1]];
        awT[0] = a_w[bp[0]];
        awT[1] = a_w[bp[1]];
        bbwT[0] = bbw[bp[0]];
        bbwT[1] = bbw[bp[1]];
        FT[0] = F[bp[0]];
        FT[1] = F[bp[1]];

        /* Calculate new a(510 (or 531 MODIS)) and bb(510 (or 531 MODIS)) values */
        if (CASEII)
            iter_ab2(rho_wT, awT, bbwT, FT, eps_b, eps_a, ab);
        else
            iter_ab(rho_wT, awT, bbwT, FT, eps_b, eps_a, ab);

        /* Evaluate a490*/
        /* This is done via the empirically calculated slopes */
        a[bp[0]] = ab[0] * eps_a; /* 490 nm */
        a[bp[1]] = ab[0]; /* 510 nm (SeaWiFS) 531 (MODIS) */
        /* df = 0.; */
        for (i = 0; i < NB; i++) {
            bb[i] = ab[1] * bn[i]; /* spectral bb assuming slope */
            /* don't re-evaluate 490 and 510 (531) */
            if ((i != bp[0]) || (i != bp[1])) {
                if (rho_w[i] > 0.0) {

                    if (CASEII)
                        a[i] = iter_a2(rho_w[i], a_w[i], bbw[i], F[i], bb[i]);
                    else
                        a[i] = iter_a(rho_w[i], a_w[i], bbw[i], F[i], bb[i]);
                } else a[i] = 0.0;
            }

            /* Calculate new values of b for F */
            b[i] = bb[i] / b_tilde_p;

            /* Evaluate new F values */
            FC[i] = f_ab(a[i], b[i], i);
        }

        df = fabs(FC[bp[0]] - F[bp[0]]) + fabs(FC[bp[1]] - F[bp[1]]);

        /* Evaluate */
        iter++;
    } while ((iter <= maxit) && (df >= tol));

    if (iter > maxit)
        fail = 2;

    return fail;
}

/* Function: iter_ab */

/* Returns the new a and b values */
int iter_ab(double rho_w[], double aw[], double bbw[], double F[], double epsb, double epsa, double ab[]) {
    int result = 0;
    double x, y, z;
    double scale;

    x = F[0] * rho_w[1];
    y = epsb*x;
    z = epsa * F[1] * rho_w[0];
    scale = y - z;
    if (scale == 0.0) {
        ab[0] = -1.0;
        ab[1] = -1.0;
        result = 1;
        return result;
    }
    /* a510 (a531)  */
    ab[0] = (F[0] * F[1]*(bbw[1] * epsb - bbw[0]) + aw[0] * F[1] * rho_w[0] - aw[1] * y) / scale;
    /* bb510 (bb531) */
    ab[1] = (rho_w[0] * rho_w[1]*(aw[0] - aw[1] * epsa) - bbw[0] * x + bbw[1] * z) / scale;
    return result;
}

/* Function: iter_ab2 */

/* Returns the new a and b values for CASEII formulation */
int iter_ab2(double rho_w[], double aw[], double bbw[], double F[], double epsb, double epsa, double ab[]) {
    int result = 0;
    double y, z;
    double scale_bb, scale_a;

    z = epsa * (rho_w[0] * rho_w[1] - F[1] * rho_w[0]);
    y = epsb * (rho_w[1] * F[0] - rho_w[1] * rho_w[0]);
    scale_bb = z + y;

    z = epsa * rho_w[0]*(rho_w[1] - F[1]);
    y = epsb * rho_w[1]*(rho_w[0] - F[0]);
    scale_a = z + y;

    if (scale_a == 0.0 || scale_bb == 0.0) {
        ab[0] = -1.0;
        ab[1] = -1.0;
        result = 1;
        return result;
    }
    /* a510 (a531) */
    ab[0] = (epsb * (rho_w[0] - F[0])*(F[1] * bbw[1] - rho_w[1]*(aw[1] + bbw[1]))-(rho_w[1] - F[1])*(F[0] * bbw[0] - rho_w[0]*(aw[0] + bbw[0]))) / scale_a;

    /* bb510 (bb531) */
    ab[1] = (rho_w[0]*(epsa * (F[1] * bbw[1] - rho_w[1]*(aw[1] + bbw[1])) + rho_w[1]*(aw[0] + bbw[0])) - rho_w[1] * F[0] * bbw[0]) / scale_bb;

    return result;
}

/* Function: iter_a */

/* Solve the absorption coefficient give the backscatter */
double iter_a(double rho_w, double aw, double bbw, double F, double bb) {
    double a;
    a = F * (bb + bbw) / rho_w - aw;
    return a;
}

/* Function: iter_a2 */
/* Solve the absorption coefficient given the backscatter */

/* Case II waters */
double iter_a2(double rho_w, double aw, double bbw, double F, double bb) {
    double a;
    a = F * (bb + bbw) / rho_w - (aw + bb + bbw);
    return a;
}

/* Function: biogeochem_iter */
/* Iterates until convergence on the "measured" value of total a */
/* The nomenclature is as follows: */
/* l: lower wavelength (412 or 443) */

/* u: upper wavelength (443 or 490) */
int biogeochem_iter(float al, float au, float *adyu, float *adyl, float *aphl) {
    float adyu_upper, adyu_lower, adyu_next;
    float al_m, df;
    float aphl_m, adyl_m, last_iter;
    float TOL = 0.001;
    float MAXIT = 20;
    int iter = 0, exit_iter_flag = 0;

    /* Lower first guess for adyu = 0.01*au */
    /* Upper first guess for adyu = au - biogeochem tolerance */
    adyu_lower = 0.01 * au;
    adyu_upper = au - TOL;

    /* first pass of the functions to get initial guesses */
    biogeochem_mod(au, adyu_lower, &al_m, &adyl_m, &aphl_m);
    biogeochem_mod(au, adyu_upper, &al_m, &adyl_m, &aphl_m);
    adyu_next = adyu_upper;
    do {

        /* bisect the interval */
        last_iter = adyu_next;
        adyu_next = 0.5 * (adyu_lower + adyu_upper);
        if (fabs(last_iter - adyu_next) < 0.001) {
            exit_iter_flag = 1;
            break;
        }
        biogeochem_mod(au, adyu_next, &al_m, &adyl_m, &aphl_m);
        /* test to see if overestimate or underestimate */
        if ((al_m - al) > 0.)
            adyu_upper = adyu_next;
        if ((al_m - al) < 0.)
            adyu_lower = adyu_next;

        /* absolute difference to see if convergence met */
        df = fabs(al_m - al);
        iter++;

    } while ((iter <= MAXIT) && (df >= TOL));

    if (exit_iter_flag) {
        /* Attribute all of the absorption due to gelbstoff */
        if ((al_m - al) < 0.) {
            adyu_next = au;
            adyl_m = al;
            aphl_m = TOL;
        }
        /* Attribute all absorption to phytoplankton */
        if ((al_m - al) > 0.) {
            adyu_next = TOL;
            adyl_m = TOL;
            aphl_m = al;
        }
    }
    /* Return error if maximum number of iterations exceeded */
    if (iter > maxit)
        return 1;

    *adyu = adyu_next;
    *adyl = adyl_m;
    *aphl = aphl_m;

    return 0;
}

/* Function: biogeochem_mod */
/* This function was created because of the failure of the */
/* analytical expression, which contained in effect a first order */
/* polynomial for the spectral slope */
/* Function returns the value of a(412 or 443) */
/* which can be compared with the value of a(412 or 443) */

/* obtained from the primary IOP part of the PML IOP model */
int biogeochem_mod(float au, float adyu, float *al_m, float *adyl_m, float *aphl_m) {
    float ady_part, aph_part;

    /* Numbers from the NOMAD regressions */
    /* 412:443 */
    float A = 0.059;
    float B = 1.099;
    float C = 0.229;
    float D = 0.004;
    float E = 1.033;
    float F = -0.059;

    /* 443:490 (to overcome problems with 412 channel) */
    if (ysbpa_l == 443) {
        A = 0.081;
        B = 1.186;
        C = 0.370;
        D = 0.003;
        E = 1.001;
        F = 0.183;
    }

    /* Relationships are quadratic in log10 space */
    ady_part = A * pow(log10(adyu), 2) + B * log10(adyu) + C;
    aph_part = D * pow(log10(au - adyu), 2) + E * log10(au - adyu) + F;

    *al_m = pow(10, ady_part) + pow(10, aph_part);
    *adyl_m = pow(10, ady_part);
    *aphl_m = pow(10, aph_part);

    return 0;
}
