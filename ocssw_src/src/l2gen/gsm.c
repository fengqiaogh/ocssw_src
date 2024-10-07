/* =================================================================== */
/* module gsm.c - Garver-Seigel-Maritorena 2001 Bio-optical model      */
/*                                                                     */
/* This module contains the functions to optimize and evaluate the     */
/* Garver-Seigel-Maritorena 2001 Bio-optical model.  The optimization  */
/* is performed for each set of Rrs values within a scanline, using    */
/* the Amoeba simplex minimization technique (Numerical Recipes).      */
/*                                                                     */
/* Reference:                                                          */
/* Maritorena S., D.A. Siegel & A. Peterson, Optimization of a Semi-   */
/* Analytical Ocean Color Model for Global Scale Applications, Applied */
/* Optics,  41(15): 2705-2714, 2002.                                   */
/*                                                                     */
/* Implementation:                                                     */
/* B. Franz, NASA/OCRPG/SAIC, September 2004 (rewrite of previous)     */
/*                                                                     */
/* =================================================================== */

#include <stdio.h>
#include <math.h>
#include "l12_proto.h"
#include "amoeba.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#define NPAR   3
#define MAXITR 500

#define GSMDEFAULT 0
#define CHESAPEAKE   1

#define AMOEBA  0
#define LEVMARQ 1

static float badval = BAD_FLT;

static int32_t nbandsvis;
static float grd1 = 0.0949;
static float grd2 = 0.0794;
static float *aphstar;
static float *adgstar;
static float *bbpstar;
static float *lambda;
static float *aw;
static float *bbw;

static int32_t coefset = GSMDEFAULT;
static int32_t pixday;
static float32 pixlon;
static float32 pixlat;

static double *dFdx;
static double *dxda;
static double *dxdb;

static int32_t GSMRecNum = -1;
static double *chl;
static double *adg;
static double *bbp;
static int16 *iter;

static float merit;

struct datastruct {
    size_t n;
    double * y;
    double * sigma;
} data;

int gsm_ran(int recnum) {
    if (recnum == GSMRecNum)
        return 1;
    else
        return 0;
}


/* ------------------------------------------------------------------- */
/* get_cb_model - retrieves model coefficients for Chesapeake region   */

/* ------------------------------------------------------------------- */
void set_cb_model(float wave[], int32_t nwave, float lon, float lat, int32_t day) {
    static int firstCall = 1;
    static int *iwtab;
    static int *interp;
    static int iseason = -1;

    static int32_t ntwave = 6;
    static float twave[6] = {412.0, 443.0, 490.0, 510.0, 555.0, 670.0};
    static char season[4][7] = {"Spring", "Summer", "Fall", "Winter"};

    static float tcoefs[5][4][6 + 2] = {
        { /* Upper Bay */
            {0.02119, 0.02509, 0.01282, 0.00910, 0.00427, 0.02087, 0.01218, 0.0}, /* Spring */
            {0.02653, 0.02979, 0.01655, 0.01208, 0.00470, 0.02122, 0.01218, 0.0}, /* Summer */
            {0.02259, 0.02573, 0.01372, 0.01019, 0.00376, 0.02066, 0.01218, 0.0}, /* Fall   */
            {0.02259, 0.02573, 0.01372, 0.01019, 0.00376, 0.02066, 0.01218, 0.0}
        }, /* Winter */
        { /* Mid-Bay */
            {0.02001, 0.02212, 0.01279, 0.00974, 0.00449, 0.01588, 0.01385, 0.0}, /* Spring */
            {0.03345, 0.03900, 0.02318, 0.01664, 0.00609, 0.02285, 0.01385, 0.0}, /* Summer */
            {0.02758, 0.03080, 0.01826, 0.01438, 0.00691, 0.02112, 0.01385, 0.0}, /* Fall   */
            {0.02758, 0.03080, 0.01826, 0.01438, 0.00691, 0.02112, 0.01385, 0.0}
        }, /* Winter */
        { /* Lower Bay */
            {0.02001, 0.02212, 0.01279, 0.00974, 0.00449, 0.01588, 0.01330, 0.0}, /* Spring */
            {0.03345, 0.03900, 0.02318, 0.01664, 0.00609, 0.02285, 0.01330, 0.0}, /* Summer */
            {0.02758, 0.03080, 0.01826, 0.01438, 0.00691, 0.02112, 0.01330, 0.0}, /* Fall   */
            {0.02758, 0.03080, 0.01826, 0.01438, 0.00691, 0.02112, 0.01330, 0.0}
        }, /* Winter */
        { /* Inshore MAB */
            {0.07123, 0.08843, 0.06024, 0.04072, 0.01693, 0.03815, 0.01236, 0.0}, /* Spring */
            {0.07123, 0.08843, 0.06024, 0.04072, 0.01693, 0.03815, 0.01236, 0.0}, /* Summer */
            {0.07123, 0.08843, 0.06024, 0.04072, 0.01693, 0.03815, 0.01236, 0.0}, /* Fall   */
            {0.07123, 0.08843, 0.06024, 0.04072, 0.01693, 0.03815, 0.01236, 0.0}
        }, /* Winter */
        { /* Offshore MAB */
            {0.11331, 0.14678, 0.09832, 0.06048, 0.01920, 0.04349, 0.01646, 1.0}, /* Spring */
            {0.11331, 0.14678, 0.09832, 0.06048, 0.01920, 0.04349, 0.01646, 1.0}, /* Summer */
            {0.11331, 0.14678, 0.09832, 0.06048, 0.01920, 0.04349, 0.01646, 1.0}, /* Fall   */
            {0.11331, 0.14678, 0.09832, 0.06048, 0.01920, 0.04349, 0.01646, 1.0}
        } /* Winter */

        /* aph*412,aph*443,aph*490,aph*510,aph*555,aph*670,    Sdg,   Sbbp */
    };

    float *taphstar, adg_s, bbp_s;
    int iw, iregion;
    if (firstCall == 1) {
        firstCall = 0;
        if ((iwtab = (int *) calloc(nbandsvis, sizeof (int))) == NULL) {
            printf("-E- %s line %d : error allocating memory for iwtab in gsm:set_cb_model.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((interp = (int *) calloc(nbandsvis, sizeof (int))) == NULL) {
            printf("-E- %s line %d : error allocating memory for interp in gsm:set_cb_model.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    }

    if (iseason < 0) {

        /* first call, determine the season */
        if (day >= 79 && day < 172) iseason = 0; /* Spring */
        else if (day >= 172 && day < 265) iseason = 1; /* Summer */
        else if (day >= 265 && day < 355) iseason = 2; /* Fall   */
        else iseason = 3; /* Winter */

        printf("\nInitializing Chesapeake GSM model for %s\n", season[iseason]);

        /* determine if we need to interpolate */
        for (iw = 0; iw < nwave; iw++) {
            iwtab[iw] = windex(wave[iw], twave, ntwave);
            if (fabsf(wave[iw] - twave[iwtab[iw]]) > 0.5) {
                printf("\nGSM coefficients will be interpolated for %5.1f nm band.\n", wave[iw]);
                interp[iw] = 1;
            } else
                interp[iw] = 0;
        }
    }

    /* determine region */
    iregion = -1;
    if (lat < 37.0 && lon > (-(lat + 152) / 2.5))
        iregion = 4; /* Offshore MAB */
    else if (lat >= 37.0 && lon > ((lat - 153.31) / 1.5385))
        iregion = 4; /* Offshore MAB */
    else if (lon > -75.7)
        iregion = 3; /* Inshore MAB  */
    else if (lat < 37.6)
        iregion = 2; /* Lower Bay    */
    else if (lat < 38.6)
        iregion = 1; /* Mid Bay      */
    else
        iregion = 0; /* Upper Bay    */

    /* compute coefficients */
    taphstar = &tcoefs[iregion][iseason][0];
    adg_s = tcoefs[iregion][iseason][ntwave + 0];
    bbp_s = tcoefs[iregion][iseason][ntwave + 1];
    for (iw = 0; iw < nwave; iw++) {
        if (interp[iw])
            aphstar[iw] = linterp(twave, taphstar, ntwave, wave[iw]);
        else
            aphstar[iw] = taphstar[iwtab[iw]];
        adgstar[iw] = exp(-adg_s * (wave[iw] - 443.0));
        bbpstar[iw] = pow((443.0 / wave[iw]), bbp_s);
    }

    return;
}


/* ------------------------------------------------------------------- */
/* gsm_amb() GSM model, set-up to work with amoeba optimization        */
/*                                                                     */
/* Input:                                                              */
/* auxdata - amoeba interface structure                                */
/* par[]   - model parameters                                          */
/*                                                                     */
/* Output:                                                             */
/* function returns chi-squared, which is the quantity minimized.      */
/*                                                                     */

/* ------------------------------------------------------------------- */
double gsm_amb(FITSTRUCT *auxdata, double par[]) {

    int iw;
    float x, F, bb, ac;
    double merit = 0.0;
    int nwave = auxdata->npnts;

    if (coefset == CHESAPEAKE) {
        set_cb_model(lambda, nwave, pixlon, pixlat, pixday);
    }

    for (iw = 0; iw < nwave; iw++) {

        ac = aw [iw] + aphstar[iw] * par[0] + par[1] * adgstar[iw];
        bb = bbw[iw] + par[2] * bbpstar[iw];
        x = bb / (ac + bb);
        F = grd1 * x + grd2 * pow(x, 2);

        auxdata->yfit[iw] = F;

        merit += pow((auxdata->y[iw] - F), 2) * auxdata->wgt[iw];
    }

    auxdata->merit = merit;

    return (merit);
}


/* ------------------------------------------------------------------- */
/* fit_gsm_amb() - run amoeba optimization of GSM01 for one spectrum   */

/* ------------------------------------------------------------------- */
int fit_gsm_amb(double Rrs[], double wts[], int32_t npts, double fitparms[],
        double Rrs_fit[], int16 *itercnt) {
    static float tol = 1.e-6; /* fractional change in chisqr   */
    static FITSTRUCT auxdata; /* amoeba interface structure    */
    static double init[NPAR * (NPAR + 1)]; /* initial simplex               */
    static double p0[] = {.1, .02, .0001}; /* starting parameters (C,adg,bb)*/
    static double d0[] = {5., 1.0, 0.01}; /* characteristic length         */

    int i, j;
    short isml;
    int status = 1;

    auxdata.niter = MAXITR; /* max number of iterations allowed   */
    auxdata.nfunc = NPAR; /* number of model parameters         */
    auxdata.npnts = npts; /* number of wavelengths (Rrs values) */
    auxdata.y = Rrs; /* Input Rrs values                   */
    auxdata.wgt = wts; /* Input weights on Rrs values        */
    auxdata.yfit = Rrs_fit; /* Output model predicted Rrs values  */

    /* initialize simplex with first guess model parameters */
    for (j = 0; j < NPAR + 1; j++) for (i = 0; i < NPAR; i++) init[j * NPAR + i] = p0[i];
    for (i = 0; i < NPAR; i++) {
        init[NPAR + i * (NPAR + 1)] += d0[i];
        fitparms[i] = 0.0;
    }

    /* run optimization */
    isml = amoeba(init, &auxdata, gsm_amb, tol);

    /* check convergence and record parameter results */
    if (auxdata.niter > MAXITR)
        /* failed to converge */
        status = 0;

    for (i = 0; i < NPAR; i++)
        fitparms[i] = init[NPAR * isml + i];

    *itercnt = auxdata.niter;

    return (status);
}


/* ------------------------------------------------------------------- */
/* run_gsm_amb() - returns requested GSM product for full scanline     */

/* ------------------------------------------------------------------- */
void run_gsm_amb(l2str *l2rec) {

    int32_t ip, iw;
    double model [NPAR];
    double *Rrs;
    double *wts;
    double *Rrs_fit;
    int16 itercnt;
    int16 bndcnt;
    int status;

    l1str *l1rec = l2rec->l1rec;

    if ((Rrs = (double *) calloc(l1rec->l1file->nbands, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for Rrs in gsm:run_gsm_amb.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((wts = (double *) calloc(l1rec->l1file->nbands, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for wts in gsm:run_gsm_amb.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((Rrs_fit = (double *) calloc(l1rec->l1file->nbands, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for Rrs_fit in gsm:run_gsm_amb.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    /* if this is a new scanline, fit model for every pixel of the scanline */

    for (ip = 0; ip < l1rec->npix; ip++) {

        status = 1;
        bndcnt = 0;
        chl [ip] = badval;
        adg [ip] = badval;
        bbp [ip] = badval;
        iter[ip] = 0;

        if (l1rec->mask[ip] == 0) {

            /* for CB model selection */
            pixlon = l1rec->lon[ip];
            pixlat = l1rec->lat[ip];
            int16_t year, day;
            double sec;
            unix2yds(l1rec->scantime, &year, &day, &sec);
            pixday = (int32_t) day;

            for (iw = 0; iw < nbandsvis; iw++) {
                wts[iw] = 1.0;
                Rrs[iw] = l2rec->Rrs[ip * l1rec->l1file->nbands + iw];
                Rrs[iw] = Rrs[iw] / (0.52 + 1.7 * Rrs[iw]);
                if (Rrs[iw] <= 0.0) status = 0;
                else bndcnt++;
            }

            /* if less than 3 valid radiances, fail pixel    */
            if (bndcnt < 3) {
                l1rec->flags[ip] |= PRODFAIL;
                continue;
            }

            /* if any negative reflectance, set warning flag */
            /* but allow processing to continue              */
            if (status == 0) {
                l1rec->flags[ip] |= PRODWARN;
            }

            /* fit model for this pixel */
            status = fit_gsm_amb(Rrs, wts, nbandsvis, model, Rrs_fit, &itercnt);

            /* save results and flag failure */
            if (status != 1 || model[0] <= 0.0) {
                l1rec->flags[ip] |= PRODFAIL;
            } else {
                chl[ip] = model[0];
                adg[ip] = model[1];
                bbp[ip] = model[2];
            }
            iter[ip] = itercnt;
        } else {
            l1rec->flags[ip] |= PRODFAIL;
        }
    }

    GSMRecNum = l1rec->iscan;

    free(Rrs);
    free(wts);
    free(Rrs_fit);
    return;
}


/* ------------------------------------------------------------------- */
/* gsm_lm_f() GSM01 model, set-up to work with L.M. optimization       */
/*                                                                     */
/* Input:                                                              */
/* par[]   - model parameters                                          */
/* data    - data to fit                                               */
/* f       - gsl_multifit_function                                     */
/*                                                                     */

/* ------------------------------------------------------------------- */
int gsm_lm_f(const gsl_vector *par, void *data, gsl_vector *f) {
    static int first = 1;

    int iw;
    float x, F, bb, ac;
    size_t nwave = ((struct datastruct *) data)->n;
    double *y = ((struct datastruct *) data)->y;
    double *sigma = ((struct datastruct *) data)->sigma;

    if (first == 1) {

        if ((dFdx = (double *) calloc(nwave, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for dFdx in gsm:gsm_lm_f.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((dxda = (double *) calloc(nwave, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for dxda in gsm:gsm_lm_f.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((dxdb = (double *) calloc(nwave, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for dxdb in gsm:gsm_lm_f.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        first = 0;
    }
    if (coefset == CHESAPEAKE) {
        set_cb_model(lambda, nwave, pixlon, pixlat, pixday);
    }

    merit = 0;
    for (iw = 0; iw < nwave; iw++) {
        ac = aw [iw] +
                aphstar[iw] * gsl_vector_get(par, 0) + gsl_vector_get(par, 1) * adgstar[iw];
        bb = bbw[iw] + gsl_vector_get(par, 2) * bbpstar[iw];
        x = bb / (ac + bb);
        F = grd1 * x + grd2 * pow(x, 2);

        dFdx[iw] = grd1 + 2 * grd2*x;
        dxda[iw] = -x * x / bb;
        dxdb[iw] = x / bb + dxda[iw];

        gsl_vector_set(f, iw, (F - y[iw]) / sigma[iw]);

        merit += pow((F - y[iw]) / sigma[iw], 2);
    }

    return GSL_SUCCESS;
}


/* ------------------------------------------------------------------- */
/* gsm_lm_df() GSM01 model, set-up to work with L.M. optimization      */
/*                                                                     */
/* Input:                                                              */
/* par[]   - model parameters                                          */
/* data    - data to fit                                               */
/* J       - Jacobian                                                  */
/*                                                                     */

/* ------------------------------------------------------------------- */
int gsm_lm_df(const gsl_vector *par, void *data, gsl_matrix *J) {
    int iw;
    size_t nwave = ((struct datastruct *) data)->n;
    double *sigma = ((struct datastruct *) data)->sigma;

    if (coefset == CHESAPEAKE) {
        set_cb_model(lambda, nwave, pixlon, pixlat, pixday);
    }

    for (iw = 0; iw < nwave; iw++) {

        /* Jacobian matrix J(iw,j) = dfi / dxj, */
        gsl_matrix_set(J, iw, 0, dFdx[iw] * dxda[iw] * aphstar[iw] / sigma[iw]);
        gsl_matrix_set(J, iw, 1, dFdx[iw] * dxda[iw] * adgstar[iw] / sigma[iw]);
        gsl_matrix_set(J, iw, 2, dFdx[iw] * dxdb[iw] * bbpstar[iw] / sigma[iw]);
    }

    return GSL_SUCCESS;
}

int gsm_lm_fdf(const gsl_vector *par, void *data,
        gsl_vector *f, gsl_matrix *J) {
    gsm_lm_f(par, data, f);
    gsm_lm_df(par, data, J);

    return GSL_SUCCESS;
}


/* ---------------------------------------------------------------------- */
/* run_gsm_lm() - returns requested GSM product for full scanline         */

/* ---------------------------------------------------------------------- */
void run_gsm_lm(l2str *l2rec) {
    static int firstCall = TRUE;

    int32_t ip, iw;
    double *Rrs;
    double *wts;
    double *sigma;
    int16 itercnt;
    int16 bndcnt;
    int status;
    l1str *l1rec = l2rec->l1rec;

    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;

    static gsl_multifit_function_fdf f;
    double x_init[3] = {0.0, 0.0, 0.0};
    gsl_vector_view x = gsl_vector_view_array(x_init, NPAR);
    if ((Rrs = (double *) calloc(l1rec->l1file->nbands, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for Rrs in gsm:run_gsm_amb.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((wts = (double *) calloc(l1rec->l1file->nbands, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for wts in gsm:run_gsm_amb.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((sigma = (double *) calloc(l1rec->l1file->nbands, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for sigma in gsm:run_gsm_amb.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    if (firstCall == TRUE) {

        firstCall = FALSE;

        /* Set up data structure */
        data.n = nbandsvis;
        data.y = Rrs;
        data.sigma = sigma;

        /* Set up multifit function structure */
        f.f = &gsm_lm_f;
        f.df = &gsm_lm_df;
        f.fdf = &gsm_lm_fdf;
        f.n = nbandsvis;
        f.p = NPAR;
        f.params = &data;
    }

    T = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc(T, nbandsvis, NPAR);

    /* if this is a new scanline, fit model for every pixel of the scanline */

    for (ip = 0; ip < l1rec->npix; ip++) {

        status = 1;
        bndcnt = 0;
        chl [ip] = badval;
        adg [ip] = badval;
        bbp [ip] = badval;
        iter[ip] = 0;

        if (l1rec->mask[ip] == 0) {

            /* for CB model selection */
            pixlon = l1rec->lon[ip];
            pixlat = l1rec->lat[ip];
            int16_t year, day;
            double sec;
            unix2yds(l1rec->scantime, &year, &day, &sec);
            pixday = (int32_t) day;


            for (iw = 0; iw < nbandsvis; iw++) {
                wts[iw] = 1.0;
                sigma[iw] = 1 / sqrt(wts[iw]);
                Rrs[iw] = l2rec->Rrs[ip * l1rec->l1file->nbands + iw];
                Rrs[iw] = Rrs[iw] / (0.52 + 1.7 * Rrs[iw]);
                if (Rrs[iw] <= 0.0) status = 0;
                else bndcnt++;
            }

            /* if less than 3 valid radiances, fail pixel    */
            if (bndcnt < 3) {
                l1rec->flags[ip] |= PRODFAIL;
                continue;
            }

            /* if any negative reflectance, set warning flag */
            /* but allow processing to continue              */
            if (status == 0) {
                l1rec->flags[ip] |= PRODWARN;
            }

            /* fit model for this pixel */
            gsl_multifit_fdfsolver_set(s, &f, &x.vector);
            itercnt = 0;
            do {
                itercnt++;
                status = gsl_multifit_fdfsolver_iterate(s);

                //printf ("status = %s\n", gsl_strerror (status));

                if (status)
                    break;

                status = gsl_multifit_test_delta(s->dx, s->x, 1e-4, 1e-4);

            } while (status == GSL_CONTINUE && itercnt < MAXITR);


            /* save results and flag failure */
            if (itercnt == MAXITR || gsl_vector_get(s->x, 0) <= 0.0) {
                l1rec->flags[ip] |= PRODFAIL;
            } else {
                chl[ip] = gsl_vector_get(s->x, 0);
                adg[ip] = gsl_vector_get(s->x, 1);
                bbp[ip] = gsl_vector_get(s->x, 2);
            }
            iter[ip] = itercnt;
        } else {
            l1rec->flags[ip] |= PRODFAIL;
        }
    }

    GSMRecNum = l1rec->iscan;
    gsl_multifit_fdfsolver_free(s);

    free(Rrs);
    free(wts);
    free(sigma);
    return;
}


/* ---------------------------------------------------------------------- */
/* run_gsm - runs optimization using requested method                     */

/* ---------------------------------------------------------------------- */
void run_gsm(l2str *l2rec) {
    static int firstCall = 1;
    static int last_npix;

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;
    int32_t nbands = l1file->nbands;

    /* return requested product, return -1 for abs/scatter in NIR */

    if (firstCall) {
        firstCall = 0;

        int32_t iw;
        int32_t taphn = -1;
        float *taphw = NULL;
        float *taphs = NULL;
        float adg_s = -1;
        float bbp_s = -1;

        int16_t year, day;
        double sec;
        unix2yds(l2rec->l1rec->scantime, &year, &day, &sec);
        pixday = day;

        if ((aphstar = (float *) calloc(nbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for gsm.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((adgstar = (float *) calloc(nbands, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for gsm.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((bbpstar = (float *) calloc(nbands, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for gsm.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((lambda = (float *) calloc(nbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for gsm.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((aw = (float *) calloc(nbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for gsm.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((bbw = (float *) calloc(nbands, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for gsm.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((chl = (double *) calloc(l1rec->npix, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for gsm.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((adg = (double *) calloc(l1rec->npix, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for gsm.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((bbp = (double *) calloc(l1rec->npix, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for gsm.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((iter = (int16 *) calloc(l1rec->npix, sizeof (int16))) == NULL) {
            printf("-E- %s line %d : error allocating memory for gsm.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        last_npix = l1rec->npix;

        /* set model coefficient set */
        coefset = input->gsm_opt;

        /* build array of sensor wavelengths */
        for (iw = 0; iw < l1rec->l1file->nbands; iw++) {
            lambda[iw] = l1rec->l1file->fwave[iw];
        }

        /* limit to visible wavelengths */
        nbandsvis = rdsensorinfo(l1rec->l1file->sensorID, l1_input->evalmask, "NbandsVIS", NULL);

        /* set static coefficients */
        if (coefset == GSMDEFAULT) {
            adg_s = input->gsm_adg_s;
            bbp_s = input->gsm_bbp_s;
            taphw = input->gsm_aphw; /* aphstar table wavelengths  */
            taphs = input->gsm_aphs; /* aphstar table coefficients */
            taphn = 0; /* aphstar table size         */
            for (taphn = 0; taphn < nbandsvis; taphn++)
                if (taphw[taphn] < 0) break;
            for (iw = 0; iw < nbandsvis; iw++) {
                aphstar[iw] = linterp(taphw, taphs, taphn, lambda[iw]);
                adgstar[iw] = exp(-adg_s * (lambda[iw] - 443.0));
                bbpstar[iw] = pow((443.0 / lambda[iw]), bbp_s);
            }
        }

        /* save associated irradiances and water absorption/backscatter     */
        /* use band-averaged vales if no nLw out-of-band correction applied */

        if (l1_input->outband_opt >= 2) {
            for (iw = 0; iw < nbandsvis; iw++) {
                aw [iw] = aw_spectra(lambda[iw], BANDW);
                bbw[iw] = bbw_spectra(lambda[iw], BANDW);
            }
        } else {
            for (iw = 0; iw < nbandsvis; iw++) {
                aw [iw] = l1rec->l1file->aw[iw];
                bbw[iw] = l1rec->l1file->bbw[iw];
            }
        }
    }

    // see if the number of pixels got bigger
    if (l1rec->npix > last_npix) {
        free(chl);
        if ((chl = (double *) calloc(l1rec->npix, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for gsm.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        free(adg);
        if ((adg = (double *) calloc(l1rec->npix, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for gsm.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        free(bbp);
        if ((bbp = (double *) calloc(l1rec->npix, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for gsm.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        free(iter);
        if ((iter = (int16 *) calloc(l1rec->npix, sizeof (int16))) == NULL) {
            printf("-E- %s line %d : error allocating memory for gsm.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        last_npix = l1rec->npix;
    }

    switch (input->gsm_fit) {
    case AMOEBA:
        run_gsm_amb(l2rec);
        break;
    case LEVMARQ:
        run_gsm_lm(l2rec);
        break;
    default:
        printf("%s Line %d: Unknown optimization method for GSM %d\n",
                __FILE__, __LINE__, input->gsm_fit);
        exit(1);
        break;
    }
}


/* ------------------------------------------------------------------- */
/* get_gsm() - returns requested GSM product for full scanline         */

/* ------------------------------------------------------------------- */
void get_gsm(l2str *l2rec, l2prodstr *p, float prod[]) {
    int prodID = p->cat_ix;
    int band = p->prod_ix;
    int32_t ip, iw;
 
    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;
 
    if (band >= 0) {
        if (l1file->fwave[band] > 700.0) {
            for (ip = 0; ip < l1rec->npix; ip++)
                prod[ip] = -1.0;
            return;
        }
    }

    iw = p->prod_ix;

    if (!gsm_ran(l1rec->iscan))
        run_gsm(l2rec);

    for (ip = 0; ip < l1rec->npix; ip++) {

        switch (prodID) {

        case CAT_chl_gsm:
            prod[ip] = (float) chl[ip];
            break;

        case CAT_adg_gsm:
            if (band < 0)
                prod[ip] = (float) adg[ip];
            else
                prod[ip] = (float) adg[ip] * adgstar[iw];
            break;

        case CAT_bbp_gsm:
            if (band < 0)
                prod[ip] = (float) bbp[ip];
            else
                prod[ip] = (float) bbp[ip] * bbpstar[iw];
            break;

        case CAT_aph_gsm:
            prod[ip] = (float) aphstar[iw] * chl[ip];
            break;

        case CAT_a_gsm:
            prod[ip] = (float) aw[iw] + aphstar[iw] * chl[ip] + adg[ip] * adgstar[iw];
            break;

        case CAT_bb_gsm:
            prod[ip] = (float) bbw[iw] + bbp[ip] * bbpstar[iw];
            break;

        default:
            printf("-E- %s line %d : erroneous product ID %d passed to gsm.\n",
                    __FILE__, __LINE__, prodID);
            exit(1);
        }
    }

    return;
}


/* ------------------------------------------------------------------- */
/* get_gsm_iter() - returns iteration count                          */

/* ------------------------------------------------------------------- */
int16 *get_iter_gsm(l2str *l2rec) {
    coefset = input->gsm_opt;
    if (!gsm_ran(l2rec->l1rec->iscan))
        run_gsm(l2rec);

    return iter;
}

/* Interface to convl12() to return GSM iops */

void iops_gsm(l2str *l2rec) {
    int32_t ib, ip, ipb;
    int32_t nbands = l2rec->l1rec->l1file->nbands;

    coefset = input->gsm_opt;
    if (!gsm_ran(l2rec->l1rec->iscan))
        run_gsm(l2rec);

    for (ip = 0; ip < l2rec->l1rec->npix; ip++)
        for (ib = 0; ib < nbands; ib++) {
            ipb = ip * nbands + ib;
            l2rec->a [ipb] = (float) aw[ib] + aphstar[ib] * chl[ip] + adg[ip] * adgstar[ib];
            l2rec->bb[ipb] = (float) bbw[ib] + bbp[ip] * bbpstar[ib];
        }

    return;
}

