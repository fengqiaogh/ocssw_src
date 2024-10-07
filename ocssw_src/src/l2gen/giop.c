/* =================================================================== */
/* module giop.c: generic iop model                                    */
/*                                                                     */
/* (May 2003)                                                          */
/* The module contains the functions to optimize and evaluate a user   */
/* controlled IOP model, which currently defaults to the GSM 2002.     */
/*                                                                     */
/* (May 2015)                                                          */
/* Raman scattering correction now applied to Rrs                      */
/*                                                                     */
/* (Dec 2015)                                                          */
/* Includes iterative adaptive SIOP matrix inversion that iterates     */
/* through a number of IOP spectral shapes.                            */
/*                                                                     */
/* (Oct 2023)                                                          */
/*  Output parameter uncertainties computed per McKinna et al (2019)   */
/*  bbp line height model per McKinna et al (2021)                     */
/*  bbp chl-based model per Huot et al (2008)                          */
/*                                                                     */
/* Implementation:                                                     */
/* B. Franz, NASA/OBPG/SAIC, May 2008                                  */
/*                                                                     */
/* Reference:                                                          */
/* Werdell et al. (2013) Generalized ocean color inversion model for   */
/* retrieving marine inherent optical properties, Applied Optics, 52,  */
/* 2019-2037.                                                          */
/*                                                                     */
/* Notes:                                                              */
/* This code was written with the intent that it may become a base for */
/* multiple algorithms (to be writtenn as wrappers).  As such, it can't*/
/* be assumed that the starting config is the only config (e.g.,       */
/* number of wavelengths to fit), which leads to some inefficiencies.  */
/* =================================================================== */

#include <stdio.h>
#include <math.h>
#include "l12_proto.h"
#include "giop.h"
#include "amoeba.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>


typedef double realtype;

typedef struct fit_data_str {
    double *y;
    double *w;
    giopstr *g;
} fitstr;

static int32_t LastRecNum = -1;
static float badval = BAD_FLT;
static float bbp_s_max = 5.0;
static float bbp_s_min = 0.0;
static float adg_s_max = 0.025;
static float adg_s_min = 0.01;
static float chl_max = 10.0;
static float chl_min = 0.03;

static float *aw;
static float *bbw;
static float *foq;

static float *aph1;
static float *adg1;
static float *bbp1;

static float *anap1;
static float *acdom1;
static float *bbph1;
static float *bbnap1;


// these need to be passed from run_giop if giop is to
// support simultaneous products (e.g., a gsm wrapper)

static int16 *iter;
static int16 *iopf;
static float *a;
static float *bb;
static float *aph;
static float *adg;
static float *bbp;

static float *a_unc;
static float *bb_unc;

static float *anap;
static float *acdom;
static float *bbph;
static float *bbnap;

static float *chl;
static float *aph_s;
static float *adg_s;
static float *uadg_s;
static float *bbp_s;
static float *ubbp_s;
static float *rrsdiff;
static float *mRrs;
static float **fit_par;
static float *chisqr;
static int16 max_npar;

static int *siop_num;
static int allocateRrsRaman = 0;

void freeArray(void **a, int32_t m) {
    int i;
    for (i = 0; i < m; ++i) {
        free(a[i]);
    }
    free(a);
}

void freeDArray(double **a, int32_t m) {
    int i;
    for (i = 0; i < m; ++i) {
        free(a[i]);
    }
    free(a);
}

/* ------------------------------------------------------------------- */
/* giop_int_tab_file - reads a table of eigenvectors and interpolates  */
/* to input wavelengths and returns a 2-D array of vectors.            */
/* Function reads in asscoicated uncertainty vectors from              */
/* companion file (if avaiable)                                        */

/* ------------------------------------------------------------------- */
void giop_int_tab_file(char *file, char *ufile, int nx, float *x, float **y, float **uy) {

    float *table;
    float *xtab;
    float *ytab;
    float *utable;
    float *uxtab;
    float *uytab;
    int ncol, nrow, uncol, unrow;
    int i, ivec;

    printf("\nLoading basis vector from %s\n", file);


    nrow = table_row_count(file);
    if (nrow <= 0) {
        printf("-E- %s line %d: error opening (%s) file", __FILE__, __LINE__, file);
        exit(1);
    }

    ncol = table_column_count(file);
    table = table_read_r4(file, ncol, nrow);

    // Interpolate to input wavelengths (x)

    xtab = &table[nrow * 0];

    for (ivec = 0; ivec < ncol - 1; ivec++) {
        ytab = &table[nrow * (ivec + 1)];
        for (i = 0; i < nx; i++) {
            y[ivec][i] = linterp(xtab, ytab, nrow, x[i]);
            uy[ivec][i] = 0.0;
        }
    }

    table_free_r4(table);

    //Check if uncertainty file provided
    if (!ufile) {
        printf("\nNo uncertainty basis vector file ascociated with %s\n", file);
    } else {

        unrow = table_row_count(ufile);
        if (unrow <= 0) {
            printf("-E- %s line %d: error opening (%s) file", __FILE__, __LINE__, ufile);
            exit(1);
        }

        uncol = table_column_count(ufile);
        utable = table_read_r4(ufile, uncol, unrow);

        // Interpolate to input wavelengths (x)

        uxtab = &utable[unrow * 0];

        for (ivec = 0; ivec < uncol - 1; ivec++) {
            uytab = &utable[unrow * (ivec + 1)];
            for (i = 0; i < nx; i++) {
                uy[ivec][i] = linterp(uxtab, uytab, unrow, x[i]);
            }
        }

        table_free_r4(utable);

    }

}


/* ------------------------------------------------------------------- */
/* giop_ctl_start - set starting values for optimization               */

/* ------------------------------------------------------------------- */
void giop_ctl_start(giopstr *g, float chl) {

    int ipar, ivec;

    // set starting params based on model elements

    ipar = 0;

    for (ivec = 0; ivec < g->aph_nvec; ivec++) {
        g->par[ipar] = chl / g->aph_nvec;
        g->len[ipar] = 0.5;
        ipar++;
    }

    for (ivec = 0; ivec < g->adg_nvec; ivec++) {
        g->par[ipar] = 0.01;
        g->len[ipar] = 0.1;
        ipar++;
    }

    switch (g->bbp_opt) {
    case BBPLASFIX:
    case BBPQAAFIX:
        break;
    case BBPLAS:
        g->par[ipar] = 1.0;
        g->len[ipar] = 0.01;
        ipar++;
        break;
    default:
        for (ivec = 0; ivec < g->bbp_nvec; ivec++) {
            g->par[ipar] = 0.001;
            g->len[ipar] = 0.01;
            ipar++;
        }
        break;
    }

    //note: for SVDSIOP ipar > npar
    if (ipar != g->npar && g->fit_opt != SVDSIOP) {
        printf("-E- %s Line %d: number of GIOP fit parameters (%d) does not equal number of initialized parameter (%d)\n",
                __FILE__, __LINE__, g->npar, ipar);
        exit(1);
    }
}

/* ------------------------------------------------------------------- */
/* giop_ctl_init - initialize model control structure                  */
/*                                                                     */
/* Provides default values or user-supplied parameters for all static  */
/* model components.  These may later be over-riden by dynamic model   */
/* parameters (e.g., band-ratio based bbp_s).                          */
/*                                                                     */

/* ------------------------------------------------------------------- */
void giop_ctl_init(giopstr *g, int nwave, float wave[],
        float aw[], float bbw[]) {

    int iw, iwx;

    float def_aph_w = 443.0;
    float def_aph_s = 0.5;
    float def_adg_w = 443.0;
    float def_adg_s = 0.02061;
    float def_bbp_w = 443.0;
    float def_bbp_s = 1.03373;
    float def_grd[] = {0.0949, 0.0794};

    // determine indices of wavelengths to be used (default to all vis)

    if (input->giop_wave[0] > 0) {
        for (iw = 0; iw < nwave; iw++) {
            if (input->giop_wave[iw] <= 0) break;
            iwx = windex(input->giop_wave[iw], wave, nwave);
            g->bindx[iw] = iwx;
            //windex(input->giop_wave[iw],wave,nwave);
        }
        g->nwave = iw;
    } else {
        g->nwave = MIN(windex(671., wave, nwave) + 1, nwave);
        for (iw = 0; iw < g->nwave; iw++)
            g->bindx[iw] = iw;
    }

    // set wavelengths and pure-water values

    for (iw = 0; iw < g->nwave; iw++) {
        g->wave[iw] = wave[g->bindx[iw]];
        g->aw [iw] = aw [g->bindx[iw]];
        g->bbw [iw] = bbw [g->bindx[iw]];
    }

    //rrs uncertainties options
    
    if (input->giop_rrs_unc_opt >= 0) {
        g->urrs_opt = input->giop_rrs_unc_opt;
    } else {
        g->urrs_opt = URRSNONE;
    }

    //Input Rrs
    /*
    if (input->giop_rrs_unc[0] > 0) {
        g->wt_opt = 1;
        for (iw = 0; iw < nwave; iw++) {
            if (input->giop_rrs_unc[iw] <= 0) break;
            g->wts[iw] = 1.0 / pow(input->giop_rrs_unc[iw], 2);
        }
        if (iw != g->nwave) {
            printf("-E- %s line %d: number of giop_rrs_unc (%d) must equal number of giop_wave (%d)",
                    __FILE__, __LINE__, iw, g->nwave);
            exit(1);
        }
    } else {
        g->wt_opt = 0;
        for (iw = 0; iw < g->nwave; iw++)
            g->wts[iw] = 1.0;
    }
    */

    // maximum iterations

    if (input->giop_maxiter > 0)
        g->maxiter = input->giop_maxiter;
    else
        g->maxiter = 500;

    // fitting method

    if (input->giop_fit_opt > 0)
        g->fit_opt = input->giop_fit_opt;
    else
        g->fit_opt = LEVMARQ;

    // Rrs to bb/(a+bb) method

    if (input->giop_rrs_opt > 0)
        g->rrs_opt = input->giop_rrs_opt;
    else
        g->rrs_opt = RRSGRD;


    // set coefficients of Gordon quadratic

    if (input->giop_grd[0] > -999.0) {
        g->grd[0] = input->giop_grd[0];
        g->grd[1] = input->giop_grd[1];
    } else {
        g->grd[0] = def_grd[0];
        g->grd[1] = def_grd[1];
    }

    // default basis vectors

    strcpy(g->aph_tab_file, input->giop_aph_file);
    strcpy(g->adg_tab_file, input->giop_adg_file);
    strcpy(g->bbp_tab_file, input->giop_bbp_file);
    strcpy(g->acdom_tab_file, input->giop_acdom_file);
    strcpy(g->anap_tab_file, input->giop_anap_file);
    strcpy(g->bbph_tab_file, input->giop_bbph_file);
    strcpy(g->bbnap_tab_file, input->giop_bbnap_file);

    //set defaults
    g->aph_opt = APHTAB;
    g->adg_opt = ADGTAB;
    g->bbp_opt = BBPTAB;
    g->acdom_opt = ACDOMNONE;
    g->anap_opt = ANAPNONE;
    g->bbnap_opt = BBPHNONE;
    g->bbnap_opt = BBNAPNONE;

    // aphstar function

    if (input->giop_aph_opt > 0)
        g->aph_opt = input->giop_aph_opt;

    if (input->giop_aph_w > 0.0)
        g->aph_w = input->giop_aph_w;
    else
        g->aph_w = def_aph_w;

    if (input->giop_aph_s > -999.0)
        g->aph_s = input->giop_aph_s;
    else
        g->aph_s = def_aph_s;


    // adgstar function

    //if (input->giop_adg_opt > 0) 
    //    g->adg_opt = input->giop_adg_opt; 
    if (input->giop_adg_opt != 1)
        g->adg_opt = input->giop_adg_opt;
    else
        g->adg_opt = ADGS;

    if (input->giop_adg_w > 0.0)
        g->adg_w = input->giop_adg_w;
    else
        g->adg_w = def_adg_w;

    if (input->giop_adg_s > -999.0) {
        g->adg_s = input->giop_adg_s;
        g->uadg_s = input->giop_uadg_s;
    } else {
        g->adg_s = def_adg_s;
        g->uadg_s = 0.0;
    }

    //acdom star function
    if (input->giop_acdom_opt != 1)
        g->acdom_opt = input->giop_acdom_opt;
    else
        g->acdom_opt = ACDOMNONE;

    //anap star function
    if (input->giop_anap_opt != 1)
        g->anap_opt = input->giop_anap_opt;
    else
        g->anap_opt = ANAPNONE;

    // bbpstar function

    //if (input->giop_bbp_opt > 0) 
    //    g->bbp_opt = input->giop_bbp_opt;
    if (input->giop_bbp_opt != 1)
        g->bbp_opt = input->giop_bbp_opt;
    else
        g->bbp_opt = BBPS;

    if (input->giop_bbp_w > 0.0)
        g->bbp_w = input->giop_bbp_w;
    else
        g->bbp_w = def_bbp_w;

    if (input->giop_bbp_s > -999.0) {
        g->bbp_s = input->giop_bbp_s;
        g->ubbp_s = input->giop_ubbp_s;
    } else {
        g->bbp_s = def_bbp_s;
        g->ubbp_s = 0.0;
    }

    //bbph star function
    if (input->giop_bbph_opt != 1)
        g->bbph_opt = input->giop_bbph_opt;
    else
        g->bbph_opt = BBPHNONE;

    //bbnap star function
    if (input->giop_bbnap_opt != 1)
        g->bbnap_opt = input->giop_bbnap_opt;
    else
        g->bbnap_opt = BBNAPNONE;


    // set number of vectors

    switch (g->aph_opt) {
    case APHTAB:
        g->aph_nvec = table_column_count(g->aph_tab_file) - 1;
        break;
    default:
        g->aph_nvec = 1;
        break;
    }

    switch (g->adg_opt) {
    case ADGTAB:
        g->adg_nvec = table_column_count(g->adg_tab_file) - 1;
        break;
    default:
        g->adg_nvec = 1;
        break;
    }

    switch (g->acdom_opt) {
    case ACDOMTAB:
        g->acdom_nvec = table_column_count(g->acdom_tab_file) - 1;
        break;
    default:
        g->acdom_nvec = 1;
        break;
    }

    switch (g->anap_opt) {
    case ANAPTAB:
        g->anap_nvec = table_column_count(g->anap_tab_file) - 1;
        break;
    default:
        g->anap_nvec = 1;
        break;
    }

    switch (g->bbp_opt) {
    case BBPTAB:
        g->bbp_nvec = table_column_count(g->bbp_tab_file) - 1;
        break;
    case BBPLASFIX:
    case BBPQAAFIX:
    case BBPLH:
    case BBPCHL:
        g->bbp_nvec = 1;
        break;
    default:
        g->bbp_nvec = 1;
        break;
    }

    switch (g->bbph_opt) {
    case BBPHTAB:
        g->bbph_nvec = table_column_count(g->bbph_tab_file) - 1;
        break;
    default:
        g->bbph_nvec = 1;
        break;
    }

    switch (g->bbnap_opt) {
    case BBNAPTAB:
        g->bbnap_nvec = table_column_count(g->bbnap_tab_file) - 1;
        break;
    default:
        g->bbnap_nvec = 1;
        break;
    }

    // total number of parameters to be optimized
    if (g->fit_opt == SVDSIOP) {
        //only fitting for chl, cdom and nap
        g->npar = 3;
    } else {
        g->npar = g->aph_nvec + g->adg_nvec + g->bbp_nvec;
    }


    // allocate space for vectors (one element per sensor wavelength)

    g->aph_tab_nw = nwave;
    g->aph_tab_w = (float *) calloc(nwave, sizeof (float));
    g->aph_tab_s = (float **) allocate2d_float(g->aph_tab_nw, g->aph_nvec);
    g->uaph_tab_s = (float **) allocate2d_float(g->aph_tab_nw, g->aph_nvec);

    g->adg_tab_nw = nwave;
    g->adg_tab_w = (float *) calloc(nwave, sizeof (float));
    g->adg_tab_s = (float **) allocate2d_float(g->adg_tab_nw, g->adg_nvec);
    g->uadg_tab_s = (float **) allocate2d_float(g->adg_tab_nw, g->adg_nvec);

    g->acdom_tab_nw = nwave;
    g->acdom_tab_w = (float *) calloc(nwave, sizeof (float));
    g->acdom_tab_s = (float **) allocate2d_float(g->acdom_tab_nw, g->acdom_nvec);
    g->uacdom_tab_s = (float **) allocate2d_float(g->acdom_tab_nw, g->acdom_nvec);

    g->anap_tab_nw = nwave;
    g->anap_tab_w = (float *) calloc(nwave, sizeof (float));
    g->anap_tab_s = (float **) allocate2d_float(g->anap_tab_nw, g->anap_nvec);
    g->uanap_tab_s = (float **) allocate2d_float(g->anap_tab_nw, g->anap_nvec);

    g->bbp_tab_nw = nwave;
    g->bbp_tab_w = (float *) calloc(nwave, sizeof (float));
    g->bbp_tab_s = (float **) allocate2d_float(g->bbp_tab_nw, g->bbp_nvec);
    g->ubbp_tab_s = (float **) allocate2d_float(g->bbp_tab_nw, g->bbp_nvec);

    g->bbph_tab_nw = nwave;
    g->bbph_tab_w = (float *) calloc(nwave, sizeof (float));
    g->bbph_tab_s = (float **) allocate2d_float(g->bbph_tab_nw, g->bbph_nvec);
    g->ubbph_tab_s = (float **) allocate2d_float(g->bbph_tab_nw, g->bbph_nvec);

    g->bbnap_tab_nw = nwave;
    g->bbnap_tab_w = (float *) calloc(nwave, sizeof (float));
    g->bbnap_tab_s = (float **) allocate2d_float(g->bbnap_tab_nw, g->bbnap_nvec);
    g->ubbnap_tab_s = (float **) allocate2d_float(g->bbnap_tab_nw, g->bbnap_nvec);

    // set vector wavelengths (same as ALL sensor for now)

    for (iw = 0; iw < nwave; iw++) {
        g->aph_tab_w[iw] = wave[iw];
        g->adg_tab_w[iw] = wave[iw];
        g->acdom_tab_w[iw] = wave[iw];
        g->anap_tab_w[iw] = wave[iw];
        g->bbp_tab_w[iw] = wave[iw];
        g->bbph_tab_w[iw] = wave[iw];
        g->bbnap_tab_w[iw] = wave[iw];
    }

    // load tabulated vectors if requested

    switch (g->aph_opt) {
    case APHTAB:
        giop_int_tab_file(g->aph_tab_file, g->uaph_tab_file, nwave, wave, g->aph_tab_s, g->uaph_tab_s);
        break;
    }

    switch (g->adg_opt) {
    case ADGTAB:
        giop_int_tab_file(g->adg_tab_file, g->uadg_tab_file,nwave, wave, g->adg_tab_s, g->uadg_tab_s);
        break;
    }

    switch (g->bbp_opt) {
    case BBPTAB:
        giop_int_tab_file(g->bbp_tab_file, g->ubbp_tab_file, nwave, wave, g->bbp_tab_s, g->ubbp_tab_s);
        break;
    }
    switch (g->acdom_opt) {
    case ACDOMTAB:
        giop_int_tab_file(g->acdom_tab_file,g->uacdom_tab_file, nwave, wave, g->acdom_tab_s, g->uacdom_tab_s);
        break;
    }
    switch (g->anap_opt) {
    case ANAPTAB:
        giop_int_tab_file(g->anap_tab_file, g->uanap_tab_file, nwave, wave, g->anap_tab_s, g->uanap_tab_s);
        break;
    }
    switch (g->bbph_opt) {
    case BBPHTAB:
        giop_int_tab_file(g->bbph_tab_file, g->ubbph_tab_file, nwave, wave, g->bbph_tab_s, g->ubbph_tab_s);
        break;
    }
    switch (g->bbnap_opt) {
    case BBNAPTAB:
        giop_int_tab_file(g->bbnap_tab_file, g->ubbnap_tab_file, nwave, wave, g->bbnap_tab_s, g->ubbnap_tab_s);
        break;
    }

}


/* ------------------------------------------------------------------- */
/*                                                                     */

/* ------------------------------------------------------------------- */
int giop_ran(int recnum) {
    if (recnum == LastRecNum)
        return 1;
    else
        return 0;
}

/* ------------------------------------------------------------------- */
/*                                                                     */

/* ------------------------------------------------------------------- */
void giop_model(giopstr *g, double par[], int nwave, float wave[], float aw[], float bbw[],
        float foq[], float aph[], float adg[], float bbp[],
        double rrs[], double **dfdpar, double **parstar) {
    //    double rrs[],double dfdpar[NBANDS][GIOPMAXPAR],double parstar[NBANDS][GIOPMAXPAR])
    int iw, iwtab, ivec, ipar;
    float bb, a;
    float x;
    float *aphstar;
    float *adgstar;
    float *bbpstar;
    float *uaphstar;
    float *uadgstar;
    float *ubbpstar;
    float dfdx;
    float dxda;
    float dxdb;

    // note: input wavelength set can vary between optimization and evaluation calls
    if ((aphstar = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:giop_model.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((adgstar = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:giop_model.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((bbpstar = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:giop_model.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((uaphstar = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:giop_model.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((uadgstar = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:giop_model.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((ubbpstar = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:giop_model.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    for (iw = 0; iw < nwave; iw++) {

        // evaluate model for each IOP component at all wavelengths
        // also compute and store associated uncertainty, based on fitted paraneter
        // uncertainties as stored in the second half of the par array
        // note that we're using a simple summation for uncertainties (for now)

        ipar = 0;

        switch (g->aph_opt) {
        case APHGAUSS:
            aphstar[ipar] = exp(-1. * pow(wave[iw] - g->aph_w, 2) / (2 * pow(g->aph_s, 2)));
            uaphstar[ipar] = 0;
            aph[iw] = par[ipar] * aphstar[ipar];
            aph[iw + nwave] = sqrt(pow(par[ipar + g->npar]*aphstar[ipar] ,2) 
                    + pow(par[ipar]*uaphstar[ipar],2));
            ipar++;
            break;
        default:
            aph[iw] = 0;
            aph[iw + nwave] = 0;
            for (ivec = 0; ivec < g->aph_nvec; ivec++) {
                iwtab = windex(wave[iw], g->aph_tab_w, g->aph_tab_nw);
                aphstar[ipar] = g->aph_tab_s[ivec][iwtab];
                uaphstar[ipar] = g->uaph_tab_s[ivec][iwtab];
                aph[iw] += par[ipar] * aphstar[ipar];
                aph[iw + nwave] += sqrt(pow(par[ipar + g->npar]*aphstar[ipar],2)
                        + pow(par[ipar]*uaphstar[ipar],2));
                ipar++;
            }
            break;
        }

        switch (g->adg_opt) {
        case ADGS:
        case ADGSQAA:
        case ADGSOBPG:
            adgstar[ipar] = exp(-g->adg_s * (wave[iw] - g->adg_w));
            uadgstar[ipar] = sqrt(pow(g->uadg_s*(g->adg_w - wave[iw])*exp(-g->adg_s * (wave[iw] - g->adg_w)) ,2));
            adg[iw] = par[ipar] * adgstar[ipar];
            adg[iw + nwave] = sqrt(pow(par[ipar + g->npar] * adgstar[ipar],2) + pow(par[ipar]*uadgstar[ipar],2));
            ipar++;
            break;
        default:
            adg[iw] = 0;
            adg[iw + nwave] = 0;
            for (ivec = 0; ivec < g->adg_nvec; ivec++) {
                iwtab = windex(wave[iw], g->adg_tab_w, g->adg_tab_nw);
                adgstar[ipar] = g->adg_tab_s[ivec][iwtab];
                uadgstar[ipar] = g->uadg_tab_s[ivec][iwtab];
                adg[iw] += par[ipar] * adgstar[ipar];
                adg[iw + nwave] += sqrt(pow(par[ipar + g->npar] * adgstar[ipar],2)
                        + pow(par[ipar]*uadgstar[ipar],2));
                ipar++;
            }
            break;
        }

        switch (g->bbp_opt) {
        case BBPLASFIX:
        case BBPQAAFIX:
            iwtab = windex(wave[iw], g->bbp_tab_w, g->bbp_tab_nw);
            bbpstar[ipar] = g->bbp_tab_s[0][iwtab];
            bbp[iw] = bbpstar[ipar];
            bbp[iw + nwave] = 0.0;
            break;
        case BBPLH:
        case BBPCHL:
            iwtab = windex(wave[iw], g->bbp_tab_w, g->bbp_tab_nw);
            bbpstar[ipar] = g->bbp_tab_s[0][iwtab];
            ubbpstar[ipar] = g->ubbp_tab_s[0][iwtab];
            bbp[iw] = bbpstar[ipar];
            bbp[iw + nwave] = ubbpstar[ipar];
            break;
        case BBPS:
        case BBPSLAS:
        case BBPSHAL:
        case BBPSQAA:
        case BBPSCIOTTI:
        case BBPSMM01:
            bbpstar[ipar] = pow((g->bbp_w / wave[iw]), g->bbp_s);
            ubbpstar[ipar] = sqrt(pow(g->ubbp_s*pow(g->bbp_w / wave[iw], g->bbp_s)*log(g->bbp_w / wave[iw]),2));
            bbp[iw] = par[ipar] * bbpstar[ipar];
            bbp[iw + nwave] = sqrt(pow(par[ipar + g->npar] * bbpstar[ipar],2) + 
                    pow(par[ipar]*ubbpstar[ipar],2));
            ipar++;
            break;
        default:
            bbp[iw] = 0;
            bbp[iw + nwave] = 0;
            for (ivec = 0; ivec < g->bbp_nvec; ivec++) {
                iwtab = windex(wave[iw], g->bbp_tab_w, g->bbp_tab_nw);
                bbpstar[ipar] = g->bbp_tab_s[ivec][iwtab];
                ubbpstar[ipar] = g->ubbp_tab_s[ivec][iwtab];
                bbp[iw] += par[ipar] * bbpstar[ipar];
                bbp[iw + nwave] += sqrt(pow(par[ipar + g->npar] * bbpstar[ipar],2) +
                        pow(par[ipar]*ubbpstar[ipar],2));
                ipar++;
            }
            break;
        }

        a = aw [iw] + aph[iw] + adg[iw];
        bb = bbw[iw] + bbp[iw];
        x = bb / (a + bb);

        switch (g->rrs_opt) {
        case RRSGRD:
            rrs[iw] = g->grd[0] * x + g->grd[1] * pow(x, 2);
            dfdx = g->grd[0] + 2 * g->grd[1] * x;
            break;
        case RRSFOQ:
            rrs[iw] = foq[iw] * x;
            dfdx = foq[iw];
            break;
        }

        if (dfdpar != NULL) {

            ipar = 0;

            dxda = -x * x / bb;
            dxdb = x / bb + dxda;

            switch (g->aph_opt) {
            default:
                for (ivec = 0; ivec < g->aph_nvec; ivec++) {
                    dfdpar[iw][ipar] = dfdx * dxda * aphstar[ipar];
                    ipar++;
                }
                break;
            }

            switch (g->adg_opt) {
            default:
                for (ivec = 0; ivec < g->adg_nvec; ivec++) {
                    dfdpar[iw][ipar] = dfdx * dxda * adgstar[ipar];
                    ipar++;
                }
                break;
            }

            switch (g->bbp_opt) {
            case BBPLASFIX:
            case BBPQAAFIX:
            case BBPLH:
            case BBPCHL:
                break;
            default:
                for (ivec = 0; ivec < g->bbp_nvec; ivec++) {
                    dfdpar[iw][ipar] = dfdx * dxdb * bbpstar[ipar];
                    ipar++;
                }
                break;
            }
        }

        if (parstar != NULL) {

            ipar = 0;

            switch (g->aph_opt) {
            default:
                for (ivec = 0; ivec < g->aph_nvec; ivec++) {
                    parstar[iw][ipar] = aphstar[ipar];
                    ipar++;
                }
                break;
            }

            switch (g->adg_opt) {
            default:
                for (ivec = 0; ivec < g->adg_nvec; ivec++) {
                    parstar[iw][ipar] = adgstar[ipar];
                    ipar++;
                }
                break;
            }

            switch (g->bbp_opt) {
            case BBPLASFIX:
            case BBPQAAFIX:
            case BBPLH:
            case BBPCHL:
                break;
            default:
                for (ivec = 0; ivec < g->bbp_nvec; ivec++) {
                    parstar[iw][ipar] = bbpstar[ipar];
                    ipar++;
                }
                break;
            }
        }

    }

    free(aphstar);
    free(adgstar);
    free(bbpstar);
    free(uaphstar);
    free(uadgstar);
    free(ubbpstar);
    return;
}
/* ------------------------------------------------------------------- */
/*                                                                     */

/* ------------------------------------------------------------------- */
void giop_model_iterate(giopstr *g, double par[], int nwave, float wave[], float aw[], float bbw[],
        float foq[], float aph[], float adg[], float bbp[], float acdom[], float anap[], float bbph[], float bbnap[],
        double rrs[], double **dfdpar, double **parstar) {
    //    double rrs[],double dfdpar[NBANDS][GIOPMAXPAR],double parstar[NBANDS][GIOPMAXPAR])
    int iw, iwtab, idx443;

    int16 isiop = g->siopIdx;

    float bb, a;
    float x;
    float acdom443;

    float *aphstar;
    float *acdomstar;
    float *anapstar;
    float *bbphstar;
    float *bbnapstar;

    float dfdx;
    float dxda;
    float dxdb;


    /* note: input wavelength set can vary between optimization and evaluation calls*/
    if ((aphstar = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:giop_model.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((acdomstar = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:giop_model.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((anapstar = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:giop_model.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((bbphstar = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:giop_model.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((bbnapstar = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:giop_model.\n",
                __FILE__, __LINE__);
        exit(1);
    }


    /*Sanity check: do the input tables have the same number of columns?*/
    switch (g->fit_opt) {
    case SVDSIOP:
        if (g->aph_nvec != g->acdom_nvec || g->aph_nvec != g->anap_nvec
                || g->aph_nvec != g->bbph_nvec || g->aph_nvec != g->bbnap_nvec
                || g->acdom_nvec != g->anap_nvec || g->acdom_nvec != g->bbph_nvec
                || g->acdom_nvec != g->bbnap_nvec || g->anap_nvec != g->bbph_nvec
                || g->anap_nvec != g->bbnap_nvec || g->bbph_nvec != g->bbnap_nvec) {

            printf("-E- %s line %d: SIOP tables must have same number of columns.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        break;
    }

    idx443 = windex(443., g->acdom_tab_w, g->acdom_tab_nw);

    for (iw = 0; iw < nwave; iw++) {

        // evaluate model for each IOP component at all wavelengths
        // also compute and store associated uncertainty, based on fitted paraneter
        // uncertainties as stored in the second half of the par array
        // note that we're using a simple summation for uncertainties (for now)

        iwtab = windex(wave[iw], g->aph_tab_w, g->aph_tab_nw);
        aphstar[iw] = g->aph_tab_s[isiop][iwtab];
        aph[iw] = par[0] * aphstar[iw];
        aph[iw + nwave] = par[g->npar] * aphstar[iw];

        /*Esnure that acdom443 is equal to 1.0*/
        iwtab = windex(wave[iw], g->acdom_tab_w, g->acdom_tab_nw);
        acdomstar[iw] = g->acdom_tab_s[isiop][iwtab];
        acdom443 = g->acdom_tab_s[isiop][idx443];
        acdom[iw] = par[1]*(acdomstar[iw] / acdom443);
        acdom[iw + nwave] = par[1 + g->npar]*(aphstar[iw] / acdom443);

        iwtab = windex(wave[iw], g->anap_tab_w, g->anap_tab_nw);
        anapstar[iw] = g->anap_tab_s[isiop][iwtab];
        anap[iw] = par[2] * anapstar[iw];
        anap[iw + nwave] = par[2 + g->npar] * aphstar[iw];

        adg[iw] = par[1] * acdomstar[iw] + par[2] * anapstar[iw];
        adg[iw + nwave] = par[1 + g->npar] * acdomstar[iw] + par[2 + g->npar] * anapstar[iw];

        iwtab = windex(wave[iw], g->bbph_tab_w, g->bbph_tab_nw);
        bbphstar[iw] = g->bbph_tab_s[isiop][iwtab];
        bbph[iw] = par[0] * bbphstar[iw];
        bbph[iw + nwave] = par[g->npar] * bbphstar[iw];

        iwtab = windex(wave[iw], g->bbnap_tab_w, g->bbnap_tab_nw);
        bbnapstar[iw] = g->bbnap_tab_s[isiop][iwtab];
        bbnap[iw] = par[2] * bbnapstar[iw];
        bbnap[iw + nwave] = par[2 + g->npar] * bbnapstar[iw];

        bbp[iw] = par[0] * bbphstar[iw] + par[2] * bbnapstar[iw];
        bbp[iw + nwave] = par[g->npar] * bbphstar[iw] + par[2 + g->npar] * bbnapstar[iw];

        a = aw [iw] + aph[iw] + adg[iw];
        bb = bbw[iw] + bbp[iw];
        x = bb / (a + bb);

        switch (g->rrs_opt) {
        case RRSGRD:
            rrs[iw] = g->grd[0] * x + g->grd[1] * pow(x, 2);
            dfdx = g->grd[0] + 2 * g->grd[1] * x;
            break;
        case RRSFOQ:
            rrs[iw] = foq[iw] * x;
            dfdx = foq[iw];
            break;
        }

        if (dfdpar != NULL) {

            dxda = -x * x / bb;
            dxdb = x / bb + dxda;

            dfdpar[iw][0] = dfdx * dxda * aphstar[iw];
            dfdpar[iw][1] = dfdx * dxda * acdomstar[iw];
            dfdpar[iw][2] = dfdx * dxdb * anapstar[iw];
            dfdpar[iw][3] = dfdx * dxdb * bbphstar[iw];
            dfdpar[iw][4] = dfdx * dxdb * bbnapstar[iw];


        }

        if (parstar != NULL) {

            parstar[iw][0] = aphstar[iw];
            parstar[iw][1] = acdomstar[iw];
            parstar[iw][2] = anapstar[iw];
            parstar[iw][3] = bbphstar[iw];
            parstar[iw][4] = bbnapstar[iw];

        }

    }

    free(aphstar);
    free(acdomstar);
    free(anapstar);
    free(bbphstar);
    free(bbnapstar);

    return;
}

/* ------------------------------------------------------------------- */
/*                                                                     */

/* ------------------------------------------------------------------- */
double giop_amb(FITSTRUCT *ambdata, double par[]) {
    int iw;

    giopstr *g = (giopstr *) (ambdata->meta);

    giop_model(g, par, g->nwave, g->wave, g->aw, g->bbw, g->foq, aph1, adg1, bbp1, ambdata->yfit, NULL, NULL);

    ambdata->merit = 0.0;
    for (iw = 0; iw < g->nwave; iw++) {
        ambdata->merit += pow((ambdata->y[iw] - ambdata->yfit[iw]), 2) * ambdata->wgt[iw];
    }

    return (ambdata->merit);
}


/* ------------------------------------------------------------------- */
/*                                                                     */

/* ------------------------------------------------------------------- */
int fit_giop_amb(giopstr *g, double Rrs[], double wts[], double par[],
        double Rrs_fit[], int16 *itercnt) {
    int status = 0;


    static float tol = 1.e-6; /* fractional change in chisqr */
    static FITSTRUCT ambdata; /* amoeba interface structure  */
    static double *init; /* initial simplex    */

    static int firstCall = 1;

    int i, j;
    short isml;

    ambdata.niter = g->maxiter; /* max number of iterations          */
    ambdata.nfunc = g->npar; /* number of model parameters        */
    ambdata.npnts = g->nwave; /* number of wavelengths (Rrs)       */
    ambdata.y = Rrs; /* Input Rrs values (subsurface)     */
    ambdata.wgt = wts; /* Input weights on Rrs values       */
    ambdata.yfit = Rrs_fit; /* Output model predicted Rrs values */
    ambdata.meta = g;

    if (firstCall == 1) {
        firstCall = 0;
        if ((init = (double *) calloc(g->nwave * (g->nwave + 1), sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP:gio_model.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

    }
    /* initialize simplex with first guess model parameters */
    for (j = 0; j < g->npar + 1; j++)
        for (i = 0; i < g->npar; i++)
            init[j * g->npar + i] = g->par[i];

    for (i = 0; i < g->npar; i++) {
        init[g->npar + i * (g->npar + 1)] += g->len[i];
        par[i] = 0.0;
    }

    /* run optimization */
    isml = amoeba(init, &ambdata, giop_amb, tol);

    /* check convergence and record parameter results */
    if (ambdata.niter >= g->maxiter)
        status = 1;

    for (i = 0; i < g->npar; i++) {
        par[i] = init[g->npar * isml + i];
    }

    *itercnt = ambdata.niter;

    return (status);
}



/* ---------------------------------------------------------------------- */
/* wrapper function for L-M fit to GIOP model                             */

/* ---------------------------------------------------------------------- */
int giop_lm_fdf(const gsl_vector *parv, void *data, gsl_vector *f, gsl_matrix *J) {
    double *y = ((fitstr *) data)->y;
    double *w = ((fitstr *) data)->w;
    float *aw = ((fitstr *) data)->g->aw;
    float *bbw = ((fitstr *) data)->g->bbw;
    float *foq = ((fitstr *) data)->g->foq;
    float *wave = ((fitstr *) data)->g->wave;
    size_t nwave = ((fitstr *) data)->g->nwave;
    size_t npar = ((fitstr *) data)->g->npar;
    giopstr *g = ((fitstr *) data)->g;

    double *par;
    double *yfit;
    double **dydpar;

    double sigma;

    int iw, ipar;

    if ((par = (double *) calloc(nwave, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:giop_lm_fdf.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    /* extract model params from vector */
    for (ipar = 0; ipar < npar; ipar++) {
        par[ipar] = gsl_vector_get(parv, ipar);
    }
    if ((yfit = (double *) calloc(nwave, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:giop_lm_fdf.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((dydpar = (double **) calloc(nwave, sizeof (double *))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:giop_lm_fdf.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    for (iw = 0; iw < nwave; iw++)
        if ((dydpar[iw] = (double *) calloc(nwave, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP:giop_lm_fdf.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

    /* run model */
    giop_model(g, par, nwave, wave, aw, bbw, foq, aph1, adg1, bbp1, yfit, dydpar, NULL);

    /* evaluate and store for lm solver */
    for (iw = 0; iw < nwave; iw++) {

        /* silly, but we need sigma and not sigma-squared */
        sigma = sqrt(1. / w[iw]);

        /* Function to be minimized */
        if (f != NULL)
            gsl_vector_set(f, iw, (yfit[iw] - y[iw]) / sigma);

        /* Jacobian matrix */
        if (J != NULL)
            for (ipar = 0; ipar < npar; ipar++)
                gsl_matrix_set(J, iw, ipar, dydpar[iw][ipar] / sigma);
    }

    free(yfit);
    freeDArray(dydpar, nwave);
    free(par);
    return GSL_SUCCESS;
}

int giop_lm_f(const gsl_vector *parv, void *data, gsl_vector *f) {
    return (giop_lm_fdf(parv, data, f, NULL));
}

int giop_lm_df(const gsl_vector *parv, void *data, gsl_matrix *J) {
    return (giop_lm_fdf(parv, data, NULL, J));
}

/* ---------------------------------------------------------------------- */
/*calc_uadg_s() - calculate uncertainty in model ed adg_s slope           */

/* ---------------------------------------------------------------------- */

float calc_uadg_s(giopstr *g, float Rrs1, float Rrs2, float uRrs1, float uRrs2, float covRrs1Rrs2) {

    float dsdr1,dsdr2;
    float uadg_s;
    
    switch (g->adg_opt) {
    case ADGSQAA:
        dsdr1 = -0.002 / (0.36*Rrs2 + 1.2*Rrs1 + Rrs1*Rrs1/Rrs2);
        dsdr2 = (0.002*Rrs1) / pow(Rrs1 + 0.6*Rrs2,2);
        uadg_s = pow(pow(dsdr1*uRrs1,2) + pow(dsdr2*uRrs2,2) + 2*dsdr1*dsdr2*covRrs1Rrs2, 0.5);
        break;
    case ADGSOBPG:
        dsdr1 = 0.0038/(log(10)*Rrs1);
        dsdr2 = -0.0038/(log(10)*Rrs2);
        uadg_s = pow(pow(dsdr1*uRrs1,2) + pow(dsdr2*uRrs2,2) +2*dsdr1*dsdr2*covRrs1Rrs2, 0.5);
        break;
    } 
    return uadg_s; 
}

/* ---------------------------------------------------------------------- */
/*calc_sbbp_unc() - calculate uncertainty in the modelled bbp_s slope     */

/* ---------------------------------------------------------------------- */

float calc_ubbp_s(giopstr *g, float Rrs1, float Rrs2, float uRrs1, float uRrs2, float covRrs1Rrs2) {

    float v,dvdr1,dvdr2,dsdr1,dsdr2,dsdv;
    float dsdchl;
    float ubbp_s;
    
    switch (g->bbp_opt) {
    case BBPSHAL:
        dsdr1 = 0.8/Rrs1;
        dsdr2 = -0.8*Rrs1/pow(Rrs2,2);
        ubbp_s = pow( pow(dsdr1*uRrs1,2) + pow(dsdr2*uRrs2,2) + 2*dsdr1*dsdr2*covRrs1Rrs2, 0.5);
        break;
    case BBPSQAA:
        v = -0.9*(Rrs1/Rrs2);
        dsdv = -2.4*exp(v);
        dvdr1 = -0.9/Rrs1;
        dvdr2 = 0.9*(Rrs1/pow(Rrs2,2));
        dsdr1 = dsdv*dvdr1;
        dsdr2 = dsdv*dvdr2;
        ubbp_s = pow(pow(dsdr1*uRrs1,2) + pow(dsdr2*uRrs2,2) + 2*dsdr1*dsdr2*covRrs1Rrs2, 0.5);
        break;
    case BBPSCIOTTI:
        dsdchl = -0.768/(log(10)*g->chl); 
        ubbp_s = pow(pow( dsdchl*g->uchl,2) ,0.5);
        break;
    case BBPSMM01:
        dsdchl = -0.5/(log(10)*g->chl);
        ubbp_s = pow(pow( dsdchl*g->uchl,2) ,0.5);
        break;
    }
        
    return ubbp_s;
} 

/* ---------------------------------------------------------------------- */
/*calc_pinv() - calculate the Moore-Penrose pseudo inverse of a matrix    */

/* ---------------------------------------------------------------------- */

gsl_matrix* calc_pinv(gsl_matrix *A, const realtype rcond) {

    gsl_matrix *V, *Sigma_pinv, *U, *A_pinv;
    gsl_matrix *_tmp_mat = NULL;
    gsl_vector *_tmp_vec;
    gsl_vector *u;
    realtype x, cutoff;
    size_t i, j;
    unsigned int n = A->size1;
    unsigned int m = A->size2;
    bool was_swapped = false;


    if (m > n) {
            /* libgsl SVD can only handle the case m <= n - transpose matrix */
            was_swapped = true;
            _tmp_mat = gsl_matrix_alloc(m, n);
            gsl_matrix_transpose_memcpy(_tmp_mat, A);
            A = _tmp_mat;
            i = m;
            m = n;
            n = i;
    }

    /* do SVD */
    V = gsl_matrix_alloc(m, m);
    u = gsl_vector_alloc(m);
    _tmp_vec = gsl_vector_alloc(m);
    gsl_linalg_SV_decomp(A, V, u, _tmp_vec);
    gsl_vector_free(_tmp_vec);

    /* compute Σ⁻¹ */
    Sigma_pinv = gsl_matrix_alloc(m, n);
    gsl_matrix_set_zero(Sigma_pinv);
    cutoff = rcond * gsl_vector_max(u);

    for (i = 0; i < m; ++i) {
            if (gsl_vector_get(u, i) > cutoff) {
                    x = 1. / gsl_vector_get(u, i);
            }
            else {
                    x = 0.;
            }
            gsl_matrix_set(Sigma_pinv, i, i, x);
    }

    /* libgsl SVD yields "thin" SVD - pad to full matrix by adding zeros */
    U = gsl_matrix_alloc(n, n);
    gsl_matrix_set_zero(U);

    for (i = 0; i < n; ++i) {
            for (j = 0; j < m; ++j) {
                    gsl_matrix_set(U, i, j, gsl_matrix_get(A, i, j));
            }
    }

    if (_tmp_mat != NULL) {
            gsl_matrix_free(_tmp_mat);
    }

    /* two dot products to obtain pseudoinverse */
    _tmp_mat = gsl_matrix_alloc(m, n);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., V, Sigma_pinv, 0., _tmp_mat);

    if (was_swapped) {
            A_pinv = gsl_matrix_alloc(n, m);
            gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., U, _tmp_mat, 0., A_pinv);
    }
    else {
            A_pinv = gsl_matrix_alloc(m, n);
            gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., _tmp_mat, U, 0., A_pinv);
    }

    gsl_matrix_free(_tmp_mat);
    gsl_matrix_free(U);
    gsl_matrix_free(Sigma_pinv);
    gsl_vector_free(u);
    gsl_matrix_free(V);

    return A_pinv;
}


/* ---------------------------------------------------------------------- */
/*par_unc_calc() - compute uncertainties in the GIOP model free parameters*/
/*                                                                        */
/* Reference:                                                             */
/* McKinna et al. (2019) AApproach for Propagating Radiometric Data       */
/* Uncertainties Through NASA Ocean Color Algorithms, Front. Earth Sci.,  */ 
/* 7, doi: 10.3389/feart.2019.00176.                                      */
/*                                                                        */
/* Implementation: L. McKinna, GO2Q/NASA GSFC, October 2023               */

/* ---------------------------------------------------------------------- */
void calc_par_unc(giopstr *g, double par[], double uRrs[], double upar[]) {
    
    int iw, ix, iy;

    size_t nwave = g->nwave;
    int32_t m = g->nwave;
    int32_t n = g->npar;
    
    double *Rrs_f;
    double **dydpar;
    
    gsl_matrix *covRrs = gsl_matrix_alloc(m,m);
    gsl_matrix *Jac = gsl_matrix_alloc(m,n); 
    gsl_matrix *invJac = gsl_matrix_alloc(n,m);  
    gsl_matrix *tmpMatrix = gsl_matrix_alloc(m,n);
    gsl_matrix *covPar = gsl_matrix_alloc(n,n);
   
    if ((Rrs_f = (double *) calloc(nwave, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:giop_lm_fdf.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    
    if ((dydpar = (double **) calloc(nwave, sizeof (double *))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:giop_lm_fdf.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    
    for (iw = 0; iw < nwave; iw++) {
        if ((dydpar[iw] = (double *) calloc(nwave, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP:giop_lm_fdf.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    }

    
    /* run model */
    switch (g->fit_opt) {
    case SVDSIOP:
        giop_model_iterate(g, par, g->nwave, g->wave, g->aw, g->bbw, g->foq, aph1, adg1, bbp1, acdom1, anap1, bbph1, bbnap1, Rrs_f, dydpar, NULL);
        break;
    default:
        giop_model(g, par, g->nwave, g->wave, g->aw, g->bbw, g->foq, aph1, adg1, bbp1, Rrs_f, dydpar, NULL);
        break;
    }
  
    /* evaluate and store for lm solver */
    for (ix = 0; ix < g->nwave; ix++) {
        for (iy = 0; iy < g->npar; iy++) {
            gsl_matrix_set(Jac, ix, iy, dydpar[ix][iy]);
        }
    }
    
    //Set the error-covariance matrix with Rrs variance diagonal terms
    for (ix = 0; ix < g->nwave; ix++) {
        for (iy = 0; iy < g->nwave; iy++) {
            if (ix == iy) {
               gsl_matrix_set(covRrs, ix, iy, pow(uRrs[g->bindx[ix]],2));  
            } else {
               gsl_matrix_set(covRrs, ix, iy, 0.0);  
            }
        } 
    }
    
    invJac = calc_pinv(Jac,1E-15);
    
    gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, covRrs, invJac, 0.0, tmpMatrix);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, invJac, tmpMatrix, 0.0, covPar);
    
    for (ix = 0; ix < g->npar; ix++){
        for (iy = 0; iy < g->npar; iy++){
            if (ix == iy) {
                upar[ix] = pow(gsl_matrix_get(covPar,ix,iy),0.5);
            }
        }
    }
    
    free(Rrs_f);
    freeDArray(dydpar, nwave);
    gsl_matrix_free(Jac);
    gsl_matrix_free(covPar); 
    gsl_matrix_free(invJac);
    gsl_matrix_free(covRrs);
    gsl_matrix_free(tmpMatrix);
    
}

/* ---------------------------------------------------------------------- */
/* fit_giop_lm() - runs Levenburg-Marquart optimization for one pixel     */

/* ---------------------------------------------------------------------- */
int fit_giop_lm(giopstr *g, double Rrs[], double wts[], double par[], double *chi, int16 *itercnt) {
    int status = 0;
    int iter;

    static fitstr data;

    size_t npar = g->npar;

    static gsl_multifit_function_fdf func;
    const gsl_multifit_fdfsolver_type *t;
    gsl_multifit_fdfsolver *s;
    gsl_vector_view x;

    gsl_matrix *J = gsl_matrix_alloc(g->nwave, npar);
    gsl_matrix *cov = gsl_matrix_alloc(npar, npar);

    double sum, dof, c;
    int ipar;

    /* Set up data structure */
    data.y = Rrs;
    data.w = wts;
    data.g = g;

    /* Set up multifit function structure */
    func.f = &giop_lm_f;
    func.df = &giop_lm_df;
    func.fdf = &giop_lm_fdf;
    func.n = g->nwave;
    func.p = g->npar;
    func.params = &data;

    /* Allocate solver space */
    t = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc(t, g->nwave, g->npar);

    /* Set start params */
    x = gsl_vector_view_array(par, npar);
    gsl_multifit_fdfsolver_set(s, &func, &x.vector);

    /* Fit model for this pixel */
    status = 1;
    for (iter = 0; iter < g->maxiter; iter++) {
        gsl_multifit_fdfsolver_iterate(s);
        if (gsl_multifit_test_delta(s->dx, s->x, 1e-4, 1e-4) == GSL_SUCCESS) {
            status = 0;
            break;
        }
    }
    *itercnt = iter;

    /* Compute covariance matrix  from Jacobian */
    gsl_multifit_fdfsolver_jac(s, J);
    gsl_multifit_covar(J, 0.0, cov);

    /* Compute chi-square */
    sum = gsl_blas_dnrm2(s->f); // quadrature sum of minimization function
    dof = func.n - func.p; // degrees of freedom
    *chi = pow(sum, 2.0) / dof; // chi-square per dof

    if (g->wt_opt == 0)
        c = sum / sqrt(dof); // use variance in fit to estimate Rrs uncertainty
    else
        c = 1.0; // Rrs uncertainty included in sum 

    /* Retrieve fit params & error estimates */
    for (ipar = 0; ipar < npar; ipar++) {
        par[ipar] = gsl_vector_get(s->x, ipar);
        par[ipar + npar] = sqrt(gsl_matrix_get(cov, ipar, ipar)) * c;
        //printf("%d %lf %lf\n",ipar,par[ipar],par[ipar+npar]);
    }

    gsl_multifit_fdfsolver_free(s);
    gsl_matrix_free(cov);
    gsl_matrix_free(J);

    return (status);
}


/* ---------------------------------------------------------------------- */
/* fit_giop_svd() - constrained matrix solution                           */

/* ---------------------------------------------------------------------- */
int fit_giop_svd(giopstr *g, double rrs[], double wts[], double par[]) {
    int status = 0; /* init to success */

    double *rrs_fit;
    double **parstar;

    size_t nwave = g->nwave;
    size_t npar = g->npar;

    int iw, ipar;

    double a0, a1, a2;
    double u0, u1, u;

    gsl_matrix *A = gsl_matrix_alloc(nwave, npar);
    gsl_matrix *V = gsl_matrix_alloc(npar, npar);
    gsl_vector *S = gsl_vector_alloc(npar);
    gsl_vector *W = gsl_vector_alloc(npar);
    gsl_vector *x = gsl_vector_alloc(npar);
    gsl_vector *b = gsl_vector_alloc(nwave);

    if ((rrs_fit = (double *) calloc(nwave, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:fit_giop_svd.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((parstar = (double **) calloc(nwave, sizeof (double *))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:fit_giop_svd.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    for (iw = 0; iw < nwave; iw++)
        if ((parstar[iw] = (double *) calloc(nwave, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP:fit_giop_svd.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

    /* run model to get parstar wavelength-dependence terms */
    giop_model(g, par, g->nwave, g->wave, g->aw, g->bbw, g->foq, aph1, adg1, bbp1, rrs_fit, NULL, parstar);

    /* clear fit parameters */
    for (ipar = 0; ipar < g->npar; ipar++)
        par[ipar] = BAD_FLT;

    /* build A matrix and b vector for A x = b */

    for (iw = 0; iw < g->nwave; iw++) {

        /* get u = bb/(a+bb) */
        switch (g->rrs_opt) {
        case RRSGRD:
            a2 = g->grd[1];
            a1 = g->grd[0];
            a0 = -rrs[iw];
            if ((gsl_poly_solve_quadratic(a2, a1, a0, &u0, &u1) == 2) && u1 > 0.0)
                u = u1;
            else {
                status = 1;
                goto cleanup;
            }
            break;
        case RRSFOQ:
            u = rrs[iw] / g->foq[iw];
            break;
        }

        gsl_vector_set(b, iw, -(u * g->aw[iw] + (u - 1.0) * g->bbw[iw]));

        for (ipar = 0; ipar < g->npar; ipar++) {
            if (ipar < (g->aph_nvec + g->adg_nvec))
                gsl_matrix_set(A, iw, ipar, parstar[iw][ipar] * u); // a
            else
                gsl_matrix_set(A, iw, ipar, parstar[iw][ipar]*(u - 1.0)); // bb
        }
    }

    /* solve A x = b for x */
    status = gsl_linalg_SV_decomp(A, V, S, W);
    status = gsl_linalg_SV_solve(A, V, S, b, x);

    /* extract fitted parameters */
    if (status == 0) {
        for (ipar = 0; ipar < g->npar; ipar++) {
            par[ipar] = gsl_vector_get(x, ipar);
        }
    }

cleanup:

    gsl_matrix_free(A);
    gsl_matrix_free(V);
    gsl_vector_free(S);
    gsl_vector_free(x);
    gsl_vector_free(b);

    free(rrs_fit);
    freeDArray(parstar, nwave);
    return (status);
}

/* ---------------------------------------------------------------------- */
/* fit_giop_svd_siop() - adaptive linear matrix inversion solution using  */
/*                       specific inherent optical properties (SIOPS) and */
/*                       SVD to solve system of equations                 */
/*                                                                        */
/* Reference:                                                             */
/* Brando et al. (2012) Adaptive semianalytical inversion of ocean color  */
/* radiometry in optically complex waters, Applied Optics, 51(15),        */
/* 2808-2833.                                                             */
/*                                                                        */
/* Implementation: L. McKinna, SAIC/NASA GSFC, November 2015              */
/*                                                                        */

/* ---------------------------------------------------------------------- */
int fit_giop_svd_siop(giopstr *g, double rrs[], double wts[], double par[], double *chi) {

    int16 status = 0; /* init to success */

    double *rrs_fit;
    double **parstar;
    double *parArr;
    double *parFit;
    double *rmse;
    int *badSolution;


    size_t nwave = g->nwave;
    size_t npar = g->npar; //Number of parameters fixed at 3
    size_t nvec = g->aph_nvec;

    int iw, ipar, iv, smlIdx;

    double a0, a1, a2;
    double u0, u1, u;

    double diffSq, diffSqSum, smlRmse, sumRrs;

    gsl_matrix *A = gsl_matrix_alloc(nwave, npar);
    gsl_matrix *V = gsl_matrix_alloc(npar, npar);
    gsl_vector *S = gsl_vector_alloc(npar);
    gsl_vector *W = gsl_vector_alloc(npar);
    gsl_vector *x = gsl_vector_alloc(npar);
    gsl_vector *b = gsl_vector_alloc(nwave);


    if ((rrs_fit = (double *) calloc(nwave, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:fit_giop_svd.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((parstar = (double **) calloc(nwave, sizeof (double *))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:fit_giop_svd.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    for (iw = 0; iw < nwave; iw++) {
        if ((parstar[iw] = (double *) calloc(5, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP:fit_giop_svd.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    }
    if ((parArr = (double *) calloc(npar * nvec, sizeof (double *))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:fit_giop_svd.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((rmse = (double *) calloc(nvec, sizeof (double *))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:fit_giop_svd.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((badSolution = (int *) calloc(nvec, sizeof (int *))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:fit_giop_svd.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((parFit = (double *) calloc(npar, sizeof (double *))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:fit_giop_svd.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    /* intialize fit parameters */
    for (ipar = 0; ipar < g->npar; ipar++)
        par[ipar] = BAD_FLT;

    for (ipar = 0; ipar < 3; ipar++)
        parFit[ipar] = BAD_FLT;


    /* clear fit parameters */
    for (iv = 0; iv < nvec; iv++)
        badSolution[iv] = BAD_INT;


    /*Iterate over siop combinations*/
    for (iv = 0; iv < nvec; iv++) {

        g->siopIdx = iv;

        /* run model to get parstar wavelength-dependence terms */
        giop_model_iterate(g, parFit, g->nwave, g->wave, g->aw, g->bbw, g->foq, aph1, adg1, bbp1, acdom1, anap1, bbph1, bbnap1, rrs_fit, NULL, parstar);

        /* build A matrix and b vector for A x = b */
        for (iw = 0; iw < g->nwave; iw++) {

            /* get u = bb/(a+bb) */
            switch (g->rrs_opt) {
            case RRSGRD:
                a2 = g->grd[1];
                a1 = g->grd[0];
                a0 = -rrs[iw];
                if ((gsl_poly_solve_quadratic(a2, a1, a0, &u0, &u1) == 2) && u1 > 0.0)
                    u = u1;
                else {
                    status = 1;
                    goto cleanup;
                }
                break;
            case RRSFOQ:
                u = rrs[iw] / g->foq[iw];
                break;
            }

            gsl_vector_set(b, iw, -(g->aw[iw] * u + g->bbw[iw]*(u - 1.0)));

            // SIOP implementation Aij in eq 12 of Brando et al 2012
            // A is 3 components + water	
            // i=0 (phy): ipar=0 for aphy* and ipar 3 for bbphy*
            // i=1 (CDOM): ipar=1 for aCDOM* 
            // i=2 (NAP): ipar=2 for anap* and ipar 4 for bbnap*    
            gsl_matrix_set(A, iw, 0, parstar[iw][0] * u + parstar[iw][3]*(u - 1.0)); // a+ bb
            gsl_matrix_set(A, iw, 1, parstar[iw][1] * u); // ONLY a
            gsl_matrix_set(A, iw, 2, parstar[iw][2] * u + parstar[iw][4]*(u - 1.0)); // a+ bb 

        }

        /* solve A x = b for x */
        /*SVD decomposition*/
        status = gsl_linalg_SV_decomp(A, V, S, W);
        status = gsl_linalg_SV_solve(A, V, S, b, x);


        for (ipar = 0; ipar < npar; ipar++) {

            /*Test for non-physical (negative) concentrations or matrix inversion*/
            /*failure. For such circumstances set to BAD_FLT*/
            if (gsl_vector_get(x, ipar) < 0.0 || status != 0) {
                badSolution[iv] = 1;
                for (ipar = 0; ipar < npar; ipar++) {
                    parArr[iv * npar + ipar] = BAD_FLT;
                }
                break;

            } else if (gsl_vector_get(x, ipar) >= 0.0) {
                badSolution[iv] = 0;
                parArr[iv * npar + ipar] = gsl_vector_get(x, ipar);
            }
        }
    }

    /*Find optimal set of parameters*/
    g->siopIdx = 0;
    diffSq = 0;
    diffSqSum = 0;

    /*Compute the root mean square error (rmse) metric for all SIOP combinations. */
    /*The smallest rmse valid values is returned as the  optimal solution.*/
    for (iv = 0; iv < nvec; iv++) {

        g->siopIdx = iv;

        for (ipar = 0; ipar < npar; ipar++) {
            //if (badSolution[iv]) {
            //    parFit[ipar] = BAD_FLT;
            //} else {
            parFit[ipar] = parArr[iv * npar + ipar];
            //}    
        }

        /*Calculate fitted Rrs*/
        giop_model_iterate(g, parFit, g->nwave, g->wave, g->aw, g->bbw, g->foq, aph1, adg1, bbp1, acdom1, anap1, bbph1, bbnap1, rrs_fit, NULL, parstar);

        /*reset the sum of squares values*/
        diffSq = 0;
        diffSqSum = 0;
        sumRrs = 0;

        for (iw = 0; iw < g->nwave; iw++) {
            diffSq = pow((rrs_fit[iw] - rrs[iw]), 2);
            diffSqSum += diffSq;
            sumRrs += rrs[iw];
        }


        if (badSolution[iv]) {
            rmse[iv] = 999;
        } else {
            /*Save value to relative RMSE array*/
            rmse[iv] = pow((diffSqSum / (g->nwave - 1)), 0.5); // sumRrs/(g->nwave - 1) );
        }
    }

    /*Find the smallest relative RMSE value*/
    /*Initialize with zeroth values*/
    smlRmse = rmse[0];
    smlIdx = 0;

    for (iv = 0; iv < nvec; iv++) {
        if (rmse[iv] < smlRmse) {
            smlRmse = rmse[iv];
            smlIdx = iv;
        }
    }


    //Save index of best SIOP combination to giop record
    g->siopIdx = smlIdx;

    /*Get optimal concentrations values from array of all concentrations*/
    /*If there is still a negative solution, set it to zero*/
    if (rmse[smlIdx] != 999 && badSolution[smlIdx] != 1) {

        for (ipar = 0; ipar < npar; ipar++) {
            par[ipar] = parArr[smlIdx * npar + ipar];
        }
        *chi = smlRmse;
        status = 0;


    } else {
        //No solution status
        for (ipar = 0; ipar < npar; ipar++) {
            par[ipar] = BAD_FLT;
        }
        *chi = BAD_FLT;
        status = -99;

    }

cleanup:
    gsl_matrix_free(A);
    gsl_matrix_free(V);
    gsl_vector_free(S);
    gsl_vector_free(x);
    gsl_vector_free(b);


    free(rrs_fit);
    free(parArr);
    free(rmse);
    free(badSolution);
    freeDArray(parstar, nwave);
    return (status);
}

/* ------------------------------------------------------------------- */
/* giop_chl() - returns magnitude of aph* as chl                       */

/* ------------------------------------------------------------------- */
float32 giop_chl(giopstr *g, int16 iopf, double par[], float *uchl) {
    float32 chl = badval;
    int ipar;

    float uchl2_temp = badval;
    *uchl = badval;

    switch (g->aph_opt) {
    case APHGAUSS:
        break;
    default:
        if ((iopf & IOPF_FAILED) == 0) {
            chl = 0.0;
            uchl2_temp = 0.0;
            for (ipar = 0; ipar < g->aph_nvec; ipar++) {
                chl += par[ipar];
                uchl2_temp += pow(par[ipar + g->npar],2);
            }
            *uchl = pow(uchl2_temp,0.5);
        }
        break;
    }


    return (chl);
}

/* ---------------------------------------------------------------------- */
/* Convert Rrs[0+] to Rrs[0-]                                             */

/* ---------------------------------------------------------------------- */
float rrs_above_to_below(float Rrs) {
    return (Rrs / (0.52 + 1.7 * Rrs));

}

/* ---------------------------------------------------------------------- */
/* Uncertainty estimate convert Rrs[0+] to Rrs[0-]                         */

/* ---------------------------------------------------------------------- */
float rrs_above_to_below_unc(float Rrs, float uRrs) {
    float dydr;
    dydr = 0.52 / pow(0.52 + 1.7*Rrs, 2.);
    return pow(pow(dydr*uRrs,2.),0.5);

}

/* ---------------------------------------------------------------------- */
/* Convert Rrs[0-] to Rrs[0+]                                             */

/* ---------------------------------------------------------------------- */
float rrs_below_to_above(float rrs_s) {
    return ( (0.52 * rrs_s) / (1.0 - 1.7 * rrs_s));

}

/* ---------------------------------------------------------------------- */
/* Uncertainty estimate convert Rrs[0-] to Rrs[0+]                        */

/* ---------------------------------------------------------------------- */
float rrs_below_to_above_unc(float rrs_s, float urrs_s) {
    float dydr;
    dydr = 0.52/pow(1. - 1.7*rrs_s,2.);
    return pow(pow(dydr*urrs_s,2.), 0.5);

}

/* ---------------------------------------------------------------------- */
/* get_bbp_lh computes spectral backscattering coefficient using the      */
/* reflectance line height model of McKinna et al (2021)                  */
/* -----------------------------------------------------------------------*/
int get_bbp_lh(l2str *l2rec, giopstr *g, int ipb2, float tab_wave[], 
        float tab_bbp[], float tab_ubbp[], int tab_nwave) {
    
    int iw;
    int i1, i2, i3;
    float Rrs1, Rrs2, Rrs3;
    float uRrs1, uRrs2, uRrs3, covRrs1Rrs2, covRrs1Rrs3, covRrs2Rrs3;
    float LH, uLH;
    float *uRrs2_new;
    uRrs2_new = (float*) calloc(1, sizeof(float)); 
    
    float bbp555_lh = badval;
    float ubbp555_lh = 0;
    
    float *wave = l2rec->l1rec->l1file->fwave;
    int32 nwave = l2rec->l1rec->l1file->nbands;
    
    static int status = -1;
    
    //intialize arrays
    for (iw = 0; iw < tab_nwave; iw++) {
        g->bbp_tab_w[iw] = wave[iw];
        tab_ubbp[iw] = badval;
    }   
    
    //LH backscattering coefficients
    float alh[2] = {-2.5770 , 281.27}; 
    //LH backscattering uncertainties and covariance
    float ualh[3] = {2.4819E-2, 20.777, 0.24852};
    
    i1 = windex(490.0, wave, nwave);
    i2 = windex(555.0, wave, nwave);
    i3 = windex(670.0, wave, nwave);

    Rrs1 = l2rec->Rrs[ipb2 + i1];
    Rrs2 = l2rec->Rrs[ipb2 + i2];
    Rrs3 = l2rec->Rrs[ipb2 + i3]; 
    
    uRrs1 = g->uRrs_a[i1];
    uRrs2 = g->uRrs_a[i2];
    uRrs3 = g->uRrs_a[i3];

    covRrs1Rrs2 = 0.0;
    covRrs1Rrs3 = 0.0;
    covRrs2Rrs3 = 0.0;

    Rrs2 = conv_rrs_to_555(Rrs2, l2rec->l1rec->l1file->fwave[i2], uRrs2, uRrs2_new);
    //conv_rrs_to_555(float Rrs, float wave, float uRrs_in, float *uRrs_out ) 
    //uRrs2 = uconv_rrs_to_555(Rrs2, uRrs2, l2rec->l1rec->l1file->fwave[i2]);

    LH = Rrs2 - (0.64*Rrs1 + 0.36*Rrs3  );
    uLH = sqrt( pow(0.64*uRrs1,2) + pow(*uRrs2_new,2) + 
            pow(0.36*uRrs3,2) - 2*0.64*covRrs1Rrs2  - 
            2*0.64*0.36*covRrs1Rrs3 - 0.36*covRrs2Rrs3);

    //Compute bbp(555) and uncertainty using an empirical model
    bbp555_lh = pow(10.,alh[0] + alh[1]*LH);
    
    ubbp555_lh = sqrt(pow(log(10)*ualh[0]*bbp555_lh,2) + 
            pow(LH*log(10)*ualh[1]*bbp555_lh,2) + 
            pow(uLH*alh[1]*log(10)*bbp555_lh,2) +
            2*ualh[2]*pow(log(10)*bbp555_lh,2));
    
    if (bbp555_lh < 0) {
        return(0);
    } else { 
        status = 1;
    }
    
    //Compute bbp_s according to Lee et al (2002)
    // update power-law exponent based on QAA band-ratio relationship
    i1 = windex(443.0, wave, nwave);
    i2 = windex(550.0, wave, nwave);
    //QAA bbp_s relationship derived from NOMAD data. 
    //No Raman scattering correction applied to Rrs.
    Rrs1 = l2rec->Rrs[ipb2 + i1];
    Rrs2 = l2rec->Rrs[ipb2 + i2];
    uRrs1 = g->uRrs_a[i1];
    uRrs2 = g->uRrs_a[i2];
    covRrs1Rrs2 = 0.0;                                          

    if (Rrs1 > 0.0 && Rrs2 > 0.0) {
        Rrs1 = rrs_above_to_below(Rrs1);
        Rrs2 = rrs_above_to_below(Rrs2);
        uRrs1 = rrs_above_to_below_unc(Rrs1, uRrs1); 
        uRrs2 = rrs_above_to_below_unc(Rrs2, uRrs2);
        g->bbp_s = MAX(MIN(2.0 * (1.0 - 1.2 * exp(-0.9 * Rrs1 / Rrs2)), bbp_s_max), bbp_s_min);
        if (g->bbp_s != bbp_s_min || g->bbp_s != bbp_s_max) {
            g->ubbp_s = calc_ubbp_s(g, Rrs1, Rrs2, uRrs1, *uRrs2_new, covRrs1Rrs2);
        } else {
            g->ubbp_s = 0.0;
        }
    } else if (Rrs2 > 0.0) {
        g->bbp_s = bbp_s_min;
    } else {
        g->bbp_s = BAD_FLT;
        status = 0;
    }
    
    if (status == 0) {
        for (iw = 0; iw < tab_nwave; iw++) {
            g->bbp_tab_w[iw] = wave[iw];
            tab_bbp[iw] = 0;
            tab_ubbp[iw] = 0;
        }
    } else {
        for (iw = 0; iw < tab_nwave; iw++) {
            g->bbp_tab_w[iw] = wave[iw];
            tab_bbp[iw] = bbp555_lh*pow(555./tab_wave[iw], g->bbp_s);
            tab_ubbp[iw] = sqrt(
                    pow(bbp555_lh*g->ubbp_s*pow(555./tab_wave[iw], g->bbp_s)*log(555./tab_wave[iw]),2)
                    + pow(ubbp555_lh*pow(555./tab_wave[iw], g->bbp_s),2));
        }
    }
    
    free(uRrs2_new);
    return status;
}

/* ----------------------------------------------------------------------  */
/* get_bbp_chl computes spectral backscattering coefficient using the      */
/* chlorophyll-dependent model of Huot et al (2008)                        */
/* ----------------------------------------------------------------------- */
int get_bbp_chl(l2str *l2rec, giopstr *g, int ipb2, float tab_wave[], 
        float tab_bbp[], float tab_ubbp[], int tab_nwave) {
    
    int iw;
    int i1, i2;
    float Rrs1, Rrs2;
    float uRrs1, uRrs2,covRrs1Rrs2;
   
    float bbp555_chl = badval;
    float ubbp555_chl = 0;
    
    float *wave = l2rec->l1rec->l1file->fwave;
    int32 nwave = l2rec->l1rec->l1file->nbands;
    
    static int status = -1;
    
    //intialize arrays
    for (iw = 0; iw < tab_nwave; iw++) {
        g->bbp_tab_w[iw] = wave[iw];
        tab_ubbp[iw] = badval;
    }   
    
    //Huot backscattering coefficients
    float achl[2] = {0, 0}; 
    //Huot backscattering uncertainties and covariance
    float uachl[2] = {1E-4, 0.02};
    
    i1 = windex(550.0, wave, nwave);
    
    achl[0] = 2.267E-3 - 5.058E-6 * (l2rec->l1rec->l1file->fwave[i1] - 550.);
    achl[1] = 0.565 - 4.86E-4 *  (l2rec->l1rec->l1file->fwave[i1] - 550.);

    //Compute bbp(555) and uncertainty estimate using an empirical model
    bbp555_chl = achl[0]*pow(g->chl,achl[1]);
    ubbp555_chl = sqrt( g->uchl*pow(achl[0]*achl[1]*pow(g->chl,(achl[1]-1)),2) +
            pow(uachl[0]*pow(g->chl,achl[1]),2) +
            pow(uachl[1]* log(g->chl)*achl[0]*pow(g->chl,achl[1]),2));
    
    
    if (bbp555_chl < 0) {
        return(0);
    } else { 
        status = 1;
    }
    
    //Compute bbp_s according to Lee et al (2002)
    // update power-law exponent based on QAA band-ratio relationship
    i1 = windex(443.0, wave, nwave);
    i2 = windex(550.0, wave, nwave);
    //QAA bbp_s relationship derived from NOMAD data. 
    //No Raman scattering correction applied to Rrs.
    Rrs1 = l2rec->Rrs[ipb2 + i1];
    Rrs2 = l2rec->Rrs[ipb2 + i2];
    uRrs1 = g->uRrs_a[i1];
    uRrs2 = g->uRrs_a[i2];
    covRrs1Rrs2 = 0.0;                                          

    if (Rrs1 > 0.0 && Rrs2 > 0.0) {
        Rrs1 = rrs_above_to_below(Rrs1);
        Rrs2 = rrs_above_to_below(Rrs2);
        uRrs1 = rrs_above_to_below_unc(Rrs1, uRrs1); 
        uRrs2 = rrs_above_to_below_unc(Rrs2, uRrs2);
        g->bbp_s = MAX(MIN(2.0 * (1.0 - 1.2 * exp(-0.9 * Rrs1 / Rrs2)), bbp_s_max), bbp_s_min);
        if (g->bbp_s != bbp_s_min || g->bbp_s != bbp_s_max) {
            g->ubbp_s = calc_ubbp_s(g, Rrs1, Rrs2, uRrs1, uRrs2, covRrs1Rrs2);
        } else {
            g->ubbp_s = 0.0;
        }
    } else if (Rrs2 > 0.0) {
        g->bbp_s = bbp_s_min;
    } else {
        g->bbp_s = BAD_FLT;
        status = 0;
    }
    
    if (status == 0) {
        for (iw = 0; iw < tab_nwave; iw++) {
            g->bbp_tab_w[iw] = wave[iw];
            tab_bbp[iw] = 0;
            tab_ubbp[iw] = 0;
        }
    } else {
        for (iw = 0; iw < tab_nwave; iw++) {
            g->bbp_tab_w[iw] = wave[iw];
            tab_bbp[iw] = bbp555_chl*pow(555./tab_wave[iw], g->bbp_s);
            tab_ubbp[iw] = sqrt(
                    pow(bbp555_chl*g->ubbp_s*pow(555./tab_wave[iw], g->bbp_s)*log(555./tab_wave[iw]),2)
                    + pow(ubbp555_chl*pow(555./tab_wave[iw], g->bbp_s),2));
        }
    }
        
    return status;
}

/* ---------------------------------------------------------------------- */
/* run_giop - runs optimization using requested method                    */

/* ---------------------------------------------------------------------- */
void run_giop(l2str *l2rec) {
    static int firstCall = 1;
    static giopstr giopctl;
    static giopstr *g = &giopctl;

    double *par;
    double *upar;
    double *Rrs_a; /* above water, per fit band */
    double *uRrs_a; /* uncertainty above water, per fit band */
    double *Rrs_b; /* below water, per fit band */
    double *uRrs_b; /* uncertainty below water, per fit band */   
    double *rrs_diff; /*observed-modeled rrs diff , per fit band */ 
    double *Rrs_f; /* modeled Rrs, per fit band */
    double *wts; /* weights, per fit band     */
    double *uRrs; 

    float Rrs1, Rrs2;
    /*uncertainties in Rrs1, Rrs2, and covariance*/
    float uRrs1, uRrs2, covRrs1Rrs2; 
    int32 i1, i2;

    int16 itercnt;
    int16 bndcnt;
    int16 status;

    int32 ipar, ip, iw, ib, ipb, ipb2, ierr;
    double chi = BAD_FLT;

    l1str *l1rec = l2rec->l1rec;
    uncertainty_t *uncertainty=l1rec->uncertainty;

    float *wave = l1rec->l1file->fwave;
    int32 nwave = l1rec->l1file->nbands;
    static int32 npix;


    float aph_norm = BAD_FLT;
    float dudtau, dtaudchl,dudchl,dvdchl,daphdchl;

    if (firstCall) {
        npix = l1rec->npix;

        firstCall = 0;

        // initialize control structure (to get npar)
        if ((bbw = calloc(nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP:run_giop.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((aw = calloc(nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP:run_giop.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((foq = calloc(nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP:run_giop.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((aph1 = calloc(2 * nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP:run_giop.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((acdom1 = calloc(2 * nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP:run_giop.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((anap1 = calloc(2 * nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP:run_giop.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((adg1 = calloc(2 * nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP:run_giop.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((bbph1 = calloc(2 * nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP:run_giop.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((bbnap1 = calloc(2 * nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP:run_giop.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((bbp1 = calloc(2 * nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP:run_giop.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((g->wave = (float *) calloc(nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((g->bindx = (int *) calloc(nwave, sizeof (int))) == NULL) {
            printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((g->aw = (float *) calloc(nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((g->bbw = (float *) calloc(nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((g->wts = (float *) calloc(nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((g->foq = (float *) calloc(nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((g->par = (double *) calloc(nwave, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((g->len = (double *) calloc(nwave, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((g->Rrs_a = (float *) calloc(nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((g->uRrs_a = (float *) calloc(nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory in init_iop_flag.\n",
                    __FILE__, __LINE__);
            exit(1);
        }


        giop_ctl_init(g, nwave, wave, aw, bbw);

        // allocate static storage for one scanline

        if ((iter = calloc(npix, sizeof (int16))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((iopf = calloc(npix, sizeof (int16))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((mRrs = calloc(npix * nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((a = calloc(npix * nwave * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((a_unc = calloc(npix * nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((aph = calloc(npix * nwave * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((acdom = calloc(npix * nwave * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((anap = calloc(npix * nwave * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((adg = calloc(npix * nwave * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((bbph = calloc(npix * nwave * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((bbnap = calloc(npix * nwave * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((bbp = calloc(npix * nwave * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((bb = calloc(npix * nwave * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((bb_unc = calloc(npix * nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((chl = calloc(npix * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((rrsdiff = calloc(npix, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((aph_s = calloc(npix, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((adg_s = calloc(npix, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((uadg_s = calloc(npix, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((bbp_s = calloc(npix, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((ubbp_s = calloc(npix, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((fit_par = allocate2d_float(npix, g->npar)) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((siop_num = calloc(npix, sizeof (int))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        max_npar = g->npar;
        if ((chisqr = calloc(npix, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        /*Note: Calculation of Raman scattering contribution to Rrs is at present*/
        /*only supported in l2gen. Where Raman Rrs is not calculated (i.e. l3gen smi)*/
        /*set Rrs_raman to zeros */
        if (l2rec->Rrs_raman == NULL) {

            printf("\n");
            printf("No Raman scattering correction applied to Rrs. \n");
            printf("\n");

            allocateRrsRaman = 1;
            l2rec->Rrs_raman = (float*) allocateMemory(npix *
                    l1rec->l1file->nbands * sizeof (float), "Rrs_ram");
        }

    }

    if ((par = calloc(2 * nwave, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:run_giop.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    if ((upar = calloc(g->npar, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP:run_giop.\n",
                __FILE__, __LINE__);
        exit(1);
    }



    // reallocate if npix changes - l3gen-ism...
    if (l1rec->npix > npix) {
        npix = l1rec->npix;
        free(iter);
        free(iopf);
        free(mRrs);
        free(a);
        free(a_unc);
        free(aph);
        free(acdom);
        free(anap);
        free(adg);
        free(bbph);
        free(bbnap);
        free(bbp);
        free(bb);
        free(bb_unc);
        free(chl);
        free(rrsdiff);
        free(aph_s);
        free(adg_s);
        free(uadg_s);
        free(bbp_s);
        free(ubbp_s);
        free(fit_par);
        free(siop_num);
        free(chisqr);
        if (allocateRrsRaman) {
            free(l2rec->Rrs_raman);
            l2rec->Rrs_raman = (float*) allocateMemory(npix *
                    l1rec->l1file->nbands * sizeof (float), "Rrs_ram");
        }

        if ((iter = calloc(npix, sizeof (int16))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((iopf = calloc(npix, sizeof (int16))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((mRrs = calloc(npix * nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((a = calloc(npix * nwave * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((a_unc = calloc(npix * nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((aph = calloc(npix * nwave * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((acdom = calloc(npix * nwave * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((anap = calloc(npix * nwave * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((adg = calloc(npix * nwave * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((bbph = calloc(npix * nwave * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((bbnap = calloc(npix * nwave * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((bbp = calloc(npix * nwave * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((bb = calloc(npix * nwave * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((bb_unc = calloc(npix * nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((chl = calloc(npix * 2, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((rrsdiff = calloc(npix, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((aph_s = calloc(npix, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((adg_s = calloc(npix, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((uadg_s = calloc(npix, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((bbp_s = calloc(npix, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((ubbp_s = calloc(npix, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((fit_par = allocate2d_float(npix, g->npar)) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((siop_num = calloc(npix, sizeof (int))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        max_npar = g->npar;
        if ((chisqr = calloc(npix, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for GIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

    }

    if ((Rrs_a = calloc(nwave, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((uRrs_a = calloc(nwave, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((Rrs_b = calloc(nwave, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((uRrs_b = calloc(nwave, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((rrs_diff = calloc(nwave, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((Rrs_f = calloc(nwave, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((wts = calloc(nwave, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((uRrs = calloc(nwave, sizeof (double))) == NULL) {
        printf("-E- %s line %d : error allocating memory for GIOP.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    // re-initialize control structure when giop_iterate=1 (under evaluation)
    // PJW 9 Jan 2013
    if (input->giop_iterate > 0) {
        g->adg_s = input->giop_adg_s;
        g->bbp_s = input->giop_bbp_s;
    }


    for (ip = 0; ip < l1rec->npix; ip++) {

        ipb = ip*nwave;
        ipb2 = ip * l1rec->l1file->nbands;
        ierr = l1rec->npix * nwave + ipb;

        // initialize output arrays and counters

        for (ib = 0; ib < nwave; ib++) {
            a [ipb + ib] = badval;
            a_unc [ipb + ib] = 0.0;
            bb [ipb + ib] = badval;
            bb_unc [ipb + ib] = 0.0;
            aph [ipb + ib] = badval;
            adg [ipb + ib] = badval;
            acdom [ipb + ib] = badval;
            anap[ipb + ib] = badval;
            bbp [ipb + ib] = badval;
            bbph[ipb + ib] = badval;
            bbnap [ipb + ib] = badval;
            mRrs[ipb + ib] = badval;
            a [ierr + ib] = 0.0;
            bb [ierr + ib] = 0.0;
            aph [ierr + ib] = 0.0;
            adg [ierr + ib] = 0.0;
            acdom [ierr + ib] = 0.0;
            anap [ierr + ib] = 0.0;
            bbp [ierr + ib] = 0.0;
            bbph [ierr + ib] = 0.0;
            bbnap [ierr + ib] = 0.0;
        }
        chl [ip] = badval;
        iter[ip] = 0;
        iopf[ip] = 0;
        status = 0;
        bndcnt = 0;
        itercnt = 0;
        rrsdiff[ip] = badval;
        chisqr[ip] = badval;

        aph_s[ip] = badval;
        adg_s[ip] = badval;
        uadg_s[ip] = 0.0;
        bbp_s[ip] = badval;
        ubbp_s[ip] = 0.0;
        siop_num[ip] = BAD_INT;

        for (ipar = 0; ipar < g->npar; ipar++)
            fit_par[ip][ipar] = badval;


        // flag and skip if pixel already masked

        if (l1rec->mask[ip]) {
            iopf[ip] |= IOPF_ISMASKED;
            iopf[ip] |= IOPF_FAILED;
            continue;
        }

        // set values that depend on sea state

        for (iw = 0; iw < nwave; iw++) {
            aw [iw] = l1rec->sw_a [ipb2 + iw];
            bbw[iw] = l1rec->sw_bb[ipb2 + iw];
        }
        for (iw = 0; iw < g->nwave; iw++) {
            g->aw [iw] = aw [g->bindx[iw]];
            g->bbw[iw] = bbw [g->bindx[iw]];
        }

        //define GIOP Rrs and accsociated uncertainties 
        for (iw = 0; iw < nwave; iw++) {
            g->Rrs_a[iw] = l2rec->Rrs[ipb2 + iw];

            switch (g->urrs_opt) {
            case URRSNONE:
                g->uRrs_a[iw] =  0.0;
                break; 
            case URRSCALC:
                if (uncertainty) {
                    g->uRrs_a[iw] = l2rec->Rrs_unc[ipb2 + iw];
                } else {
                    g->uRrs_a[iw] = 0.0;
                }
                break;
            case URRSREL:
                g->uRrs_a[iw] = g->Rrs_a[iw]*input->giop_rrs_unc[iw];
                break;
            case URRSABS:
                g->uRrs_a[iw] = input->giop_rrs_unc[iw];
                break;
            default:
                g->uRrs_a[iw] = 0.0;
                break;
            }
        }
        
        
        //define giop model fit weights for select bands
        for (iw = 0; iw < g->nwave; iw++) {
            ib = g->bindx[iw];
            switch (g->urrs_opt) {
            case URRSCALC:
                if (uncertainty) {
                    g->wts[iw] = 1./pow(g->uRrs_a[g->bindx[iw]], 2);
                } else {
                    g->wts[iw] = 1.0;
                }
                break;
            case URRSREL:
            case URRSABS:
                g->wt_opt = 1;
                g->wts[iw] = 1./pow(g->uRrs_a[g->bindx[iw]], 2);
                break;
            default:
                g->wt_opt = 0;
                g->wts[iw] = 1.0;
                break;
            }
        }

        // set dynamic, non-optimized model parameters
        //

        aph_s[ip] = g->aph_s;
        adg_s[ip] = g->adg_s;
        bbp_s[ip] = g->bbp_s;
        uadg_s[ip] = g->uadg_s;
        ubbp_s[ip] = g->ubbp_s;
        

        // get starting chlorophyll from default algorithm

        g->chl = MAX(MIN(l2rec->chl[ip], chl_max), chl_min);

        //Get starting chlorophyll uncertainties for default algorithm
        //If not computed, set to zero
        if (uncertainty) {
            g->uchl = l2rec->chl_unc[ip];
        } else {
            g->uchl = 0.0;
        }

        switch (g->rrs_opt) {
        case RRSFOQ:
            foqint_morel(input->fqfile, wave, nwave, 0.0, 0.0, 0.0, g->chl, foq);
            for (iw = 0; iw < g->nwave; iw++) {
                g->foq[iw] = foq[g->bindx[iw]];
            }
            break;
        }

        // aph function

        switch (g->aph_opt) {
        case APHBRICAUD:
            // replaces default tabbed values with Bricaud chl-based values
            for (iw = 0; iw < nwave; iw++) {
                g->aph_tab_w[iw] = wave[iw];
                // get the aph* normalization factor - the 0.055 is a "typical"
                // aph* value for 443nm
                aph_norm = 0.055 / get_aphstar(443., BANDW, APHBRICAUD, g->chl);
                g->aph_tab_s[0][iw] = aph_norm * get_aphstar(wave[iw], BANDW, APHBRICAUD, g->chl);

                //Estimate aph* Bricaud standard uncertainty
                dudtau = -0.055/pow(get_aphstar(443., BANDW, APHBRICAUD, g->chl),2);
                dtaudchl = get_aphstar_pderiv(443., BANDW, APHBRICAUD, g->chl);
                dudchl = dudtau*dtaudchl;
                dvdchl = get_aphstar_pderiv(wave[iw], BANDW, APHBRICAUD, g->chl);
                daphdchl = dudchl*g->aph_tab_s[0][iw]  + dvdchl*aph_norm;
                g->uaph_tab_s[0][iw] = sqrt(pow(daphdchl*g->uchl,2));    
            }
            g->aph_tab_nw = nwave;
            g->aph_nvec = 1;
            break;
        case APHCIOTTI:
            // replaces default tabbed values with Ciotti size-fraction-based values
            for (iw = 0; iw < nwave; iw++) {
                g->aph_tab_w[iw] = wave[iw];
                g->aph_tab_s[0][iw] = get_aphstar(wave[iw], BANDW, APHCIOTTI, g->aph_s);
                g->uaph_tab_s[0][iw] = 0;
            }
            g->aph_tab_nw = nwave;
            g->aph_nvec = 1;
            break;
        }

        // adg function

        switch (g->adg_opt) {
        case ADGSQAA:
            // update exponential based on QAA band-ratio relationship
            i1 = windex(443.0, wave, nwave);
            i2 = windex(550.0, wave, nwave);
            //QAA Sdg relationship derived from NOMAD data. 
            //No Raman scattering correction applied to Rrs.
            Rrs1 = l2rec->Rrs[ipb2 + i1];
            Rrs2 = l2rec->Rrs[ipb2 + i2];
            uRrs1 = g->uRrs_a[i1];
            uRrs2 = g->uRrs_a[i2];
            covRrs1Rrs2 = 0.0;
            if (Rrs1 > 0.0 && Rrs2 > 0.0) {
                Rrs1 = rrs_above_to_below(Rrs1);
                Rrs2 = rrs_above_to_below(Rrs2);
                uRrs1 = rrs_above_to_below_unc(Rrs1,uRrs2);
                uRrs2 = rrs_above_to_below_unc(Rrs2,uRrs2);
                covRrs1Rrs2 = 0.0;
                g->adg_s = MAX(MIN(0.015 + 0.002 / (0.6 + Rrs1 / Rrs2), adg_s_max), adg_s_min);
                if (g->adg_s != adg_s_min || g->adg_s != adg_s_max) {
                    g->uadg_s = calc_uadg_s(g, Rrs1, Rrs2, uRrs1, uRrs2, covRrs1Rrs2);
                } else {
                    g->uadg_s = 0.0;
                }
            } else if (Rrs2 > 0.0) {
                g->adg_s = adg_s_min;
                g->uadg_s = 0.0;
            } else {
                g->adg_s = BAD_FLT;
                iopf[ip] |= IOPF_BADRRS;
                status = 1;
            }
            break;
        case ADGSOBPG:
            // update exponential based on OBPG band-ratio relationship
            i1 = windex(412.0, wave, nwave);
            i2 = windex(550.0, wave, nwave);
            //QAA Sdg relationship derived from NOMAD data. 
            //No Raman scattering correction applied to Rrs.
            Rrs1 = l2rec->Rrs[ipb2 + i1];
            Rrs2 = l2rec->Rrs[ipb2 + i2];
            uRrs1 = g->uRrs_a[i1];
            uRrs2 = g->uRrs_a[i2];
            covRrs1Rrs2 = 0.0;
            if (Rrs1 > 0.0 && Rrs2 > 0.0) {
                g->adg_s = MAX(MIN(0.015 + 0.0038 * log10(Rrs1 / Rrs2), adg_s_max), adg_s_min);
            if (g->adg_s != adg_s_min || g->adg_s != adg_s_max) {
                    g->uadg_s = calc_uadg_s(g, Rrs1, Rrs2, uRrs1, uRrs2, covRrs1Rrs2);
                } else {
                    g->uadg_s = 0.0;
                }
            } else if (Rrs2 > 0.0) {
                g->adg_s = adg_s_min;
                g->uadg_s = 0.0;
            } else {
                g->adg_s = BAD_FLT;
                iopf[ip] |= IOPF_BADRRS;
                status = 1;
            }
            break;
        }

        // bbp function

        switch (g->bbp_opt) {
        case BBPSHAL:
            // update power-law exponent based on HAL band-ratio relationship
            i1 = windex(490.0, wave, nwave);
            i2 = windex(550.0, wave, nwave);
            Rrs1 = l2rec->Rrs[ipb2 + i1] - l2rec->Rrs_raman[ipb2 + i1]; //CHECK IF RAMAN SHOULD BE HERE!!
            Rrs2 = l2rec->Rrs[ipb2 + i2] - l2rec->Rrs_raman[ipb2 + i2];
            uRrs1 = g->uRrs_a[i1];
            uRrs2 = g->uRrs_a[i2];
            covRrs1Rrs2 = 0.0;
            if (Rrs1 > 0.0 && Rrs2 > 0.0) {
                g->bbp_s = MAX(MIN(0.8 * (Rrs1 / Rrs2) + 0.2, bbp_s_max), bbp_s_min);
            if (g->bbp_s != bbp_s_min || g->bbp_s != bbp_s_max) {
                    g->ubbp_s = calc_ubbp_s(g, Rrs1, Rrs2, uRrs1, uRrs2, covRrs1Rrs2);
                } else {
                    g->ubbp_s = 0.0;
                }
            } else if (Rrs2 > 0.0) {
                g->bbp_s = bbp_s_min;
            } else {
                g->bbp_s = BAD_FLT;
                iopf[ip] |= IOPF_BADRRS;
                status = 1;
            }
            break;
        case BBPSQAA:
            // update power-law exponent based on QAA band-ratio relationship
            i1 = windex(443.0, wave, nwave);
            i2 = windex(550.0, wave, nwave);
            //QAA bbp_s relationship derived from NOMAD data. 
            //No Raman scattering correction applied to Rrs.
            Rrs1 = l2rec->Rrs[ipb2 + i1];
            Rrs2 = l2rec->Rrs[ipb2 + i2];
            uRrs1 = g->uRrs_a[i1];
            uRrs2 = g->uRrs_a[i2];
            covRrs1Rrs2 = 0.0;                                          
            
            if (Rrs1 > 0.0 && Rrs2 > 0.0) {
                Rrs1 = rrs_above_to_below(Rrs1);
                Rrs2 = rrs_above_to_below(Rrs2);
                uRrs1 = rrs_above_to_below_unc(Rrs1, uRrs1); 
                uRrs2 = rrs_above_to_below_unc(Rrs2, uRrs2);
                g->bbp_s = MAX(MIN(2.0 * (1.0 - 1.2 * exp(-0.9 * Rrs1 / Rrs2)), bbp_s_max), bbp_s_min);
            if (g->bbp_s != bbp_s_min || g->bbp_s != bbp_s_max) {
                    g->ubbp_s = calc_ubbp_s(g, Rrs1, Rrs2, uRrs1, uRrs2, covRrs1Rrs2);
                } else {
                    g->ubbp_s = 0.0;
                }
            } else if (Rrs2 > 0.0) {
                g->bbp_s = bbp_s_min;
            } else {
                g->bbp_s = BAD_FLT;
                iopf[ip] |= IOPF_BADRRS;
                status = 1;
            }
            break;
        case BBPSCIOTTI:
            // update power-law exponent based on Ciotti chl relationship
            if (l2rec->chl[ip] > 0.0) {
                g->bbp_s = MAX(MIN(1.0 - 0.768 * log10(g->chl), bbp_s_max), bbp_s_min);
            if (g->bbp_s != bbp_s_min || g->bbp_s != bbp_s_max) {
                    g->ubbp_s = calc_ubbp_s(g, 0, 0, 0, 0, 0);
                } else {
                    g->ubbp_s = 0.0;
                }
            } else {
                g->bbp_s = BAD_FLT;
                iopf[ip] |= IOPF_BADRRS;
                status = 1;
            }
            break;
        case BBPSMM01:
            // update power-law exponent based on Morel chl relationship
            if (l2rec->chl[ip] > 0.0) {
                g->bbp_s = MAX(MIN(0.5 * (0.3 - log10(g->chl)), bbp_s_max), bbp_s_min);
            if (g->bbp_s != bbp_s_min || g->bbp_s != bbp_s_max) {
                    g->ubbp_s = calc_ubbp_s(g, 0, 0, 0, 0, 0);
                } else {
                    g->ubbp_s = 0.0;
                }
            } else {
                g->bbp_s = BAD_FLT;
                iopf[ip] |= IOPF_BADRRS;
                status = 1;
            }
            break;
        case BBPSLAS:
            // update power-law exponent based on Loisel & Stramski model
            if ((g->bbp_s = get_bbp_las_eta(l2rec, ip)) == BAD_FLT) {
                iopf[ip] |= IOPF_BADRRS;
                status = 1;
            }
            g->ubbp_s = 0.0;
            printf("GIOP message: bbp slope uncertainties not computed during interface to LAS model.\n");
            printf("              ubbp_s uncertainty values set to: %f.\n", g->ubbp_s);
            g->bbp_nvec = 1;
            break;
        case BBPLAS:
        case BBPLASFIX:
            // replaces default tabbed values with results from Loisel & Stramski model
            if (get_bbp_las(l2rec, ip, g->bbp_tab_w, g->bbp_tab_s[0], g->bbp_tab_nw) == 0) {
                iopf[ip] |= IOPF_BADRRS;
                status = 1;
            }
            g->ubbp_s = 0.0;
            printf("GIOP message: bbp slope uncertainties not computed during interface to LAS model.\n");
            printf("              ubbp_s uncertainty values set to: %f.\n", g->ubbp_s);
            g->bbp_nvec = 1;
            break;
        case BBPQAAFIX:
            // replaces default tabbed values with results from QAA model
            if (get_bbp_qaa(l2rec, ip, g->bbp_tab_w, g->bbp_tab_s[0], g->bbp_tab_nw) == 0) {
                iopf[ip] |= IOPF_BADRRS;
                status = 1;
            }
            g->ubbp_s = 0.0;
            printf("GIOP message: bbp slope uncertainties not computed during interface to QAA model.\n");
            printf("              ubbp_s uncertainty values set to: %f.\n", g->ubbp_s);
            g->bbp_nvec = 1;
            break;
        case BBPLH:
            if (get_bbp_lh(l2rec, g, ipb2, g->bbp_tab_w, g->bbp_tab_s[0], g->ubbp_tab_s[0], g->bbp_tab_nw) == 0) {
                iopf[ip] |= IOPF_BADRRS;
                status = 1;
            }
            g->bbp_nvec = 1;
            break;
        case BBPCHL:
            if (get_bbp_chl(l2rec, g, ipb2, g->bbp_tab_w, g->bbp_tab_s[0], g->ubbp_tab_s[0], g->bbp_tab_nw) == 0) {
                iopf[ip] |= IOPF_BADRRS;
                status = 1;
            }
            g->bbp_nvec = 1;
            break;
        }

        // save updated model config, flag and skip if unable to compute

        aph_s[ip] = g->aph_s;
        adg_s[ip] = g->adg_s;
        bbp_s[ip] = g->bbp_s;
        uadg_s[ip] = g->uadg_s;
        ubbp_s[ip] = g->ubbp_s;

        if (status != 0) {
            iopf[ip] |= IOPF_BADRRS;
            continue;
        }

        // convert to subsurface reflectance

        for (iw = 0; iw < g->nwave; iw++) {
            ib = g->bindx[iw];
            Rrs_a[iw] = l2rec->Rrs[ipb2 + ib] - l2rec->Rrs_raman[ipb2 + ib];
            Rrs_b[iw] = rrs_above_to_below(Rrs_a[iw]);
            uRrs_b[iw] = rrs_above_to_below_unc(Rrs_a[iw],g->uRrs_a[iw]);
            
            wts [iw] = g->wts[iw]; // /rrs_above_to_below(Rrs_a[iw])*sqrt(l2rec->nobs[ip]); 
            //wts[iw] = 1.0 / pow(uRrs_b[iw], 2);  // testing with fitting weights
            if (Rrs_b[iw] > 0.0) {
                bndcnt++;
            }
        }

        // if less than npar valid reflectances, flag and skip pixel */

        if (bndcnt < g->npar) {
            iopf[ip] |= IOPF_BADRRS;
            continue;
        }

        // initialize model parameters

        giop_ctl_start(g, g->chl);

        for (ipar = 0; ipar < g->npar; ipar++) {
            par[ipar] = g->par[ipar]; // model params
            par[ipar + g->npar] = 0.0; // model param errors
            upar[ipar] = 0.0; // uncertainty model param errors
        }

        // run model optimization for this pixel

        switch (g->fit_opt) {
        case AMOEBA:
            status = fit_giop_amb(g, Rrs_b, wts, par, Rrs_f, &itercnt);
            break;
        case LEVMARQ:
            status = fit_giop_lm(g, Rrs_b, wts, par, &chi, &itercnt);
            break;
        case SVDFIT:
            status = fit_giop_svd(g, Rrs_b, wts, par);
            break;
        case SVDSIOP:
            status = fit_giop_svd_siop(g, Rrs_b, wts, par, &chi);
            /*If an optimal svd_siop solution not reached*/
            if (status == -99) {
                siop_num[ip] = -1;
            }
            break;
        default:
            printf("%s Line %d: Unknown optimization method for GIOP %d\n",
                    __FILE__, __LINE__, g->fit_opt);
            exit(1);
            break;
        }

        //Compute fit parameter uncertainties
        
        if (status == 0) {
            calc_par_unc(g, par, uRrs_b, upar);
        }
         
        //if no solution, flag as IOP prod failure

        if (status != 0) {
            iopf[ip] |= IOPF_FAILED;
            if (itercnt >= g->maxiter)
                iopf[ip] |= IOPF_MAXITER;
            continue;
        }

        //save parameter uncertainties to par array
        
        for (ipar = 0; ipar < g->npar; ipar++) {
            if (isfinite(par[ipar + g->npar]))
                //estimate combined model misfit uncertainty + radiometric uncertainty
                par[ipar+g->npar] = sqrt(pow(upar[ipar],2) +pow(par[ipar+g->npar],2));
        }

        // save final params 
        
        for (ipar = 0; ipar < g->npar; ipar++) {
            if (isfinite(par[ipar]))
                fit_par[ip][ipar] = par[ipar];
        }

        chisqr[ip] = chi;

        // evaluate model at fitted bands

        switch (g->fit_opt) {
        case SVDSIOP:
            giop_model_iterate(g, par, g->nwave, g->wave, g->aw, g->bbw, g->foq, aph1, adg1, bbp1, acdom1, anap1, bbph1, bbnap1, Rrs_f, NULL, NULL);
            break;
        default:
            giop_model(g, par, g->nwave, g->wave, g->aw, g->bbw, g->foq, aph1, adg1, bbp1, Rrs_f, NULL, NULL);
            break;
        }


        // bogus evaluation, flag and skip

        for (iw = 0; iw < g->nwave; iw++) {
            if (!isfinite(Rrs_f[iw])) {
                iopf[ip] |= IOPF_NAN;
                break;
            }
        }
        if (iopf[ip] & IOPF_NAN) {
            iopf[ip] |= IOPF_FAILED;
            continue;
        }

        // check goodness of fit

        rrsdiff[ip] = 0.0;
        for (iw = 0; iw < g->nwave; iw++) {
            if (g->wave[iw] >= 400 && g->wave[iw] <= 600) {
                if (fabs(Rrs_b[iw]) > 1e-7 && fabs(Rrs_f[iw] - Rrs_b[iw]) > 1e-5)
                    rrsdiff[ip] += fabs(Rrs_f[iw] - Rrs_b[iw]) / fabs(Rrs_b[iw]);
            }
        }
        rrsdiff[ip] /= g->nwave;
        if (rrsdiff[ip] > input->giop_rrs_diff) {
            iopf[ip] |= IOPF_RRSDIFF;
        }

        //if (status == 0) {
        //    calc_par_unc(g, par, rrs_diff, upar);
        //}

        // store in static globals

        switch (g->fit_opt) {
        case SVDSIOP:
            giop_model_iterate(g, par, nwave, wave, aw, bbw, foq, aph1, adg1, bbp1, acdom1, anap1, bbph1, bbnap1, Rrs_f, NULL, NULL);

            mRrs[ipb + ib] = rrs_below_to_above(Rrs_f[ib]) + l2rec->Rrs_raman[ipb2 + ib];

            for (ib = 0; ib < nwave; ib++) {

                if (isfinite(aph1[ib])) {
                    aph[ipb + ib] = aph1[ib];
                    aph[ierr + ib] = aph1[ib + nwave];
                } else {
                    iopf[ip] |= IOPF_NAN;
                }
                if (isfinite(acdom1[ib])) {
                    acdom[ipb + ib] = acdom1[ib];
                    acdom[ierr + ib] = acdom1[ib + nwave];
                } else {
                    iopf[ip] |= IOPF_NAN;
                }
                if (isfinite(anap1[ib])) {
                    anap[ipb + ib] = anap1[ib];
                    anap[ierr + ib] = anap1[ib + nwave];
                } else {
                    iopf[ip] |= IOPF_NAN;
                }
                if (isfinite(aph1[ib]) && isfinite(acdom1[ib]) && isfinite(anap1[ib])) {
                    adg[ipb + ib] = acdom1[ib] + anap1[ib];
                    adg[ierr + ib] = pow(pow(adg1[ib + nwave],2) + pow(acdom1[ib + nwave],2),0.5);
                    a[ipb + ib] = aw[ib] + aph1[ib] + acdom1[ib] + anap1[ib];
                    a[ierr + ib] = pow( pow(aph1[ib + nwave],2) + pow(adg1[ib + nwave],2) + pow(acdom1[ib + nwave],2) ,0.5 );
                } else {
                    iopf[ip] |= IOPF_NAN;
                }
                if (isfinite(bbnap1[ib])) {
                    bbnap [ipb + ib] = bbnap1[ib];
                    bbnap [ierr + ib] = bbnap1[ib + nwave];
                } else {
                    iopf[ip] |= IOPF_NAN;
                }
                if (isfinite(bbph1[ib])) {
                    bbph [ipb + ib] = bbph1[ib];
                    bbph [ierr + ib] = bbph1[ib + nwave];
                } else {
                    iopf[ip] |= IOPF_NAN;
                }
                if (isfinite(bbph1[ib]) && isfinite(bbnap1[ib])) {
                    bbp[ipb + ib] = bbph1[ib] + bbnap1[ib];
                    bbp[ierr + ib] = pow(pow(bbph1[ib + nwave],2) + pow(bbnap1[ib + nwave],2),0.5);
                    bb [ipb + ib] = bbw[ib] + bbph1[ib] + bbnap1[ib];
                    bb [ierr + ib] = pow(pow(bbph1[ib + nwave],2) + pow(bbnap[ib + nwave],2),0.5);
                } else {
                    iopf[ip] |= IOPF_NAN;
                }
            }

            /* check IOP ranges and set flags */
            set_iop_flag(wave, nwave, &a[ipb], &aph[ipb], &adg[ipb],
                    &bb[ipb], &bbp[ipb], &iopf[ip]);

            /*Optimal SIOP combination (for iterative aLMI solution)*/
            siop_num[ip] = 1 + g->siopIdx;

            /* aLMI compute chlorophyll */
            chl [ip] = par[0];

            break;

        default:
            //Fit model to full visible bandset
            giop_model(g, par, nwave, wave, aw, bbw, foq, aph1, adg1, bbp1, Rrs_f, NULL, NULL);

            for (ib = 0; ib < nwave; ib++) {

                mRrs[ipb + ib] = rrs_below_to_above(Rrs_f[ib]) + l2rec->Rrs_raman[ipb2 + ib];

                if (isfinite(aph1[ib])) {
                    aph[ipb + ib] = aph1[ib];
                    aph[ierr + ib] = aph1[ib + nwave];
                } else {
                    iopf[ip] |= IOPF_NAN;
                }
                if (isfinite(adg1[ib])) {
                    adg[ipb + ib] = adg1[ib];
                    adg[ierr + ib] = adg1[ib + nwave];
                } else {
                    iopf[ip] |= IOPF_NAN;
                }
                if (isfinite(aph1[ib]) && isfinite(adg1[ib])) {
                    a[ipb + ib] = aw[ib] + aph1[ib] + adg1[ib];
                    a[ierr + ib] = pow( pow(aph1[ib + nwave],2) + pow(adg1[ib + nwave],2), 0.5 );
                    a_unc[ipb + ib] = pow( pow(aph1[ib + nwave],2) + pow(adg1[ib + nwave],2), 0.5 );

                } else {
                    iopf[ip] |= IOPF_NAN;
                }
                if (isfinite(bbp1[ib])) {
                    bbp[ipb + ib] = bbp1[ib];
                    bbp[ierr + ib] = bbp1[ib + nwave];
                    bb [ipb + ib] = bbw[ib] + bbp1[ib];
                    bb [ierr + ib] = bbp1[ib + nwave];
                    bb_unc[ipb + ib] =  bbp1[ib + nwave];
                } else {
                    iopf[ip] |= IOPF_NAN;
                }
            }
            // check IOP ranges and set flags

            set_iop_flag(wave, nwave, &a[ipb], &aph[ipb], &adg[ipb],
                    &bb[ipb], &bbp[ipb], &iopf[ip]);

            // compute chlorophyll

            chl [ip] = giop_chl(g, iopf[ip], par, &chl[ip + l1rec->npix]);
            break;
        }


        iter[ip] = itercnt;
    }



    // fail pixels where any flags were set

    for (ip = 0; ip < l1rec->npix; ip++)
        if (iopf[ip] != 0) l1rec->flags[ip] |= PRODFAIL;


    LastRecNum = l1rec->iscan;
    free(Rrs_a);
    free(Rrs_b);
    free(Rrs_f);
    free(wts);
    free(par);
    free(upar);
    free(uRrs_a);
    free(uRrs_b);
}


/*----------------------------------------------------------------------*/
/* extract_band_3d() - utility function to extract selected band subset */

/*----------------------------------------------------------------------*/
static void extract_band_3d(float *in_buf, float *out_buf, int numPixels, int numBands) {
    float * out_ptr = out_buf;
    for (int pix = 0; pix < numPixels; pix++) {
        float *in_ptr = in_buf + pix * numBands;
        for (int band_3d = 0; band_3d < input->nwavelengths_3d; band_3d++) {
            int band = input->wavelength_3d_index[band_3d];
            *out_ptr = in_ptr[band];
            out_ptr++; 
        }
    }
}

/*----------------------------------------------------------------------*/
/* extract_band_3d_unc() - utility function to extract selected band    */
/*                         subset of spectral product uncertainties     */

/*----------------------------------------------------------------------*/
static void extract_band_3d_unc(float *in_buf, float *out_buf, int numPixels, int numBands) {
    float * out_ptr = out_buf;
    for (int pix = 0; pix < numPixels; pix++) {
        float *in_ptr = in_buf + (pix * numBands + numPixels*numBands) ;
        for (int band_3d = 0; band_3d < input->nwavelengths_3d; band_3d++) {
            int band = input->wavelength_3d_index[band_3d];
            *out_ptr = in_ptr[band];
            out_ptr++; 
        }
    }
}

/* ------------------------------------------------------------------- */
/* get_giop() - returns requested GIOP product for full scanline       */

/* ------------------------------------------------------------------- */
void get_giop(l2str *l2rec, l2prodstr *p, float prod[]) {
    int prodID = p->cat_ix;
    int ib = p->prod_ix;

    int32_t ip, ipb, ierr;
    l1str *l1rec = l2rec->l1rec;

    /*Before running GIOP, check if valid output products selected*/
    /*Note: aCDOM, aNAP, bbPh, bbNAP are only valid choices for giop_fit_opt SVDSIOP*/
    switch (input->giop_fit_opt) {
    case SVDSIOP:
        break;
    default:
        switch (prodID) {
        case CAT_acdom_giop:
        case CAT_anap_giop:
        case CAT_bbph_giop:
        case CAT_bbnap_giop:
        case CAT_acdom_unc_giop:
        case CAT_anap_unc_giop:
        case CAT_bbph_unc_giop:
        case CAT_bbnap_unc_giop:
            printf("-E- %s line %d : products acdom, anap, bbph and bbnap are only applicable with giop_fit_opt=SVDSIOP.\n",
                    __FILE__, __LINE__);
            exit(1);
            break;
        }
    }

    if (!giop_ran(l1rec->iscan))
        run_giop(l2rec);

    if(p->rank == 3) {
        switch (prodID) {

        case CAT_aph_giop:
            extract_band_3d(aph, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;

        case CAT_adg_giop:
            extract_band_3d(adg, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;

        case CAT_acdom_giop:
            extract_band_3d(acdom, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;

        case CAT_anap_giop:
            extract_band_3d(anap, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;

        case CAT_bbp_giop:
            extract_band_3d(bbp, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;

        case CAT_a_giop:
            extract_band_3d(a, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;

        case CAT_bb_giop:
            extract_band_3d(bb, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;

        case CAT_bbph_giop:
            extract_band_3d(bbph, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;

        case CAT_bbnap_giop:
            extract_band_3d(bbph, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;

        case CAT_a_unc_giop:
            extract_band_3d_unc(a, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;

        case CAT_aph_unc_giop:
            extract_band_3d_unc(aph, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;

        case CAT_adg_unc_giop:
            extract_band_3d_unc(adg, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;
            
        case CAT_acdom_unc_giop:
            extract_band_3d_unc(acdom, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;
        
        case CAT_anap_unc_giop:
            extract_band_3d_unc(anap, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;
        
        case CAT_bb_unc_giop:
            extract_band_3d_unc(bb, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;

        case CAT_bbp_unc_giop:
            extract_band_3d_unc(bbp, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;
        
        case CAT_bbph_unc_giop:
            extract_band_3d_unc(bbph, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;
        
        case CAT_bbnap_unc_giop:
            extract_band_3d_unc(bbnap, prod, l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
            break;         

        default:
            printf("-E- %s line %d : product ID %d passed to GIOP does not support 3D data sets.\n",
                    __FILE__, __LINE__, prodID);
            exit(1);
        } // switch prodID

    } else {

        for (ip = 0; ip < l1rec->npix; ip++) {

            ipb = ip * l1rec->l1file->nbands + ib;
            ierr = l1rec->npix * l1rec->l1file->nbands + ipb;

            switch (prodID) {

            case CAT_mRrs_giop:
                prod[ip] = (float) mRrs[ipb];
                break;

            case CAT_aph_giop:
                prod[ip] = (float) aph[ipb];
                break;

            case CAT_adg_giop:
                prod[ip] = (float) adg[ipb];
                break;

            case CAT_bbp_giop:
                prod[ip] = (float) bbp[ipb];
                break;

            case CAT_a_giop:
                prod[ip] = (float) a[ipb];
                break;

            case CAT_bb_giop:
                prod[ip] = (float) bb[ipb];
                break;

            case CAT_acdom_giop:
                prod[ip] = (float) acdom[ipb];
                break;

            case CAT_anap_giop:
                prod[ip] = (float) anap[ipb];
                break;

            case CAT_bbph_giop:
                prod[ip] = (float) bbph[ipb];
                break;

            case CAT_bbnap_giop:
                prod[ip] = (float) bbnap[ipb];
                break;

            case CAT_chl_giop:
                prod[ip] = (float) chl[ip];
                break;

            case CAT_opt_siop_giop:
                prod[ip] = (int) siop_num[ip];
                break;

            case CAT_aph_unc_giop:
                prod[ip] = (float) aph[ierr];
                break;

            case CAT_adg_unc_giop:
                prod[ip] = (float) adg[ierr];
                break;

            case CAT_acdom_unc_giop:
                prod[ip] = (float) adg[ierr];
                break;

            case CAT_anap_unc_giop:
                prod[ip] = (float) adg[ierr];
                break;

            case CAT_bbp_unc_giop:
                prod[ip] = (float) bbp[ierr];
                break;

            case CAT_bbph_unc_giop:
                prod[ip] = (float) bbph[ierr];
                break;

            case CAT_bbnap_unc_giop:
                prod[ip] = (float) bbnap[ierr];
                break;

            case CAT_a_unc_giop:
                prod[ip] = (float) a[ierr];
                break;

            case CAT_bb_unc_giop:
                prod[ip] = (float) bb[ierr];
                break;

            case CAT_chl_unc_giop:
                prod[ip] = (float) chl[ip + l1rec->npix];
                break;

            case CAT_aphs_giop:
                prod[ip] = (float) aph_s[ip];
                break;

            case CAT_adgs_giop:
                prod[ip] = (float) adg_s[ip];
                break;
            
            case CAT_adgs_unc_giop:
                prod[ip] = (float) uadg_s[ip];
                break;

            case CAT_bbps_giop:
                prod[ip] = (float) bbp_s[ip];
                break;

            case CAT_bbps_unc_giop:
                prod[ip] = (float) ubbp_s[ip];
                break;

            case CAT_rrsdiff_giop:
                prod[ip] = (float) rrsdiff[ip];
                break;

            case CAT_chisqr_giop:
                prod[ip] = (float) chisqr[ip];
                break;

            case CAT_fitpar_giop:
                if (ib >= max_npar) {
                    printf("-E- %s line %d : output request for GIOP fit parameter %d exceeds number of fit parameters %d.\n",
                            __FILE__, __LINE__, ib, max_npar);
                    exit(1);
                }
                prod[ip] = (float) fit_par[ip][ib];
                break;

            default:
                printf("-E- %s line %d : erroneous product ID %d passed to GIOP.\n",
                        __FILE__, __LINE__, prodID);
                exit(1);
            }
        } // for ip
    } // if rank != 3

    return;
}


/* ------------------------------------------------------------------- */
/* get_iter_giop() - returns iteration count                           */

/* ------------------------------------------------------------------- */
int16 *get_iter_giop(l2str *l2rec) {
    if (!giop_ran(l2rec->l1rec->iscan))
        run_giop(l2rec);

    return iter;
}


/* ------------------------------------------------------------------- */
/* get_flags_giop() - returns iteration count                          */

/* ------------------------------------------------------------------- */
int16 *get_flags_giop(l2str *l2rec) {
    if (!giop_ran(l2rec->l1rec->iscan))
        run_giop(l2rec);

    return iopf;
}


/* ------------------------------------------------------------------- */
/* Interface to convl12() to return GIOP iops                          */

/* ------------------------------------------------------------------- */
void iops_giop(l2str *l2rec) {
    int32_t ib, ip, ipb, ipb2;

    int32_t nbands = l2rec->l1rec->l1file->nbands;
    int32_t npix = l2rec->l1rec->npix;

    if (!giop_ran(l2rec->l1rec->iscan))
        run_giop(l2rec);

    for (ip = 0; ip < npix; ip++) {
        for (ib = 0; ib < nbands; ib++) {
            ipb2 = ip * nbands + ib;
            ipb = ip * nbands + ib;
            l2rec->a [ipb2] = (float) a[ipb];
            l2rec->bb[ipb2] = (float) bb[ipb];
        }
    }

    return;
}


// PJW 9 Jan 2013
/* ------------------------------------------------------------------- */
/* give external access to local *chl                          */

/* ------------------------------------------------------------------- */
float* giop_get_chl_pointer() {
    return chl;
}
/* ------------------------------------------------------------------- */
/* give external access to local *adg                          */

/* ------------------------------------------------------------------- */
float* giop_get_adg_pointer() {
    return adg;
}
/* ------------------------------------------------------------------- */
/* give external access to local *bbp                          */

/* ------------------------------------------------------------------- */
float* giop_get_bbp_pointer() {
    return bbp;
}
/* ------------------------------------------------------------------- */
/* give external access to local *aph                          */

/* ------------------------------------------------------------------- */
float* giop_get_aph_pointer() {
    return aph;
}
/* ------------------------------------------------------------------- */
/* give external access to local **fit_par                         */

/* ------------------------------------------------------------------- */
float** giop_get_fitpar_pointer() {
    return fit_par;
}

/* ------------------------------------------------------------------- */
/* give external access to local *bbp_s                          */

/* ------------------------------------------------------------------- */
float* giop_get_bbp_s_pointer() {
    return bbp_s;
}

/* ------------------------------------------------------------------- */
float* giop_get_ubbp_s_pointer() {
    return ubbp_s;
}

/* ------------------------------------------------------------------- */
float* giop_get_a_unc_pointer() {
    return a_unc;
}

/* ------------------------------------------------------------------- */
float* giop_get_bb_unc_pointer() {
    return bb_unc;
}