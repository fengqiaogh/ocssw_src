#ifndef tjsm_pml_iop_tables 
#define tjsm_pml_iop_tables
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE 180

/* These are declared here as an extern so they will be useable */
/* by the rest of the code (static was previously used incorrectly)*/
extern float *lambda, lc[MAX_BANDS], *a_w, *b_w;

/* Geophysical (GOP) variables */
extern int32_t nband, ch_n, sp_n, od_n;
extern float *ch_lev, *ac[MAX_BANDS], *bc[MAX_BANDS];
extern float *sp_lev, *as[MAX_BANDS], *bs[MAX_BANDS];
extern float *od_lev, *od[MAX_BANDS];

/* IOP variables */
extern int32_t th_s_n, th_v_n, dphi_n;
extern float *th_s_lev, *th_v_lev, *dphi_lev;
extern int32_t ap_n, bp_n;
extern float *ap_lev, *bp_lev;

extern float *refen; /* The pointer for refen.*/

/* 490:510 ratio - step 1 of the model */
extern int bp[2], maxit;
extern float eps_a_init, eps_a_init_modis, init_chl, tol;
extern float b_tilde_w, b_tilde_p;
extern float scat_a, scat_b, scat_c, scat_n, scat_l, scat_l_modis;

/* Gelbstoff and pigment parameters */
extern float eps_y_412_443, eps_p_412_443;
extern float ysbpa_0, ysbpa_s, ysbpa_l;

/* Bright pixel externals */
extern float tol_n, tol_b, n_init;
extern int bp_base, bp_1, bp_2, max_iter, n_bands;
extern double lc1, lc2, log_lc1, e_init;
extern double b_low_init, b_high_init, b_init, delta_b_init, min_db, max_db;
extern double spm_max, spm_min, rst_spm, rst_db, rst_n, rst_n_step;
extern double n_min, n_max, iter_scale, clim_spm;

/* Function declarations */
int load_work_tab(char *configname, int sensorID);
void load_config(char *configname);
float geo2iop(float *levels, float *iopv[MAX_BANDS], int band, float value, int size);
float interp_l(float *x, float u, int n);
float interp(float *x, float u, int n);
int setgeom(float sun_theta, float sen_theta, float dphi);
double f_ab(double a, double b, int band);
double fint(double a, double b, int band);

/* These functions are used in the bright pixel code */
float r_ab(float a, float b, int band);
float sed_ref(float spm, int band);
float chl_ref(float chl, int band);
float iop_ref(float conc, int band, int iop);

#endif
