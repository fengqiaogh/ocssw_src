#include <math.h>
#include "pml_iop.h"
#include "pml_iop_config.h"
#include "pml_iop_tables.h"
#include "pml_iop_calculate.h"
#include "pml_bright.h"
#include "l12_parms.h"
#include "genutils.h"

/* Module: pml_iop_tables.c */
/* Authors: Gerald Moore and Tim Smyth */
/* Date: 06/03/13 */
/* Version: 1.5 */
/* Description: Module containing functions to read in the various */
/* tables defined in the config file and to perform interpolations */

#define MAX_LINE 180

/* set the endian-ness of the PML iop table binary file (1=little endian) */
#define PML_FILE_ENDIAN 1


float *lambda, lc[MAX_BANDS], *a_w, *b_w;
/* Geophysical (GOP) variables */
int32_t nband, ch_n, sp_n, od_n;
float *ch_lev, *ac[MAX_BANDS], *bc[MAX_BANDS];
float *sp_lev, *as[MAX_BANDS], *bs[MAX_BANDS];
float *od_lev, *od[MAX_BANDS];

/* IOP variables */
int32_t th_s_n, th_v_n, dphi_n;
float *th_s_lev, *th_v_lev, *dphi_lev;
int32_t ap_n, bp_n;
float *ap_lev, *bp_lev;

float *refen; /* The pointer for refen.*/
static int32_t refind = 0l; /* current index */
static int32_t refind_max = 0l;

/* Constants for interior conversions */
/* loaded by config routine */
int bp[2], maxit;
float b_tilde_w, b_tilde_p, init_chl, eps_a_init, eps_a_init_modis, eps_bb_init;
float scat_a, scat_b, scat_c, scat_n, scat_l, scat_l_modis;
float tol;

/* Gelbstoff and pigment parameters */
float eps_y_412_443, eps_p_412_443;
float ysbpa_0, ysbpa_s, ysbpa_l;

/* Bright pixel variables */
int bp_base, bp_1, bp_2, max_iter, n_bands;
float tol_n, tol_b, n_init;
double lc1, lc2, log_lc1, e_init;
double b_low_init, b_high_init, b_init, delta_b_init, min_db, max_db;
double spm_max, spm_min, rst_spm, rst_db, rst_n, rst_n_step;
double n_min, n_max, iter_scale, clim_spm;

/* Function: load_work_tab */

/* Loads the sensor Look-Up Tables into memory  */
int load_work_tab(char *configfname, int sensorID) {
    char *tmp_str;
    char fname[FILENAME_MAX];

    int i, j;
    int32_t t_nbands, v_len;
    float *t_lambda;
    char created[MAX_LINE];
    char sensor_name[MAX_LINE];

    FILE *table;

    /* Geophysical LUT */
    if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
        printf("OCDATAROOT environment variable is not defined.\n");
        exit(1);
    }

    /* sensor specific LUTs used - only available for MODIS and SeaWiFS */
    switch (sensorID) {
    case SEAWIFS:
        sprintf(sensor_name, "seawifs");
        t_nbands = 8;
        break;
    case MODISA:
    case MODIST:
        sprintf(sensor_name, "modis");
        t_nbands = 9;
        break;
    default:
        printf("Sensor not either SeaWiFS or MODIS and the LUTs do not exist\n");
        exit(1);
        break;
    }

    sprintf(fname, "%s/%s/iop/pml/%s", tmp_str, sensor_name, get_cfg_s("gop_table", configfname));
    if ((table = fopen(fname, "r")) == NULL) {
        printf("Error opening %s\n", fname);
        exit(1);
    }

    /* start by reading bands and alocating variables */
    /* t_nbands=(int32_t)get_cfg_i("n_bands",configfname); */
    fread_swap(PML_FILE_ENDIAN, &nband, sizeof (int32_t), 1, table);
    /* verify correct size of table */
    if (nband != t_nbands) {
        printf("Table band mismatch expected %d read %d\n", t_nbands, nband);
        printf("Using Geophysical Look-up Table: %s \n", fname);
        exit(1);
    }
    /* Allocate the arrays according to the number of wavelengths */
    lambda = calloc(nband, sizeof (float));
    a_w = calloc(nband, sizeof (float));
    b_w = calloc(nband, sizeof (float));
    /* Read in the wavelengths and aw,bw (in that order) */
    fread_swap(PML_FILE_ENDIAN, lambda, sizeof (float), nband, table);
    fread_swap(PML_FILE_ENDIAN, a_w, sizeof (float), nband, table);
    fread_swap(PML_FILE_ENDIAN, b_w, sizeof (float), nband, table);

    /* Chlorophyll */
    fread_swap(PML_FILE_ENDIAN, &ch_n, sizeof (int32_t), 1, table);
    ch_lev = calloc(ch_n, sizeof (float));
    fread_swap(PML_FILE_ENDIAN, ch_lev, sizeof (float), ch_n, table);

    for (i = 0; i < nband; i++) {
        ac[i] = calloc(ch_n, sizeof (float));
        bc[i] = calloc(ch_n, sizeof (float));
    }
    for (i = 0; i < ch_n; i++) {
        for (j = 0; j < nband; j++) {
            fread_swap(PML_FILE_ENDIAN, &ac[j][i], sizeof (float), 1, table);
        }
    }
    for (i = 0; i < ch_n; i++) {
        for (j = 0; j < nband; j++) {
            fread_swap(PML_FILE_ENDIAN, &bc[j][i], sizeof (float), 1, table);
        }
    }

    /* SPM */
    fread_swap(PML_FILE_ENDIAN, &sp_n, sizeof (int32_t), 1, table);
    sp_lev = calloc(sp_n, sizeof (float));
    fread_swap(PML_FILE_ENDIAN, sp_lev, sizeof (float), sp_n, table);

    for (i = 0; i < nband; i++) {
        as[i] = calloc(sp_n, sizeof (float));
        bs[i] = calloc(sp_n, sizeof (float));
    }
    for (i = 0; i < sp_n; i++) {
        for (j = 0; j < nband; j++) {
            fread_swap(PML_FILE_ENDIAN, &as[j][i], sizeof (float), 1, table);
        }
    }
    for (i = 0; i < sp_n; i++) {
        for (j = 0; j < nband; j++) {
            fread_swap(PML_FILE_ENDIAN, &bs[j][i], sizeof (float), 1, table);
        }
    }

    /* Gelbstoff */
    fread_swap(PML_FILE_ENDIAN, &od_n, sizeof (int32_t), 1, table);
    od_lev = calloc(od_n, sizeof (float));
    fread_swap(PML_FILE_ENDIAN, od_lev, sizeof (float), od_n, table);

    for (i = 0; i < nband; i++) {
        od[i] = calloc(od_n, sizeof (float));
    }
    for (i = 0; i < od_n; i++) {
        for (j = 0; j < nband; j++) {
            fread_swap(PML_FILE_ENDIAN, &od[j][i], sizeof (float), 1, table);
        }
    }

    fread_swap(PML_FILE_ENDIAN, &v_len, sizeof (int32_t), 1, table);
    fread_swap(PML_FILE_ENDIAN, created, sizeof (char), v_len, table);
    created[v_len] = '\0';

    fclose(table);

    /* IOP tables */
    /* IOP header file */
    sprintf(fname, "%s/%s/iop/pml/%s", tmp_str, sensor_name, get_cfg_s("iop_F_head", configfname));
    if ((table = fopen(fname, "r")) == NULL) {
        printf("Error opening %s\n", fname);
        exit(1);
    }

    /* start by reading bands and alocating variables */
    /*    t_nbands=(int32_t)get_cfg_i("n_bands",configfname); */
    fread_swap(PML_FILE_ENDIAN, &nband, sizeof (int32_t), 1, table);
    /* verify correct size of table */
    if (nband != t_nbands) {
        printf("Table band mismatch expected %d read %d\n", t_nbands, nband);
        printf("Using IOP  header Table: %s \n", fname);
        exit(1);
    }
    /* First wavelengths */
    t_lambda = calloc(nband, sizeof (float));
    fread_swap(PML_FILE_ENDIAN, t_lambda, sizeof (float), nband, table);

    /* verify that these match */
    for (i = 0; i < nband; i++) {
        if (fabsf(t_lambda[i] - lambda[i]) > 15.0) {
            printf("Error IOP table wavelength mismatch expected %f got %f\n", lambda[i], t_lambda[i]);
            exit(1);
        }
    }
    free(t_lambda);
    /* Solar zenith angle */
    fread_swap(PML_FILE_ENDIAN, &th_s_n, sizeof (int32_t), 1, table);
    th_s_lev = calloc(th_s_n, sizeof (float));
    fread_swap(PML_FILE_ENDIAN, th_s_lev, sizeof (float), th_s_n, table);
    /* view zenith angle */
    fread_swap(PML_FILE_ENDIAN, &th_v_n, sizeof (int32_t), 1, table);
    th_v_lev = calloc(th_v_n, sizeof (float));
    fread_swap(PML_FILE_ENDIAN, th_v_lev, sizeof (float), th_v_n, table);
    /* azimuth angle difference (view and satellite)*/
    fread_swap(PML_FILE_ENDIAN, &dphi_n, sizeof (int32_t), 1, table);
    dphi_lev = calloc(dphi_n, sizeof (float));
    fread_swap(PML_FILE_ENDIAN, dphi_lev, sizeof (float), dphi_n, table);
    /* absorption (a) */
    fread_swap(PML_FILE_ENDIAN, &ap_n, sizeof (int32_t), 1, table);
    ap_lev = calloc(ap_n, sizeof (float));
    fread_swap(PML_FILE_ENDIAN, ap_lev, sizeof (float), ap_n, table);
    /* scatter (b) */
    fread_swap(PML_FILE_ENDIAN, &bp_n, sizeof (int32_t), 1, table);
    bp_lev = calloc(bp_n, sizeof (float));
    fread_swap(PML_FILE_ENDIAN, bp_lev, sizeof (float), bp_n, table);

    fclose(table);

    /* Master IOP table for the geometry */
    sprintf(fname, "%s/%s/iop/pml/%s", tmp_str, sensor_name, get_cfg_s("iop_F_table", configfname));
    if ((table = fopen(fname, "r")) == NULL) {
        printf("Error opening %s\n", fname);
        exit(1);
    }

    /* Allocate memory for refen array */
    refind = th_s_n * th_v_n * dphi_n * ap_n * bp_n*nband;
    refind_max = refind - (ap_n * bp_n * nband);
    refen = (float *) calloc((refind), sizeof (float));
    if (refen == NULL) {
        fprintf(stderr, "load_work_tab: memory allocation failure for refen\n");
        exit(1);
    }

    fread_swap(PML_FILE_ENDIAN, refen, sizeof (float), refind, table);

    /* Check that table limits are not exceeded */
    if (feof(table)) {
        printf("IOP master table error: table too short!\n");
        exit(1);
    }
    /* Try next byte as constitancy check */
    i = fgetc(table);
    if (!feof(table)) {
        printf("IOP master table error: table too long!\n");
        exit(1);
    }

    fclose(table);
    refind = 0l;
    return refind;
}

/* Function: load_config */
/* Stores the config table values such that they are available to */

/* the main calling routine */
void load_config(char *configfname) {
    bp[0] = get_cfg_i("low_band", configfname);
    bp[1] = get_cfg_i("high_band", configfname);
    b_tilde_p = get_cfg_f("b_tilde_p", configfname);
    b_tilde_w = get_cfg_f("b_tilde_w", configfname);

    /* Step 1: 490:510 ratio absorption and backscatter */
    eps_a_init = get_cfg_f("eps_a_init", configfname);
    eps_bb_init = get_cfg_f("eps_bb_init", configfname);
    /* Additional information if processing a MODIS pass */
    eps_a_init_modis = get_cfg_f("eps_a_init_modis", configfname);
    scat_l_modis = get_cfg_f("scat_l_modis", configfname);
    scat_a = get_cfg_f("scat_a", configfname);
    scat_b = get_cfg_f("scat_b", configfname);
    scat_c = get_cfg_f("scat_c", configfname);
    scat_n = get_cfg_f("scat_n", configfname);
    scat_l = get_cfg_f("scat_l", configfname);

    /* Iterations to get the F parameter correctly */
    init_chl = get_cfg_f("init_chl", configfname);
    tol = get_cfg_f("iop_tol", configfname);
    maxit = get_cfg_i("iop_maxit", configfname);

    /* Step 2: 412:443 ratios to get the gelbstoff and pigment */
    eps_y_412_443 = get_cfg_f("eps_y_412_443", configfname);
    eps_p_412_443 = get_cfg_f("eps_p_412_443", configfname);

    ysbpa_0 = get_cfg_f("YSBPA_0", configfname);
    ysbpa_s = get_cfg_f("YSBPA_S", configfname);
    ysbpa_l = get_cfg_f("YSBPA_l", configfname);

    /* Additional variables for the bright pixel code */
    bp_base = get_cfg_i("bp_base", configfname); /* base band */
    bp_1 = get_cfg_i("bp_1", configfname); /* first NIR */
    bp_2 = get_cfg_i("bp_2", configfname); /* second(lower) NIR */

    /* convergence criterion */
    tol_b = get_cfg_f("bp_tol_b", configfname); /* Tolerence for converge */
    tol_n = get_cfg_f("bp_tol_n", configfname); /* Tolerence for converge */
    max_iter = get_cfg_i("bp_iter", configfname); /* Max iterations */

    n_init = get_cfg_f("bp_n", configfname); /* inititial Angstrom exponent */
    n_min = get_cfg_f("n_min", configfname); /* min. Angstrom exponent */
    n_max = get_cfg_f("n_max", configfname); /* max. Angstrom exponent */
    b_low_init = get_cfg_f("low_spm_init", configfname);
    b_high_init = get_cfg_f("high_spm_init", configfname);
    b_init = get_cfg_f("b_init", configfname);

    delta_b_init = get_cfg_f("delta_b_init", configfname);
    min_db = get_cfg_f("min_db", configfname);
    max_db = get_cfg_f("max_db", configfname);
    spm_min = get_cfg_f("spm_min", configfname);
    spm_max = get_cfg_f("spm_max", configfname);
    iter_scale = get_cfg_f("iter_scale", configfname);
    rst_n = get_cfg_f("rst_n", configfname);
    rst_n_step = get_cfg_f("rst_n_step", configfname);
    rst_spm = get_cfg_f("rst_spm", configfname);
    rst_db = get_cfg_f("rst_db", configfname);
    clim_spm = get_cfg_f("clim_spm", configfname);

}

/* Function: geo2iop */

/* Returns the IOP values for a given geophysical input */
float geo2iop(float *levels, float *iopv[MAX_BANDS], int band, float value, int size) {

    float res, ind;
    /* Conversion of band number to account for an 11 band IOP table */
    ind = interp_l(levels, value, size);
    res = iopv[band][(int) ind]*(floor(ind) + 1.0 - ind) + iopv[band][(int) ind + 1]*(ind - floor(ind));

    return res;
}

/* Function: interp */

/* Returns an interpolated value */
float interp(float *x, float u, int n) {
    int s, i;
    s = 0;
    if (u > x[n - 1]) return n - 1;
    for (i = 0; i < n; i++) {
        if (x[i] >= u) {
            s = i - 1;
            break;
        }
    }
    if (s < 0) s = 0;
    return ((u - x[s]) / (x[s + 1] - x[s]) + s);
}

/* Function: interp_l */
/* Returns the log interpolated value */

/* Fast version for geophysical variables */
float interp_l(float *x, float u, int n) {
    int s, i;
    s = 0;
    if (u > x[n - 1]) return n - 1;
    for (i = 0; i < n; i++) {
        if (x[i] >= u) {
            s = i - 1;
            break;
        }
    }
    if (s < 0) s = 0;
    if (s == 0) return ((u - x[s]) / (x[s + 1] - x[s]) + s);
    return ((log(u) - log(x[s])) / (log(x[s + 1]) - log(x[s])) + s);
}

/* Function: setgeom */

/* Interpolates from the LUT to get correct value for geometry */
int setgeom(float sun_theta, float sen_theta, float dphi) {
    int status = 0;
    int32_t th_s_ent, th_v_ent, dphi_ent;

    /* Page in IOP table */
    /* check that the angle matches index */
    th_s_ent = (int32_t) floor(interp(th_s_lev, sun_theta, th_s_n) + 0.5);
    th_v_ent = (floor) (interp(th_v_lev, sen_theta, th_v_n) + 0.5);
    if (dphi < 0.0) dphi = dphi + OEL_PI * 2.0;
    if (dphi > OEL_PI * 2.0) dphi = dphi - OEL_PI * 2.0;
    dphi_ent = (int32_t) floor(interp(dphi_lev, dphi, dphi_n) + 0.5);
    /* IOP data */
    refind = (int) (th_s_ent * th_v_n * dphi_n * nband * bp_n * ap_n)
            +(th_v_ent * dphi_n * nband * bp_n * ap_n)
            +(dphi_ent * nband * bp_n * ap_n);

    /* Geometry beyond table limits */
    if (refind > refind_max) {
        status = 1;
        return status;
    }
    return status;
}

/* Function: f_ab */
/* Interpolates to get a value f/Q */

/* for a pair of a and b */
double f_ab(double a, double b, int band) {
    double ain, bin;
    double res;
    ain = interp_l(ap_lev, a, ap_n);
    bin = interp_l(bp_lev, b, bp_n);
    res = fint(ain, bin, band);

    return res;
}

/* Function: fint */

/* Interpolates for a pair of doubles and a given waveband */
double fint(double a, double b, int band) {
    double ral, rah, res;
    int al, bl, bh, ah;
    if (a > 0) al = (int) floor(a);
    else al = 0;
    if (b > 0) bl = (int) floor(b);
    else bl = 0;
    bh = bl + 1;
    ah = al + 1;

    if ((bh > 15) || (ah > 15)) return 0.0;
    ral = refen[refind + band * bp_n * ap_n + bl * ap_n + ah]*(a - (double) al)
            + refen[refind + band * bp_n * ap_n + bl * ap_n + al]*((double) ah - a);
    rah = refen[refind + band * bp_n * ap_n + bh * ap_n + ah]*(a - (double) al)
            + refen[refind + band * bp_n * ap_n + bh * ap_n + al]*((double) ah - a);
    res = rah * (b - (double) bl) + ral * ((double) bh - b);

    return res;
}

/* Function: iop_ref */

/* Depending on the iop binary flag returns sediment / chlorophyll conc */
float iop_ref(float conc, int band, int iop) {
    /* Sediment is set for 0, chlorophyll for 1 .. */
    if (iop == 0) return sed_ref(conc, band);
    if (iop == 1) return chl_ref(conc, band);
    printf("Call error in iop_ref (type = %d)\n", iop);
    exit(1);
}

/* Function: r_ab */

/* Given value of a and b return reflectance (from LUT) */
float r_ab(float a, float b, int band) {
    double F;
    float res;
    F = f_ab(a, b, band);
    res = F * (b_w[band] * b_tilde_w + b * b_tilde_p) / (a_w[band] + a);

    return res;
}

/* Function: sed_ref */

/* Returns reflectance for a given sediment value */
float sed_ref(float spm, int band) {
    float res;
    if (spm < spm_min) spm = spm_min;
    else if (spm > spm_max) spm = spm_max;
    res = r_ab(geo2iop(sp_lev, as, band, spm, sp_n), geo2iop(sp_lev, bs, band, spm, sp_n), band);

    return res;
}

/* Function: chl_ref */

/* Returns reflectance for a given chlorophyll value */
float chl_ref(float chl, int band) {
    float res;
    /* Need to debug this (or at least make nomenclature change */
    if (chl < spm_min) chl = spm_min;
    else if (chl > spm_max) chl = spm_max;
    res = r_ab(geo2iop(ch_lev, ac, band, chl, ch_n), geo2iop(ch_lev, bc, band, chl, ch_n), band);

    return res;
}

