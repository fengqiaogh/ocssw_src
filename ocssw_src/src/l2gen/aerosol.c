/* ======================================================================================== */
/* module aerosol.c  - functions to facilitate aerosol model selection and application      */
/*                                                                                          */
/* Description:                                                                             */
/*                                                                                          */
/* This code replaces the original set of fortran subroutines developed by M.Wang, H.Gordon,*/
/* and others (e.g., rho_a_sub_quad, linear_abc, funct_eps, load_aer, load_ss11) as well as */
/* the original get_aerosol() developed for MSl12.                                          */
/*                                                                                          */
/* The functions herein read and interpolate the aerosol model tables, which are now stored */
/* as individual HDF files per model.  Whfere sensor wavelengths differ from tabulated model */
/* wavelengths, interpolation is performed. Efficiencies are gained by remembering and      */
/* reusing computational results when applicable.                                           */
/*                                                                                          */
/* Primary Function:                                                                        */
/* ----------------                                                                         */
/* aerosol() - compute aerosol reflectance using specified algorithm (main interface)       */
/*                                                                                          */
/* Secondary Functions:                                                                     */
/* -------------------                                                                      */
/* wangaer() - compute aerosol reflectance using Gordon & Wang 1994 algorithm               */
/* fixedaot() - compute aerosol reflectance for fixed aot(lambda)                           */
/* fixedaer() - compute aerosol reflectance for fixed aerosol model                         */
/* get_angstrom() - compute angstrom coefficient (interface to l2_hdf_generic)              */
/* diff_tran() - compute Rayleigh-aerosol diffuse trans for selected model pair, both paths */
/*                                                                                          */
/* Supporting Functions:                                                                    */
/* --------------------                                                                     */
/* load_aermod() - loads the entire aerosol model table for the specified model list        */
/* ss_to_ms_coef() - for one model, return coefficients of function relating single         */
/*                   scattering to multiple scattering at the input geometry.               */
/* rhoas_to_rhoa() - SS aerosol reflectance to MS aerosol reflectance                       */
/* rhoa_to_rhoas() - MS aerosol reflectance to SS aerosol reflectance                       */
/* model_select_wang() - M. Wang aerosol model selection process                 .          */
/* model_select_angst() - Select model pair based on input Angstrom coefficient             */
/* model_phase() - return phase function at model wavelengths at the input geometry.        */
/* model_epsilon() - return model epsilon at input wavelengths for input geometry.          */
/* model_transmittance() - compute path Rayleigh-aerosol diffuse trans for specified model  */
/* model_taua() - compute AOT at input wavelengths using specified model                    */
/* aeroob - out-of-band water-vapor correction                                              */
/*                                                                                          */
/*                                                                                          */
/* Written By:  B. Franz, SAIC, NASA/GSFC Ocean Biology Processing Group, Summer 2004.      */
/* W. Robinson, SAIC, 24 MAr 2017, modularized and enhanced for band-dependent              */
/*  geometry                                                                                */
/*                                                                                          */
/* ======================================================================================== */

#include <float.h>

#include "l12_proto.h"

#include <allocate2d.h>
#include <allocate3d.h>
#include <allocate5d.h>

#define MAXMODEL MAXAERMOD
//#define MAXSOLZ  33
//#define MAXSENZ  35
//#define MAXPHI   19
//#define MAXSCATT 75
//#define DTNTHETA 33

static int32_t first_solz;
static int32_t first_senz;
static int32_t first_phi;
static int32_t first_scatt;
static int32_t first_dtnwave;
static int32_t first_dtntheta;
static int32_t first_npc;
static int32_t first_ntau_870;

static float pi = PI;
static double radeg = RADEG;
static float p0 = STDPR;

static int have_ms = 0;
static int have_rh = 0;
static int have_sd = 0;
static int use_rh = 0;
static int use_netcdf = 0;
static int use_pca_lut = 0;

static int32_t Nbands;
static int32_t Maxband; /* must be >= NBANDS */

float *noise_global;

typedef struct aermod_struct {
    char name[32];
    float rh;
    int sd;

    /* angstrom exponent (nbands+1)*/
    float *angstrom;

    /* single-scattering albedo(nbands+1), extinction coefficient(nbands+1), phase function */
    float *albedo;
    float *extc;
    float **phase;

    /* quadratic coefficients for SS to MS relationship */
    float *acost;
    float *bcost;
    float *ccost;

    /* cubic coefficients for ms_epsilon atmospheric correction ..ZA */
    float *ams_all;
    float *bms_all;
    float *cms_all;
    float *dms_all;
    float *ems_all;


    /* Rayleigh-aerosol diffuse transmittance coeffs */
    float **dtran_a;
    float **dtran_b;

    /* derived quantities */
    float **lnphase;
    float **d2phase;

    /* PCA table variables */
    float *****pc_rhoa;         // pc_rhoa(tau_870, solz, phi, senz, pc)
    float **pc_components_rhoa; // pc_components_rhoa(pc, wave)
    float *pc_mean_rhoa;        // pc_mean_rhoa(wave)
    float ***pc_td;             // pc_td(tau_870, solz, pc)
    float **pc_components_td;   // pc_components_td(pc, wave)
    float *pc_mean_td;          // pc_mean_td(wave)
    float *tau_870;             // tau_870(tau_870)
    float *pc;                  // pc(pc)

} aermodstr;
    // sorting in ascending way elements of chi-squared --  will need clean-up
struct str { float value;int index;};
// 
static int cmp(const void *a,const void *b)
    {
        struct str *a1 = (struct str *)a;
        struct str *a2 = (struct str*)b;
        if((*a1).value<(*a2).value)return -1;
        else if((*a1).value>(*a2).value)return 1;
        else return 0;
    }
aermodstr* alloc_aermodstr(int nbands, int nscatt, int nphi, int nsolz, int nsenz, int ntheta, int npc, int ntau_870) {
    aermodstr *model;

    // make sure model starts out all zeros/null
    model = (aermodstr *) calloc(1, sizeof(aermodstr));
    model->angstrom = (float*) malloc(nbands * sizeof (float));
    model->extc = (float*) malloc(nbands * sizeof (float));
    if(use_pca_lut) {
        // pc_rhoa(tau_870, solz, phi, senz, pc)
        model->pc_rhoa = allocate5d_float(ntau_870, nsolz, nphi, nsenz, npc);

        // pc_components_rhoa(pc, wave)
        model->pc_components_rhoa = allocate2d_float(npc, nbands); 
        
        // pc_mean_rhoa(wave)
        model->pc_mean_rhoa = (float*) malloc(nbands * sizeof(float));
        
        // pc_td(tau_870, solz, pc)
        model->pc_td = allocate3d_float(ntau_870,nsenz , npc);//nsolz
        
        // pc_components_td(pc, wave)
        model->pc_components_td = allocate2d_float(npc, nbands);   
        
        // pc_mean_td(wave)
        model->pc_mean_td = (float*) malloc(nbands * sizeof(float));
        
        // tau_870(tau_870)
        model->tau_870 = (float*) malloc(ntau_870 * sizeof(float));             

        // pc(pc)
        model->pc = (float*) malloc(npc * sizeof(float));
    } else {
        if(nscatt > 0) {
            model->albedo = (float*) malloc(nbands * sizeof (float));
            model->phase = allocate2d_float(nbands, nscatt);
            model->lnphase = allocate2d_float(nbands, nscatt);
            model->d2phase = allocate2d_float(nbands, nscatt);
            model->acost = (float*) malloc(nbands * nsolz * nphi * nsenz * sizeof (float));
            model->bcost = (float*) malloc(nbands * nsolz * nphi * nsenz * sizeof (float));
            model->ccost = (float*) malloc(nbands * nsolz * nphi * nsenz * sizeof (float));
        }
        model->dtran_a = allocate2d_float(nbands, ntheta);
        model->dtran_b = allocate2d_float(nbands, ntheta);
        model->ams_all = (float*) malloc(nbands * nsolz * nphi * nsenz * sizeof (float));
        model->bms_all = (float*) malloc(nbands * nsolz * nphi * nsenz * sizeof (float));
        model->cms_all = (float*) malloc(nbands * nsolz * nphi * nsenz * sizeof (float));
        model->dms_all = (float*) malloc(nbands * nsolz * nphi * nsenz * sizeof (float));
        model->ems_all = (float*) malloc(nbands * nsolz * nphi * nsenz * sizeof (float));
    }
    return model;
}

typedef struct aermodtab_struct {
    int32_t sensorID;

    /* table dimensions */
    int32_t nwave;
    int32_t nmodel;
    int32_t nsolz;
    int32_t nsenz;
    int32_t nphi;
    int32_t nscatt;
    int32_t dtran_nwave;
    int32_t dtran_ntheta;
    int32_t npc;
    int32_t ntau_870;

    /* table spectral bands and angles */
    float *wave;
    float *solz;
    float *senz;
    float *phi;
    float *scatt;

    /* diffuse transmittance spectral bands and angles */
    float *dtran_wave;
    float *dtran_theta;
    float *dtran_airmass;

    aermodstr **model;

} aermodtabstr;

typedef struct alphaT_struct {
    int32_t modnum;
    float angstrom;
} alphaTstr;

typedef struct epsilonT_struct {
    int32_t modnum;
    float eps_obs;
} epsilonTstr;

/* aerosol table */
static aermodtabstr *aertab = NULL;

/* structure for carrying the geometry information */
typedef struct geom_strdef {
    int gmult; /* band offset multiplier: 0 for nominal geometry
                 1 for band-dependent geometry  */
    float *senz;
    float *solz;
    float *phi;
    float *csolz; /* cosine of solz */
    float *csenz; /* cosine of senz */
    float *airmass;
    float *airmass_plp;
    float *airmass_sph;
} geom_str;

geom_str geom;

/* global variable declarations */
static int loaded = 0;
static int interpol = 0;
static int32_t *iwatab;
static int32_t *iwdtab;

static int32_t iwnir_s = -1;
static int32_t iwnir_l = -1;

static float mu0;
static float mu;
static float airmass;

static int32_t evalmask = 0;
static int32_t aer_opt = 0;
static float airmass_plp;
static float airmass_sph;

int cmpfunc(const void * a, const void * b) {
    if (*(double*) a > *(double*) b) return 1;
    else if (*(double*) a < *(double*) b) return -1;
    else return 0;
}



/* ---------------------------------------------------------------------------------------- */
/* first_deriv() - returns first derivative (dy/dx) of 1st or last array indices using a    */
/*                 4-pt Lagrangian interpolation.  NOTE: It is assumed that 4 points exist. */

/* ---------------------------------------------------------------------------------------- */
float first_deriv(float x[], float y[], int n) {
    float a1, a2, a3, a4, a5, a6, d1;

    if (n == 0) {

        a1 = x[0] - x[1];
        a2 = x[0] - x[2];
        a3 = x[0] - x[3];
        a4 = x[1] - x[2];
        a5 = x[1] - x[3];
        a6 = x[2] - x[3];

        d1 = y[0]*(1.0 / a1 + 1.0 / a2 + 1.0 / a3)
                - a2 * a3 * y[1] / (a1 * a4 * a5)
                + a1 * a3 * y[2] / (a2 * a4 * a6)
                - a1 * a2 * y[3] / (a3 * a5 * a6);

    } else {

        a1 = x[n - 1] - x[n - 4];
        a2 = x[n - 1] - x[n - 3];
        a3 = x[n - 1] - x[n - 2];
        a4 = x[n - 2] - x[n - 4];
        a5 = x[n - 2] - x[n - 3];
        a6 = x[n - 3] - x[n - 4];

        d1 = -a2 * a3 * y[n - 4] / (a6 * a4 * a1)
                + a1 * a3 * y[n - 3] / (a6 * a5 * a2)
                - a1 * a2 * y[n - 2] / (a4 * a5 * a3)
                + y[n - 1]*(1.0 / a1 + 1.0 / a2 + 1.0 / a3);
    }

    return (d1);
}


static void read_dimension_size(const char* file_name, int nc_id, const char* dim_name, int32_t *length) {
    int dim_id;
    size_t dim_length;
    int status;

    status = nc_inq_dimid (nc_id, dim_name, &dim_id);
    if (status != NC_NOERR) {
        printf("-E- %s:  Error looking for dimension %s from %s.\n", __FILE__, dim_name, file_name);
        exit(1);
    }
    status = nc_inq_dimlen(nc_id, dim_id, &dim_length);
    if (status != NC_NOERR) {
        printf("-E- %s:  Error reading dimension %s from %s.\n", __FILE__, dim_name, file_name);
        exit(1);
    }
    *length = (int32_t) dim_length;
}


static void read_lut_variable(const char* file, int nc_id, const char* var_name, float* data) {
    int var_id;
    int status;

    status = nc_inq_varid(nc_id, var_name, &var_id);
    if (status != NC_NOERR) {
        printf("-E- %s:  Error looking for variable %s from %s.\n", __FILE__, var_name, file);
        exit(1);
    }
    status = nc_get_var_float(nc_id, var_id, data);
    if (status != NC_NOERR) {
        printf("-E- %s:  Error reading variable %s from %s.\n", __FILE__, var_name, file);
        exit(1);
    }
}

// static void read_lut_variable_short(const char* file, int nc_id, const char* var_name, short* data) {
//     int var_id;
//     int status;

//     status = nc_inq_varid(nc_id, var_name, &var_id);
//     if (status != NC_NOERR) {
//         printf("-E- %s:  Error looking for variable %s from %s.\n", __FILE__, var_name, file);
//         exit(1);
//     }
//     status = nc_get_var_short(nc_id, var_id, data);
//     if (status != NC_NOERR) {
//         printf("-E- %s:  Error reading variable %s from %s.\n", __FILE__, var_name, file);
//         exit(1);
//     }
// }

/* ---------------------------------------------------------------------------------------- */
/* load_aermod() - loads the entire aerosol model table for the specified model list        */

/* ---------------------------------------------------------------------------------------- */
int load_aermod(int32_t sensorID, float wave[], int32_t nwave, char *aermodfile, char models[MAXAERMOD][32], int32_t nmodels) {
    int nc_id;
    int status;

    int32_t aer_nwave, nsolz, nsenz, nphi, nscatt, dtran_nwave, dtran_ntheta, npc, ntau_870;

    float d1phase1;
    float d1phaseN;
    float rh;
    int16_t sd;

    char file [FILENAME_MAX] = "";
    char path [FILENAME_MAX] = "";

    int iw, im, is, iwbase, i;
    static int firstCall = 1;

    if (firstCall == 1) {
        if ((iwatab = (int32_t *) calloc(nwave, sizeof (int32_t))) == NULL) {
            printf("Unable to allocate space for iwatab.\n");
            exit(1);
        }
        if ((iwdtab = (int32_t *) calloc(nwave, sizeof (int32_t))) == NULL) {
            printf("Unable to allocate space for iwdtab.\n");
            exit(1);
        }
        firstCall = 0;
    }

    printf("Loading aerosol models from %s\n", aermodfile);

    for (im = 0; im < nmodels + 1; im++) {

        if (im < nmodels) {

            // try netCDF first
            strcpy(file, path);
            strcat(file, aermodfile);
            strcat(file, "_");
            strcat(file, models[im]);
            strcat(file, ".nc");
            use_netcdf = 1;

            // if not try HDF4
            if(access(file, R_OK) == -1) {
                strcpy(file, path);
                strcat(file, aermodfile);
                strcat(file, "_");
                strcat(file, models[im]);
                strcat(file, ".hdf");
                use_netcdf = 0;
            }

        } else {

            // first try netCDF
            strcpy(file, path);
            strcat(file, aermodfile);
            strcat(file, "_default.nc");
            use_netcdf = 1;

            // if not try HDF4
            if(access(file, R_OK) == -1) {
                strcpy(file, path);
                strcat(file, aermodfile);
                strcat(file, "_default.hdf");
                use_netcdf = 0;
            }
        }

        status = nc_open(file, NC_NOWRITE, &nc_id);
        if (status != NC_NOERR) {
            printf("-E- %s:  Error opening file %s.\n", __FILE__, file);
            exit(1);
        }

        use_pca_lut = 0;
        nscatt = 0;
        dtran_nwave = 0;
        dtran_ntheta = 0;
        npc = 0;
        ntau_870 = 0;

        /* read dimensions which should be constant between models */
        if(use_netcdf) {
            //see if it is a PCA LUT
            size_t attlen;
            status = nc_inq_attlen(nc_id, NC_GLOBAL, "title", &attlen);
            if (status != NC_NOERR) {
                printf("-E- %s: could not find title in %s.\n", __FILE__, file);
                exit(1);
            }
            char* title = malloc(sizeof(char) * attlen);
            nc_get_att_text(nc_id, NC_GLOBAL, "title", title);
            if(strcasestr(title, "pca aerosol model data")) {
                use_pca_lut = 1;
            }
            free(title);

            if(use_pca_lut) {
                read_dimension_size(file, nc_id, "wavelength", &aer_nwave);
                read_dimension_size(file, nc_id, "solar_zenith", &nsolz);
                read_dimension_size(file, nc_id, "sensor_zenith", &nsenz);
                read_dimension_size(file, nc_id, "relative_azimuth", &nphi);
                read_dimension_size(file, nc_id, "principal_component", &npc);
                read_dimension_size(file, nc_id, "aerosol_optical_thickness", &ntau_870);
            } else {
                read_dimension_size(file, nc_id, "nwave", &aer_nwave);
                read_dimension_size(file, nc_id, "nsolz", &nsolz);
                read_dimension_size(file, nc_id, "nsenz", &nsenz);
                read_dimension_size(file, nc_id, "nphi", &nphi);
                read_dimension_size(file, nc_id, "dtran_nwave", &dtran_nwave);
                read_dimension_size(file, nc_id, "dtran_ntheta", &dtran_ntheta);
            }
         } else {
            read_dimension_size(file, nc_id, "nwave", &aer_nwave);
            read_dimension_size(file, nc_id, "nsolz", &nsolz);
            read_dimension_size(file, nc_id, "nsenz", &nsenz);
            read_dimension_size(file, nc_id, "nphi", &nphi);
            read_dimension_size(file, nc_id, "nscatt", &nscatt);
            read_dimension_size(file, nc_id, "dtran_nwave", &dtran_nwave);
            read_dimension_size(file, nc_id, "dtran_ntheta", &dtran_ntheta);
         }
        // if(aer_nwave != nwave) {
        //     printf("-E- %s:%d -  Error nwave dimension = %d not equal to nwave = %d in file = %s.\n", 
        //         __FILE__, __LINE__, aer_nwave, nwave, file);
        //     exit(1);
        // }

        if (im == 0) {
            printf("Number of Wavelengths                          %d\n", aer_nwave);
            printf("Number of Solar Zenith Angles                  %d\n", nsolz);
            printf("Number of View Zenith Angles                   %d\n", nsenz);
            printf("Number of Relative Azimuth Angles              %d\n", nphi);
            if(use_pca_lut) {
                printf("Number of Principal Components                 %d\n", npc);
                printf("Number of tau 870 Angles                       %d\n", ntau_870);
            } else {
                printf("Number of Scattering Angles                    %d\n", nscatt);
                printf("Number of Diffuse Transmittance Wavelengths    %d\n", dtran_nwave);
                printf("Number of Diffuse Transmittance Zenith Angles  %d\n", dtran_ntheta);
            }

            first_solz = nsolz;
            first_senz = nsenz;
            first_phi = nphi;
            first_scatt = nscatt;
            first_dtnwave = dtran_nwave;
            first_dtntheta = dtran_ntheta;
            first_npc = npc;
            first_ntau_870 = ntau_870;

            // allocate the aerosol table
            if ((aertab = (aermodtabstr *) calloc(1, sizeof (aermodtabstr))) == NULL) {
                printf("Unable to allocate space for aerosol table.\n");
                exit(1);
            }

            aertab->nmodel = nmodels;
            aertab->nwave = aer_nwave;
            aertab->nsolz = nsolz;
            aertab->nsenz = nsenz;
            aertab->nphi = nphi;
            aertab->nscatt = nscatt;
            aertab->dtran_nwave = dtran_nwave;
            aertab->dtran_ntheta = dtran_ntheta;

            if(use_pca_lut){
                aertab->npc = npc;
                aertab->ntau_870 = ntau_870;
            }

            aertab->wave = (float *) malloc(aer_nwave * sizeof (float));
            aertab->solz = (float *) malloc(nsolz * sizeof (float));
            aertab->senz = (float *) malloc(nsenz * sizeof (float));
            aertab->phi = (float *) malloc(nphi * sizeof (float));
            if(use_netcdf)
                aertab->scatt = NULL;
            else
                aertab->scatt = (float *) malloc(nscatt * sizeof (float));

            if(!use_pca_lut) {
                aertab->dtran_wave = (float *) malloc(dtran_nwave * sizeof (float));
                aertab->dtran_theta = (float *) malloc(dtran_ntheta * sizeof (float));
                aertab->dtran_airmass = (float *) malloc(dtran_ntheta * sizeof (float));
            }

            // allocate the model tables
            if ((aertab->model = (aermodstr **) calloc(1, (nmodels + 1) * sizeof (aermodstr*))) == NULL) {
                printf("Unable to allocate space for %d aerosol models.\n", nmodels + 1);
                exit(1);
            }
            for (i = 0; i < nmodels + 1; i++) {
                if ((aertab->model[i] = alloc_aermodstr(aer_nwave, nscatt, nphi, nsolz, nsenz, dtran_ntheta, npc, ntau_870)) == NULL) {
                    printf("Unable to allocate space for aerosol model %d.\n", im);
                    exit(1);
                }
            }

            /* read SDSes which are constant between models */

            if(!use_netcdf)
                read_lut_variable(file, nc_id, "scatt", aertab->scatt);
            if(!use_netcdf){
                read_lut_variable(file, nc_id, "wave", aertab->wave);
                read_lut_variable(file, nc_id, "solz", aertab->solz);
                read_lut_variable(file, nc_id, "senz", aertab->senz);
                read_lut_variable(file, nc_id, "phi", aertab->phi);
            }
            else{
                if(use_pca_lut){
                read_lut_variable(file, nc_id, "wavelength", aertab->wave);
                read_lut_variable(file, nc_id, "solar_zenith", aertab->solz);
                read_lut_variable(file, nc_id, "sensor_zenith", aertab->senz);
                read_lut_variable(file, nc_id, "relative_azimuth", aertab->phi);
                }
                else{
                    read_lut_variable(file, nc_id, "wave", aertab->wave);
                    read_lut_variable(file, nc_id, "solz", aertab->solz);
                    read_lut_variable(file, nc_id, "senz",aertab->senz);
                    read_lut_variable(file, nc_id, "phi", aertab->phi);
                }
            }
            if(!use_pca_lut) {
                read_lut_variable(file, nc_id, "dtran_wave", aertab->dtran_wave);
                read_lut_variable(file, nc_id, "dtran_theta", aertab->dtran_theta);
            }
        } else {
            /*  check that all the aerosol files at least have the same
                main dimensions  */
            if ((aertab->nsolz != first_solz) || (aertab->nsenz != first_senz) || 
                (aertab->nphi != first_phi) || (aertab->nscatt != first_scatt) ||
                (aertab->dtran_nwave != first_dtnwave) || (aertab->dtran_ntheta != first_dtntheta) ||
                (use_pca_lut && aertab->ntau_870 != first_ntau_870) || (use_pca_lut && aertab->npc != first_npc)) {
                printf("-E- %s, %d:  Error, Aerosol table %s\n",
                        __FILE__, __LINE__, file);
                printf("    has different dimensions from previous tables\n");
                exit(1);
            }
        }

        if (im < nmodels)
            strncpy(aertab->model[im]->name, models[im], 32);
        else
            strncpy(aertab->model[im]->name, "default", 32);

        status = nc_get_att_float(nc_id, NC_GLOBAL, "RelativeHumidity", &rh);
        if(status != NC_NOERR) {
            status = nc_get_att_float(nc_id, NC_GLOBAL, "Relative Humidity", &rh);
        }
        if(status != NC_NOERR) {
            status = nc_get_att_float(nc_id, NC_GLOBAL, "relative_humidity", &rh);
        }
        if (status == NC_NOERR) {
            if(use_pca_lut)
                rh=rh*100;
            aertab->model[im]->rh = rh;
            have_rh = 1;
        } else {
            aertab->model[im]->rh = -1.0;
            have_rh = 0;
        }

        if(use_netcdf) {
            float fmf;
            status = nc_get_att_float(nc_id, NC_GLOBAL, "AerosolFMF", &fmf);
            if(status != NC_NOERR) {
                status = nc_get_att_float(nc_id, NC_GLOBAL, "fine_mode_fraction", &fmf);
                if(status!=NC_NOERR){
                    printf("-E- %s, %d:  Error, Aerosol table %s does not have AerosolFMF\n",
                            __FILE__, __LINE__, file);
                    exit(1);
                }
            }
            sd = 100 - fmf*100;
            aertab->model[im]->sd = sd;
        } else {
            status = nc_get_att_short(nc_id, NC_GLOBAL, "Size Distribution", &sd);
            if (status == NC_NOERR) {
                aertab->model[im]->sd = sd;
                have_sd = 1;
            } else {
                aertab->model[im]->sd = -1;
                have_sd = 0;
            }
        }

        if(!use_netcdf)
            read_lut_variable(file, nc_id, "extc", aertab->model[im]->extc);
        else{
            if(use_pca_lut)
                read_lut_variable(file, nc_id, "extinction_coefficient", aertab->model[im]->extc);
            else
                read_lut_variable(file, nc_id, "extc", aertab->model[im]->extc);
        }

        if(use_pca_lut) {
            read_lut_variable(file, nc_id, "aerosol_reflectance_score", aertab->model[im]->pc_rhoa[0][0][0][0]);
            read_lut_variable(file, nc_id, "aerosol_reflectance_component", aertab->model[im]->pc_components_rhoa[0]);
            read_lut_variable(file, nc_id, "aerosol_reflectance_mean", aertab->model[im]->pc_mean_rhoa);
            read_lut_variable(file, nc_id, "diffuse_transmittance_score", aertab->model[im]->pc_td[0][0]);
            read_lut_variable(file, nc_id, "diffuse_transmittance_component", aertab->model[im]->pc_components_td[0]);
            read_lut_variable(file, nc_id, "diffuse_transmittance_mean", aertab->model[im]->pc_mean_td);
            read_lut_variable(file, nc_id, "angstrom", aertab->model[im]->angstrom);
            read_lut_variable(file, nc_id, "aerosol_optical_thickness", aertab->model[im]->tau_870);
        } else {
            if(!use_netcdf) {
                read_lut_variable(file, nc_id, "albedo", aertab->model[im]->albedo);
                read_lut_variable(file, nc_id, "phase", aertab->model[im]->phase[0]);
                read_lut_variable(file, nc_id, "acost", aertab->model[im]->acost);
                read_lut_variable(file, nc_id, "bcost", aertab->model[im]->bcost);
                read_lut_variable(file, nc_id, "ccost", aertab->model[im]->ccost);
            }

            // check for multi-scattering epsilon tables
            if(nc_inq_varid(nc_id, "ams_all", &status) == NC_NOERR) {
                have_ms = 1;
                read_lut_variable(file, nc_id, "ams_all", aertab->model[im]->ams_all);
                read_lut_variable(file, nc_id, "bms_all", aertab->model[im]->bms_all);
                read_lut_variable(file, nc_id, "cms_all", aertab->model[im]->cms_all);
                read_lut_variable(file, nc_id, "dms_all", aertab->model[im]->dms_all);
                read_lut_variable(file, nc_id, "ems_all", aertab->model[im]->ems_all);
            }

            read_lut_variable(file, nc_id, "dtran_a", aertab->model[im]->dtran_a[0]);
            read_lut_variable(file, nc_id, "dtran_b", aertab->model[im]->dtran_b[0]);
        }

        nc_close(nc_id);

        if(!use_pca_lut) {
            /* compute angstrom exponent for each model wavelength relative to max wavelength */
            iwbase = windex(865, aertab->wave, aertab->nwave);
            for (iw = 0; iw < aertab->nwave; iw++) {
                if (iw != iwbase)
                    aertab->model[im]->angstrom[iw] = -log(aertab->model[im]->extc[iw] / aertab->model[im]->extc[iwbase]) /
                    log(aertab->wave[iw] / aertab->wave[iwbase]);
            }
            aertab->model[im]->angstrom[iwbase] = aertab->model[im]->angstrom[iwbase - 1];

            /* precompute log of phase function and 2nd derivative (for cubic interp) */
            if(aertab->scatt) {
                for (iw = 0; iw < aertab->nwave; iw++) {
                    for (is = 0; is < aertab->nscatt; is++) {
                        aertab->model[im]->lnphase[iw][is] = log(aertab->model[im]->phase[iw][is]);
                    }
                    d1phase1 = first_deriv(aertab->scatt, &aertab->model[im]->lnphase[iw][0], 0);
                    d1phaseN = first_deriv(aertab->scatt, &aertab->model[im]->lnphase[iw][0], aertab->nscatt);
                    spline(aertab->scatt,
                            &aertab->model[im]->lnphase[iw][0],
                            aertab->nscatt,
                            d1phase1,
                            d1phaseN,
                            &aertab->model[im]->d2phase[iw][0]);
                }
            }
        }

    }

    aertab->nmodel = nmodels;
    aertab->sensorID = sensorID;

    /* Require number of table wavelengths to equal number of sensor wavelengths */
    if (aertab->nwave != nwave) {
        printf("Number of aerosol LUT wavelengths (%d) is not equal to number of sensor wavelengths (%d).\n", 
            aertab->nwave,nwave);
        exit(1);
    }

    if(!use_pca_lut) {
        /* precompute airmass for diffuse transmittance */
        for (is = 0; is < aertab->dtran_ntheta; is++) {
            aertab->dtran_airmass[is] = 1.0 / cos(aertab->dtran_theta[is] / radeg);
        }

        if (aertab->dtran_nwave != nwave) {
            printf("Number of aerosol diffue trans LUT wavelengths (%d) is not equal to number of sensor wavelengths (%d).\n", 
                aertab->dtran_nwave,nwave);
            exit(1);
        }

        /* map sensor wavelengths to table wavelengths */
        printf("Wavelengths - Sensor\tAerosol model\tDiffuse transmittance\n");
        for (iw = 0; iw < nwave; iw++) {
            iwatab[iw] = windex(wave[iw], aertab->wave, aertab->nwave);
            iwdtab[iw] = windex(wave[iw], aertab->dtran_wave, aertab->dtran_nwave);
            printf("\t      %6.2f\t    %6.2f\t   %6.2f\n", wave[iw],aertab->wave[iwatab[iw]], aertab->dtran_wave[iwdtab[iw]]);
        }

    }

    loaded = 1;

    return (0);
}



#define INDEX(iw,isol,iphi,isen) (iw*aertab->nsolz*aertab->nphi*aertab->nsenz + isol*aertab->nphi*aertab->nsenz + iphi*aertab->nsenz + isen)

/* ---------------------------------------------------------------------------------------- */
/* ss_to_ms_coef() - for one model, return coefficients of function relating single         */
/*                   scattering to multiple scattering at the input geometry.               */
/*                                                                                          */
/* This is effectively a C version of M. Wangs linear_a_b_c.f.  The program optimizes for   */
/* multiple calls at the same geometry by computing for all models on the first call with a */
/* new geometry.  It returns pointers to the internal static arrays of coefficients.        */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */
/* W. Robinson, SAIC  adapt to band-ependent geometry                                       */

/* ---------------------------------------------------------------------------------------- */
void ss_to_ms_coef(int modnum, geom_str *geom, float **a, float **b, float **c) {
    static float lastsolz = -999.;
    static float lastsenz = -999.;
    static float lastphi = -999.;

    static int computed[MAXMODEL];

    static float *a_coef[MAXMODEL];
    static float *b_coef[MAXMODEL];
    static float *c_coef[MAXMODEL];

    static float *p, *q, *r;
    static float as000, as100, as010, as110, as001, as011, as101, as111;
    static float ai000, ai100, ai010, ai110, ai001, ai011, ai101, ai111;
    static float ac000, ac100, ac010, ac110, ac001, ac011, ac101, ac111;

    static int *isolz1, *isolz2;
    static int *isenz1, *isenz2;
    static int *iphi1, *iphi2;

    static float *p_ar, *q_ar, *r_ar;
    static float p_cnst, q_cnst, r_cnst;
    static int *isolz1_ar, *isolz2_ar, isolz1_cnst, isolz2_cnst;
    static int *isenz1_ar, *isenz2_ar, isenz1_cnst, isenz2_cnst;
    static int *iphi1_ar, *iphi2_ar, iphi1_cnst, iphi2_cnst;

    float aphi;
    float px, qx, rx;
    int im, iw, i, ig;
    static int firstCall = 1;
    static int gmult;

    if (firstCall == 1) {
        firstCall = 0;
        for (i = 0; i < MAXMODEL; i++) {
            if ((a_coef[i] = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for a_coef.\n");
                exit(1);
            }
            if ((b_coef[i] = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for rhoa.\n");
                exit(1);
            }
            if ((c_coef[i] = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for rhoa.\n");
                exit(1);
            }
        }
        /*  set up indicies, weights for band-dependent or nominal geometry */
        if ((geom->gmult == 0) || (interpol == 1)) {
            gmult = 0;
            p = &p_cnst;
            q = &q_cnst;
            r = &r_cnst;
            isolz1 = &isolz1_cnst;
            isolz2 = &isolz2_cnst;
            isenz1 = &isenz1_cnst;
            isenz2 = &isenz2_cnst;
            iphi1 = &iphi1_cnst;
            iphi2 = &iphi2_cnst;
        } else
            gmult = 1;
        {
            if (((p_ar = (float *) malloc(aertab->nwave * sizeof (float)))
                    == NULL) ||
                    ((q_ar = (float *) malloc(aertab->nwave * sizeof (float)))
                    == NULL) ||
                    ((r_ar = (float *) malloc(aertab->nwave * sizeof (float)))
                    == NULL)) {
                printf("Unable to allocate space for p, q, r weights.\n");
                exit(1);
            }
            if (((isolz1_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL) ||
                    ((isenz1_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL) ||
                    ((iphi1_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL)) {
                printf("Unable to allocate space for interp indicies 1.\n");
                exit(1);
            }
            if (((isolz2_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL) ||
                    ((isenz2_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL) ||
                    ((iphi2_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL)) {
                printf("Unable to allocate space for interp indicies 2.\n");
                exit(1);
            }
            p = p_ar;
            q = q_ar;
            r = r_ar;
            isolz1 = isolz1_ar;
            isolz2 = isolz2_ar;
            isenz1 = isenz1_ar;
            isenz2 = isenz2_ar;
            iphi1 = iphi1_ar;
            iphi2 = iphi2_ar;
        }
    }

    if ((geom->solz[0] != lastsolz) || (geom->senz[0] != lastsenz) ||
            (geom->phi[0] != lastphi)) {
        for (im = 0; im < aertab->nmodel; im++)
            computed[im] = 0;

        for (iw = 0; iw < aertab->nwave; iw++) {
            ig = iw * gmult;
            /* find bracketing solar indices */
            for (i = 0; i < aertab->nsolz; i++) {
                if (geom->solz[ig] < aertab->solz[i])
                    break;
            }
            isolz1[iw] = MAX(i - 1, 0);
            isolz2[iw] = MIN(i, aertab->nsolz - 1);
            if (isolz2[iw] != isolz1[iw])
                r[iw] = (geom->solz[ig] - aertab->solz[isolz1[iw]]) /
                (aertab->solz[isolz2[iw]] - aertab->solz[isolz1[iw]]);
            else
                r[iw] = 0.0;

            /* find bracketing view indices */
            for (i = 0; i < aertab->nsenz; i++) {
                if (geom->senz[ig] < aertab->senz[i])
                    break;
            }
            isenz1[iw] = MAX(i - 1, 0);
            isenz2[iw] = MIN(i, aertab->nsenz - 1);
            if (isenz2[iw] != isenz1[iw])
                p[iw] = (geom->senz[ig] - aertab->senz[isenz1[iw]]) /
                (aertab->senz[isenz2[iw]] - aertab->senz[isenz1[iw]]);
            else
                p[iw] = 0.0;

            /* find bracketing azimuth indices */
            aphi = fabs(geom->phi[ig]);
            for (i = 0; i < aertab->nphi; i++) {
                if (aphi < aertab->phi[i])
                    break;
            }
            iphi1[iw] = MAX(i - 1, 0);
            iphi2[iw] = MIN(i, aertab->nphi - 1);
            if (iphi2[iw] != iphi1[iw])
                q[iw] = (aphi - aertab->phi[iphi1[iw]]) /
                (aertab->phi[iphi2[iw]] - aertab->phi[iphi1[iw]]);
            else
                q[iw] = 0.0;
            if (gmult == 0) break;
        }

        /* remember last geometry */
        lastsolz = geom->solz[0];
        lastsenz = geom->senz[0];
        lastphi = geom->phi[0];
    }

    if (!computed[modnum]) {
        im = modnum;
        computed[modnum] = 1;

        for (iw = 0; iw < aertab->nwave; iw++) {
            ig = iw * gmult;
            px = p[ig];
            qx = q[ig];
            rx = r[ig];
            if (isolz2[ig] == 0) {
                as000 = aertab->model[im]->acost[INDEX(iw, isolz1[ig], 0, isenz1[ig])];
                as100 = aertab->model[im]->acost[INDEX(iw, isolz1[ig], 0, isenz2[ig])];
                as001 = aertab->model[im]->acost[INDEX(iw, isolz2[ig], 0, isenz1[ig])];
                as101 = aertab->model[im]->acost[INDEX(iw, isolz2[ig], 0, isenz2[ig])];

                ai000 = aertab->model[im]->bcost[INDEX(iw, isolz1[ig], 0, isenz1[ig])];
                ai100 = aertab->model[im]->bcost[INDEX(iw, isolz1[ig], 0, isenz2[ig])];
                ai001 = aertab->model[im]->bcost[INDEX(iw, isolz2[ig], 0, isenz1[ig])];
                ai101 = aertab->model[im]->bcost[INDEX(iw, isolz2[ig], 0, isenz2[ig])];

                ac000 = aertab->model[im]->ccost[INDEX(iw, isolz1[ig], 0, isenz1[ig])];
                ac100 = aertab->model[im]->ccost[INDEX(iw, isolz1[ig], 0, isenz2[ig])];
                ac001 = aertab->model[im]->ccost[INDEX(iw, isolz2[ig], 0, isenz1[ig])];
                ac101 = aertab->model[im]->ccost[INDEX(iw, isolz2[ig], 0, isenz2[ig])];

                a_coef[im][iw] = (1. - px)*(1. - rx) * as000 + px * rx * as101
                        + (1. - px) * rx * as001 + px * (1. - rx) * as100;

                b_coef[im][iw] = (1. - px)*(1. - rx) * ai000 + px * rx * ai101
                        + (1. - px) * rx * ai001 + px * (1. - qx)*(1. - rx) * ai100;

                c_coef[im][iw] = (1. - px)*(1. - rx) * ac000 + px * rx * ac101
                        + (1. - px) * rx * ac001 + px * (1. - qx)*(1. - rx) * ac100;
            } else {
            	as000 = aertab->model[im]->acost[INDEX(iw, isolz1[ig], iphi1[ig], isenz1[ig])];
                as100 = aertab->model[im]->acost[INDEX(iw, isolz1[ig], iphi1[ig], isenz2[ig])];
                as010 = aertab->model[im]->acost[INDEX(iw, isolz1[ig], iphi2[ig], isenz1[ig])];
                as110 = aertab->model[im]->acost[INDEX(iw, isolz1[ig], iphi2[ig], isenz2[ig])];
                as001 = aertab->model[im]->acost[INDEX(iw, isolz2[ig], iphi1[ig], isenz1[ig])];
                as011 = aertab->model[im]->acost[INDEX(iw, isolz2[ig], iphi2[ig], isenz1[ig])];
                as101 = aertab->model[im]->acost[INDEX(iw, isolz2[ig], iphi1[ig], isenz2[ig])];
                as111 = aertab->model[im]->acost[INDEX(iw, isolz2[ig], iphi2[ig], isenz2[ig])];

                ai000 = aertab->model[im]->bcost[INDEX(iw, isolz1[ig], iphi1[ig], isenz1[ig])];
                ai100 = aertab->model[im]->bcost[INDEX(iw, isolz1[ig], iphi1[ig], isenz2[ig])];
                ai010 = aertab->model[im]->bcost[INDEX(iw, isolz1[ig], iphi2[ig], isenz1[ig])];
                ai110 = aertab->model[im]->bcost[INDEX(iw, isolz1[ig], iphi2[ig], isenz2[ig])];
                ai001 = aertab->model[im]->bcost[INDEX(iw, isolz2[ig], iphi1[ig], isenz1[ig])];
                ai011 = aertab->model[im]->bcost[INDEX(iw, isolz2[ig], iphi2[ig], isenz1[ig])];
                ai101 = aertab->model[im]->bcost[INDEX(iw, isolz2[ig], iphi1[ig], isenz2[ig])];
                ai111 = aertab->model[im]->bcost[INDEX(iw, isolz2[ig], iphi2[ig], isenz2[ig])];

                ac000 = aertab->model[im]->ccost[INDEX(iw, isolz1[ig], iphi1[ig], isenz1[ig])];
                ac100 = aertab->model[im]->ccost[INDEX(iw, isolz1[ig], iphi1[ig], isenz2[ig])];
                ac010 = aertab->model[im]->ccost[INDEX(iw, isolz1[ig], iphi2[ig], isenz1[ig])];
                ac110 = aertab->model[im]->ccost[INDEX(iw, isolz1[ig], iphi2[ig], isenz2[ig])];
                ac001 = aertab->model[im]->ccost[INDEX(iw, isolz2[ig], iphi1[ig], isenz1[ig])];
                ac011 = aertab->model[im]->ccost[INDEX(iw, isolz2[ig], iphi2[ig], isenz1[ig])];
                ac101 = aertab->model[im]->ccost[INDEX(iw, isolz2[ig], iphi1[ig], isenz2[ig])];
                ac111 = aertab->model[im]->ccost[INDEX(iw, isolz2[ig], iphi2[ig], isenz2[ig])];

                a_coef[im][iw] = (1. - px)*(1. - qx)*(1. - rx) * as000 + px * qx * rx * as111
                        + px * (1. - qx) * rx * as101 + (1. - px) * qx * (1. - rx) * as010
                        + px * qx * (1. - rx) * as110 + (1. - px)*(1. - qx) * rx * as001
                        + (1. - px) * qx * rx * as011 + px * (1. - qx)*(1. - rx) * as100;

                b_coef[im][iw] = (1. - px)*(1. - qx)*(1. - rx) * ai000 + px * qx * rx * ai111
                        + px * (1. - qx) * rx * ai101 + (1. - px) * qx * (1. - rx) * ai010
                        + px * qx * (1. - rx) * ai110 + (1. - px)*(1. - qx) * rx * ai001
                        + (1. - px) * qx * rx * ai011 + px * (1. - qx)*(1. - rx) * ai100;

                c_coef[im][iw] = (1. - px)*(1. - qx)*(1. - rx) * ac000 + px * qx * rx * ac111
                        + px * (1. - qx) * rx * ac101 + (1. - px) * qx * (1. - rx) * ac010
                        + px * qx * (1. - rx) * ac110 + (1. - px)*(1. - qx) * rx * ac001
                        + (1. - px) * qx * rx * ac011 + px * (1. - qx)*(1. - rx) * ac100;
            }
        }
    }

    /* return pointers to coeffs for this geometry */
    *a = &a_coef[modnum][0];
    *b = &b_coef[modnum][0];
    *c = &c_coef[modnum][0];

    return;
}

/* ---------------------------------------------------------------------------------------- */
/* fresnel_coef() - computes Fresnel reflectance coefficient for specified index of refr.   */

/* ---------------------------------------------------------------------------------------- */
float fresnel_coef(float mu, float index) {
    float sq, r2, q1;

    sq = sqrt(pow(index, 2.0) - 1.0 + pow(mu, 2.0));
    r2 = pow((mu - sq) / (mu + sq), 2.0);
    q1 = (1.0 - pow(mu, 2.0) - mu * sq) / (1.0 - pow(mu, 2.0) + mu * sq);

    return (r2 * (q1 * q1 + 1.0) / 2.0);
}



/* ---------------------------------------------------------------------------------------- */
/* ms_eps_coef() -   for a given model, returns  ams, bms, cms, dms and ems coefficients to */
/*                   compute ms_reflectance at the input geometry.                          */
/*                   Also, the program optimizes for multiple calls at the same geometry    */
/*                   by computing for all models on the first call with a new geometry.     */
/*                   It returns pointers to the internal static arrays of coefficients.     */
/*                                                                                          */
/* Z. Ahmad July 08, 2014                                                                   */
/* W. Robinson, SAIC  23 Mar 2017  add band-dependent geometry to code                      */
/* M. Zhang 06/19/2019   fixed the interpolation part                                      */

/* ---------------------------------------------------------------------------------------- */



void ms_eps_coef(int modnum, int32_t iwnir_l, float wave[], geom_str *geom,
        float **a, float **b, float **c, float **d, float **e)
 {
    static float lastsolz = -999.;
    static float lastsenz = -999.;
    static float lastphi = -999.;

    static int computed[MAXMODEL];

    static float *a_coef[MAXMODEL];
    static float *b_coef[MAXMODEL];
    static float *c_coef[MAXMODEL];
    static float *d_coef[MAXMODEL];
    static float *e_coef[MAXMODEL];

    static float *p, *q, *r;
    static float as000, as100, as010, as110, as001, as011, as101, as111;
    static float ai000, ai100, ai010, ai110, ai001, ai011, ai101, ai111;
    static float ac000, ac100, ac010, ac110, ac001, ac011, ac101, ac111;
    static float ad000, ad100, ad010, ad110, ad001, ad011, ad101, ad111;
    static float ae000, ae100, ae010, ae110, ae001, ae011, ae101, ae111;

    static int *isolz1, *isolz2;
    static int *isenz1, *isenz2;
    static int *iphi1, *iphi2;

    static float *p_ar, *q_ar, *r_ar;
    static float p_cnst, q_cnst, r_cnst;
    static int *isolz1_ar, *isolz2_ar, isolz1_cnst, isolz2_cnst;
    static int *isenz1_ar, *isenz2_ar, isenz1_cnst, isenz2_cnst;
    static int *iphi1_ar, *iphi2_ar, iphi1_cnst, iphi2_cnst;
    static int gmult;

    float aphi;
    float px, qx, rx;
    int im, iw, i, ig;
    static int firstCall = 1;

    if (firstCall == 1) {
        firstCall = 0;
        for (i = 0; i < MAXMODEL; i++) {
            if ((a_coef[i] = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for a_coef.\n");
                exit(1);
            }
            if ((b_coef[i] = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for b_coef.\n");
                exit(1);
            }
            if ((c_coef[i] = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for c_coef.\n");
                exit(1);
            }
            if ((d_coef[i] = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for d_coef.\n");
                exit(1);
            }

            if ((e_coef[i] = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for e_coef.\n");
                exit(1);
            }

        }
        /*  set up indicies, weights for band-dependent or nominal geometry */
        if ((geom->gmult == 0) || (interpol == 1)) {
            gmult = 0;
            p = &p_cnst;
            q = &q_cnst;
            r = &r_cnst;
            isolz1 = &isolz1_cnst;
            isolz2 = &isolz2_cnst;
            isenz1 = &isenz1_cnst;
            isenz2 = &isenz2_cnst;
            iphi1 = &iphi1_cnst;
            iphi2 = &iphi2_cnst;
        } else {
            gmult = 1;
            if (((p_ar = (float *) malloc(aertab->nwave * sizeof (float)))
                    == NULL) ||
                    ((q_ar = (float *) malloc(aertab->nwave * sizeof (float)))
                    == NULL) ||
                    ((r_ar = (float *) malloc(aertab->nwave * sizeof (float)))
                    == NULL)) {
                printf("Unable to allocate space for p, q, r weights.\n");
                exit(1);
            }
            if (((isolz1_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL) ||
                    ((isenz1_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL) ||
                    ((iphi1_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL)) {
                printf("Unable to allocate space for interp indicies 1.\n");
                exit(1);
            }
            if (((isolz2_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL) ||
                    ((isenz2_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL) ||
                    ((iphi2_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                    == NULL)) {
                printf("Unable to allocate space for interp indicies 2.\n");
                exit(1);
            }
            p = p_ar;
            q = q_ar;
            r = r_ar;
            isolz1 = isolz1_ar;
            isolz2 = isolz2_ar;
            isenz1 = isenz1_ar;
            isenz2 = isenz2_ar;
            iphi1 = iphi1_ar;
            iphi2 = iphi2_ar;
        }
    }

    if (geom->solz[0] != lastsolz || geom->senz[0] != lastsenz ||
            geom->phi[0] != lastphi) {
        for (im = 0; im < aertab->nmodel; im++)
            computed[im] = 0;

        for (iw = 0; iw < aertab->nwave; iw++) {
            ig = iw * gmult;
            /* find bracketing solar indices */
            for (i = 0; i < aertab->nsolz; i++) {
                if (geom->solz[ig] < aertab->solz[i])
                    break;
            }
            isolz1[iw] = MAX(i - 1, 0);
            isolz2[iw] = MIN(i, aertab->nsolz - 1);
            if (isolz2[iw] != isolz1[iw])
                r[iw] = (geom->solz[ig] - aertab->solz[isolz1[iw]]) /
                (aertab->solz[isolz2[iw]] - aertab->solz[isolz1[iw]]);
            else
                r[iw] = 0.0;

            /* find bracketing view indices */
            for (i = 0; i < aertab->nsenz; i++) {
                if (geom->senz[ig] < aertab->senz[i])
                    break;
            }
            isenz1[iw] = MAX(i - 1, 0);
            isenz2[iw] = MIN(i, aertab->nsenz - 1);
            if (isenz2[iw] != isenz1[iw])
                p[iw] = (geom->senz[ig] - aertab->senz[isenz1[iw]]) /
                (aertab->senz[isenz2[iw]] - aertab->senz[isenz1[iw]]);
            else
                p[iw] = 0.0;

            /* find bracketing azimuth indices */
            aphi = fabs(geom->phi[ig]);
            for (i = 0; i < aertab->nphi; i++) {
                if (aphi < aertab->phi[i])
                    break;
            }
            iphi1[iw] = MAX(i - 1, 0);
            iphi2[iw] = MIN(i, aertab->nphi - 1);
            if (iphi2[iw] != iphi1[iw])
                q[iw] = (aphi - aertab->phi[iphi1[iw]]) /
                (aertab->phi[iphi2[iw]] - aertab->phi[iphi1[iw]]);
            else
                q[iw] = 0.0;
            if (gmult == 0) break;
        }

        /* remember last geometry */
        lastsolz = geom->solz[0];
        lastsenz = geom->senz[0];
        lastphi = geom->phi[0];

    }

    im = modnum;

    if (!computed[modnum]) {
        im = modnum;
        computed[modnum] = 1;

        for (iw = 0; iw < aertab->nwave; iw++) {
            ig = iw * gmult;
            px = p[ig];
            qx = q[ig];
            rx = r[ig];
            if (isolz2[ig] == 0) {
                as000 = aertab->model[im]->ams_all[INDEX(iw, isolz1[ig], 0, isenz1[ig])];
                as100 = aertab->model[im]->ams_all[INDEX(iw, isolz1[ig], 0, isenz2[ig])];
                as001 = aertab->model[im]->ams_all[INDEX(iw, isolz2[ig], 0, isenz1[ig])];
                as101 = aertab->model[im]->ams_all[INDEX(iw, isolz2[ig], 0, isenz2[ig])];

                ai000 = aertab->model[im]->bms_all[INDEX(iw, isolz1[ig], 0, isenz1[ig])];
                ai100 = aertab->model[im]->bms_all[INDEX(iw, isolz1[ig], 0, isenz2[ig])];
                ai001 = aertab->model[im]->bms_all[INDEX(iw, isolz2[ig], 0, isenz1[ig])];
                ai101 = aertab->model[im]->bms_all[INDEX(iw, isolz2[ig], 0, isenz2[ig])];

                ac000 = aertab->model[im]->cms_all[INDEX(iw, isolz1[ig], 0, isenz1[ig])];
                ac100 = aertab->model[im]->cms_all[INDEX(iw, isolz1[ig], 0, isenz2[ig])];
                ac001 = aertab->model[im]->cms_all[INDEX(iw, isolz2[ig], 0, isenz1[ig])];
                ac101 = aertab->model[im]->cms_all[INDEX(iw, isolz2[ig], 0, isenz2[ig])];

                ad000 = aertab->model[im]->dms_all[INDEX(iw, isolz1[ig], 0, isenz1[ig])];
                ad100 = aertab->model[im]->dms_all[INDEX(iw, isolz1[ig], 0, isenz2[ig])];
                ad001 = aertab->model[im]->dms_all[INDEX(iw, isolz2[ig], 0, isenz1[ig])];
                ad101 = aertab->model[im]->dms_all[INDEX(iw, isolz2[ig], 0, isenz2[ig])];

                ae000 = aertab->model[im]->ems_all[INDEX(iw, isolz1[ig], 0, isenz1[ig])];
                ae100 = aertab->model[im]->ems_all[INDEX(iw, isolz1[ig], 0, isenz2[ig])];
                ae001 = aertab->model[im]->ems_all[INDEX(iw, isolz2[ig], 0, isenz1[ig])];
                ae101 = aertab->model[im]->ems_all[INDEX(iw, isolz2[ig], 0, isenz2[ig])];

                a_coef[im][iw] = (1. - px)*(1. - rx) * as000 + px * rx * as101
                        + (1. - px) * rx * as001 + px * (1. - rx) * as100;

                b_coef[im][iw] = (1. - px)*(1. - rx) * ai000 + px * rx * ai101
                        + (1. - px) * rx * ai001 + px * (1. - qx)*(1. - rx) * ai100;

                c_coef[im][iw] = (1. - px)*(1. - rx) * ac000 + px * rx * ac101
                        + (1. - px) * rx * ac001 + px * (1. - qx)*(1. - rx) * ac100;

                d_coef[im][iw] = (1. - px)*(1. - rx) * ad000 + px * rx * ad101
                        + (1. - px) * rx * ad001 + px * (1. - qx)*(1. - rx) * ad100;

                e_coef[im][iw] = (1. - px)*(1. - rx) * ae000 + px * rx * ae101
                        + (1. - px) * rx * ae001 + px * (1. - qx)*(1. - rx) * ae100;
            } else {
                /*       printf("coeffs: as000,ai000,ac000,ad000,ae000\n");   */

            as000 = aertab->model[im]->ams_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz1[ig])];
            as100 = aertab->model[im]->ams_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz1[ig])];
            as010 = aertab->model[im]->ams_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz1[ig])];
            as110 = aertab->model[im]->ams_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz1[ig])];
            as001 = aertab->model[im]->ams_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz2[ig])];
            as011 = aertab->model[im]->ams_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz2[ig])];
            as101 = aertab->model[im]->ams_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz2[ig])];
            as111 = aertab->model[im]->ams_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz2[ig])];

            ai000 = aertab->model[im]->bms_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz1[ig])];
            ai100 = aertab->model[im]->bms_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz1[ig])];
            ai010 = aertab->model[im]->bms_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz1[ig])];
            ai110 = aertab->model[im]->bms_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz1[ig])];
            ai001 = aertab->model[im]->bms_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz2[ig])];
            ai011 = aertab->model[im]->bms_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz2[ig])];
            ai101 = aertab->model[im]->bms_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz2[ig])];
            ai111 = aertab->model[im]->bms_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz2[ig])];

            ac000 = aertab->model[im]->cms_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz1[ig])];
            ac100 = aertab->model[im]->cms_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz1[ig])];
            ac010 = aertab->model[im]->cms_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz1[ig])];
            ac110 = aertab->model[im]->cms_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz1[ig])];
            ac001 = aertab->model[im]->cms_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz2[ig])];
            ac011 = aertab->model[im]->cms_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz2[ig])];
            ac101 = aertab->model[im]->cms_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz2[ig])];
            ac111 = aertab->model[im]->cms_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz2[ig])];

            ad000 = aertab->model[im]->dms_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz1[ig])];
            ad100 = aertab->model[im]->dms_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz1[ig])];
            ad010 = aertab->model[im]->dms_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz1[ig])];
            ad110 = aertab->model[im]->dms_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz1[ig])];
            ad001 = aertab->model[im]->dms_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz2[ig])];
            ad011 = aertab->model[im]->dms_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz2[ig])];
            ad101 = aertab->model[im]->dms_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz2[ig])];
            ad111 = aertab->model[im]->dms_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz2[ig])];

            ae000 = aertab->model[im]->ems_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz1[ig])];
            ae100 = aertab->model[im]->ems_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz1[ig])];
            ae010 = aertab->model[im]->ems_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz1[ig])];
            ae110 = aertab->model[im]->ems_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz1[ig])];
            ae001 = aertab->model[im]->ems_all[INDEX(iw, isolz1[ig], iphi1[ig], isenz2[ig])];
            ae011 = aertab->model[im]->ems_all[INDEX(iw, isolz1[ig], iphi2[ig], isenz2[ig])];
            ae101 = aertab->model[im]->ems_all[INDEX(iw, isolz2[ig], iphi1[ig], isenz2[ig])];
            ae111 = aertab->model[im]->ems_all[INDEX(iw, isolz2[ig], iphi2[ig], isenz2[ig])];


           a_coef[im][iw] = (1. - rx)*(1. - qx)*(1. - px) * as000 + rx * qx * px * as111
                            + rx * (1. - qx) * px * as101 + (1. - rx) * qx * (1. - px) * as010
                            + rx * qx * (1. - px) * as110 + (1. - rx)*(1. - qx) * px * as001
                            + (1. - rx) * qx * px * as011 + rx * (1. - qx)*(1. - px) * as100;

           b_coef[im][iw] = (1. - rx)*(1. - qx)*(1. - px) * ai000 + rx * qx * px * ai111
                            + rx * (1. - qx) * px * ai101 + (1. - rx) * qx * (1. - px) * ai010
                            + rx * qx * (1. - px) * ai110 + (1. - rx)*(1. - qx) * px * ai001
                            + (1. - rx) * qx * px * ai011 + rx * (1. - qx)*(1. - px) * ai100;

           c_coef[im][iw] = (1. - rx)*(1. - qx)*(1. - px) * ac000 + rx * qx * px * ac111
                            + rx * (1. - qx) * px * ac101 + (1. - rx) * qx * (1. - px) * ac010
                            + rx * qx * (1. - px) * ac110 + (1. - rx)*(1. - qx) * px * ac001
                            + (1. - rx) * qx * px * ac011 + rx * (1. - qx)*(1. - px) * ac100;

           d_coef[im][iw] = (1. - rx)*(1. - qx)*(1. - px) * ad000 + rx * qx * px * ad111
                            + rx * (1. - qx) * px * ad101 + (1. - rx) * qx * (1. - px) * ad010
                            + rx * qx * (1. - px) * ad110 + (1. - rx)*(1. - qx) * px * ad001
                            + (1. - rx) * qx * px * ad011 + rx * (1. - qx)*(1. - px) * ad100;

           e_coef[im][iw] = (1. - rx)*(1. - qx)*(1. - px) * ae000 + rx * qx * px * ae111
                            + rx * (1. - qx) * px * ae101 + (1. - rx) * qx * (1. - px) * ae010
                            + rx * qx * (1. - px) * ae110 + (1. - rx)*(1. - qx) * px * ae001
                            + (1. - rx) * qx * px * ae011 + rx * (1. - qx)*(1. - px) * ae100;

            }
        }
    }

    // root finding is quadratic, but table includes cubic term, make sure it's zero
    if (fabs(d_coef[modnum][iwnir_l]) > 1e-9) {
        printf("non zero cubic term found in longest NIR wavelength of aerosol table. Zia!!\n");
        exit(1);
    }

    /* return pointers to coeffs for this geometry */
    *a = &a_coef[modnum][0];
    *b = &b_coef[modnum][0];
    *c = &c_coef[modnum][0];
    *d = &d_coef[modnum][0];
    *e = &e_coef[modnum][0];


    return;
}

/* Interpolate pc_rhoa using geometry           */
/* pc_rhoa:  ntau_870*npc                       */
/*  M.Zhang, Oct. 2022                          */

void get_pc_rhoa(int modnum, geom_str *geom, float **pc_rhoa)
 {
    static float lastsolz = -999.;
    static float lastsenz = -999.;
    static float lastphi = -999.;

    static int computed[MAXMODEL];

    static float ***pc_coef;    //MAXMODEL*ntau_870*npc
    static int ntau_870,npc;

    static float *p, *q, *r;
    static float as000, as100, as010, as110, as001, as011, as101, as111;

    static int *isolz1, *isolz2;
    static int *isenz1, *isenz2;
    static int *iphi1, *iphi2;

    static float *p_ar, *q_ar, *r_ar;
    static int *isolz1_ar, *isolz2_ar;
    static int *isenz1_ar, *isenz2_ar;
    static int *iphi1_ar, *iphi2_ar;
    static int gmult;

    float aphi;
    float px, qx, rx;
    int im, iw, i, ig, itau,ipc;
    static int firstCall = 1;

    if (firstCall == 1) {
        firstCall = 0;
        ntau_870=aertab->ntau_870;
        npc=aertab->npc;

        pc_coef=(float ***)calloc(MAXMODEL,sizeof(float **));

        for (i = 0; i < MAXMODEL; i++) {
            if ((pc_coef[i] = (float **) calloc(ntau_870, sizeof (float*))) == NULL) {
                printf("Unable to allocate space for pc_coef.\n");
                exit(1);
            }
            for(itau=0;itau<ntau_870;itau++){
                if ((pc_coef[i][itau] = (float *) calloc(aertab->npc, sizeof (float))) == NULL) {
                    printf("Unable to allocate space for pc_coef.\n");
                    exit(1);
                }
            }
        }
        /*  set up indicies, weights for band-dependent or nominal geometry */
        if ((geom->gmult == 0) || (interpol == 1)) {
            gmult = 0;
        } else {
            gmult = 1;
        }
        if (((p_ar = (float *) malloc(aertab->nwave * sizeof (float)))
                == NULL) ||
                ((q_ar = (float *) malloc(aertab->nwave * sizeof (float)))
                        == NULL) ||
                        ((r_ar = (float *) malloc(aertab->nwave * sizeof (float)))
                                == NULL)) {
            printf("Unable to allocate space for p, q, r weights.\n");
            exit(1);
        }
        if (((isolz1_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                == NULL) ||
                ((isenz1_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                        == NULL) ||
                        ((iphi1_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                                == NULL)) {
            printf("Unable to allocate space for interp indicies 1.\n");
            exit(1);
        }
        if (((isolz2_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                == NULL) ||
                ((isenz2_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                        == NULL) ||
                        ((iphi2_ar = (int *) malloc(aertab->nwave * sizeof (int)))
                                == NULL)) {
            printf("Unable to allocate space for interp indicies 2.\n");
            exit(1);
        }
        p = p_ar;
        q = q_ar;
        r = r_ar;
        isolz1 = isolz1_ar;
        isolz2 = isolz2_ar;
        isenz1 = isenz1_ar;
        isenz2 = isenz2_ar;
        iphi1 = iphi1_ar;
        iphi2 = iphi2_ar;
    }

    if (geom->solz[0] != lastsolz || geom->senz[0] != lastsenz ||
            geom->phi[0] != lastphi) {
        for (im = 0; im < aertab->nmodel; im++)
            computed[im] = 0;

        for (iw = 0; iw < aertab->nwave; iw++) {
            ig = iw * gmult;
            /* find bracketing solar indices */
            for (i = 0; i < aertab->nsolz; i++) {
                if (geom->solz[ig] < aertab->solz[i])
                    break;
            }
            isolz1[iw] = MAX(i - 1, 0);
            isolz2[iw] = MIN(i, aertab->nsolz - 1);
            if (isolz2[iw] != isolz1[iw])
                r[iw] = (geom->solz[ig] - aertab->solz[isolz1[iw]]) /
                (aertab->solz[isolz2[iw]] - aertab->solz[isolz1[iw]]);
            else
                r[iw] = 0.0;

            /* find bracketing view indices */
            for (i = 0; i < aertab->nsenz; i++) {
                if (geom->senz[ig] < aertab->senz[i])
                    break;
            }
            isenz1[iw] = MAX(i - 1, 0);
            isenz2[iw] = MIN(i, aertab->nsenz - 1);
            if (isenz2[iw] != isenz1[iw])
                p[iw] = (geom->senz[ig] - aertab->senz[isenz1[iw]]) /
                (aertab->senz[isenz2[iw]] - aertab->senz[isenz1[iw]]);
            else
                p[iw] = 0.0;

            /* find bracketing azimuth indices */
            aphi = fabs(geom->phi[ig]);
            for (i = 0; i < aertab->nphi; i++) {
                if (aphi < aertab->phi[i])
                    break;
            }
            iphi1[iw] = MAX(i - 1, 0);
            iphi2[iw] = MIN(i, aertab->nphi - 1);
            if (iphi2[iw] != iphi1[iw])
                q[iw] = (aphi - aertab->phi[iphi1[iw]]) /
                (aertab->phi[iphi2[iw]] - aertab->phi[iphi1[iw]]);
            else
                q[iw] = 0.0;
            if (gmult == 0) break;
        }

        /* remember last geometry */
        lastsolz = geom->solz[0];
        lastsenz = geom->senz[0];
        lastphi = geom->phi[0];

    }

    im = modnum;

    if (!computed[modnum]) {
        im = modnum;
        computed[modnum] = 1;

        iw=0;
        ig = iw * gmult;
        px = p[ig];
        qx = q[ig];
        rx = r[ig];
        if (isolz2[ig] == 0) {

            for(itau=0;itau<ntau_870;itau++)
                for(ipc=0;ipc<npc;ipc++){
                    as000 = aertab->model[im]->pc_rhoa[itau][isolz1[ig]][0][isenz1[ig]][ipc];
                    as100 = aertab->model[im]->pc_rhoa[itau][isolz1[ig]][0][isenz2[ig]][ipc];
                    as001 = aertab->model[im]->pc_rhoa[itau][isolz2[ig]][0][isenz1[ig]][ipc];
                    as101 = aertab->model[im]->pc_rhoa[itau][isolz2[ig]][0][isenz2[ig]][ipc];

                    pc_coef[im][itau][ipc] = (1. - px)*(1. - rx) * as000 + px * rx * as101
                            + (1. - px) * rx * as001 + px * (1. - rx) * as100;
                }
        } else {
            /*       printf("coeffs: as000,ai000,ac000,ad000,ae000\n");   */

            for(itau=0;itau<ntau_870;itau++)
                for(ipc=0;ipc<npc;ipc++){
                    as000 = aertab->model[im]->pc_rhoa[itau][isolz1[ig]][iphi1[ig]][isenz1[ig]][ipc];
                    as100 = aertab->model[im]->pc_rhoa[itau][isolz2[ig]][iphi1[ig]][isenz1[ig]][ipc];
                    as010 = aertab->model[im]->pc_rhoa[itau][isolz1[ig]][iphi2[ig]][isenz1[ig]][ipc];
                    as110 = aertab->model[im]->pc_rhoa[itau][isolz2[ig]][iphi2[ig]][isenz1[ig]][ipc];
                    as001 = aertab->model[im]->pc_rhoa[itau][isolz1[ig]][iphi1[ig]][isenz2[ig]][ipc];
                    as011 = aertab->model[im]->pc_rhoa[itau][isolz1[ig]][iphi2[ig]][isenz2[ig]][ipc];
                    as101 = aertab->model[im]->pc_rhoa[itau][isolz2[ig]][iphi1[ig]][isenz2[ig]][ipc];
                    as111 = aertab->model[im]->pc_rhoa[itau][isolz2[ig]][iphi2[ig]][isenz2[ig]][ipc];

                    pc_coef[im][itau][ipc] = (1. - rx)*(1. - qx)*(1. - px) * as000 + rx * qx * px * as111
                            + rx * (1. - qx) * px * as101 + (1. - rx) * qx * (1. - px) * as010
                            + rx * qx * (1. - px) * as110 + (1. - rx)*(1. - qx) * px * as001
                            + (1. - rx) * qx * px * as011 + rx * (1. - qx)*(1. - px) * as100;
                }
        }
    }

    /* return pointers to coeffs for this geometry */
  //  *pc_rhoa=&pc_coef[modnum][0][0];
    for(itau=0;itau<ntau_870;itau++)
        for(ipc=0;ipc<npc;ipc++){
            pc_rhoa[itau][ipc]=pc_coef[modnum][itau][ipc];
        }

    /*FILE *fp=fopen("/accounts/mzhang11/Rsdata/OCI/20220321/pc.txt","w");

    for(ipc=0;ipc<npc;ipc++){
       fprintf(fp,"%f %f\n",aertab->model[modnum]->pc_rhoa[0][0][0][0][ipc],aertab->model[modnum]->pc_components_rhoa[ipc][176]);
    }

    fclose(fp);*/

    return;
}
/* ---------------------------------------------------------------------------------------- */
/* model_select_ahmad() - select two aerosol models whose epsilon values bracket the        */
/*                        the observed ms epsilon, eps_obs                                  */
/*                                                                                          */
/* Z Ahmad July 2014.                                                                       */

/* ---------------------------------------------------------------------------------------- */
int comp_epsilonT(epsilonTstr *x, epsilonTstr *y) {
    return (x->eps_obs < y->eps_obs ? -1 : 1);
}

void model_select_ahmad(int32_t nmodels, int32_t *mindx, float eps_pred[], float eps_obs, int32_t *modmin,
        int32_t *modmax, float *modrat) {
    static epsilonTstr epsilonT[MAXAERMOD];

    int im, im1, im2;

    // set-up table to keep model epsilon and model number pairs together

    for (im = 0; im < nmodels; im++) {
        epsilonT[im].eps_obs = eps_pred[im];
        epsilonT[im].modnum = im;
        /*     printf("%d %7.4f %7.4f\n",im,eps_pred[im],eps_obs);              */
    }

    // sort into ascending model epsilon order

    qsort(epsilonT, nmodels, sizeof (epsilonTstr),
            (int (*)(const void *, const void *)) comp_epsilonT);

    // find bounding epsilon indices in table

    for (im = 0; im < nmodels; im++) {
        if (eps_obs < epsilonT[im].eps_obs)
            break;
    }

    //if(im==0) //no lower bounding by M. Zhang
    //{
    //	*modmin=-1;
    //	return;
    //}


    im1 = MAX(MIN(im - 1, nmodels - 1), 0);
    im2 = MAX(MIN(im, nmodels - 1), 0);


    // convert table indices to model indices of the input order
    // compute model weighting

    *modmin = epsilonT[im1].modnum;
    *modmax = epsilonT[im2].modnum;
    *modrat = (eps_obs - epsilonT[im1].eps_obs) / (epsilonT[im2].eps_obs - epsilonT[im1].eps_obs);
    /*    printf("%d %d %7.4f %7.4f %7.4f\n",im1,im2,eps_obs,epsilonT[im1].eps_obs,epsilonT[im2].eps_obs); */

    // If eps_obs is higer or lower than epsilon from table, then set the weight to 1
    if (*modmin == *modmax)
        *modrat = 1;

    return;
}


/*------------------------------------------------------------------------------------------ */
/* comp_rhoa_ms_eps() -  for a given model, computes the AOT and aerosol reflectance         */
/*                                                                                           */
/* Z ahmad July2014                                                                          */

/*------------------------------------------------------------------------------------------ */
int comp_rhoa_ms_eps(int32_t nwave, float wave[], geom_str *geom,
        float tau_iwnir_l, int32_t modl, float tau_pred[], float rho_pred[],float derv_rhoa[], float derv_taua_tmp[])
 {
    float *ac, *bc, *cc, *dc, *ec;
    float ext_modl[nwave];
    float lg_tau_pred[nwave];
    float lg_rho_pred[nwave];
    int iw, iwtab;

    /* get the coefficients for lg_rho vs lg_aot  */

    // Zia's function ---> something is wrong
    ms_eps_coef(modl, iwnir_l, wave, geom, &ac, &bc, &cc, &dc, &ec);
    // ms_eps_coef_cal(modl,nwave,solz,senz,phi,&ac,&bc,&cc,&dc,&ec);


    /* get the extinction coefficients and compute AOT at all wavelengths */

    for (iw = 0; iw < nwave; iw++) {
        iwtab = iwatab[iw];
        ext_modl[iw] = aertab->model[modl]->extc[iwtab];
    }

    /*   printf("tau_pred[iw],tau_iwnir_l\n");    */
    for (iw = 0; iw < nwave; iw++) {
        tau_pred[iw] = (ext_modl[iw] / ext_modl[iwnir_l]) * tau_iwnir_l;
        lg_tau_pred[iw] = log(tau_pred[iw]);

        if(derv_taua_tmp){

            derv_taua_tmp[iw] = ext_modl[iw] / ext_modl[iwnir_l];
            derv_rhoa[iw] = 1 / tau_pred[iw]
                    * (ext_modl[iw] / ext_modl[iwnir_l]);
        }
    }

    /* compute rho_pred */

    for (iw = 0; iw < nwave; iw++) {
        iwtab = iwatab[iw];
        lg_rho_pred[iw] = ac[iwtab] +
                bc[iwtab] * lg_tau_pred[iw] +
                cc[iwtab] * pow(lg_tau_pred[iw], 2) +
                dc[iwtab] * pow(lg_tau_pred[iw], 3) +
                ec[iwtab] * pow(lg_tau_pred[iw], 4);
        rho_pred[iw] = exp(lg_rho_pred[iw]);

        if(derv_taua_tmp){
            derv_rhoa[iw] *= (bc[iwtab] + 2 * cc[iwtab] * lg_tau_pred[iw]
                    + 3 * dc[iwtab] * pow(lg_tau_pred[iw], 2)
                    + 4 * ec[iwtab] * pow(lg_tau_pred[iw], 3));
            derv_rhoa[iw] *= rho_pred[iw];
        }
    }


    return (0);
}

/*------------------------------------------------------------------------------------------ */
/* comp_rhoa_pca() -  for a given model, computes the AOT and aerosol reflectance         */
/*                                                                                           */
/* M. Zhang,  Oct. 2022                                                                          */

/*------------------------------------------------------------------------------------------ */
int comp_rhoa_pca(int32_t nwave, float wave[], geom_str *geom,
        float tau_iwnir_l, int32_t modl, float tau_pred[], float rho_pred[],float derv_rhoa[], float derv_taua_tmp[])
 {
    static int firstcall=1;
    static float **pc_rhoa;
    static int npc, ntau_870;
    aermodstr *aermod=aertab->model[modl];

    int iw,itau,ipc, itau1,itau2;
    float rhoa1,rhoa2;

    if(firstcall){
        firstcall=0;
        ntau_870=aertab->ntau_870;
        npc=aertab->npc;
        pc_rhoa=(float **)malloc(aertab->ntau_870*sizeof(float *));

        for(itau=0;itau<ntau_870;itau++)
            pc_rhoa[itau]=(float *)malloc(aertab->npc*sizeof(float));
    }

  /*  if(tau_iwnir_l<aermod->tau_870[0]|| tau_iwnir_l>aermod->tau_870[ntau_870-1]){
        printf("tau_870 is out of the range in LUTs\n");
        return -1;
    }*/


    for(itau=0;itau<ntau_870;itau++){
        if(tau_iwnir_l<aermod->tau_870[itau])
            break;
    }
    itau1=MAX(itau-1,0);
    itau2=MIN(itau,ntau_870-1);

    if(itau1==0)
        itau2=1;
    if(itau2==ntau_870-1)
        itau1=itau2-1;

    get_pc_rhoa(modl,geom,pc_rhoa);

    for(iw=0;iw<nwave;iw++){
        rhoa1=0.;
        rhoa2=0.;

        for(ipc=0;ipc<npc;ipc++){
            rhoa1+=pc_rhoa[itau1][ipc]*aermod->pc_components_rhoa[ipc][iw];
            rhoa2+=pc_rhoa[itau2][ipc]*aermod->pc_components_rhoa[ipc][iw];
        }
        rhoa1=rhoa1+aermod->pc_mean_rhoa[iw];
        rhoa2=rhoa2+aermod->pc_mean_rhoa[iw];
        rho_pred[iw]=rhoa1+(rhoa2-rhoa1)*(tau_iwnir_l-aermod->tau_870[itau1])/(aermod->tau_870[itau2]-aermod->tau_870[itau1]);
        tau_pred[iw]=tau_iwnir_l*aermod->extc[iw]/aermod->extc[iwnir_l];

        if(derv_taua_tmp){
            derv_taua_tmp[iw]=aermod->extc[iw]/aermod->extc[iwnir_l];
            derv_rhoa    [iw]=(rhoa2-rhoa1)/(aermod->tau_870[itau2]-aermod->tau_870[itau1]);
        }
    }

    return (0);
}

/*------------------------------------------------------------------------------------------ */
/* comp_rhoa_ms_eps() -  for a given model, computes the AOT and aerosol reflectance         */
/*                                                                                           */
/* Z ahmad July2014                                                                          */
/* Modified by Amir Ibrahim January 2015 for linear space coefficients                       */

/*------------------------------------------------------------------------------------------ */
int comp_rhoa_ms_eps_lin(int32_t nwave, float wave[], geom_str *geom,
        float tau_iwnir_l, int32_t modl, float tau_pred[], float rho_pred[])
 {
    float *ac, *bc, *cc, *dc, *ec;
    float ext_modl[nwave];
    float lg_tau_pred[nwave];
    float lg_rho_pred[nwave];
    int iw, iwtab;

    /* get the coefficients for lg_rho vs lg_aot  */

    // Zia's function ---> something is wrong
    ms_eps_coef(modl, iwnir_l, wave, geom, &ac, &bc, &cc, &dc, &ec);
    // ms_eps_coef_cal(modl,nwave,geom,&ac,&bc,&cc,&dc,&ec);


    /* get the extinction coefficients and compute AOT at all wavelengths */

    for (iw = 0; iw < nwave; iw++) {
        iwtab = iwatab[iw];
        ext_modl[iw] = aertab->model[modl]->extc[iwtab];
    }

    /*   printf("tau_pred[iw],tau_iwnir_l\n");    */
    for (iw = 0; iw < nwave; iw++) {
        tau_pred[iw] = (ext_modl[iw] / ext_modl[iwnir_l]) * tau_iwnir_l;
        lg_tau_pred[iw] = (tau_pred[iw]);
    }

    /* compute rho_pred */

    for (iw = 0; iw < nwave; iw++) {
        iwtab = iwatab[iw];
        lg_rho_pred[iw] = ac[iwtab] +
                bc[iwtab] * lg_tau_pred[iw] +
                cc[iwtab] * pow(lg_tau_pred[iw], 2) +
                dc[iwtab] * pow(lg_tau_pred[iw], 3) +
                ec[iwtab] * pow(lg_tau_pred[iw], 4);
        rho_pred[iw] = (lg_rho_pred[iw]);
    }


    return (0);
}

/*------------------------------------------------------------------------------------------*/
/* ahmad_atm_corr_pca() - compute aerosol reflectance at all wavelengths using pca LUTs     */
/*                                                                                          */
/* M. Zhang, Oct. 2022                                                                       */

/*----------------------------------------------------------------------------------------  */

int ahmad_atm_corr_pca(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s, int32_t iwnir_l,
        int32_t nmodels, int32_t mindx[],
        geom_str *geom, float wv, float rhoa[],
        int32_t *modmin, int32_t *modmax, float *modrat, float *epsnir,
        float tau_pred_max[], float tau_pred_min[], float rho_pred_max[], float rho_pred_min[],
        float tau_aer[], float rho_aer[],int ip,uncertainty_t *uncertainty) {

    static int firstcall=1;
    static float **pc_rhoa;
    static int npc, ntau_870;
    aermodstr *aermod;
    float tau_iwnir_l[nmodels];
    float rho_iwnir_s_pred[nmodels];
    float eps_pred[nmodels];
    int im, modl;
    int im1, im2;
    float mwt;

    float derv_taua_rhoa[nmodels]; //derivative of taua[aer_l] to rha[aer_l]
    float derv_eps_rhoa_l[nmodels];
    float *derv_rhoa_min=NULL,*derv_rhoa_max=NULL;  //derivative of rhoa[]to taua[aer_l] for aermodmin and aermodmax
    float *derv_taua_min=NULL,*derv_taua_max=NULL;  //derivative of taua[]to taua[aer_l] for aermodmin and aermodmax
    float derv_mwt_rhoa_s=0., derv_mwt_rhoa_l=0.,derv_mwt_taua_l=0.,derv_mwt_rhow_l=0.;

    float eps_obs;//,deps_obs;

    float derv_eps_Lrc_s,derv_eps_Lrc_l,derv_eps_taua_l,derv_eps_rhow_l,derv_Lg_taua_l;
    float derv_eps_mod[4][nmodels]; //derivative of modeled eps to 0: rhorc_s, 1: rhorc_l, 2: taua_l and 3: t_sen*t_sol*rhow[aer_l]

    if(uncertainty){
        derv_eps_Lrc_s = uncertainty->derv_eps_Lrc_s;
        derv_eps_Lrc_l = uncertainty->derv_eps_Lrc_l;
        derv_eps_taua_l = uncertainty->derv_eps_taua_l;
        derv_eps_rhow_l = uncertainty->derv_eps_rhow_l;
        derv_Lg_taua_l = uncertainty->derv_Lg_taua[iwnir_l];

        derv_rhoa_min=(float *)malloc(nwave*sizeof(float));
        derv_rhoa_max=(float *)malloc(nwave*sizeof(float));
        derv_taua_min=(float *)malloc(nwave*sizeof(float));
        derv_taua_max=(float *)malloc(nwave*sizeof(float));
    }

    int iw,itau,ipc;

    int status = 0.0;
    static float *rhoa_modl_s,*rhoa_modl_l;

    if(firstcall){
        firstcall=0;
        ntau_870=aertab->ntau_870;
        npc=aertab->npc;
        rhoa_modl_s=(float *)malloc(aertab->ntau_870*sizeof(float));
        rhoa_modl_l=(float *)malloc(aertab->ntau_870*sizeof(float));
        pc_rhoa=(float **)malloc(aertab->ntau_870*sizeof(float *));

        for(itau=0;itau<aertab->ntau_870;itau++)
            pc_rhoa[itau]=(float *)malloc(aertab->npc*sizeof(float));
    }


    /* compute the observed epsilon */

    eps_obs = rhoa[iwnir_s] / rhoa[iwnir_l];

    /*             printf("rho_869,rho_748,eps_obs\n");                    */
    /*    printf("%10.5f %10.5f %10.5f\n",rhoa[iwnir_l],rhoa[iwnir_s],eps_obs);   */


    // compute MS epsilon for all nmodels.  note that nmodels may be only a subset of
    // the full model suite.  mindx maps from subset index to full index in suite


    for (im = 0; im < nmodels; im++) {

        modl = mindx[im]; // index in full model suite, needed for coefficient look-up
        aermod=aertab->model[modl];

        get_pc_rhoa(modl,geom,pc_rhoa);

        for(itau=0;itau<ntau_870;itau++){
            rhoa_modl_l[itau]=0;
            rhoa_modl_s[itau]=0;

            for(ipc=0;ipc<npc;ipc++){
                rhoa_modl_s[itau]+=pc_rhoa[itau][ipc]*aermod->pc_components_rhoa[ipc][iwnir_s];//
                rhoa_modl_l[itau]+=pc_rhoa[itau][ipc]*aermod->pc_components_rhoa[ipc][iwnir_l];
            }
            rhoa_modl_s[itau]+=aermod->pc_mean_rhoa[iwnir_s];
            rhoa_modl_l[itau]+=aermod->pc_mean_rhoa[iwnir_l];
        }

        /* compute model epsilon */
        if(rhoa[iwnir_l]<rhoa_modl_l[0])
            itau=1;
        else if (rhoa[iwnir_l]>rhoa_modl_l[ntau_870-1])
            itau=ntau_870-1;
        else{
            for(itau=1;itau<ntau_870;itau++){
                if(rhoa[iwnir_l]<rhoa_modl_l[itau] && rhoa[iwnir_l]>rhoa_modl_l[itau-1])
                    break;
            }
        }
        tau_iwnir_l[im]=aermod->tau_870[itau]+(aermod->tau_870[itau]-aermod->tau_870[itau-1])*(rhoa[iwnir_l]-rhoa_modl_l[itau])/(rhoa_modl_l[itau]-rhoa_modl_l[itau-1]);
        rho_iwnir_s_pred[im]=rhoa_modl_s[itau]+(rhoa_modl_s[itau]-rhoa_modl_s[itau-1])*(tau_iwnir_l[im]-aermod->tau_870[itau])/(aermod->tau_870[itau]-aermod->tau_870[itau-1]);

       // rho_iwnir_s_pred[im]=exp(rho_iwnir_s_pred[im]);
        eps_pred[im] = rho_iwnir_s_pred[im] / rhoa[iwnir_l];

        ///!!!! TBD: should include the tlw[nir] in the future  zhang
        if(uncertainty){
            derv_taua_rhoa[im]=(aermod->tau_870[itau]-aermod->tau_870[itau-1])/(rhoa_modl_l[itau]-rhoa_modl_l[itau-1]);
            derv_eps_rhoa_l[im]=derv_taua_rhoa[im]*((rhoa_modl_s[itau]-rhoa_modl_s[itau-1])/(aermod->tau_870[itau]-aermod->tau_870[itau-1]));

            derv_eps_rhoa_l[im]=derv_eps_rhoa_l[im]/rhoa[iwnir_l]-rho_iwnir_s_pred[im]/rhoa[iwnir_l]/rhoa[iwnir_l];
            derv_eps_mod[1][im]=derv_eps_rhoa_l[im];
            derv_eps_mod[2][im] = -derv_eps_rhoa_l[im] * uncertainty->derv_Lg_taua[iwnir_l];
            derv_eps_mod[3][im] = -derv_eps_rhoa_l[im];
        }
    }

    *epsnir = eps_obs;


    /* now, find the bounding models, but skip this if we only have one */
    /* model (as when vicariously calibrating)                          */

    if (nmodels > 1) {

        /* locate two model_epsilons that bracket the observed epsilon */
        /* this will return the model numbers for the subset of models */

        model_select_ahmad(nmodels, mindx, eps_pred, eps_obs, &im1, &im2, &mwt);

        if(uncertainty){
            if (im1 == im2) {
                for(iw=iwnir_s;iw<=iwnir_l;iw++)
                    uncertainty->derv_modrat_rhorc[iw-iwnir_s]=0.;
                uncertainty->derv_modrat_taua_l=0.;
            }
            else {
                derv_mwt_rhoa_s = 1 / (eps_pred[im2] - eps_pred[im1])
                        * derv_eps_Lrc_s;
                //derv_mwt_taua_s = 1 / (eps_pred[im2] - eps_pred[im1])
                //        * derv_eps_taua_s;

                derv_mwt_rhoa_l = 1 / (eps_pred[im2] - eps_pred[im1])
                        * derv_eps_Lrc_l;
                derv_mwt_rhoa_l -= ((mwt / (eps_pred[im2] - eps_pred[im1])
                        * derv_eps_mod[1][im2]));
                derv_mwt_rhoa_l += ((eps_obs - eps_pred[im2])
                        / pow(eps_pred[im2] - eps_pred[im1], 2)
                        * derv_eps_mod[1][im1]);

                derv_mwt_taua_l = 1 / (eps_pred[im2] - eps_pred[im1])
                        * derv_eps_taua_l;
                derv_mwt_taua_l -= ((mwt / (eps_pred[im2] - eps_pred[im1])
                        * derv_eps_mod[2][im2]));
                derv_mwt_taua_l += ((eps_obs - eps_pred[im2])
                        / pow(eps_pred[im2] - eps_pred[im1], 2)
                        * derv_eps_mod[2][im1]);

                derv_mwt_rhow_l = 1 / (eps_pred[im2] - eps_pred[im1])
                        * derv_eps_rhow_l;
                derv_mwt_rhow_l -= ((mwt / (eps_pred[im2] - eps_pred[im1])
                        * derv_eps_mod[3][im2]));
                derv_mwt_rhow_l += ((eps_obs - eps_pred[im2])
                        / pow(eps_pred[im2] - eps_pred[im1], 2)
                        * derv_eps_mod[3][im1]);

                uncertainty->derv_modrat_rhorc[0] = derv_mwt_rhoa_s;
                uncertainty->derv_modrat_rhorc[iwnir_l-iwnir_s] = derv_mwt_rhoa_l;
                uncertainty->derv_modrat_taua_l = derv_mwt_taua_l;
                uncertainty->derv_modrat_rhow_l = derv_mwt_rhow_l;
            }
        }
    } else {

        /* single model case (allows fixed model by limiting model suite) */

        im1 = 0;
        im2 = 0;
        mwt = 0.0;
    }

    /* map selection to full suite for subsequent processing and return */

    *modmin = mindx[im1];
    *modmax = mindx[im2];
    *modrat = mwt;

    /* compute tau_pred and rho_predicted */

    if(comp_rhoa_pca(nwave, wave, geom, tau_iwnir_l[im1], *modmin, tau_pred_min, rho_pred_min,derv_rhoa_min,derv_taua_min) !=0)
        return -1;
    if(comp_rhoa_pca(nwave, wave, geom, tau_iwnir_l[im2], *modmax, tau_pred_max, rho_pred_max,derv_rhoa_max,derv_taua_max) !=0)
        return -1;

    /* compute weighted tau_aer and rho_aer */

    for (iw = 0; iw < nwave; iw++) {
        tau_aer[iw] = (1.0 - mwt) * tau_pred_min[iw] + mwt * tau_pred_max[iw];
        rho_aer[iw] = (1.0 - mwt) * rho_pred_min[iw] + mwt * rho_pred_max[iw];


        if(uncertainty){
            uncertainty->derv_La_rhorc[iw][iwnir_l-iwnir_s] = (1.0 - mwt) * derv_rhoa_min[iw]
                    * derv_taua_rhoa[im1]
                    + mwt * derv_rhoa_max[iw] * derv_taua_rhoa[im2]
                    + (rho_pred_max[iw] - rho_pred_min[iw]) * derv_mwt_rhoa_l;
            uncertainty->derv_La_taua_l[iw] = -(1.0 - mwt) * derv_rhoa_min[iw]
                    * derv_taua_rhoa[im1] * derv_Lg_taua_l
                    - mwt * derv_rhoa_max[iw] * derv_taua_rhoa[im2]
                            * derv_Lg_taua_l
                    + (rho_pred_max[iw] - rho_pred_min[iw]) * derv_mwt_taua_l;
            uncertainty->derv_La_rhorc[iw][0] = (rho_pred_max[iw] - rho_pred_min[iw])
                    * derv_mwt_rhoa_s;
            uncertainty->derv_La_rhow_l[iw] = -(1.0 - mwt) * derv_rhoa_min[iw]
                    * derv_taua_rhoa[im1]
                    - mwt * derv_rhoa_max[iw] * derv_taua_rhoa[im2]
                    + (rho_pred_max[iw] - rho_pred_min[iw]) * derv_mwt_rhow_l;
           // uncertainty->derv_taua_s[iw] = (rho_pred_max[iw] - rho_pred_min[iw])
            //        * derv_mwt_taua_s;

            uncertainty->derv_taua_rhorc[iw][iwnir_l-iwnir_s]= (1.0 - mwt) * derv_taua_min[iw]
                    * derv_taua_rhoa[im1]
                    + mwt * derv_taua_max[iw] * derv_taua_rhoa[im2]
                    + (tau_pred_max[iw] - tau_pred_min[iw]) * derv_mwt_rhoa_l;
            uncertainty->derv_taua_taua_l[iw] = -(1.0 - mwt) * derv_taua_min[iw]
                    * derv_taua_rhoa[im1] * derv_Lg_taua_l
                    - mwt * derv_taua_max[iw] * derv_taua_rhoa[im2]
                            * derv_Lg_taua_l
                    + (tau_pred_max[iw] - tau_pred_min[iw]) * derv_mwt_taua_l;
            uncertainty->derv_taua_rhorc[iw][0] = (tau_pred_max[iw] - tau_pred_min[iw])
                    * derv_mwt_rhoa_s;
            uncertainty->derv_taua_rhow_l[iw] = -(1.0 - mwt) * derv_taua_min[iw]
                    * derv_taua_rhoa[im1]
                    - mwt * derv_taua_max[iw] * derv_taua_rhoa[im2]
                    + (tau_pred_max[iw] - tau_pred_min[iw]) * derv_mwt_rhow_l;
            //uncertainty->derv_taua_taua_s[iw] = (tau_pred_max[iw] - tau_pred_min[iw])
            //        * derv_mwt_taua_s;

           // derv_taua[0][iw][0] = 0.;
            uncertainty->derv_taua_min_rhorc_l[iw] = derv_taua_min[iw] * derv_taua_rhoa[im1];
            uncertainty->derv_taua_min_taua_l[iw] = -derv_taua_min[iw] * derv_taua_rhoa[im1]
                    * derv_Lg_taua_l;
            uncertainty->derv_taua_min_rhow_l[iw] = -derv_taua_min[iw] * derv_taua_rhoa[im1];
           // derv_taua[0][iw][4] = 0.;

           // derv_taua[1][iw][0] = 0.;
            uncertainty->derv_taua_max_rhorc_l[iw]= derv_taua_max[iw] * derv_taua_rhoa[im2];
            uncertainty->derv_taua_max_taua_l [iw]= -derv_taua_max[iw] * derv_taua_rhoa[im2]
                    * derv_Lg_taua_l;
            uncertainty->derv_taua_max_rhow_l [iw]= -derv_taua_max[iw] * derv_taua_rhoa[im2];
            //derv_taua[1][iw][4] = 0.;
        }
    }

    if(uncertainty){
        free(derv_rhoa_min);
        free(derv_rhoa_max);
        free(derv_taua_min);
        free(derv_taua_max);
    }

    return (status);
}
/*------------------------------------------------------------------------------------------*/
/* ahmad_atm_corr() - compute aerosol reflectance at all wavelengths                        */
/*                                                                                          */
/* Z Ahmad. July 2014                                                                       */

/*----------------------------------------------------------------------------------------  */
int ahmad_atm_corr(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s, int32_t iwnir_l,
        int32_t nmodels, int32_t mindx[],
        geom_str *geom, float wv, float rhoa[],
        int32_t *modmin, int32_t *modmax, float *modrat, float *epsnir,
        float tau_pred_max[], float tau_pred_min[], float rho_pred_max[], float rho_pred_min[],
        float tau_aer[], float rho_aer[],int ip,uncertainty_t *uncertainty) {


    int status = 0.0;
    if(use_pca_lut){
        status=ahmad_atm_corr_pca(sensorID,wave,nwave,iwnir_s,iwnir_l,nmodels,mindx,geom,wv,rhoa,
                modmin,modmax,modrat,epsnir,tau_pred_max,tau_pred_min,rho_pred_max,rho_pred_min,tau_aer,
                rho_aer,ip,uncertainty);
        return status;
    }
    float *ac, *bc, *cc, *dc, *ec;
    float ext_iwnir_l, ext_iwnir_s;
    double ax, bx, cx, fx;
    float tau_iwnir_l[nmodels];
    float lg_tau_iwnir_s[nmodels], tau_iwnir_s[nmodels];
    float lg_rho_iwnir_s_pred[nmodels], rho_iwnir_s_pred[nmodels];
    float eps_pred[nmodels];
    int im, modl;
    int im1, im2;
    float mwt;
    float lg_tau_iwnir_l;
    int iwtab_l, iwtab_s;

    float derv_taua_rhoa[nmodels]; //derivative of taua[aer_l] to rhoa[aer_l]
    float derv_eps_rhoa_l[nmodels];
    float *derv_rhoa_min=NULL,*derv_rhoa_max=NULL;  //derivative of rhoa[]to taua[aer_l] for aermodmin and aermodmax
    float *derv_taua_min=NULL,*derv_taua_max=NULL;  //derivative of taua[]to taua[aer_l] for aermodmin and aermodmax
    float derv_mwt_rhoa_s=0., derv_mwt_rhoa_l=0.,derv_mwt_taua_l=0.,derv_mwt_rhow_l=0.;

    float eps_obs;//,deps_obs;

    float derv_eps_Lrc_s,derv_eps_Lrc_l,derv_eps_taua_l,derv_eps_rhow_l,derv_Lg_taua_l;
    float derv_eps_mod[4][nmodels]; //derivative of modeled eps to 0: rhorc_s, 1: rhorc_l, 2: taua_l and 3: t_sen*t_sol*rhow[aer_l]

    if(uncertainty){
        derv_eps_Lrc_s = uncertainty->derv_eps_Lrc_s;
        derv_eps_Lrc_l = uncertainty->derv_eps_Lrc_l;
        derv_eps_taua_l = uncertainty->derv_eps_taua_l;
        derv_eps_rhow_l = uncertainty->derv_eps_rhow_l;
        derv_Lg_taua_l = uncertainty->derv_Lg_taua[iwnir_l];

        derv_rhoa_min=(float *)malloc(nwave*sizeof(float));
        derv_rhoa_max=(float *)malloc(nwave*sizeof(float));
        derv_taua_min=(float *)malloc(nwave*sizeof(float));
        derv_taua_max=(float *)malloc(nwave*sizeof(float));
    }

    int iw;

    /* compute the observed epsilon */

    eps_obs = rhoa[iwnir_s] / rhoa[iwnir_l];

    /*             printf("rho_869,rho_748,eps_obs\n");                    */
    /*    printf("%10.5f %10.5f %10.5f\n",rhoa[iwnir_l],rhoa[iwnir_s],eps_obs);   */


    // compute MS epsilon for all nmodels.  note that nmodels may be only a subset of
    // the full model suite.  mindx maps from subset index to full index in suite


    for (im = 0; im < nmodels; im++) {

        modl = mindx[im]; // index in full model suite, needed for coefficient look-up

        /* compute AOT at longest aerosol wavelength (iwnir_l) */

        // Zia's function ---> something is wrong
        //ms_eps_coef(modl,iwnir_l,wave,solz,senz,phi,&ac,&bc,&cc,&dc,&ec);

        //        ms_eps_coef_cal(modl,nwave,solz,senz,phi,&ac,&bc,&cc,&dc,&ec);
        ms_eps_coef(modl, iwnir_l, wave, geom, &ac, &bc, &cc, &dc, &ec);
        iwtab_l = iwatab[iwnir_l];
        iwtab_s = iwatab[iwnir_s];

        ax = (double) ac[iwtab_l] - log((double) rhoa[iwnir_l]);
        bx = (double) bc[iwtab_l];
        cx = (double) cc[iwtab_l];

        fx = bx * bx - 4.0 * ax*cx;
        if (fx > 0.0 && cx != 0.0) {
            lg_tau_iwnir_l = 0.5 * (-bx + sqrt(fx)) / cx;
            if(uncertainty)
                derv_taua_rhoa[im] = 1 / rhoa[iwnir_l] / sqrt(fx);

            tau_iwnir_l[im] = exp(lg_tau_iwnir_l);
            if(uncertainty)
                derv_taua_rhoa[im]*=tau_iwnir_l[im];
        } else {
            status = 1;
            break;    // TODO: need to think about it
        }

        /* compute AOT at shortest aerosol wavelength (iwnir_s) */

        ext_iwnir_l = aertab->model[modl]->extc[iwtab_l];
        ext_iwnir_s = aertab->model[modl]->extc[iwtab_s];

        tau_iwnir_s[im] = (ext_iwnir_s / ext_iwnir_l) * tau_iwnir_l[im];
        lg_tau_iwnir_s[im] = log(tau_iwnir_s[im]);

        if(uncertainty)
            derv_eps_rhoa_l[im]=1/tau_iwnir_s[im]*(ext_iwnir_s / ext_iwnir_l)*derv_taua_rhoa[im];


        /* compute reflectance at (iwnir_s) */

        lg_rho_iwnir_s_pred[im] = ac[iwtab_s] + bc[iwtab_s] * lg_tau_iwnir_s[im] +
                cc[iwtab_s] * pow(lg_tau_iwnir_s[im], 2) +
                dc[iwtab_s] * pow(lg_tau_iwnir_s[im], 3) +
                ec[iwtab_s] * pow(lg_tau_iwnir_s[im], 4);

        if(uncertainty)
            derv_eps_rhoa_l[im]=( bc[iwtab_s]+ 2*cc[iwtab_s]*lg_tau_iwnir_s[im]+3*dc[iwtab_s]*pow(lg_tau_iwnir_s[im], 2)+4*ec[iwtab_s] * pow(lg_tau_iwnir_s[im], 3) )*derv_eps_rhoa_l[im];

        rho_iwnir_s_pred[im] = exp(lg_rho_iwnir_s_pred[im]);

        if(uncertainty)
            derv_eps_rhoa_l[im]=rho_iwnir_s_pred[im]*derv_eps_rhoa_l[im];

        /* compute model epsilon */

        eps_pred[im] = rho_iwnir_s_pred[im] / rhoa[iwnir_l];

        ///!!!! TBD: should include the tlw[nir] in the future  zhang
        if(uncertainty){
            derv_eps_rhoa_l[im]=derv_eps_rhoa_l[im]/rhoa[iwnir_l]-rho_iwnir_s_pred[im]/rhoa[iwnir_l]/rhoa[iwnir_l];
            derv_eps_mod[1][im]=derv_eps_rhoa_l[im];
            derv_eps_mod[2][im] = -derv_eps_rhoa_l[im] * uncertainty->derv_Lg_taua[iwnir_l];
            derv_eps_mod[3][im] = -derv_eps_rhoa_l[im];
        }
    }


    *epsnir = eps_obs;


    /* now, find the bounding models, but skip this if we only have one */
    /* model (as when vicariously calibrating)                          */

    if (nmodels > 1) {

        /* locate two model_epsilons that bracket the observed epsilon */
        /* this will return the model numbers for the subset of models */

        model_select_ahmad(nmodels, mindx, eps_pred, eps_obs, &im1, &im2, &mwt);

        if(uncertainty){
            if (im1 == im2) {
                for(iw=iwnir_s;iw<=iwnir_l;iw++)
                    uncertainty->derv_modrat_rhorc[iw-iwnir_s]=0.;
                uncertainty->derv_modrat_taua_l=0.;
            }
            else {
                derv_mwt_rhoa_s = 1 / (eps_pred[im2] - eps_pred[im1])
                        * derv_eps_Lrc_s;
                //derv_mwt_taua_s = 1 / (eps_pred[im2] - eps_pred[im1])
                //        * derv_eps_taua_s;

                derv_mwt_rhoa_l = 1 / (eps_pred[im2] - eps_pred[im1])
                        * derv_eps_Lrc_l;
                derv_mwt_rhoa_l -= ((mwt / (eps_pred[im2] - eps_pred[im1])
                        * derv_eps_mod[1][im2]));
                derv_mwt_rhoa_l += ((eps_obs - eps_pred[im2])
                        / pow(eps_pred[im2] - eps_pred[im1], 2)
                        * derv_eps_mod[1][im1]);

                derv_mwt_taua_l = 1 / (eps_pred[im2] - eps_pred[im1])
                        * derv_eps_taua_l;
                derv_mwt_taua_l -= ((mwt / (eps_pred[im2] - eps_pred[im1])
                        * derv_eps_mod[2][im2]));
                derv_mwt_taua_l += ((eps_obs - eps_pred[im2])
                        / pow(eps_pred[im2] - eps_pred[im1], 2)
                        * derv_eps_mod[2][im1]);

                derv_mwt_rhow_l = 1 / (eps_pred[im2] - eps_pred[im1])
                        * derv_eps_rhow_l;
                derv_mwt_rhow_l -= ((mwt / (eps_pred[im2] - eps_pred[im1])
                        * derv_eps_mod[3][im2]));
                derv_mwt_rhow_l += ((eps_obs - eps_pred[im2])
                        / pow(eps_pred[im2] - eps_pred[im1], 2)
                        * derv_eps_mod[3][im1]);

                uncertainty->derv_modrat_rhorc[0] = derv_mwt_rhoa_s;
                uncertainty->derv_modrat_rhorc[iwnir_l-iwnir_s] = derv_mwt_rhoa_l;
                uncertainty->derv_modrat_taua_l = derv_mwt_taua_l;
                uncertainty->derv_modrat_rhow_l = derv_mwt_rhow_l;
            }
        }
    } else {

        /* single model case (allows fixed model by limiting model suite) */

        im1 = 0;
        im2 = 0;
        mwt = 0.0;
    }

    /* map selection to full suite for subsequent processing and return */

    *modmin = mindx[im1];
    *modmax = mindx[im2];
    *modrat = mwt;

    /* compute tau_pred and rho_predicted */

    comp_rhoa_ms_eps(nwave, wave, geom, tau_iwnir_l[im1], *modmin, tau_pred_min, rho_pred_min,derv_rhoa_min,derv_taua_min);
    comp_rhoa_ms_eps(nwave, wave, geom, tau_iwnir_l[im2], *modmax, tau_pred_max, rho_pred_max,derv_rhoa_max,derv_taua_max);

    /* compute weighted tau_aer and rho_aer */

    for (iw = 0; iw < nwave; iw++) {
        tau_aer[iw] = (1.0 - mwt) * tau_pred_min[iw] + mwt * tau_pred_max[iw];
        rho_aer[iw] = (1.0 - mwt) * rho_pred_min[iw] + mwt * rho_pred_max[iw];


        if(uncertainty){
            uncertainty->derv_La_rhorc[iw][iwnir_l-iwnir_s] = (1.0 - mwt) * derv_rhoa_min[iw]
                    * derv_taua_rhoa[im1]
                    + mwt * derv_rhoa_max[iw] * derv_taua_rhoa[im2]
                    + (rho_pred_max[iw] - rho_pred_min[iw]) * derv_mwt_rhoa_l;
            uncertainty->derv_La_taua_l[iw] = -(1.0 - mwt) * derv_rhoa_min[iw]
                    * derv_taua_rhoa[im1] * derv_Lg_taua_l
                    - mwt * derv_rhoa_max[iw] * derv_taua_rhoa[im2]
                            * derv_Lg_taua_l
                    + (rho_pred_max[iw] - rho_pred_min[iw]) * derv_mwt_taua_l;
            uncertainty->derv_La_rhorc[iw][0] = (rho_pred_max[iw] - rho_pred_min[iw])
                    * derv_mwt_rhoa_s;
            uncertainty->derv_La_rhow_l[iw] = -(1.0 - mwt) * derv_rhoa_min[iw]
                    * derv_taua_rhoa[im1]
                    - mwt * derv_rhoa_max[iw] * derv_taua_rhoa[im2]
                    + (rho_pred_max[iw] - rho_pred_min[iw]) * derv_mwt_rhow_l;
           // uncertainty->derv_taua_s[iw] = (rho_pred_max[iw] - rho_pred_min[iw])
            //        * derv_mwt_taua_s;

            uncertainty->derv_taua_rhorc[iw][iwnir_l-iwnir_s]= (1.0 - mwt) * derv_taua_min[iw]
                    * derv_taua_rhoa[im1]
                    + mwt * derv_taua_max[iw] * derv_taua_rhoa[im2]
                    + (tau_pred_max[iw] - tau_pred_min[iw]) * derv_mwt_rhoa_l;
            uncertainty->derv_taua_taua_l[iw] = -(1.0 - mwt) * derv_taua_min[iw]
                    * derv_taua_rhoa[im1] * derv_Lg_taua_l
                    - mwt * derv_taua_max[iw] * derv_taua_rhoa[im2]
                            * derv_Lg_taua_l
                    + (tau_pred_max[iw] - tau_pred_min[iw]) * derv_mwt_taua_l;
            uncertainty->derv_taua_rhorc[iw][0] = (tau_pred_max[iw] - tau_pred_min[iw])
                    * derv_mwt_rhoa_s;
            uncertainty->derv_taua_rhow_l[iw] = -(1.0 - mwt) * derv_taua_min[iw]
                    * derv_taua_rhoa[im1]
                    - mwt * derv_taua_max[iw] * derv_taua_rhoa[im2]
                    + (tau_pred_max[iw] - tau_pred_min[iw]) * derv_mwt_rhow_l;
            //uncertainty->derv_taua_taua_s[iw] = (tau_pred_max[iw] - tau_pred_min[iw])
            //        * derv_mwt_taua_s;

           // derv_taua[0][iw][0] = 0.;
            uncertainty->derv_taua_min_rhorc_l[iw] = derv_taua_min[iw] * derv_taua_rhoa[im1];
            uncertainty->derv_taua_min_taua_l[iw] = -derv_taua_min[iw] * derv_taua_rhoa[im1]
                    * derv_Lg_taua_l;
            uncertainty->derv_taua_min_rhow_l[iw] = -derv_taua_min[iw] * derv_taua_rhoa[im1];
           // derv_taua[0][iw][4] = 0.;

           // derv_taua[1][iw][0] = 0.;
            uncertainty->derv_taua_max_rhorc_l[iw]= derv_taua_max[iw] * derv_taua_rhoa[im2];
            uncertainty->derv_taua_max_taua_l [iw]= -derv_taua_max[iw] * derv_taua_rhoa[im2]
                    * derv_Lg_taua_l;
            uncertainty->derv_taua_max_rhow_l [iw]= -derv_taua_max[iw] * derv_taua_rhoa[im2];
            //derv_taua[1][iw][4] = 0.;
        }
    }

    if(uncertainty){
        free(derv_rhoa_min);
        free(derv_rhoa_max);
        free(derv_taua_min);
        free(derv_taua_max);
    }

    return (status);
}

/*------------------------------------------------------------------------------------------*/
/* ahmad_atm_corr_lin() - compute aerosol reflectance at all wavelengths in linear space    */
/*                                                                                          */
/* Z Ahmad. July 2014                                                                       */
/* Modified to linear spaced coefficient by Amir Ibrahim January 2017                       */

/*----------------------------------------------------------------------------------------  */
int ahmad_atm_corr_lin(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s, int32_t iwnir_l,
        int32_t nmodels, int32_t mindx[],
        geom_str *geom, float wv, float rhoa[],
        int32_t *modmin, int32_t *modmax, float *modrat, float *epsnir,
        float tau_pred_max[], float tau_pred_min[], float rho_pred_max[], float rho_pred_min[],
        float tau_aer[], float rho_aer[]) {


    float *ac, *bc, *cc, *dc, *ec;
    float ext_iwnir_l, ext_iwnir_s;
    double ax, bx, cx, fx;
    float tau_iwnir_l[nmodels];
    float lg_tau_iwnir_s[nmodels], tau_iwnir_s[nmodels];
    float lg_rho_iwnir_s_pred[nmodels], rho_iwnir_s_pred[nmodels];
    float eps_pred[nmodels];
    int im, modl;
    int im1, im2;
    float mwt;
    float lg_tau_iwnir_l;
    int iwtab_l, iwtab_s;


    static float eps_obs;

    int iw;

    int status = 0.0;


    /* compute the observed epsilon */

    eps_obs = rhoa[iwnir_s] / rhoa[iwnir_l];

    /*             printf("rho_869,rho_748,eps_obs\n");                    */
    /*    printf("%10.5f %10.5f %10.5f\n",rhoa[iwnir_l],rhoa[iwnir_s],eps_obs);   */


    // compute MS epsilon for all nmodels.  note that nmodels may be only a subset of
    // the full model suite.  mindx maps from subset index to full index in suite


    for (im = 0; im < nmodels; im++) {

        modl = mindx[im]; // index in full model suite, needed for coefficient look-up

        /* compute AOT at longest aerosol wavelength (iwnir_l) */

        // Zia's function ---> something is wrong
        //ms_eps_coef(modl,iwnir_l,wave,solz,senz,phi,&ac,&bc,&cc,&dc,&ec);

        //        ms_eps_coef_cal(modl,nwave,solz,senz,phi,&ac,&bc,&cc,&dc,&ec);
        ms_eps_coef(modl, iwnir_l, wave, geom, &ac, &bc, &cc, &dc, &ec);
        iwtab_l = iwatab[iwnir_l];
        iwtab_s = iwatab[iwnir_s];

        ax = (double) ac[iwtab_l]-(double) rhoa[iwnir_l];
        bx = (double) bc[iwtab_l];
        cx = (double) cc[iwtab_l];

        fx = bx * bx - 4.0 * ax*cx;
        if (fx > 0.0 && cx != 0.0) {
            lg_tau_iwnir_l = 0.5 * (-bx + sqrt(fx)) / cx;
            tau_iwnir_l[im] = (lg_tau_iwnir_l);
        } else {
            status = 1;
            break;
        }

        /* compute AOT at shortest aerosol wavelength (iwnir_s) */

        ext_iwnir_l = aertab->model[modl]->extc[iwtab_l];
        ext_iwnir_s = aertab->model[modl]->extc[iwtab_s];

        tau_iwnir_s[im] = (ext_iwnir_s / ext_iwnir_l) * tau_iwnir_l[im];
        lg_tau_iwnir_s[im] = (tau_iwnir_s[im]);

        if (ax > 1e5) {
            printf("\nErroneous aerosol option, %d\n", aer_opt);
            printf("You are using a log-space LUT. The aerosol LUT coefficients need to be in linear-space. Use aer_opt=-18 instead. \n");
            exit(FATAL_ERROR);
        }


        /* compute reflectance at (iwnir_s) */

        lg_rho_iwnir_s_pred[im] = ac[iwtab_s] + bc[iwtab_s] * lg_tau_iwnir_s[im] +
                cc[iwtab_s] * pow(lg_tau_iwnir_s[im], 2) +
                dc[iwtab_s] * pow(lg_tau_iwnir_s[im], 3) +
                ec[iwtab_s] * pow(lg_tau_iwnir_s[im], 4);

        rho_iwnir_s_pred[im] = (lg_rho_iwnir_s_pred[im]);

        /* compute model epsilon */

        eps_pred[im] = rho_iwnir_s_pred[im] / rhoa[iwnir_l];

    }


    *epsnir = eps_obs;


    /* now, find the bounding models, but skip this if we only have one */
    /* model (as when vicariously calibrating)                          */

    if (nmodels > 1) {

        /* locate two model_epsilons that bracket the observed epsilon */
        /* this will return the model numbers for the subset of models */

        model_select_ahmad(nmodels, mindx, eps_pred, eps_obs, &im1, &im2, &mwt);

    } else {

        /* single model case (allows fixed model by limiting model suite) */

        im1 = 0;
        im2 = 0;
        mwt = 0.0;
    }

    /* map selection to full suite for subsequent processing and return */

    *modmin = mindx[im1];
    *modmax = mindx[im2];
    *modrat = mwt;

    /* compute tau_pred and rho_predicted */

    comp_rhoa_ms_eps_lin(nwave, wave, geom, tau_iwnir_l[im1], *modmin, tau_pred_min, rho_pred_min);
    comp_rhoa_ms_eps_lin(nwave, wave, geom, tau_iwnir_l[im2], *modmax, tau_pred_max, rho_pred_max);

    /* compute weighted tau_aer and rho_aer */

    for (iw = 0; iw < nwave; iw++) {
        tau_aer[iw] = (1.0 - mwt) * tau_pred_min[iw] + mwt * tau_pred_max[iw];
        rho_aer[iw] = (1.0 - mwt) * rho_pred_min[iw] + mwt * rho_pred_max[iw];
    }


    return (status);
}

/* ---------------------------------------------------------------------------------------- */
/* model_phase() - return phase function at model wavelengths at the input geometry.        */
/*                                                                                          */
/* This is effectively a C version of load_ss.f by M. Wang.  The program optimizes for      */
/* multiple calls at the same geometry by computing for all models on the first call with a */
/* new geometry.  It returns a pointer to the internal static array for the requested model.*/
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */
/*  W. Robinson, SAIC modified to use band-dependent viewing geometry, if avail             */

/* ---------------------------------------------------------------------------------------- */
float *model_phase(int modnum, geom_str *geom) {
    static float nw = 1.334;

    static int computed[MAXMODEL];
    static float lastsolz = -999.;
    static float lastsenz = -999.;
    static float lastphi = -999.;
    static float *phase[MAXMODEL];
    static float *fres1, *fres2;
    static float *scatt1, *scatt2;
    static float scatt1_cnst, scatt2_cnst, fres1_cnst, fres2_cnst;
    static float *scatt1_ar, *scatt2_ar, *fres1_ar, *fres2_ar;
    static int firstCall = 1, gmult = 1;

    float phase1, phase2;
    int im, iw, ig;

    if (firstCall == 1) {
        firstCall = 0;
        for (im = 0; im < MAXMODEL; im++) {
            if ((phase[im] = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for phase.\n");
                exit(1);
            }
        }
        /* set up geometry based scattering and fresnel */
        if ((geom->gmult == 0) || (interpol == 1)) {
            gmult = 0;
            scatt1 = &scatt1_cnst;
            scatt2 = &scatt2_cnst;
            fres1 = &fres1_cnst;
            fres2 = &fres2_cnst;
        } else {
            gmult = 1;
            if ((scatt1_ar =
                    (float *) malloc(aertab->nwave * sizeof (float))) == NULL) {
                printf("Unable to allocate space for scatt1.\n");
                exit(1);
            }
            if ((scatt2_ar =
                    (float *) malloc(aertab->nwave * sizeof (float))) == NULL) {
                printf("Unable to allocate space for scatt2.\n");
                exit(1);
            }
            if ((fres1_ar =
                    (float *) malloc(aertab->nwave * sizeof (float))) == NULL) {
                printf("Unable to allocate space for fres1.\n");
                exit(1);
            }
            if ((fres2_ar =
                    (float *) malloc(aertab->nwave * sizeof (float))) == NULL) {
                printf("Unable to allocate space for fres2.\n");
                exit(1);
            }
            scatt1 = scatt1_ar;
            scatt2 = scatt2_ar;
            fres1 = fres1_ar;
            fres2 = fres2_ar;
        }
    }

    /* recalculate only if geometry changes */

    if ((geom->solz[0] != lastsolz) || (geom->senz[0] != lastsenz) ||
            (geom->phi[0] != lastphi)) {

        float temp;

        for (im = 0; im < aertab->nmodel; im++)
            computed[im] = 0;

        /* determine scattering angles (direct and surface reflected) */
        for (iw = 0; iw < aertab->nwave; iw++) {
            ig = gmult * iw;
            temp = sqrt((1.0 - geom->csenz[ig] * geom->csenz[ig]) *
                    (1.0 - geom->csolz[ig] * geom->csolz[ig])) *
                    cos(geom->phi[ig] / radeg);
            scatt1[iw] = acos(
                    MAX(-geom->csenz[ig] * geom->csolz[ig] + temp, -1.0)) * radeg;
            scatt2[iw] = acos(
                    MIN(geom->csenz[ig] * geom->csolz[ig] + temp, 1.0)) * radeg;

            /* compute Fresnel coefficients */
            fres1[iw] = fresnel_coef(geom->csenz[ig], nw);
            fres2[iw] = fresnel_coef(geom->csolz[ig], nw);
            if (gmult == 0) break;
        }

        lastsolz = geom->solz[0];
        lastsenz = geom->senz[0];
        lastphi = geom->phi[0];
    }

    if (!computed[modnum]) {

        im = modnum;
        computed[modnum] = 1;

        /* compute phase function for this geometry, all models */
        for (iw = 0; iw < aertab->nwave; iw++) {
            ig = gmult * iw;
            splint(aertab->scatt,
                    &aertab->model[im]->lnphase[iw][0],
                    &aertab->model[im]->d2phase[iw][0],
                    aertab->nscatt, scatt1[ig], &phase1);
            splint(aertab->scatt,
                    &aertab->model[im]->lnphase[iw][0],
                    &aertab->model[im]->d2phase[iw][0],
                    aertab->nscatt, scatt2[ig], &phase2);
            //          incident diffuse   reflected   diff  dir
            phase[im][iw] = exp(phase1) +
                    exp(phase2)*(fres1[ig] + fres2[ig]);
        }
    }

    return (&phase[modnum][0]);
}


/* ---------------------------------------------------------------------------------------- */
/* aeroob_cf() - out-of-band water-vapor scale factor                                       */

/* ---------------------------------------------------------------------------------------- */
float aeroob_cf(int modnum, geom_str *geom) {
    static int firstCall = 1;
    static int iw1;
    static int iw2;

    float *phase;
    float rhoas1, rhoas2;
    float eps;
    float cf;

    if (firstCall) {
        iw1 = windex(765, aertab->wave, aertab->nwave);
        iw2 = windex(865, aertab->wave, aertab->nwave);
        if (iw1 == iw2) iw1--;
        firstCall = 0;
    }

    phase = model_phase(modnum, geom);
    rhoas1 = aertab->model[modnum]->albedo[iw1] * phase[iw1] * aertab->model[modnum]->extc[iw1];
    rhoas2 = aertab->model[modnum]->albedo[iw2] * phase[iw2] * aertab->model[modnum]->extc[iw2];
    eps = rhoas1 / rhoas2;
    cf = log(eps) / (aertab->wave[iw2] - aertab->wave[iw1]);

    return (cf);
}


/* ---------------------------------------------------------------------------------------- */
/* aeroob - out-of-band water-vapor correction                                              */

/* ---------------------------------------------------------------------------------------- */
float aeroob(int32_t sensorID, int32_t iw, float airmass, float cf, float wv) {
    static float *a01;
    static float *a02;
    static float *a03;
    static float *a04;
    static float *a05;
    static float *a06;
    static float *a07;
    static float *a08;
    static float *a09;
    static float *a10;
    static float *a11;
    static float *a12;

    static int firstCall = 1;

    float f;
    f = 1.0;
    if ((input->gas_opt & GAS_TRANS_TBL_BIT) == 0) {
        if (firstCall) {
            firstCall = 0;
            rdsensorinfo(sensorID, evalmask, "oobwv01", (void **) &a01); /* coeff #1 per sensor wave */
            rdsensorinfo(sensorID, evalmask, "oobwv02", (void **) &a02);
            rdsensorinfo(sensorID, evalmask, "oobwv03", (void **) &a03);
            rdsensorinfo(sensorID, evalmask, "oobwv04", (void **) &a04);
            rdsensorinfo(sensorID, evalmask, "oobwv05", (void **) &a05);
            rdsensorinfo(sensorID, evalmask, "oobwv06", (void **) &a06);
            rdsensorinfo(sensorID, evalmask, "oobwv07", (void **) &a07);
            rdsensorinfo(sensorID, evalmask, "oobwv08", (void **) &a08);
            rdsensorinfo(sensorID, evalmask, "oobwv09", (void **) &a09);
            rdsensorinfo(sensorID, evalmask, "oobwv10", (void **) &a10);
            rdsensorinfo(sensorID, evalmask, "oobwv11", (void **) &a11);
            rdsensorinfo(sensorID, evalmask, "oobwv12", (void **) &a12);
            printf("\nLoading water-vapor correction coefficients.\n");
        }

        f = (a01[iw] + a02[iw] * airmass + cf * (a03[iw] + a04[iw] * airmass))
                + (a05[iw] + a06[iw] * airmass + cf * (a07[iw] + a08[iw] * airmass)) * wv
                + (a09[iw] + a10[iw] * airmass + cf * (a11[iw] + a12[iw] * airmass)) * wv*wv;
    }
    return (f);
}


/* ---------------------------------------------------------------------------------------- */
/* rhoa_to_rhoas() - MS aerosol reflectance to SS aerosol reflectance                       */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */

/* ---------------------------------------------------------------------------------------- */
int rhoa_to_rhoas(int32_t sensorID, int modnum, geom_str *geom, float wv,
        float rhoa[], float wave[], int32_t nwave, int iw1, int iw2, float rhoas[]) {
    float *ac, *bc, *cc;
    double a, b, c;
    double f;
    int iw, ig;
    int iwtab;
    float cf;
    int status = 0;

    ss_to_ms_coef(modnum, geom, &ac, &bc, &cc);

    cf = aeroob_cf(modnum, geom);

    for (iw = iw1; iw <= iw2; iw++) {
        ig = iw * geom->gmult;
        if (rhoa[iw] < 1.e-20)
            rhoas[iw] = rhoa[iw];
        else {
            iwtab = iwatab[iw];
            a = (double) ac[iwtab];
            b = (double) bc[iwtab];
            c = (double) cc[iwtab];
            f = b * b - 4 * c * (a - log((double) rhoa[iw]));
            if (f > 0.00001) { // this was 0.0, but small values caused segfault (BAF, 9/2014)
                if (fabs(c) > 1.e-20) {
                    rhoas[iw] = exp(0.5 * (-b + sqrt(f)) / c);
                } else if (fabs(a) > 1.e-20 && fabs(b) > 1.e-20) {
                    rhoas[iw] = pow(rhoa[iw] / a, 1. / b);
                } else {
                    status = 1;
                    break;
                }
                rhoas[iw] = rhoas[iw] / aeroob(sensorID, iw, geom->airmass[ig], cf, wv);
                if (!isfinite(rhoas[iw]) || rhoas[iw] < 1.e-20) {
                    status = 1;
                    break;
                }
            } else {
                status = 1;
                break;
            }
        }
    }

    // return input values and failure status if any wavelengths failed

    if (status != 0) {
        for (iw = iw1; iw <= iw2; iw++)
            rhoas[iw] = rhoa[iw];
    }

    return (status);
}



/* ---------------------------------------------------------------------------------------- */
/* rhoas_to_rhoa() - SS aerosol reflectance to MS aerosol reflectance                       */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */

/* ---------------------------------------------------------------------------------------- */
void rhoas_to_rhoa(int32_t sensorID, int modnum, geom_str *geom, float wv,
        float rhoas[], float wave[], int32_t nwave, int iw1, int iw2, float rhoa[]) {
    float *ac, *bc, *cc;
    float a, b, c;
    float lnrhoas;
    int iw, ig;
    int iwtab;
    float cf;

    ss_to_ms_coef(modnum, geom, &ac, &bc, &cc);

    cf = aeroob_cf(modnum, geom);

    for (iw = iw1; iw <= iw2; iw++) {
        ig = iw * geom->gmult;

        /* these changes ensure that rhoa1 -> rhoas -> rhoa2 == rhoa1 */
        /* but, that changes everything, slightly (tau)               */

        if (rhoas[iw] < 1.e-20)
            rhoa[iw] = rhoas[iw];
        else {
            iwtab = iwatab[iw];
            a = ac[iwtab];
            b = bc[iwtab];
            c = cc[iwtab];
            lnrhoas = log(rhoas[iw] * aeroob(sensorID, iw, geom->airmass[ig], cf, wv));
            rhoa[iw] = exp(a + b * lnrhoas + c * lnrhoas * lnrhoas);
        }

        /*
          iwtab = iwatab[iw];
          a = ac[iwtab];
          b = bc[iwtab];
          c = cc[iwtab];
          lnrhoas = log(rhoas[iw]);
          rhoa[iw] = exp(a + b*lnrhoas + c*lnrhoas*lnrhoas)
         * aeroob(sensorID,iw, geom->airmass[ig],cf,wv);
         */
    }

    return;
}


/* ---------------------------------------------------------------------------------------- */
/* model_epsilon() - return model epsilon at input wavelengths for input geometry.          */
/*                                                                                          */
/* If the input wavelengths are not equal to the model wavelengths, the model epsilon will  */
/* be interpolated to the input wavelengths.  It is assumed that the longest (last) input   */
/* wavelength is equivalent to the longest wavelength of the model table.  Hence,           */
/* the function should always be called with the full array of input sensor wavelengths.    */
/*                                                                                          */
/* The program optimizes for multiple calls at the same geometry by computing for all       */
/* models on the first call with a new geometry.  It returns a pointer to the internal      */
/* static arrays of epsilon for the requested model.                               .        */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */

/* ---------------------------------------------------------------------------------------- */
float *model_epsilon(int modnum, int32_t iwnir_l, float wave[], int32_t nwave, geom_str *geom) {
    static float lastsolz = -999.;
    static float lastsenz = -999.;
    static float lastphi = -999.;
    static int32_t lastiwl = -999;
    static int lastmod = -999;
    static float *epsilon[MAXMODEL];
    static int firstCall = 1;
    int i;
    float maxwave;
    /* recalculate only if geometry changes */

    if (firstCall == 1) {
        firstCall = 0;
        maxwave = MAX(aertab->nwave, nwave);
        for (i = 0; i < MAXMODEL; i++) {
            if ((epsilon[i] = (float *) calloc(maxwave, sizeof (float))) == NULL) {
                printf("Unable to allocate space for epsilon.\n");
                exit(1);
            }

        }
    }

    /* recalculate only if geometry changes */

    if (modnum != lastmod || geom->solz[0] != lastsolz ||
            geom->senz[0] != lastsenz || geom->phi[0] != lastphi ||
            iwnir_l != lastiwl) {

        int iwnir = iwatab[iwnir_l];
        float *phase;
        float *lneps;
        float rhoas1, rhoas2;
        int im, iw, iwtab;

        if ((lneps = (float *) calloc(aertab->nwave, sizeof (float))) == NULL) {
            printf("Unable to allocate space for lneps.\n");
            exit(1);
        }

        im = modnum;
        phase = model_phase(im, geom);
        for (iw = 0; iw < aertab->nwave; iw++) {
            rhoas1 = aertab->model[im]->albedo[iw] * phase[iw] * aertab->model[im]->extc[iw];
            rhoas2 = aertab->model[im]->albedo[iwnir] * phase[iwnir] * aertab->model[im]->extc[iwnir];
            epsilon[im][iw] = rhoas1 / rhoas2;
            if (interpol)
                lneps[iw] = log(epsilon[im][iw]);
        }
        if (interpol) {
            for (iw = 0; iw < nwave; iw++) {
                iwtab = iwatab[iw];
                if (aertab->wave[iwtab] != wave[iw] && wave[iw] > 0)
                    epsilon[im][iw] = exp(linterp(aertab->wave, lneps, aertab->nwave, wave[iw]));
                else
                    epsilon[im][iw] = exp(lneps[iwtab]);
            }
        }

        lastsolz = geom->solz[0];
        lastsenz = geom->senz[0];
        lastphi = geom->phi[0];
        lastiwl = iwnir_l;
        lastmod = modnum;
        free(lneps);

    }

    return (&epsilon[modnum][0]);
}


/* ---------------------------------------------------------------------------------------- */
/* model_select_wang() - M. Wang aerosol model selection process                 .          */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */

/* ---------------------------------------------------------------------------------------- */
int model_select_wang(int32_t sensorID, float wave[], int32_t nwave,
        int32_t nmodel, int32_t mindx[], geom_str *geom, float wv,
        float rhoa[], int32_t iwnir_s, int32_t iwnir_l, int32_t *modmin,
        int32_t *modmax, float *modrat, float *epsnir) {
    float *rhoas;
    float eps_ret [MAXMODEL];
    float eps_mod [MAXMODEL];
    float eps_err [MAXMODEL];
    int imod [MAXMODEL];
    int nmod = nmodel;
    float eps_ave;
    float *eps;
    float err_m;
    float err_p;
    int jm, im, iim;
    int eps_flg = 0;
    float wt;
    float tot_err;
    int itmp;

    *modmin = -1;
    *modmax = -1;
    *modrat = 0.0;
    *epsnir = 0.0;

    if ((rhoas = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoas.\n");
        exit(1);
    }

    /* get theoretical and retrieved epsilon for each model, and */
    /* compute average retrieved epsilon                         */

    eps_ave = 0.0;
    for (jm = 0; jm < nmod; jm++) {

        im = mindx[jm];

        /* retrieve epsilon in NIR assuming each model */
        rhoa_to_rhoas(sensorID, im, geom, wv, rhoa, wave, nwave, iwnir_s, iwnir_l, rhoas);

        if (rhoas[iwnir_l] > 0.0000001)
            eps_ret[jm] = rhoas[iwnir_s] / rhoas[iwnir_l];
        else
            eps_ret[jm] = 0;

        /* get model epsilon for each model at this geometry */
        eps = model_epsilon(im, iwnir_l, wave, nwave, geom);
        eps_mod[jm] = eps[iwnir_s];

        eps_ave += eps_ret[jm];
    }
    if (isfinite(eps_ave))
        eps_ave /= nmod;
    else
        eps_ave = 1.0;


    /* determine average retrieved epsilon for the four retrievals which most */
    /* closely match the theoretical model epsilons. the model set is reduced */
    /* from the full suite to the final four by iterative outlier rejection.  */

    while (nmod > 4) {

        /* compute differences between retrieved and model epsilon */
        for (im = 0; im < nmodel; im++) {
            imod[im] = im;
            eps_err[im] = eps_ave - eps_mod[im];
        }

        /* sort model indices by smallest to largest absolute differences */
        for (im = 0; im < nmodel - 1; im++) {
            for (iim = im + 1; iim < nmodel; iim++)
                if (fabs(eps_err[imod[im]]) > fabs(eps_err[imod[iim]])) {
                    itmp = imod[im ];
                    imod[im ] = imod[iim];
                    imod[iim] = itmp;
                }
        }

        /* recompute average retrieved epsilon over the n-2 closest models  */
        /* averaging is done as a weighted mean with wt=1-|eps_err|/tot_err */
        /* note that the sum of the weights is equal to n-1                 */

        nmod = nmod - 2;

        tot_err = 0.0;
        for (iim = 0; iim < nmod; iim++) {
            im = imod[iim];
            tot_err += fabs(eps_err[im]);
        }

        eps_ave = 0.0;
        for (iim = 0; iim < nmod; iim++) {
            im = imod[iim];
            wt = 1.0 - fabs(eps_err[im]) / tot_err;
            eps_ave += eps_ret[im] * wt;
        }
        eps_ave /= (nmod - 1);
    }

    /* now select the two models which bracket eps_ave  */
    err_m = 0 - FLT_MAX;
    err_p = FLT_MAX;
    for (im = 0; im < nmodel; im++) {
        eps_err[im] = eps_ave - eps_mod[im];
        if (eps_err[im] >= 0.0) {
            if (eps_err[im] < err_p) {
                err_p = eps_err[im];
                *modmin = im;
            }
        } else {
            if (eps_err[im] > err_m) {
                err_m = eps_err[im];
                *modmax = im;
            }
        }
    }

    /* compute interpolation ratio */
    if (*modmin < 0) {
        /* no lower-bounding model */
        *modmin = *modmax;
        *modrat = 0.0;
        eps_flg = -1;
    } else if (*modmax < 0) {
        /* no upper-bounding model */
        *modmax = *modmin;
        *modrat = 0.0;
        eps_flg = 1;
    } else
        *modrat = (eps_ave - eps_mod[*modmin]) / (eps_mod[*modmax] - eps_mod[*modmin]);

    *modmin = mindx[*modmin];
    *modmax = mindx[*modmax];

    /* return retrieved epsilon */
    *epsnir = eps_ave;

    free(rhoas);

    return (eps_flg);
}


/* ---------------------------------------------------------------------------------------- */
/* model_select_angst() - Select model pair based on input Angstrom coefficient             */
/*                                                                                          */
/* B. Franz, 1 August 2004.                                                                 */
/* ---------------------------------------------------------------------------------------- */
int compalphaT(alphaTstr *x, alphaTstr *y) {
    return (x->angstrom < y->angstrom ? -1 : 1);
}

void model_select_angstrom(float angstrom, int32_t *modmin, int32_t *modmax, float *modrat) {
    static alphaTstr alphaT[MAXAERMOD];
    static int firstCall = 1;

    int im, im1, im2;

    if (firstCall) {

        int ib = windex(520, aertab->wave, aertab->nwave);

        /* grab angstrom coefficients and sort in ascending order */
        for (im = 0; im < aertab->nmodel; im++) {
            alphaT[im].angstrom = aertab->model[im]->angstrom[ib];
            alphaT[im].modnum = im;
        }
        qsort(alphaT, aertab->nmodel, sizeof (alphaTstr),
                (int (*)(const void *, const void *)) compalphaT);

        firstCall = 0;
    }

    for (im = 0; im < aertab->nmodel; im++) {
        if (angstrom < alphaT[im].angstrom)
            break;
    }
    im1 = MAX(MIN(im - 1, aertab->nmodel - 1), 0);
    im2 = MAX(MIN(im, aertab->nmodel - 1), 0);

    *modmin = alphaT[im1].modnum;
    *modmax = alphaT[im2].modnum;

    if (im1 == im2)
        *modrat = 1.0;
    else
        *modrat = (angstrom - alphaT[im1].angstrom) /
        (alphaT[im2].angstrom - alphaT[im1].angstrom);

    return;

}

/* ---------------------------------------------------------------------------------------- */
/* model_taua() - compute AOT at input wavelengths using specified model                    */
/*  Note by W. Robinson - the aot is determined for SENSOR wavelength nir_l
    using the mu, mu0 for the geometry for that band.  To get the aot at other
    TABLE bands, it starts with the aot(nir_l).  So, how would the geometry
    change get properly translated to the other table bands?  For now, I'm
    not accounting for that translation.  Possible method would be to (for
    interpol=0 only) compute the aot(nir_l, geom(sensor band)) and use
    that aot to get the aot(geom(sensor band))                                              */

/* ---------------------------------------------------------------------------------------- */
void model_taua(int32_t sensorID, int modnum, float wave[], int32_t nwave,
        int32_t iwnir_l, float rhoa[], geom_str *geom, float wv, float taua[]) {
    float *aot;
    float *lnaot;
    float *rhoas;

    int iwnir = iwatab[iwnir_l];
    float *phase = model_phase(modnum, geom);
    int iw, iwg, iwtab;
    float maxwave;

    iwg = iwnir * geom->gmult; /* index for geometry */

    maxwave = MAX(aertab->nwave, nwave);

    if ((rhoas = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoas.\n");
        exit(1);
    }
    if ((aot = (float *) calloc(maxwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for aot.\n");
        exit(1);
    }
    if ((lnaot = (float *) calloc(maxwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for lnaot.\n");
        exit(1);
    }

    /* get SS aerosol reflectance at longest sensor wavelength */
    rhoa_to_rhoas(sensorID, modnum, geom, wv, rhoa, wave, nwave, iwnir_l, iwnir_l, rhoas);

    /* get aerosol optical thickness at longest sensor wavelength */
    aot[iwnir] = rhoas[iwnir_l]*(4.0 * geom->csolz[iwg] * geom->csenz[iwg])
            / (phase[iwnir] * aertab->model[modnum]->albedo[iwnir]);

    /* get aerosol optical thickness at all other table wavelengths */
    for (iw = 0; iw < aertab->nwave; iw++) {
        /* note to self: i actually have aot(865) for czcs at this point */
        aot[iw] = aot[iwnir] * aertab->model[modnum]->extc[iw] / aertab->model[modnum]->extc[iwnir];
        if (interpol)
            lnaot[iw] = log(aot[iw]);
    }

    /* interpolate aot from table to sensor wavelengths */
    if (interpol) {
        for (iw = 0; iw < nwave; iw++) {
            iwtab = iwatab[iw];
            if (aertab->wave[iwtab] != wave[iw] && wave[iw] > 0)
                taua[iw] = exp(linterp(aertab->wave, lnaot, aertab->nwave, wave[iw]));
            else
                taua[iw] = aot[iwtab];
        }
    } else
        for (iw = 0; iw < nwave; iw++)
            taua[iw] = aot[iw];

    free(rhoas);
    free(aot);
    free(lnaot);
    return;
}

/*------------------------------------------------------------------------------------------ */
/* model_taua_mseps_pca() -  for a given model and rhoa[iwnir_l], computes the AOT using PCA  LUTs */

/*------------------------------------------------------------------------------------------ */
int model_taua_mseps_pca(int32_t modl, float wave[], int32_t nwave, int32_t iwnir_l,
             float rhoa[], geom_str *geom, float tau_pred[])
{
    static int firstcall=1;
    static float **pc_rhoa;
    static int npc, ntau_870;
    aermodstr *aermod=aertab->model[modl];

    int iw,itau,ipc;
    static float *rhoa_l;

    if(firstcall){
        firstcall=0;
        ntau_870=aertab->ntau_870;
        npc=aertab->npc;
        pc_rhoa=(float **)malloc(aertab->ntau_870*sizeof(float *));

        for(itau=0;itau<ntau_870;itau++)
            pc_rhoa[itau]=(float *)malloc(aertab->npc*sizeof(float));
        rhoa_l=(float *)malloc(ntau_870*sizeof(float));
    }

    get_pc_rhoa(modl,geom,pc_rhoa);

    for(itau=0;itau<ntau_870;itau++){
        rhoa_l[itau]=0;
        for(ipc=0;ipc<npc;ipc++){
            rhoa_l[itau]+=pc_rhoa[itau][ipc]*aermod->pc_components_rhoa[ipc][iwnir_l];
            rhoa_l[itau]+=aermod->pc_mean_rhoa[iwnir_l];
        }
    }

    for(itau=0;itau<ntau_870;itau++){
        if(rhoa[iwnir_l]<rhoa_l[itau] && rhoa[iwnir_l]>rhoa_l[itau-1])
            break;
    }
    tau_pred[iwnir_l]=aermod->tau_870[itau]+(aermod->tau_870[itau]-aermod->tau_870[itau-1])*(rhoa[iwnir_l]-rhoa_l[itau])/(rhoa_l[itau]-rhoa_l[itau-1]);

    for(iw=0;iw<nwave;iw++)
        tau_pred[iw]=tau_pred[iwnir_l]*aermod->extc[iw]/aermod->extc[iwnir_l];

    return 0;

}

/*------------------------------------------------------------------------------------------ */
/* model_taua_mseps() -  for a given model and rhoa[iwnir_l], computes the AOT using MSEPS   */
/*                       coefficients relating rhoa to taua (in log space)                   */
/*------------------------------------------------------------------------------------------ */
int model_taua_mseps(int32_t modl, float wave[], int32_t nwave, int32_t iwnir_l, 
		     float rhoa[], geom_str *geom, float tau_pred[])
 {
    float *ac, *bc, *cc, *dc, *ec;
    double ax, bx, cx, fx;
    float ext_modl[nwave];
    float tau_iwnir_l;
    int iw, iwtab;
    int iwtab_l = iwatab[iwnir_l];

    /* get the coefficients for lg_rho vs lg_aot  */
    ms_eps_coef(modl, iwnir_l, wave, geom, &ac, &bc, &cc, &dc, &ec);


    ax = (double) ac[iwtab_l] - log((double) rhoa[iwnir_l]);
    bx = (double) bc[iwtab_l];
    cx = (double) cc[iwtab_l];

    fx = bx * bx - 4.0 * ax*cx;
    if (fx > 0.0 && cx != 0.0) {
        tau_iwnir_l = exp(0.5*(-bx + sqrt(fx)) / cx);
    } else {
        tau_iwnir_l = 0.0;
    }

    /* get the extinction coefficients and compute AOT at all wavelengths */
    for (iw = 0; iw < nwave; iw++) {
        iwtab = iwatab[iw];
        ext_modl[iw] = aertab->model[modl]->extc[iwtab];
    }
    for (iw = 0; iw < nwave; iw++) {
        tau_pred[iw] = (ext_modl[iw] / ext_modl[iwnir_l]) * tau_iwnir_l;
    }

    return (0);
}

/* ---------------------------------------------------------------------------------------- */
/* model_transmittance() - compute path Rayleigh-aerosol diffuse trans for specified model  */
/*                                                                                          */
/*  W. Robinson, SAIC, 24 Mar 2017, modify for band-dependent geometry                      */
/*  M. Zhang,    SAIC,  July 2021,  adding the uncertainty propagation                      */

/* ---------------------------------------------------------------------------------------- */
void model_transmittance(int modnum, float wave[], int32_t nwave,
        float *theta, int gmult, float taua[], float dtran[], float dt[]) {

    int i1, i2, i, ig;
    int iw, iwtab;
    float a1, b1;
    float a2, b2;
    float x1, x2;
    float y1, y2;
    float xbar;
    float wt;

    /* disable interpolation, since we now allow LUT wavelengths to differ from sensor nominal wavelengths
       BAF, 6/2022

    if (firstCall) {
        float taur1, um1;
        float taur2, um2;
        firstCall = 0;

        intexp = (float *) malloc(nwave * sizeof (float));
        inttst = (int *) malloc(nwave * sizeof (int));

        for (iw = 0; iw < nwave; iw++) {
            intexp[iw] = 1.0;
            inttst[iw] = 0;
            iwtab = iwdtab[iw];
            if (fabs(wave[iw] - aertab->dtran_wave[iwtab]) > 0.51) {
                um1 = aertab->dtran_wave[iwtab] / 1000.0;
                um2 = wave[iw] / 1000.0;
                taur1 = 0.008569 * pow(um1, -4)*(1.0 + (0.0113 * pow(um1, -2))+(0.00013 * pow(um1, -4)));
                taur2 = 0.008569 * pow(um2, -4)*(1.0 + (0.0113 * pow(um2, -2))+(0.00013 * pow(um2, -4)));
                intexp[iw] = taur2 / taur1;
                inttst[iw] = 1;
                printf("Interpolating diffuse transmittance for %d from %f by %f\n",
                        (int) wave[iw], aertab->dtran_wave[iwtab], intexp[iw]);
            }
        }
    }
    */

    /* use coefficients of nearest wavelength */
    for (iw = 0; iw < nwave; iw++) {
        if ((iw == 0) || (gmult != 0)) {
            ig = iw * gmult;
            /* find bracketing zenith angle indices */
            for (i = 0; i < aertab->dtran_ntheta; i++) {
                if (theta[ig] < aertab->dtran_theta[i])
                    break;
            }
            if (i == aertab->dtran_ntheta) {
                i1 = i - 1;
                i2 = i - 1;
                wt = 0.0;
            } else {
                i1 = MIN(MAX(i - 1, 0), aertab->dtran_ntheta - 2);
                i2 = i1 + 1;
                x1 = aertab->dtran_airmass[i1];
                x2 = aertab->dtran_airmass[i2];
                xbar = 1.0 / cos(theta[ig] / radeg);
                wt = (xbar - x1) / (x2 - x1);
            }
        }
        iwtab = iwdtab[iw];

        a1 = aertab->model[modnum]->dtran_a[iwtab][i1];
        b1 = aertab->model[modnum]->dtran_b[iwtab][i1];

        a2 = aertab->model[modnum]->dtran_a[iwtab][i2];
        b2 = aertab->model[modnum]->dtran_b[iwtab][i2];

	/* disable interpolation, BAF, 6/2022
        if (inttst[iw]) {
            a1 = pow(a1, intexp[iw]);
            a2 = pow(a2, intexp[iw]);
        }
	*/

        y1 = a1 * exp(-b1 * taua[iw]);
        y2 = a2 * exp(-b2 * taua[iw]);

        dtran[iw] = MAX(MIN((1.0 - wt) * y1 + wt*y2, 1.0), 1e-5);

        if(dt){
            if( fabs(dtran[iw]-1.0)<0.0000001 || fabs(dtran[iw]-1e-5)<0.0000001 )
                dt[iw]=0;
            else
                dt[iw]=-(1.0-wt)*y1*b1-wt*y2*b2;
        }
    }

    return;
}
void model_transmittance_pca(int modnum, float wave[], int32_t nwave,
        float *theta, int gmult, float taua[], float dtran[], float dt[]) {

    static int firstcall=1;
    static float *pc_td, *solz;
    static float *deriv_pc; //derivative of pc_td to taua[npc]
    static int npc, ntau_870,nsolz,nir_base;
    aermodstr *aermod=aertab->model[modnum];
    float a00,a01,a10,a11;

    int iw,itau,itheta,ipc,itau1,itau2,itheta1,itheta2;
    float r,p,deriv_r=0.;

    if(firstcall){
        firstcall=0;
        nir_base=windex(input->aer_wave_base,wave,nwave);
        solz=aertab->senz;
        nsolz=aertab->nsenz;
        ntau_870=aertab->ntau_870;
        npc=aertab->npc;
        pc_td=(float *)malloc(aertab->npc*sizeof(float));
        if(dt)
            deriv_pc=(float *)malloc(npc*sizeof(float));
    }

    /// interpolate the pc_td using tau_870 and *theta

    if(taua[nir_base]<aermod->tau_870[0]){
        itau1=0;
        itau2=1;
    }else if(taua[nir_base]>=aermod->tau_870[ntau_870-1]){
        itau1=ntau_870-2;
        itau2=ntau_870-1;
    }else{
        for(itau=0;itau<ntau_870;itau++)
            if(taua[nir_base]<aermod->tau_870[itau])
                break;
        itau1=itau-1;
        itau2=itau;
    }
    if(itau1!=itau2){
        if(dt)
            deriv_r=1/(aermod->tau_870[itau2]-aermod->tau_870[itau1]);
        r=(taua[nir_base]-aermod->tau_870[itau1])/(aermod->tau_870[itau2]-aermod->tau_870[itau1]);
    }
    else
        r=0.0;

    for(itheta=0;itheta<nsolz;itheta++)
        if(*theta<solz[itheta])
            break;
    itheta1=MAX(itheta-1,0);
    itheta2=MIN(itheta,nsolz-1);
    if(itheta1!=itheta2)
        p=(*theta-solz[itheta1])/(solz[itheta2]-solz[itheta1]);
    else
        p=0.0;

    for(ipc=0;ipc<npc;ipc++){
        a00=aermod->pc_td[itau1][itheta1][ipc];
        a01=aermod->pc_td[itau1][itheta2][ipc];
        a10=aermod->pc_td[itau2][itheta1][ipc];
        a11=aermod->pc_td[itau2][itheta2][ipc];

        pc_td[ipc]=a00*(1-r)*(1-p)+a01*(1-r)*p+a10*r*(1-p)+a11*r*p;
        if(dt)
            deriv_pc[ipc]=(-a00*(1-p)-a01*p+a10*(1-p)+a11*p)*deriv_r;
    }

    for(iw=0;iw<nwave;iw++){
        dtran[iw]=0.;
        if(dt)
            dt[iw]=0.;
        for(ipc=0;ipc<npc;ipc++){
            dtran[iw]+=pc_td[ipc]*aermod->pc_components_td[ipc][iw];
            if(dt)
                dt[iw]+=aermod->pc_components_td[ipc][iw]*deriv_pc[ipc];
        }
        dtran[iw]+=aermod->pc_mean_td[iw];
        dtran[iw]=exp(dtran[iw]);
        if(dt)
            dt[iw]*=dtran[iw];
    }
    return;
}
void diff_tran_pca(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_l,
        geom_str *geom, float wv, float pr, float taur[], int32_t modmin,
        int32_t modmax, float modrat, float rhoa[], float taua[], float tsol[],
        float tsen[], float tauamin[], float tauamax[], int taua_opt, int ip,uncertainty_t *uncertainty) {

    int iw, gmult, ig,inir,ib;
    float *tsolmin;
    float *tsolmax;
    float *tsenmin;
    float *tsenmax;
    float *dtmin=NULL,*dtmax=NULL;// derivative of tsol[iw] or tsen [iw] with respect to taua[nir_l]
    float tmp;
    static int firstcall=1;
    static int32_t nbands_ac, *acbands_index;


    if(firstcall){
        firstcall=0;
        nbands_ac=input->nbands_ac;
        acbands_index=input->acbands_index;
    }


    if ((tsolmin = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsolmin.\n");
        exit(1);
    }
    if ((tsolmax = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsolmax.\n");
        exit(1);
    }
    if ((tsenmin = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsenmin.\n");
        exit(1);
    }
    if ((tsenmax = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsenmax.\n");
        exit(1);
    }

    if(uncertainty){
        iwnir_s=bindex_get(input->aer_wave_short);
        //dmodrat=uncertainty->dmodrat;
        if ((dtmin = (float *) calloc(nwave, sizeof(float))) == NULL) {
            printf("Unable to allocate space for dtmin.\n");
            exit(1);
        }
        if ((dtmax = (float *) calloc(nwave, sizeof(float))) == NULL) {
            printf("Unable to allocate space for dtmax.\n");
            exit(1);
        }
    }


    gmult = geom->gmult;
    if (interpol == 1) gmult = 0; /* geom used for model wavelengths */

    /* get AOT per band for each model, if not already computed */
    if (taua_opt == 0) {
        // using single scattering approximation space (Gordan and Wang)
        model_taua(sensorID, modmin, wave, nwave, iwnir_l, rhoa, geom, wv, tauamin);
        model_taua(sensorID, modmax, wave, nwave, iwnir_l, rhoa, geom, wv, tauamax);
    //printf("%d %d %d %f %f %f\n",taua_opt,modmin, iwnir_l, rhoa[iwnir_l], tauamin[0], tauamin[iwnir_l]);
    } else if (taua_opt == 2) {
        // using multi-scattering relationship between rhoa and taua (Ahmad)
        model_taua_mseps_pca(modmin, wave, nwave, iwnir_l, rhoa, geom, tauamin);
        model_taua_mseps_pca(modmax, wave, nwave, iwnir_l, rhoa, geom, tauamax);
    //printf("%d %d %d %f %f %f\n",taua_opt,modmin, iwnir_l, rhoa[iwnir_l], tauamin[0], tauamin[iwnir_l]);
    }

    /* get diff trans sun to ground, per band for each model */
    model_transmittance_pca(modmin, wave, nwave, geom->solz, gmult, tauamin, tsolmin, dtmin);
    model_transmittance_pca(modmax, wave, nwave, geom->solz, gmult, tauamax, tsolmax, dtmax);

    for (iw = 0; iw < nwave; iw++) {
        ig = iw * geom->gmult; /* geom used for sensor wavelengths */
        tsol[iw] = tsolmin[iw]*(1.0 - modrat) + tsolmax[iw] * modrat;
        //uncertainty->dt_sol[ip*nwave+iw]=sqrt( pow((1.0-modrat)*dtmin[iw],2) + pow(modrat*dtmax[iw],2) + pow((tsolmax[iw]-tsolmin[iw])*dmodrat[ip],2) );

        if(uncertainty){
            for(ib=0;ib<nbands_ac;ib++){
                inir=acbands_index[ib]-iwnir_s;
                uncertainty->derv_tsol_rhorc[iw][inir]=(tsolmax[iw]-tsolmin[iw])*uncertainty->derv_modrat_rhorc[inir];
                if(acbands_index[ib]==iwnir_l){
                    uncertainty->derv_tsol_rhorc[iw][inir]+=(1.0 - modrat)*dtmin[iw]*uncertainty->derv_taua_min_rhorc_l[iwnir_l];
                    uncertainty->derv_tsol_rhorc[iw][inir]+=  modrat      *dtmax[iw]*uncertainty->derv_taua_max_rhorc_l[iwnir_l];
                }
            }
            uncertainty->derv_tsol_taua_l[iw]=(tsolmax[iw]-tsolmin[iw])*uncertainty->derv_modrat_taua_l;
            uncertainty->derv_tsol_taua_l[iw]+=(1.0 - modrat)*dtmin[iw]*uncertainty->derv_taua_min_taua_l[iw];
            uncertainty->derv_tsol_taua_l[iw]+=  modrat      *dtmax[iw]*uncertainty->derv_taua_max_taua_l[iw];

            uncertainty->derv_tsol_rhow_l[iw]=(tsolmax[iw]-tsolmin[iw])*uncertainty->derv_modrat_rhow_l;
            uncertainty->derv_tsol_rhow_l[iw]+=(1.0 - modrat)*dtmin[iw]*uncertainty->derv_taua_min_rhow_l[iw];
            uncertainty->derv_tsol_rhow_l[iw]+=modrat        *dtmax[iw]*uncertainty->derv_taua_max_rhow_l[iw];
        }

        /* correct for pressure difference from standard pressure */
        double tmp_pressure_diff = exp(-0.5 * taur[iw] / geom->csolz[ig]*(pr / p0 - 1));
        tsol[iw] = tsol[iw] * tmp_pressure_diff;

        //uncertainty->dt_sol[ip*nwave+iw]*=tmp_pressure_diff;

        if(uncertainty){
            for(ib=0;ib<nbands_ac;ib++){
                inir=acbands_index[ib]-iwnir_s;
                uncertainty->derv_tsol_rhorc[iw][inir]*=tmp_pressure_diff;
            }
            uncertainty->derv_tsol_taua_l[iw]*=tmp_pressure_diff;
            uncertainty->derv_tsol_rhow_l[iw]*=tmp_pressure_diff;
        }

        if ((evalmask & TRANSSPHER) != 0) {
            /* correct for airmass difference, plane-parallel to spherical atmosphere */
            tmp=tsol[iw];
            tsol[iw] = pow(tsol[iw], geom->airmass_sph[ig] / geom->airmass_plp[ig]);

            if(uncertainty)
                uncertainty->dt_sol[ip*nwave+iw]*= ( geom->airmass_sph[ig] / geom->airmass_plp[ig] *pow(tmp,geom->airmass_sph[ig] / geom->airmass_plp[ig]-1) );
        }
    }

    /* get diff trans ground to sensor, per band for each model */
    model_transmittance_pca(modmin, wave, nwave, geom->senz, gmult, tauamin, tsenmin, dtmin);
    model_transmittance_pca(modmax, wave, nwave, geom->senz, gmult, tauamax, tsenmax, dtmax);

    /* interpolate and pressure correct */
    for (iw = 0; iw < nwave; iw++) {
        ig = iw * geom->gmult; /* geom used for sensor wavelengths */
        taua[iw] = tauamin[iw]*(1.0 - modrat) + tauamax[iw] * modrat;
        tsen[iw] = tsenmin[iw]*(1.0 - modrat) + tsenmax[iw] * modrat;

        if(uncertainty){
            for(ib=0;ib<nbands_ac;ib++){
                inir=acbands_index[ib]-iwnir_s;
                uncertainty->derv_tsen_rhorc[iw][inir]=(tsenmax[iw]-tsenmin[iw])*uncertainty->derv_modrat_rhorc[inir];
                if(acbands_index[ib]==iwnir_l){
                    uncertainty->derv_tsen_rhorc[iw][inir]+=(1.0 - modrat)*dtmin[iw]*uncertainty->derv_taua_min_rhorc_l[iwnir_l];
                    uncertainty->derv_tsen_rhorc[iw][inir]+=  modrat      *dtmax[iw]*uncertainty->derv_taua_max_rhorc_l[iwnir_l];
                }
            }
            uncertainty->derv_tsen_taua_l[iw]=(tsenmax[iw]-tsenmin[iw])*uncertainty->derv_modrat_taua_l;
            uncertainty->derv_tsen_taua_l[iw]+=(1.0 - modrat)*dtmin[iw]*uncertainty->derv_taua_min_taua_l[iw];
            uncertainty->derv_tsen_taua_l[iw]+=  modrat      *dtmax[iw]*uncertainty->derv_taua_max_taua_l[iw];

            uncertainty->derv_tsen_rhow_l[iw]=(tsenmax[iw]-tsenmin[iw])*uncertainty->derv_modrat_rhow_l;
            uncertainty->derv_tsen_rhow_l[iw]+=(1.0 - modrat)*dtmin[iw]*uncertainty->derv_taua_min_rhow_l[iw];
            uncertainty->derv_tsen_rhow_l[iw]+=modrat        *dtmax[iw]*uncertainty->derv_taua_max_rhow_l[iw];
        }
        /* correct for pressure difference from standard pressure */
        double tmp_pressure_diff = exp(-0.5 * taur[iw] / geom->csenz[ig] *(pr / p0 - 1));
        tsen[iw] = tsen[iw] * tmp_pressure_diff;

        if(uncertainty){
            for(ib=0;ib<nbands_ac;ib++){
                inir=acbands_index[ib]-iwnir_s;
                uncertainty->derv_tsen_rhorc[iw][inir]*=tmp_pressure_diff;
            }
            uncertainty->derv_tsen_taua_l[iw]*=tmp_pressure_diff;
            uncertainty->derv_tsen_rhow_l[iw]*=tmp_pressure_diff;
        }

        if ((evalmask & TRANSSPHER) != 0) {
            /* correct for airmass difference, plane-parallel to spherical atmosphere */
            tmp=tsen[iw];
            tsen[iw] = pow(tsen[iw], geom->airmass_sph[ig] / geom->airmass_plp[ig]);
            if(uncertainty)
                uncertainty->dt_sen[ip*nwave+iw]*= ( geom->airmass_sph[ig] / geom->airmass_plp[ig] *pow(tmp,geom->airmass_sph[ig] / geom->airmass_plp[ig]-1) );
        }
    }

    free(tsolmin);
    free(tsolmax);
    free(tsenmin);
    free(tsenmax);
    if(uncertainty){
        free(dtmin);
        free(dtmax);
    }

    return;
}

/* ---------------------------------------------------------------------------------------- */
/* diff_tran() - compute Rayleigh-aerosol diffuse trans for selected model pair, both paths */
/* ---------------------------------------------------------------------------------------- */
void diff_tran(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_l,
        geom_str *geom, float wv, float pr, float taur[], int32_t modmin,
        int32_t modmax, float modrat, float rhoa[], float taua[], float tsol[],
        float tsen[], float tauamin[], float tauamax[], int taua_opt, int ip,uncertainty_t *uncertainty) {

    int iw, gmult, ig;
    float *tsolmin;
    float *tsolmax;
    float *tsenmin;
    float *tsenmax;
    float *dtmin=NULL,*dtmax=NULL;// derivative of tsol or tsen with respect to the corresponding taua
    float tmp;

    int iwnir_s;
    int inir;

    if(use_pca_lut){
        diff_tran_pca(sensorID,wave,nwave,iwnir_l,geom,wv,pr,taur,modmin,modmax,modrat,rhoa,taua,tsol,tsen,\
                tauamin,tauamax,taua_opt,ip,uncertainty);
        return;
    }

    if ((tsolmin = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsolmin.\n");
        exit(1);
    }
    if ((tsolmax = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsolmax.\n");
        exit(1);
    }
    if ((tsenmin = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsenmin.\n");
        exit(1);
    }
    if ((tsenmax = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsenmax.\n");
        exit(1);
    }

    if(uncertainty){
        iwnir_s=bindex_get(input->aer_wave_short);
        //dmodrat=uncertainty->dmodrat;
        if ((dtmin = (float *) calloc(nwave, sizeof(float))) == NULL) {
            printf("Unable to allocate space for dtmin.\n");
            exit(1);
        }
        if ((dtmax = (float *) calloc(nwave, sizeof(float))) == NULL) {
            printf("Unable to allocate space for dtmax.\n");
            exit(1);
        }
    }


    gmult = geom->gmult;
    if (interpol == 1) gmult = 0; /* geom used for model wavelengths */

    /* get AOT per band for each model, if not already computed */
    if (taua_opt == 0) {
        // using single scattering approximation space (Gordan and Wang) 
        model_taua(sensorID, modmin, wave, nwave, iwnir_l, rhoa, geom, wv, tauamin);
        model_taua(sensorID, modmax, wave, nwave, iwnir_l, rhoa, geom, wv, tauamax);
	//printf("%d %d %d %f %f %f\n",taua_opt,modmin, iwnir_l, rhoa[iwnir_l], tauamin[0], tauamin[iwnir_l]);
    } else if (taua_opt == 2) {
        // using multi-scattering relationship between rhoa and taua (Ahmad)
        model_taua_mseps(modmin, wave, nwave, iwnir_l, rhoa, geom, tauamin);
        model_taua_mseps(modmax, wave, nwave, iwnir_l, rhoa, geom, tauamax);
	//printf("%d %d %d %f %f %f\n",taua_opt,modmin, iwnir_l, rhoa[iwnir_l], tauamin[0], tauamin[iwnir_l]);
    }

    /* get diff trans sun to ground, per band for each model */
    model_transmittance(modmin, wave, nwave, geom->solz, gmult, tauamin, tsolmin, dtmin);
    model_transmittance(modmax, wave, nwave, geom->solz, gmult, tauamax, tsolmax, dtmax);

    for (iw = 0; iw < nwave; iw++) {
    	ig = iw * geom->gmult; /* geom used for sensor wavelengths */
    	tsol[iw] = tsolmin[iw]*(1.0 - modrat) + tsolmax[iw] * modrat;
    	//uncertainty->dt_sol[ip*nwave+iw]=sqrt( pow((1.0-modrat)*dtmin[iw],2) + pow(modrat*dtmax[iw],2) + pow((tsolmax[iw]-tsolmin[iw])*dmodrat[ip],2) );

    	if(uncertainty){
    	    for(inir=iwnir_s;inir<=iwnir_l;inir++){
    	        uncertainty->derv_tsol_rhorc[iw][inir-iwnir_s]=(tsolmax[iw]-tsolmin[iw])*uncertainty->derv_modrat_rhorc[inir-iwnir_s];
    	        if(inir==iwnir_l){
    	            uncertainty->derv_tsol_rhorc[iw][inir-iwnir_s]+=(1.0 - modrat)*dtmin[iw]*uncertainty->derv_taua_min_rhorc_l[iw];
    	            uncertainty->derv_tsol_rhorc[iw][inir-iwnir_s]+=  modrat      *dtmax[iw]*uncertainty->derv_taua_max_rhorc_l[iw];
    	        }
    	    }
            uncertainty->derv_tsol_taua_l[iw]=(tsolmax[iw]-tsolmin[iw])*uncertainty->derv_modrat_taua_l;
            uncertainty->derv_tsol_taua_l[iw]+=(1.0 - modrat)*dtmin[iw]*uncertainty->derv_taua_min_taua_l[iw];
            uncertainty->derv_tsol_taua_l[iw]+=  modrat      *dtmax[iw]*uncertainty->derv_taua_max_taua_l[iw];

            uncertainty->derv_tsol_rhow_l[iw]=(tsolmax[iw]-tsolmin[iw])*uncertainty->derv_modrat_rhow_l;
            uncertainty->derv_tsol_rhow_l[iw]+=(1.0 - modrat)*dtmin[iw]*uncertainty->derv_taua_min_rhow_l[iw];
            uncertainty->derv_tsol_rhow_l[iw]+=modrat        *dtmax[iw]*uncertainty->derv_taua_max_rhow_l[iw];
    	}

    	/* correct for pressure difference from standard pressure */
        double tmp_pressure_diff = exp(-0.5 * taur[iw] / geom->csolz[ig]*(pr / p0 - 1));
    	tsol[iw] = tsol[iw] * tmp_pressure_diff;

    	//uncertainty->dt_sol[ip*nwave+iw]*=tmp_pressure_diff;

    	if(uncertainty){
    	    for(inir=iwnir_s;inir<=iwnir_l;inir++)
    	        uncertainty->derv_tsol_rhorc[iw][inir-iwnir_s]*=tmp_pressure_diff;
    	    uncertainty->derv_tsol_taua_l[iw]*=tmp_pressure_diff;
    	    uncertainty->derv_tsol_rhow_l[iw]*=tmp_pressure_diff;
    	}

    	if ((evalmask & TRANSSPHER) != 0) {
    		/* correct for airmass difference, plane-parallel to spherical atmosphere */
    		tmp=tsol[iw];
    		tsol[iw] = pow(tsol[iw], geom->airmass_sph[ig] / geom->airmass_plp[ig]);

    		if(uncertainty)
    		    uncertainty->dt_sol[ip*nwave+iw]*= ( geom->airmass_sph[ig] / geom->airmass_plp[ig] *pow(tmp,geom->airmass_sph[ig] / geom->airmass_plp[ig]-1) );
    	}
    }

    /* get diff trans ground to sensor, per band for each model */
    model_transmittance(modmin, wave, nwave, geom->senz, gmult, tauamin, tsenmin, dtmin);
    model_transmittance(modmax, wave, nwave, geom->senz, gmult, tauamax, tsenmax, dtmax);

    /* interpolate and pressure correct */
    for (iw = 0; iw < nwave; iw++) {
        ig = iw * geom->gmult; /* geom used for sensor wavelengths */
        taua[iw] = tauamin[iw]*(1.0 - modrat) + tauamax[iw] * modrat;
        tsen[iw] = tsenmin[iw]*(1.0 - modrat) + tsenmax[iw] * modrat;
   
        if(uncertainty){
            for(inir=iwnir_s;inir<=iwnir_l;inir++){
                uncertainty->derv_tsen_rhorc[iw][inir-iwnir_s]=(tsenmax[iw]-tsenmin[iw])*uncertainty->derv_modrat_rhorc[inir-iwnir_s];
                if(inir==iwnir_l){
                    uncertainty->derv_tsen_rhorc[iw][inir-iwnir_s]+=(1.0 - modrat)*dtmin[iw]*uncertainty->derv_taua_min_rhorc_l[iw];
                    uncertainty->derv_tsen_rhorc[iw][inir-iwnir_s]+=  modrat      *dtmax[iw]*uncertainty->derv_taua_max_rhorc_l[iw];
                }
            }
            uncertainty->derv_tsen_taua_l[iw]=(tsenmax[iw]-tsenmin[iw])*uncertainty->derv_modrat_taua_l;
            uncertainty->derv_tsen_taua_l[iw]+=(1.0 - modrat)*dtmin[iw]*uncertainty->derv_taua_min_taua_l[iw];
            uncertainty->derv_tsen_taua_l[iw]+=  modrat      *dtmax[iw]*uncertainty->derv_taua_max_taua_l[iw];

            uncertainty->derv_tsen_rhow_l[iw]=(tsenmax[iw]-tsenmin[iw])*uncertainty->derv_modrat_rhow_l;
            uncertainty->derv_tsen_rhow_l[iw]+=(1.0 - modrat)*dtmin[iw]*uncertainty->derv_taua_min_rhow_l[iw];
            uncertainty->derv_tsen_rhow_l[iw]+=modrat        *dtmax[iw]*uncertainty->derv_taua_max_rhow_l[iw];
        }
        /* correct for pressure difference from standard pressure */
        double tmp_pressure_diff = exp(-0.5 * taur[iw] / geom->csenz[ig] *(pr / p0 - 1));
        tsen[iw] = tsen[iw] * tmp_pressure_diff;
     
        if(uncertainty){
            for(inir=iwnir_s;inir<=iwnir_l;inir++)
                uncertainty->derv_tsen_rhorc[iw][inir-iwnir_s]*=tmp_pressure_diff;
            uncertainty->derv_tsen_taua_l[iw]*=tmp_pressure_diff;
            uncertainty->derv_tsen_rhow_l[iw]*=tmp_pressure_diff;
        }

        if ((evalmask & TRANSSPHER) != 0) {
            /* correct for airmass difference, plane-parallel to spherical atmosphere */
            tmp=tsen[iw];
            tsen[iw] = pow(tsen[iw], geom->airmass_sph[ig] / geom->airmass_plp[ig]);
            if(uncertainty)
                uncertainty->dt_sen[ip*nwave+iw]*= ( geom->airmass_sph[ig] / geom->airmass_plp[ig] *pow(tmp,geom->airmass_sph[ig] / geom->airmass_plp[ig]-1) );
        }
    }

    free(tsolmin);
    free(tsolmax);
    free(tsenmin);
    free(tsenmax);
    if(uncertainty){
        free(dtmin);
        free(dtmax);
    }

    return;
}

/* ---------------------------------------------------------------------------------------- */
/* smaer_pca() - spectral matching approach based on PCA LUTs                            */
/*
 * M. Zhang, Dec. 2022.                                                                     */

/* ---------------------------------------------------------------------------------------- */

int smaer_pca(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s,
        int32_t iwnir_l, int32_t nmodels, int32_t mindx[], geom_str *geom, float wv, float rhoa[], int32_t *modmin, int32_t *modmax,
        float32 *modrat, float *tau_pred_min, float *tau_pred_max, int32_t ip, float chi[], float *mbac_w,uncertainty_t *uncertainty) {
    static float *tau_iwnir_l;//[nmodels]
    static int firstcall=1, nir_l;
    static float **pc_rhoa;
    static int npc, ntau_870;
    aermodstr *aermod;
    static int wave_base_index; //index of the nir_l in the mbac_wave

    int iw, ib, itau,ipc,itemp;

    int status = 0;

    static float *rho_all_wav_pred,*rho_pred_min,*rho_pred_max;

    int im, modl;

    float diff_2 = 0.0,diff_1=0.0;
    float mbac_wsum = 0.0;
    static float *chi_old;

    float *noise=&noise_global[ip*nwave];

    static float **derv_chi_rhorc;//dimension: [nmodels][nbands_ac], deivative of chi to rhorc
    static float *derv_chi_taua_l; //dimension: [nmodels], deivative of chi to taua_l
    static float **derv_rhoa_rhorc_l;// derivative of modeled rhoa[nwave] to rhorc[iwnir_l]
    float *derv_modrat_rhorc;// dimension: [nbands_ac], deivative of modrat to rhorc[iwnir_s to iwnir_l]
    float derv_modrat_taua_l=0.; // deivative of modrat to taua [iwnir_l]
    float derv_modrat_rhow_l=0.;
    float *derv_taua_min_rhorc_l;
    float *derv_taua_min_taua_l;
    float *derv_taua_min_rhow_l;
    float *derv_taua_max_rhorc_l;
    float *derv_taua_max_taua_l;
    float *derv_taua_max_rhow_l;
    static float *derv_taua_rhorc_l;//[nmodels], derivative of modeled taua[nir_l] to rhorc_l
    static float *derv_chi_obs_rhorc; //[nbands_ac], derivative of chi_obs to rhorc at nbands_ac
    float derv_chi_obs_taua_l; // derivative of chi_obs to taua_l

    static float *derv_rhoa_min=NULL,*derv_rhoa_max=NULL;  //derivative of rhoa[]to taua[aer_l] for aermodmin and aermodmax
    static float *derv_taua_min=NULL,*derv_taua_max=NULL;  //derivative of taua[]to taua[aer_l] for aermodmin and aermodmax
    

    float derv_temp_rhorc;
    float derv_modrat_chi0,derv_modrat_chi1;

    static float **rhoa_modl; //[ntau_870][nwave]
    struct str chi_struct[nmodels];

    static int32_t nbands_ac, *acbands_index;

    if(firstcall){
        firstcall=0;
        nir_l=bindex_get(input->aer_wave_base);

        if(nir_l!=iwnir_l){
            printf("aer_wave_base is different from the aer_wave_long, which may break the uncertainty for mbac\n");
            exit(1);
        }

        ntau_870=aertab->ntau_870;
        npc=aertab->npc;
        rhoa_modl=(float **)malloc(aertab->ntau_870*sizeof(float*));
        pc_rhoa=(float **)malloc(aertab->ntau_870*sizeof(float *));

        for(itau=0;itau<aertab->ntau_870;itau++){
            pc_rhoa[itau]=(float *)malloc(aertab->npc*sizeof(float));
            rhoa_modl[itau]=(float *)malloc(nwave*sizeof(float));
        }
        rho_all_wav_pred = (float*) malloc(nwave * sizeof (float));
        rho_pred_min = (float*) malloc(nwave * sizeof (float));
        rho_pred_max = (float*) malloc(nwave * sizeof (float));
        chi_old=(float*) malloc((nmodels) * sizeof (float));
        tau_iwnir_l=(float*) malloc((nmodels) * sizeof (float));

        nbands_ac=input->nbands_ac;
        acbands_index=input->acbands_index;

        float *tempwave=(float *)malloc(nbands_ac*sizeof(float));
        for(iw=0;iw<nbands_ac;iw++)
            tempwave[iw]=input->mbac_wave[iw];

        wave_base_index=windex(input->aer_wave_base*1.,tempwave,nbands_ac);
        free(tempwave);

        if(uncertainty){
            derv_chi_taua_l =(float *)malloc(nmodels*sizeof(float));
            derv_chi_rhorc  =(float **)malloc(nmodels*sizeof(float *));
            derv_rhoa_rhorc_l=(float **)malloc(nmodels*sizeof(float*));
            derv_taua_rhorc_l=(float *)malloc(nmodels*sizeof(float));

            for(im=0;im<nmodels;im++){
                derv_chi_rhorc[im]   =(float *)malloc(nbands_ac*sizeof(float));
                derv_rhoa_rhorc_l[im]=(float *)malloc(nwave*sizeof(float));
            }
            derv_chi_obs_rhorc=(float *)malloc(nbands_ac*sizeof(float));
            derv_rhoa_min=(float *)malloc(nwave*sizeof(float));
            derv_rhoa_max=(float *)malloc(nwave*sizeof(float));
            derv_taua_min=(float *)malloc(nwave*sizeof(float));
            derv_taua_max=(float *)malloc(nwave*sizeof(float));
        }
    }

    if(uncertainty){
        derv_modrat_rhorc=uncertainty->derv_modrat_rhorc;
        derv_taua_min_rhorc_l=uncertainty->derv_taua_min_rhorc_l;
        derv_taua_min_taua_l=uncertainty->derv_taua_min_taua_l;
        derv_taua_min_rhow_l=uncertainty->derv_taua_min_rhow_l;
        derv_taua_max_rhorc_l=uncertainty->derv_taua_max_rhorc_l;
        derv_taua_max_taua_l=uncertainty->derv_taua_max_taua_l;
        derv_taua_max_rhow_l=uncertainty->derv_taua_max_rhow_l;
    }

    float chi_obs=0.0;
    if(uncertainty){
        for(iw = 0; iw <nbands_ac; iw++)
            derv_chi_obs_rhorc[iw]=0.;
        derv_chi_obs_taua_l=0.;
    }

    for (ib = 0; ib <nbands_ac; ib++) {

        iw=acbands_index[ib];
        if(iw==nir_l)
            continue;
        chi_obs += rhoa[iw]/rhoa[nir_l];
        mbac_wsum += mbac_w[iw];
        if(uncertainty){
            derv_chi_obs_rhorc[ib]=1/rhoa[nir_l];
            derv_chi_obs_rhorc[wave_base_index]+=(-rhoa[iw]/rhoa[nir_l]/rhoa[nir_l]);
            derv_chi_obs_taua_l+=(-1/rhoa[nir_l]*uncertainty->derv_Lg_taua[iw]+rhoa[iw]/rhoa[nir_l]/rhoa[nir_l]*uncertainty->derv_Lg_taua[nir_l]);
        }
    }
    chi_obs/=mbac_wsum;
    if(uncertainty){
        derv_chi_obs_taua_l/=mbac_wsum;
        for (ib = 0; ib <nbands_ac; ib++)
            derv_chi_obs_rhorc[ib]/=mbac_wsum;
    }

    for (im = 0; im < nmodels; im++) {

        modl = mindx[im];
        aermod=aertab->model[modl];
        //ext_iwnir_l=aermod->extc[nir_l];

        get_pc_rhoa(modl,geom,pc_rhoa);

        for(itau=0;itau<ntau_870;itau++)
            for (ib = 0; ib <nbands_ac; ib++){

                iw=acbands_index[ib];
                rhoa_modl[itau][iw]=0;

                for(ipc=0;ipc<npc;ipc++)
                    rhoa_modl[itau][iw]+=pc_rhoa[itau][ipc]*aermod->pc_components_rhoa[ipc][iw];

                rhoa_modl[itau][iw]+=aermod->pc_mean_rhoa[iw];
                //rhoa_modl[itau][iw]=exp(rhoa_modl[itau][iw]);
            }


        /* compute model epsilon */
        if(rhoa[nir_l]<=rhoa_modl[0][nir_l])
            itau=1;
        else if (rhoa[nir_l]>= rhoa_modl[ntau_870-1][nir_l])
            itau=ntau_870-1;
        else{
            for(itau=1;itau<ntau_870;itau++){
                if(rhoa[nir_l]<rhoa_modl[itau][nir_l] && rhoa[nir_l]>=rhoa_modl[itau-1][nir_l])
                    break;
            }
        }
        tau_iwnir_l[im]=aermod->tau_870[itau]+(aermod->tau_870[itau]-aermod->tau_870[itau-1])*(rhoa[nir_l]-rhoa_modl[itau][nir_l])/(rhoa_modl[itau][nir_l]-rhoa_modl[itau-1][nir_l]);

        if(nmodels<=1){
            if(comp_rhoa_pca(nwave, wave, geom, tau_iwnir_l[im], modl, tau_pred_min, rhoa,derv_rhoa_min,derv_taua_min) !=0)
                return -1;
            for(ib=0;ib<nwave;ib++)
                tau_pred_max[ib]=tau_pred_min[ib];
            *modmin=modl;*modmax=modl;*modrat=1.0;
            return 0;
        }

        if(uncertainty)
            derv_taua_rhorc_l[im]=(aermod->tau_870[itau]-aermod->tau_870[itau-1])/(rhoa_modl[itau][nir_l]-rhoa_modl[itau-1][nir_l]);

        /* compute reflectance at all wavelength */
        for (ib = 0; ib <nbands_ac; ib++){
            iw=acbands_index[ib];
            rho_all_wav_pred[iw]=rhoa_modl[itau][iw]+(rhoa_modl[itau][iw]-rhoa_modl[itau-1][iw])*(tau_iwnir_l[im]-aermod->tau_870[itau])/(aermod->tau_870[itau]-aermod->tau_870[itau-1]);

          //  ext_coef = aermod->extc[iw];
            if(uncertainty){
                derv_rhoa_rhorc_l[im][iw]=derv_taua_rhorc_l[im]*(rhoa_modl[itau][iw]-rhoa_modl[itau-1][iw])/(aermod->tau_870[itau]-aermod->tau_870[itau-1]);
               // derv_taua_rhorc_l[im][iw]=(ext_coef / ext_iwnir_l)*derv_temp_rhorc;
            }
        }
        if(uncertainty){
            derv_chi_rhorc[im][wave_base_index]=0.;
            derv_chi_taua_l[im]=0.;
        }
        mbac_wsum=0.;
        for (ib = 0; ib <nbands_ac; ib++){
            iw=acbands_index[ib];
            if(iw==nir_l)
                continue;

            diff_2 += rho_all_wav_pred[iw]/rhoa[nir_l];
            diff_1 +=pow((rhoa[iw] - rho_all_wav_pred[iw])/noise[iw], 2)*mbac_w[iw];//
            mbac_wsum += mbac_w[iw];

            if(uncertainty){
                derv_temp_rhorc=(derv_rhoa_rhorc_l[im][iw]/rhoa[nir_l]-rho_all_wav_pred[iw]/rhoa[nir_l]/rhoa[nir_l]);
                derv_chi_rhorc [im][wave_base_index]+=derv_temp_rhorc;
                derv_chi_taua_l[im]+=(-derv_temp_rhorc*uncertainty->derv_Lg_taua[nir_l]);
            }
        }

        chi[im] = diff_2 / mbac_wsum;
        chi_old[im]=diff_1/mbac_wsum;
        if(uncertainty){
            derv_chi_rhorc[im][wave_base_index]/=mbac_wsum;
            derv_chi_taua_l [im]/=mbac_wsum;
        }

        diff_2 = 0.0;
        diff_1 = 0.0;
        mbac_wsum = 0.0;
        chi_struct[im].value=chi_old[im];
        chi_struct[im].index=im;
    }

    qsort(chi_struct,nmodels,sizeof(chi_struct[0]),cmp);

    /* for (im = 0; im < nmodels; im++)
           if(chi_obs<chi_struct[im].value)
               break;
     *modmin = MAX(MIN(im - 1, nmodels - 1), 0);
     *modmax = MAX(MIN(im, nmodels - 1), 0);*/

    *modrat=chi_struct[0].value/(chi_struct[0].value+chi_struct[1].value);
    chi_old[0]=chi_struct[0].value*(1-*modrat)+chi_struct[1].value*(*modrat);  //temporary keeping the chi2 value from two closest aerosol models

    /*chi_struc is switched back to the ratio after selecting the model using the chi based on difference*/
    for(im=0;im<nmodels;im++)
        chi_struct[im].value=chi[chi_struct[im].index];


    /*if( (chi_obs-chi_struct[0].value)*(chi_obs-chi_struct[1].value)>0){
     *modrat=1.0;    //need to flag
     *modmin=chi_struct[0].index;
     *modmax=*modmin;
       }
       else*/{
        if(chi_obs>chi_struct[1].value){
            im=chi_struct[0].index;
            diff_1=chi_struct[0].value;

            chi_struct[0].value=chi_struct[1].value;
            chi_struct[0].index=chi_struct[1].index;

            chi_struct[1].value=diff_1;
            chi_struct[1].index=im;
        }
        *modmin=chi_struct[0].index;
        *modmax=chi_struct[1].index;

        *modrat=(chi_obs-chi_struct[0].value)/(chi_struct[1].value-chi_struct[0].value);
        if(uncertainty){
            derv_modrat_chi0= (chi_obs-chi_struct[1].value)/pow(chi_struct[1].value-chi_struct[0].value,2);
            derv_modrat_chi1=-(chi_obs-chi_struct[0].value)/pow(chi_struct[1].value-chi_struct[0].value,2);

            for (ib = 0; ib <nbands_ac; ib++){
                iw=acbands_index[ib];
                if(iw!=nir_l){
                    derv_modrat_rhorc[ib]=derv_chi_obs_rhorc[ib]/(chi_struct[1].value-chi_struct[0].value);
                    derv_modrat_rhow_l+=-derv_modrat_rhorc[ib]*uncertainty->ratio_rhow[ib];
                }
            }

            derv_modrat_rhorc[wave_base_index]=derv_modrat_chi0*derv_chi_rhorc[*modmin][wave_base_index]+derv_modrat_chi1*derv_chi_rhorc[*modmax][wave_base_index];
            derv_modrat_rhorc[wave_base_index]+=derv_chi_obs_rhorc[wave_base_index]/(chi_struct[1].value-chi_struct[0].value);
            derv_modrat_rhow_l+=-derv_modrat_rhorc[wave_base_index];

            derv_modrat_taua_l =derv_modrat_chi0*derv_chi_taua_l [*modmin]+derv_modrat_chi1*derv_chi_taua_l [*modmax];
            derv_modrat_taua_l +=derv_chi_obs_taua_l/(chi_struct[1].value-chi_struct[0].value);

            uncertainty->derv_modrat_rhow_l=derv_modrat_rhow_l;
            uncertainty->derv_modrat_taua_l=derv_modrat_taua_l;
        }
       }

       if(comp_rhoa_pca(nwave, wave, geom, tau_iwnir_l[*modmin], mindx[*modmin], tau_pred_min, rho_pred_min,derv_rhoa_min,derv_taua_min) !=0)
           return -1;
       if(comp_rhoa_pca(nwave, wave, geom, tau_iwnir_l[*modmax], mindx[*modmax], tau_pred_max, rho_pred_max,derv_rhoa_max,derv_taua_max) !=0)
           return -1;

       for (iw = 0; iw <nwave; iw++) {

           rhoa[iw]=rho_pred_min[iw]*(1-*modrat)+rho_pred_max[iw]*(*modrat);

           if(uncertainty){
               for(itemp=0;itemp<nbands_ac;itemp++){
                   ib=acbands_index[itemp];
                   if(ib!=nir_l){
                       uncertainty->derv_La_rhorc[iw][itemp]=(rho_pred_max[iw]-rho_pred_min[iw])*derv_modrat_rhorc[itemp];
                       uncertainty->derv_La_rhow_l[iw]+=-uncertainty->derv_La_rhorc[iw][itemp]*uncertainty->ratio_rhow[itemp];
                   }
               }
               uncertainty->derv_La_rhorc[iw][wave_base_index]=(rho_pred_max[iw]-rho_pred_min[iw])*derv_modrat_rhorc[wave_base_index];
               uncertainty->derv_La_rhorc[iw][wave_base_index]+= (1-*modrat)*derv_rhoa_min[iw]*derv_taua_rhorc_l[*modmin];
               uncertainty->derv_La_rhorc[iw][wave_base_index]+= (*modrat)  *derv_rhoa_max[iw]*derv_taua_rhorc_l[*modmax];
               uncertainty->derv_La_rhow_l[iw]+=-uncertainty->derv_La_rhorc[iw][wave_base_index];

               uncertainty->derv_La_taua_l[iw]=(rho_pred_max[iw]-rho_pred_min[iw])*derv_modrat_taua_l;
               uncertainty->derv_La_taua_l[iw]+=(1-*modrat)*(-derv_rhoa_min[iw]*derv_taua_rhorc_l[*modmin]*uncertainty->derv_Lg_taua[iw]);
               uncertainty->derv_La_taua_l[iw]+=(*modrat)  *(-derv_rhoa_max[iw]*derv_taua_rhorc_l[*modmax]*uncertainty->derv_Lg_taua[iw]);

               derv_taua_min_rhorc_l[iw]=derv_taua_min[iw]*derv_taua_rhorc_l[*modmin];
               derv_taua_max_rhorc_l[iw]=derv_taua_max[iw]*derv_taua_rhorc_l[*modmax];
               derv_taua_min_taua_l[iw]=-derv_taua_min_rhorc_l[iw]*uncertainty->derv_Lg_taua[iw];
               derv_taua_max_taua_l[iw]=-derv_taua_max_rhorc_l[iw]*uncertainty->derv_Lg_taua[iw];
               derv_taua_min_rhow_l[iw]=-derv_taua_min_rhorc_l[iw];
               derv_taua_max_rhow_l[iw]=-derv_taua_max_rhorc_l[iw];

           }
       }

       *modmin=mindx[*modmin];
       *modmax=mindx[*modmax];
       chi[0]=chi_old[0];



       return (status);
}



/* ---------------------------------------------------------------------------------------- */
/* smaer() - compute aerosol reflectance using MSEPS approach of Ahmad
 * perform spectral matching in reflectance space                                           */
/* output the reflectance of the two bracketing aerosol models                              */
/*
 * Amir Ibrahim, Jul 2016.                                                                  */

/* ---------------------------------------------------------------------------------------- */

int smaer(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s,
        int32_t iwnir_l, int32_t nmodels, int32_t mindx[], geom_str *geom, float wv, float rhoa[], int32_t *modmin, int32_t *modmax,
        float32 *modrat, float *tau_pred_min, float *tau_pred_max, int32_t ip, float chi[], float *mbac_w,uncertainty_t *uncertainty) {

    if (!have_ms && !use_pca_lut) {
        printf("\nThe multi-scattering spectral matching atmospheric correction method requires\n");
        printf("ams_all, bms_all, cms_all, dms_all, ems in the aerosol model tables.\n");
        exit(1);
    }

    float *ac, *bc, *cc, *dc, *ec;
    float ext_iwnir_l,tau_iwnir_l;
    double ax, bx, cx, fx;
    static int firstcall=1, nir_l;

    int iw, i;

    int status = 0;

    if(use_pca_lut){
        status=smaer_pca(sensorID,wave,nwave,iwnir_s,iwnir_l,nmodels,mindx,geom,wv,rhoa,
                modmin,modmax,modrat,tau_pred_min,tau_pred_max,ip,chi,mbac_w,uncertainty);
        return status;
    }

    static float **tau_all_wav, **rho_all_wav_pred;

    float ext_coef, lg_taua;

    int iwtab, im, modl;

    float diff_2 = 0.0;
    float mbac_wsum = 0.0;

    float *noise;
    int iwtab_l;
    int nbands_ac=iwnir_l-iwnir_s+1;

    float **derv_chi_rhorc;//dimension: [nmodels][nbands_ac], deivative of chi to rhorc
    float *derv_chi_taua_l; //dimension: [nmodels], deivative of chi to taua_l
    float **derv_rhoa_rhorc_l;// derivative of modeled rhoa[nwave] to rhorc[iwnir_l]
    float *derv_modrat_rhorc;// dimension: [nbands_ac], deivative of modrat to rhorc[iwnir_s to iwnir_l]
    float derv_modrat_taua_l=0.; // deivative of modrat to taua [iwnir_l]
    float derv_modrat_rhow_l=0.;
    float *derv_taua_min_rhorc_l;
    float *derv_taua_min_taua_l;
    float *derv_taua_min_rhow_l;
    float *derv_taua_max_rhorc_l;
    float *derv_taua_max_taua_l;
    float *derv_taua_max_rhow_l;
    float **derv_taua_rhorc_l;//[nmodels][nwave], derivative of modeled taua[nwave] to rhorc_l

    float derv_temp_rhorc;
    float derv_modrat_chi0,derv_modrat_chi1;

    if(uncertainty){
        derv_modrat_rhorc=uncertainty->derv_modrat_rhorc;
        derv_taua_min_rhorc_l=uncertainty->derv_taua_min_rhorc_l;
        derv_taua_min_taua_l=uncertainty->derv_taua_min_taua_l;
        derv_taua_min_rhow_l=uncertainty->derv_taua_min_rhow_l;
        derv_taua_max_rhorc_l=uncertainty->derv_taua_max_rhorc_l;
        derv_taua_max_taua_l=uncertainty->derv_taua_max_taua_l;
        derv_taua_max_rhow_l=uncertainty->derv_taua_max_rhow_l;

        derv_chi_taua_l =(float *)malloc(nmodels*sizeof(float));
        derv_chi_rhorc  =(float **)malloc(nmodels*sizeof(float *));
        derv_rhoa_rhorc_l=(float **)malloc(nmodels*sizeof(float*));
        derv_taua_rhorc_l=(float **)malloc(nmodels*sizeof(float*));

        for(im=0;im<nmodels;im++){
            derv_chi_rhorc[im]   =(float *)malloc(nbands_ac*sizeof(float));
            derv_rhoa_rhorc_l[im]=(float *)malloc(nwave*sizeof(float));
            derv_taua_rhorc_l[im]=(float *)malloc(nwave*sizeof(float));
        }
        for(im=0;im<nmodels;im++)
            for(iw=0;iw<nbands_ac;iw++)
                derv_chi_rhorc[im][iw]=0.;
    }



    if(firstcall){
        firstcall=0;
        nir_l=bindex_get(input->aer_wave_base);

        tau_all_wav = (float**) malloc(nwave * sizeof (float*));
        rho_all_wav_pred = (float**) malloc(nwave * sizeof (float*));

        for (i = 0; i < nwave; i++) {
            tau_all_wav[i] = (float*) malloc((nmodels) * sizeof (float));
            rho_all_wav_pred[i] = (float*) malloc((nmodels) * sizeof (float));
        }
    }

    noise= &noise_global[ip*nwave];

    // sorting in ascending way elements of chi-squared --  will need clean-up
    struct str { float value;int index;};
    struct str chi_struct[nmodels];


    if(sensorID==OCIS || sensorID==OCI){
        //nir_l=bindex_get(870);
       /* for(iw=iwnir_s;iw<=iwnir_l;iw++){
            if( (wave[iw]>750 && wave[iw]<865) )
                mbac_w[iw]=1.0;
            //if( (iw>=229 && iw<=233) || (iw>=235 &&iw<=237))
               // mbac_w[iw]=0.;

        }*/
    }
    iwtab_l=iwatab[nir_l];


    for (im = 0; im < nmodels; im++) {

        modl = mindx[im];

        ms_eps_coef(modl, nir_l, wave, geom, &ac, &bc, &cc, &dc, &ec);

        ax = (double) ac[iwtab_l] - log((double) rhoa[nir_l]);
        bx = (double) bc[iwtab_l];
        cx = (double) cc[iwtab_l];

        fx = bx * bx - 4.0 * ax*cx;
        if (fx > 0.0 && cx != 0.0) {
            tau_iwnir_l = 0.5 * (-bx + sqrt(fx)) / cx;
            if(uncertainty)
                derv_temp_rhorc= 1 / rhoa[nir_l] / sqrt(fx);

            tau_iwnir_l = exp(tau_iwnir_l);
            if(uncertainty)
                derv_temp_rhorc *= tau_iwnir_l;
        } else {
            status = 1;
            return(status);
            //break;
        }

        /* compute reflectance at all wavelength */

        ext_iwnir_l = aertab->model[modl]->extc[iwtab_l];
        for (iw = 0; iw <nwave ; iw++) {
            iwtab = iwatab[iw];
            ext_coef = aertab->model[modl]->extc[iwtab];
            tau_all_wav[iw][im] = (ext_coef / ext_iwnir_l) * tau_iwnir_l;
            lg_taua = log(tau_all_wav[iw][im]);
            if(uncertainty){
                derv_taua_rhorc_l[im][iw]=(ext_coef / ext_iwnir_l)*derv_temp_rhorc;
                derv_rhoa_rhorc_l[im][iw]=1/tau_all_wav[iw][im]*(ext_coef / ext_iwnir_l)*derv_temp_rhorc;
            }

            lg_taua = ac[iwtab] + bc[iwtab] * lg_taua +
                    cc[iwtab] * pow(lg_taua, 2) +
                    dc[iwtab] * pow(lg_taua, 3) +
                    ec[iwtab] * pow(lg_taua, 4);

            if(uncertainty)
                derv_rhoa_rhorc_l[im][iw] *=( bc[iwtab]+ 2*cc[iwtab]*lg_taua+3*dc[iwtab]*pow(lg_taua, 2)+4*ec[iwtab] * pow(lg_taua, 3) );

            rho_all_wav_pred[iw][im] = exp(lg_taua);

            if(uncertainty){
                derv_rhoa_rhorc_l[im][iw]*= rho_all_wav_pred[iw][im];
               // derv_rhoa_taua_l [im][iw] =
            }
        }

        if(uncertainty){
            derv_chi_rhorc[im][iwnir_l-iwnir_s]=0.;
            derv_chi_taua_l[im]=0.;
        }
        for (iw = iwnir_s; iw <= iwnir_l; iw++) {
            if(mbac_w[iw]==0.|| iw==nir_l)
                continue;
            //noise[iw]=uncertainty->dsensor[ip*nwave+iw];
            diff_2 += pow((rhoa[iw] - rho_all_wav_pred[iw][im])/noise[iw], 2)*mbac_w[iw];//
            mbac_wsum += mbac_w[iw];

            if(uncertainty){
                /*derv_chi_rhorc[im][iw-iwnir_s]=2*fabs(rhoa[iw] - rho_all_wav_pred[iw][im])/noise[iw]/noise[iw]*mbac_w[iw];
                  derv_chi_rhorc [im][iwnir_l-iwnir_s]+=fabs(-derv_chi_rhorc[im][iw-iwnir_s])*derv_rhoa_rhorc_l[im][iw];*/
                if(rhoa[iw] >=rho_all_wav_pred[iw][im]){
                    derv_chi_rhorc[im][iw-iwnir_s]=2*(rhoa[iw] - rho_all_wav_pred[iw][im])/noise[iw]/noise[iw]*mbac_w[iw];
                    derv_chi_rhorc [im][iwnir_l-iwnir_s]+=(-derv_chi_rhorc[im][iw-iwnir_s])*derv_rhoa_rhorc_l[im][iw];
                    derv_chi_taua_l[im]+=(derv_chi_rhorc[im][iw-iwnir_s])*(-uncertainty->derv_Lg_taua[iw]+ derv_rhoa_rhorc_l[im][iw]*uncertainty->derv_Lg_taua[iwnir_l]);
                }
                else {
                    derv_chi_rhorc[im][iw-iwnir_s]=2*(-rhoa[iw] + rho_all_wav_pred[iw][im])/noise[iw]/noise[iw]*mbac_w[iw];
                    derv_chi_rhorc [im][iwnir_l-iwnir_s]+=(-derv_chi_rhorc[im][iw-iwnir_s])*derv_rhoa_rhorc_l[im][iw];
                    derv_chi_taua_l[im]+=(derv_chi_rhorc[im][iw-iwnir_s])*(uncertainty->derv_Lg_taua[iw]-derv_rhoa_rhorc_l[im][iw]*uncertainty->derv_Lg_taua[iwnir_l]);
                }
            }
        }

        chi[im] = (diff_2 / mbac_wsum);
        if(uncertainty){
            for (iw = iwnir_s; iw <= iwnir_l; iw++){
                derv_chi_rhorc[im][iw-iwnir_s]/=mbac_wsum;
                derv_chi_taua_l [im]/=mbac_wsum;
            }
        }


        diff_2 = 0.0;
        mbac_wsum = 0.0;
        chi_struct[im].value=chi[im];
        chi_struct[im].index=im;

    }

    qsort(chi_struct,nmodels,sizeof(chi_struct[0]),cmp);

    *modmin=chi_struct[0].index;
    *modmax=chi_struct[1].index;
   // *modrat=1/chi_struct[1].value/(1/chi_struct[0].value+1/chi_struct[1].value); //in case one of the chi=0.0;
    *modrat=chi_struct[0].value/(chi_struct[0].value+chi_struct[1].value);
    //*modrat=-chi_struct[0].value/(chi_struct[1].value-chi_struct[0].value);
    chi[0]=chi_struct[0].value*(1-*modrat)+chi_struct[1].value*(*modrat);  //temporary keeping the chi2 value from two closest aerosol models

    if(uncertainty){
        derv_modrat_chi0= chi_struct[1].value/pow(chi_struct[0].value+chi_struct[1].value,2);
        derv_modrat_chi1=-chi_struct[0].value/pow(chi_struct[0].value+chi_struct[1].value,2);

        //derv_modrat_chi0= -chi_struct[1].value/pow(chi_struct[0].value-chi_struct[1].value,2);
        //derv_modrat_chi1=  chi_struct[0].value/pow(chi_struct[0].value-chi_struct[1].value,2);

        for(iw=iwnir_s;iw<=iwnir_l;iw++){
            derv_modrat_rhorc[iw-iwnir_s]=derv_modrat_chi0*derv_chi_rhorc[*modmin][iw-iwnir_s]+derv_modrat_chi1*derv_chi_rhorc[*modmax][iw-iwnir_s];
            derv_modrat_rhow_l+=-derv_modrat_rhorc[iw-iwnir_s]*uncertainty->ratio_rhow[iw-iwnir_s];
        }
        derv_modrat_taua_l =derv_modrat_chi0*derv_chi_taua_l [*modmin]+derv_modrat_chi1*derv_chi_taua_l [*modmax];

        uncertainty->derv_modrat_rhow_l=derv_modrat_rhow_l;
        uncertainty->derv_modrat_taua_l=derv_modrat_taua_l;
    }

    for (iw = 0; iw <nwave; iw++) {

        tau_pred_min[iw] = tau_all_wav[iw][*modmin];
        tau_pred_max[iw] = tau_all_wav[iw][*modmax];
        rhoa[iw]=rho_all_wav_pred[iw][*modmin]*(1-*modrat)+rho_all_wav_pred[iw][*modmax]*(*modrat);

        if(uncertainty){
            for(i=iwnir_s;i<iwnir_l;i++){
                uncertainty->derv_La_rhorc[iw][i-iwnir_s]=(rho_all_wav_pred[iw][*modmax]-rho_all_wav_pred[iw][*modmin])*derv_modrat_rhorc[i-iwnir_s];
                uncertainty->derv_La_rhow_l[iw]+=-uncertainty->derv_La_rhorc[iw][i-iwnir_s]*uncertainty->ratio_rhow[i-iwnir_s];

            }
            uncertainty->derv_La_rhorc[iw][iwnir_l-iwnir_s]=(rho_all_wav_pred[iw][*modmax]-rho_all_wav_pred[iw][*modmin])*derv_modrat_rhorc[iwnir_l-iwnir_s];
            uncertainty->derv_La_rhorc[iw][iwnir_l-iwnir_s]+= (1-*modrat)*derv_rhoa_rhorc_l[*modmin][iw];
            uncertainty->derv_La_rhorc[iw][iwnir_l-iwnir_s]+= (*modrat)  *derv_rhoa_rhorc_l[*modmax][iw];
            uncertainty->derv_La_rhow_l[iw]+=-uncertainty->derv_La_rhorc[iw][iwnir_l-iwnir_s];

            uncertainty->derv_La_taua_l[iw]=(rho_all_wav_pred[iw][*modmax]-rho_all_wav_pred[iw][*modmin])*derv_modrat_taua_l;
            uncertainty->derv_La_taua_l[iw]+=(1-*modrat)*(-derv_rhoa_rhorc_l[*modmin][iw]*uncertainty->derv_Lg_taua[iw]);
            uncertainty->derv_La_taua_l[iw]+=(*modrat)  *(-derv_rhoa_rhorc_l[*modmax][iw]*uncertainty->derv_Lg_taua[iw]);

            derv_taua_min_rhorc_l[iw]=derv_taua_rhorc_l[*modmin][iw];
            derv_taua_max_rhorc_l[iw]=derv_taua_rhorc_l[*modmax][iw];
            derv_taua_min_taua_l[iw]=-derv_taua_rhorc_l[*modmin][iw]*uncertainty->derv_Lg_taua[iw];
            derv_taua_max_taua_l[iw]=-derv_taua_rhorc_l[*modmax][iw]*uncertainty->derv_Lg_taua[iw];
            derv_taua_min_rhow_l[iw]=-derv_taua_rhorc_l[*modmin][iw];
            derv_taua_max_rhow_l[iw]=-derv_taua_rhorc_l[*modmax][iw];

        }
    }

    *modmin=mindx[chi_struct[0].index];
    *modmax=mindx[chi_struct[1].index];

   /* for (i = 0; i < nwave; i++) {
        free(tau_all_wav[i]);
        free(rho_all_wav_pred[i]);
    }
    free(tau_all_wav);
    free(rho_all_wav_pred);*/

    if(uncertainty){
        for(im=0;im<nmodels;im++){
            free(derv_chi_rhorc[im]);
            free(derv_rhoa_rhorc_l[im]);
            free(derv_taua_rhorc_l[im]);
        }
        free(derv_chi_rhorc);
        free(derv_chi_taua_l);
        free(derv_rhoa_rhorc_l);
        free(derv_taua_rhorc_l);
    }

    return (status);
}

/* ---------------------------------------------------------------------------------------- */
/* ahmadaer() - compute aerosol reflectance using MSEPS approach of Ahmad                   */
/*                                                                                          */
/* Z. Ahmad, August 2014.                                                                   */
/* M. Zhang,    SAIC,  July 2021,  adding the uncertainty propagation                       */

/* ---------------------------------------------------------------------------------------- */
int ahmadaer(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s, int32_t iwnir_l,
        int32_t nmodels, int32_t mindx[],
        geom_str *geom, float wv, float rhoa[],
        int32_t *modmin, int32_t *modmax, float *modrat, float *epsnir,
        float tau_pred_min[], float tau_pred_max[],int ip,uncertainty_t *uncertainty) {
    int iw;

    float rho_pred_min[nwave], rho_pred_max[nwave];
    float rho_aer[nwave], tau_aer[nwave];

    if (!have_ms && !use_pca_lut) {
        printf("\nThe multi-scattering epsilon atmospheric correction method requires\n");
        printf("ams_all, bms_all, cms_all, dms_all, ems in the aerosol model tables.\n");
        exit(1);
    }
    if (aer_opt == AERRHMSEPS_lin) {
        if (ahmad_atm_corr_lin(sensorID, wave, nwave, iwnir_s, iwnir_l, nmodels, mindx, geom, wv,
                rhoa, modmin, modmax, modrat, epsnir,
                tau_pred_max, tau_pred_min, rho_pred_max, rho_pred_min, tau_aer, rho_aer) != 0)
            return (1);
    } else {
        /* use the ms_epsilon method to get rhoa */
        if (ahmad_atm_corr(sensorID, wave, nwave, iwnir_s, iwnir_l, nmodels, mindx, geom, wv,
                rhoa, modmin, modmax, modrat, epsnir,
                tau_pred_max, tau_pred_min, rho_pred_max, rho_pred_min, tau_aer, rho_aer,ip,uncertainty) != 0)
            return (1);
    }

    for (iw = 0; iw < nwave; iw++) {
        rhoa[iw] = rho_aer[iw];
    }

    return (0);
}


/* ---------------------------------------------------------------------------------------- */
/* wangaer() - compute aerosol reflectance using Gordon & Wang 1994 algorithm               */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */

/* ---------------------------------------------------------------------------------------- */

int wangaer(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s,
        int32_t iwnir_l, int32_t nmodels, int32_t mindx[], geom_str *geom,
        float wv, float rhoa[], int32_t *modmin, int32_t *modmax,
        float *modrat, float *epsnir, float tauamin[], float tauamax[]) {
    int modflg;
    float *epsmin;
    float *epsmax;

    float epsmin1;
    float epsmax1;

    float *rhoasmin;
    float *rhoasmax;
    float *rhoamin;
    float *rhoamax;

    float cc = 0.0;
    int iw;

    if ((rhoasmin = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoasmin.\n");
        exit(1);
    }
    if ((rhoasmax = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoasmax.\n");
        exit(1);
    }
    if ((rhoamin = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoamin.\n");
        exit(1);
    }
    if ((rhoamax = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoasmax.\n");
        exit(1);
    }

    /* find upper and lower-bounding models */
    modflg = model_select_wang(sensorID, wave, nwave, nmodels, mindx, geom, wv,
            rhoa, iwnir_s, iwnir_l, modmin, modmax, modrat, epsnir);

    /* if no lower-bounding aerosol model, set-up for extrapolation */
    if (modflg < 0)
        cc = log(*epsnir) / (wave[iwnir_l] - wave[iwnir_s]);

    /* get model epsilon for each bounding model, all wavelengths */
    epsmin = model_epsilon(*modmin, iwnir_l, wave, nwave, geom);
    epsmax = model_epsilon(*modmax, iwnir_l, wave, nwave, geom);

    /* get SS aerosol reflectance at longest wavelength for the two models */
    if (rhoa_to_rhoas(sensorID, *modmin, geom, wv, rhoa, wave, nwave, iwnir_l, iwnir_l, rhoasmin) != 0) {
        free(rhoamin);
        free(rhoasmin);
        free(rhoamax);
        free(rhoasmax);
        return (1);
    }
    if (rhoa_to_rhoas(sensorID, *modmax, geom, wv, rhoa, wave, nwave, iwnir_l, iwnir_l, rhoasmax) != 0) {
        free(rhoamin);
        free(rhoasmin);
        free(rhoamax);
        free(rhoasmax);
        return (1);
    }
    /* compute SS aerosol reflectance in all bands */
    for (iw = 0; iw < nwave; iw++) {

        epsmin1 = epsmin[iw];
        epsmax1 = epsmax[iw];

        if (modflg < 0) {
            epsmin1 = exp(cc * (wave[iwnir_l] - wave[iw]));
            epsmax1 = epsmin1;
        }

        rhoasmin[iw] = rhoasmin[iwnir_l] * epsmin1;
        rhoasmax[iw] = rhoasmax[iwnir_l] * epsmax1;
    }

    /* compute MS aerosol reflectance in visible bands */
    rhoas_to_rhoa(sensorID, *modmin, geom, wv, rhoasmin, wave, nwave, 0, nwave - 1, rhoamin);
    rhoas_to_rhoa(sensorID, *modmax, geom, wv, rhoasmax, wave, nwave, 0, nwave - 1, rhoamax);

    /* interpolate between upper and lower-bounding models */
    for (iw = 0; iw < nwave; iw++) {
        rhoa[iw] = rhoamin[iw]*(1.0 - (*modrat)) + rhoamax[iw]*(*modrat);
    }

    model_taua(sensorID, *modmin, wave, nwave, iwnir_l, rhoa, geom, wv, tauamin);
    model_taua(sensorID, *modmax, wave, nwave, iwnir_l, rhoa, geom, wv, tauamax);

    free(rhoamin);
    free(rhoasmin);
    free(rhoamax);
    free(rhoasmax);

    return (0);
}

/* ---------------------------------------------------------------------------------------- */
/* model_select_franz() - Franz aerosol model selection process.                            */
/*                                                                                          */
/* B. Franz, 1 February 2009.                                                               */

/* ---------------------------------------------------------------------------------------- */

typedef struct rhoaT_struct {
    int32_t modnum;
    float rhoa;
    float eps;
} rhoaTstr;

int comp_rhoaT(rhoaTstr *x, rhoaTstr *y) {
    return (x->rhoa < y->rhoa ? -1 : 1);
}

int model_select_franz(int32_t sensorID, float wave[], int32_t nwave,
        int32_t nmodel, int32_t mindx[], geom_str *geom, float wv,
        float rhoa[], int32_t iwnir_s, int32_t iwnir_l, int32_t *modmin,
        int32_t *modmax, float *modrat, float *epsnir) {

    float *rhoas;
    float *rhoa_tmp;
    rhoaTstr rhoa_tab[MAXMODEL];

    float *eps;
    int jm, im;
    int jm1, jm2;
    float wt;

    *modmin = -1;
    *modmax = -1;
    *modrat = 0.0;
    *epsnir = 0.0;

    if ((rhoas = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoas.\n");
        exit(1);
    }
    if ((rhoa_tmp = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoas_tmp.\n");
        exit(1);
    }

    // predict MS aerosol reflectance assuming each model

    for (jm = 0; jm < nmodel; jm++) {

        im = mindx[jm];

        // get model epsilon at this geometry
        eps = model_epsilon(im, iwnir_l, wave, nwave, geom);

        // get SS aerosol reflectance at iwnir_l
        if (rhoa_to_rhoas(sensorID, im, geom, wv, rhoa, wave, nwave, iwnir_l, iwnir_l, rhoas) != 0) {
            free(rhoas);
            free(rhoa_tmp);
            return (1);
        }
        // get SS aerosol reflectance at iwnir_s
        rhoas[iwnir_s] = rhoas[iwnir_l] * eps[iwnir_s];

        // get MS aerosol reflectance at iwnir_s and save
        rhoas_to_rhoa(sensorID, im, geom, wv, rhoas, wave, nwave, iwnir_s, iwnir_s, rhoa_tmp);

        rhoa_tab[jm].modnum = im;
        rhoa_tab[jm].rhoa = rhoa_tmp[iwnir_s];
        rhoa_tab[jm].eps = eps[iwnir_s];
    }

    // put results in ascending order of predicted rhoa[iwnir_s]
    qsort(rhoa_tab, nmodel, sizeof (rhoaTstr), (int (*)(const void *, const void *)) comp_rhoaT);

    // compare observed rhoa with model predictions at iwnir_s to select models
    for (jm = 0; jm < nmodel; jm++) {
        if (rhoa_tab[jm].rhoa > rhoa[iwnir_s])
            break;
    }
    if (jm == 0) {
        jm1 = 0;
        jm2 = 1;
    } else if (jm == nmodel) {
        jm1 = nmodel - 2;
        jm2 = nmodel - 1;
    } else {
        jm1 = jm - 1;
        jm2 = jm1 + 1;
    }
    wt = (rhoa[iwnir_s] - rhoa_tab[jm1].rhoa) / (rhoa_tab[jm2].rhoa - rhoa_tab[jm1].rhoa);

    *modmin = rhoa_tab[jm1].modnum;
    *modmax = rhoa_tab[jm2].modnum;
    *modrat = wt;
    *epsnir = rhoa_tab[jm1].eps * (1.0 - wt) + rhoa_tab[jm2].eps*wt;

    free(rhoas);
    free(rhoa_tmp);

    return (0);
}


/* ---------------------------------------------------------------------------------------- */
/* order_models is the sorting function used in rhaer                                       */

/* ---------------------------------------------------------------------------------------- */
static int order_models(const void *p1, const void *p2) {
    aermodstr *x = *(aermodstr **) p1;
    aermodstr *y = *(aermodstr **) p2;

    if (x->rh == y->rh) {
        if (x->sd > y->sd)
            return ( 1);
        else
            return (-1);
    } else {
        if (x->rh > y->rh)
            return ( 1);
        else
            return (-1);
    }
}


/* ---------------------------------------------------------------------------------------- */
/* rhaer() - compute aerosol reflectance using RH descrimination + desired selection scheme */
/* M. Zhang, SAIC,  July 2021,  adding the uncertainty propagation                          */

/* ---------------------------------------------------------------------------------------- */
int rhaer(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s, int32_t iwnir_l,
        geom_str *geom, float wv, float rh, float pr, float taur[], float rhoa[],
        int32_t *modmin1, int32_t *modmax1, float *modrat1, int32_t *modmin2, int32_t *modmax2, float *modrat2,
        float *eps, float taua[], float tsol[], float tsen[], int32_t ip, float *mbac_w,float *chi, uncertainty_t *uncertainty) {
    static int firstCall = 1;
    static int nrh;
    static float rhtab[MAXAERMOD];
    static int nsd;
    static int sdtab[MAXAERMOD];

    float *rhoa1;
    float *rhoa2;
    float *taua1;
    float *taua2;
    float *tsol1;
    float *tsol2;
    float *tsen1;
    float *tsen2;
    float eps1=1.0;
    float eps2=1.0;
    //float modrat1;
    //float modrat2;
    int nmodels;
    float *chi1, *chi2;

    float *tau_pred_min1;
    float *tau_pred_max1;
    float *tau_pred_min2;
    float *tau_pred_max2;

    int32_t mindx1[MAXAERMOD];
    int32_t mindx2[MAXAERMOD];
    int irh1, irh2, irh;
    //int irh3; // Third RH index --> Amir
    int isd;
    float wt;
    static uncertainty_t *uncertainty1=NULL, *uncertainty2=NULL;
    float derv_wt_rh; //derivative of wt to rh

    int iw, im, inir;

    if (firstCall) {
        firstCall = 0;
        float lastrh = -1.0;
        int lastsd = -1;
        if (!have_rh) {
            printf("-E- %s line %d: This aerosol selection method requires models with a Relative Humidity attribute and size distribution.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        if(uncertainty){
            uncertainty1 = (uncertainty_t *) malloc(sizeof(uncertainty_t));
            uncertainty2 = (uncertainty_t *) malloc(sizeof(uncertainty_t));

            if (alloc_uncertainty(uncertainty->nbands, uncertainty->nbands_ac, uncertainty->npix, uncertainty1) != 0) {
                printf("unable to allocate error record in rhaer()\n");
                exit(1);
            }
            if (alloc_uncertainty(uncertainty->nbands, uncertainty->nbands_ac, uncertainty->npix, uncertainty2) != 0) {
                printf("unable to allocate error record in rhaer()\n");
                exit(1);
            }
        }
        // need in order of rh and sd within rh
        qsort(aertab->model, aertab->nmodel, sizeof (aermodstr*), (int (*)(const void *, const void *)) order_models);

        // count the number of model humidities and the number of model size distributions
        // note that use of a single model suite will yield nrh=1, which inherently avoids RH weighting that case

        nsd = 0;
        nrh = 0;

        for (im = 0; im < aertab->nmodel; im++) {
            if (aertab->model[im]->rh != lastrh) {
                rhtab[nrh] = aertab->model[im]->rh;
                lastrh = rhtab[nrh];
                nrh++;
            }
            if (nrh == 1 && aertab->model[im]->sd != lastsd) {
                sdtab[nsd] = aertab->model[im]->sd;
                lastsd = sdtab[nsd];
                nsd++;
            }
        }
        if (nrh * nsd != aertab->nmodel) {
            printf("-E- %s line %d: number of humidities (%d) x number of size distributions (%d) must equal number of models (%d).\n",
                    __FILE__, __LINE__, nrh, nsd, aertab->nmodel);
            exit(1);
        } else {
            printf("%d aerosol models: %d humidities x %d size fractions\n", aertab->nmodel, nrh, nsd);
            for (irh = 0; irh < nrh; irh++) {
                for (isd = 0; isd < nsd; isd++) {
                    im = irh * nsd + isd;
                    printf("model %d, rh=%f, sd=%d, alpha=%f, name=%s\n",
                            im, aertab->model[im]->rh, aertab->model[im]->sd, aertab->model[im]->angstrom[0], aertab->model[im]->name);
                }
            }
        }
    }

    nmodels=nsd;

    if(uncertainty){
        init_uncertainty(uncertainty1,0);
        init_uncertainty(uncertainty2,0);
        if (cp_uncertainty(uncertainty, uncertainty1, ip) != 0) {
            printf("unable to copy the error record in rhaer()\n");
            exit(1);
        }
        if (cp_uncertainty(uncertainty, uncertainty2, ip) != 0) {
            printf("unable to copy the error record in rhaer()\n");
            exit(1);
        }
    }

    // initialize
    if ((taua1 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for taua1.\n");
        exit(1);
    }
    if ((taua2 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for taua2.\n");
        exit(1);
    }
    if ((tsol1 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsol1.\n");
        exit(1);
    }
    if ((tsol2 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsol2.\n");
        exit(1);
    }
    if ((rhoa1 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoa1.\n");
        exit(1);
    }
    if ((rhoa2 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoa2.\n");
        exit(1);
    }
    if ((tsen1 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsen1.\n");
        exit(1);
    }
    if ((tsen2 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for tsen2.\n");
        exit(1);
    }
    if ((tau_pred_min1 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for tau_pred_min1.\n");
        exit(1);
    }
    if ((tau_pred_min2 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for tau_pred_min2.\n");
        exit(1);
    }
    if ((tau_pred_max1 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for tau_pred_max1.\n");
        exit(1);
    }
    if ((tau_pred_max2 = (float *) malloc(nwave * sizeof (float))) == NULL) {
        printf("Unable to allocate space for tau_pred_max2.\n");
        exit(1);
    }

    if ((chi1 = (float *) calloc(nmodels, sizeof (float))) == NULL) {
                printf("Unable to allocate space for chi1.\n");
                exit(1);
        }
     if ((chi2 = (float *) calloc(nmodels, sizeof (float))) == NULL) {
                printf("Unable to allocate space for chi2.\n");
                exit(1);
        }


    for (iw = 0; iw < nwave; iw++) {
        taua[iw] = -1.0;
        tsol[iw] = -1.0;
        tsen[iw] = -1.0;
        rhoa1[iw] = rhoa[iw];
        rhoa2[iw] = rhoa[iw];
        rhoa [iw] = BAD_FLT;
    }


    // adjust rh for spectral matching
    if (aer_opt == AERRHSM) {
        if (rh >= 95) {
            //printf("Warning rh is greater than 95%%. Reset to 94%% rh=%f\n", rh);
            rh = 94;
        }
    }

    // find RH index and wts
    if (nrh == 1 || rhtab[0] > rh) { // actual RH < smallest model RH or only one model RH
        irh1 = 0;
        irh2 = 0;
        wt = 0.0;
        derv_wt_rh=0.0;
    } else if (rhtab[nrh - 1] < rh) { // actual RH > smallestlargest model RH
        irh1 = nrh - 1;
        irh2 = nrh - 1;
        wt = 0.0;
        derv_wt_rh=0.0;
    } else {
        for (irh = 0; irh < nrh; irh++) {
            if (rhtab[irh] > rh)
                break;
        }
        irh1 = MIN(MAX(0, irh - 1), nrh - 2);
        irh2 = irh1 + 1;
        wt = (rh - rhtab[irh1]) / (rhtab[irh2] - rhtab[irh1]);
        if(uncertainty)
            derv_wt_rh = 1.0 / (rhtab[irh2] - rhtab[irh1]);     
    }

    // set indices of active model sets

    for (im = 0; im < nsd; im++) {
        mindx1[im] = irh1 * nsd + im;
        mindx2[im] = irh2 * nsd + im;
    }

    // compute aerosol reflectances, aot, diffuse trans, eps from first model set

    /* perform spectral matching from Red to SWIR in radiance space, based on Ahmad
     * Multi-scattering coeff method --> Amir*/
    if (aer_opt == AERRHSM) {

        if (smaer(sensorID, wave, nwave, iwnir_s, iwnir_l, nmodels, mindx1,geom, wv, rhoa1, modmin1, modmax1, modrat1, tau_pred_min1, tau_pred_max1,ip, chi1, mbac_w,uncertainty1) != 0) {

            free(taua1);
            free(taua2);
            free(tsol1);
            free(tsol2);
            free(tsen1);
            free(tsen2);
            free(rhoa1);
            free(rhoa2);
            free(tau_pred_min1);
            free(tau_pred_max1);
            free(tau_pred_min2);
            free(tau_pred_max2);
            return (1);
        }
    } else if (aer_opt == AERRHMSEPS || aer_opt == AERRHMSEPS_lin) {
    	   if (ahmadaer(sensorID, wave, nwave, iwnir_s, iwnir_l, nsd, mindx1,
                geom, wv, rhoa1, modmin1, modmax1, modrat1, &eps1, tau_pred_min1, tau_pred_max1,ip,uncertainty1) != 0) {

            free(taua1);
            free(taua2);
            free(tsol1);
            free(tsol2);
            free(tsen1);
            free(tsen2);
            free(rhoa1);
            free(rhoa2);
            free(tau_pred_min1);
            free(tau_pred_max1);
            free(tau_pred_min2);
            free(tau_pred_max2);
            return (1);

        }
    } else {
        if (wangaer(sensorID, wave, nwave, iwnir_s, iwnir_l, nsd, mindx1,
                geom, wv, rhoa1, modmin1, modmax1, modrat1, &eps1, tau_pred_min1, tau_pred_max1) != 0) {
            free(taua1);
            free(taua2);
            free(tsol1);
            free(tsol2);
            free(tsen1);
            free(tsen2);
            free(rhoa1);
            free(rhoa2);
            free(tau_pred_min1);
            free(tau_pred_max1);
            free(tau_pred_min2);
            free(tau_pred_max2);
            return (1);
        }
    }

    diff_tran(sensorID, wave, nwave, iwnir_l, geom, wv, pr, taur,
            *modmin1, *modmax1, *modrat1, rhoa1, taua1, tsol1, tsen1, tau_pred_min1, tau_pred_max1, 1, ip,uncertainty1);

    // compute aerosol reflectances, aot, diffuse trans, eps from second model set (if needed)

    if (irh2 != irh1) {

        if (aer_opt == AERRHSM) {

            if (smaer(sensorID, wave, nwave, iwnir_s, iwnir_l, nmodels, mindx2,
                    geom, wv, rhoa2, modmin2, modmax2, modrat2, tau_pred_min2, tau_pred_max2, ip, chi2, mbac_w,uncertainty2) != 0) {

                free(taua1);
                free(taua2);
                free(tsol1);
                free(tsol2);
                free(tsen1);
                free(tsen2);
                free(rhoa1);
                free(rhoa2);
                free(tau_pred_min1);
                free(tau_pred_max1);
                free(tau_pred_min2);
                free(tau_pred_max2);

                return (1);
            }

        }

        else if (aer_opt == AERRHMSEPS || aer_opt == AERRHMSEPS_lin) {
            if (ahmadaer(sensorID, wave, nwave, iwnir_s, iwnir_l, nsd, mindx2,
                    geom, wv, rhoa2, modmin2, modmax2, modrat2, &eps2, tau_pred_min2, tau_pred_max2,ip,uncertainty2) != 0) {
                free(taua1);
                free(taua2);
                free(tsol1);
                free(tsol2);
                free(tsen1);
                free(tsen2);
                free(rhoa1);
                free(rhoa2);
                free(tau_pred_min1);
                free(tau_pred_max1);
                free(tau_pred_min2);
                free(tau_pred_max2);
                return (1);
            }
        } else {
            if (wangaer(sensorID, wave, nwave, iwnir_s, iwnir_l, nsd, mindx2,
                    geom, wv, rhoa2, modmin2, modmax2, modrat2, &eps2, tau_pred_min2, tau_pred_max2) != 0) {
                free(taua1);
                free(taua2);
                free(tsol1);
                free(tsol2);
                free(tsen1);
                free(tsen2);
                free(rhoa1);
                free(rhoa2);
                free(tau_pred_min1);
                free(tau_pred_max1);
                free(tau_pred_min2);
                free(tau_pred_max2);
                return (1);
            }
        }

        diff_tran(sensorID, wave, nwave, iwnir_l, geom, wv, pr, taur,
                *modmin2, *modmax2, *modrat2, rhoa2, taua2, tsol2, tsen2, tau_pred_min2, tau_pred_max2, 1,ip,uncertainty2);

        for (iw = 0; iw < nwave; iw++) {
            rhoa[iw] = rhoa1[iw]*(1 - wt) + rhoa2[iw] * wt;
            taua[iw] = taua1[iw]*(1 - wt) + taua2[iw] * wt;
            tsol[iw] = tsol1[iw]*(1 - wt) + tsol2[iw] * wt;
            tsen[iw] = tsen1[iw]*(1 - wt) + tsen2[iw] * wt;

            if(uncertainty){
                for(inir=0;inir<=iwnir_l-iwnir_s;inir++){
                    uncertainty->derv_La_rhorc[iw][inir]=(1-wt)*uncertainty1->derv_La_rhorc[iw][inir]+ wt*uncertainty2->derv_La_rhorc[iw][inir];
                    uncertainty->derv_taua_rhorc[iw][inir]=(1-wt)*uncertainty1->derv_taua_rhorc[iw][inir]+ wt*uncertainty2->derv_taua_rhorc[iw][inir];
                    uncertainty->derv_tsen_rhorc[iw][inir]=(1-wt)*uncertainty1->derv_tsen_rhorc[iw][inir]+ wt*uncertainty2->derv_tsen_rhorc[iw][inir];
                    uncertainty->derv_tsol_rhorc[iw][inir]=(1-wt)*uncertainty1->derv_tsol_rhorc[iw][inir]+ wt*uncertainty2->derv_tsol_rhorc[iw][inir];
                }
                uncertainty->derv_La_taua_l[iw]=(1-wt)*uncertainty1->derv_La_taua_l[iw]+ wt*uncertainty2->derv_La_taua_l[iw];
                uncertainty->derv_La_rhow_l[iw]=(1-wt)*uncertainty1->derv_La_rhow_l[iw]+ wt*uncertainty2->derv_La_rhow_l[iw];
                uncertainty->derv_La_rh[iw]=(rhoa2[iw]-rhoa1[iw])*derv_wt_rh;

                //uncertainty->derv_taua_rhoa_l[iw]=(1-wt)*uncertainty1->derv_taua_rhoa_l[iw]+ wt*uncertainty2->derv_taua_rhoa_l[iw];
                uncertainty->derv_taua_taua_l[iw]=(1-wt)*uncertainty1->derv_taua_taua_l[iw]+ wt*uncertainty2->derv_taua_taua_l[iw];
                //rrrec->derv_taua_rhoa_s[iw]=(1-wt)*uncertainty1->derv_taua_rhoa_s[iw]+ wt*uncertainty2->derv_taua_rhoa_s[iw];
                uncertainty->derv_taua_rhow_l[iw]=(1-wt)*uncertainty1->derv_taua_rhow_l[iw]+ wt*uncertainty2->derv_taua_rhow_l[iw];
                uncertainty->derv_taua_rh[iw]=(taua2[iw]-taua1[iw])*derv_wt_rh;
                //uncertainty->derv_taua_taua_s[iw]=(1-wt)*uncertainty1->derv_taua_taua_s[iw]+ wt*uncertainty2->derv_taua_taua_s[iw];

                // uncertainty->derv_tsol_rhoa_l[iw]=(1-wt)*uncertainty1->derv_tsol_rhoa_l[iw]+ wt*uncertainty2->derv_tsol_rhoa_l[iw];
                //uncertainty->derv_tsol_rhoa_s[iw]=(1-wt)*uncertainty1->derv_tsol_rhoa_s[iw]+ wt*uncertainty2->derv_tsol_rhoa_s[iw];
                uncertainty->derv_tsol_taua_l[iw]=(1-wt)*uncertainty1->derv_tsol_taua_l[iw]+ wt*uncertainty2->derv_tsol_taua_l[iw];
                uncertainty->derv_tsol_rhow_l[iw]=(1-wt)*uncertainty1->derv_tsol_rhow_l[iw]+ wt*uncertainty2->derv_tsol_rhow_l[iw];
                uncertainty->derv_tsol_rh[iw]=(tsol2[iw]-tsol1[iw])*derv_wt_rh;
                //uncertainty->derv_tsol_taua_s[iw]=(1-wt)*uncertainty1->derv_tsol_taua_s[iw]+ wt*uncertainty2->derv_tsol_taua_s[iw];

                //uncertainty->derv_tsen_rhoa_l[iw]=(1-wt)*uncertainty1->derv_tsen_rhoa_l[iw]+ wt*uncertainty2->derv_tsen_rhoa_l[iw];
                //uncertainty->derv_tsen_rhoa_s[iw]=(1-wt)*uncertainty1->derv_tsen_rhoa_s[iw]+ wt*uncertainty2->derv_tsen_rhoa_s[iw];
                uncertainty->derv_tsen_taua_l[iw]=(1-wt)*uncertainty1->derv_tsen_taua_l[iw]+ wt*uncertainty2->derv_tsen_taua_l[iw];
                uncertainty->derv_tsen_rhow_l[iw]=(1-wt)*uncertainty1->derv_tsen_rhow_l[iw]+ wt*uncertainty2->derv_tsen_rhow_l[iw];
                uncertainty->derv_tsen_rh[iw]=(tsen2[iw]-tsen1[iw])*derv_wt_rh;
                //uncertainty->derv_tsen_taua_s[iw]=(1-wt)*uncertainty1->derv_tsen_taua_s[iw]+ wt*uncertainty2->derv_tsen_taua_s[iw];
            }
        }
        *eps = eps1 * (1 - wt) + eps2*wt;
        *chi =chi1[0]*(1-wt) +chi2[0]*wt;

    } else {

        for (iw = 0; iw < nwave; iw++) {
            rhoa[iw] = rhoa1[iw];
            taua[iw] = taua1[iw];
            tsol[iw] = tsol1[iw];
            tsen[iw] = tsen1[iw];
            *chi=chi1[0];

            if(uncertainty){
                for(inir=0;inir<=iwnir_l-iwnir_s;inir++){
                    uncertainty->derv_La_rhorc[iw][inir]=uncertainty1->derv_La_rhorc[iw][inir];
                    uncertainty->derv_taua_rhorc[iw][inir]=uncertainty1->derv_taua_rhorc[iw][inir];
                    uncertainty->derv_tsen_rhorc[iw][inir]=uncertainty1->derv_tsen_rhorc[iw][inir];
                    uncertainty->derv_tsol_rhorc[iw][inir]=uncertainty1->derv_tsol_rhorc[iw][inir];
                }
                uncertainty->derv_La_taua_l[iw]=uncertainty1->derv_La_taua_l[iw];
                //uncertainty->derv_rhoa_l[iw]=uncertainty1->derv_rhoa_l[iw];
                uncertainty->derv_La_rhow_l[iw]=uncertainty1->derv_La_rhow_l[iw];
                //uncertainty->derv_rhoa_s[iw]=uncertainty1->derv_rhoa_s[iw];
                //uncertainty->derv_taua_s[iw]=uncertainty1->derv_taua_s[iw];

                // uncertainty->derv_taua_rhoa_l[iw]=uncertainty1->derv_taua_rhoa_l[iw];
                uncertainty->derv_taua_taua_l[iw]=uncertainty1->derv_taua_taua_l[iw];
                // uncertainty->derv_taua_rhoa_s[iw]=uncertainty1->derv_taua_rhoa_s[iw];
                uncertainty->derv_taua_rhow_l[iw]=uncertainty1->derv_taua_rhow_l[iw];
                //uncertainty->derv_taua_taua_s[iw]=uncertainty1->derv_taua_taua_s[iw];

                // uncertainty->derv_tsol_rhoa_l[iw]=uncertainty1->derv_tsol_rhoa_l[iw];
                // uncertainty->derv_tsol_rhoa_s[iw]=uncertainty1->derv_tsol_rhoa_s[iw];
                uncertainty->derv_tsol_taua_l[iw]=uncertainty1->derv_tsol_taua_l[iw];
                uncertainty->derv_tsol_rhow_l[iw]=uncertainty1->derv_tsol_rhow_l[iw];
                //uncertainty->derv_tsol_taua_s[iw]=uncertainty1->derv_tsol_taua_s[iw];

                //uncertainty->derv_tsen_rhoa_l[iw]=uncertainty1->derv_tsen_rhoa_l[iw];
                //uncertainty->derv_tsen_rhoa_s[iw]=uncertainty1->derv_tsen_rhoa_s[iw];
                uncertainty->derv_tsen_taua_l[iw]=uncertainty1->derv_tsen_taua_l[iw];
                uncertainty->derv_tsen_rhow_l[iw]=uncertainty1->derv_tsen_rhow_l[iw];
                //uncertainty->derv_tsen_taua_s[iw]=uncertainty1->derv_tsen_taua_s[iw];
            }
        }
        *eps = eps1;
    }

    free(taua1);
    free(taua2);
    free(tsol1);
    free(tsol2);
    free(tsen1);
    free(tsen2);
    free(rhoa1);
    free(rhoa2);
    free(tau_pred_min1);
    free(tau_pred_max1);
    free(tau_pred_min2);
    free(tau_pred_max2);
    free(chi1);
    free(chi2);

    return (0);
}


/* ---------------------------------------------------------------------------------------- */
/* fixedaer() - compute aerosol reflectance for fixed aerosol model                         */
/*                                                                                          */
/* B. Franz, August 2004.                                                                   */

/* ---------------------------------------------------------------------------------------- */
int fixedaer(int32_t sensorID, int32_t modnum, float wave[], int32_t nwave, int32_t iwnir_s, int32_t iwnir_l,
        char models[MAXAERMOD][32], int32_t nmodels,
        geom_str *geom, float wv, float rhoa[], float *epsnir) {
    float *eps;
    float *rhoas;
    int iw;

    if ((rhoas = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoas.\n");
        exit(1);
    }

    if (rhoa[iwnir_l] < 0.0) {
        *epsnir = BAD_FLT;
        //        for (iw=0; iw<nwave; iw++)
        //            rhoas[iw] = BAD_FLT;
        free(rhoas);
        return (1);
    }

    /* get model epsilon for all wavelengths at this geometry */
    eps = model_epsilon(modnum, iwnir_l, wave, nwave, geom);

    /* get SS aerosol reflectance at longest wavelength */
    if (rhoa_to_rhoas(sensorID, modnum, geom, wv, rhoa, wave, nwave, iwnir_l, iwnir_l, rhoas) != 0) {
        printf("Error getting rhoas\n");
        free(rhoas);
        return (1);
    }

    /* compute SS aerosol reflectance in visible bands */
    for (iw = 0; iw < nwave; iw++) {
        rhoas[iw] = rhoas[iwnir_l] * eps[iw];
    }

    /* compute MS aerosol reflectance in visible bands */
    rhoas_to_rhoa(sensorID, modnum, geom, wv, rhoas, wave, nwave, 0, nwave - 1, rhoa);

    if (iwnir_s == iwnir_l)
        *epsnir = eps[iwnir_l - 1];
    else
        *epsnir = eps[iwnir_s];

    free(rhoas);

    return (0);
}


/* ---------------------------------------------------------------------------------------- */
/* fixedmodpair() - compute aerosol reflectance for fixed model pair                        */
/*                                                                                          */
/* B. Franz, August 2004.                                                                   */

/* ---------------------------------------------------------------------------------------- */
int fixedmodpair(int32_t sensorID, float wave[], int32_t nwave,
        int32_t iwnir_s, int32_t iwnir_l, geom_str *geom, float wv,
        int32_t modmin, int32_t modmax, float modrat, float rhoa[], float *eps) {
    float *rhoa1;
    float *rhoa2;
    float eps1;
    float eps2;
    int iw;
    int status;

    if ((rhoa1 = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoa1.\n");
        exit(1);
    }
    if ((rhoa2 = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoa2.\n");
        exit(1);
    }

    if (modmin < 0 || modmin >= input->naermodels ||
            modmax < 0 || modmax >= input->naermodels ||
            modrat < 0.0 || modrat > 1.0) {
        printf("Bogus input for fixed model pair: %d %d %f\n", modmin + 1, modmax + 1, modrat);
        exit(1);
    }


    if (rhoa[iwnir_l] > input->rhoamin) {

        rhoa2[iwnir_l] = rhoa1[iwnir_l] = rhoa[iwnir_l];

        /* get aerosol reflectance in visible for fixed model, and epsilon */
        status = fixedaer(sensorID, modmin, wave, nwave, iwnir_s, iwnir_l,
                input->aermodels, input->naermodels, geom, wv, rhoa1, &eps1);
        status = fixedaer(sensorID, modmax, wave, nwave, iwnir_s, iwnir_l,
                input->aermodels, input->naermodels, geom, wv, rhoa2, &eps2);

        /* convert aerosol relectance to radiance */
        if (status == 0) {
            for (iw = 0; iw < nwave; iw++) {
                if (iw != iwnir_l) // without this check tLw-La may go slight negative
                    rhoa[iw] = MAX((1.0 - modrat) * rhoa1[iw] + modrat * rhoa2[iw], 0.0);
            }
            *eps = (1.0 - modrat) * eps1 + modrat*eps2;
        }

    } else if (rhoa[iwnir_l] > -(input->rhoamin)) {

        /* if input NIR is near zero, assume a small white aerosol signal */
        *eps = 1.0;
        for (iw = 0; iw < nwave; iw++) {
            rhoa[iw] = MAX(rhoa[iwnir_l], 1e-6);
        }

        status = 0;

    } else {

        /* if input NIR is very negative, fail the case */
        status = 1;
    }

    free(rhoa1);
    free(rhoa2);

    return (status);
}


/* ---------------------------------------------------------------------------------------- */
/* fixedaot() - compute aerosol reflectance for fixed aot(lambda)                           */
/*                                                                                          */
/* B. Franz, August 2004.                                                                   */

/* ---------------------------------------------------------------------------------------- */

int fixedaot(int32_t sensorID, float aot[], float wave[], int32_t nwave,
        int32_t iwnir_s, int32_t iwnir_l, geom_str *geom, float wv,
        int32_t *modmin, int32_t *modmax, float *modrat, float rhoa[],
        float *epsnir) {
    static int firstCall = 1;
    static int angst_band1 = -1;
    static int angst_band2 = -1;

    float *phase1;
    float *phase2;
    float *f1;
    float *f2;
    float *lnf1;
    float *lnf2;
    float *rhoas1;
    float *rhoas2;
    float *rhoa1;
    float *rhoa2;
    float eps1;
    float eps2;
    float angstrom;
    int ig, gmult, iw, iwtab;
    float maxwave;

    maxwave = MAX(aertab->nwave, nwave);

    if ((rhoa1 = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoa1.\n");
        exit(1);
    }
    if ((rhoa2 = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoa2.\n");
        exit(1);
    }
    if ((rhoas1 = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoas1.\n");
        exit(1);
    }
    if ((rhoas2 = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoas2.\n");
        exit(1);
    }
    if ((f1 = (float *) calloc(maxwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for f1.\n");
        exit(1);
    }
    if ((f2 = (float *) calloc(maxwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for f2.\n");
        exit(1);
    }
    if ((lnf1 = (float *) calloc(maxwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for lnf1.\n");
        exit(1);
    }
    if ((lnf2 = (float *) calloc(maxwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for lnf2.\n");
        exit(1);
    }

    if (firstCall) {
        angst_band1 = windex(520, wave, nwave);
        angst_band2 = windex(865, wave, nwave);
        firstCall = 0;
    }

    /* bail on negative input aot */
    for (iw = 0; iw < nwave; iw++) {
        if (aot[iw] < 0.0) {
            free(rhoa1);
            free(rhoa2);
            free(rhoas1);
            free(rhoas2);
            free(f1);
            free(f2);
            free(lnf1);
            free(lnf2);
            return (1);
        }
    }

    /* compute angstrom and use to select bounding models */
    if (aot[iwnir_l] > 0.0)
        angstrom = -log(aot[angst_band1] / aot[angst_band2]) /
        log(wave[angst_band1] / wave[angst_band2]);
    else
        angstrom = 0.0;

    model_select_angstrom(angstrom, modmin, modmax, modrat);


    /* get model phase function for all wavelengths at this geometry for the two models */
    phase1 = model_phase(*modmin, geom);
    phase2 = model_phase(*modmax, geom);

    gmult = (interpol == 1) ? 0 : geom->gmult;

    /* compute factor for SS approximation, set-up for interpolation */
    for (iw = 0; iw < aertab->nwave; iw++) {
        ig = gmult * iw;
        f1[iw] = aertab->model[*modmin]->albedo[iw] *
                phase1[iw] / 4.0 / geom->csolz[ig] / geom->csenz[ig];
        f2[iw] = aertab->model[*modmax]->albedo[iw] *
                phase2[iw] / 4.0 / geom->csolz[ig] / geom->csenz[ig];
        if (interpol) {
            lnf1[iw] = log(f1[iw]);
            lnf2[iw] = log(f2[iw]);
        }
    }

    /* compute SS aerosol reflectance */
    if (interpol) {
        for (iw = 0; iw < nwave; iw++) {
            iwtab = iwatab[iw];
            if (aertab->wave[iwtab] != wave[iw] && wave[iw] > 0) {
                rhoas1[iw] = aot[iw] * exp(linterp(aertab->wave, lnf1, aertab->nwave, wave[iw]));
                rhoas2[iw] = aot[iw] * exp(linterp(aertab->wave, lnf2, aertab->nwave, wave[iw]));
            } else {
                rhoas1[iw] = aot[iw] * f1[iwtab];
                rhoas2[iw] = aot[iw] * f2[iwtab];
            }
        }
    } else {
        for (iw = 0; iw < nwave; iw++) {
            iwtab = iwatab[iw];
            rhoas1[iw] = aot[iw] * f1[iwtab];
            rhoas2[iw] = aot[iw] * f2[iwtab];
        }
    }
    eps1 = rhoas1[iwnir_s] / rhoas1[iwnir_l];
    eps2 = rhoas2[iwnir_s] / rhoas2[iwnir_l];

    /* compute MS aerosol reflectance */
    rhoas_to_rhoa(sensorID, *modmin, geom, wv, rhoas1, wave, nwave, 0, nwave - 1, rhoa1);
    rhoas_to_rhoa(sensorID, *modmax, geom, wv, rhoas2, wave, nwave, 0, nwave - 1, rhoa2);

    for (iw = 0; iw < nwave; iw++) {
        rhoa[iw] = (1.0 - *modrat) * rhoa1[iw] + *modrat * rhoa2[iw];
    }
    *epsnir = (1.0 - *modrat) * eps1 + *modrat * eps2;

    free(rhoa1);
    free(rhoa2);
    free(rhoas1);
    free(rhoas2);
    free(f1);
    free(f2);
    free(lnf1);
    free(lnf2);

    return (0);
}


/* ---------------------------------------------------------------------------------------- */
/* aerosol() - compute aerosol reflectance using specified algorithm                        */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */

/* ---------------------------------------------------------------------------------------- */
int aerosol(l2str *l2rec, int32_t aer_opt_in, aestr *aerec, int32_t ip,
        float wave[], int32_t nwave, int32_t iwnir_s_in, int32_t iwnir_l_in,
        float F0_in[], float La1_in[], float La2_out[],
        float t_sol_out[], float t_sen_out[], float *eps, float taua_out[],
        int32_t *modmin, int32_t *modmax, float *modrat,
        int32_t *modmin2, int32_t *modmax2, float *modrat2, float *mbac_w) {
    static int firstCall = 1;
    static int32_t mindx[MAXAERMOD];

    int status = 1;
    float *rhoa;
    float *radref;
    float temp;
    float *aot;
    float angstrom;

    int iw, ipb, inir;

    float *F0;
    float *taur;
    float *La1;
    float *La2;
    float *t_sol;
    float *t_sen;
    float *taua;
    float *taua_pred_min;
    float *taua_pred_max;

    static int32_t last_scan=-99;
    static int32_t aer_wave_base;
    float  eps_tmp;
    float chi2=0.;

    l1str *l1rec = l2rec->l1rec;
    uncertainty_t *uncertainty=l1rec->uncertainty;

    static float *noise_temp;

    if(input->aer_opt==AERRHSM && firstCall){
        noise_global=(float *)malloc(l1rec->npix*nwave*sizeof(float));
    }

    if(input->aer_opt==AERRHSM && l1rec->iscan!=last_scan){
        if(uncertainty)
            noise_temp=uncertainty->dsensor;
        else
            noise_temp=get_uncertainty(l1rec);

        last_scan=l1rec->iscan;

    }

    int32_t sensorID = l1rec->l1file->sensorID;
    float solz = l1rec->solz [ip];
    float senz = l1rec->senz [ip];
    float wv = l1rec->wv [ip];
    float rh = l1rec->rh [ip];
    float pr = l1rec->pr [ip];
    int taua_opt;

    if (firstCall) {
        Nbands = nwave;
        Maxband = nwave + 1; /* Must be >= Nbands */
    }

    if ((rhoa = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for rhoa.\n");
        exit(1);
    }
    if ((radref = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for raderef.\n");
        exit(1);
    }
    if ((F0 = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for F0.\n");
        exit(1);
    }
    if ((taur = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for taur.\n");
        exit(1);
    }
    if ((La1 = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for La1.\n");
        exit(1);
    }
    if ((La2 = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for La2.\n");
        exit(1);
    }
    if ((t_sol = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for t_sol.\n");
        exit(1);
    }
    if ((t_sen = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for t_sen.\n");
        exit(1);
    }
    if ((taua = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for taua.\n");
        exit(1);
    }
    if ((taua_pred_min = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for taua_pred_min.\n");
        exit(1);
    }
    if ((taua_pred_max = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("Unable to allocate space for taua_pred_max.\n");
        exit(1);
    }

    /* set up the geometry structure to be passed  */
    if (l1rec->geom_per_band == NULL) {
        /* nominal geometry setup */
        geom.senz = &l1rec->senz[ip];
        geom.solz = &l1rec->solz[ip];
        geom.phi = &l1rec->delphi[ip];
        geom.csolz = &l1rec->csolz[ip];
        geom.csenz = &l1rec->csenz[ip];
        geom.gmult = 0;

        geom.airmass = (float *) malloc(sizeof (float));
        *geom.airmass = 1. / geom.csolz[0] + 1. / geom.csenz[0];
        if ((evalmask & TRANSSPHER) != 0) {
            geom.airmass_sph = (float *) malloc(sizeof (float));
            geom.airmass_plp = (float *) malloc(sizeof (float));
            *geom.airmass_sph = ky_airmass(geom.solz[0]) +
                    ky_airmass(geom.senz[0]);
            *geom.airmass_plp = pp_airmass(geom.solz[0]) +
                    pp_airmass(geom.senz[0]);
        }
    } else {
        ipb = ip * Nbands;
        geom.senz = &l1rec->geom_per_band->senz[ipb];
        geom.solz = &l1rec->geom_per_band->solz[ipb];
        geom.phi = &l1rec->geom_per_band->delphi[ipb];
        geom.csolz = &l1rec->geom_per_band->csolz[ipb];
        geom.csenz = &l1rec->geom_per_band->csenz[ipb];
        geom.gmult = 1;

        geom.airmass = (float *) malloc(Nbands * sizeof (float));
        for (iw = 0; iw < Nbands; iw++) {
            geom.airmass[iw] = 1. / geom.csolz[iw] + 1. / geom.csenz[iw];
        }
        if ((evalmask & TRANSSPHER) != 0) {
            geom.airmass_plp = (float *) malloc(Nbands * sizeof (float));
            geom.airmass_sph = (float *) malloc(Nbands * sizeof (float));
            for (iw = 0; iw < Nbands; iw++) {
                geom.airmass_plp[iw] = pp_airmass(geom.solz[iw]) +
                        pp_airmass(geom.senz[iw]);
                geom.airmass_sph[iw] = ky_airmass(geom.solz[iw]) +
                        ky_airmass(geom.senz[iw]);
            }
        }
    }

    /* set static global evaluation level */
    evalmask = l1_input->evalmask;
    aer_opt = aer_opt_in;

    /* transfer inputs per band to inputs per wavelength */
    for (iw = 0; iw < nwave; iw++) {
        F0 [iw] = F0_in [iw];
        taur[iw] = l1rec->l1file->Tau_r[iw];
        La1 [iw] = La1_in[iw];
        if (iw == iwnir_s_in) iwnir_s = iw;
        if (iw == iwnir_l_in) iwnir_l = iw;
    }

    /* compute total airmass (static global) */
    mu0 = cos(solz / radeg);
    mu = cos(senz / radeg);
    airmass = 1.0 / mu0 + 1.0 / mu;
    if ((evalmask & TRANSSPHER) != 0) {
        airmass_plp = pp_airmass(solz) + pp_airmass(senz);
        airmass_sph = ky_airmass(solz) + ky_airmass(senz);
    }

    /* initialize epsilon and diffuse transmittance to defaults */
    *eps = 1.0;
    *modmin = BAD_INT;
    *modmax = BAD_INT;
    *modrat = BAD_FLT;
    *modmin2 = BAD_INT;
    *modmax2 = BAD_INT;
    *modrat2 = BAD_FLT;
    for (iw = 0; iw < nwave; iw++) {
        taua [iw] = 0.0;
        t_sol [iw] = 1.0;
        t_sen [iw] = 1.0;
        La2 [iw] = 0.0;
        radref[iw] = pi / F0[iw] / geom.csolz[iw * geom.gmult];

        if(input->aer_opt==AERRHSM){
            if(noise_temp)
                noise_global[ip*nwave+iw]=noise_temp[ip*nwave+iw]*radref[iw];
            else
                noise_global[ip*nwave+iw]=1.0;
        }
    }

    /* load aerosol model tables */
    if (!loaded) {
        int32_t im;
        load_aermod(sensorID, wave, nwave, input->aermodfile, input->aermodels, input->naermodels);
        for (im = 0; im < aertab->nmodel; im++) mindx[im] = im;
        if (have_rh && aertab->nmodel >= 30) {
            printf("\nLimiting aerosol models based on RH.\n");
            use_rh = 1;
        }
    }
    /* do not permit use of geom_per_band if interpolation is done */
    /* WDR note that it CAN work, but currently, we don't want users to
       use this combination */
    if ((interpol == 1) && (geom.gmult == 1)) {
        fprintf(stderr, "-E- %s line %d: Interpolated aerosol tables are\n",
                __FILE__, __LINE__);
        fprintf(stderr, "           not permitted for use with band-dependent geometry, set geom_per_band=0\n");
        exit(1);
    }

    /* Change the aerosol option if need be */
    if (use_rh)
        switch (aer_opt) {
        case AERWANG: aer_opt = AERRH;
            break;
        case AERWANGNIR: aer_opt = AERRHNIR;
            break;
        case AERWANGSWIR: aer_opt = AERRHSWIR;
            break;
        case AERMUMM: aer_opt = AERRHMUMM;
            break;
        }


    if (firstCall) {
        aer_wave_base=bindex_get(input->aer_wave_base);
        if (aer_opt > 0 && aer_opt <= MAXAERMOD) {
            printf("\nUsing fixed aerosol model #%d (%s)\n", aer_opt, input->aermodels[aer_opt - 1]);
            printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
        } else {
            switch (aer_opt) {
            case AERRHSM: // Spectral Matching --> Amir
                printf("\nUsing Spectral Matching of aerosols reflectance for\n");
                printf("wavelength from %4.1f nm to %4.1f nm for model selection\n", wave[iwnir_s], wave[iwnir_l]);
                break;
            case AERWHITE:
                printf("\nUsing white-aerosol approximation\n");
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            case AERWANG:
            case AERRH:
                printf("\nUsing Gordon & Wang aerosol model selection\n");
                printf("Using bands at %4.1f and %4.1f nm for model selection\n", wave[iwnir_s], wave[iwnir_l]);
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            case AERWANGNIR:
            case AERRHNIR:
                printf("\nUsing Gordon & Wang aerosol model selection\n");
                printf("  and NIR correction with up to %d iterations\n", input->aer_iter_max);
                printf("Using bands at %4.1f and %4.1f nm for model selection\n", wave[iwnir_s], wave[iwnir_l]);
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            case AERWANGSWIR:
            case AERRHSWIR:
                printf("\nUsing Gordon & Wang aerosol model selection with NIR/SWIR switching.\n");
                printf("NIR correction with up to %d iterations\n", input->aer_iter_max);
                printf("NIR bands at %d and %d nm\n", input->aer_wave_short, input->aer_wave_long);
                printf("SWIR bands at %d and %d nm\n\n", input->aer_swir_short, input->aer_swir_long);
                break;
            case AERMUMM:
            case AERRHMUMM:
                printf("\nUsing Gordon & Wang aerosol model selection\n");
                printf("  and MUMM correction\n");
                printf("Using bands at %4.1f and %4.1f nm for model selection\n", wave[iwnir_s], wave[iwnir_l]);
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            case AERRHMSEPS:
                printf("\nUsing multi-scattering aerosol model selection.\n");
                printf("  and NIR correction with up to %d iterations\n", input->aer_iter_max);
                printf("Using bands at %4.1f and %4.1f nm for model selection\n", wave[iwnir_s], wave[iwnir_l]);
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            case AERRHMSEPS_lin:
                printf("\nUsing multi-scattering aerosol model selection in linear space.\n");
                printf("  and NIR correction with up to %d iterations\n", input->aer_iter_max);
                printf("Using bands at %4.1f and %4.1f nm for model selection\n", wave[iwnir_s], wave[iwnir_l]);
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            case FIXAOT:
                printf("\nUsing fixed, input aerosol optical thicknesses for aerosol selection.\n");
                break;
            case FIXMODPAIR:
                printf("\nUsing multiple scattering aerosols with fixed model pair\n");
                break;
            case FIXMODPAIRNIR:
                printf("\nUsing multiple scattering aerosols with fixed model pair\n");
                printf("  and NIR iteration with up to %d iterations\n", input->aer_iter_max);
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            case FIXANGSTROM:
                printf("\nUsing fixed aerosol model based on predefined Angstrom exponent\n");
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            case FIXANGSTROMNIR:
                printf("\nUsing fixed aerosol model based on predefined Angstrom exponent\n");
                printf("  and NIR iteration with up to %d iterations\n", input->aer_iter_max);
                printf("Extrapolating from %4.1f nm\n", wave[iwnir_l]);
                break;
            default:
                printf("\nErroneous aerosol option, %d\n", aer_opt);
                exit(FATAL_ERROR);
                break;
            }
        }
        firstCall = 0;
    }


    /* convert input aerosol radiances to relectance */
    for (iw = iwnir_s; iw <= iwnir_l; iw += MAX(iwnir_l - iwnir_s, 1))
        rhoa[iw] = La1[iw] * radref[iw];


    /* compute aerosol using requested method */
    /* -------------------------------------- */

    switch (aer_opt) {

    case AERWANG: case AERWANGNIR: case AERWANGSWIR: case AERMUMM:

        if (iwnir_l <= iwnir_s || wave[iwnir_s] < 600) {
            printf("Aerosol selection bands must be greater than 600nm with short wave less than long wave (%d,%d)\n", iwnir_l, iwnir_s);
            exit(1);
        }

        /* Multi-Scattering with Gordon & Wang Aerosol Selection */
        /* ----------------------------------------------------- */

        /* convert input NIR aerosol radiances to reflectance */
        for (iw = iwnir_s; iw <= iwnir_l; iw++) {
            rhoa[iw] = La1[iw] * radref[iw];
        }

        /* require sufficient signal in two NIR bands */
        if (rhoa[iwnir_s] > input->rhoamin && rhoa[iwnir_l] > input->rhoamin) {

            /* require MS epsilon to be reasonable */
            if (La1[iwnir_s] / La1[iwnir_l] > 0.1) {

                status = wangaer(sensorID, wave, nwave, iwnir_s, iwnir_l,
                        aertab->nmodel, mindx,
                        &geom, wv, rhoa, modmin, modmax, modrat, eps, taua, taua); // this taua is not used //

                if (status == 0)
                    for (iw = 0; iw < nwave; iw++)
                        La2[iw] = rhoa[iw] / radref[iw];
            }

        } else if (rhoa[iwnir_s] > -(input->rhoamin) && rhoa[iwnir_l] > -(input->rhoamin)) {

            /* if input NIR is near zero, assume a small white aerosol signal */
            *eps = 1.0;
            *modmin = aertab->nmodel;
            *modmax = aertab->nmodel;
            *modrat = 0.0;
            temp = MAX(rhoa[iwnir_l], 1e-6);
            for (iw = 0; iw < nwave; iw++) {
                rhoa[iw] = temp;
                La2 [iw] = rhoa[iw] / radref[iw];
            }

            status = 0;

        } else {

            /* if input NIR is very negative, fail the case */
            status = 1;
        }

        break;

    case AERRH:
    case AERRHNIR:
    case AERRHSWIR:
    case AERRHMUMM:
    case AERRHMSEPS:
    case AERRHMSEPS_lin:
    case AERRHSM: // spectral matching --> Amir

        if (iwnir_l <= iwnir_s || wave[iwnir_s] < 600) {
            printf("Aerosol selection bands must be greater than 600nm with short wave less than long wave");
            exit(1);
        }

        /* Multi-Scattering with RH-based selection              */
        /* ----------------------------------------------------- */

        /* convert input NIR aerosol radiances to relectance */
        for (iw = iwnir_s; iw <= iwnir_l; iw++) {
            rhoa[iw] = La1[iw] * radref[iw];
        }

        /* require sufficient signal in two NIR bands */
        if (rhoa[iwnir_s] > input->rhoamin && rhoa[iwnir_l] > input->rhoamin) {

            /* require MS epsilon to be reasonable */
            eps_tmp=rhoa[iwnir_s] / rhoa[iwnir_l];
            if(aer_opt==AERRHSM)
                eps_tmp=rhoa[iwnir_s] / rhoa[aer_wave_base];
            if (eps_tmp > 0.1 && eps_tmp < 10.0) {

                status = rhaer(sensorID, wave, nwave, iwnir_s, iwnir_l,
                        &geom, wv, rh, pr, taur, rhoa,
                        modmin, modmax, modrat, modmin2, modmax2, modrat2, eps, taua, t_sol, t_sen, ip,mbac_w,&chi2,uncertainty);
                status = 0;
                if (status == 0){
                    l2rec->chi2[ip]=chi2;
                    for (iw = 0; iw < nwave; iw++){
                        La2[iw] = rhoa[iw] / radref[iw];

                    	if(uncertainty){
                    	    uncertainty->derv_La_taua_l[iw] /= radref[iw];
                    	    uncertainty->derv_La_rhow_l[iw] /= radref[iw];
                    	    uncertainty->derv_La_rh[iw] /= radref[iw];
                    	    for(inir=0;inir<=iwnir_l-iwnir_s;inir++){
                    	        uncertainty->derv_La_rhorc[iw][inir] /= radref[iw];
                    	    }
                    	}
                    }
                }
            }

        } else if (rhoa[iwnir_s] > -(input->rhoamin) && rhoa[iwnir_l] > -(input->rhoamin)) {

            /* if input NIR is near zero, assume a small white aerosol signal */
            *eps = 1.0;
            *modmin = aertab->nmodel;
            *modmax = aertab->nmodel;
            *modrat = 0.0;
            *modmin2 = aertab->nmodel;
            *modmax2 = aertab->nmodel;
            *modrat2 = 0.0;
            temp = MAX(rhoa[iwnir_l], 1e-6);
            for (iw = 0; iw < nwave; iw++) {
                rhoa[iw] = temp;
                La2 [iw] = rhoa[iw] / radref[iw];
            }

	    // The diff_tran function needs taua. In this exception case, taua was not yet computed.  The taua_opt
	    // tells diff_tran what method to use for the computation (0=SS, 2=MS).  I used the aer_opt to decide,
	    // but we could use the global have_ms flag.  For now I want to limit the impact, as some sensors are
	    // run with the G&W SS aerosol option, but their tables "have_ms".  BAF, 6/2022.
	    
            if (aer_opt == AERRHMSEPS || aer_opt == AERRHMSEPS_lin|| aer_opt == AERRHSM)
                taua_opt = 2;
            else
                taua_opt = 0;

            diff_tran(sensorID, wave, nwave, iwnir_l, &geom, wv, pr, taur,
                     *modmin, *modmax, *modrat, rhoa, taua, t_sol, t_sen, taua_pred_min, taua_pred_max, taua_opt, ip, uncertainty);


            status = 0;

        } else {

            /* if input NIR is very negative, fail the case */
            status = 1;
        }

        break;

    case AERWHITE:

        /* White Aerosol */
        /* ------------- */

        if (La1[iwnir_l] > 0.0) {

            *eps = 1.0;
            *modmin = 0;
            *modmax = 0;
            *modrat = 0.0;

            for (iw = 0; iw < nwave; iw++) {
                La2 [iw] = *eps * F0[iw] / F0[iwnir_l] * La1[iwnir_l];
                rhoa[iw] = La2[iw] * radref[iw];
            }


            status = 0;
        }
        break;

    case FIXMODPAIR: case FIXMODPAIRNIR:

        /* Multi-Scattering with Fixed Model Pair */
        /* -------------------------------------- */

        if (aerec != NULL && aerec->mode == ON) {
            *modmin = aerec->mod_min[ip] - 1;
            *modmax = aerec->mod_max[ip] - 1;
            *modrat = aerec->mod_rat[ip];
        } else {
            *modmin = input->aermodmin - 1;
            *modmax = input->aermodmax - 1;
            *modrat = input->aermodrat;
        }

        status = fixedmodpair(sensorID, wave, nwave, iwnir_s, iwnir_l, &geom, wv,
                *modmin, *modmax, *modrat, rhoa, eps);
        if (status == 0) {
            for (iw = 0; iw < nwave; iw++)
                La2[iw] = rhoa[iw] / radref[iw];
        }

        break;

    case FIXANGSTROM: case FIXANGSTROMNIR:

        if (input->aer_angstrom > -2) {
            angstrom = input->aer_angstrom;
        } else {
            int16_t year, day;
            double sec;
            unix2yds(l2rec->l1rec->scantime, &year, &day, &sec);
            angstrom = bin_climatology(input->aerbinfile, day, l1rec->lon[ip], l1rec->lat[ip], "angstrom");
        }

        if (angstrom > -2) {

            model_select_angstrom(angstrom, modmin, modmax, modrat);

            status = fixedmodpair(sensorID, wave, nwave, iwnir_s, iwnir_l, &geom, wv,
                    *modmin, *modmax, *modrat, rhoa, eps);
        } else
            status = 1;

        if (status == 0) {
            for (iw = 0; iw < nwave; iw++)
                La2[iw] = rhoa[iw] / radref[iw];
        }

        break;

    case FIXAOT:

        /* Multi-Scattering with fixed AOT */
        /* ------------------------------- */
        if (aerec != NULL && aerec->mode == ON)
            aot = &aerec->taua[ip * Nbands];
        else
            aot = input->taua;

        status = fixedaot(sensorID, aot, wave, nwave, iwnir_s, iwnir_l, &geom, wv,
                modmin, modmax, modrat, rhoa, eps);
        if (status == 0) {
            for (iw = 0; iw < nwave; iw++)
                La2[iw] = rhoa[iw] / radref[iw];
        }

        break;

    default:

        /* Multi-Scattering with Fixed Model */
        /* --------------------------------- */

        *modmin = aer_opt - 1;
        *modmax = aer_opt - 1;
        *modrat = 0.0;

        if (*modmin < 0 || *modmin > input->naermodels - 1) {
            printf("Invalid aerosol option: %d\n", *modmin);
            exit(1);
        }

        /* convert input NIR aerosol radiance to relectance */
        rhoa[iwnir_l] = La1[iwnir_l] * radref[iwnir_l];

        /* get aerosol reflectance in visible for fixed model, and epsilon */
        status = fixedaer(sensorID, *modmin, wave, nwave, iwnir_s, iwnir_l,
                input->aermodels, input->naermodels,
                &geom, wv, rhoa, eps);

        /* convert aerosol relectance to radiance */
        if (status == 0)
            for (iw = 0; iw < nwave; iw++)
                La2[iw] = rhoa[iw] / radref[iw];
        break;
    }

    /* compute diffuse transmittance through aerosol and Rayleigh, if not yet computed */
    if (status == 0 && aer_opt != AERRHNIR && aer_opt != AERRHSWIR && aer_opt != AERRH && aer_opt != AERRHMSEPS && aer_opt != AERRHSM) {

        diff_tran(sensorID, wave, nwave, iwnir_l, &geom, wv, pr, taur,
                *modmin, *modmax, *modrat, rhoa, taua, t_sol, t_sen, taua_pred_min, taua_pred_max, 0,ip, uncertainty);
    }

    /* transfer outputs per wavelength to outputs per band */
    for (iw = 0; iw < nwave; iw++) {
        La2_out [iw] = La2 [iw];
        taua_out [iw] = taua [iw];
        t_sol_out[iw] = t_sol[iw];
        t_sen_out[iw] = t_sen[iw];
    }

    /* model numbers are reported as 1-based */
    *modmin = *modmin + 1;
    *modmax = *modmax + 1;
    *modmin2 = *modmin2 + 1;
    *modmax2 = *modmax2 + 1;

    free(rhoa);
    free(radref);
    free(F0);
    free(taur);
    free(La1);
    free(La2);
    free(t_sol);
    free(t_sen);
    free(taua);
    free(taua_pred_min);
    free(taua_pred_max);
    free(geom.airmass);
    if ((evalmask & TRANSSPHER) != 0) {
        free(geom.airmass_plp);
        free(geom.airmass_sph);
    }

    return (status);
}



/* --------------------------------------------------------------- */
/* get_angstrom.c - compute angstrom coefficient                   */
/*                                                                 */
/* Inputs:                                                         */
/*     l2rec - level-2 structure containing one complete scan      */
/*             after atmospheric correction.                       */
/*     band  = band number 0-7                                     */
/* Outputs:                                                        */
/*     angst - angstrom coefficient                                */
/*                                                                 */
/* Written By: B. A. Franz, SAIC, August 2004                      */
/*                                                                 */

/* --------------------------------------------------------------- */
void get_angstrom(l2str *l2rec, int band, float angst[]) {
    static int firstCall = 1;
    static int32_t ib2;
    static float wave2;

    int32_t ip;
    int32_t ib1;
    float wave1;
    float aot1;
    float aot2;

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;

    if (firstCall) {
        ib2 = windex(865.0, l1file->fwave, l1file->nbands);
        wave2 = l1file->fwave[ib2];
        firstCall = 0;
    }

    if (band < 0)
        ib1 = windex(443.0, l1file->fwave, l1file->nbands);
    else
        ib1 = band;
    wave1 = l1file->fwave[ib1];

    for (ip = 0; ip < l1rec->npix; ip++) {
        aot1 = l2rec->taua[ip * l1file->nbands + ib1];
        aot2 = l2rec->taua[ip * l1file->nbands + ib2];
        if (aot1 > 0.0 && aot2 > 0.0)
            angst[ip] = -log(aot1 / aot2) / log(wave1 / wave2);
        else if (aot1 == 0.0 && aot2 == 0.0)
            angst[ip] = 0.0;
        else {
            angst[ip] = BAD_FLT;
            l1rec->flags[ip] |= PRODFAIL;
        }
        //     printf("angst = %f aot1= %f aot2= %f ib1=%d ib2=%d\n",angst[ip],aot1,aot2,ib1,ib2);
    }

    return;
}


/* --------------------------------------------------------------- */
/* get_ms_epsilon.c -                                              */

/* --------------------------------------------------------------- */
void get_ms_epsilon(l2str *l2rec, float eps[]) {
    int32_t ip, ipb1, ipb2;

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;

    for (ip = 0; ip < l1rec->npix; ip++) {
        ipb1 = ip * l1file->nbands + iwnir_s;
        ipb2 = ip * l1file->nbands + iwnir_l;
        if (l2rec->La[ipb2] > 0.0) {
            eps[ip] = l2rec->La[ipb1] / l2rec->La[ipb2]
                    * l1rec->Fo[iwnir_l] / l1rec->Fo[iwnir_s];
        } else {
            /* "should" indicate ATMFAIL, but just to be safe */
            eps[ip] = BAD_FLT;
            l1rec->flags[ip] |= PRODFAIL;
        }
    }

    return;
}
