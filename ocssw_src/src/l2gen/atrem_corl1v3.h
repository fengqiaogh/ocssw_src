/*
 * atrem.h
 *
 *  Created on: Feb 11, 2015
 *      Author: rhealy
 */

#ifndef SRC_ATREM_ATREM_H_
#define SRC_ATREM_ATREM_H_
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <timeutils.h>
#include "l12_proto.h"

#define NH2OMAX       60
#define NH2OMAXM1  (NH2OMAX-1)
#define FINSTMAX     100
#define MODELMAX      25
#define ABSCF_MAXLYR  19
#define NO3PT       5001
#define NO2PT       5001
#define NINSTRF     5000
#define NP_HI     300000
#define NP_MED     25401
#define NP_STD     28001
#define NINSTR_MAX  3001
#define NPSHIF      2600 //NPSHIF = NP_STD - NP_MED

#define VSTART   0.56   //Starting wavelength for calculations
#define VEND     3.1
#define DWAVLN   1e-4
#define DWAVNO   0.05   //Point spacing of high res. spectra (cm-1).
#define DLT_MED  0.0002 // 0.2-nm medium resolution spectrum.
#define FACDLT   2.0    //Factor to multiply DLT by to determine the
//range used in the Gaussian function
//calculation for smoothing spectra.

#define ATREM_O3    1
#define ATREM_CO2   2
#define ATREM_NO2   4
#define ATREM_CO    8
#define ATREM_CH4  16
#define ATREM_O2   32
#define ATREM_N2O  64

#define RAD_DEG    OEL_RADEG
#define CONST1      2.7725887    // CONST1=4.0*ln(2)=2.7725887

//typedef float t_array[NBANDS];

typedef struct param_table {
    /*Number of narrow channels to form broader window and absorption channels.  */
    int32_t nb1, nb2, nb3, nb4;
    /*number of points used in channel ratios for both the .94- and 1.14-um regions */
    int32_t nbp94, nb1p14;

    int32_t nh2o; /* number of water vapor values */
    int32_t nbands; /* number of spectral bands in sensor */
    /* 3-channel ratioing
     *   parameters for bands up to the 0.94-um water vapor band [index=0,1]
     *   parameters for bands up to the 1.14-um water vapor band [index=2,3]
     *                          */
    int32_t start_ndx[4]; //ist1,ist2,ist3,ist4
    int32_t end_ndx[4]; //ied1,ied2,ied3,ied4
    /*  3-channel ratioing   parameters for the 0.94-um water vapor band */
    int32_t start_p94;
    int32_t end_p94;
    /*  3-channel ratioing parameters for the 1.14-um water vapor band */
    int32_t start_1p14;
    int32_t end_1p14;
    /* Parameters for smoothing output reflectance spectra. */
    int32_t start2;
    int32_t end2;
    int32_t ncv2;
    /* number of channels of the four AVIRIS spectrometers. */
    int32_t natot, nbtot, nctot, ndtot;
    /* how often to recalculate geometry
     * dogeom = 1 - every pixel
     */
    int32_t dogeom;
    /* Atmospheric model number */
    int model;
    /* Relative weights for the four window     *
     *      channels used in channel-ratioing calculations */
    int idx450; // 450 nm wavelength index

    double wt1, wt2, wt3, wt4;
    /*
     *   delta, delta2 - resolution, in units of nm, of input
     *                   spectra and resolution of output surface reflectance
     *                   spectra. If DLT2>DLT, output spectra are smoothed
     *                   using a gaussian function.
     */
    double delta, delta2;
    /*  The "equivalent" geometrical     *
     *       factor corresponding to the total slant vapor amount      *
     *       VAP_SLANT_MDL and the column vapor amount CLMVAP.  */
    double g_vap_equiv; // This depends on zenith angle and lat/lon of measurement
    float *r0p94; //ratio for the 0.94 um H2O absorption band
    float *r1p14; //ratio for the 1.14 um H2O absorption band
    float *finst2; // some kind of smoothing factor calculated in INIT_SPECCAL only used for AVIRIS?
    /* TOTAL SUN-SURFACE-SENSOR PATH WATER VAPOR IN UNITS OF CM */
    float *vaptot;

    /* total transmittances of all gases that match the     *
     *                         resolutions of imaging spectrometers */
    float **trntbl;
    float water_vapor; // returned water vapor value from get_atrem
    int ja, jb; // indices to tran_table from get_atrem
    float f1a, f2a, f1b, f2b; // fractions for interpolation of transmittances in tran_table from get_atrem
} paramstr;

float get_atrem(float *tg_tot, float *rhot, paramstr *P);
int get_atrem_cor(l1str *l1rec, int32_t ip, float *rhot, float *tg_tot, double *tg_sol, double *tg_sen);
int init_atrem(int32_t sensorID, paramstr *P, l1str *l2rec, int32_t nbands);
int32_t rdatreminfo(int32_t sensorID, int32_t evalmask, const char *pname, void **pval);
int get_angle_limits(float **angle_limit, float **senz, float **solz, int *n_senz, int *n_solz);
float get_current_angle_limit(float insenz, float insolz, int *i, int *j, float **anglelimit, float senz[], float solz[], int n_senz, int n_solz);
int32_t hunt(float *xx, int32_t n, double x, int32_t jlo);
int init_tpvmr(int model);
int getModelNum(float lat, int32_t day);
int32_t findMatch(float *list, int32_t nbands, float elem);
void channelRatio();
void ecdf_(float *xcdf, float *ycdf, int32_t *bin_number, float *xs, int32_t *sample_size);
void kdistgasabs(float *kcdf, float *abscf, float*waveno, float *wavobs, int32_t np_hi, int32_t nlayers, int32_t nbands);
void model_adjust();
void locate_pos_(float *xx, int32_t *n1, float *x1, int32_t *jj);
void geometry();
void init_spectral_calculations();
void get_abscf_data(int levels, int bands, int sds_id, char filename[FILENAME_MAX], float* abscf, char *varname);
void get_tpvmr(size_t layers, size_t models, int sds_id,
        char filename[FILENAME_MAX], char *varname, float* var_a);
int init_tpvmr_nc(int model);
void tran_table();
void tran_smooth(float *tran_hi);
void tran_smooth_others();

struct {
    int32_t h2o, co2, o3, n2o, co, ch4, o2, no2;
} getinput1_;

struct {
    char filename[FILENAME_MAX];
    int32_t dln;
} input_l2gen_;

struct {
    float *tg_sol, *tg_sen, *tg_solo, *tg_seno; // nbands
} tran_table_l2gen_;

struct {
    float *tran_hi_sa[2], *tran_hi_sap1[2], *tran_hi_sb[2], *tran_hi_sbp1[2];
} tran_tables_;

struct {
    float *tran_med_index_sa_sol, *tran_med_index_sa_sen,
    *tran_med_index_sap1_sol, *tran_med_index_sap1_sen,
    *tran_med_index_sb_sol, *tran_med_index_sb_sen,
    *tran_med_index_sbp1_sol, *tran_med_index_sbp1_sen;
    float *tran_med_sa_sol, *tran_med_sa_sen,
    *tran_med_sap1_sol, *tran_med_sap1_sen,
    *tran_med_sb_sol, *tran_med_sb_sen,
    *tran_med_sbp1_sol, *tran_med_sbp1_sen;
    float *tran_std_sa_sol, *tran_std_sa_sen,
    *tran_std_sap1_sol, *tran_std_sap1_sen,
    *tran_std_sb_sol, *tran_std_sb_sen,
    *tran_std_sbp1_sol, *tran_std_sbp1_sen;

} tran_tables1_;

struct {
    float h[MODELMAX], t[MODELMAX], p[MODELMAX], vmr[MODELMAX];
    int32_t nb, nl, model, iaer; //iaer not used because call to ssssss routine commented out in fortran code
    float v, taer55, vrto3, sno2;
    //nb      = number of atmospheric levels in model
    //nl      = nb - 1
    //sno2    = scaling factor for no2 (builtin NO2 column amount is 5.0E+15 molecules/cm^2
    //vrto3   = column ozone amount (atm-cm) 0.28-0.55 is typical
    //v       = visibility (km)
    //iaer    = aerosol model value
    //taerr55 = aerosol optical depth at 550 nm
} getinput3_;

struct {
    float *wavobs, *fwhm;
} getinput4_;

struct {
    int32_t nbands, full_calc;
    float hsurf, dlt, dlt2;
} getinput5_;

struct {
    float wndow1, wndow2, wp94c, wndow3, wndow4, w1p14c;
} getinput6_;

struct {
    int32_t nb1, nb2, nbp94, nb3, nb4, nb1p14;
} getinput7_;

struct {
    int32_t imn, idy, iyr, ih, im, is;
} getinput8_;

struct {
    float xpss, xppp;
} getinput14_;

struct {
    float clmvap, q;
} model_adj1_;

struct {
    float hp[MODELMAX], tp[MODELMAX], pp[MODELMAX], vmrp[MODELMAX];
} model_adj2_;

struct {
    int32_t k_plane;
    float dvap_plane, dvap_layer, dp_plane, dp_layer, clmvapp;
} model_adj3_;

struct {
    int32_t k_surf;
} model_adj4_;

struct {
    float *tran_hi_others; //Transmittance of high res. data(.05cm-1)
} init_speccal1_;

struct {
    int32_t nh2o;
} init_speccal3_;

struct {
    float dp[MODELMAX], pm[MODELMAX], tm[MODELMAX], vmrm[MODELMAX];
} init_speccal5_;

struct {
    int32_t ist1, ied1, ist2, ied2, istp94, iedp94;
} init_speccal6_;

struct {
    int32_t ist3, ied3, ist4, ied4, ist1p14, ied1p14;
} init_speccal7_;

struct {
    float wt1, wt2, wt3, wt4;
    int32_t ja;
} init_speccal8_;

struct {
    int32_t ncv2, ncvhf2, ncvtt2;
    int32_t istrt2, iend2;
    float finst2[FINSTMAX];
} init_speccal10_;

struct {
    int32_t natot, nbtot, nctot, ndtot;
} init_speccal11_;

struct {
    float *wavln_med, *wavln_std;
} init_speccal12_;

struct {
    int32_t *index_med;
    float *wavln_med_index, **tran_med_index;
} init_speccal13_;

/*       Arrays for wavelength positions and FWHM of measured imaging spectrometer
         data and for smoothing the medium resolution spectrum (FWHM = 0.2 nm,
         and point spacing = 0.1 nm) to match coarser resolution spectrum
         of imaging spectrometer data
 */
struct {
    float **finstr;
    int32_t *ncvhf;
} init_speccal15_;

struct {
    float *tran_o3_std;
} init_speccal16_;

struct {
    float *tran_no2_std;
} init_speccal17_;

struct {
    float *sh2o, *vaptot, *r0p94, *r1p14, **trntbl;
    float **tran_kd, **diff_tran, *trntblo;
} tran_table1_;

struct {
    float solzni, solaz, obszni, obsphi;
    int32_t day;
} geometry1_;

struct {
    float gco2, go3, gn2o, gco, gch4, go2, ssh2o[NH2OMAX], totlo3, ggeom;
} geometry2_;

struct {
    float g_vap[MODELMAX], g_other[MODELMAX], g_vap_equiv;
} geometry3_;

struct {
    float vap_slant_mdl;
} geometry4_;

struct {
    float mu, mu0, ssh2o_s[2][NH2OMAX];
} geometry5_;

struct {
    float senzn_l2, senaz_l2, solzn_l2;
    float water_vapor;
    int32_t ja, jb;
    int32_t splitpaths;
    float f1a, f2a, f1b, f2b;
} geometry_l2gen_;

struct {
    float tpvmr[81][7];
} tpvmr_init1_;
//struct {
//    float o3cf[NO3PT];
//} o3cf_init1_;
//struct {
//    float rno2cf[NO2PT];
//} no2cf_init1_;

struct {
    float cst1, cst2, cst3, cst4, cst5, cst6;
    double rp94, r1p14;
    int32_t jac, jbc;
} debug_atrem;
#endif /* SRC_ATREM_ATREM_H_ */
