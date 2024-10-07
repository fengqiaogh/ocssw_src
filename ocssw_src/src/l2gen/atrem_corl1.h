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

#define NBANDS  1024 /* maximum number of bands */
#define NH2OMAX    60
#define NH2OMAXM1  (NH2OMAX-1)
#define FINSTMAX  100
#define MODELMAX   25

#define ATREM_O3    1
#define ATREM_CO2   2
#define ATREM_NO2   4
#define ATREM_CO    8
#define ATREM_CH4  16
#define ATREM_O2   32
#define ATREM_N2O  64

typedef float t_array[NBANDS];

typedef struct param_table {
    /*Number of narrow channels to form broader window and absorption channels.  */
    int32_t nb1, nb2, nb3, nb4;
    /*number of points used in channel ratios for both the .94- and 1.14-um regions */
    int32_t nbp94, nb1p14;

    int32_t nh2o; /* number of water vapor values */
    int32_t nobs; /* number of spectral observations - this should be nbands? number of channels? */
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
    t_array *trntbl;
    float water_vapor; // returned water vapor value from get_atrem
    int ja, jb; // indices to tran_table from get_atrem
    float f1a, f2a, f1b, f2b; // fractions for interpolation of transmittances in tran_table from get_atrem
} paramstr;

float get_atrem(float *tg_tot, float *rhot, paramstr *P);
int get_atrem_cor(l1str *l1rec, int32_t ip, float *rhot, float *tg_tot, float *tg_sol, float *tg_sen);
int init_atrem(int32_t sensorID, paramstr *P, l1str *l2rec, int32_t nbands);
int32_t rdatreminfo(int32_t sensorID, int32_t evalmask, const char *pname, void **pval);
int get_angle_limits(float **angle_limit, float **senz, float **solz, int *n_senz, int *n_solz);
float get_current_angle_limit(float insenz, float insolz, int *i, int *j, float **anglelimit, float senz[], float solz[], int n_senz, int n_solz);
int32_t hunt(float *xx, int32_t n, double x, int32_t jlo);
int init_tpvmr(int model);
int getModelNum(float lat, int32_t day);
int findMatch_(float *list, int *nobs, float *elem);
void channelRatio_();
void ecdf_(float *xcdf, float *ycdf, int32_t *bin_number, float *xs, int32_t *sample_size);
void kdistgasabs_(float *kcdf, float *abscf, float*waveno, float *wavobs, int32_t *np_hi, int32_t *nlayers, int32_t *nwave);
void model_adjust();
void locate_pos_(float *xx, int32_t *n1, float *x1, int32_t *jj);

extern void get_input_();
extern void model_adj_();
extern void geometry_();
extern void init_speccal_();
//extern void solar_irr_pc_();
extern void tran_table_();

extern struct {
    int32_t h2o, co2, o3, n2o, co, ch4, o2, no2;
} getinput1_;

extern struct {
    char filename[FILENAME_MAX];
    int32_t dln;
} input_l2gen_;

extern struct {
    float tg_sol[NBANDS], tg_sen[NBANDS], tg_solo[NBANDS], tg_seno[NBANDS];
} tran_table_l2gen_;

extern struct {
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

extern struct {
    float wavobs[NBANDS], fwhm[NBANDS];
} getinput4_;

extern struct {
    int32_t nobs, full_calc;
    float hsurf, dlt, dlt2;
} getinput5_;

extern struct {
    float wndow1, wndow2, wp94c, wndow3, wndow4, w1p14c;
} getinput6_;

extern struct {
    int32_t nb1, nb2, nbp94, nb3, nb4, nb1p14;
} getinput7_;

extern struct {
    int32_t imn, idy, iyr, ih, im, is;
} getinput8_;

extern struct {
    float xpss, xppp;
} getinput14_;

extern struct {
    float clmvap, q;
} model_adj1_;

extern struct {
    float hp[MODELMAX], tp[MODELMAX], pp[MODELMAX], vmrp[MODELMAX];
} model_adj2_;

extern struct {
    int32_t k_plane;
    float dvap_plane, dvap_layer, dp_plane, dp_layer, clmvapp;
} model_adj3_;

extern struct {
    int32_t k_surf;
} model_adj4_;

extern struct {
    int32_t nh2o;
} init_speccal3_;

extern struct {
    int32_t ist1, ied1, ist2, ied2, istp94, iedp94;
} init_speccal6_;

extern struct {
    int32_t ist3, ied3, ist4, ied4, ist1p14, ied1p14;
} init_speccal7_;

extern struct {
    float wt1, wt2, wt3, wt4;
    int32_t ja;
} init_speccal8_;

extern struct {
    float ncv2, ncvhf2, ncvtt2;
    int32_t istrt2, iend2;
    float finst2[FINSTMAX];
} init_speccal10_;

extern struct {
    int32_t natot, nbtot, nctot, ndtot;
} init_speccal11_;

extern struct {
    float sh2o, vaptot[NH2OMAX], r0p94[NH2OMAX], r1p14[NH2OMAX], trntbl[NH2OMAX][NBANDS];
} tran_table1_;

extern struct {
    float g_vap[MODELMAX], g_other[MODELMAX], g_vap_equiv;
} geometry3_;

extern struct {
    float vap_slant_mdl;
} geometry4_;

extern struct {
    float senzn_l2, senaz_l2, solzn_l2;
    float water_vapor;
    int ja, jb;
    int splitpaths;
    float f1a, f2a, f1b, f2b;
} geometry_l2gen_;

extern struct {
    float tpvmr[81][7];
} tpvmr_init1_;

#endif /* SRC_ATREM_ATREM_H_ */
