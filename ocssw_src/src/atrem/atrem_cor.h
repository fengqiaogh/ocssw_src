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

#define NBANDS  1024 /* maximum number of bands */
#define TBLMAX    60
#define TBLMAXM1  (TBLMAX-1)
#define FINSTMAX 100
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
    /*   number of channels of the four AVIRIS spectrometers. */
    int32_t natot, nbtot, nctot, ndtot;
    /* Relative weights for the four window     *
     *      channels used in channel-ratioing calculations */
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
} paramstr;

int get_atrem(double *yy, paramstr P, double *, double *, double *, double *, double *, double *, double *, double *, int32_t *, int32_t *);
int32_t hunt(float *xx, int32_t n, double x, int32_t jlo);
//DIMENSION VAPTOT(60), R0P94(60), R1P14(60), TRNTBL(1024,60)
//DIMENSION FINST2(100)
//COMMON /GETINPUT5/ NOBS,HSURF,DLT,DLT2
//COMMON /GETINPUT7/ NB1,NB2,NBP94,NB3,NB4,NB1P14
//COMMON /INIT_SPECCAL3/ NH2O
//COMMON /INIT_SPECCAL6/ IST1,IED1,IST2,IED2,ISTP94,IEDP94
//COMMON /INIT_SPECCAL7/ IST3,IED3,IST4,IED4,IST1P14,IED1P14
//COMMON /INIT_SPECCAL8/ WT1,WT2,WT3,WT4,JA
//COMMON /INIT_SPECCAL10/ NCV2,NCVHF2,NCVTT2,ISTRT2,IEND2,FINST2
//COMMON /INIT_SPECCAL11/ NATOT,NBTOT,NCTOT,NDTOT
//COMMON /TRAN_TABLE1/ SH2O,VAPTOT,R0P94,R1P14,TRNTBL
//COMMON /GEOMETRY3/ G_VAP, G_OTHER, G_VAP_EQUIV, VAP_SLANT_MDL
extern void get_input_();
extern void model_adj_();
extern void geometry_();
extern void init_speccal();
extern void solar_irr_pc();
extern void tran_table();

extern struct {
    int32_t nobs;
    float hsurf, dlt, dlt2;
} getinput5_;

extern struct {
    int32_t nb1, nb2, nbp94, nb3, nb4, nb1p14;
} getinput7_;

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
    float wt1, wt2, wt3, wt4, ja;
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
    float sh2o, vaptot[TBLMAX], r0p94[TBLMAX], r1p14[TBLMAX], trntbl[TBLMAX][NBANDS];
} tran_table1_;

extern struct {
    float g_vap, g_other, g_vap_equiv, vap_slant_mdl;
} geometry3_;
#endif /* SRC_ATREM_ATREM_H_ */
