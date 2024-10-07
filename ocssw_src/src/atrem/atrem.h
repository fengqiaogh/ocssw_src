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

#define NBANDS  1024 /* maximum number of pixels in a scan line */
#define TBLMAX    60
#define TBLMAXM1  (TBLMAX-1)
#define FINSTMAX 100

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
    double r0p94[TBLMAX]; //ratio for the 0.94 um H2O absorption band
    double r1p14[TBLMAX]; //ratio for the 1.14 um H2O absorption band
    double finst2[FINSTMAX]; // some kind of smoothing factor calculated in INIT_SPECCAL only used for AVIRIS?
    /* TOTAL SUN-SURFACE-SENSOR PATH WATER VAPOR IN UNITS OF CM */
    double vaptot[TBLMAX];
    /* total transmittances of all gases that match the     *
     *                         resolutions of imaging spectrometers */
    double trntbl[TBLMAX][NBANDS];
} paramstr;

int get_atrem(double *yy, paramstr P, double *, double *, double *, double *, double *, double *, double *, double *, int32_t *, int32_t *);
int32_t hunt(double *xx, int32_t n, double x, int32_t jlo);

#endif /* SRC_ATREM_ATREM_H_ */
