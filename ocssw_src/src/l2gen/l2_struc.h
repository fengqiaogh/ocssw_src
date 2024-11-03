#ifndef _L2_STRUC_H
#define _L2_STRUC_H

#include "input_struc.h"
#include "target_struc.h"
#include <l1.h>

typedef struct l2_struct {
    l1str *l1rec;
    int32_t length;
    char *data;

    // var[npix]
    int32_t *num_iter;
    int32_t *aermodmin;
    int32_t *aermodmax;
    int32_t *aermodmin2;
    int32_t *aermodmax2;

    float *chl;
    float *eps; // NIR aerosol reflectance ratio (single scattering)
    float *aerratio;
    float *aerratio2;
    float *aerindex;

    // var[npix*nbands]
    float *taua; // aerosol optical thickness
    float *La; // aerosol radiance1
    float *Lw; // water-leaving radiance
    float *nLw; // normalized water-leaving radiance
    float *nLw_unc;
    float *brdf; //bi-direction reflectance function
    float *Rrs; //Remote sensing reflectance
    float *Rrs_unc;
    float *chi2; //chi square from the spectral match in mbac
    float *chl_unc;
    float *covariance_matrix;// covariance matrix for Rrs
    float *outband_correction; //square bandpass correction for Rrs
    float *a; //absoprtion coefficient
    float *bb; //backscattering coefficient

    // allocated or set later
    int32_t *bindx;
    float *sst;
    float *Rrs_raman;
    tgstr *tgrec;

} l2str;

#endif
