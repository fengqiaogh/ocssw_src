/**
 *  @file liboc/qaa.h
 *  @brief Compute IOPs from rrs using Quasi-Analytic Algorithm
 *  @author Paul Martinolich
 *  @author Naval Research Laboratory, Stennis Space Center, MS
 */

#ifndef _QAA_H
#define _QAA_H

enum {
    QAA_S_PARAM = 1,
    QAA_COEFS_PARAM = 3,
    QAA_APH_CHECK = 8
};

int qaa_is_initialized(void);

int qaa_init(int i410, int i440, int i490, int i555, int i640);
int qaa_set_param(int param, ...);
int qaa_v6(int nbands, double *wavel, double *Rrs, double *aw, double *bbw,
        double *rrs, double *u, double *a, double *bb,
        unsigned char *flags);
int qaa_decomp(int nbands, double *wavel, double *rrs, double *a, double *aw,
        double *adg, double *aph, unsigned char *flags);
int qaaf_v6(int nbands, float *wavel, float *Rrs, float *aw, float *bbw,
        float *rrs, float *u, float *a, float *bb,
        unsigned char *flags);
int qaaf_decomp(int nbands, float *wavel, float *rrs, float *a, float *aw,
        float *adg, float *aph, unsigned char *flags);

#endif
