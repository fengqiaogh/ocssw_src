#ifndef _AER_STRUC_H
#define _AER_STRUC_H

typedef struct aer_struct {
    int32_t length;
    int32_t npix;
    int32_t mode;
    char *data;
    int32_t *mod_min;
    int32_t *mod_max;
    float *mod_rat;
    float *taua;
} aestr;

#endif



