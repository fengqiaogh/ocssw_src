#ifndef _READ_L3BIN_H
#define _READ_L3BIN_H

#include "l12_proto.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct l3bin_struct {
    char *file;
    int32_t sensorID;
    int32_t nbins;
    int32_t nprods;
    int32_t nrows;
    int32_t *numbin;
    int32_t *basebin;
    int32_t *bins;
    int32_t *nobs;
    int32_t *nscenes;
    float *chl;
    float *tau;
    float **data; // This is to contain nLw or Rrs
    char **prods;
    int32_t nwave;
    int hasRrs;
    float *wavelengths;
} l3binstr;

//int32_t read_l3bin (char *file, char prodlist_wt[MAXPROD][32], int32_t nprods_wt, l3binstr *l3bin);
int32_t read_l3bin(char *file, l3binstr *l3bin, int32_t nbands);

#ifdef __cplusplus
}
#endif

#endif
