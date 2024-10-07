#ifndef tjsm_pml_bright
#define tjsm_pml_bright
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Function declarations */
int bp_gen(int sensorID, float *geom, float *rho_rc, float *td, int32_t nwave, float *rho_w, float *rho_a, float *n_ret, float *b_ret, int verbose);
double spm_rho(double ref, int band);

#endif
