/** @file allocate5d.c
    @brief Utility functions for allocating and freeing five-dimensional arrays of various types.

    This file was created by allocate5d.pl and should not be edited manually.
*/

#include "allocate5d.h"

#include <stdio.h>
#include <stdlib.h>

float *****allocate5d_float(size_t nq, size_t nr, size_t nz, size_t ny, size_t nx) {
    float *x_ptr = (float*) malloc(nq * nr * nz * ny * nx * sizeof(float));
    if (x_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of data block failed.n", __FILE__, __LINE__);
        return NULL;
    }
    float **y_ptr = (float**) malloc(nq * nr * nz * ny * sizeof(float*));
    if (y_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of y array failed.n", __FILE__, __LINE__);
        return NULL;
    }
    float ***z_ptr = (float***) malloc(nq * nr * nz * sizeof(float**));
    if (z_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of z array failed.n", __FILE__, __LINE__);
        return NULL;
    }
    float ****r_ptr = (float****) malloc(nq * nr * sizeof(float***));
    if (r_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of r array failed.n", __FILE__, __LINE__);
        return NULL;
    }
    float *****q_ptr = (float*****) malloc(nq * sizeof(float****));
    if (q_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of q array failed.n", __FILE__, __LINE__);
        return NULL;
    }

	for(size_t q=0; q<nq; q++) {
		for(size_t r=0; r<nr; r++) {
			for(size_t z=0; z<nz; z++) {
				for(size_t y=0; y<ny; y++) {
					y_ptr[y] = x_ptr;
					x_ptr += nx;
				}
				z_ptr[z] = y_ptr;
				y_ptr += ny;
			}
			r_ptr[r] = z_ptr;
			z_ptr += nz;
		}
		q_ptr[q] = r_ptr;
		r_ptr += nr;
	}
    return q_ptr;
}

void free5d_float(float *****p) {
    free(p[0][0][0][0]);
    free(p[0][0][0]);
    free(p[0][0]);
    free(p[0]);
    free(p);
}

