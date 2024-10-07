/** @file allocate4d.c
    @brief Utility functions for allocating and freeing four-dimensional arrays of various types.

    This file was created by allocate4d.pl and should not be edited manually.
*/

#include "allocate4d.h"

#include <stdio.h>
#include <stdlib.h>

int ****allocate4d_int(size_t nr, size_t nz, size_t ny, size_t nx) {
    int *x_ptr = (int*) malloc(nr * nz * ny * nx * sizeof(int));
    if (x_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of data block failed.n", __FILE__, __LINE__);
        return NULL;
    }
    int **y_ptr = (int**) malloc(nr * nz * ny * sizeof(int*));
    if (y_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of y array failed.n", __FILE__, __LINE__);
        return NULL;
    }
    int ***z_ptr = (int***) malloc(nr * nz * sizeof(int**));
    if (z_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of z array failed.n", __FILE__, __LINE__);
        return NULL;
    }
    int ****r_ptr = (int****) malloc(nr * sizeof(int***));
    if (r_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of r array failed.n", __FILE__, __LINE__);
        return NULL;
    }

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
    return r_ptr;
}

void free4d_int(int ****p) {
    free(p[0][0][0]);
    free(p[0][0]);
    free(p[0]);
    free(p);
}

float ****allocate4d_float(size_t nr, size_t nz, size_t ny, size_t nx) {
    float *x_ptr = (float*) malloc(nr * nz * ny * nx * sizeof(float));
    if (x_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of data block failed.n", __FILE__, __LINE__);
        return NULL;
    }
    float **y_ptr = (float**) malloc(nr * nz * ny * sizeof(float*));
    if (y_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of y array failed.n", __FILE__, __LINE__);
        return NULL;
    }
    float ***z_ptr = (float***) malloc(nr * nz * sizeof(float**));
    if (z_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of z array failed.n", __FILE__, __LINE__);
        return NULL;
    }
    float ****r_ptr = (float****) malloc(nr * sizeof(float***));
    if (r_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of r array failed.n", __FILE__, __LINE__);
        return NULL;
    }

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
    return r_ptr;
}

void free4d_float(float ****p) {
    free(p[0][0][0]);
    free(p[0][0]);
    free(p[0]);
    free(p);
}

double ****allocate4d_double(size_t nr, size_t nz, size_t ny, size_t nx) {
    double *x_ptr = (double*) malloc(nr * nz * ny * nx * sizeof(double));
    if (x_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of data block failed.n", __FILE__, __LINE__);
        return NULL;
    }
    double **y_ptr = (double**) malloc(nr * nz * ny * sizeof(double*));
    if (y_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of y array failed.n", __FILE__, __LINE__);
        return NULL;
    }
    double ***z_ptr = (double***) malloc(nr * nz * sizeof(double**));
    if (z_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of z array failed.n", __FILE__, __LINE__);
        return NULL;
    }
    double ****r_ptr = (double****) malloc(nr * sizeof(double***));
    if (r_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of r array failed.n", __FILE__, __LINE__);
        return NULL;
    }

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
    return r_ptr;
}

void free4d_double(double ****p) {
    free(p[0][0][0]);
    free(p[0][0]);
    free(p[0]);
    free(p);
}

