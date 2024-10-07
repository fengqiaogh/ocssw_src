/** @file allocate3d.c
	@brief Utility functions for allocating and freeing three-dimensional arrays of various types.

	This file was created by allocate3d.pl and should not be edited manually.
*/

#include "allocate3d.h"

#include <stdio.h>
#include <stdlib.h>


short ***allocate3d_short(size_t nz, size_t ny, size_t nx){
    short *x_ptr = (short*) malloc(nz * ny * nx * sizeof(short));
    if (x_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of data block failed.n", __FILE__, __LINE__);
        return NULL;
    }
    short **y_ptr = (short**) malloc(nz * ny * sizeof(short*));
    if (y_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of y array failed.n", __FILE__, __LINE__);
        return NULL;
    }
    short ***z_ptr = (short***) malloc(nz * sizeof(short**));
    if (z_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of z array failed.n", __FILE__, __LINE__);
        return NULL;
    }
    for(size_t z=0; z<nz; z++) {
        for(size_t y=0; y<ny; y++) {
            y_ptr[y] = x_ptr;
            x_ptr += nx;
        }
        z_ptr[z] = y_ptr;
        y_ptr += ny;
    }
    return z_ptr;
}
void free3d_short(short ***p) {
    free(p[0][0]);
    free(p[0]);
    free(p);
}

int ***allocate3d_int(size_t nz, size_t ny, size_t nx){
    int *x_ptr = (int*) malloc(nz * ny * nx * sizeof(int));
    if (x_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of data block failed.n", __FILE__, __LINE__);
        return NULL;
    }
    int **y_ptr = (int**) malloc(nz * ny * sizeof(int*));
    if (y_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of y array failed.n", __FILE__, __LINE__);
        return NULL;
    }
    int ***z_ptr = (int***) malloc(nz * sizeof(int**));
    if (z_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of z array failed.n", __FILE__, __LINE__);
        return NULL;
    }
    for(size_t z=0; z<nz; z++) {
        for(size_t y=0; y<ny; y++) {
            y_ptr[y] = x_ptr;
            x_ptr += nx;
        }
        z_ptr[z] = y_ptr;
        y_ptr += ny;
    }
    return z_ptr;
}
void free3d_int(int ***p) {
    free(p[0][0]);
    free(p[0]);
    free(p);
}

float ***allocate3d_float(size_t nz, size_t ny, size_t nx){
    float *x_ptr = (float*) malloc(nz * ny * nx * sizeof(float));
    if (x_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of data block failed.n", __FILE__, __LINE__);
        return NULL;
    }
    float **y_ptr = (float**) malloc(nz * ny * sizeof(float*));
    if (y_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of y array failed.n", __FILE__, __LINE__);
        return NULL;
    }
    float ***z_ptr = (float***) malloc(nz * sizeof(float**));
    if (z_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of z array failed.n", __FILE__, __LINE__);
        return NULL;
    }
    for(size_t z=0; z<nz; z++) {
        for(size_t y=0; y<ny; y++) {
            y_ptr[y] = x_ptr;
            x_ptr += nx;
        }
        z_ptr[z] = y_ptr;
        y_ptr += ny;
    }
    return z_ptr;
}
void free3d_float(float ***p) {
    free(p[0][0]);
    free(p[0]);
    free(p);
}

double ***allocate3d_double(size_t nz, size_t ny, size_t nx){
    double *x_ptr = (double*) malloc(nz * ny * nx * sizeof(double));
    if (x_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of data block failed.n", __FILE__, __LINE__);
        return NULL;
    }
    double **y_ptr = (double**) malloc(nz * ny * sizeof(double*));
    if (y_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of y array failed.n", __FILE__, __LINE__);
        return NULL;
    }
    double ***z_ptr = (double***) malloc(nz * sizeof(double**));
    if (z_ptr == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation of z array failed.n", __FILE__, __LINE__);
        return NULL;
    }
    for(size_t z=0; z<nz; z++) {
        for(size_t y=0; y<ny; y++) {
            y_ptr[y] = x_ptr;
            x_ptr += nx;
        }
        z_ptr[z] = y_ptr;
        y_ptr += ny;
    }
    return z_ptr;
}
void free3d_double(double ***p) {
    free(p[0][0]);
    free(p[0]);
    free(p);
}

