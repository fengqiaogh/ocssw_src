/** @file allocate2d.c
	@brief Utility functions for allocating and freeing two-dimensional arrays of various types.

	This file was created by allocate2d.pl and should not be edited manually.
*/

#include "allocate2d.h"

#include <stdio.h>
#include <stdlib.h>


char **allocate2d_char(size_t h, size_t w){
    char **p = (char**) malloc(h * sizeof(char*));
    if (p == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation failed.n", __FILE__, __LINE__);
        return NULL;
    }
    p[0] = (char*) malloc(w * h * sizeof(char));
    if (p[0] == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation failed.n", __FILE__, __LINE__);
        free(p);
        return NULL;
    }
    for (size_t i = 1; i < h; i++) {
        p[i] = &(p[0][i * w]);
    }
    return p;
}
void free2d_char(char **p) {
    free(p[0]);
    free(p);
}

unsigned char **allocate2d_uchar(size_t h, size_t w){
    unsigned char **p = (unsigned char**) malloc(h * sizeof(unsigned char*));
    if (p == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation failed.n", __FILE__, __LINE__);
        return NULL;
    }
    p[0] = (unsigned char*) malloc(w * h * sizeof(unsigned char));
    if (p[0] == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation failed.n", __FILE__, __LINE__);
        free(p);
        return NULL;
    }
    for (size_t i = 1; i < h; i++) {
        p[i] = &(p[0][i * w]);
    }
    return p;
}
void free2d_uchar(unsigned char **p) {
    free(p[0]);
    free(p);
}

signed char **allocate2d_schar(size_t h, size_t w){
    signed char **p = (signed char**) malloc(h * sizeof(signed char*));
    if (p == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation failed.n", __FILE__, __LINE__);
        return NULL;
    }
    p[0] = (signed char*) malloc(w * h * sizeof(signed char));
    if (p[0] == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation failed.n", __FILE__, __LINE__);
        free(p);
        return NULL;
    }
    for (size_t i = 1; i < h; i++) {
        p[i] = &(p[0][i * w]);
    }
    return p;
}
void free2d_schar(signed char **p) {
    free(p[0]);
    free(p);
}

short **allocate2d_short(size_t h, size_t w){
    short **p = (short**) malloc(h * sizeof(short*));
    if (p == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation failed.n", __FILE__, __LINE__);
        return NULL;
    }
    p[0] = (short*) malloc(w * h * sizeof(short));
    if (p[0] == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation failed.n", __FILE__, __LINE__);
        free(p);
        return NULL;
    }
    for (size_t i = 1; i < h; i++) {
        p[i] = &(p[0][i * w]);
    }
    return p;
}
void free2d_short(short **p) {
    free(p[0]);
    free(p);
}

int **allocate2d_int(size_t h, size_t w){
    int **p = (int**) malloc(h * sizeof(int*));
    if (p == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation failed.n", __FILE__, __LINE__);
        return NULL;
    }
    p[0] = (int*) malloc(w * h * sizeof(int));
    if (p[0] == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation failed.n", __FILE__, __LINE__);
        free(p);
        return NULL;
    }
    for (size_t i = 1; i < h; i++) {
        p[i] = &(p[0][i * w]);
    }
    return p;
}
void free2d_int(int **p) {
    free(p[0]);
    free(p);
}

float **allocate2d_float(size_t h, size_t w){
    float **p = (float**) malloc(h * sizeof(float*));
    if (p == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation failed.n", __FILE__, __LINE__);
        return NULL;
    }
    p[0] = (float*) malloc(w * h * sizeof(float));
    if (p[0] == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation failed.n", __FILE__, __LINE__);
        free(p);
        return NULL;
    }
    for (size_t i = 1; i < h; i++) {
        p[i] = &(p[0][i * w]);
    }
    return p;
}
void free2d_float(float **p) {
    free(p[0]);
    free(p);
}

double **allocate2d_double(size_t h, size_t w){
    double **p = (double**) malloc(h * sizeof(double*));
    if (p == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation failed.n", __FILE__, __LINE__);
        return NULL;
    }
    p[0] = (double*) malloc(w * h * sizeof(double));
    if (p[0] == NULL) {
        fprintf(stderr, "-E- %s line %d: Memory allocation failed.n", __FILE__, __LINE__);
        free(p);
        return NULL;
    }
    for (size_t i = 1; i < h; i++) {
        p[i] = &(p[0][i * w]);
    }
    return p;
}
void free2d_double(double **p) {
    free(p[0]);
    free(p);
}

