#ifndef _FILTER_H
#define _FILTER_H

#include <stdint.h>
#include "l12_parms.h"

#define FDILATE    1
#define FLTMEAN    2
#define FLTMED     3 
#define FLTRMEAN   4
#define FLTRMED    5
#define FEPSMEAN   6
#define FCLEAN     7
#define FLTRIQMEAN 8
#define FSTLIGHT   9
#define FBTDETAVG 10
#define FLTRREJECT 11
#define FTEST     12

static const char * const filter_names[] = {"",
    "dilation masking on flag ",
    "mean smoothing of Lt for band ",
    "median smoothing of LT for band ",
    "mean smoothing of Lt-Lr for band ",
    "median smoothing of Lt-Lr for band ",
    "mean epsilon smoothing for band ",
    "cleaning for flag ",
    "interquartile mean smoothing of Lt-Lr for band ",
    "straylight filter on flag ",
    "replacing BT with detector average for IR band ",
    "outlier rejection of Lt-Lr for band ",
    "test filter for band "};

typedef struct filter_struct {
    int32_t func;
    int32_t band;
    int32_t nx;
    int32_t ny;
    int32_t minfill;
    char *kernel;
} filstr;

typedef struct filter_ctl_struct {
    int32_t nfilt;
    int32_t nscan;
    int32_t npix;
    filstr f[FILTMAX];
} fctlstr;

#endif
