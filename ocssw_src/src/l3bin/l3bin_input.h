#ifndef _INPUT_STR_H
#define _INPUT_STR_H

#include <stdio.h>
#include <stdint.h>
#include "clo.h"

typedef struct input_struct {
    char infile [FILENAME_MAX];
    char ofile [FILENAME_MAX];
    char pfile [FILENAME_MAX];
    char out_parm[16384];
    char tflag;
    char parms[16384];
    char pversion[16];

    int32_t syear;
    int32_t sday;
    int32_t eyear;
    int32_t eday;

    int32_t sorbit;
    int32_t eorbit;

    int32_t reduce_fac;
    char resolve[10];

    int32_t noext;

    char merged[4096];

    float loneast;
    float lonwest;
    float latnorth;
    float latsouth;

    int32_t verbose;
    int32_t unit_wgt;
    int32_t median;
    int32_t union_bins;

    uint32_t deflate;
    char oformat [20];

    char composite_prod[1024];
    char composite_scheme[1024];
    char doi[1024];
} instr;

int l3bin_input(int argc, char **argv, instr *input, const char* prog, const char* version);
//int l3bin_init_options(clo_optionList_t* list, const char* prog, const char* version);

#endif



