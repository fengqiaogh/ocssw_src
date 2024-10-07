#ifndef _INPUT_STR_H
#define _INPUT_STR_H

#include <stdio.h>
#include "clo.h"

#define DEF_FLAG "ATMFAIL,LAND,HILT,HISATZEN,STRAYLIGHT,CLDICE,COCCOLITH,LOWLW,CHLFAIL,CHLWARN,NAVWARN,ABSAER,MAXAERITER,ATMWARN,HISOLZEN,NAVFAIL,FILTER"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct input_struct {
    char infile [FILENAME_MAX];
    char ofile [FILENAME_MAX];
    char pfile [FILENAME_MAX];
    char fileuse [FILENAME_MAX];
    char flaguse[2048];
    char l3bprod[2048];
    char prodtype[32];
    char qual_prod[2048];
    char composite_prod[2048];
    char composite_scheme[2048];
    char pversion[16];
    char suite [32];
    char output_wavelengths[2048];
    char output_product_names[2048];
    char parms [4096];

    int32_t sday;
    int32_t eday;
    char resolve[4];
    int32_t rowgroup;
    int32_t meminfo;
    int32_t dcinfo;
    int32_t night;
    int32_t verbose;
    int32_t minobs;
    float deltaeqcross;
    int32_t deflate;

    float latsouth;
    float latnorth;
    float lonwest;
    float loneast;

    uint8_t qual_max;
    std::vector<std::string> files;

    int area_weighting;
    char doi[1024];
    
} instr;

int l2bin_input(int argc, char **argv, instr *input, const char* prog, const char* version);
int l2bin_init_options(clo_optionList_t* list, const char* prog, const char* version);

#ifdef __cplusplus
}
#endif
#endif
