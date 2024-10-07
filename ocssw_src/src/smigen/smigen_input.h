#ifndef _INPUT_STR_H
#define _INPUT_STR_H

#include <stdio.h>
#include <stdint.h>

typedef struct input_struct {
    char ifile [FILENAME_MAX];
    char ofile [FILENAME_MAX];
    char pfile [FILENAME_MAX];
    char palfile[FILENAME_MAX];
    char parms [4096];

    char pversion[255];
    char prod[255];
    int32_t stype;
    int32_t meas;
    float datamin;
    float datamax;
    float lonwest;
    float loneast;
    float latnorth;
    float latsouth;
    char resolution[32];
    char projection[8];
    char palette[768];
    int32_t gap_fill;
    float seam_lon;
    char proddesc[128];
    char units[64];
    char precision[4];
    uint32_t minobs;

    int32_t deflate;
    char oformat[20];
} instr;

#endif



