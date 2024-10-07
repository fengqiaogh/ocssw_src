#ifndef L2MAPGEN_INPUT_H
#define L2MAPGEN_INPUT_H

#include <stdio.h>
#include <clo.h>

typedef struct input_struct {
    char ifile [FILENAME_MAX];
    char ofile [FILENAME_MAX];
    char palfile[FILENAME_MAX];
    char palette_dir[FILENAME_MAX];
    char product_table[FILENAME_MAX];
    char flaguse[1024];
    char parms [4096];

    char prod[255];
    int32_t stype;
    float datamin;
    float datamax;
    float west;
    float east;
    float north;
    float south;
    int32_t width;
    float threshold;
    int mask;
    int quality;
    int apply_pal;
    char palette[768];
    int outmode;
} instr;

int l2mapgen_init_options(clo_optionList_t* list);
int l2mapgen_read_options(clo_optionList_t* list, int argc, char* argv[]);
int l2mapgen_input(int argc, char **argv, clo_optionList_t* list, instr *input);

#endif



