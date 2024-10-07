#ifndef L2BINMATCH_INPUT_H
#define L2BINMATCH_INPUT_H

#include <stdio.h>
#include <clo.h>
#include "l12_proto.h"
#include "input_struc.h"

//typedef struct binmatch_input_struct {
//
//  char    ifile  [FILENAME_MAX];
//  char    ofile  [FILENAME_MAX];
//  char    oformat[20];  // output file type
//  char    flaguse[1024];
//  char    input_parms  [4096];
//
//  char    l2prod[PRODSTRLEN];
//  char    def_l2prod[PRODSTRLEN];
//  int32_t   spixl;          /* starting pixel no. of the input (1-rel)  */
//  int32_t   epixl;          /* ending pixel no. of the input (1-rel)    */
//  int32_t   dpixl;          /* pixel subsampling increment              */
//  int32_t   sline;          /* starting line no. of the input (1-rel)   */
//  int32_t   eline;          /* ending line no. of the input (1-rel)     */
//  int32_t   dline;          /* line subsampling increment               */
//  int32_t band_shift_opt; /* band fill, 0=lin.interp. 1=bio.opt.band shift */
//  char   tgtfile[FILENAME_MAX];         /* input cal target file     */
//  float   vcal_depth;  /*  vcaltarget depth mask value */
//  int32_t vcal_min_nbin;  /* min # samples in bin to accept */
//  int32_t vcal_min_nscene;  /* min # scenes in bin to accept */
//  char    elevfile [FILENAME_MAX];
//  int32_t  deflate;
//
//
//} binmatch_instr;

int l2binmatch_input_init();
int l2binmatch_init_options(clo_optionList_t* list);
int l2binmatch_read_options(clo_optionList_t* list, int argc, char* argv[]);
int l2binmatch_input(int argc, char **argv, clo_optionList_t* list);

#endif



