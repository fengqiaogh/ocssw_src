#ifndef L2STAT_H
#define L2STAT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hdf.h"
#include "mfhdf.h"

#define NPARAMS  12
#define NFLAGS  32

#define NSAMP   "Pixels per Scan Line"
#define NSCANS   "Number of Scan Lines"
#define TITLE    "Title"
#define DTYPE    "Data Type"
#define PERCENTFLAGS "Flag Percentages"

#define GAC   "GAC"
#define LAC   "LAC"
#define HRPT   "HRPT"

#define MASKNAMES "Mask Names"
#define L2FLAGS  "l2_flags"

/*
 * cntl2_struct will hold requested band threshold values for gain1, gain2,
 *  	and zero pixels
 */

typedef struct cntl_struct {
    int32 param;
    float32 err_thresh;
    float32 low_thresh;
    float32 high_thresh;
} cntl_str;

typedef struct flag_struct {
    int32 flag;
    float32 err_low_thresh;
    float32 err_high_thresh;
} flag_str;

/*
 *  Following define returns an absolute value of the given number
 */

#define fltabs(x)                 (x>=0 ? x : -(x))

#endif /* L2STAT_H */
