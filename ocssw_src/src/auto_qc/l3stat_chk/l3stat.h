#ifndef L3STAT_H
#define L3STAT_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NPARAMS  11
#define NFLAGS  32
#define BUFSZ    1000
#define MAXVAL 255
#define SET 1

#define NSAMP   "Pixels per Scan Line"
#define NSCANS   "Number of Scan Lines"
#define TITLE    "Title"
#define DTYPE    "Data Type"
#define PERCENTFLAGS "Flag Percentages"

#define GAC   "GAC"
#define LAC   "LAC"
#define HRPT   "HRPT"

#define MASKNAMES "Mask Names"
#define L3FLAGS  "l3_flags"

/*
 * cntl_struct will hold requested parameter threshold values 
 */

typedef struct cntl_struct {
    int32_t param;
    float err_thresh;
    float low_thresh;
    float high_thresh;
} cntl_str;

typedef struct clim_struct {
    int32_t param;
    float thresh1H;
    float thresh1L;
    float thresh2H;
    float thresh2L;
    float thresh3H;
    float thresh3L;
    char climfile[255];
} clim_str;

/*
 *  Following define returns an absolute value of the given number
 */

#define fltabs(x)                 (x>=0 ? x : -(x))

#endif /* L3STAT_H */
