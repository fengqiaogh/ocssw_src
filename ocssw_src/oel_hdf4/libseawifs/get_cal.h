#ifndef GETCAL_H
#define GETCAL_H
#include <stdio.h>
#include <stdlib.h>

#define BANDS 8
#define DETS  4
#define GAINS 4
#define KNEE  4
#define TIME "Time"
#define RDERR  -2
#define TMERR  -3
#define BUFERR -4

#define TEMPS           "temps"
#define FPTEMPS  "fp_temps"
#define TDILIST  "TDI_list"
#define SCANMOD  "scan_mod"

#define REFYEAR  "Reference Year"
#define REFDAY  "Reference Day"
#define REFMIN  "Reference Minute"
#define ENTRY_YEAR "entry_year"
#define ENTRY_DAY "entry_day"
#define SYEAR  "syear"
#define SDAY  "sday"
#define SMSEC  "smsec"
#define EYEAR  "eyear"
#define EDAY  "eday"
#define EMSEC  "emsec"
#define GAC             "gac"

#define CAL_CLASS "Calibration"

static const char slp_flds[] = {
    "g1d1,g1d2,g1d3,g1d4,g2d1,g2d2,g2d3,g2d4,g3d1,g3d2,g3d3,g3d4,g4d1,g4d2,g4d3,g4d4"
};

static const char * const slp_names[] = {
    "B1Slopes", "B2Slopes", "B3Slopes", "B4Slopes", "B5Slopes", "B6Slopes",
    "B7Slopes", "B8Slopes"
};

static const char * const parm_names[] = {
    "B1Parms", "B2Parms", "B3Parms", "B4Parms", "B5Parms", "B6Parms",
    "B7Parms", "B8Parms"
};


/* These fields are for new exponential cal table of January 2001 */
static const char OFFSET_FLDS[] = {
    "g1offs1,g1offs2,g1offs3,g1offs4,g2offs1,g2offs2,g2offs3,g2offs4,g3offs1,g3offs2,g3offs3,g3offs4,g4offs1,g4offs2,g4offs3,g4offs4"
};

static const char TFACTOR_FLDS_NEW[] = {"t_const,t_linear_1,t_exponential_1,t_linear_2,t_exponential_2,cal_offs,inst_tcorr,inst_tref,fp_tcorr,fp_tref"};

static const char TFACTOR_FLDS[] = {"t_const,t_linear_1,t_exponential_1,t_linear_2,t_exponential_2,cal_offs"};

static const char MIRROR_FLDS[] = {"ms1_const,ms1_linear_1,ms1_exponential_1,ms1_linear_2,ms1_exponential_2,ms2_const,ms2_linear_1,ms2_exponential_1,ms2_linear_2,ms2_exponential_2"};

#endif  /* GETCAL_H */
