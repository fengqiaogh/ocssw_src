/*
 *  W. Robinson, SAIC, 26 Oct 2001 have defs for all data types
 *  to replace the string defines
 *  24 Jul 2002  W. Robinson  add time range checking to 3 check types
 *        in thr_ctl structure
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hdf.h"
#include "mfhdf.h"

#define BAND1 0
#define BAND2 1
#define BAND3 2
#define BAND4 3
#define BAND5 4
#define BAND6 5
#define BAND7 6
#define BAND8 7

#define NSAMP "Pixels per Scan Line"
#define NSCANS "Number of Scan Lines"
#define TITLE  "Title"
#define DTYPE  "Data Type"
#define GN1SAT "Gain 1 Saturated Pixels"
#define GN2SAT "Gain 2 Saturated Pixels"
#define GN1UNSAT "Gain 1 Non-Saturated Pixels"
#define GN2UNSAT "Gain 2 Non-Saturated Pixels"
#define ZEROPIX  "Zero Pixels"

#define GAC  0
#define LAC  1
#define HRPT 2
#define LUN  3
#define SOL  4
#define IGC  5
#define TDI  6

/* define label for l1a_data sds */
#define L1ADATA "l1a_data"

/* define labels for tilt datasets */
#define NTILTS  "ntilts"
#define TILT_RANGES "tilt_ranges"
#define TILT_FLAGS "tilt_flags"

/* define labels for navigation datasets */
#define ORBVEC  "orb_vec"
#define SENMAT  "sen_mat"
#define SCANELL  "scan_ell"
#define SUNREF  "sun_ref"
#define NFLAG           "nflag"

/* gain_struct has all the information about controls for checking
 *      gain1 and gain2 data
 */

typedef struct cntl1_struct {
    int32 band;
    float32 threshold;
} cntl1_str;

typedef struct cntl2_struct {
    int32 band;
    float32 err_thresh;
    float32 cnt_thresh;
} cntl2_str;

/* WDR new structure for some of the controls */
typedef struct thr_ctl_struct {
    /*  for doing the time range check  */
    int16 trng_chk_do; /* check file for being in possible bad time range */
    /*  instrument analog telem info */
    int16 rpt_inst_ana; /* 1 to report failure for inst_ana checks */
    int16 inst_ana_do[32]; /* perform this check - the control is there */
    float32 inst_ana_lo[32]; /* low threshold */
    float32 inst_ana_hi[32]; /* high threshold */
    float32 inst_ana_pct[32]; /* % acceptable outside thresholds */
    /*  Gain check info */
    int16 rpt_gainv_chk; /* 1 to report failure for gainv checks */
    int16 gainv_chk_do[8]; /* check the gain value for GAC, LAC, SOL, LUN */
    float32 gainv_chk_pct[8]; /* % not at right gain acceptable  */
    /*  TDI check info  */
    int16 rpt_tdi_vchk; /* 1 to report failure for tdi checks */
    int16 tdiv_chk_do[8]; /* check the tdi for GAC, LAC, SOL, LUN */
    float32 tdiv_chk_pct[8]; /* % not at right tdi setting that is acceptable */
} thr_ctl_def;

/*  WDR have a convinent container for nav and time info  */
typedef struct nav_t_struct {
    float32 *orb_vec;
    float32 *sun_ref;
    float32 *sen_mat;
    float32 *scal_ell;
    int32_t *msec;
    int32 *nflag;
} nav_t_str;

#define fltabs(x)                 (x>=0 ? x : -(x))
