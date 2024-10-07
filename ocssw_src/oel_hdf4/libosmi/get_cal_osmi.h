#ifndef GETCAL_H
#define GETCAL_H
#include <stdio.h>
#include <stdlib.h>

#define BANDS 8
#define DETS  4
#define GAINS 4
#define KNEE  4
#define TIME "CalTime"
#define RDERR  -2
#define TMERR  -3
#define BUFERR -4

#define TEMPS   "temps"
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

static const char gain_flds[] = {
    "egain0,egain1,egain2,egain3,egain4,egain5,egain6,egain7,eoffset"
};

static const char mirror_flds[] = {
    "mirror0,mirror1,mirror2,mirror3,mirror4,mirror5,mirror6,mirror7"
};
static const char t_const_flds[] = {
    "t_const0,t_const1,t_const2,t_const3,t_const4,t_const5,t_const6,t_const7"
};
static const char t_linear_flds[] = {
    "t_linear0,t_linear1,t_linear2,t_linear3,t_linear4,t_linear5,t_linear6,t_linear7"
};
static const char t_quadratic_flds[] = {
    "t_quadratic0,t_quadratic1,t_quadratic2,t_quadratic3,t_quadratic4,t_quadratic5,t_quadratic6,t_quadratic7"
};
static const char temp_flds[] = {
    "tempc0,tempc1,tempc2,tempc3,tempc4,tempc5,tempc6,tempc7"
};
static const char slope_flds[] = {
    "slope0_00,slope0_08,slope0_16,slope0_24,slope0_32,slope0_40,slope0_48,slope0_56,slope0_64,slope0_72,slope0_80,slope0_88,slope1_00,slope1_08,slope1_16,slope1_24,slope1_32,slope1_40,slope1_48,slope1_56,slope1_64,slope1_72,slope1_80,slope1_88,slope2_00,slope2_08,slope2_16,slope2_24,slope2_32,slope2_40,slope2_48,slope2_56,slope2_64,slope2_72,slope2_80,slope2_88,slope3_00,slope3_08,slope3_16,slope3_24,slope3_32,slope3_40,slope3_48,slope3_56,slope3_64,slope3_72,slope3_80,slope3_88,slope4_00,slope4_08,slope4_16,slope4_24,slope4_32,slope4_40,slope4_48,slope4_56,slope4_64,slope4_72,slope4_80,slope4_88,slope5_00,slope5_08,slope5_16,slope5_24,slope5_32,slope5_40,slope5_48,slope5_56,slope5_64,slope5_72,slope5_80,slope5_88,slope6_00,slope6_08,slope6_16,slope6_24,slope6_32,slope6_40,slope6_48,slope6_56,slope6_64,slope6_72,slope6_80,slope6_88,slope7_00,slope7_08,slope7_16,slope7_24,slope7_32,slope7_40,slope7_48,slope7_56,slope7_64,slope7_72,slope7_80,slope7_88"
};
static const char dc_flds[] = {
    "dc0_00,dc0_08,dc0_16,dc0_24,dc0_32,dc0_40,dc0_48,dc0_56,dc0_64,dc0_72,dc0_80,dc0_88,dc1_00,dc1_08,dc1_16,dc1_24,dc1_32,dc1_40,dc1_48,dc1_56,dc1_64,dc1_72,dc1_80,dc1_88,dc2_00,dc2_08,dc2_16,dc2_24,dc2_32,dc2_40,dc2_48,dc2_56,dc2_64,dc2_72,dc2_80,dc2_88,dc3_00,dc3_08,dc3_16,dc3_24,dc3_32,dc3_40,dc3_48,dc3_56,dc3_64,dc3_72,dc3_80,dc3_88,dc4_00,dc4_08,dc4_16,dc4_24,dc4_32,dc4_40,dc4_48,dc4_56,dc4_64,dc4_72,dc4_80,dc4_88,dc5_00,dc5_08,dc5_16,dc5_24,dc5_32,dc5_40,dc5_48,dc5_56,dc5_64,dc5_72,dc5_80,dc5_88,dc6_00,dc6_08,dc6_16,dc6_24,dc6_32,dc6_40,dc6_48,dc6_56,dc6_64,dc6_72,dc6_80,dc6_88,dc7_00,dc7_08,dc7_16,dc7_24,dc7_32,dc7_40,dc7_48,dc7_56,dc7_64,dc7_72,dc7_80,dc7_88"
};
static const char sm_flds[] = {
    "sm0,sm12,sm24,sm36,sm48,sm60,sm72,sm84,sm96,sm108,sm120,sm132,sm144,sm156,sm168,sm180,sm192,sm204,sm216,sm228,sm240,sm252,sm264,sm276,sm288,sm300,sm312,sm324,sm336,sm348,sm360,sm372,sm384,sm396,sm408,sm420,sm432,sm444,sm456,sm468,sm480,sm492,sm504,sm516,sm528,sm540,sm552,sm564,sm576,sm588,sm600,sm612,sm624,sm636,sm648,sm660,sm672,sm684,sm696,sm708,sm720,sm732,sm744,sm756,sm768,sm780,sm792,sm804,sm816,sm828,sm840,sm852,sm864,sm876,sm888,sm900,sm912,sm924,sm936,sm948,sm960,sm972,sm984,sm996,sm1008,sm1020,sm1032"
};


#endif  /* GETCAL_H */
