#ifndef _SST_FLAGS_H
#define _SST_FLAGS_H

/* flag bit settings */
/* avhrr needs 3 bits that modis doesn't and since avhrr doesn't have */
/* a 4 um product, we can share those bits (40, 80, and 400) */
#define SSTF_ISMASKED  0x0001    /* Pixel aready masked	 M1B3 */
#define SSTF_BTBAD  0x0002    /* Brights are bad		 M1B4 */
#define SSTF_BTRANGE  0x0004    /* Brights are out-of-range	 M1B1 */
#define SSTF_BTDIFF  0x0008    /* Brights are too different	  */
#define SSTF_SSTRANGE  0x0010    /* SST outside valid range	 M2B4 */
#define SSTF_SSTREFDIFF  0x0020    /* Different from reference	 M1B8 */
#define SSTF_SST4DIFF  0x0040    /* Different from SST4	  */
#define SSTF_SST3DIFF  0x0040    /* Different from SST3	  */
#define SSTF_SUNLIGHT    0x0040           /* Stray sunlight test (AVHRR) M2B2 */
#define SSTF_SST4VDIFF  0x0080    /* Very different from SST4	  */
#define SSTF_SST3VDIFF  0x0080    /* Very different from SST3	  */
#define SSTF_ASCEND      0x0080           /* ASCENDING           (AVHRR) M2B6 */
#define SSTF_BTNONUNIF  0x0100    /* BT window nonuniform	 M1B5 */
#define SSTF_BTVNONUNIF  0x0200    /* BT window very nonuniform	 M1B6 */
#define SSTF_BT4REFDIFF  0x0400    /* 4um BT diff relative to ref  */
#define SSTF_GLINT       0x0400           /* Glint_coef>glintmax (AVHRR) M2B8 */
#define SSTF_REDNONUNIF  0x0800    /* Red band spatial nonuniform  */
#define SSTF_HISENZ  0x1000    /* Sensor zenith high		 M1B7 */
#define SSTF_VHISENZ  0x2000    /* Sensor zenith very high	 M2B1 */
#define SSTF_SSTREFVDIFF 0x4000    /* Very diff from reference	  */
#define SSTF_CLOUD  0x8000    /* Cloud test (tree test)	 M1B2 */

/* This set would match original pathfinder processing */
//static char *avhrr_sst_flag_lname[NSSTFLAGS] = {"BTRANGE", 
//                                      "CLOUD",
//                                      "ISMASKED",
//                                      "BTBAD",
//                                      "BTNONUNIF",
//                                      "BTVNONUNIF",
//                                      "HISENZ",
//                                      "SSTREFDIFF",
//                                      "VHISENZ",
//                                      "SUNLIGHT",
//                                      "SPARE",
//                                      "SSTRANGE",
//                                      "SPARE",
//                                      "ASCEND",
//                                      "SPARE",
//				      "GLINT"};
/* This set fits the avhrr bits into the same slots as the modis sst flags */
/* with some 4 um specific bits shared with avhrr only tests */
static const char * const sst_flag_lname[NSSTFLAGS] = {"ISMASKED",
    "BTBAD",
    "BTRANGE",
    "BTDIFF",
    "SSTRANGE",
    "SSTREFDIFF",
    "SST4DIFF",
    "SST4VDIFF",
    "BTNONUNIF",
    "BTVNONUNIF",
    "BT4REFDIFF",
    "REDNONUNIF",
    "HISENZ",
    "VHISENZ",
    "SSTREFVDIFF",
    "CLOUD"};
static const char * const avhrr_sst_flag_lname[NSSTFLAGS] = {"ISMASKED",
    "BTBAD",
    "BTRANGE",
    "BTDIFF",
    "SSTRANGE",
    "SSTREFDIFF",
    "SUNLIGHT",
    "ASCEND",
    "BTNONUNIF",
    "BTVNONUNIF",
    "GLINT",
    "REDNONUNIF",
    "HISENZ",
    "VHISENZ",
    "SSTREFVDIFF",
    "CLOUD"};
static const char * const viirs_sst_flag_lname[NSSTFLAGS] = {"ISMASKED",
    "BTBAD",
    "BTRANGE",
    "BTDIFF",
    "SSTRANGE",
    "SSTREFDIFF",
    "SST3DIFF",
    "SST3VDIFF",
    "BTNONUNIF",
    "BTVNONUNIF",
    "SPARE",
    "REDNONUNIF",
    "HISENZ",
    "VHISENZ",
    "SSTREFVDIFF",
    "CLOUD"};
static const char * const qual_sst_flag_lname[NQSSTFLAGS] = {"BEST",
    "GOOD",
    "QUESTIONABLE",
    "BAD",
    "NOTPROCESSED"};
#endif
