#ifndef GETCAL_H
#define GETCAL_H
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "mfhdf.h"

#define BANDS 8
#define DETS  4
#define GAINS 4
#define KNEE  4
#define TIME "Time"
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

#endif  /* GETCAL_H */
