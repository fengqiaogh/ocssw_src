/*
  $Header: /app/shared/RCS/irix-5.2/seawifsd/src/hdfio/Shared.V4.4/L012_Util/util/hdf/l2_flags.h,v 1.1 1995/12/07 16:20:23 seawifsd Exp seawifsd $
  $Log: l2_flags.h,v $
  Revision 1.1  1995/12/07 16:20:23  seawifsd
  Initial revision

 */
#ifndef L2_FLAGS_SW_H_
#define L2_FLAGS_SW_H_

#define TOTAL_FLAG_NAMES 16
#define LEN_FLAG_NAMES  19
#define LEN_FLAG_NAMES_P1 (LEN_FLAG_NAMES + 1)
#define MAXLEN_FLAG_NAMES (TOTAL_FLAG_NAMES * LEN_FLAG_NAMES_P1)

static const char * const l2_flags_name[][3] = {
    {"EPSILON1", "f01_name", "atmospheric correction algorithm failure"},
    {"LAND1", "f02_name", "land"},
    {"ANCIL1", "f03_name", "missing ancillary data"},
    {"SUNGLINT1", "f04_name", "Sun glint"},
    {"HIGHLT1", "f05_name", "total radiance greater than knee value"},
    {"SATZEN1", "f06_name", "large spacecraft zenith angle"},
    {"COASTZ1", "f07_name", "shallow water"},
    {"NEGLW1", "f08_name", "negative water-leaving radiance"},
    {"STRAYLIGHT1", "f09_name", "stray light"},
    {"CLDICE1", "f10_name", "cloud and ice"},
    {"COCCOLITH1", "f11_name", "coccolithophores"},
    {"TURBIDW1", "f12_name", "turbid, case-2 water"},
    {"SOLZEN1", "f13_name", "large solar zenith angle"},
    {"HIGHTAU1", "f14_name", "high aerosol concentration"},
    {"LOWLW1", "f15_name", "low water-leaving radiance at 555 nm"},
    {"CHLOR1", "f16_name", "chlorophyll algorithm failure"}
};

#define MASK_EPSILON1  1
#define MASK_LAND1  2
#define MASK_ANCIL1  4
#define MASK_SUNGLINT1  8
#define MASK_HIGHLT1  16
#define MASK_SATZEN1  32
#define MASK_COSTZ1  64
#define MASK_NEGLW1  128
#define MASK_STRAYLIGHT1 256
#define MASK_CLDICE1  512
#define MASK_COCCOLITH1  1024
#define MASK_TURBIDW1  2048
#define MASK_SOLZEN1  4096
#define MASK_HIGHTAU1  8192
#define MASK_LOWLW1  16384
#define MASK_CHLOR1  32768

#include "l2_flags_proto.h"

#endif /* L2_FLAGS_SW_H_ */
