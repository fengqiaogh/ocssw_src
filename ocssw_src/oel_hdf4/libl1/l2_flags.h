#ifndef _L2_FLAGS_H
#define _L2_FLAGS_H

#include <filehandle.h>

#ifdef __cplusplus
extern "C" {
#endif


#define ATMFAIL             1
#define LAND                2
#define PRODWARN            4
#define HIGLINT             8
#define HILT               16
#define HISATZEN           32
#define COASTZ             64
#define SPARE1            128
#define STRAYLIGHT        256
#define CLOUD             512
#define COCCOLITH        1024
#define TURBIDW          2048
#define HISOLZEN         4096
#define SPARE2          8192
#define LOWLW           16384
#define CHLFAIL         32768
#define NAVWARN         65536
#define ABSAER         131072
#define SPARE3         262144
#define MAXAERITER     524288
#define MODGLINT      1048576
#define CHLWARN       2097152
#define ATMWARN       4194304
#define OPSHAL        8388608
#define SEAICE       16777216
#define NAVFAIL      33554432
#define FILTER       67108864
#define SPARE5      134217728
#define BOWTIEDEL   268435456
#define HIPOL       536870912
#define PRODFAIL   1073741824
#define SPARE7     2147483648

// assign the georegion flag to one of the spare flags
#define GEOREGION SPARE1


static const char * const l2_flag_lname[L1_NFLAGS] = {"ATMFAIL",
    "LAND",
    "PRODWARN",
    "HIGLINT",
    "HILT",
    "HISATZEN",
    "COASTZ",
    "SPARE",
    "STRAYLIGHT",
    "CLDICE",
    "COCCOLITH",
    "TURBIDW",
    "HISOLZEN",
    "SPARE",
    "LOWLW",
    "CHLFAIL",
    "NAVWARN",
    "ABSAER",
    "SPARE",
    "MAXAERITER",
    "MODGLINT",
    "CHLWARN",
    "ATMWARN",
    "OPSHAL",
    "SEAICE",
    "NAVFAIL",
    "FILTER",
    "SPARE",
    "BOWTIEDEL",
    "HIPOL",
    "PRODFAIL",
    "SPARE"};

static const char * const l2_flag_sname[L1_NFLAGS] = {"f01_name",
    "f02_name",
    "f03_name",
    "f04_name",
    "f05_name",
    "f06_name",
    "f07_name",
    "f08_name",
    "f09_name",
    "f10_name",
    "f11_name",
    "f12_name",
    "f13_name",
    "f14_name",
    "f15_name",
    "f16_name",
    "f17_name",
    "f18_name",
    "f19_name",
    "f20_name",
    "f21_name",
    "f22_name",
    "f23_name",
    "f24_name",
    "f25_name",
    "f26_name",
    "f27_name",
    "f28_name",
    "f29_name",
    "f30_name",
    "f31_name",
    "f32_name"};


#ifdef __cplusplus
}
#endif

#endif
