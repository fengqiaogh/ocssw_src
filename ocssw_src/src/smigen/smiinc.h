#ifndef SMIINC_H
#define SMIINC_H

#define L3M_SOFTNM_VAL  "smigen"
#define L3M_SOFTVER_VAL "5.20"

#define SMI_MIN(a,b)   (a) > (b) ? (b) : (a)
#define SMI_MAX(a,b)   (a) > (b) ? (a) : (b)
#define L3M_MEASURES     7
#define SMI_NX          4096
#define SMI_NY          2048
#define SMI_BUFSIZE     250000
#define NUM_SMI_ARGS    6
#define SMI_LAT_NORTH   90.
#define SMI_LAT_SOUTH   -90.
#define SMI_LON_WEST    -180.
#define SMI_LON_EAST    180.
#define SMI_MAX_STR     1024
#define SMI_MAX_INFILES  30
#define SMI_MAX_STR_SHORT 50
#endif


