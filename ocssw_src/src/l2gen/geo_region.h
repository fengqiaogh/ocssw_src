#ifndef GEO_REGION_H
#define GEO_REGION_H
#include "l2_flags.h"
#define GEOREGION SPARE1
#ifdef __cplusplus
extern "C" {
#endif
int get_georegion(float lat, float lon);
void close_georegion_file();
#ifdef __cplusplus
}
#endif

#endif