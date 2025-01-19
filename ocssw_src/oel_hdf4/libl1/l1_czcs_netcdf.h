/*
 *  W. Robinson, SAIC, 10 Dec 2004  new for CZCS
 */
#ifndef L1_CZCS_NETCDF_H
#define L1_CZCS_NETCDF_H

#include <stdint.h>
#include "filehandle.h"
#include "l1.h"

#ifdef __cplusplus
extern "C" {
#endif

int32_t openl1_czcs_netcdf(filehandle *file);
int32_t readl1_czcs_netcdf(filehandle *file, int32_t recnum, l1str *l1rec);
int32_t closel1_czcs_netcdf(filehandle *file);

#ifdef __cplusplus
}
#endif

#endif
