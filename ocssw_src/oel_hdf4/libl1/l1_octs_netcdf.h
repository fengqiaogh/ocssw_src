#ifndef  _L1_OCTS_NETCDF_H
#define  _L1_OCTS_NETCDF_H

#include <stdint.h>
#include "filehandle.h"
#include "l1.h"

#ifdef __cplusplus
extern "C" {
#endif

int32_t openl1_octs_netcdf(filehandle *l1file);
int32_t readl1_octs_netcdf(filehandle *l1file, int32_t recnum, l1str *l1rec);
int32_t closel1_octs_netcdf(filehandle *l1file);

#ifdef __cplusplus
}
#endif

#endif

