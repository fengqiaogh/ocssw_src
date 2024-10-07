#ifndef L1_SEAWIFS_NETCDF_H
#define L1_SEAWIFS_NETCDF_H

#include <stdint.h>
#include "filehandle.h"
#include "l1.h"

#ifdef __cplusplus
extern "C" {
#endif

int openl1_seawifs_netcdf(filehandle *file);
int readl1_seawifs_netcdf(filehandle *file, int32_t recnum, l1str *l1rec);
int readl1_lonlat_seawifs_netcdf(filehandle *file, int32_t recnum, l1str *l1rec);
int closel1_seawifs_netcdf(filehandle *file);

#ifdef __cplusplus
}
#endif


#endif
