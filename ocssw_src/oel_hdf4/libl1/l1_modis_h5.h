#ifndef _L1_MODIS_H5_H
#define _L1_MODIS_H5_H

#include <stdint.h>
#include "l1.h"

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus

int closel1_modis_h5();
int openl1_modis_h5(filehandle *l1file);
int readl1_modis_h5(filehandle *l1file, int32_t recnum, l1str *l1rec, int lonlat);
int readl1_lonlat_modis_h5(filehandle *l1file, int32_t recnum, l1str *l1rec);

#ifdef __cplusplus
}
#endif  // __cplusplus

#endif  // _L1_MODIS_H5_H
