#ifndef  _L1_MODIS_H
#define  _L1_MODIS_H

#include <stdint.h>
#include "l1.h"

int closel1_modis();
int openl1_modis(filehandle *l1file);
int readl1_modis(filehandle *l1file, int32_t recnum, l1str *l1rec, int lonlat);
int readl1_lonlat_modis(filehandle *l1file, int32_t recnum, l1str *l1rec);

#endif
