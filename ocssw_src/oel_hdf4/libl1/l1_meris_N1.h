#ifndef  _L1_MERIS_N1_H
#define  _L1_MERIS_N1_H

#include <stdint.h>
#include "l1.h"

int closel1_meris_N1(filehandle *l1file);
int openl1_meris_N1(filehandle *l1file);
int readl1_meris_N1(filehandle *l1file, int32_t recnum, l1str *l1rec);
int readl1_lonlat_meris_N1(filehandle *l1file, int32_t recnum, l1str *l1rec);

#endif
