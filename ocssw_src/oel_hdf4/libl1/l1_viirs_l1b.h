#ifndef  _L1_VIIRS_L1B_H
#define  _L1_VIIRS_L1B_H

#include <stdint.h>
#include "l1.h"

int closel1_viirs_l1b();
int openl1_viirs_l1b(filehandle *l1file);
int readl1_viirs_l1b(filehandle *l1file, int32_t recnum, l1str *l1rec, int lonlat);

#endif
