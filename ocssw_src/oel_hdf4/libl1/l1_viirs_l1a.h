#ifndef  _L1_VIIRS_L1A_H
#define  _L1_VIIRS_L1A_H

#include <stdint.h>
#include "l1.h"

int closel1_viirs_l1a(filehandle *l1file);
int openl1_viirs_l1a(filehandle *l1file);
int readl1_viirs_l1a(filehandle *l1file, int32_t recnum, l1str *l1rec);

#endif
