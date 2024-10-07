#ifndef  _L1_MERIS_CC_H
#define  _L1_MERIS_CC_H

#include <stdint.h>
#include "l1.h"

int closel1_meris_CC(filehandle *l1file);
int openl1_meris_CC(filehandle *l1file);
int readl1_meris_CC(filehandle *l1file, int32_t recnum, l1str *l1rec);

#endif
