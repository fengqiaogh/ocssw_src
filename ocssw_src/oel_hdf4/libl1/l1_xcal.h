#ifndef  _L1_XCAL_H
#define  _L1_XCAL_H

#include <stdint.h>
#include "l1.h"

int closel1_xcal(filehandle *l1file);
int openl1_xcal(filehandle *l1file);
int readl1_xcal(filehandle *l1file, int32_t recnum, l1str *l1rec);

#endif
