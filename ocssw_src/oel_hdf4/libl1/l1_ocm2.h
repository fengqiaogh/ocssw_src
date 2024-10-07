#ifndef  _L1_OCM2_H
#define  _L1_OCM2_H

#include <stdint.h>
#include "l1.h"

int closel1_ocm2(filehandle *l1file);
int openl1_ocm2(filehandle *l1file);
int readl1_ocm2(filehandle *l1file, int32_t recnum, l1str *l1rec);

#endif
