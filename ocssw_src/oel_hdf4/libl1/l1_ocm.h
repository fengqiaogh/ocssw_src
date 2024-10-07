#ifndef  _L1_OCM_H
#define  _L1_OCM_H

#include <stdint.h>
#include "l1.h"

int closel1_ocm(filehandle *l1file);
int openl1_ocm(filehandle *l1file);
int readl1_ocm(filehandle *l1file, int32_t recnum, l1str *l1rec);

#endif
