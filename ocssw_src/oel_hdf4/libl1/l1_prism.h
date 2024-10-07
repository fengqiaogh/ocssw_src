#ifndef  _L1_PRISM_H
#define  _L1_PRISM_H

#include <stdint.h>
#include "l1.h"

int closel1_prism(filehandle *l1file);
int openl1_prism(filehandle *l1file);
int readl1_prism(filehandle *l1file, int32_t recnum, l1str *l1rec);

#endif
