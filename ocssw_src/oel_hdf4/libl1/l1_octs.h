#ifndef  _L1_OCTS_H
#define  _L1_OCTS_H

#include <stdint.h>
#include "l1.h"


int openl1_octs(filehandle *l1file);
int readl1_octs(filehandle *l1file, int32_t recnum, l1str *l1rec);
int closel1_octs(filehandle *l1file);

#endif

