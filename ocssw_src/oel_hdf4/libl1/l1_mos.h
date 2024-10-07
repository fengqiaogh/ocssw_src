#ifndef  _L1_MOS_H
#define  _L1_MOS_H

#include <stdint.h>
#include "l1.h"

int openl1_read_mos(filehandle *l1file);
int readl1_mos(filehandle *l1file, int32_t recnum, l1str *l1rec);
int closel1_mos(filehandle *l1file);

#endif
