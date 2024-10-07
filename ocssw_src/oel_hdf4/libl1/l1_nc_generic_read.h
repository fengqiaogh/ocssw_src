#ifndef  _L1_NC_GENERIC_READ_H
#define  _L1_NC_GENERIC_READ_H


#include <stdint.h>
#include "l1.h"

int closel1_nc_generic(filehandle *l1file);
int openl1_nc_generic(filehandle *l1file);
int readl1_nc_generic(filehandle *l1file, int32_t recnum, l1str *l1rec);

#endif
