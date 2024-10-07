#ifndef  _L1_HDF_GENERIC_READ_H
#define  _L1_HDF_GENERIC_READ_H

#include <stdint.h>
#include "l1.h"

int closel1_hdf_g(filehandle *l1file);
int openl1_read_hdf_g(filehandle *l1file);
int readl1_hdf_g(filehandle *l1file, int32_t recnum, l1str *l1rec);

#endif
