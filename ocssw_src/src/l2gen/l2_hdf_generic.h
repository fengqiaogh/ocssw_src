#ifndef  _L2_HDF_GENERIC_H
#define  _L2_HDF_GENERIC_H

#include <stdint.h>
#include "l2_struc.h"
#include "filehandle.h"

int closel2_hdf(filehandle *l2file);
int openl2_hdf(filehandle *l2file);
int writel2_hdf(filehandle *l2file, int32_t recnum, l2str *l2rec);

#endif
