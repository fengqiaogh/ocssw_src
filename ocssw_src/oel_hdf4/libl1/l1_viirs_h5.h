#ifndef  _L1_VIIRS_H5_H
#define  _L1_VIIRS_H5_H

#include <stdint.h>
#include "l1.h"
#include "h5io.h"

int closel1_viirs_h5(filehandle *l1file);
int openl1_viirs_h5(filehandle *l1file);
int readl1_viirs_h5(filehandle *l1file, int32_t recnum, l1str *l1rec, int lonlat);

#endif
