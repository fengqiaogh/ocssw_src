#ifndef  _L1_VIIRS_L1B_H
#define  _L1_VIIRS_L1B_H

#include <stdint.h>
#include "l1.h"
#ifdef __cplusplus
extern "C" {
#endif

int closel1_viirs_l1b();
/**
 *
 * @param l1file input l1file struct
 * @return success flag
 */
int openl1_viirs_l1b(filehandle *l1file);
int readl1_viirs_l1b(filehandle *l1file, int32_t recnum, l1str *l1rec, int lonlat);

#ifdef __cplusplus
    }
#endif

#endif
