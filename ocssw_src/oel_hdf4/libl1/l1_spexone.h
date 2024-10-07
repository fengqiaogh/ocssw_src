#ifndef L1_SPEXONE_H
#define L1_SPEXONE_H

#include <stdint.h>
#include "l1.h"

#ifdef __cplusplus
extern "C" {
#endif

// Open, read, close
int openl1_spexone(filehandle *file);
int readl1_spexone(filehandle *file, int32_t recnum, l1str *l1rec);
int closel1_spexone(filehandle *file);

#ifdef __cplusplus
}
#endif
    
#endif