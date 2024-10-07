#ifndef L1_L1C_ANC_H
#define L1_L1C_ANC_H

#include <stdint.h>
#include "l1.h"

#ifdef __cplusplus
extern "C" {
#endif

// Open, read, close
int openl1_l1c_anc(filehandle *file);
int readl1_l1c_anc(filehandle *file, int32_t recnum, l1str *l1rec);
int closel1_l1c_anc(filehandle *file);

#ifdef __cplusplus
}
#endif
    
#endif