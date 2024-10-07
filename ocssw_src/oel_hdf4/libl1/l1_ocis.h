#ifndef L1_OCIS_H
#define L1_OCIS_H

#ifdef __cplusplus
extern "C" {
#endif
    
#include <stdint.h>
#include "l1.h"

// Open, read, close
int openl1_ocis(filehandle *file);
int readl1_ocis(filehandle *file, int32_t recnum, l1str *l1rec);
int closel1_ocis(filehandle *file);

    
#ifdef __cplusplus
}
#endif
    
#endif