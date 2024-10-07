#ifndef L1_OCI_H
#define L1_OCI_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include "l1.h"


// Open, read, close
int openl1_oci(filehandle *file);
int readl1_oci(filehandle *file, int32_t recnum, l1str *l1rec, int lonlat);
int closel1_oci(filehandle *file);
    
#ifdef __cplusplus
}
#endif
    
#endif