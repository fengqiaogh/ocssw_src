/*
 * l1_safe.h
 */

#ifndef L1_SAFE_H_
#define L1_SAFE_H_

#include <stdint.h>
#include "l1.h"

#ifdef __cplusplus
extern "C" {
#endif

int closel1_safe(filehandle *file);
int openl1_safe(filehandle *file);
int readl1_safe(filehandle *file, int32_t recnum, l1str *l1rec, int lonlat);

#ifdef __cplusplus
}
#endif

#endif /* L1_SAFE_H_ */

