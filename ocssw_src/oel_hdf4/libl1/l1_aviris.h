#ifndef  _L1_AVIRIS_H
#define  _L1_AVIRIS_H

#include <stdint.h>
#include "l1.h"
#include "l1_aviris_struc.h"

#ifdef __cplusplus
extern "C" {
#endif

aviris_t* createPrivateData_aviris(int numBands);
void freePrivateData_aviris(aviris_t* data);

int closel1_aviris(filehandle *l1file);
int openl1_aviris(filehandle *l1file);
int readl1_aviris(filehandle *l1file, int32_t recnum, l1str *l1rec);
int closel1_aviris_nc(filehandle *l1file);
int openl1_aviris_nc(filehandle *l1file);
int readl1_aviris_nc(filehandle *l1file, int32_t recnum, l1str *l1rec);

#ifdef __cplusplus
}
#endif

#endif
