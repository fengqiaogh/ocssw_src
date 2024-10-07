#ifndef L1_SEAWIFS_H
#define L1_SEAWIFS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

int openl1_seawifs(filehandle *file);
int readl1_seawifs(filehandle *file, int32_t recnum, l1str *l1rec);
int readl1_lonlat_seawifs(filehandle *file, int32_t recnum, l1str *l1rec);
int closel1_seawifs(filehandle *file);

#ifdef __cplusplus
}
#endif

#endif
