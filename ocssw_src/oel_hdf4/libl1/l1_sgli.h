/*
 * l1_sgli.h
 */

#ifndef L1_SGLI_H_
#define L1_SGLI_H_

#include <stdint.h>
#include "l1.h"

int closel1_sgli(filehandle *file);
int openl1_sgli(filehandle *file);
int readl1_sgli(filehandle *file, int32_t recnum, l1str *l1rec);

#endif /* L1_SGLI_H_ */
