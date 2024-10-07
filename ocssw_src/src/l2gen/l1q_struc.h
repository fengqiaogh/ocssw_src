#ifndef _L1Q_STRUC
#define _L1Q_STRUC

#include "l1.h"

typedef struct l1q_struct {
    int32_t nq;
    int32_t cscan;
    l1str r[NQMAX];
} l1qstr;

#endif
