#include "l12_proto.h"

static float *nw = NULL;
static float *aw = NULL;
static float *bbw = NULL;
static int32_t nbands;

void seawater_set(l1str *l1rec) {
    nw = l1rec->sw_n;
    aw = l1rec->sw_a;
    bbw = l1rec->sw_bb;
    nbands = l1rec->l1file->nbands;

    return;
}

float seawater_get_n(int32_t ip, int32_t ib) {
    return (nw[ip * nbands + ib]);
}

float seawater_get_a(int32_t ip, int32_t ib) {
    return (aw[ip * nbands + ib]);
}

float seawater_get_bb(int32_t ip, int32_t ib) {
    return (bbw[ip * nbands + ib]);
}

