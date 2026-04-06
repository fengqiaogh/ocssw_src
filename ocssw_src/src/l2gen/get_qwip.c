#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"

const static double p[] = {-8.399885e-9, 1.715532e-5, -1.301670e-2, 4.357838e0, -5.449532e2};

static float* avw = NULL;

void get_qwip(l2str* l2rec, float prod[]) {
    int32_t nbands = l2rec->l1rec->l1file->nbands;
    int32_t npix = l2rec->l1rec->npix;
    int ib1 = bindex_get(490);
    int ib2 = bindex_get(665);
    static int first = 1;
    if (first) {
        if (ib1 < 0 || ib2 < 0) {
            printf("QWIP: incompatible sensor wavelengths for this algorithm\n");
            exit(1);
        }
        avw = calloc(npix, sizeof(float));
        first = 0;
    }
    float* Rrs = l2rec->Rrs;
    get_avw(l2rec, avw);
    for (int32_t ip = 0; ip < npix; ip++) {
        float Rrs1 = Rrs[ib1 + ip * nbands];
        float Rrs2 = Rrs[ib2 + ip * nbands];
        if (Rrs1 == BAD_FLT || Rrs2 == BAD_FLT) {
            prod[ip] = BAD_FLT;
            continue;
        }
        float ndi = (Rrs2 - Rrs1) / (Rrs2 + Rrs1);

        if (avw[ip] == BAD_FLT) {
            prod[ip] = BAD_FLT;
            continue;
        }
        float power = 1;
        float ndi_predicted = 0;
        for (int ipower = 0; ipower < 5; ipower++) {
            ndi_predicted += power * p[4 - ipower];
            power *= avw[ip];
        }
        prod[ip] = ndi - ndi_predicted;
    }
};
