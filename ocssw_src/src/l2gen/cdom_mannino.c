#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"

static float badval = BAD_FLT;

void cdom_mannino(l2str *l2rec, int prodnum, float prod[]) {
    static int firstCall = 1;
    static int ib1 = -1;
    static int ib2 = -1;

    static float b_ag412[] = {-2.784, -1.146, 1.008};
    static float b_sg275[] = {-3.325, 0.3, -0.252};
    static float b_sg300[] = {-3.679, 0.168, -0.134};

    l1str *l1rec = l2rec->l1rec;
    int32_t nbands = l1rec->l1file->nbands;

    float *wave = l1rec->l1file->fwave;
    float *b;
    float *Rrs, Rrs1, Rrs2;
    float x1, x2;
    int32_t ip;

    if (firstCall) {
        firstCall = 0;
        ib1 = bindex_get(443);
        ib2 = bindex_get(545);
        if (ib2 < 0) ib2 = bindex_get(550);
        if (ib2 < 0) ib2 = bindex_get(555);
        if (ib2 < 0) ib2 = bindex_get(560);
        if (ib2 < 0) ib2 = bindex_get(565);
        if (ib1 < 0 || ib2 < 0) {
            printf("-E- %s line %d: required bands not available for Stramski POC\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    }

    switch (prodnum) {
    case CAT_ag_412_mlrc:
        b = b_ag412;
        break;
    case CAT_Sg_275_295_mlrc:
        b = b_sg275;
        break;
    case CAT_Sg_300_600_mlrc:
        b = b_sg300;
        break;
    default:
        printf("Error: %s : Unknown product specifier: %d\n", __FILE__, prodnum);
        exit(1);
        break;
    }

    for (ip = 0; ip < l1rec->npix; ip++) {

        prod[ip] = badval;

        Rrs = &l2rec->Rrs[ip * nbands];
        Rrs1 = Rrs[ib1];
        Rrs2 = Rrs[ib2];

        if (Rrs1 > 0.0 && Rrs2 > 0.0) {

            Rrs2 = conv_rrs_to_555(Rrs2, wave[ib2], -99,NULL);

            x1 = log(Rrs1);
            x2 = log(Rrs2);

            prod[ip] = exp(b[0] + b[1] * x1 + b[2] * x2);

        } else {
            l1rec->flags[ip] |= PRODFAIL;
        }

    }

    return;
}


