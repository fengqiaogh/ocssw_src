#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"

static int32_t npix = -1;
static int32_t nbands = -1;
static vcstr vrec;

int alloc_vcal(int32_t npix, int32_t nbands, vcstr *rec);

void vgain(l2str *l2rec, int band, float vgain[]) {
    int32_t ip, ipb;

    for (ip = 0; ip < npix; ip++) {
        ipb = ip * nbands + band;
        if (vrec.vLt[ipb] > 0.0 && l2rec->l1rec->Lt[ipb] > 0.0) {
            vgain[ip] = vrec.vLt[ipb] / (l2rec->l1rec->Lt[ipb] / l1_input->gain[band]);
        } else
            vgain[ip] = -1.0;
    }

    return;
}

void vcal(l2str *l2rec, l2prodstr *p, float prod[]) {
    static int firstCall = 1;
    static int32_t lastScan = -1;

    int32_t ip, ipb;
    int32_t band = p->prod_ix;

    if (firstCall) {
        npix = l2rec->l1rec->npix;
        nbands = l2rec->l1rec->l1file->nbands;
        if (alloc_vcal(npix, nbands, &vrec) == 0) {
            printf("-E- %s: Unable to allocate vcal record.\n", __FILE__);
            exit(1);
        }
        firstCall = 0;
    }

    if (l2rec->l1rec->iscan != lastScan) {
        convl21(l2rec, l2rec->tgrec, 0, npix - 1, vrec.vLt, &vrec);
        lastScan = l2rec->l1rec->iscan;
    }

    switch (p->cat_ix) {
    case CAT_vgain:
        vgain(l2rec, band, prod);
        break;
    case CAT_vLt:
        for (ip = 0; ip < npix; ip++) {
            ipb = ip * nbands + band;
            prod[ip] = vrec.vLt[ipb];
        }
        break;
    case CAT_vtLw:
        for (ip = 0; ip < npix; ip++) {
            ipb = ip * nbands + band;
            prod[ip] = vrec.tLw[ipb];
        }
        break;
    case CAT_vLw:
        for (ip = 0; ip < npix; ip++) {
            ipb = ip * nbands + band;
            prod[ip] = vrec.Lw[ipb];
        }
        break;
    case CAT_vnLw:
        for (ip = 0; ip < npix; ip++) {
            ipb = ip * nbands + band;
            prod[ip] = vrec.nLw[ipb];
        }
        break;
    case CAT_vbsat:
        for (ip = 0; ip < npix; ip++) {
            ipb = ip * nbands + band;
            prod[ip] = vrec.brdfsat[ipb];
        }
        break;
    case CAT_vbtgt:
        for (ip = 0; ip < npix; ip++) {
            ipb = ip * nbands + band;
            prod[ip] = vrec.brdftgt[ipb];
        }
        break;
    }

    return;
}
