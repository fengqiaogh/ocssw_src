#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"

static float badval = BAD_FLT;
static float minval = 0.0;
static float maxval = 10.0 * 1000;

static void init_poc_bands(const int sensorId, const int band, int *bandIndex1, int *bandIndex2) {
    *bandIndex1 = bindex_get(band);

    if ((sensorId == MODISA || sensorId == MODIST || sensorId == HICO) && bandIndex2 >= 0) {
        *bandIndex2 = bindex_get(547);
        return; // Don't need to check other bands
    } else {
        *bandIndex2 = bindex_get(545);
    }

    if (*bandIndex2 < 0) 
        *bandIndex2 = bindex_get(550);
    if (*bandIndex2 < 0) 
        *bandIndex2 = bindex_get(555);
    if (*bandIndex2 < 0) 
        *bandIndex2 = bindex_get(560);

    if (bandIndex1 < 0 || bandIndex2 < 0) {
        printf("-E- %s line %d: required bands not available for Stramski POC\n",
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
}

float poc_stramski_443(l2str* l2rec, int currPixel) {
    static int firstCall = 1;
    static float a = 203.2;
    static float b = -1.034;
    static int ib1 = -1;
    static int ib2 = -1;
    float *Rrs = &l2rec->Rrs[currPixel * l2rec->l1rec->l1file->nbands];
    float *wave = l2rec->l1rec->l1file->fwave;

    float poc = badval;
    float Rrs1 = 0.0;
    float Rrs2 = 0.0;

    if (firstCall) {
        firstCall = 0;
        init_poc_bands(l2rec->l1rec->l1file->sensorID, 443, &ib1, &ib2);
    }

    Rrs1 = Rrs[ib1];
    Rrs2 = Rrs[ib2];

    if (Rrs1 >= 0.0 && Rrs2 > 0.0) {
        Rrs2 = conv_rrs_to_555(Rrs2, wave[ib2],-99,NULL);
        poc = a * pow(Rrs1 / Rrs2, b);
    }

    return (poc);
}

/*------------------------------------------------*/
/*Standard uncertainty for poc_stramksi_443       */
/*esimated using analytical propagation           */
/*------------------------------------------------*/

float unc_poc_stramski_443(l2str *l2rec, int ipb) {

    l1str *l1rec = l2rec->l1rec;
    uncertainty_t *uncertainty=l1rec->uncertainty;
    float *wave = l2rec->l1rec->l1file->fwave;

    static int firstCall = 1;
    static float a = 203.2;
    static float b = -1.034;
    static int ib1 = -1;
    static int ib2 = -1;

    float upoc = badval;

    float Rrs1 = 0.0;
    float Rrs2 = 0.0;
    float uRrs1 = 0.0;
    float uRrs2 = 0.0;
    float *uRrs2_new;
    uRrs2_new = calloc(1, sizeof(float)); 
    float dRatdRrs1;
    float dRatdRrs2;
    float dPocdRat;
    float dPocdRrs1;
    float dPocdRrs2;
    float rat;

    if (firstCall) {
        firstCall = 0;
        init_poc_bands(l2rec->l1rec->l1file->sensorID, 443, &ib1, &ib2);
    }

    if (uncertainty) {
        Rrs1 = l2rec->Rrs[ipb + ib1];
        Rrs2 = l2rec->Rrs[ipb + ib2];
        uRrs1 = l2rec->Rrs_unc[ipb + ib1];
        uRrs2 = l2rec->Rrs_unc[ipb + ib2];
    }

    
    if (Rrs1 >= 0.0 && Rrs2 > 0.0) {
        Rrs2 = conv_rrs_to_555(Rrs2, wave[ib2],uRrs2,uRrs2_new);
        rat = Rrs1/Rrs2;
        dPocdRat = a*b*pow(rat,b-1.);
        dRatdRrs1 = 1./Rrs2;
        dRatdRrs2 = -Rrs1/pow(Rrs2,2.);
        dPocdRrs1 = dPocdRat * dRatdRrs1;
        dPocdRrs2 = dPocdRat * dRatdRrs2;

        upoc = sqrt ( pow(dPocdRrs1*uRrs1,2) + pow(dPocdRrs2*(*uRrs2_new),2));

    }
    
    free(uRrs2_new);
    return (upoc);
}



float poc_stramski_490(l2str* l2rec, int currPixel) {
    static int firstCall = 1;
    static float a = 308.3;
    static float b = -1.639;
    static int ib1 = -1;
    static int ib2 = -1;
    float *Rrs = &l2rec->Rrs[currPixel * l2rec->l1rec->l1file->nbands];
    float *wave = l2rec->l1rec->l1file->fwave;

    float poc = badval;
    float Rrs1 = 0.0;
    float Rrs2 = 0.0;

    if (firstCall) {
        firstCall = 0;
        init_poc_bands(l2rec->l1rec->l1file->sensorID, 490, &ib1, &ib2);
    }

    Rrs1 = Rrs[ib1];
    Rrs2 = Rrs[ib2];

    if (Rrs1 >= 0.0 && Rrs2 > 0.0) {
        Rrs2 = conv_rrs_to_555(Rrs2, wave[ib2],-99,NULL);
        poc = a * pow(Rrs1 / Rrs2, b);
    }

    return (poc);
}

/*------------------------------------------------*/
/*Standard uncertainty for poc_stramksi_490       */
/*esimated using analytical propagation           */
/*------------------------------------------------*/

float unc_poc_stramski_490(l2str *l2rec, int ipb) {

    l1str *l1rec = l2rec->l1rec;
    uncertainty_t *uncertainty=l1rec->uncertainty;
    float *wave = l2rec->l1rec->l1file->fwave;

    static int firstCall = 1;
    static float a = 308.3;
    static float b = -1.639;
    static int ib1 = -1;
    static int ib2 = -1;

    float upoc = badval;
    float Rrs1 = 0.0;
    float Rrs2 = 0.0;
    float uRrs1 = 0.0;
    float uRrs2 = 0.0;
    float *uRrs2_new;
    uRrs2_new = (float*) calloc(1, sizeof(float)); 
    float dRatdRrs1;
    float dRatdRrs2;
    float dPocdRat;
    float dPocdRrs1;
    float dPocdRrs2;
    float rat;


    if (firstCall) {
        firstCall = 0;
        init_poc_bands(l2rec->l1rec->l1file->sensorID, 490, &ib1, &ib2);
    }

    if (uncertainty) {
        Rrs1 = l2rec->Rrs[ipb + ib1];
        Rrs2 = l2rec->Rrs[ipb + ib2];
        uRrs1 = l2rec->Rrs_unc[ipb + ib1];
        uRrs2 = l2rec->Rrs_unc[ipb + ib2];
    }

    if (Rrs1 >= 0.0 && Rrs2 > 0.0) {
        Rrs2 = conv_rrs_to_555(Rrs2, wave[ib2],uRrs2,uRrs2_new);

        rat = Rrs1/Rrs2;
        dPocdRat = a*b*pow(rat,b-1.);
        dRatdRrs1 = 1./Rrs2;
        dRatdRrs2 = -Rrs1/pow(Rrs2,2.);
        dPocdRrs1 = dPocdRat * dRatdRrs1;
        dPocdRrs2 = dPocdRat * dRatdRrs2;

        upoc = sqrt ( pow(dPocdRrs1*uRrs1,2) + pow(dPocdRrs2*(*uRrs2_new),2));
    }

    free(uRrs2_new);
    return (upoc);
}


void get_poc(l2str *l2rec, l2prodstr *p, float prod[]) {
    int32_t ip, ipb;
    l1str *l1rec = l2rec->l1rec;
    int32_t nbands = l1rec->l1file->nbands;

    for (ip = 0; ip < l1rec->npix; ip++) {

        prod[ip] = badval;
        ipb = nbands * ip;

        switch (p->cat_ix) {
            case CAT_poc_stramski_443:
                prod[ip] = poc_stramski_443(l2rec, ip);
                break;

            case CAT_poc_stramski_490:
                prod[ip] = poc_stramski_490(l2rec, ip);
                break;

            case CAT_poc_unc_stramski_443:
                prod[ip] = unc_poc_stramski_443(l2rec, ipb);
                break;

            case CAT_poc_unc_stramski_490:
                prod[ip] = unc_poc_stramski_490(l2rec, ipb);
                break;
            case CAT_poc_stramski_hybrid:
                prod[ip] = poc_stramski_hybrid(&l2rec->Rrs[ip * nbands], l1rec->l1file->sensorID);
                if (prod[ip] < p->min || prod[ip] > p->max)
                    prod[ip] = BAD_FLT;
                break;
            default:
                printf("Error: %s : Unknown product specifier: %d\n", __FILE__, p->cat_ix);
                exit(1);
                break;
        };

        if (prod[ip] == badval)
            l1rec->flags[ip] |= PRODFAIL;
        else if (prod[ip] < minval || prod[ip] > maxval)
            l1rec->flags[ip] |= PRODWARN;

    }

}

