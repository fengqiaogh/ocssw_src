#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"
#include "chl.h"


float get_chl_ocx(l2str *l2rec, float Rrs[]);

float chl_oc2(l2str *l2rec, float Rrs[]) {
    static int32_t *w = NULL;
    static float *a = NULL;
    static int ib1 = -1;
    static int ib2 = -1;

    float rat;
    float Rrs1, Rrs2;
    float chl = chlbad;

    if (w == NULL) {
        w = input->chloc2w;
        a = input->chloc2c;
        if (w[0] < 0 || w[1] < 0) {
            printf("chl_oc2: algorithm coefficients not provided for this sensor.\n");
            exit(1);
        }
        ib1 = bindex_get(w[0]);
        ib2 = bindex_get(w[1]);
        if (ib1 < 0 || ib2 < 0) {
            printf("chl_oc2: incompatible sensor wavelengths for this algorithm\n");
            exit(1);
        }
    }

    Rrs1 = Rrs[ib1];
    Rrs2 = Rrs[ib2];

    if (Rrs1 > 0.0 && Rrs2 > 0.0 && Rrs1 / Rrs2 <= 10.0) {
        rat = Rrs1 / Rrs2;
        if (rat > minrat && rat < maxrat) {
            rat = log10(rat);
            chl = (float)
                    pow(10.0, (a[0] + rat * (a[1] + rat * (a[2] + rat * (a[3] + rat * a[4])))));
            chl = (chl > chlmin ? chl : chlmin);
            chl = (chl < chlmax ? chl : chlmax);
        }
    }

    return (chl);
}

float chl_oc3(l2str *l2rec, float Rrs[]) {
    static int32_t *w = NULL;
    static float *a = NULL;
    static int ib1 = -1;
    static int ib2 = -1;
    static int ib3 = -1;

    float rat, minRrs;
    float Rrs1, Rrs2, Rrs3;
    float chl = chlbad;
    float *dRrs, *dchl;
    uncertainty_t *uncertainty=l2rec->l1rec->uncertainty;

    if(uncertainty){
        dRrs=uncertainty->dRrs;
        dchl=&uncertainty->dchl;
    }

    if (w == NULL) {
        w = input->chloc3w;
        a = input->chloc3c;
        if (w[0] < 0 || w[1] < 0 || w[2] < 0) {
            printf("chl_oc3: algorithm coefficients not provided for this sensor.\n");
            exit(1);
        }
        ib1 = bindex_get(w[0]);
        ib2 = bindex_get(w[1]);
        ib3 = bindex_get(w[2]);
        if (ib1 < 0 || ib2 < 0 || ib3 < 0) {
            printf("chl_oc3: incompatible sensor wavelengths for this algorithm\n");
            exit(1);
        }
    }

    Rrs1 = Rrs[ib1];
    Rrs2 = Rrs[ib2];
    Rrs3 = Rrs[ib3];

    minRrs = MIN(Rrs1, Rrs2);

    if (Rrs3 > 0.0 && Rrs2 > 0.0 && minRrs > -0.001) {
        rat = MAX(Rrs1, Rrs2) / Rrs3;
        if (rat > minrat && rat < maxrat) {
            rat = log10(rat);
            chl = (float)
                    pow(10.0, (a[0] + rat * (a[1] + rat * (a[2] + rat * (a[3] + rat * a[4])))));
  
            if(uncertainty){
                if(Rrs1>Rrs2){
                    *dchl= a[1]+2*a[2]*rat+3*a[3]*rat*rat+4*a[4]*rat*rat*rat;
                    *dchl=chl*chl*pow(*dchl,2)*pow(dRrs[ib1]/Rrs1-dRrs[ib3]/Rrs[ib3],2);
                    *dchl=sqrt(*dchl);
                }
                else{
                    *dchl= a[1]+2*a[2]*rat+3*a[3]*rat*rat+4*a[4]*rat*rat*rat;
                    *dchl=chl*chl*pow(*dchl,2)*pow(dRrs[ib2]/Rrs2-dRrs[ib3]/Rrs[ib3],2);
                    *dchl=sqrt(*dchl);
                }
            }

            chl = (chl > chlmin ? chl : chlmin);
            chl = (chl < chlmax ? chl : chlmax);

            if(uncertainty){
                if(fabs(chl-chlmin)<0.0000001 || fabs(chl-chlmax)<0.0000001)
                    *dchl=0;
            }
        }
    }

    return (chl);
}

float chl_oc3c(l2str *l2rec, float Rrs[]) {
    static float a[] = {0.2515, -2.3798, 1.5823, -0.6372, -0.5692};
    static int ib1 = -1;
    static int ib2 = -1;
    static int ib3 = -1;

    float rat, minRrs;
    float Rrs1, Rrs2, Rrs3;
    float chl = chlbad;

    if (ib1 < 0) {
        ib1 = bindex_get(443);
        ib2 = bindex_get(490);
        ib3 = bindex_get_555();

        if (ib1 < 0 || ib2 < 0 || ib3 < 0) {
            printf("chl_oc3: incompatible sensor wavelengths for this algorithm\n");
            exit(1);
        }
    }

    Rrs1 = Rrs[ib1];
    Rrs2 = Rrs[ib2];
    Rrs3 = Rrs[ib3];

    minRrs = MIN(Rrs1, Rrs2);

    if (Rrs3 > 0.0 && Rrs2 > 0.0 && minRrs > -0.001) {
        Rrs3 = conv_rrs_to_555(Rrs3, l2rec->l1rec->l1file->fwave[ib3],-99,NULL);
        rat = MAX(Rrs1, Rrs2) / Rrs3;
        if (rat > minrat && rat < maxrat) {
            rat = log10(rat);
            chl = (float)
                    pow(10.0, (a[0] + rat * (a[1] + rat * (a[2] + rat * (a[3] + rat * a[4])))));
            chl = (chl > chlmin ? chl : chlmin);
            chl = (chl < chlmax ? chl : chlmax);
        }
    }

    return (chl);
}

float chl_oc4(l2str *l2rec, float Rrs[]) {
    static int32_t *w = NULL;
    static float *a = NULL;
    static int ib1 = -1;
    static int ib2 = -1;
    static int ib3 = -1;
    static int ib4 = -1;

    float rat, minRrs;
    float Rrs1, Rrs2, Rrs3, Rrs4;
    float chl = chlbad;

    if (w == NULL) {
        w = input->chloc4w;
        a = input->chloc4c;
        if (w[0] < 0 || w[1] < 0 || w[2] < 0 || w[3] < 0) {
            printf("chl_oc4: algorithm coefficients not provided for this sensor.\n");
            exit(1);
        }
        ib1 = bindex_get(w[0]);
        ib2 = bindex_get(w[1]);
        ib3 = bindex_get(w[2]);
        ib4 = bindex_get(w[3]);
        if (ib1 < 0 || ib2 < 0 || ib3 < 0 || ib4 < 0) {
            printf("chl_oc4: incompatible sensor wavelengths for this algorithm\n");
            exit(1);
        }
    }

    Rrs1 = Rrs[ib1];
    Rrs2 = Rrs[ib2];
    Rrs3 = Rrs[ib3];
    Rrs4 = Rrs[ib4];

    minRrs = MIN(Rrs1, Rrs2);

    if (Rrs4 > 0.0 && Rrs3 > 0.0 && (Rrs2 > 0.0 || Rrs1 * Rrs2 > 0.0) && minRrs > -0.001) {
        rat = MAX(MAX(Rrs1, Rrs2), Rrs3) / Rrs4;
        if (rat > minrat && rat < maxrat) {
            rat = log10(rat);
            chl = (float)
                    pow(10.0, (a[0] + rat * (a[1] + rat * (a[2] + rat * (a[3] + rat * a[4])))));
            chl = (chl > chlmin ? chl : chlmin);
            chl = (chl < chlmax ? chl : chlmax);
        }
    }

    return (chl);
}

float chl_hu(l2str *l2rec, float Rrs[]) {
    static float w[] = {443, 555, 670};
    static float a[] = {-0.4287, 230.47};
    static int ib1 = -1;
    static int ib2 = -1;
    static int ib3 = -1;

    float ci;
    float Rrs1, Rrs2, Rrs3;
    float chl = chlbad;

    if (ib1 == -1) {

        ib1 = bindex_get(443);
        ib2 = bindex_get_555();
        ib3 = bindex_get(670);
        if (ib3 < 0) ib3 = bindex_get(665);
        if (ib3 < 0) ib3 = bindex_get(655);
        if (ib3 < 0) ib3 = bindex_get(620); // for OCM2: need to add band correction 
        if (ib1 < 0 || ib2 < 0 || ib3 < 0) {
            printf("chl_hu: incompatible sensor wavelengths for this algorithm\n");
            printf("chl_hu: %d %d %d\n", ib1, ib2, ib3);
            exit(1);
        } else {
            printf("chl_hu: using %7.2f %7.2f %7.2f\n", l2rec->l1rec->l1file->fwave[ib1], l2rec->l1rec->l1file->fwave[ib2], l2rec->l1rec->l1file->fwave[ib3]);
        }
    }

    Rrs1 = Rrs[ib1];
    Rrs2 = Rrs[ib2];
    Rrs3 = Rrs[ib3];

    // If all Rrs are bad then the pixel is masked, so return badval. For 
    // any other case, we want to return a valid value so that OCI can decide.

    if (Rrs3 > BAD_FLT + 1 && Rrs2 > BAD_FLT + 1 && Rrs1 > BAD_FLT + 1) {

        // The chl index (ci) is negative line height; ci=0 will yield chl > 0.3,
        // where the algorithm is not considered valid, and not used in the
        // merged algorithm. 

        ci = 0;

        // We require that the red channel Rrs is valid (may be slightly negative
        // in clear water due to noise), and that Rrs in blue and green be positive.
        // Cases of negative blue radiance are likely high chl anyway.

        if (Rrs3 > BAD_FLT + 1 && Rrs2 > 0.0 && Rrs1 > 0.0) {
            // shift Rrs2 to 555
            Rrs2 = conv_rrs_to_555(Rrs2, l2rec->l1rec->l1file->fwave[ib2],-99, NULL);
            // compute index
            ci = MIN(Rrs2 - (Rrs1 + (w[1] - w[0]) / (w[2] - w[0])*(Rrs3 - Rrs1)), 0.0);
            // index should be negative in algorithm-validity range
        }
        chl = (float) pow(10.0, a[0] + a[1] * ci);
        chl = (chl > chlmin ? chl : chlmin);
        chl = (chl < chlmax ? chl : chlmax);
    }

    return (chl);
}

float chl_oci(l2str *l2rec, float Rrs[]) {
    static float t1 = 0.25;
    static float t2 = 0.35;

    float chl1 = chlbad;
    float chl2 = chlbad;
    float chl = chlbad;

    float dchl1,dchl2;

    uncertainty_t *uncertainty=l2rec->l1rec->uncertainty;
    float *dchl;
    if(uncertainty)
        dchl=&uncertainty->dchl;

    chl1 = chl_hu(l2rec, Rrs);
    if(uncertainty)
        dchl1=*dchl;

    if (chl1 <= t1)
        chl = chl1;
    else {
        chl2 = get_chl_ocx(l2rec, Rrs);
        if(uncertainty)
            dchl2=*dchl;
        if (chl2 > 0.0) {
            if (chl1 >= t2)
                chl = chl2;
            else {
                chl = chl1 * (t2 - chl1) / (t2 - t1)
                        + chl2 * (chl1 - t1) / (t2 - t1);
                if(uncertainty)
                    *dchl= sqrt(pow((t2+chl2-2*chl1)*dchl1,2) + pow((chl1-t1)*dchl2,2) )/(t2-t1);
            }
        }
    }

    return (chl);
}

float chl_cdr(l2str *l2rec, float Rrs[]) {
    static float chl_lo = 0.35;
    static float chl_hi = 20.0;

    static float c_modisa[6] = {1.030e+00, 7.668e-02, 4.152e-01, 5.335e-01, -5.040e-01, 1.265e-01};
    static float c_meris [6] = {1.0, 0.0, 0.0, 0.0, 0.0};
    static float c_null [6] = {1.0, 0.0, 0.0, 0.0, 0.0};

    float chl1 = chlbad;
    float chl2 = chlbad;
    float lchl = chlbad;
    float chl = chlbad;
    float *c;
    float rat;

    chl1 = get_default_chl(l2rec, Rrs);

    if (chl1 > 0.0) {
        switch (l2rec->l1rec->l1file->sensorID) {
        case MODISA:
            c = c_modisa;
            break;
        case MERIS:
            c = c_meris;
            break;
        case SEAWIFS:
        case OCTS:
        case OCM1:
        case OCM2:
        case MOS:
        case HICO:
        case MODIST:
        case CZCS:
        case OSMI:
        case VIIRSN:
        case VIIRSJ1:
        case VIIRSJ2:
        case OCRVC:
        case GOCI:
            c = c_null;
            break;
        default:
            printf("%s Line %d: need a default chlorophyll algorithm for this sensor\n",
                    __FILE__, __LINE__);
            exit(1);
            break;
        }

        chl = MAX(MIN(chl_hi, chl1), chl_lo);
        lchl = log10(chl);
        rat = c[0] + lchl * (c[1] + lchl * c[2] + lchl * (c[3] + lchl * (c[4] + lchl * c[5])));

        chl2 = chl1 / rat;
    }

    return (chl2);
}

/*======================================================================*
 *                                                                      *
 * Implementation of ABI Chlorophyll (Shanmugam, 2011)                  *
 * A new bio-optical algorithm for the remote sensing of algal blooms   *
 * in complex ocean waters                                              *
 *                                                                      *
 * Journal of Geophysical Research, 116(C4), C04016.                    *
 * doi:10.1029/2010JC006796                                             *
 *======================================================================*/

float chl_abi(l2str *l2rec, float nLw[]) {
    float chl1, chl2, chl3, chl4, chl5, chl6, chl7, Z;
    float nLw1, nLw2, nLw3;
    float chl = chlbad;
    static int ib1 = -1, ib2 = -1, ib3 = -1;

    if (ib1 == -1) {
        ib1 = bindex_get(443);

        ib2 = bindex_get(488);
        if (ib2 < 0) ib2 = bindex_get(490);

        ib3 = bindex_get(547);
        if (ib3 < 0) ib3 = bindex_get(555);

        if (ib1 < 0 || ib2 < 0 || ib3 < 0) {
            printf("chl_abi is not compatible with this sensor\n");
            exit(1);
        } else {
            printf("Calculating chl_abi using %7.2f %7.2f %7.2f\n",
                    l2rec->l1rec->l1file->fwave[ib1],
                    l2rec->l1rec->l1file->fwave[ib2],
                    l2rec->l1rec->l1file->fwave[ib3]);
        }
    }




    nLw1 = nLw[ib1]; //443nm
    nLw2 = nLw[ib2]; //488nm
    nLw3 = nLw[ib3]; //547nm
    if (nLw1 >-10.0) { //condition to avoid negative nLw
        chl1 = ((nLw2 / nLw3) - nLw1) / ((nLw2 / nLw3) + nLw1);
        chl2 = pow(10.0, chl1);
        chl3 = (((nLw1 * nLw2) / 47.0)*((nLw1 * nLw2) / (nLw3 * nLw3)));
        chl4 = chl2*chl3;
        chl5 = 0.1403 * (pow(chl4, (-0.572)));
        chl6 = pow(chl5, 1.056);
        Z = pow(98, 1.056);
        chl7 = chl6 / (chl6 + Z);
        chl = 125.0 * chl7;
    }
    return (chl);
}

float get_default_chl(l2str *l2rec, float Rrs[]) {

    switch (l2rec->l1rec->l1file->sensorID) {
    case MOS:
        return (chl_oc4(l2rec, Rrs));
    case OSMI:
        return (chl_oc3(l2rec, Rrs));
    case L5TM:
    case L7ETMP:
    case MISR:
        return (chl_oc2(l2rec, Rrs));
    default:
        return (chl_oci(l2rec, Rrs));
    }
}

float get_chl_ocx(l2str *l2rec, float Rrs[]) {
    float chl;

    chl = chlbad;

    switch (l2rec->l1rec->l1file->sensorID) {
    case SEAWIFS:
    case OCTS:
    case OCM1:
    case OCM2:
    case MOS:
    case MERIS:
    case HICO:
    case OCIA:
    case OCI:
    case OCIS:
    case HAWKEYE:
    case AVIRIS:
    case OLCIS3A:
    case OLCIS3B:
        chl = chl_oc4(l2rec, Rrs);
        break;
    case MODIST:
    case MODISA:
    case CZCS:
    case OSMI:
    case VIIRSN:
    case VIIRSJ1:
    case VIIRSJ2:
    case OCRVC:
    case GOCI:
    case OLIL8:
    case OLIL9:
    case PRISM:
    case SGLI:
    case MSIS2A:
    case MSIS2B:
        chl = chl_oc3(l2rec,Rrs);
        break;
    case L5TM:
    case L7ETMP:
    case MISR:
        chl = chl_oc2(l2rec,Rrs);
        break;
    default:
        printf("%s Line %d: need a default chlorophyll algorithm for this sensor\n",
                __FILE__,__LINE__);
        exit(1);
        break;
    }

    return (chl);
}

void get_chl(l2str *l2rec, int prodnum, float prod[]) {
    int32_t ip;
    int32_t ipb;

    int32_t nbands = l2rec->l1rec->l1file->nbands;
    int32_t npix = l2rec->l1rec->npix;

    /*                                                      */
    /* Compute desired products at each pixel               */
    /*                                                      */
    for (ip = 0; ip < npix; ip++) {

        ipb = ip*nbands;

        switch (prodnum) {

        case DEFAULT_CHL:
            prod[ip] = get_default_chl(l2rec, &l2rec->Rrs[ipb]);
            break;
        case CAT_chl_oc2:
            prod[ip] = chl_oc2(l2rec, &l2rec->Rrs[ipb]);
            break;
        case CAT_chl_oc3:
            prod[ip] = chl_oc3(l2rec, &l2rec->Rrs[ipb]);
            break;
        case CAT_chl_oc3c:
            prod[ip] = chl_oc3c(l2rec, &l2rec->Rrs[ipb]);
            break;
        case CAT_chl_oc4:
            prod[ip] = chl_oc4(l2rec, &l2rec->Rrs[ipb]);
            break;
        case CAT_chl_hu:
            prod[ip] = chl_hu(l2rec, &l2rec->Rrs[ipb]);
            break;
        case CAT_chl_oci:
            prod[ip] = chl_oci(l2rec, &l2rec->Rrs[ipb]);
            break;
        case CAT_chl_cdr:
            prod[ip] = chl_cdr(l2rec, &l2rec->Rrs[ipb]);
            break;
        case CAT_chl_abi:
            prod[ip] = chl_abi(l2rec, &l2rec->nLw[ipb]);
            break;
        default:
            printf("Error: %s : Unknown product specifier: %d\n", __FILE__, prodnum);
            exit(FATAL_ERROR);
            break;
        }

        if (prod[ip] == chlbad)
            l2rec->l1rec->flags[ip] |= PRODFAIL;
    }
}

