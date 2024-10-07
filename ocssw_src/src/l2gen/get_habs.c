/*
 * get_hab.c
 *
 * Harmful Algal Blooms
 *  Created on: Aug 31, 2015
 *      Author: Rick Healy (richard.healy@nasa.gov)
 */
#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"
#include <sensorInfo.h>
#include "mph_flags.h"
/*
 * Harmful Algal Bloom Indexes and related functions
 * for flagging products and cloud masking
 * R. Healy 9/1/2015 (richard.healy@nasa.gov)
 *
 * Cyanobacteria Index
 * Wynne,T.T.; Stumpf, R.P.; 2013. Spatial and Temporal Patterns in the Seasonal Distribution of
            Toxic Cyanobacteria in Western Lake Erie from 2002–2014,Toxins 2015, 7, 1649-1663; doi:10.3390/toxins7051649

   Maximum Chlorophyll Index
    C.E. Binding ⁎, T.A. Greenberg, R.P. Bukata, The MERIS Maximum Chlorophyll Index; its merits and limitations for inland water
                algal bloom monitoring, JGLR-00579
 *
 */

static int numFlagMPHPixels = -1;
static int numFlagHABSPixels = -1;
static uint8_t *flags_mph = NULL;
static uint8_t *flags_habs = NULL;
static int32_t mask = LAND;

void allocateFlagsMPH(int numPix) {
    if((numFlagMPHPixels != numPix) || !flags_mph) {
        numFlagMPHPixels = numPix;
        if(flags_mph)
            free(flags_mph);
        flags_mph = (uint8_t*) malloc(numPix);
    }
}
void allocateFlagsHABS(int numPix) {
    if((numFlagHABSPixels != numPix) || !flags_habs) {
        numFlagHABSPixels = numPix;
        if(flags_habs)
            free(flags_habs);
        flags_habs = (uint8_t*) malloc(numPix);
    }
}

void habs_meris_ci_corr(float rhos443, float rhos490, float rhos560,
        float rhos620, float rhos665, float rhos709, float rhos865,
        float rhos885, float *ci) {

    float kd, kd_709, ss_560;

    kd = (0.7 * ((((rhos620 + rhos665) / 2.) - rhos885) / (((rhos443 + rhos490) / 2.) - rhos885)));

    //calculate Kd using rhos_709 instead of rhos_620 & rhos_665
    kd_709 = (0.7 * ((rhos709 - rhos885) / (((rhos443 + rhos490) / 2.) - rhos885)));

    //709 switching modification for scum conditions
    if (kd_709 > kd) kd = kd_709;

    //expect a 560 peak with cyano blooms
    ss_560 = rhos560 - rhos443 + (rhos443 - rhos620)*(560 - 443) / (620 - 443);

    //identify over-corrected pixels that would result in negative Kd's
    if (((((rhos620 + rhos665) / 2.) - rhos885) > 0) &&
            ((((rhos443 + rhos490) / 2.) - rhos885) > 0)) {
        if ((kd < 0.25) &&
            ((rhos865 <= rhos490) ||
            (rhos865 <= rhos665) ||
            (rhos865 <= rhos709)) &&
            (ss_560 < 0.01)) {
            *ci = 0;
        }
    } else {
        if (rhos885 < 0.005) {
            *ci = 0;
        }
    }
}

void get_habs_ci(l2str *l2rec, l2prodstr *p, float ci[]) {

    int ib0, ib1, ib2, ib3, ib4, ib5, ib6, ib7, ib8;
    float wav0, wav1, wav2, wav3, wav4, wav5, wav6, wav7, wav8, fac, w681_fac;
    int nonci;

    int ip, ipb;
    static productInfo_t* ci_product_info;

    switch (l2rec->l1rec->l1file->sensorID) {
    case MODISA:
    case MODIST:
        fac = 1.3;
        w681_fac = 1.0;
        break;
    case OLCIS3A:
    case OLCIS3B:        
        fac = 1.052;
        w681_fac = 0.99;
        break;
    default:
        fac = 1.0;
        w681_fac = 1.0;        
    }
    switch (p->cat_ix) {

    case CAT_CI_stumpf:
    case CAT_CI_cyano:
    case CAT_CI_noncyano:
        // Cyanobacteria Index
        // Wynne, Stumpf algorithm 2013
        switch (l2rec->l1rec->l1file->sensorID) {
        case MODISA:
        case MODIST:
            wav0 = 547;
            wav1 = 667;
            wav2 = 678;
            wav3 = 748;
            break;
        default:
            wav0 = 620;
            wav1 = 665;
            wav2 = 681;
            wav3 = 709;
            wav4 = 443;
            wav5 = 490;
            wav6 = 560;
            wav7 = 865;
            wav8 = 885;
            ib4 = bindex_get(wav4);
            ib5 = bindex_get(wav5);
            ib6 = bindex_get(wav6);
            ib7 = bindex_get(wav7);
            ib8 = bindex_get(wav8);
        }
        ib0 = bindex_get(wav0);
        break;

    case CAT_MCI_stumpf:
        if (l2rec->l1rec->l1file->sensorID != MERIS && 
                l2rec->l1rec->l1file->sensorID != OLCIS3A && 
                l2rec->l1rec->l1file->sensorID != OLCIS3B) {
            printf("MCI not supported for this sensor (%s).\n",
                    sensorId2SensorName(l2rec->l1rec->l1file->sensorID));
            exit(EXIT_FAILURE);
        }
        wav1 = 681;
        wav2 = 709;
        wav3 = 754;
        break;

    default:
        printf("HABS_CI: Hmm, something's really messed up.\n");
        exit(EXIT_FAILURE);
    }

    if (ci_product_info == NULL) {
        ci_product_info = allocateProductInfo();
        findProductInfo("CI_cyano", l2rec->l1rec->l1file->sensorID, ci_product_info);
    }

    ib1 = bindex_get(wav1);
    ib2 = bindex_get(wav2);
    ib3 = bindex_get(wav3);

    if (ib1 < 0 || ib2 < 0 || ib3 < 0) {
        printf("(M)CI_stumpf: incompatible sensor wavelengths for this algorithm\n");
        exit(EXIT_FAILURE);
    }

    for (ip = 0; ip < l2rec->l1rec->npix; ip++) {

        ipb = l2rec->l1rec->l1file->nbands*ip;
        //  check for near coastal and inland waters (height == 0 and depth < shallow_water_depth)
        //  and valid data for processing CI
        if ((l2rec->l1rec->height[ip] == 0 && l2rec->l1rec->dem[ip] < -1 * input->shallow_water_depth) ||
                (l2rec->l1rec->flags[ip] & mask) != 0 ||
                l2rec->l1rec->rhos[ipb + ib1] <= 0.0 ||
                l2rec->l1rec->rhos[ipb + ib2] <= 0.0 ||
                l2rec->l1rec->rhos[ipb + ib3] <= 0.0) {
            ci[ip] = BAD_FLT;
            l2rec->l1rec->flags[ip] |= PRODFAIL;
        } else {
            switch (p->cat_ix) {

            case CAT_CI_stumpf:
            case CAT_CI_cyano:
            case CAT_CI_noncyano:
                ci[ip] = fac * ((l2rec->l1rec->rhos[ipb + ib3] - l2rec->l1rec->rhos[ipb + ib1])*(wav2 - wav1) / (wav3 - wav1)
                        - (l2rec->l1rec->rhos[ipb + ib2]*w681_fac - l2rec->l1rec->rhos[ipb + ib1]));
           
                //following corrections currently only applicable for MERIS/OLCI
                if (l2rec->l1rec->l1file->sensorID == MERIS || l2rec->l1rec->l1file->sensorID == OLCIS3A 
                        || l2rec->l1rec->l1file->sensorID == OLCIS3B) {

                    //turbidity correction based on ss620
                    if (l2rec->l1rec->rhos[ipb + ib0] - l2rec->l1rec->rhos[ipb + ib6]
                        + (l2rec->l1rec->rhos[ipb + ib6] - l2rec->l1rec->rhos[ipb + ib1]) * (wav0 - wav6) / (wav1 - wav6) > 0) {
                            ci[ip] = 0;
                    }

                    habs_meris_ci_corr(l2rec->l1rec->rhos[ipb + ib4], l2rec->l1rec->rhos[ipb + ib5],
                                       l2rec->l1rec->rhos[ipb + ib6], l2rec->l1rec->rhos[ipb + ib0],
                                       l2rec->l1rec->rhos[ipb + ib1], l2rec->l1rec->rhos[ipb + ib3],
                                       l2rec->l1rec->rhos[ipb + ib7], l2rec->l1rec->rhos[ipb + ib8], &ci[ip]);

                    if (p->cat_ix == CAT_CI_cyano || p->cat_ix == CAT_CI_noncyano) {
                        if (l2rec->l1rec->rhos[ipb + ib1]
                                - l2rec->l1rec->rhos[ipb + ib0]
                                + (l2rec->l1rec->rhos[ipb + ib0]
                                - l2rec->l1rec->rhos[ipb + ib2]*w681_fac)
                                *(wav1 - wav0) / (wav2 - wav0) >= 0) {
                            nonci = 0;
                        } else {
                            nonci = 1;
                        }

                        if (l2rec->l1rec->rhos[ipb + ib0] <= 0 && p->cat_ix == CAT_CI_noncyano) {
                            ci[ip] = BAD_FLT;
                            l2rec->l1rec->flags[ip] |= PRODFAIL;
                        } else {
                            if (p->cat_ix == CAT_CI_noncyano) {

                                if (nonci == 0) ci[ip] = 0;
                            } else {
                                if (nonci == 1) ci[ip] = 0;
                            }
                        }
                    }
                }
                break;
            case CAT_MCI_stumpf:
                ci[ip] = fac * (l2rec->l1rec->rhos[ipb + ib2] - l2rec->l1rec->rhos[ipb + ib1]*w681_fac
                        - (l2rec->l1rec->rhos[ipb + ib3] - l2rec->l1rec->rhos[ipb + ib1]*w681_fac)*(wav2 - wav1) / (wav3 - wav1));
                break;
            default:
                ci[ip] = BAD_FLT;
                break;
            }

            // set non-detect levels of CI to the minimum valid value in product.xml definition
            if (ci[ip] < ci_product_info->validMin) {
                ci[ip] = ci_product_info->validMin;
            }
        }
    }

    flags_habs = get_flags_habs(l2rec);
    for (ip = 0; ip < l2rec->l1rec->npix; ip++) {
        if (flags_habs[ip] != 0){
            ci[ip] = BAD_FLT;
            l2rec->l1rec->flags[ip] |= PRODFAIL;
        }
    }
    if (l2rec->l1rec->iscan == l2rec->l1rec->l1file->nscan) {
        freeProductInfo(ci_product_info);
    }
}
/*
 * Maximum Peak Height of chlorophyll for MERIS
 *
 * Updated: J. Scott (6/1/2018) joel.scott@nasa.gov; removed Rrs check, fixed mistyped coefficient, removed negative chl filter
 *
 * R. Healy (9/1/2015) richard.healy@nasa.gov
 *
 * Mark William Matthews , Daniel Odermatt
 * Remote Sensing of Environment 156 (2015) 374–382
 *
 *
 */
void get_habs_mph(l2str *l2rec, l2prodstr *p, float chl_mph[]) {
    float wav6, wav7, wav8, wav9, wav10, wav14;
    int ip, ipb, ib6, ib7, ib8, ib9, ib10, ib14;
    float Rmax0, Rmax1, wavmax0, wavmax1, ndvi;
    float sipf, sicf, bair, mph0, mph1;
    float *rhos = l2rec->l1rec->rhos;

    if ((l2rec->l1rec->l1file->sensorID != MERIS) && 
            (l2rec->l1rec->l1file->sensorID != OLCIS3A) &&
            (l2rec->l1rec->l1file->sensorID != OLCIS3B)) {
        printf("MPH not supported for this sensor (%s).\n",
                sensorId2SensorName(l2rec->l1rec->l1file->sensorID));
        exit(EXIT_FAILURE);
    }

    wav6 = 620;
    wav7 = 664;
    wav8 = 681;
    wav9 = 709;
    wav10 = 753;
    wav14 = 885;

    ib6 = bindex_get(wav6);
    ib7 = bindex_get(wav7);
    ib8 = bindex_get(wav8);
    ib9 = bindex_get(wav9);
    ib10 = bindex_get(wav10);
    ib14 = bindex_get(wav14);

    if (ib6 < 0 || ib7 < 0 || ib8 < 0 || ib9 < 0 || ib10 < 0 || ib14 < 0) {
        printf("MPH_stumpf: incompatible sensor wavelengths for this algorithm\n");
        exit(EXIT_FAILURE);
    }

    for (ip = 0; ip < l2rec->l1rec->npix; ip++) {
        ipb = l2rec->l1rec->l1file->nbands*ip;

        //if (l2rec->Rrs[ipb + ib6] <= 0.0) {
            //            if(l2rec->Rrs[ipb+ib6] <= 0.0 || l2rec->Rrs[ipb+ib7] <= 0.0 || l2rec->Rrs[ipb+ib8] <= 0.0
            //            || l2rec->Rrs[ipb+ib9] <= 0.0) {
            //chl_mph[ip] = BAD_FLT;
            //l2rec->l1rec->flags[ip] |= PRODFAIL;
        //} else {
        if (rhos[ipb + ib8] > rhos[ipb + ib9]) {
            wavmax0 = wav8;
            Rmax0 = rhos[ipb + ib8];
        } else {
            wavmax0 = wav9;
            Rmax0 = rhos[ipb + ib9];
        }
        if (Rmax0 > rhos[ipb + ib10]) {
            wavmax1 = wavmax0;
            Rmax1 = Rmax0;
        } else {
            wavmax1 = wav10;
            Rmax1 = rhos[ipb + ib10];
        }

        //sun-induced phycocyanin absorption fluorescence
        sipf = rhos[ipb + ib7] - rhos[ipb + ib6] - (rhos[ipb + ib8] - rhos[ipb + ib6]) * (664 - 619) / (681 - 619);
        //sun induced chlorophyll fluorescence
        sicf = rhos[ipb + ib8] - rhos[ipb + ib7] - (rhos[ipb + ib9] - rhos[ipb + ib7]) * (681 - 664) / (709 - 664);
        //normalised difference vegetation index
        ndvi = (rhos[ipb + ib14] - rhos[ipb + ib7]) / (rhos[ipb + ib14] + rhos[ipb + ib7]);
        //backscatter and absorption induced reflectance
        bair = rhos[ipb + ib9] - rhos[ipb + ib7] - (rhos[ipb + ib14] - rhos[ipb + ib7]) * (709 - 664) / (885 - 664);
        mph0 = Rmax0 - rhos[ipb + ib7] - (rhos[ipb + ib14] - rhos[ipb + ib7]) * (wavmax0 - 664) / (885 - 664);
        mph1 = Rmax1 - rhos[ipb + ib7] - (rhos[ipb + ib14] - rhos[ipb + ib7]) * (wavmax1 - 664) / (885 - 664);

        if (wavmax1 != wav10) {

            if (sicf >= 0 || sipf <= 0 || bair <= 0.002) {
                chl_mph[ip] = 5.24e9 * pow(mph0, 4) - 1.95e8 * pow(mph0, 3) + 2.46e6 * pow(mph0, 2) + 4.02e3 * mph0 + 1.97;
            } else {
                chl_mph[ip] = 22.44 * exp(35.79 * mph1);
            }
        } else {
            if (mph1 >= 0.02 || ndvi >= 0.2) {

                if (sicf < 0 && sipf > 0) {
                    chl_mph[ip] = 22.44 * exp(35.79 * mph1);

                } else {
                    chl_mph[ip] = BAD_FLT;
                }
            } else {
                chl_mph[ip] = 5.24e9 * pow(mph0, 4) - 1.95e8 * pow(mph0, 3) + 2.46e6 * pow(mph0, 2) + 4.02e3 * mph0 + 1.97;
            }
        }
        //}
        //if (chl_mph[ip] < 0.) chl_mph[ip] = 0.;
    }

}

uint8_t* get_flags_habs_mph(l2str *l2rec) {
    //MPH - Maximum Peak Height
    // from "Improved algorithm for routine monitoring of cyanobacteria and
    // eutrophication in inland and near-coastal waters"
    // Matthews and Odermatt, Remote Sensing of Environment (doi:10.1016/j.rse.2014.10.010)
    //

    float wav6, wav7, wav8, wav9, wav10, wav14;
    int ip, ipb, ib6, ib7, ib8, ib9, ib10, ib14;
    float Rmax0, Rmax1, wavmax0, wavmax1, ndvi;
    float sipf, sicf, bair, mph1, chl_mph;
    float *rhos = l2rec->l1rec->rhos;
    static float thresh = 350;

    if ((l2rec->l1rec->l1file->sensorID != MERIS) && 
            (l2rec->l1rec->l1file->sensorID != OLCIS3A) &&
            (l2rec->l1rec->l1file->sensorID != OLCIS3B)) {
        printf("MPH not supported for this sensor (%s).\n",
                sensorId2SensorName(l2rec->l1rec->l1file->sensorID));
        exit(EXIT_FAILURE);
    }

    allocateFlagsMPH(l2rec->l1rec->npix);

    wav6 = 620;
    wav7 = 664;
    wav8 = 681;
    wav9 = 709;
    wav10 = 753;
    wav14 = 885;

    ib6 = bindex_get(wav6);
    ib7 = bindex_get(wav7);
    ib8 = bindex_get(wav8);
    ib9 = bindex_get(wav9);
    ib10 = bindex_get(wav10);
    ib14 = bindex_get(wav14);

    if (ib6 < 0 || ib7 < 0 || ib8 < 0 || ib9 < 0 || ib10 < 0 || ib14 < 0) {
        printf("MPH_stumpf: incompatible sensor wavelengths for this algorithm\n");
        exit(EXIT_FAILURE);
    }

    for (ip = 0; ip < l2rec->l1rec->npix; ip++) {

        flags_mph[ip] = 0;
        ipb = l2rec->l1rec->l1file->nbands*ip;

        //     if(l2rec->rhos[ipb+ib6] <= 0.0 ) {
        if (l2rec->l1rec->rhos[ipb + ib6] <= 0.0 ||
                l2rec->l1rec->rhos[ipb + ib7] <= 0.0 ||
                l2rec->l1rec->rhos[ipb + ib8] <= 0.0 ||
                l2rec->l1rec->rhos[ipb + ib9] <= 0.0 ||
                (l2rec->l1rec->flags[ip] & LAND) != 0 ||
                (l2rec->l1rec->flags[ip] & NAVFAIL) != 0) {
            l2rec->l1rec->flags[ip] |= PRODFAIL;
            flags_mph[ip] |= MPH_BADINPUT;
        } else {
            if (rhos[ipb + ib8] > rhos[ipb + ib9]) {
                wavmax0 = wav8;
                Rmax0 = rhos[ipb + ib8];
            } else {
                wavmax0 = wav9;
                Rmax0 = rhos[ipb + ib9];
            }
            if (Rmax0 > rhos[ipb + ib10]) {
                wavmax1 = wavmax0;
                Rmax1 = Rmax0;
            } else {
                wavmax1 = wav10;
                Rmax1 = rhos[ipb + ib10];
            }

            //sun-induced phycocyanin absorption fluorescence
            sipf = rhos[ipb + ib7] - rhos[ipb + ib6] - (rhos[ipb + ib8] - rhos[ipb + ib6]) * (664 - 619) / (681 - 619);
            //sun induced chlorophyll fluorescence
            sicf = rhos[ipb + ib8] - rhos[ipb + ib7] - (rhos[ipb + ib9] - rhos[ipb + ib7]) * (681 - 664) / (709 - 664);
            //normalised difference vegetation index
            ndvi = (rhos[ipb + ib14] - rhos[ipb + ib7]) / (rhos[ipb + ib14] + rhos[ipb + ib7]);
            //backscatter and absorption induced reflectance
            bair = rhos[ipb + ib9] - rhos[ipb + ib7] - (rhos[ipb + ib14] - rhos[ipb + ib7]) * (709 - 664) / (885 - 664);
            //mph0 = Rmax0 - rhos[ipb+ib7] - ( rhos[ipb+ib14] - rhos[ipb+ib7]) * (wavmax0 - 664)/(885 - 664);
            mph1 = Rmax1 - rhos[ipb + ib7] - (rhos[ipb + ib14] - rhos[ipb + ib7]) * (wavmax1 - 664) / (885 - 664);

            if (wavmax1 != wav10) {

                if (sicf < 0 && sipf > 0 && bair > 0.002) {
                    flags_mph[ip] |= MPH_CYANO;
                    chl_mph = 22.44 * exp(35.79 * mph1);
                    if (chl_mph > thresh)
                        flags_mph[ip] |= MPH_FLOAT;
                }
            } else {
                if (mph1 >= 0.02 || ndvi >= 0.2) {
                    flags_mph[ip] |= MPH_FLOAT;

                    if (sicf < 0 && sipf > 0) {
                        flags_mph[ip] |= MPH_CYANO;
                        chl_mph = 22.44 * exp(35.79 * mph1);
                        if (chl_mph > thresh)
                            flags_mph[ip] |= MPH_FLOAT;

                    }
                } else {
                    flags_mph[ip] |= MPH_ADJ;
                }
            }

        }

    }
    return flags_mph;
}

uint8_t* get_flags_habs_meris(l2str *l2rec) {
    // Cloud Masking for MERIS & OLCI

    int ib443, ib490, ib510, ib560, ib620, ib665, ib681, ib709, ib754, ib865, ib885;
    int ip, ipb;
    float *rhos = l2rec->l1rec->rhos;
    float mdsi, cv, mean, sum, sdev, w681_fac;
    int i, n=7;

    allocateFlagsHABS(l2rec->l1rec->npix);

    if (l2rec->l1rec->l1file->sensorID != MERIS && 
            l2rec->l1rec->l1file->sensorID != OLCIS3A &&
            l2rec->l1rec->l1file->sensorID != OLCIS3B) {
        printf("HABS flags not supported for this sensor (%s).\n",
                sensorId2SensorName(l2rec->l1rec->l1file->sensorID));
        exit(EXIT_FAILURE);
    }

    ib443 = bindex_get(443);
    ib490 = bindex_get(490);
    ib510 = bindex_get(510);
    ib560 = bindex_get(560);
    ib620 = bindex_get(620);
    ib665 = bindex_get(665);
    ib681 = bindex_get(681);
    ib709 = bindex_get(709);
    ib754 = bindex_get(754);
    ib865 = bindex_get(865);
    ib885 = bindex_get(885);

    if (ib443 < 0 || ib490 < 0 || ib510 < 0 || ib560 < 0 || ib620 < 0 || ib665 < 0 || ib681 < 0 || ib709 < 0 || ib754 < 0 || ib865 < 0 || ib885 < 0) {
        printf("get_flags_habs: incompatible sensor wavelengths for this algorithm\n");
        exit(EXIT_FAILURE);
    }
    
    if (l2rec->l1rec->l1file->sensorID == OLCIS3A ||
        l2rec->l1rec->l1file->sensorID == OLCIS3B) {
        w681_fac = 0.99;
    } else {
        w681_fac = 1.0;
    }

    for (ip = 0; ip < l2rec->l1rec->npix; ip++) {
        flags_habs[ip] = 0;
        ipb = l2rec->l1rec->l1file->nbands*ip;

        // if the eval cloud masking is used, just use the cloud mask
        // otherwise, run the get_cldmask function
        if ((l1_input->evalmask | 2) == 2){
            flags_habs[ip] = (uint8_t) l2rec->l1rec->cloud[ip];
        } else {
            if (get_cldmask(l2rec->l1rec, ip)) {
                flags_habs[ip] |= HABS_CLOUD;
            }
        }

        if (rhos[ipb + ib443] >= 0.0 &&
            rhos[ipb + ib490] >= 0.0 &&
            rhos[ipb + ib510] >= 0.0 &&
            rhos[ipb + ib560] >= 0.0 &&
            rhos[ipb + ib620] >= 0.0 &&
            rhos[ipb + ib665] >= 0.0 &&
            rhos[ipb + ib681] >= 0.0 &&
            rhos[ipb + ib709] >= 0.0 &&
            rhos[ipb + ib754] >= 0.0 &&
            rhos[ipb + ib885] >= 0.0) {

            if (rhos[ipb + ib885] > rhos[ipb + ib620] &&
                    rhos[ipb + ib885] > rhos[ipb + ib709] &&
                    rhos[ipb + ib885] > rhos[ipb + ib754] &&
                    rhos[ipb + ib885] > 0.01) {
                flags_habs[ip] |= HABS_NONWTR;
            }

            // test dry lake condition: (rhos_620 > rhos_560)
            if ((rhos[ipb + ib620] > rhos[ipb + ib560]) &&
                (rhos[ipb + ib560] > 0.15) &&
                (rhos[ipb + ib885] > 0.15)){
                flags_habs[ip] |= HABS_NONWTR;
            }

            // test snow/ice condition - meris differential snow index
            mdsi = (rhos[ipb + ib865] - rhos[ipb + ib885]) / (rhos[ipb + ib865] + rhos[ipb + ib885]);

            //  test snow/ice condition - exclude potential bloom conditions (based on coefficient of variation in visible bands)
            int b[7] = {ib443, ib490, ib510, ib560, ib620, ib665, ib681};
            sum = 0;
            for(i = 0; i < n; ++i) {
                sum += rhos[ipb + b[i]];
            }
            mean = sum / n;

            sdev = 0;
            for(i = 0; i < n; ++i) {
                sdev += pow((rhos[ipb + b[i]] - mean),2);
            }
            cv = sqrt(sdev / (n - 1)) / mean;

            if ((mdsi > 0.01) && (rhos[ipb + ib885] > 0.15) && (cv < 0.1)) {
                flags_habs[ip] |= HABS_SNOWICE;
            }
            //adjacency flagging
            float mci = 1.0 * (rhos[ipb + ib709] - rhos[ipb + ib681]*w681_fac
                + (rhos[ipb + ib681]*w681_fac - rhos[ipb + ib754])*(709 - 681) / (754 - 681));
            float ci = 1.0 * ((rhos[ipb + ib709] - rhos[ipb + ib665])*(681 - 665) / (709 - 665)
                - (rhos[ipb + ib681]*w681_fac - rhos[ipb + ib665]));

            habs_meris_ci_corr(rhos[ipb + ib443], rhos[ipb + ib490],
                    rhos[ipb + ib560], rhos[ipb + ib620],
                    rhos[ipb + ib665], rhos[ipb + ib709],
                    rhos[ipb + ib865], rhos[ipb + ib885], &ci);

            if ((mci < 0) && (ci > 0)) {
                flags_habs[ip] |= HABS_ADJ;
            }
        } else {
            if (rhos[ipb + ib620] < 0.0 ||
                rhos[ipb + ib665] < 0.0 ||
                rhos[ipb + ib681] < 0.0 ||
                rhos[ipb + ib709] < 0.0) {
                flags_habs[ip] |= HABS_BADINPUT;
            }
        }
    }
    return flags_habs;
}

uint8_t* get_flags_habs_modis(l2str *l2rec) {
    // Cloud Masking for MODIS

    int ib469, ib555, ib645, ib667, ib859, ib1240, ib2130;
    int ip, ipb;
    float *rhos = l2rec->l1rec->rhos, *Rrs = l2rec->Rrs, *cloud_albedo = l2rec->l1rec->cloud_albedo;
    float ftemp, ftemp2, ftemp3;
    float cloudthr = 0.027;

    allocateFlagsHABS(l2rec->l1rec->npix);

    if (l2rec->l1rec->l1file->sensorID != MODISA && l2rec->l1rec->l1file->sensorID != MODIST) {
        printf("HABS flags not supported for this sensor (%s).\n",
                sensorId2SensorName(l2rec->l1rec->l1file->sensorID));
        exit(EXIT_FAILURE);
    }

    ib469 = bindex_get(469);
    ib555 = bindex_get(555);
    ib645 = bindex_get(645);
    ib667 = bindex_get(667);
    ib859 = bindex_get(859);
    ib1240 = bindex_get(1240);
    ib2130 = bindex_get(2130);

    if (ib469 < 0 || ib555 < 0 || ib645 < 0 || ib667 < 0 || ib859 < 0 || ib1240 < 0 || ib2130 < 0) {
        printf("get_flags_habs: incompatible sensor wavelengths for this algorithm\n");
        exit(EXIT_FAILURE);
    }

    for (ip = 0; ip < l2rec->l1rec->npix; ip++) {
        flags_habs[ip] = 0;
        ipb = l2rec->l1rec->l1file->nbands*ip;
        // TODO: make consistent with get_cloudmask_modis and use it instead of
        //       this duplication of code...

        if (l2rec->chl[ip] < 0)
            ftemp = Rrs[ipb + ib667];
        else
            ftemp = Rrs[ipb + ib667]*(0.45 + l2rec->chl[ip]*0.005) / 4.3;
        //      first correct for turbid water
        if (Rrs[ipb + ib667] < 0.0) ftemp = 0.0;
        ftemp2 = cloud_albedo[ip] - ftemp;

        if (ftemp2 > 0.027) flags_habs[ip] |= HABS_CLOUD;
        //        non-water check  1240 is bright relative to 859 and the combination is bright
        //        this may hit glint by accident, need to be checked.

        if (rhos[ipb + ib1240] / rhos[ipb + ib859] > 0.5 && (rhos[ipb + ib1240] + rhos[ipb + ib2130]) > 0.10) flags_habs[ip] |= HABS_CLOUD;

        //        now try to correct for glint
        //        region check was thrown out {IF (region = "OM") cloudthr = 0.04} rjh 11/2/2015

        ftemp = rhos[ipb + ib645] - rhos[ipb + ib555] + (rhos[ipb + ib555] - rhos[ipb + ib859])*(645.0 - 555.0) / (859.0 - 555.0);
        ftemp2 = cloud_albedo[ip] + ftemp;
        if (ftemp2 < cloudthr) flags_habs[ip] = 0;
        if (rhos[ipb + ib859] / rhos[ipb + ib1240] > 4.0) flags_habs[ip] = 0;

        //     scum areas

        if ((rhos[ipb + ib859] - rhos[ipb + ib469]) > 0.01 && cloud_albedo[ip] < 0.30) flags_habs[ip] = 0;
        if ((rhos[ipb + ib859] - rhos[ipb + ib645]) > 0.01 && cloud_albedo[ip] < 0.15) flags_habs[ip] = 0;
        if (rhos[ipb + ib1240] < 0.2)
            ftemp2 = ftemp2 - (rhos[ipb + ib859] - rhos[ipb + ib1240]) * fabs(rhos[ipb + ib859] - rhos[ipb + ib1240]) / cloudthr;

        ftemp3 = ftemp2;
        if (ftemp2 < cloudthr * 2) {
            if ((rhos[ipb + ib555] - rhos[ipb + ib1240]) > (rhos[ipb + ib469] - rhos[ipb + ib1240])) {
                ftemp3 = ftemp2 - (rhos[ipb + ib555] - rhos[ipb + ib1240]);
            } else {
                ftemp3 = ftemp2 - (rhos[ipb + ib469] - rhos[ipb + ib1240]);
            }
        }

        if (ftemp3 < cloudthr) flags_habs[ip] = 0;

        if (rhos[ipb + ib555] >= 0 && rhos[ipb + ib1240] >= 0.0 && rhos[ipb + ib1240] > rhos[ipb + ib555])
            flags_habs[ip] |= HABS_NONWTR;
    }
    return flags_habs;
}

uint8_t* get_flags_habs(l2str *l2rec) {
    static int32_t currentLine = -1;

    if (currentLine == l2rec->l1rec->iscan )
        return flags_habs;
    currentLine = l2rec->l1rec->iscan;
    switch (l2rec->l1rec->l1file->sensorID) {
    case MERIS:
    case OLCIS3A:
    case OLCIS3B:
        return get_flags_habs_meris(l2rec);
        break;
    case MODISA:
    case MODIST:
        return get_flags_habs_modis(l2rec);
        break;
    default:
        printf("HABS flags not supported for this sensor (%s).\n",
                sensorId2SensorName(l2rec->l1rec->l1file->sensorID));
        exit(EXIT_FAILURE);
    }
    return NULL;
}
