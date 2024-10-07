#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "l12_proto.h"
#include "l1_misr.h"
#include "geo_region.h"
/* ---------------------------------------------------------- */
/* Converts a sensor-generic level-1b record to level-2       */
/*                                                            */
/* B. A. Franz, GSC, SIMBIOS Project, March 1998              */
/* W. Robinson, SAIC  15 Dec 2006  fix Western, Eastern most  */
/*     long for CZCS                                          */

/* ---------------------------------------------------------- */

int convl12(l1str *l1rec, l2str *l2rec, int32_t spix, int32_t epix,
        aestr *aerec) {
    int32_t ip; /* Pixel index        */
    int32_t status; /* 0=OK, 1=bad        */

    /*                                                      */
    /* Clear the L2 record                                  */
    /*                                                      */
    init_l2(l2rec, l1rec->l1file->nbands);

    /* Point L2 rec to computed SST, if requested (before atmcor) */
    if (input->proc_sst)
        l2rec->sst = get_sst(l2rec);
    else
        l2rec->sst = NULL;

    if (l1rec->l1file->sensorID == MISR) {
      misr_t *private_data = l1rec->l1file->private_data;
      int32_t block = l1rec->iscan / 128;
      if (private_data->multipleInput == 1) block /= 9;
      if(//private_data->isOceanBlock[block] == 0 ||
         (block+1) < private_data->startBlock   ||
         (block+1) > private_data->endBlock)
        return 0;
    }

    /*                                                      */
    /* Loop through each pixel and do atmospheric correction*/
    /*                                                      */
    for (ip = spix; ip <= epix; ip++) {
        // setting the georegion
        if(input->georegionfile[0]){
            float lat = l1rec->lat[ip];
            float lon = l1rec->lon[ip];
            if(get_georegion(lat, lon))
                l2rec->l1rec->flags[ip] |=GEOREGION;
        }
        /* ------------------------------------------------ */
        /* Ocean processing                                 */
        /* ------------------------------------------------ */
        if ((input->proc_ocean != 0) &&
            !l1rec->mask[ip] &&
            l1rec->solz[ip] < SOLZNIGHT) {
            if (l1rec->is_l2){
                /* Lt values are reflectances: skip atmocor, but calc chl */
                int nbands = l1rec->l1file->nbands;
                for (int ib = 0; ib < nbands; ib++) {
                    int ipb = ip*nbands+ib;
                    l2rec->Rrs[ipb] = l1rec->Lt[ipb];
                    l2rec->nLw[ipb] = l2rec->Rrs[ipb]*l1rec->l1file->Fobar[ib];
                }

                // l2rec->chl[ip] = get_default_chl(l2rec, &l2rec->Rrs[ip * nbands]);
                l2rec->chl[ip] = get_default_chl(l2rec, l2rec->Rrs);
            } else if (input->atmocor) {
                /* set aerosol values from input rec, if supplied */
//                if (input->aer_opt == AERSOA || input->aer_opt == AERSMA)
//                    status = run_soa_sma(l2rec, ip);
//                else
                    status = atmocor2(l2rec, aerec, ip);

                /*                                                           */
                /* If the atmospheric correction failed, flag and mask. Else,*/
                /* set flags which depend on complete atmospheric correction.*/
                /*                                                           */
                if (status != 0) {
                    l2rec->l1rec->flags[ip] |= ATMFAIL;
                    l2rec->l1rec->mask[ip] = 1;
                } else {
                    setflagbits_l2(l2rec, ip);
                }

            }

        } // if ocean

    } // for ip


    /* Load L2 rec with inherent optical properties */
    if (l1rec->is_l2 || (input->iop_opt > 0 && (input->proc_ocean != 0) && input->atmocor)){
        get_iops(l2rec, input->iop_opt);
    }

    return (0);
}



/* --------------------------------------------------------------- */
/* get_iops.c - load IOP (a & bb) fields in L2 rec.                */
/*                                                                 */
/* Inputs:                                                         */
/*     l2rec - level-2 structure containing one complete scan      */
/*             after atmospheric correction.                       */
/* Outputs:                                                        */
/*     iop_opt - algorithm selector                                */
/*                                                                 */
/* Written By: B. Franz, NASA/OBPG/SAIC, 25 Feb 2005               */
/*                                                                 */

/* --------------------------------------------------------------- */
void get_iops(l2str *l2rec, int32_t iop_opt) {
    switch (iop_opt) {
    case IOPCARDER:
        iops_carder(l2rec);
        break;
    case IOPGSM:
        iops_gsm(l2rec);
        break;
    case IOPQAA:
        iops_qaa(l2rec);
        break;
    case IOPPML:
        iops_pml(l2rec);
        break;
    case IOPLAS:
        iops_las(l2rec);
        break;
    case IOPNIWA:
        iops_niwa(l2rec);
        break;
    case IOPGIOP:
        iops_giop(l2rec);
        break;
    case IOPSWIM:
        iops_swim(l2rec);
        break;
    default:
        break;
    }

    return;
}
