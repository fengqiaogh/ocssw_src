/*
 * get_Cphyt.c
 *
 * This algorithm returns the concentration of the phytoplankton carbon (Cphyt) in mg m^-3, calculated using an
 * empirical relationship derived from field data between analytical measurements of Cphyt and particulate
 * backscattering coefficient.
 *
 * Reference:
 * Graff, J.R., Westberry, T.K., Milligan, A.J., Brown, M.B., Dall’Olmo, G., Dongen-Vogels, V.v., Reifel, K.M., & Behrenfeld, M.J. (2015).
 * Analytical phytoplankton carbon measurements spanning diverse ecosystems.
 * Deep Sea Research Part I: Oceanographic Research Papers, 102, 16-25
 *
 *  Created on: Mar 4, 2021
 *  Author: Minwei Zhang
 */

#include "l12_proto.h"

void get_Cphyt(l2str *l2rec, float cphyt[])
{
    float bbp_470;
    float a1=12128., b1=0.59;

    int32_t ip;
    int32_t ipb;

    static int32_t firstcall=1,ib443, nbands;
    int32_t npix;
    float *bbp,*bbp_s;
    l1str *l1rec=l2rec->l1rec;

    if(firstcall){
        firstcall=0;
        ib443=bindex_get(443);
        if (ib443 < 0) {
            printf("get_Cphyt: incompatible sensor wavelengths (no 443 for bbp and bbp_s).\n");
            exit(1);
        }
        nbands=l1rec->l1file->nbands;
    }

    npix = l1rec->npix;
    if (!giop_ran(l1rec->iscan)){
        run_giop(l2rec);

    }
    bbp=giop_get_bbp_pointer();
    bbp_s=giop_get_bbp_s_pointer();
    /*                                                      */
    /* Compute desired products at each pixel               */
    /*                                                      */
    for (ip = 0; ip < npix; ip++) {
        cphyt[ip]=BAD_FLT;

        if(l2rec->l1rec->mask[ip] || l2rec->l1rec->solz[ip] >= SOLZNIGHT)
            continue;

        ipb=ip*nbands+ib443;

        if(bbp[ipb]<0 || bbp[ipb]>0.15 || bbp_s[ip]<-3 || bbp_s[ip]>3)
            continue;

        if(input->cphyt_opt==1){
            bbp_470=bbp[ipb]*pow(443./470.,bbp_s[ip]);
            cphyt[ip]=a1*bbp_470+b1;
        }else if(input->cphyt_opt==2){
            cphyt[ip]=13000.*bbp[ipb]-0.00035;
        }else{
            printf("get_Cphyt: unrecognized cphyt_opt=%d\n", input->cphyt_opt);
            exit(1);
        }

    }

}


/* ---------------------------------------------------------------------- */
/* get_Cphyt_unc() computes uncertaintis for Cphyt models                 */
/*                                                                        */
/* Author: L. McKinna (GO2Q)                                              */
/* Created: Oct 2023                                                      */
/* -----------------------------------------------------------------------*/


void get_Cphyt_unc(l2str *l2rec, float cphyt_unc[]) {

    float bbp_470, ubbp_470;
    //Model coefficients
    float a1=12128.; 
    //float b1=0.59;
    float a2 = 13000.; 
    //float b2 = 0.00035;

    /*Uncertainties unreported. Update in next revision*/
    float ua1 = 0.0; 
    float ub1 = 0.0;
    float ua2 = 0.0; 
    float ub2 = 0.0;

    float dyd1, dyd2, dcda, dcdb, dcdbbp;

    int32_t ip;
    int32_t ipb;
    int32_t ierr;

    static int32_t firstcall=1,ib443, nbands;
    int32_t npix;
    float *bbp,*bbp_s, *ubbp_s;
    l1str *l1rec=l2rec->l1rec;


    if(firstcall){
        firstcall=0;
        ib443=bindex_get(443);
        if (ib443 < 0) {
            printf("get_Cphyt: incompatible sensor wavelengths (no 443 for bbp and bbp_s).\n");
            exit(1);
        }
        nbands=l1rec->l1file->nbands;
    }

    npix = l1rec->npix;
    if (!giop_ran(l1rec->iscan)){
        run_giop(l2rec);

    }
    bbp=giop_get_bbp_pointer();
    bbp_s=giop_get_bbp_s_pointer();
    ubbp_s=giop_get_ubbp_s_pointer();
    /*                                                      */
    /* Compute desired products at each pixel               */
    /*                                                      */
    for (ip = 0; ip < npix; ip++) {

        ipb = ip*nbands;
        ierr = l1rec->npix * nbands + ipb;
        
        cphyt_unc[ip] = BAD_FLT;

        if(l2rec->l1rec->mask[ip] || l2rec->l1rec->solz[ip] >= SOLZNIGHT)
            continue;

        ipb=ip*nbands+ib443;
        ierr = npix*nbands + ipb;

        if(bbp[ipb]<0 || bbp[ipb]>0.15 || bbp_s[ip]<-3 || bbp_s[ip]>3)
            continue;

        if(input->cphyt_opt==1){
            bbp_470=bbp[ipb]*pow(443./470.,bbp_s[ip]);
            dyd1 = pow(470./443.,bbp_s[ip]);
            dyd2 = bbp_470 * log(443./470);
            ubbp_470 = sqrt(pow(dyd1*bbp[ierr],2.) + pow(dyd2*ubbp_s[ip],2.));

            dcda = bbp_470;
            dcdb = 0.0;
            dcdbbp = a1;

            cphyt_unc[ip] =  sqrt(pow(dcda*ua1,2.) + pow(dcdb*ub1,2.) + pow(dcdbbp*ubbp_470,2.));

        }else if(input->cphyt_opt==2){
            dcda = bbp[ipb];
            dcdb = 0.0;
            dcdbbp = a2;

            cphyt_unc[ip] =  sqrt(pow(dcda*ua2,2.) + pow(dcdb*ub2,2.) + pow(dcdbbp*bbp[ierr],2.));

        }else{
            printf("get_Cphyt: unrecognized cphyt_opt=%d\n", input->cphyt_opt);
            exit(1);
        }

    }

}
