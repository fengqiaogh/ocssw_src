/*
 * alloc_uncertainty.c
 *
 *  Created on: Aug 8, 2019
 *      Author: mzhang11
 */

#include <uncertainty.h>
#include <l1.h>

//#include <stdio.h>
//#include <stdlib.h>
//#include "filehandle.h"


/* --------------------------------------------------------- */
/* alloc_uncertainty() - allocates 1 uncertainty record to hold data for */
/*              a single scan of "npix" pixels.              */

/* --------------------------------------------------------- */
int alloc_uncertainty(int32_t nbands, int32_t nbands_ac, int32_t npix, uncertainty_t *uncertainty) {

    char *p;
    int32_t len;

    uncertainty->nbands=nbands;
    uncertainty->nbands_ac=nbands_ac;
    uncertainty->npix=npix;

    /*                                                      */
    /* allocate data block as contiguous bytes              */
    /*                                                      */
    len =  11 * sizeof (float)*npix
         + 11 * sizeof (float)*npix * nbands
         + 16 * sizeof (float)* nbands
         + 2 * sizeof (float)* nbands*nbands
         + 4 * sizeof (float)* nbands_ac
         + 13 * sizeof (float)* nbands*nbands_ac
         + 1  * sizeof (int)*nbands;


    /* Force to 4-byte increments for good measure */
    len = (len / 4 + 1)*4;

    if ((p = (char *) malloc(len)) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: Memory allocation failure.\n",
                __FILE__, __LINE__);
        return (0);
    }

    /* Note: positional allocation is in order of datatype size, to
       ensure that all 4-byte words start on 4-byte boundaries. Some
       machines seem to have trouble if this is not done. */

    uncertainty->data = p;
    uncertainty->length=len;

    uncertainty->dwv     =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->dmw     =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->dzw     =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->dws     =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->dwd     =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->doz     =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->dpr     =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->drh     =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->dno2_tropo =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->dno2_strat =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->derv_gc_ws =(float *) p;  p+=sizeof(float)*npix;

    uncertainty->dsensor =(float *) p;  p+=sizeof(float)*npix*nbands;
    uncertainty->dsyst   =(float *) p;  p+=sizeof(float)*npix*nbands;
    uncertainty->dvc     =(float *) p;  p+=sizeof(float)*npix*nbands;
    uncertainty->dtg_sol =(float *) p;  p+=sizeof(float)*npix*nbands;
    uncertainty->dtg_sen =(float *) p;  p+=sizeof(float)*npix*nbands;
    uncertainty->dt_sol  =(float *) p;  p+=sizeof(float)*npix*nbands;
    uncertainty->dt_sen  =(float *) p;  p+=sizeof(float)*npix*nbands;
    uncertainty->dtaua   =(float *) p;  p+=sizeof(float)*npix*nbands;
    uncertainty->derv_polcor_Lt    =(float *) p;  p+=sizeof(float)*npix*nbands;
    uncertainty->derv_polcor_tgsen =(float *) p;  p+=sizeof(float)*npix*nbands;
    uncertainty->derv_polcor_tgsol =(float *) p;  p+=sizeof(float)*npix*nbands;
   
    uncertainty->derv_La_rh    =(float  *)p; p+=sizeof(float)*nbands; 
    uncertainty->derv_taua_rh    =(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_tsen_rh    =(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_tsol_rh    =(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_taua_min_rhorc_l=(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_taua_min_taua_l=(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_taua_min_rhow_l=(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_taua_max_rhorc_l=(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_taua_max_taua_l=(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_taua_max_rhow_l=(float  *)p; p+=sizeof(float)*nbands; 
    uncertainty->drhown_nir   =(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->dbrdf        =(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_rhog_taua  = (float*)p; p += sizeof(float) * nbands;
    uncertainty->derv_rhog_rhorc_l  = (float*)p; p += sizeof(float) * nbands;
    uncertainty->derv_TLg_gc     = (float*)p; p += sizeof(float) * nbands;
    uncertainty->derv_chl_rrs    = (float*)p; p += sizeof(float) * nbands;

    uncertainty->corr_coef_rhot      = (float*)p; p += sizeof(float) * nbands * nbands;
    uncertainty->covariance_matrix   = (float*)p; p += sizeof(float) * nbands * nbands;
    

    uncertainty->derv_modrat_rhorc =(float * )p; p+=sizeof(float)*nbands_ac;
    uncertainty->derv_modrat_taua  =(float * )p; p+=sizeof(float)*nbands_ac;
    uncertainty->ratio_rhow        =(float * )p; p+=sizeof(float)*nbands_ac;
    uncertainty->derv_rhownir_chl  =(float * )p; p+=sizeof(float)*nbands_ac;

    uncertainty->derv_La_rhorc    =(float *) p; p+=sizeof(float)*nbands*nbands_ac;
    uncertainty->derv_taua_rhorc  =(float *) p; p+=sizeof(float)*nbands*nbands_ac;
    uncertainty->derv_tsen_rhorc  =(float *) p; p+=sizeof(float)*nbands*nbands_ac;
    uncertainty->derv_tsol_rhorc  =(float *) p; p+=sizeof(float)*nbands*nbands_ac;
    uncertainty->derv_La_taua_l   =(float  *)p; p+=sizeof(float)*nbands*nbands_ac;
    uncertainty->derv_taua_taua_l =(float  *)p; p+=sizeof(float)*nbands*nbands_ac;
    uncertainty->derv_tsol_taua   =(float  *)p; p+=sizeof(float)*nbands*nbands_ac;
    uncertainty->derv_tsen_taua   =(float  *)p; p+=sizeof(float)*nbands*nbands_ac;
    uncertainty->derv_rhownir_rrs =(float  *)p; p+=sizeof(float)*nbands*nbands_ac;
    uncertainty->derv_La_rhow     =(float  *)p; p+=sizeof(float)*nbands*nbands_ac;
    uncertainty->derv_taua_rhow   =(float  *)p; p+=sizeof(float)*nbands*nbands_ac;
    uncertainty->derv_tsen_rhow   =(float  *)p; p+=sizeof(float)*nbands*nbands_ac;
    uncertainty->derv_tsol_rhow   =(float  *)p; p+=sizeof(float)*nbands*nbands_ac;

    uncertainty->bindx_rhownir    =(int *)p; p+=sizeof(int)*nbands;

    if ((len - (int32_t) (p - uncertainty->data)) < 0) {
        printf("%s Line %d: bad allocation on error record\n", __FILE__, __LINE__);
        exit(1);
    }
    if (want_verbose)
        printf("Allocated %d bytes in error record.\n", (int) (p - uncertainty->data));
   
    return (0);
}

/*
 *
 *ifscan=1: initializing uncertainty for a scan
 *ifscan=0:  initializing uncertainty for one pixel
 */
void init_uncertainty(uncertainty_t *uncertainty, int ifscan){
    int32_t ip, ib, ipb,inir;

    int32_t npix = uncertainty->npix;
    int32_t nbands = uncertainty->nbands;
    int32_t nbands_ac = uncertainty->nbands_ac;

    if (ifscan) {
        for (ip = 0; ip < npix; ip++) {
            uncertainty->dwv[ip] = 0.;
            uncertainty->dmw[ip] = 0.;
            uncertainty->dzw[ip] = 0.;
            uncertainty->dws[ip] = 0.;
            uncertainty->dwd[ip] = 0.;
            uncertainty->doz[ip] = 0.;
            uncertainty->dpr[ip] = 0.;
            uncertainty->drh[ip] = 0.;
            uncertainty->dno2_tropo[ip] = 0.;
            uncertainty->dno2_strat[ip] = 0.;
            uncertainty->derv_gc_ws[ip] = 0.;

            for (ib = 0; ib < nbands; ib++) {
                ipb = ip * nbands + ib;
                uncertainty->dsensor[ipb] = 0.0;
                uncertainty->dsyst[ipb] = 0.0;
                uncertainty->dvc[ipb] = 0.0;

                uncertainty->dtg_sol[ipb] = 0.0;
                uncertainty->dtg_sen[ipb] = 0.0;
                uncertainty->dt_sol[ipb] = 0.0;
                uncertainty->dt_sen[ipb] = 0.0;
                uncertainty->dtaua[ipb] = 0.0;
                uncertainty->derv_polcor_Lt[ipb] = 0.0;
                uncertainty->derv_polcor_tgsol[ipb] = 0.0;
                uncertainty->derv_polcor_tgsen[ipb] = 0.0;

                if (ip == 0)
                    uncertainty->bindx_rhownir[ib] = 0;
            }
        }

        for (ib = 0; ib < nbands * nbands; ib++)
            uncertainty->corr_coef_rhot[ib] = 0.;

        return;
    }

    uncertainty->derv_modrat_taua_l = 0.;
    uncertainty->derv_modrat_rhow_l = 0.0;
    uncertainty->derv_eps_rhorc_s = 0.;
    uncertainty->derv_eps_rhorc_l = 0.;
    //uncertainty->derv_eps_taua_l = 0;
   // uncertainty->derv_eps_taua_s = 0;
   // uncertainty->derv_eps_rhow_l = 0;
    uncertainty->dchl = 0;
    uncertainty->dkd490 = 0;
    uncertainty->aer_l_glint=0;

    for (inir = 0; inir < nbands_ac; inir++) {
        uncertainty->derv_modrat_rhorc[inir] = 0.;
        uncertainty->derv_modrat_taua[inir] =0.;
        uncertainty->ratio_rhow[inir] = 0.;
        uncertainty->derv_rhownir_chl[inir]=0.;
    }

    for (ib = 0; ib < nbands; ib++) {
        uncertainty->derv_taua_min_rhorc_l[ib] = 0.;
        uncertainty->derv_taua_max_rhorc_l[ib] = 0.;
        uncertainty->derv_taua_min_taua_l[ib] = 0.;
        uncertainty->derv_taua_max_taua_l[ib] = 0.;
        uncertainty->derv_chl_rrs        [ib] = 0.;

        for (inir = 0; inir < nbands_ac; inir++) {
            ipb=ib*nbands_ac+inir;

            uncertainty->derv_La_rhorc[ipb] = 0.;
            uncertainty->derv_taua_rhorc[ipb] = 0.;
            uncertainty->derv_tsen_rhorc[ipb] = 0.;
            uncertainty->derv_tsol_rhorc[ipb] = 0.;

            uncertainty->derv_La_taua_l[ipb]   = 0.0;
            uncertainty->derv_taua_taua_l[ipb] = 0.0;

            uncertainty->derv_tsen_taua[ipb] = 0.0;
            uncertainty->derv_tsol_taua[ipb] = 0.0;

            uncertainty->derv_rhownir_rrs[ipb] = 0.;

            uncertainty->derv_La_rhow  [ipb] = 0.0;
            uncertainty->derv_taua_rhow[ipb] = 0.0;
            uncertainty->derv_tsen_rhow[ipb] = 0.0;
            uncertainty->derv_tsol_rhow[ipb] = 0.0;
        }

       
        uncertainty->derv_taua_min_rhow_l[ib] = 0.0;
        uncertainty->derv_taua_max_rhow_l[ib] = 0.0;

        uncertainty->derv_La_rh[ib] = 0.0;
        uncertainty->derv_taua_rh[ib] = 0.0;
        uncertainty->derv_tsen_rh[ib] = 0.0;
        uncertainty->derv_tsol_rh[ib] = 0.0;

        uncertainty->drhown_nir[ib] = 0.0;
        uncertainty->dbrdf[ib] = 0.0;

        uncertainty->derv_rhog_taua[ib]=0.0;
        uncertainty->derv_rhog_rhorc_l[ib]=0.0;
        uncertainty->derv_TLg_gc[ib]   =0.;

         for (inir = 0; inir < nbands; inir++)
            uncertainty->covariance_matrix[ib*nbands+inir] = 0.;

    }
       
}

void free_uncertainty(uncertainty_t *uncertainty) {
    /*for(i=0;i<uncertainty->nbands;i++){
        free(uncertainty->derv_La_rhorc[i]);
        free(uncertainty->derv_taua_rhorc[i]);
        free(uncertainty->derv_tsen_rhorc[i]);
        free(uncertainty->derv_tsol_rhorc[i]);
    }
    free(uncertainty->derv_La_rhorc);
    free(uncertainty->derv_taua_rhorc);
    free(uncertainty->derv_tsen_rhorc);
    free(uncertainty->derv_tsol_rhorc);

    free(uncertainty->derv_modrat_rhorc);
    free(uncertainty->ratio_rhow);
    free(uncertainty->corr_coef_rhot);*/

    free((void *) uncertainty->data);

}
