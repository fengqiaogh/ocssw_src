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
    int32_t i;

    uncertainty->nbands=nbands;
    uncertainty->nbands_ac=nbands_ac;
    uncertainty->npix=npix;

    /*                                                      */
    /* allocate data block as contiguous bytes              */
    /*                                                      */
    len =  11 * sizeof (float)*npix
         + 10 * sizeof (float)*npix * nbands
         + 21 * sizeof (float)* nbands
         + 1 * sizeof (float)* nbands*nbands;


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

    uncertainty->dsensor =(float *) p;  p+=sizeof(float)*npix*nbands;
    uncertainty->dsyst   =(float *) p;  p+=sizeof(float)*npix*nbands;
    uncertainty->dvc     =(float *) p;  p+=sizeof(float)*npix*nbands;

    uncertainty->dwv     =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->dmw     =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->dzw     =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->dws     =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->dwd     =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->doz     =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->dpr     =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->drh     =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->dmodrat =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->dno2_tropo =(float *) p;  p+=sizeof(float)*npix;
    uncertainty->dno2_strat =(float *) p;  p+=sizeof(float)*npix;


    uncertainty->dLr     =(float *) p;  p+=sizeof(float)*npix*nbands;
    uncertainty->dtg_sol =(float *) p;  p+=sizeof(float)*npix*nbands;
    uncertainty->dtg_sen =(float *) p;  p+=sizeof(float)*npix*nbands;
    uncertainty->dt_sol  =(float *) p;  p+=sizeof(float)*npix*nbands;
    uncertainty->dt_sen  =(float *) p;  p+=sizeof(float)*npix*nbands;
    uncertainty->dTLg    =(float *) p;  p+=sizeof(float)*npix*nbands;

    uncertainty->derv_pol =(float *) p;  p+=sizeof(float)*npix*nbands;

    uncertainty->derv_La_taua_l=(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_La_rhow_l=(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_La_rh    =(float  *)p; p+=sizeof(float)*nbands;

    uncertainty->derv_taua_taua_l=(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_taua_rhow_l=(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_taua_rh    =(float  *)p; p+=sizeof(float)*nbands;

    uncertainty->derv_tsen_taua_l=(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_tsen_rhow_l=(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_tsen_rh    =(float  *)p; p+=sizeof(float)*nbands;

    uncertainty->derv_tsol_taua_l=(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_tsol_rhow_l=(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_tsol_rh    =(float  *)p; p+=sizeof(float)*nbands;

    uncertainty->derv_taua_min_rhorc_l=(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_taua_min_taua_l=(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_taua_min_rhow_l=(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_taua_max_rhorc_l=(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_taua_max_taua_l=(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->derv_taua_max_rhow_l=(float  *)p; p+=sizeof(float)*nbands;


    uncertainty->derv_Lg_taua=(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->drhown_nir   =(float  *)p; p+=sizeof(float)*nbands;
    uncertainty->dbrdf        =(float  *)p; p+=sizeof(float)*nbands;

    uncertainty->pixel_covariance  =(float * )p; p+=sizeof(float)*nbands*nbands;


    if ((len - (int32_t) (p - uncertainty->data)) < 0) {
        printf("%s Line %d: bad allocation on error record\n", __FILE__, __LINE__);
        exit(1);
    }
    if (want_verbose)
        printf("Allocated %d bytes in error record.\n", (int) (p - uncertainty->data));

    uncertainty->derv_La_rhorc  =(float **)malloc(nbands*sizeof(float *));
    uncertainty->derv_taua_rhorc=(float **)malloc(nbands*sizeof(float *));
    uncertainty->derv_tsen_rhorc=(float **)malloc(nbands*sizeof(float *));
    uncertainty->derv_tsol_rhorc=(float **)malloc(nbands*sizeof(float *));

    for(i=0;i<nbands;i++){
        uncertainty->derv_La_rhorc[i]=(float *)malloc(nbands_ac*sizeof(float));
        uncertainty->derv_taua_rhorc[i]=(float *)malloc(nbands_ac*sizeof(float));
        uncertainty->derv_tsen_rhorc[i]=(float *)malloc(nbands_ac*sizeof(float));
        uncertainty->derv_tsol_rhorc[i]=(float *)malloc(nbands_ac*sizeof(float));
    }
    uncertainty->derv_modrat_rhorc=(float *)malloc(nbands_ac*sizeof(float));
    uncertainty->ratio_rhow       =(float *)malloc(nbands_ac*sizeof(float));

    uncertainty->corr_coef_rhot   =(float *)malloc(nbands*nbands*sizeof(float));

    return (0);
}

int cp_uncertainty(uncertainty_t *oldrec, uncertainty_t *newrec, int32_t ip)
{
    int32_t nbands=oldrec->nbands;
    int32_t nbands_ac=oldrec->nbands_ac;

    int32_t ipb = ip * nbands;
    memcpy(&newrec->dsensor[ipb], &oldrec->dsensor[ipb], sizeof(float)*nbands);
    memcpy(&newrec->dsyst[ipb], &oldrec->dsyst[ipb], sizeof(float)*nbands);
    memcpy(&newrec->dLr[ipb], &oldrec->dLr[ipb], sizeof(float)*nbands);
    memcpy(&newrec->dtg_sol[ipb], &oldrec->dtg_sol[ipb], sizeof(float)*nbands);
    memcpy(&newrec->dtg_sen[ipb], &oldrec->dtg_sen[ipb], sizeof(float)*nbands);
    memcpy(&newrec->dt_sol[ipb], &oldrec->dt_sol[ipb], sizeof(float)*nbands);
    memcpy(&newrec->dt_sen[ipb], &oldrec->dt_sen[ipb], sizeof(float)*nbands);
    memcpy(&newrec->dTLg[ipb], &oldrec->dTLg[ipb], sizeof(float)*nbands);

    memcpy(newrec->ratio_rhow, oldrec->ratio_rhow, sizeof(float)*nbands_ac);

    newrec->dwv[ip] = oldrec->dwv[ip];
    newrec->dws[ip] = oldrec->dws[ip];
    newrec->doz[ip] = oldrec->doz[ip];
    newrec->dpr[ip] = oldrec->dpr[ip];
    newrec->drh[ip] = oldrec->drh[ip];

    newrec->derv_eps_Lrc_l=oldrec->derv_eps_Lrc_l;
    newrec->derv_eps_Lrc_s=oldrec->derv_eps_Lrc_s;
    newrec->derv_eps_taua_l=oldrec->derv_eps_taua_l;
    newrec->derv_eps_taua_s=oldrec->derv_eps_taua_s;
    newrec->derv_eps_rhow_l=oldrec->derv_eps_rhow_l;

    return 0;
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

    if(ifscan){

        for(ip=0;ip<npix;ip++){
            uncertainty->dwv[ip]=0.;
            uncertainty->dmw[ip]=0.;
            uncertainty->dzw[ip]=0.;
            uncertainty->dws[ip]=0.;
            uncertainty->dwd[ip]=0.;
            uncertainty->doz[ip]=0.;
            uncertainty->dpr[ip]=0.;
            uncertainty->drh[ip]=0.;
            uncertainty->dmodrat[ip]=0.;
            uncertainty->dno2_tropo[ip]=0.;
            uncertainty->dno2_strat[ip]=0.;

            for(ib=0;ib<nbands;ib++){
                ipb=ip*nbands+ib;
                uncertainty->dsensor[ipb]=0.0;
                uncertainty->dsyst  [ipb]=0.0;
                uncertainty->dvc    [ipb]=0.0;

                uncertainty->dLr    [ipb]=0.0;
                uncertainty->dtg_sol[ipb]=0.0;
                uncertainty->dtg_sen [ipb]=0.0;
                uncertainty->dt_sol  [ipb]=0.0;
                uncertainty->dt_sen  [ipb]=0.0;
                uncertainty->dTLg    [ipb]=0.0;
                uncertainty->derv_pol[ipb]=0.0;
            }
        }
    }
    else{

        uncertainty->derv_modrat_taua_l=0.;
        uncertainty->derv_eps_Lrc_s = 0.;
        uncertainty->derv_eps_Lrc_l = 0.;
        uncertainty->derv_eps_taua_l = 0;
        uncertainty->derv_eps_rhow_l = 0;
        uncertainty->derv_eps_taua_s = 0;
        uncertainty->dchl=0;

        for(inir=0;inir<nbands_ac;inir++){
            uncertainty->derv_modrat_rhorc[inir]=0.;
            uncertainty->ratio_rhow[inir]=0.;
        }

        for(ib=0;ib<nbands;ib++){
            uncertainty->derv_Lg_taua[ib]=0.;

            uncertainty->derv_taua_min_rhorc_l[ib]=0.;
            uncertainty->derv_taua_max_rhorc_l[ib]=0.;
            uncertainty->derv_taua_min_taua_l[ib]=0.;
            uncertainty->derv_taua_max_taua_l[ib]=0.;

            uncertainty->derv_La_taua_l[ib]=0.0;
            uncertainty->derv_taua_taua_l[ib]=0.0;
            uncertainty->derv_tsen_taua_l[ib]=0.0;
            uncertainty->derv_tsol_taua_l[ib]=0.0;

            for(inir=0;inir<nbands_ac;inir++){
                uncertainty->derv_La_rhorc[ib][inir]=0.;
                uncertainty->derv_taua_rhorc[ib][inir]=0.;
                uncertainty->derv_tsen_rhorc[ib][inir]=0.;
                uncertainty->derv_tsol_rhorc[ib][inir]=0.;
            }

            uncertainty->derv_La_rhow_l[ib]=0.0;
            uncertainty->derv_taua_rhow_l[ib]=0.0;
            uncertainty->derv_tsen_rhow_l[ib]=0.0;
            uncertainty->derv_tsol_rhow_l[ib]=0.0;
            uncertainty->derv_taua_min_rhow_l [ib]=0.0;
            uncertainty->derv_taua_max_rhow_l [ib]=0.0;

            uncertainty->derv_La_rh[ib]=0.0;
            uncertainty->derv_taua_rh[ib]=0.0;
            uncertainty->derv_tsen_rh[ib]=0.0;
            uncertainty->derv_tsol_rh[ib]=0.0;

            uncertainty->drhown_nir [ib]=0.0;
            uncertainty->dbrdf  [ib]=0.0;
        }
    }
    for(ib=0;ib<nbands*nbands;ib++)
        uncertainty->pixel_covariance[ib]=0.;
}

void free_uncertainty(uncertainty_t *uncertainty) {
    int i;

    for(i=0;i<uncertainty->nbands;i++){
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
    free(uncertainty->corr_coef_rhot);

    free((void *) uncertainty->data);

}
