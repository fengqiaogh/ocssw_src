/*
 * uncertainty.h
 *
 *  Created on: Mar 8, 2021
 *      Author: mzhang11
 */

#ifndef OEL_HDF4_LIBL1_UNCERTAINTY_H_
#define OEL_HDF4_LIBL1_UNCERTAINTY_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct uncertainty_struct{

    char *data;
    int32_t length;

    int32_t nbands;
    int32_t nbands_ac;
    int32_t npix;
   
    /* uncertainty attributes*/
    float *dsensor;
    float *dsyst;
    float *dvc;

    float *drh;
    float *doz;
    float *dpr;
    float *dwv;
    float *dmw;
    float *dzw;
    float *dws;
    float *dwd;
    float *dno2_tropo;
    float *dno2_strat;

    /* uncertainty in temporary products during L1->L2 */
    float *dtg_sol;
    float *dtg_sen;
    float *dt_sol;
    float *dt_sen;
    float *dtaua;

    /* derivative of polcor to Lt,tg */
    float *derv_polcor_Lt;
    float *derv_polcor_tgsol;
    float *derv_polcor_tgsen;

    float *derv_La_rhorc;   //derivative of La[wave] to rhorc[nbands_ac]=rhot-rhor  1st dimesion:nwave, 2nd: nbands_ac
    float *derv_La_taua_l;     //derivative of La[wave] to taua[nbands_ac], which is calculated from last iteration
    float *derv_La_rhow;   //derivative of La[wave] to t_sen.t_sol.rhow[nbands_ac]
    float *derv_La_rh ;    //derivative of La[wave] to rh

    float *derv_taua_rhorc;   //derivative of taua[wave] to rhorc[nbands_ac]=rhot-rhor  1st dimesion:nwave, 2nd: nbands_ac
    float *derv_taua_taua_l;   //derivative of taua[wave] to taua[nbands_ac], which is calculated from last iteration
    float *derv_taua_rhow;   //derivative of taua[wave] to t_sen.t_sol.rhow[nbands_ac]
    float *derv_taua_rh ;      //derivative of taua[wave] to rh

    float *derv_tsen_rhorc;   //derivative of tsen[wave] to rhorc[nbands_ac]=rhot-rhor  1st dimesion:nwave, 2nd: nbands_ac
    float *derv_tsen_taua;   //derivative of tsen[wave] to taua[nbands_ac], which is calculated from last iteration
    float *derv_tsen_rhow;   //derivative of tsen[wave] to t_sen.t_sol.rhow[nbands_ac]
    float *derv_tsen_rh ;      //derivative of tsen[wave] to rh

    float *derv_tsol_rhorc;   //derivative of tsol[wave] to rhorc[nbands_ac]=rhot-rhor  1st dimesion:nwave, 2nd: nbands_ac
    float *derv_tsol_taua;   //derivative of tsol[wave] to taua[nbands_ac], which is calculated from last iteration
    float *derv_tsol_rhow;   //derivative of tsol[wave] to t_sen.t_sol.rhow[nbands_ac]
    float *derv_tsol_rh ;      //derivative of tsol[wave] to rh

    float *derv_rhog_taua; //derivative of glint reflectance rhog[wave] to either taua at the corresponding bands or aer_l
    float *derv_rhog_rhorc_l; //derivative of glint reflectance rhog[wave] to rhorc[aer_l]
    float *drhown_nir;    //uncertainty in rhown(NIR)
    float *dbrdf;         //uncertainty in brdf correction

    float derv_eps_rhorc_l;//derivative of observed eps to rhorc[aer_l]
    float derv_eps_rhorc_s;//derivative of observed eps to rhorc[aer_s]
   // float derv_eps_taua_l;//derivative of observed eps to taua[aer_l]
    //float derv_eps_taua_s;//derivative of observed eps to taua[aer_s]
    float *derv_eps_rhow;//derivative of observed eps to t_sen*t_sol*rhow[nbands_ac]
    float *ratio_rhow;  //ratio of tLw[aer_s to aer_l].to tLw[aer_l]

    float *derv_modrat_rhorc;//derivative of modrat to rhorc[nbands_ac]
    float *derv_modrat_taua; //derivative of modrat to taua [nbands_ac], which is calculated from last iteration
    float derv_modrat_taua_l;//derivative of mwt to taua[aer_l], which is calculated from last iteration
    
    float derv_modrat_rhow_l;//derivative of mwt to t_sen*t_sol*rhow[aer_l]
    float dchl;           //uncertainty in chla
    float dkd490;         //uncertainty in kd490

    float *derv_taua_min_rhorc_l;   //derivative of tauamin[nwave] to rhorc[aer_l]=rhot-rhor
    float *derv_taua_min_taua_l;   //derivative of tauamin[nwave] to taua[aer_l], which is calculated from last iteration
    float *derv_taua_min_rhow_l;   //derivative of tauamin[nwave] to t_sen.t_sol.rhow[aer_l]
    float *derv_taua_max_rhorc_l;   //derivative of tauamax[nwave] to rhorc[aer_l]=rhot-rhor
    float *derv_taua_max_taua_l;   //derivative of tauamax[nwave] to taua[aer_l], which is calculated from last iteration
    float *derv_taua_max_rhow_l;   //derivative of tauamax[nwave] to t_sen.t_sol.rhow[aer_l]

    float *dRrs; // just a pointer pointing to Rrs_unc[ip*nbands] in l2rec, don't need to allocate memory
    
    float *covariance_matrix;//// just a pointer to covariance_matrix at current pixel ip, don't need to allocate memory

    float *corr_coef_rhot;//correlation coefficient matrix for rhot (nwave*nwave)

    int *acbands_index;// index for the AC bands
    int aer_l_glint;   //1: taua(aer_l) is used for TLg; 0: taua at corresponding band is used for TLg
    float *derv_gc_ws;  // derivative of glint coefficient to wind speed, dimension npixl
    float *derv_TLg_gc;  // derivative of TLg[nwave] to glint coefficient 

    float *derv_rhownir_rrs; // derivative of rhown_nir[nbands_ac] to rrs[nwave]
    float *derv_rhownir_chl; // derivative of rhown_nir[nbands_ac] to chl
    
    float *derv_chl_rrs;  //derivative of chl to rrs[nwave]

    int nbands_rhownir;  // No. of bands used to calculate rhownir
    int *bindx_rhownir;   // index of the bands used to calculate rhownir

} uncertainty_t;

int alloc_uncertainty(int32_t nbands, int32_t nbands_ac, int32_t npix, uncertainty_t *uncertainty);
void init_uncertainty(uncertainty_t *uncertainty, int ifscan);
void free_uncertainty(uncertainty_t *uncertainty);

#ifdef __cplusplus
}
#endif

#endif /* OEL_HDF4_LIBL1_UNCERTAINTY_H_ */
