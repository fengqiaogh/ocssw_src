#include "l12_proto.h"

/* ============================================================================

c       This corrects the sunglint contamination reflectances in the
c       SeaWiFS 8 bands using the Cox & Munk model for the ocean
c       sun glitter radiance distribution as function of the sea surface
c       wind.
c
c Inputs:
c       num_iter, I, --- iteration number in the atmospheric corrections.
c       nband, I, --- number of bands for the sensor (e.g., 8 for SeaWiFS).
c       glint_coef, R, --- glitter radiance (F0=1) from Cox & Munk.
c       air_mass, R, --- airmass value.
c       mu0, R, --- cosine of the solar zenith angle.
c       taur(nband), R, --- Rayleigh optical thicknesses.
c       taua(nband), R, --- retrieved aerosol optical thicknesses.
c       La(nband), R, --- aerosol reflectances at 8 SeaWiFS bands.
c Output:
c       TLg(nband), R, --- sunglint radiances at 8 SeaWiFS bands.
c Credits:
c       Written by
c        Menghua Wang
c        UMBC, Code 970.2, NASA/GSFC, Greenbelt, MD
c        10-7-1999
c
c       Modification History
c       - Modified to return glint radiance at TOA rather than to 
c         correct an input reflectance. Bryan Franz, 14 October 1999.
c       - Added low rhoa enhancements from M. Wang, 28 November 1999.
c	- Added a check and necessary correction to make sure there is 
c	  no over-correction, i.e., there is no pixel lost due to sun 
c	  glint corrections.  Menghua Wang, 2-18-00
c ============================================================================= */

static float taua_est(float x) {
    return (-0.8 - 0.4 * log(x));
};

void glint_rad(int32_t iter, int32_t nband, int32_t nir_s, int32_t nir_l,
        float *glint_coef, float airmass,
        float mu0, float F0[], float taur[], float taua[], float La[], float TLg[], uncertainty_t *uncertainty) {
    static float glint_min = GLINT_MIN;
    static float taua_min = 0.08;
    static float taua_ave = 0.1;
    static float rhoa_min = 0.01;
    static float rhoa_min2 = 0.008;
    static int32_t iter_max = 2;
    static float rfac = 0.8;
    float *derv_rhog_taua=NULL, *derv_TLg_gc=NULL,*derv_rhog_rhorc_l;
    int aer_base=bindex_get(input->aer_wave_base);
    float temp_derv;
    float derv_taua_ave_gc=0.,derv_taua_ave_rhorc_l=0.;

    if (uncertainty) {
        derv_rhog_taua = uncertainty->derv_rhog_taua;
        derv_rhog_rhorc_l = uncertainty->derv_rhog_rhorc_l;
        derv_TLg_gc    = uncertainty->derv_TLg_gc;
    }

    int32_t ib;
    float taua_c,multiply_taua_c;
    float refl_test;
    float taua_ave2;
    float fac;
    float derv_taua_c;


    /* If the number of iterations exceeds the maximum, we just return.  This assumes   */
    /* that the calling routine saved the TLg from the previous iteration. If the glint */
    /* coefficient is very low, return 0.                                               */

    if (iter > iter_max)
        return;

    else if (glint_coef[aer_base] <= glint_min)
        for (ib = 0; ib < nband; ib++) {
            TLg[ib] = 0.0;
            if (uncertainty){
                derv_rhog_taua[ib] = 0.0;
                derv_rhog_rhorc_l[ib] = 0.0;
            }              
        } else {

        refl_test = OEL_PI / mu0 * (La[nir_l] / F0[nir_l] - glint_coef[nir_l] * exp(-(taur[nir_l] + taua_ave) * airmass));

        if (refl_test <= rhoa_min){
            taua_ave2 = taua_est(10. * MAX(refl_test, 0.0001));  // !!!!! This part may need to include rhorc[nir_l]
            if(uncertainty){
                if (refl_test > 0.0001) {
                    temp_derv=-0.4 / (10 * refl_test)* OEL_PI / mu0;

                    derv_taua_ave_rhorc_l= 1/F0[nir_l];
                    derv_taua_ave_gc =-exp(-(taur[nir_l] + taua_ave) * airmass);

                    derv_taua_ave_gc *= temp_derv;
                    derv_taua_ave_rhorc_l*=temp_derv;
                }
            }
        } else {
            taua_ave2 = taua_ave;
        }

        for (ib = 0; ib < nband; ib++) {

            if(uncertainty){
                derv_rhog_taua[ib] = 0.0;
                temp_derv=0.;
            }
                
            if (iter <= 1) {
            	taua_c = taua_ave2;
            	derv_taua_c=0.0;
                temp_derv=derv_taua_ave_gc;

            } else if (taua[nir_l] <= taua_min) {
                derv_taua_ave_rhorc_l=0.;
            	taua_c = taua_est(taua[nir_l]);
                
            	if(uncertainty) {
                    derv_taua_c=-0.4/taua[nir_l];
                    uncertainty->aer_l_glint=1;
            	}
            } else {
                derv_taua_ave_rhorc_l=0.;
            	taua_c = taua[ib];
            	derv_taua_c=1;
            }

            /* Check for over-correction */
            if (ib == 0)
                refl_test = OEL_PI / mu0 * (La[nir_l] / F0[nir_l] - glint_coef[nir_l] * exp(-(taur[nir_l] + taua_c) * airmass));

            if (refl_test <= rhoa_min2)
                multiply_taua_c = 1.5;
            else
                multiply_taua_c = 1.0;

            TLg[ib] = F0[ib] * glint_coef[ib] * exp(-(taur[ib] + multiply_taua_c * taua_c) * airmass);
            if (uncertainty) {
                derv_rhog_taua[ib]    = -TLg[ib] * multiply_taua_c * airmass * derv_taua_c;
                derv_rhog_rhorc_l[ib] = -TLg[ib] * multiply_taua_c * airmass * derv_taua_ave_rhorc_l;
                derv_TLg_gc[ib] = temp_derv * (-multiply_taua_c * TLg[ib] * airmass);
                derv_TLg_gc[ib] += F0[ib] * exp(-(taur[ib] + multiply_taua_c * taua_c) * airmass);
            }
        }

        /* Make sure there is no over-correction */
        if (La[nir_l] > 0.0 && La[nir_s] > 0.0) {
            fac = MAX(TLg[nir_l] / La[nir_l], TLg[nir_s] / La[nir_s]);
            if (fac >= rfac) {
                for (ib = 0; ib < nband; ib++) {
                    TLg[ib] = rfac * TLg[ib] / fac;
                    if(uncertainty){
                        derv_rhog_taua[ib]*=(rfac/fac);
                        derv_TLg_gc[ib]*=(rfac/fac);
                        derv_rhog_rhorc_l[ib]*=(rfac/fac);
                    }
                }

            }

        }

    }
}






