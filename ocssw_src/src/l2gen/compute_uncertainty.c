#include "l12_proto.h"
#include "atrem_corl1.h"
#include "glint.h"
#if !defined(MACINTOSH)
__attribute__((target_clones("avx2", "default")))
#endif
void compute_uncertainty(l2str* __restrict__ l2rec, int32_t ip, int aer_l, int rhownir_corr, int want_glintcorr, int proc_uncertainty,
                         float* __restrict__  radref, float* __restrict__  delta_Lt, float* __restrict__  derv_rhownir_chl, float* __restrict__  derv_rhownir_rrs,
                         float* __restrict__  tLw, float* __restrict__  last_derv_taua_rhorc, float ** __restrict__ derv_taua_rhorc, float* __restrict__ brdf, float* __restrict__ covariance_matrix) {
    static int firstcall = 1;
    int32_t ib, ix1, ix2, ix3, iw=0, j, i, ipb,inir,iy;
    static int32_t nwave, nbands_ac, nbands_rhownir, nbands_uncertainty;
    float tmp, tmp_derv1,tmp_derv2, tmp_derv3, cov, * __restrict__ t_sol, *  __restrict__ t_sen,* __restrict__ dt_sen,*  __restrict__ dt_sol, * __restrict__ dtaua;
    static float * __restrict__ corr_coef_rhot;
    static int32_t * __restrict__ acbands_index, * __restrict__ bindex_rhownir, * __restrict__ bindex_uncertainty;
    static float* __restrict__ Fobar;
    l1str * __restrict__ l1rec=l2rec->l1rec;
    uncertainty_t* __restrict__ uncertainty=l1rec->uncertainty;
    
    static float* __restrict__ derv_nlw_rh; //derivative of nlw[nwave] to rh
    static float** __restrict__ derv_tsen_rhorc;
    static float** __restrict__ derv_tsol_rhorc;
    static float** __restrict__ derv_nlw_rhorc;
    static float* __restrict__ derv_nlw_chl;
    static float * __restrict__ derv_nlw_rhorc_corrected;
    static float * covariance_matrix_test;
    static float *__restrict__  corr_coef_rhot_corrected;
    static float *__restrict__  corr_coef_rhot_diag;
    static float *__restrict__  derv_nlw_rhorc_corr_coef_rhot_product;
    static float *__restrict__  derv_nlw_rhorc_rrs_cov_product;
    static float *__restrict__  derv_nlw_rhorc_rrs_product;
    double fsol=l1rec->fsol;
    double * __restrict__ tg;
    float mu0=cos(l1rec->solz[ip]/OEL_RADEG);
    int nbands_covariance=0; // No. of bands included in calclating the covariance matrix
    static int *bindex_covaraince; // index of bands included in calclating the covariance matrix

    if (firstcall) {
        firstcall = 0;
        Fobar = l2rec->l1rec->l1file->Fobar;
        nwave = l2rec->l1rec->l1file->nbands;
        nbands_ac = input->nbands_ac;
        nbands_uncertainty = nbands_ac + 1;
        nbands_rhownir = uncertainty->nbands_rhownir;
        bindex_rhownir = uncertainty->bindx_rhownir;

        if (proc_uncertainty == 1) {
            bindex_covaraince = (int*)malloc((nbands_rhownir + 1) * sizeof(int));
            for (ib = 1; ib < nbands_rhownir + 1; ib++)
                bindex_covaraince[ib] = bindex_rhownir[ib - 1];
        } 

        acbands_index = input->acbands_index;
        corr_coef_rhot = uncertainty->corr_coef_rhot;

        bindex_uncertainty = (int32_t*)malloc(nbands_uncertainty * sizeof(int32_t));
        for (ib = 0; ib < nbands_ac; ib++)
            bindex_uncertainty[ib + 1] = acbands_index[ib];

        if ((derv_nlw_rh = (float *)malloc(nwave*sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to derv_rh\n");
            exit(FATAL_ERROR);
        }
        if ((derv_nlw_chl = (float *)malloc(nwave*sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to derv_nlw_chl\n");
            exit(FATAL_ERROR);
        }
        if ((derv_nlw_rhorc = (float **)calloc(nwave, sizeof(float *))) == NULL) {
            printf("-E- : Error allocating memory to derv_nlw_rhorc\n");
            exit(FATAL_ERROR);
        }
        if ((derv_tsol_rhorc = (float **)malloc(nwave*sizeof(float *))) == NULL) {
            printf("-E- : Error allocating memory to derv_tsol_rhorc\n");
            exit(FATAL_ERROR);
        }
        if ((derv_tsen_rhorc = (float **)malloc(nwave*sizeof(float *))) == NULL) {
            printf("-E- : Error allocating memory to derv_tsen_rhorc\n");
            exit(FATAL_ERROR);
        }
        for (ib = 0; ib < nwave; ib++) {
            if ((derv_nlw_rhorc[ib] = (float *)calloc(nwave, sizeof(float))) == NULL) {
                printf("-E- : Error allocating memory to derv_nlw_rhorc[%d]\n", ib);
                exit(FATAL_ERROR);
            }
            if ((derv_tsen_rhorc[ib] = (float *)malloc(nwave*sizeof(float))) == NULL) {
                printf("-E- : Error allocating memory to derv_tsen_rhorc[%d]\n", ib);
                exit(FATAL_ERROR);
            }
            if ((derv_tsol_rhorc[ib] = (float *)malloc(nwave*sizeof(float))) == NULL) {
                printf("-E- : Error allocating memory to derv_tsol_rhorc[%d]\n", ib);
                exit(FATAL_ERROR);
            }
        }
       // preallocate
        derv_nlw_rhorc_corrected =  (float *) malloc(nwave * nbands_uncertainty *sizeof(float));
        corr_coef_rhot_corrected =  (float *) malloc(nbands_uncertainty * nbands_uncertainty *sizeof(float));
        corr_coef_rhot_diag =(float *) malloc(nwave * nbands_uncertainty *sizeof(float));
        derv_nlw_rhorc_corr_coef_rhot_product =  (float *) malloc(nwave * nbands_uncertainty *sizeof(float));
        covariance_matrix_test = (float *) malloc(nwave * nwave*sizeof(float));
        derv_nlw_rhorc_rrs_product = (float *) malloc(nwave * nbands_rhownir *sizeof(float));
        derv_nlw_rhorc_rrs_cov_product = (float *) malloc(nwave * nbands_rhownir *sizeof(float));
    }

    ipb=ip*nwave;
    t_sol = &l2rec->l1rec->t_sol[ipb];
    t_sen = &l2rec->l1rec->t_sen[ipb];
    dt_sen= &uncertainty->dt_sen[ipb];
    dt_sol= &uncertainty->dt_sol[ipb];
    dtaua= &uncertainty->dtaua[ipb];
    tg=&l1rec->tg[ipb];

    if (proc_uncertainty == 1)
        nbands_covariance = nbands_rhownir + 1;  // the additional 1 means the band itself;

    /* calculate the derivative of t_sol,t_sen, taua, nLw*/

    ix2 = nbands_rhownir;
    if (proc_uncertainty == 2)
        ix2 = nwave;

    /* only calculate the covariance at the bands used to calculate chl within the iteration to account for rhownir*/
    for (ix3=0; ix3 < ix2; ix3++) {
        ib=bindex_rhownir[ix3];
        if (proc_uncertainty == 2)
            ib = ix3;
        derv_nlw_chl[ib] = 0.;
        bindex_uncertainty[0] = ib;
        for (ix1 = 0; ix1 < nbands_uncertainty; ix1++) {
            iw = bindex_uncertainty[ix1];
            derv_nlw_rhorc[ib][iw] = 0.;
            derv_tsen_rhorc[ib][iw] = 0.;
            derv_tsol_rhorc[ib][iw] = 0.;
            derv_taua_rhorc[ib][iw] = 0.;

            /* make sure band corresponding to index 0 is not same as the bands in nbands_rhownir and
             * nbands_ac*/
            if (ix1 != 0 && ib == iw)
                bindex_uncertainty[0] = BAD_INT;
        }

        /* Derivative of nLw[ib] to rhorc[ib]             */
        derv_nlw_rhorc[ib][ib] = 1 / t_sol[ib] / t_sen[ib] / radref[ib];

        /* Derivative of nLw,tsen, tsol, and taua to rhorc[nbands_ac]             */
        for (inir = 0; inir < nbands_ac; inir++) {
            i = acbands_index[inir];
            ipb = ib * nbands_ac + inir;

            /* contribution from La to nLw[ib]   */
            derv_nlw_rhorc[ib][i] = -uncertainty->derv_La_rhorc[ipb] / (t_sol[ib] * t_sen[ib]);

            derv_tsen_rhorc[ib][i] = uncertainty->derv_tsen_rhorc[ipb];
            derv_tsol_rhorc[ib][i] = uncertainty->derv_tsol_rhorc[ipb];
            derv_taua_rhorc[ib][i] = uncertainty->derv_taua_rhorc[ipb];

            if (want_glintcorr) {
                /* contribution through TLg[ib] */
                j = ib;
                if (uncertainty->aer_l_glint)
                    j = aer_l;
                tmp = -1 / (t_sol[ib] * t_sen[ib] * radref[ib]);
                derv_nlw_rhorc[ib][i] +=
                    (tmp * uncertainty->derv_rhog_taua[ib] * last_derv_taua_rhorc[j * nbands_ac + inir]);

                /* contribution from rhorc[aer_l] through rhog[ib] */
                if (i == aer_l && uncertainty->derv_rhog_rhorc_l[ib] != 0.)
                    derv_nlw_rhorc[ib][i] += (tmp * uncertainty->derv_rhog_rhorc_l[ib]);

                /*contribution through rhog[ac_bands]*/
                for (iw = 0; iw < nbands_ac; iw++) {
                    ix1 = acbands_index[iw];
                    j = ix1;
                    if (uncertainty->aer_l_glint)
                        j = aer_l;
                    tmp = uncertainty->derv_rhog_taua[ix1] * last_derv_taua_rhorc[j * nbands_ac + inir];
                    if (tmp != 0) {
                        ix1 = ib * nbands_ac + iw;
                        derv_nlw_rhorc[ib][i] +=
                            (uncertainty->derv_La_rhorc[ix1] / (t_sol[ib] * t_sen[ib]) * tmp);

                        derv_tsen_rhorc[ib][i] += (-uncertainty->derv_tsen_rhorc[ix1] * tmp);
                        derv_tsol_rhorc[ib][i] += (-uncertainty->derv_tsol_rhorc[ix1] * tmp);
                        derv_taua_rhorc[ib][i] += (-uncertainty->derv_taua_rhorc[ix1] * tmp);
                    }
                }

                if (i == aer_l) {
                    for (iw = 0; iw < nbands_ac; iw++) {
                        tmp = uncertainty->derv_rhog_rhorc_l[acbands_index[iw]];
                        if (tmp != 0.) {
                            ix1 = ib * nbands_ac + iw;
                            derv_nlw_rhorc[ib][i] +=
                                (uncertainty->derv_La_rhorc[ix1] / (t_sol[ib] * t_sen[ib]) * tmp);
                            derv_tsen_rhorc[ib][i] += (-uncertainty->derv_tsen_rhorc[ix1] * tmp);
                            derv_tsol_rhorc[ib][i] += (-uncertainty->derv_tsol_rhorc[ix1] * tmp);
                            derv_taua_rhorc[ib][i] += (-uncertainty->derv_taua_rhorc[ix1] * tmp);
                        }
                    }
                }
            }
        }

        /* the contribution from tsol[ib], tsen[ib] to nLw[ib]*/
        tmp_derv1 = -tLw[ib] / (t_sol[ib] * t_sen[ib] * t_sen[ib]);
        tmp_derv2 = -tLw[ib] / (t_sol[ib] * t_sol[ib] * t_sen[ib]);
        for (ix1 = 0; ix1 < nbands_ac; ix1++) {
            i = acbands_index[ix1];
            derv_nlw_rhorc[ib][i] += (tmp_derv1 * derv_tsen_rhorc[ib][i]);
            derv_nlw_rhorc[ib][i] += (tmp_derv2 * derv_tsol_rhorc[ib][i]);
            if (rhownir_corr) {
                derv_nlw_chl[ib] += (tmp_derv1 * derv_tsen_rhorc[ib][i] * derv_rhownir_chl[ix1]);
                derv_nlw_chl[ib] += (tmp_derv2 * derv_tsol_rhorc[ib][i] * derv_rhownir_chl[ix1]);
            }
        }

        /* the derivative to rhorc calculated so far is to rhot/(tg. polcor)                 */
        /* now need to transfer the derivative to rhot which are used in calculating polcor   */

        for (ix1 = 0; ix1 < nbands_uncertainty; ix1++) {
            iw = bindex_uncertainty[ix1];
            if (iw == BAD_INT)
                continue;

            ipb = ip * nwave + iw;

            /* tmp is the derivative of rhot/(tg.polcor) to rhot  */
            tmp = l1rec->polcor[ipb];
            tmp = 1 / (tmp * tg[iw]) -
                  (l1rec->Lt[ipb] / (tmp * tmp * tg[iw]) * uncertainty->derv_polcor_Lt[ipb]);

            derv_nlw_rhorc[ib][iw] *= tmp;
            derv_tsol_rhorc[ib][iw] *= tmp;
            derv_tsen_rhorc[ib][iw] *= tmp;
            derv_taua_rhorc[ib][iw] *= tmp;

            derv_nlw_rhorc[ib][iw] *= (brdf[ib] / mu0 / fsol);
        }

        if (rhownir_corr) {
            for (inir = 0; inir < nbands_ac; inir++) {
                ipb = ib * nbands_ac + inir;
                derv_nlw_chl[ib] +=
                    uncertainty->derv_La_rhorc[ipb] / (t_sol[ib] * t_sen[ib]) * derv_rhownir_chl[inir];
            }
            derv_nlw_chl[ib] *= (brdf[ib] / mu0 / fsol);
        }

        /* Derivative of nLw to rh             */
        derv_nlw_rh[ib] = -uncertainty->derv_La_rh[ib] / (t_sol[ib] * t_sen[ib]);
        derv_nlw_rh[ib] += (-tLw[ib] / (t_sol[ib] * t_sen[ib] * t_sen[ib]) * uncertainty->derv_tsen_rh[ib]);
        derv_nlw_rh[ib] += (-tLw[ib] / (t_sol[ib] * t_sol[ib] * t_sen[ib]) * uncertainty->derv_tsol_rh[ib]);
    }

    /* calculate dt_sen, dt_sol, dtaua*/
    if (proc_uncertainty == 2) {
        for (ib = 0; ib < nwave; ib++) {
            bindex_uncertainty[0] = ib;
            dt_sen[ib] = 0.;
            dt_sol[ib] = 0.;
            dtaua[ib] = 0.;

            for (ix1 = 0; ix1 < nbands_uncertainty; ix1++) {
                iw = bindex_uncertainty[ix1];
                if (iw == BAD_INT || derv_tsen_rhorc[ib][iw] == 0.)
                    continue;

                tmp = radref[iw] * delta_Lt[iw];
                dt_sen[ib] += pow(derv_tsen_rhorc[ib][iw] * tmp, 2);
                dt_sol[ib] += pow(derv_tsol_rhorc[ib][iw] * tmp, 2);
                dtaua[ib] += pow(derv_taua_rhorc[ib][iw] * tmp, 2);

                for (ix2 = ix1 + 1; ix2 < nbands_uncertainty; ix2++) {
                    j = bindex_uncertainty[ix2];
                    if (derv_tsen_rhorc[ib][j] == 0.)
                        continue;

                    tmp = 2 * corr_coef_rhot[iw * nwave + j] * radref[iw] * delta_Lt[iw] * radref[j] *
                          delta_Lt[j];
                    dt_sen[ib] += tmp * derv_tsen_rhorc[ib][iw] * derv_tsen_rhorc[ib][j];
                    dt_sol[ib] += tmp * derv_tsol_rhorc[ib][iw] * derv_tsol_rhorc[ib][j];
                    dtaua[ib] += tmp * derv_taua_rhorc[ib][iw] * derv_taua_rhorc[ib][j];
                }
            }

            dt_sen[ib] += pow(uncertainty->derv_tsen_rh[ib] * uncertainty->drh[ip], 2);
            dt_sol[ib] += pow(uncertainty->derv_tsol_rh[ib] * uncertainty->drh[ip], 2);
            dtaua[ib] += pow(uncertainty->derv_taua_rh[ib] * uncertainty->drh[ip], 2);

            if (rhownir_corr) {
                /* contribution from chl through rhownir */
                tmp_derv1 = 0.;
                tmp_derv2 = 0.;
                tmp_derv3 = 0.;
                for (inir = 0; inir < nbands_ac; inir++) {
                    iw = acbands_index[inir];

                    tmp_derv1 -= derv_tsen_rhorc[ib][iw] * derv_rhownir_chl[inir];
                    tmp_derv2 -= derv_tsol_rhorc[ib][iw] * derv_rhownir_chl[inir];
                    tmp_derv3 -= derv_taua_rhorc[ib][iw] * derv_rhownir_chl[inir];
                }
                tmp = uncertainty->dchl * uncertainty->dchl;
                dt_sen[ib] += tmp_derv1 * tmp_derv1 * tmp;
                dt_sol[ib] += tmp_derv2 * tmp_derv2 * tmp;
                dtaua[ib] += tmp_derv3 * tmp_derv3 * tmp;

                float tmp_tsen, tmp_tsol, tmp_taua;
                /* contribution from Rrs through rhownir */
                for (j = 0; j < nbands_rhownir; j++) {
                    ix1 = bindex_rhownir[j];
                    tmp_tsen = 0.;
                    tmp_tsol = 0.;
                    tmp_taua = 0.;
                    for (i = 0; i < nbands_ac; i++) {
                        tmp_tsen -= derv_tsen_rhorc[ib][acbands_index[i]] * derv_rhownir_rrs[i * nwave + ix1];
                        tmp_tsol -= derv_tsol_rhorc[ib][acbands_index[i]] * derv_rhownir_rrs[i * nwave + ix1];
                        tmp_taua -= derv_taua_rhorc[ib][acbands_index[i]] * derv_rhownir_rrs[i * nwave + ix1];
                    }

                    if (tmp_tsen == 0.)
                        continue;

                    for (i = 0; i < nbands_rhownir; i++) {
                        ix2 = bindex_rhownir[i];

                        tmp_derv1 = 0.;
                        tmp_derv2 = 0.;
                        tmp_derv3 = 0.;
                        for (int k = 0; k < nbands_ac; k++) {
                            tmp_derv1 -=
                                derv_tsen_rhorc[iw][acbands_index[k]] * derv_rhownir_rrs[k * nwave + ix2];
                            tmp_derv2 -=
                                derv_tsol_rhorc[iw][acbands_index[k]] * derv_rhownir_rrs[k * nwave + ix2];
                            tmp_derv3 -=
                                derv_taua_rhorc[iw][acbands_index[k]] * derv_rhownir_rrs[k * nwave + ix2];
                        }
                        if (tmp_derv1 == 0.)
                            continue;

                        if (ix1 <= ix2)
                            tmp = covariance_matrix[ix1 * nwave + ix2];
                        else
                            tmp = covariance_matrix[ix2 * nwave + ix1];

                        dt_sen[ib] += tmp * tmp_derv1 * tmp_tsen;
                        dt_sol[ib] += tmp * tmp_derv2 * tmp_tsol;
                        dtaua[ib] += tmp * tmp_derv3 * tmp_taua;
                    }
                }
            }

            if (dt_sen[ib] < 0)
                dt_sen[ib] = 0.;
            dt_sen[ib] = sqrt(dt_sen[ib]);

            if (dt_sol[ib] < 0)
                dt_sol[ib] = 0.;
            dt_sol[ib] = sqrt(dt_sol[ib]);

            if (dtaua[ib] < 0)
                dtaua[ib] = 0.;
            dtaua[ib] = sqrt(dtaua[ib]);
        }
    }
    // precompute intermidiate matrixes
    for (ib = 0; ib < nwave; ib++) {
        bindex_uncertainty[0] = ib;
        for ( ix1 = 0; ix1 < nbands_uncertainty; ix1++) {
            ix2 = bindex_uncertainty[ix1];
            derv_nlw_rhorc_corrected[ib * nbands_uncertainty + ix1] = derv_nlw_rhorc[ib][ix2];
        }
        for (ix1 = 1; ix1 < nbands_uncertainty; ix1++) {
            ix2 = bindex_uncertainty[ix1];
            if (ix2 == ib) {
                derv_nlw_rhorc_corrected[ib * nbands_uncertainty] = 0;
                break;
            }
        }
        for (ix1 = 0; ix1 < nbands_uncertainty; ix1++) {
            ix2 = bindex_uncertainty[ix1];
            derv_nlw_rhorc_corrected[ib * nbands_uncertainty + ix1] *= delta_Lt[ix2] * radref[ix2];
        }
    }
   
    for ( ix1 = 0; ix1 < nbands_uncertainty; ix1++) {
        for (iy = 0; iy < nbands_uncertainty; iy++) {
            int ind1 = bindex_uncertainty[ix1];
            int ind2 = bindex_uncertainty[iy];
            if (ix1 == 0 || iy == 0) {
                corr_coef_rhot_corrected[ix1 * nbands_uncertainty + iy] = 0;
                continue;
            }
            corr_coef_rhot_corrected[ix1 * nbands_uncertainty + iy] = corr_coef_rhot[ind1 * nwave + ind2];
        }
    }

    for (ib = 0; ib < nwave; ib++) {
        corr_coef_rhot_diag[ib * nbands_uncertainty] = 0;
        derv_nlw_rhorc_corr_coef_rhot_product[nbands_uncertainty * ib ] = 0;
        for (ix1 = 1; ix1 < nbands_uncertainty; ix1++) {
            derv_nlw_rhorc_corr_coef_rhot_product[nbands_uncertainty * ib + ix1] = 0;
            ix2 = bindex_uncertainty[ix1];
            corr_coef_rhot_diag[ib * nbands_uncertainty + ix1] = corr_coef_rhot[ib * nwave + ix2];
        }
    }

    // compute the first  product
    for (ib = 0; ib < nwave; ib++) {
        for ( iy = 1; iy < nbands_uncertainty; iy++) {
            for ( ix1 = 1; ix1 < nbands_uncertainty; ix1++) {
                derv_nlw_rhorc_corr_coef_rhot_product[nbands_uncertainty * ib + ix1] +=
                    derv_nlw_rhorc_corrected[ib * nbands_uncertainty + iy] *
                    corr_coef_rhot_corrected[iy * nbands_uncertainty + ix1];
            }
        }
    }

    for (ib = 0; ib < nwave; ib++) {

        if (proc_uncertainty == 1) {
            bindex_covaraince[0] = ib;
        } else //if (proc_uncertainty == 2)
            nbands_covariance = nwave - ib;

        for (ix1 = 0; ix1 < nbands_covariance; ix1++){
           if (proc_uncertainty == 1)
                iw = bindex_covaraince[ix1];
            else if (proc_uncertainty == 2)
                iw = ib + ix1;
            covariance_matrix_test[ib * nwave + iw] = 0;
        }
    }

    static float *__restrict__  temp_array = NULL;
    if (temp_array == NULL) {
        temp_array = (float *)malloc(sizeof(float) * nwave * nbands_uncertainty);
    }
    // transpose derv_nlw_rhorc_corr_coef_rhot_product and derv_nlw_rhorc_corrected using temp_array
    for (ib = 0; ib < nwave; ib++) {
        for (ix1 = 0; ix1 < nbands_uncertainty; ix1++) {
            temp_array[ix1 * nwave + ib] = derv_nlw_rhorc_corr_coef_rhot_product[nbands_uncertainty * ib + ix1];
        }
    }
    // copy temp_array to derv_nlw_rhorc_corr_coef_rhot_product
    memcpy(derv_nlw_rhorc_corr_coef_rhot_product, temp_array, sizeof(float) * nwave * nbands_uncertainty);

    for ( ib = 0; ib < nwave; ib++) {
        for (ix1 = 0; ix1 < nbands_uncertainty; ix1++) {
            temp_array[ix1 * nwave + ib] = derv_nlw_rhorc_corrected[nbands_uncertainty * ib + ix1];
        }
    }
    memcpy(derv_nlw_rhorc_corrected, temp_array, sizeof(float) * nwave * nbands_uncertainty);


    for (ib = 0; ib < nwave; ib++) {
        for (ix1 = 0; ix1 < nbands_uncertainty; ix1++) {
            temp_array[ix1 * nwave + ib] = corr_coef_rhot_diag[nbands_uncertainty * ib + ix1];
        }
    }
    memcpy(corr_coef_rhot_diag, temp_array, sizeof(float) * nwave * nbands_uncertainty);


    // compute the covariance matrix (the second product)
    for (int ix = 1; ix < nbands_uncertainty; ix++) {
    for (ib = 0; ib < nwave; ib++) {
        if (proc_uncertainty == 1) {
            bindex_covaraince[0] = ib;
            for (ix1 = 1; ix1 < nbands_covariance; ix1++) {
                if (ib == bindex_covaraince[ix1]) {
                    bindex_covaraince[0] = BAD_INT;
                    break;
                }
            }
        } else //if (proc_uncertainty == 2)
            nbands_covariance = nwave - ib;

        for (ix1 = 0; ix1 < nbands_covariance; ix1++) {
            if (proc_uncertainty == 1) {
                iw = bindex_covaraince[ix1];
                if (iw == BAD_INT)
                    continue;
            } else //if (proc_uncertainty == 2)
                iw = ib + ix1;

            covariance_matrix_test[ib * nwave + iw] +=
                derv_nlw_rhorc_corr_coef_rhot_product[ix * nwave + ib] * derv_nlw_rhorc_corrected[ix * nwave + iw];
            }
        }
    }

    for (ix1 = 0; ix1 < nbands_uncertainty; ix1++) {
        for (ib = 0; ib < nwave; ib++) {
            temp_array[ix1 * nwave + ib] = corr_coef_rhot_diag[ix1 * nwave + ib] * derv_nlw_rhorc_corrected[ib];
        }
    }

    // compute the diagonal term
    for (int ix = 0; ix < nbands_uncertainty; ix++) {
        for (ib = 0; ib < nwave; ib++) {
            if (proc_uncertainty == 1) {
                bindex_covaraince[0] = ib;
                for (ix1 = 1; ix1 < nbands_covariance; ix1++) {
                    if (ib == bindex_covaraince[ix1]){
                        bindex_covaraince[0] = BAD_INT;
                        break;
                    }
                }
            } else //if (proc_uncertainty == 2)
                nbands_covariance = nwave - ib;

            for (ix1 = 0; ix1 < nbands_covariance; ix1++) {
                if (proc_uncertainty == 1) {
                    iw = bindex_covaraince[ix1];
                    if (iw == BAD_INT)
                        continue;
                } else //if (proc_uncertainty == 2)
                    iw = ib + ix1;
            
                if (ix == 0) {
                    covariance_matrix_test[ib * nwave + iw] += derv_nlw_rhorc_corrected[ib] *
                                                               derv_nlw_rhorc_corrected[iw] *
                                                               corr_coef_rhot[ib * nwave + iw];
                } else {
                    covariance_matrix_test[ib * nwave + iw] +=
                                                               derv_nlw_rhorc_corrected[iw + ix * nwave] *
                                                               temp_array[ib + ix * nwave];
                    covariance_matrix_test[ib * nwave + iw] +=
                                                               derv_nlw_rhorc_corrected[ib + ix * nwave] *
                                                               temp_array[iw + ix * nwave];
                }
            }
        }
    }

    // precompute derv_nlw_rhorc_rrs_product
    static float * __restrict__  delta_cov_nir_array = NULL;
    if (rhownir_corr) {
        for (ib = 0; ib < nwave; ib++) {
            for (int j_nir = 0; j_nir < nbands_rhownir; j_nir++) {
                derv_nlw_rhorc_rrs_cov_product[ib * nbands_rhownir + j_nir] = 0;
                int ind = bindex_rhownir[j_nir];
                tmp = 0.;
                for (int i_ac = 0; i_ac < nbands_ac; i_ac++) {
                    tmp += -derv_nlw_rhorc[ib][acbands_index[i_ac]] * derv_rhownir_rrs[i_ac * nwave + ind];
                }
                derv_nlw_rhorc_rrs_product[nbands_rhownir * ib + j_nir] = tmp;
            }
        }

        static float * __restrict__  cov_temp_array = NULL;
        if (cov_temp_array == NULL) {
            cov_temp_array = (float *)malloc(sizeof(float) * nbands_rhownir * nbands_rhownir);
        }
        for (int j_nir = 0; j_nir < nbands_rhownir; j_nir++) {
            for (int i_nir = 0; i_nir < nbands_rhownir; i_nir++) {
                int ind1 = bindex_rhownir[j_nir];
                int ind2 = bindex_rhownir[i_nir];
                cov_temp_array[j_nir * nbands_rhownir + i_nir] = covariance_matrix[ind1 * nwave + ind2];
            }
        }

        for (ib = 0; ib < nwave; ib++) {
            for (int j_nir = 0; j_nir < nbands_rhownir; j_nir++) {
                for (int i_nir = 0; i_nir < nbands_rhownir; i_nir++) {
                    derv_nlw_rhorc_rrs_cov_product[ib * nbands_rhownir + i_nir] +=
                        derv_nlw_rhorc_rrs_product[ib * nbands_rhownir + j_nir] *
                        cov_temp_array[j_nir * nbands_rhownir + i_nir];
                }
            }
        }

        static float * __restrict__  temp_array2 = NULL;
        if (temp_array2 == NULL) {
            temp_array2 = (float *)malloc(sizeof(float) * nwave * nbands_rhownir);
        }
        // transpose derv_nlw_rhorc_rrs_cov_product and derv_nlw_rhorc_rrs_product using temp_array2
        for ( ib = 0; ib < nwave; ib++) {
            for (int j_nir = 0; j_nir < nbands_rhownir; j_nir++) {
                temp_array2[j_nir * nwave + ib] = derv_nlw_rhorc_rrs_cov_product[nbands_rhownir * ib + j_nir];
            }
        }
        memcpy(derv_nlw_rhorc_rrs_cov_product, temp_array2, sizeof(float) * nwave * nbands_rhownir);
        for ( ib = 0; ib < nwave; ib++) {
            for (int j_nir = 0; j_nir < nbands_rhownir; j_nir++) {
                temp_array2[j_nir * nwave + ib] = derv_nlw_rhorc_rrs_product[nbands_rhownir * ib + j_nir];
            }
        }

        memcpy(derv_nlw_rhorc_rrs_product, temp_array2, sizeof(float) * nwave * nbands_rhownir);

        if (delta_cov_nir_array == NULL) {
            delta_cov_nir_array = (float *)malloc(sizeof(float) * nwave * nwave);
        }
        for ( ib = 0; ib < nwave; ib++) {
            if (proc_uncertainty == 1) {
                bindex_covaraince[0] = ib;
                for (ix1 = 1; ix1 < nbands_covariance; ix1++) {
                    if (ib == bindex_covaraince[ix1]) {
                        bindex_covaraince[0] = BAD_INT;
                        break;
                    }
                }
            } else //if (proc_uncertainty == 2)
                nbands_covariance = nwave - ib;

            for (ix1 = 0; ix1 < nbands_covariance; ix1++) {
                if (proc_uncertainty == 1) {
                    iw = bindex_covaraince[ix1];
                    if (iw == BAD_INT)
                        continue;
                } else //if (proc_uncertainty == 2)
                    iw = ib + ix1;

                delta_cov_nir_array[ib * nwave + iw] = 0;
            }
        }

        for (int j_nir = 0; j_nir < nbands_rhownir; j_nir++) {
            for ( ib = 0; ib < nwave; ib++) {
                if (proc_uncertainty == 1) {
                    bindex_covaraince[0] = ib;
                    for (ix1 = 1; ix1 < nbands_covariance; ix1++) {
                        if (ib == bindex_covaraince[ix1]) {
                            bindex_covaraince[0] = BAD_INT;
                            break;
                        }
                    }
                } else //if (proc_uncertainty == 2)
                    nbands_covariance = nwave - ib;

                for (int ix = 0; ix < nbands_covariance; ix++) {
                    if (proc_uncertainty == 1) {
                        iw = bindex_covaraince[ix];
                        if (iw == BAD_INT)
                            continue;
                    } else //if (proc_uncertainty == 2)
                        iw = ib + ix;

                    delta_cov_nir_array[ib * nwave + iw] +=
                        derv_nlw_rhorc_rrs_cov_product[j_nir * nwave + ib] *
                        derv_nlw_rhorc_rrs_product[j_nir * nwave + iw];
                }
            }
        }
    }
    static float *__restrict__ glint_corr_temp = NULL;
    if (want_glintcorr) {
        if (glint_corr_temp == NULL) {
            glint_corr_temp = (float *)malloc(sizeof(float) * nwave);
        }
        for (ib = 0; ib < nwave; ib++) {
            glint_corr_temp[ib] = -1 / t_sol[ib] / t_sen[ib];
            for (int i_ac = 0; i_ac < nbands_ac; i_ac++) {
                int index = acbands_index[i_ac];
                glint_corr_temp[ib] -= derv_nlw_rhorc[ib][index] * radref[index];
            }
        }
    }
    /* calculate cov(Rrs[ib],Rrs[iw])   */
    ix2 = nbands_rhownir;
    if (proc_uncertainty == 2)
        ix2 = nwave;

    /* only calculate the covariance at the bands used to calculate chl within the iteration to account for rhownir*/
    for (ix3 = 0; ix3 < ix2; ix3++) {
        ib = bindex_rhownir[ix3];
        if (proc_uncertainty == 2)
            ib = ix3;

        if (proc_uncertainty == 1) {
            bindex_covaraince[0] = ib;
            for (ix1 = 1; ix1 < nbands_covariance; ix1++) {
                if (ib == bindex_covaraince[ix1]) {
                    bindex_covaraince[0] = BAD_INT;
                    break;
                }
            }
        } else //if (proc_uncertainty == 2)
            nbands_covariance = nwave - ib;

        for (ix1 = 0; ix1 < nbands_covariance; ix1++) {
            if (proc_uncertainty == 1) {
                iw = bindex_covaraince[ix1];
                if (iw == BAD_INT)
                    continue;
            } else //if (proc_uncertainty == 2)
                iw = ib + ix1;

            cov = covariance_matrix_test[ib * nwave + iw];
            cov += derv_nlw_rh[ib] * derv_nlw_rh[iw] * uncertainty->drh[ip] * uncertainty->drh[ip];

            if (rhownir_corr) {
                /* contribution from chl through rhownir */
                cov += derv_nlw_chl[ib] * derv_nlw_chl[iw] * uncertainty->dchl * uncertainty->dchl;

                /* contribution from Rrs through rhownir */
                cov += delta_cov_nir_array[ib * nwave + iw];
            }
            /* add the contribution from glint coefficient resulted from uncertianty in wind speed*/
            if (want_glintcorr) {
                float derv_ws[2];  //
                derv_ws[0] = glint_corr_temp[ib];
                derv_ws[1] = glint_corr_temp[iw];

                derv_ws[0] *= (uncertainty->derv_TLg_gc[ib] * uncertainty->derv_gc_ws[ip]);
                derv_ws[1] *= (uncertainty->derv_TLg_gc[iw] * uncertainty->derv_gc_ws[ip]);

                cov += derv_ws[0] * derv_ws[1] * uncertainty->dws[ip] * uncertainty->dws[ip];
            }

            covariance_matrix_test[ib * nwave + iw] = cov / Fobar[ib] / Fobar[iw];
            covariance_matrix_test[iw * nwave + ib] = covariance_matrix_test[ib * nwave + iw];
        }
    }
    memcpy(covariance_matrix, covariance_matrix_test, sizeof(float) * nwave * nwave);
}