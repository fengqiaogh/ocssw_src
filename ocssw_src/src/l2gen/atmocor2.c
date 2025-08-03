#include "l12_proto.h"
#include "atrem_corl1.h"
#include "glint.h"

/* --------------------------------------------------------------------------------------- */
/* atmocor2() - atmospheric correction, converts Lt -> nLw                                 */
/*                                                                                         */
/* C-version, September 2004, B. Franz                                                     */

/* --------------------------------------------------------------------------------------- */

int atmocor2(l2str *l2rec, aestr *aerec, int32_t ip) {
    static int firstCall = 1;
    static int want_nirRrs = 0;
    static int want_nirLw = 0;
    static int want_mumm = 0;
    static int want_ramp = 1;
    static int32_t aer_iter_max = 1;
    static int32_t aer_iter_min = 1;

    static float pi = OEL_PI;
    static float radeg = OEL_RADEG;
    static float badval = BAD_FLT;
    static float badchl = BAD_FLT;
    static float p0 = STDPR;
    static float df = 0.33;
    static float cbot = 0.7;
    static float ctop = 1.3;
    static float seed_chl = 0.0;
    static float seed_green = 0.0;
    static float seed_red = 0.0;
    static float nir_chg = 0.02;
    static float glint_min = GLINT_MIN;
    static float cslp;
    static float cint;

    static int32_t green;
    static int32_t red;
    static int32_t nir_s;
    static int32_t nir_l;
    static int32_t swir_s;
    static int32_t swir_l;
    static int32_t aer_s;
    static int32_t aer_l;
    static int32_t daer;
    static int32_t aer_base;
    static int32_t nwvis;
    static int32_t nwave;
    static float *wave;
    static float *gc;
    float glint_coef_q,glint_coef_u;
    double nw;
    int32_t nWaveCovariance = 1;

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;
    uncertainty_t *uncertainty=l1rec->uncertainty;
    int32_t proc_uncertainty = input->proc_uncertainty;

    int32_t sensorID = l1file->sensorID;
    int32_t brdf_opt = input->brdf_opt;
    int32_t aer_opt = input->aer_opt;
    int32_t glint_opt = input->glint_opt;
    int32_t cirrus_opt = input->cirrus_opt;

    float *Fo = l1rec->Fo;
    float *Fobar = l1file->Fobar;
    float fsol = l1rec->fsol;
    float solz = l1rec->solz [ip];
    float senz = l1rec->senz [ip];
    float delphi = l1rec->delphi[ip];
    float ws = l1rec->ws [ip];

    int32_t *aermodmin = &l2rec->aermodmin[ip];
    int32_t *aermodmax = &l2rec->aermodmax[ip];
    float *aermodrat = &l2rec->aerratio [ip];
    int32_t *aermodmin2 = &l2rec->aermodmin2[ip];
    int32_t *aermodmax2 = &l2rec->aermodmax2[ip];
    float *aermodrat2 = &l2rec->aerratio2 [ip];
    float *eps = &l2rec->eps [ip];

    float *TLg = &l1rec->TLg [ip * l1file->nbands];
    float *Lw = &l2rec->Lw [ip * l1file->nbands];
    float *nLw = &l2rec->nLw [ip * l1file->nbands];
    float *La = &l2rec->La [ip * l1file->nbands];
    float *taua = &l2rec->taua [ip * l1file->nbands];
    float *t_sol = &l1rec->t_sol[ip * l1file->nbands];
    float *t_sen = &l1rec->t_sen[ip * l1file->nbands];
    float *brdf = &l2rec->brdf [ip * l1file->nbands];
    float *Rrs = &l2rec->Rrs [ip * l1file->nbands];

    float *Rrs_unc;

    double *tg_sol = &l1rec->tg_sol[ip * l1file->nbands];
    double *tg= &l1rec->tg[ip * l1file->nbands];;  // double way absorption

    float *dsensor=NULL, *dLr=NULL, *dtaua=NULL, *dtg_sol=NULL;
    float *dtg_sen=NULL, *dt_sol=NULL, *dt_sen=NULL, *dvc=NULL, *last_dtaua=NULL;
    float dLt;

    float last_dtaua_aer_l, *derv_taua_l, **derv_rhorc, *derv_rhow_l, *derv_rh;

    float *derv_Lg_taua=NULL; //derivative of TLg[wave] to corresponding taua[nir_l]
    float *drhown_nir=NULL;
    float *dchl;
    float *covariance_matrix;
    static float *corr_coef_rhot;
    int dim;
    float *F1, *F1_temp, *F2, **COV, **corr_nir;
    static float *delta_Lt = NULL;

    float *taur,*radref;
    float *tLw, *dtLw;
    float *rhown_nir;
    float *tLw_nir, *dtLw_nir, *last_tLw_nir, *dlast_tLw_nir;
    int32_t ib, ipb,inir,iw, i, j, ix1, ix2;
    int32_t status;
    int32_t iter, iter_max, iter_min, last_iter, iter_reset;
    float mu, mu0;
    int32_t nneg;
    float airmass;
    float *Ltemp, *dLtemp;
    float chl;
    int want_glintcorr;
    float refl_nir = 0, last_refl_nir = 0;
    float tindx;
    float Ka = 0.8; /* 1.375 um channel transmittance for water vapor (<=1)  Gao et al. 1998 JGR */
    float *mbac_w;
    static int nbands_ac=2;
    static int *acbands_index=NULL;

    if (firstCall == 1) {
        firstCall = 0;
        nwave = l1file->nbands;
        if ((wave = (float *) calloc(nwave, sizeof (float))) == NULL) {
            printf("-E- : Error allocating memory to wave\n");
            exit(FATAL_ERROR);
        }
        gc=(float *)malloc(nwave*sizeof(float));
        for (ib = 0; ib < nwave; ib++) {
            wave[ib] = l1file->fwave[ib];
        }
        if (sensorID == CZCS) {
            want_ramp = 0;
            seed_chl = 0.01;
            seed_green = 0.003;
            seed_red = 0.00125;
        }
        if(!input->acbands_index)
            input->acbands_index=(int *)malloc(nbands_ac*sizeof(int));
        acbands_index=input->acbands_index;
        nwvis = rdsensorinfo(l1file->sensorID, l1_input->evalmask, "NbandsVIS", NULL);
        nir_s = bindex_get(input->aer_wave_short);
        nir_l = bindex_get(input->aer_wave_long);
        aer_base=bindex_get(input->aer_wave_base);
        if (nir_s < 0 || nir_l < 0) {
            printf("Aerosol selection bands %d and %d not available for this sensor\n",
                    input->aer_wave_short, input->aer_wave_long);
            exit(1);
        }
        if (nir_l < nir_s) {
            printf("Invalid aerosol selection bands: long (%d) must be greater than short (%d).\n",
                   input->aer_wave_long, input->aer_wave_short);
            exit(1);
        }
        if (wave[nir_s] < 600) {
            printf("Aerosol selection band(s) must be greater than 600nm");
            exit(1);
        }

        aer_s = nir_s;
        aer_l = nir_l;
        daer = MAX(nir_l - nir_s, 1);
        cslp = 1. / (ctop - cbot);
        cint = -cslp * cbot;
        printf("Aerosol selection bands %d and %d\n", l1file->iwave[aer_s], l1file->iwave[aer_l]);

        switch (aer_opt) {
            case AERRHNIR:
            case FIXANGSTROMNIR:
            case FIXMODPAIRNIR:
                want_nirLw = 1;
                aer_iter_min = 1;
                aer_iter_max = input->aer_iter_max;
                printf("NIR correction enabled.\n");
                break;
            case AERRHMUMM:
                want_mumm = 1;
                aer_iter_min = 3;
                aer_iter_max = input->aer_iter_max;
                printf("MUMM correction enabled.\n");
                break;
            case AERRHSWIR:
                want_nirLw = 1;
                aer_iter_min = 1;
                aer_iter_max = input->aer_iter_max;
                swir_s = bindex_get(input->aer_swir_short);
                swir_l = bindex_get(input->aer_swir_long);
                if (swir_s < 0 || swir_l < 0) {
                    printf("Aerosol selection bands %d and %d not available for this sensor\n",
                           input->aer_swir_short, input->aer_swir_long);
                    exit(1);
                }
                printf("NIR/SWIR switching correction enabled.\n");
                break;
            case AERRHMSEPS:
                want_nirLw = 1;
                aer_iter_min = 1;
                aer_iter_max = input->aer_iter_max;
                input->nbands_ac=nbands_ac;
                acbands_index[0]=windex(input->aer_wave_short,wave,nwave);
                acbands_index[1]=windex(input->aer_wave_long,wave,nwave);
                printf("NIR correction enabled --> for multi-scattering epsilon.\n");
                break;
            case AERRHSM:
            want_nirLw = 1; //This needs to be turned on, but a new SWIR water model is needed for it to work
                aer_iter_min = 1;
                aer_iter_max = input->aer_iter_max;
            daer=1;
            nbands_ac=input->nbands_ac;

            if(input->mbac_wave[0]<input->aer_wave_short){
                printf("%s line %d: the first mbac wavelength %d shouldn't be shorter than aer_wave_short %d \n", __FILE__, __LINE__,input->mbac_wave[0],input->aer_wave_short);
                    exit(1);
                }
            for(ib=0;ib<nbands_ac;ib++){
                acbands_index[ib]=0;
                for(iw=aer_s;iw<nwave;iw++){
                    if(input->mbac_wave[ib]==l1file->iwave[iw]){
                        acbands_index[ib]=1;
                            break;
                        }
                    }
                if(!acbands_index[ib]){
                    printf("%s line %d: the band %d specified in mbac_wave doesn't exist \n", __FILE__, __LINE__,input->mbac_wave[ib]);
                        exit(1);
                    }
                    acbands_index[ib]=windex(input->mbac_wave[ib],wave,nwave);
                }

                printf("NIR correction enabled --> for spectral matching.\n");
                break;
            default:
                if (input->aer_rrs_short >= 0.0 && input->aer_rrs_long >= 0.0) {
                    want_nirRrs = 1;
                    aer_iter_min = 3;
                    aer_iter_max = input->aer_iter_max;
                    printf("NIR correction via input Rrs enabled.\n");
                }
                break;
        }
        if (input->aer_iter_max < 1)
            want_nirLw = 0;
        if ( want_nirLw || want_nirRrs) {
            if ((red = bindex_get(670)) < 0) {
                if ((red = bindex_get(680)) < 0)
                    if ((red = bindex_get(620)) < 0)
                        if ((red = bindex_get(765)) < 0)
                            if ((red = bindex_get(655)) < 0)
                                if ((red = bindex_get(664)) < 0) /* added for MSI */ {
                                    printf("%s line %d: can't find red band\n", __FILE__, __LINE__);
                                    exit(1);
                                }
            }
            if ((green = bindex_get(550)) < 0) {
                if ((green = bindex_get(555)) < 0)
                    if ((green = bindex_get(560)) < 0)
                        if ((green = bindex_get(565)) < 0) {
                            printf("%s line %d: can't find green band\n", __FILE__, __LINE__);
                            exit(1);
                        }
            }
        }

        if(uncertainty){
            if ((delta_Lt = (float *)calloc(nwave, sizeof(float))) == NULL) {
                printf("-E- : Error allocating memory to delta_Lt\n");
                exit(FATAL_ERROR);
            }
            corr_coef_rhot = uncertainty->corr_coef_rhot;
        }
    }

    if ((taur = (float *)calloc(nwave, sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to taur\n");
        exit(FATAL_ERROR);
    }
    if ((tLw = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to tLw\n");
        exit(FATAL_ERROR);
    }
    if ((rhown_nir = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to rhown_nir\n");
        exit(FATAL_ERROR);
    }
    if ((tLw_nir = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to tLw_nir\n");
        exit(FATAL_ERROR);
    }
    if ((last_tLw_nir = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to last_tLw_nir\n");
        exit(FATAL_ERROR);
    }
    if ((Ltemp = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to Ltemp\n");
        exit(FATAL_ERROR);
    }
    if ((mbac_w = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to mbac_w\n");
        exit(FATAL_ERROR);
    }
    if ((radref = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to radref\n");
        exit(FATAL_ERROR);
    }

    if (uncertainty) {
        init_uncertainty(uncertainty,0);
        derv_Lg_taua=uncertainty->derv_Lg_taua;
        drhown_nir=uncertainty->drhown_nir;
        dchl=&uncertainty->dchl;
        Rrs_unc = &l2rec->Rrs_unc[ip * l1file->nbands];
        uncertainty->dRrs=Rrs_unc;
        if (proc_uncertainty == 2) {
            covariance_matrix = &l2rec->covariance_matrix[ip * l1file->nbands * l1file->nbands];
            uncertainty->covaraince_matrix = covariance_matrix;
        } else
            covariance_matrix = uncertainty->pixel_covariance;

        dsensor = &uncertainty->dsensor[ip * nwave];
        dLr = &uncertainty->dLr[ip * nwave];
        dtg_sol = &uncertainty->dtg_sol[ip * nwave];
        dtg_sen = &uncertainty->dtg_sen[ip * nwave];
        dt_sol = &uncertainty->dt_sol[ip * nwave];
        dt_sen = &uncertainty->dt_sen[ip * nwave];
        dvc = &uncertainty->dvc[ip * nwave];
        //dbrdf = &uncertainty->dbrdf[ip * nwave];

        if ((last_dtaua = (float *) calloc(nwave, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to Lunc\n");
            exit(FATAL_ERROR);
        }
        if ((dtaua = (float *) calloc(nwave, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to Lunc\n");
            exit(FATAL_ERROR);
        }

        if ((dLtemp = (float *) calloc(nwave, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to dLtemp\n");
            exit(FATAL_ERROR);
        }
        if ((dlast_tLw_nir = (float *) calloc(nwave, sizeof(float))) == NULL) {

            printf("-E- : Error allocating memory to last_tLw_nir\n");
            exit(FATAL_ERROR);
        }
        if ((dtLw = (float *) calloc(nwave, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to tLw\n");
            exit(FATAL_ERROR);
        }
        if ((dtLw_nir = (float *) calloc(nwave, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to dtLw_nir\n");
            exit(FATAL_ERROR);
        }
        if ((derv_taua_l = (float *)calloc(nwave, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to derv_taua_l\n");
            exit(FATAL_ERROR);
        }
        if ((derv_rhow_l = (float *)calloc(nwave, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to derv_rhow_l\n");
            exit(FATAL_ERROR);
        }
        if ((derv_rh = (float *)calloc(nwave, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to derv_rh\n");
            exit(FATAL_ERROR);
        }
        if ((derv_rhorc = (float **)calloc(nwave, sizeof(float *))) == NULL) {
            printf("-E- : Error allocating memory to derv_rhorc\n");
            exit(FATAL_ERROR);
        }
        for (ib = 0; ib < nwave; ib++) {
            if ((derv_rhorc[ib] = (float *)calloc(nbands_ac, sizeof(float))) == NULL) {
                printf("-E- : Error allocating memory to derv_rhorc[%d]\n", ib);
                exit(FATAL_ERROR);
            }
        }

        /// used for calculation cov(Rrs1,Rrs2)=F1.cov(X1,X2).F2

        dim = nbands_ac + 4;
        F1 = (float *)malloc(dim * sizeof(float));
        F1_temp = (float *)malloc(dim * sizeof(float));
        F2 = (float *)malloc(dim * sizeof(float));
        COV = (float **)malloc(dim * sizeof(float *));
        for (ib = 0; ib < dim; ib++)
            COV[ib] = (float *)malloc(dim * sizeof(float));

        for (ib = 0; ib < dim; ib++)
            for (iw = 0; iw < dim; iw++)
                COV[ib][iw] = 0.;

        corr_nir = (float **)malloc(nwave * sizeof(float *));
        for (ib = 0; ib < nwave; ib++)
            corr_nir[ib] = (float *)malloc(nbands_ac * sizeof(float));
        for (ib = 0; ib < nwave; ib++)
            for (iw = 0; iw < nbands_ac; iw++)
                corr_nir[ib][iw] = 0;

        for (iw = 0; iw < nwave; iw++) {
            for (ib = 0; ib < nbands_ac; ib++) {
                if ((acbands_index[ib]) >= iw)
                    corr_nir[iw][ib] = corr_coef_rhot[iw * nwave + acbands_index[ib]];
                else
                    corr_nir[iw][ib] = corr_coef_rhot[acbands_index[ib] * nwave + iw];
            }
        }
    }

    /* Initialize output values. If any input radiances are negative, just return. */

    status = 0;
    iter = 0;
    nneg = 0;

    for (ib = 0; ib < nwave; ib++) {
        ipb = ip * nwave + ib;

        //t_sol [ib]  = 1.0;    leave them as rayleigh only
        //t_sen [ib]  = 1.0;
        La [ib] = badval;
        tLw [ib] = badval;
        Lw [ib] = badval;
        nLw [ib] = badval;
        taua [ib] = badval;
        rhown_nir[ib]=0.;
        if (glint_opt != 2) TLg [ib] = 0.0;
        brdf [ib] = 1.0;
        Rrs [ib] = badval;

        l2rec->eps[ip] = badval;

        if (l2rec->l1rec->Lt[ipb] <= 0.0)
            if (wave[ib] < 1000) nneg++; /* don't test SWIR */

        if(aer_opt==AERRHSM)
            mbac_w[ib]=0.0;
        if (uncertainty) {
            dtaua[ib] = 0.;
            //derv_Lg_taua[ib] = 0.;
            last_dtaua[ib] = 0;
            // Rrs_unc[ib] = 0.;
            dvc[ib] = dvc[ib] / tg[ib] / l1rec->polcor[ipb];
        }
    }

    /* If any expected channels are negative */
    if (nneg > 0) {
        free(taur);
        free(tLw);
        free(rhown_nir);
        free(tLw_nir);
        free(last_tLw_nir);
        free(Ltemp);
        free(mbac_w);
        free(radref);

        if(uncertainty){
            free(dtLw);
            free(dtLw_nir);
            free(dlast_tLw_nir);
            free(dLtemp);
            free(last_dtaua);
            free(dtaua);
            free(derv_taua_l);
            free(derv_rhow_l);
            free(derv_rh);
            for (ib = 0; ib < nwave; ib++)
                free(derv_rhorc[ib]);
            free(derv_rhorc);

            free(F1);
            free(F2);
            free(F1_temp);
            for (ib = 0; ib < dim; ib++)
                free(COV[ib]);
            free(COV);

            for (ib = 0; ib < nwave; ib++)
                free(corr_nir[ib]);
            free(corr_nir);
        }

        status = 1;
        return (status);
    }

    mu0 = cos(solz / radeg);
    mu = cos(senz / radeg);
    airmass = 1.0 / mu0 + 1.0 / mu;

    /* Remove pre-computed atmospheric effects */
    /* --------------------------------------- */
    for (ib = 0; ib < nwave; ib++) {

        ipb = ip * nwave + ib;

        /* Pressure correct the Rayleigh optical thickness */
        taur[ib] = l1rec->pr[ip] / p0 * l1file->Tau_r[ib];

        /* Copy TOA radiance to temp var, eliminate missing bands */
        Ltemp[ib] = l2rec->l1rec->Lt[ipb];

        /* Correct for ozone absorption.  We correct for inbound and outbound here, */
        /* then we put the inbound back when computing Lw.                          */
        Ltemp[ib] = Ltemp[ib] / tg[ib];

        /* Do Cirrus correction - subtract off cirrus reflectance from Ltemp */
        /* Ka is the 1.375 um transmittance of water vapor above cirrus clouds */
        /* Add cirrus_opt to input */
        if (cirrus_opt) Ltemp[ib] -= l1rec->rho_cirrus[ip] / Ka * Fo[ib] * mu0 / pi;

        /*  Apply polarization correction */
        Ltemp[ib] /= l1rec->polcor[ipb];

        /* Remove whitecap radiance */
        Ltemp[ib] -= l1rec->tLf[ipb];

        /* Subtract the Rayleigh contribution for this geometry. */
        Ltemp[ib] -= l1rec->Lr[ipb];

        /* If selected, correct for Ding and Gordon O2 absorption for the aerosol component (replace O2 losses) */
        if (input->oxaband_opt == 1) {
            Ltemp[ib] /= l1rec->t_o2[ipb];
        }


        /* calculate error in Ltemp */
        if (uncertainty) {
            dLt = (1 - l2rec->l1rec->Lt[ipb] * uncertainty->derv_pol[ipb] / l1rec->polcor[ipb]);
            dLt *= (1 / l1rec->tg_sol[ipb] / l1rec->tg_sen[ipb]  / l1rec->polcor[ipb]);
            dLt = pow(dLt * dsensor[ib], 2);

            dLt += pow( l2rec->l1rec->Lt[ipb] / l1rec->tg_sol[ipb]  / pow(l1rec->tg_sen[ipb], 2) / l1rec->polcor[ipb] * dtg_sen[ib], 2);
            dLt += pow( l2rec->l1rec->Lt[ipb] / l1rec->tg_sen[ipb]  / pow(l1rec->tg_sol[ipb], 2) / l1rec->polcor[ipb] * dtg_sol[ib], 2);

            dLtemp[ib] = dLt + dLr[ib] * dLr[ib];
        }
        radref[ib]=pi/(Fo[ib] * mu0);
    }


    /* Compute BRDF correction, if not dependent on radiance */
    if (brdf_opt < FOQMOREL && brdf_opt > NOBRDF)
        ocbrdf(l2rec, ip, brdf_opt, wave, nwvis, solz, senz, delphi, ws, -1.0, NULL, NULL, brdf);

    /* Initialize iteration loop */
    chl = seed_chl;
    iter = 0;
    last_iter = 0;
    iter_max = aer_iter_max;
    iter_min = aer_iter_min;
    iter_reset = 0;
    last_refl_nir = 100.;
    want_glintcorr = 0;

    for (ib = 0; ib < nwave; ib++) {
        gc[ib] = l1rec->glint_coef[ip];
        if (input->glint_opt == 3) {
            nw = get_nw(wave[ib]);
            getglint_iqu(senz, solz, delphi, ws, 0., &gc[ib], &glint_coef_q, &glint_coef_u, nw);
        }   
    }
    /*  new glint_opt usage - a 2 will use the simple TLg from atmocor1 */
    if ((glint_opt == 1 || glint_opt==3) && gc[aer_base] > glint_min) {
        iter_max = MAX(2, iter_max);
        iter_min = MAX(2, iter_min);
        want_glintcorr = 1;
    }
    if (glint_opt == 2) {
        want_glintcorr = 1;
    }

    if (aer_opt == AERRHSWIR) {
        tindx_shi(l2rec, ip, &tindx);
        if (tindx >= 1.05) {
            iter_max = 1;
            aer_s = swir_s;
            aer_l = swir_l;
            want_nirLw = 0;
        } else {
            aer_s = nir_s;
            aer_l = nir_l;
            want_nirLw = 1;
        }
        daer = MAX(aer_l - aer_s, 1);
    }
NIRSWIR:

    if (want_nirLw || want_nirRrs) {
        for (ib = 0; ib < nwave; ib++) {
            last_tLw_nir[ib] = 0.0;
            tLw_nir[ib] = 0.0;
            Rrs[ib] = 0.0;

            if(uncertainty){
                dtLw_nir[ib] = 0.;
                dlast_tLw_nir[ib] = 0.;
            }
        }
        Rrs[green] = seed_green;
        Rrs[red ] = seed_red;
    }
    

    /* -------------------------------------------------------------------- */
    /* Begin iterations for aerosol with corrections for non-zero nLw(NIR) */
    /* -------------------------------------------------------------------- */

    if (aer_opt==AERRHSM ) {//&& tLw_nir[aer_s]>0.
        for(ib=0;ib<nbands_ac;ib++)
            mbac_w[acbands_index[ib]]=1.;
    }

    if (uncertainty) {
        for (ib = 0; ib < nwave; ib++)
            delta_Lt[ib] = sqrt(dLtemp[ib] + dvc[ib] * dvc[ib]);
    }

    while (!last_iter) {
        iter++;
        status = 0;

        if (uncertainty)
            last_dtaua_aer_l = dtaua[aer_l];

        /* Initialize tLw as surface + aerosol radiance */
        for (ib = 0; ib < nwave; ib++) {
            tLw[ib] = Ltemp[ib];
            if (uncertainty){
                dtLw[ib] = sqrt(dLtemp[ib]);
                last_dtaua[ib] = dtaua[ib];
            }
        }

        /*  Compute and subtract glint radiance */
        if (want_glintcorr) {

            if (glint_opt == 1 || glint_opt==3)
                glint_rad(iter, nwave, aer_s, aer_base, gc, airmass, mu0, Fo, taur, taua, tLw, TLg, uncertainty);

            for (ib = 0; ib < nwave; ib++) {
                tLw[ib] -= TLg[ib];
            }
        }

        /* Adjust for non-zero NIR water-leaving radiances using input Rrs */
        if (want_nirRrs) {

            rhown_nir[aer_s] = pi * input->aer_rrs_short;
            rhown_nir[aer_l] = pi * input->aer_rrs_long;

            for (ib = aer_s; ib <= aer_l; ib += daer) {

                /* Convert NIR reflectance to TOA W-L radiance */
                tLw_nir[ib] = rhown_nir[ib] / pi * Fo[ib] * mu0 * t_sol[ib] * t_sen[ib] / brdf[ib];

                /* Avoid over-subtraction */
                if (tLw_nir[ib] > tLw[ib] && tLw[ib] > 0.0)
                    tLw_nir[ib] = tLw[ib];

                /* Remove estimated NIR water-leaving radiance */
                tLw[ib] -= tLw_nir[ib];
            }
        }


        /* Adjust for non-zero NIR water-leaving radiances using MUMM */
        if (want_mumm) {

            get_rhown_mumm(l2rec, ip, aer_s, aer_l, rhown_nir);

            for (ib = aer_s; ib <= aer_l; ib += daer) {

                /* Convert NIR reflectance to TOA W-L radiance */
                tLw_nir[ib] = rhown_nir[ib] / pi * Fo[ib] * mu0 * t_sol[ib] * t_sen[ib] / brdf[ib];

                /* Remove estimated NIR water-leaving radiance */
                tLw[ib] -= tLw_nir[ib];
            }
        }


        /* Adjust for non-zero NIR water-leaving radiances using IOP model */
        if (want_nirLw) {

            ipb = ip*nwave;
            get_rhown_eval(input->fqfile, Rrs, wave, aer_s, aer_l, nwave, &l1rec->sw_a_avg[ipb], &l1rec->sw_bb_avg[ipb], chl, solz, senz, delphi, rhown_nir, l2rec,ip);

            for (ib = aer_s; ib <= aer_l; ib += daer) {

                /* Convert NIR reflectance to TOA W-L radiance */
                tLw_nir[ib] = rhown_nir[ib] / pi * Fo[ib] * mu0 * t_sol[ib] * t_sen[ib] / brdf[ib];
                if(uncertainty)
                    dtLw_nir[ib] = Fo[ib] * mu0 / (pi * brdf[ib]) *
                                   sqrt(pow(rhown_nir[ib] * t_sol[ib] * dt_sen[ib], 2) +
                                        pow(rhown_nir[ib] * t_sen[ib] * dt_sol[ib], 2) +
                                        pow(t_sen[ib] * t_sol[ib] * drhown_nir[ib], 2));

                /* Iteration damping */
                tLw_nir[ib] = (1.0 - df) * tLw_nir[ib] + df * last_tLw_nir[ib];

                if(uncertainty)
                    dtLw_nir[ib] = sqrt( pow((1. - df) * dtLw_nir[ib], 2) + pow(df * dlast_tLw_nir[ib], 2));

                /* Ramp-up ?*/
                if (want_ramp) {
                    if (chl > 0.0 && chl <= cbot) {
                        tLw_nir[ib] = 0.0;
                        if(uncertainty)
                            dtLw_nir[ib] = 0.0;
                    }
                    else if ((chl > cbot) && (chl < ctop)) {
                        tLw_nir[ib] *= (cslp * chl + cint);
                        if(uncertainty)
                            dtLw_nir[ib] *= (cslp * chl + cint); //error in chl should be included in the next
                    }
                }

                /* Remove estimated NIR water-leaving radiance */
                tLw[ib] -= tLw_nir[ib];
            }

            if(uncertainty && tLw_nir[aer_l] > 0.){
                for (ib = 0; ib <nbands_ac; ib ++){ 
                    inir=acbands_index[ib];
                    uncertainty->ratio_rhow[ib] = tLw_nir[inir] / tLw_nir[aer_l] * (Fo[aer_l] / Fo[inir]);
                }
            }
        }

        if (uncertainty) {
            if (aer_opt == AERRHMSEPS && fabs(uncertainty->derv_Lg_taua[aer_s]) > 0.) {
                uncertainty->derv_eps_taua_s = -1/ (Ltemp[aer_l] - TLg[aer_l] - tLw_nir[aer_l]) * derv_Lg_taua[aer_s];
                uncertainty->derv_eps_taua_s *= (Fo[aer_l] / Fo[aer_s]);

                uncertainty->derv_eps_taua_l = (Ltemp[aer_s] - TLg[aer_s] - tLw_nir[aer_s]) /
                                               pow(Ltemp[aer_l] - TLg[aer_l] - tLw_nir[aer_l], 2) *
                                               derv_Lg_taua[aer_l];
                uncertainty->derv_eps_taua_l *= (Fo[aer_l] / Fo[aer_s]);
                uncertainty->derv_eps_taua_l+=uncertainty->derv_eps_taua_s;
            }

            for (ib = 0; ib < nwave; ib++)
                uncertainty->derv_Lg_taua[ib] *= radref[ib];
        }

        /*  Compute the aerosol contribution */
        if (status == 0) {
            if (aer_opt != AERNULL)
                status = aerosol(l2rec, aer_opt, aerec, ip, wave, nwave, aer_s, aer_l, Fo, tLw,
                    La, t_sol, t_sen, eps, taua, aermodmin, aermodmax, aermodrat,
                    aermodmin2, aermodmax2, aermodrat2, mbac_w);
            else {
                for (ib = 0; ib < nwave; ib++) {

                    ipb = ip * nwave + ib;

                    t_sol [ib] = 1.0;
                    t_sen [ib] = 1.0;
                    La [ib] = 0.0;
                    taua [ib] = 0.0;
                    *eps = 1.0;
                }
            }
        }


        /* Subtract aerosols and normalize */
        if (status == 0) {

            for (ib = 0; ib < nwave; ib++) {

                /* subtract aerosol and normalize */
                tLw[ib] = tLw[ib] - La[ib];
                Lw [ib] = tLw[ib] / t_sen[ib] * tg_sol[ib];
                nLw[ib] = Lw [ib] / t_sol[ib] / tg_sol[ib] / mu0 / fsol * brdf[ib];

                if(uncertainty){

                    dt_sen[ib] = pow(uncertainty->derv_tsen_taua_l[ib] * last_dtaua_aer_l, 2);

                    for(inir=0;inir<nbands_ac;inir++){
                        i=acbands_index[inir];
                        dt_sen[ib] += pow(uncertainty->derv_tsen_rhorc[ib][inir] * sqrt(dLtemp[i]) * radref[i], 2);
                    }

                    dt_sen[ib] += pow(uncertainty->derv_tsen_rhow_l[ib] * dtLw_nir[aer_l] *  radref[aer_l], 2);
                    dt_sen[ib] += pow(uncertainty->derv_tsen_rh[ib] * uncertainty->drh[ip],2);

                    /* vicarious calibration contribution */
                    for(inir=0;inir<nbands_ac;inir++){
                        i=acbands_index[inir];
                        dt_sen[ib] += pow(uncertainty->derv_tsen_rhorc[ib][inir] *radref[i]* dvc[i], 2);
                    }

                    //TBD,only works for mseps right now, need to tune for mbac.
                    for(inir=0;inir<nbands_ac;inir++){
                        i=acbands_index[inir];
                        dt_sen[ib] +=2*corr_nir[i][nbands_ac - 1] *(uncertainty->derv_tsen_rhorc[ib][inir] * radref[i] * dvc[i] *uncertainty->derv_tsen_rhorc[ib][nbands_ac - 1] * radref[aer_l] * dvc[aer_l]);
                    }

                    if(dt_sen[ib]<0)
                        dt_sen[ib]=0.;
                    dt_sen[ib] = sqrt(dt_sen[ib]);
                    /*                                    */

                    dt_sol[ib] = pow(uncertainty->derv_tsol_taua_l[ib] * last_dtaua_aer_l, 2);

                    for(inir=0;inir<nbands_ac;inir++){
                        i=acbands_index[inir];
                        dt_sol[ib] += pow(uncertainty->derv_tsol_rhorc[ib][inir] * sqrt(dLtemp[i]) * radref[i], 2);
                    }

                    dt_sol[ib] += pow(uncertainty->derv_tsol_rhow_l[ib] * dtLw_nir[aer_l] * radref[aer_l], 2);
                    dt_sol[ib] += pow(uncertainty->derv_tsol_rh[ib] * uncertainty->drh[ip], 2);

                    /* vicarious calibration contribution */
                    for(inir=0;inir<nbands_ac;inir++){
                        i=acbands_index[inir];
                        dt_sol[ib] += pow(uncertainty->derv_tsol_rhorc[ib][inir] * radref[i]* dvc[i], 2);
                    }

                    for(inir=0;inir<nbands_ac;inir++){
                        i=acbands_index[inir];
                        dt_sol[ib] +=
                            2 * corr_nir[i][nbands_ac - 1] *
                            (uncertainty->derv_tsol_rhorc[ib][inir] * radref[i] * dvc[i] *
                             uncertainty->derv_tsol_rhorc[ib][nbands_ac - 1] * radref[aer_l] * dvc[aer_l]);
                    }

                    if(dt_sol[ib]<0)
                        dt_sol[ib]=0;
                    dt_sol[ib] = sqrt(dt_sol[ib]);
                    /*                                    */

                    derv_taua_l[ib] = -derv_Lg_taua[ib] / (t_sol[ib] * t_sen[ib]) * (1 / radref[ib]);
                    derv_taua_l[ib] += -uncertainty->derv_La_taua_l[ib] / (t_sol[ib] * t_sen[ib]);
                    derv_taua_l[ib] +=
                        (-tLw[ib] / (t_sol[ib] * t_sen[ib] * t_sen[ib]) * uncertainty->derv_tsen_taua_l[ib]);
                    derv_taua_l[ib] +=
                        (-tLw[ib] / (t_sol[ib] * t_sol[ib] * t_sen[ib]) * uncertainty->derv_tsol_taua_l[ib]);

                    for (inir = 0; inir < nbands_ac; inir++) {
                        derv_rhorc[ib][inir] =
                            -uncertainty->derv_La_rhorc[ib][inir] / (t_sol[ib] * t_sen[ib]);
                        derv_rhorc[ib][inir] += (-tLw[ib] / (t_sol[ib] * t_sen[ib] * t_sen[ib]) *
                                                 uncertainty->derv_tsen_rhorc[ib][inir]);
                        derv_rhorc[ib][inir] += (-tLw[ib] / (t_sol[ib] * t_sol[ib] * t_sen[ib]) *
                                                 uncertainty->derv_tsol_rhorc[ib][inir]);
                    }

                    derv_rhow_l[ib] = -uncertainty->derv_La_rhow_l[ib] / (t_sol[ib] * t_sen[ib]);
                    derv_rhow_l[ib] +=
                        (-tLw[ib] / (t_sol[ib] * t_sen[ib] * t_sen[ib]) * uncertainty->derv_tsen_rhow_l[ib]);
                    derv_rhow_l[ib] +=
                        (-tLw[ib] / (t_sol[ib] * t_sol[ib] * t_sen[ib]) * uncertainty->derv_tsol_rhow_l[ib]);

                    derv_rh[ib] = -uncertainty->derv_La_rh[ib] / (t_sol[ib] * t_sen[ib]);
                    derv_rh[ib] +=
                        (-tLw[ib] / (t_sol[ib] * t_sen[ib] * t_sen[ib]) * uncertainty->derv_tsen_rh[ib]);
                    derv_rh[ib] +=
                        (-tLw[ib] / (t_sol[ib] * t_sol[ib] * t_sen[ib]) * uncertainty->derv_tsol_rh[ib]);

                    /*  end of vicarious calibration   */

                    /*   end of using a different way to calculate dnLw  */

                    dtaua[ib] = pow(uncertainty->derv_taua_taua_l[ib] * last_dtaua_aer_l, 2);
                    for (inir = 0; inir < nbands_ac; inir++) {
                        i = acbands_index[inir];
                        dtaua[ib] +=
                            pow(uncertainty->derv_taua_rhorc[ib][inir] * sqrt(dLtemp[i]) * radref[i], 2);
                    }

                    dtaua[ib] += pow(uncertainty->derv_taua_rhow_l[ib] * dtLw_nir[aer_l] * radref[aer_l], 2);
                    dtaua[ib] += pow(uncertainty->derv_taua_rh[ib] * uncertainty->drh[ip], 2);

                    /* vicarious calibration contribution */
                    for (inir = 0; inir < nbands_ac; inir++) {
                        i = acbands_index[inir];
                        dtaua[ib] += pow(uncertainty->derv_taua_rhorc[ib][inir] * radref[i] * dvc[i], 2);
                    }
                    dtaua[ib] +=
                        2 * corr_nir[aer_s][nbands_ac - 1] *
                        (uncertainty->derv_taua_rhorc[ib][nbands_ac - 1] * radref[aer_l] * dvc[aer_l] *
                         uncertainty->derv_taua_rhorc[ib][0] * radref[aer_s] * dvc[aer_s]);

                    /*                                    */
                    if (dtaua[ib] < 0)
                        dtaua[ib] = 0;

                    dtaua[ib] = sqrt(dtaua[ib]);
                }
            }

            /*calculate covariance matrix for Rrs                                         */
            /* cov(Rrs1, Rrs2)=F1.Cov(X1,X2).F2                                           */
            /*  F1: matrix of partial derivative of Rrs1 to X1(matrix)  dimension: 1*dim  */
            /*  F2: matrix of partial derivative of Rrs2 to X2(matrix)  dimension: dim*1  */
            /*  COV(X1,X2):variance-covariance matrix of (X1,X2)        dimension: dim*dim */
            /*  X1(rhorc1, rhorc_NIR, rh, rhow_l,taua_l)                                   */
            /*  X2(rhorc2, rhorc_NIR, rh, rhow_l,taua_l)                                   */

            if (uncertainty) {
                for (ib = 0; ib < nwave; ib++) {
                    if (proc_uncertainty == 2)
                        nWaveCovariance = nwave - ib;
                    for (iw = ib; iw < ib + nWaveCovariance; iw++) {
                        /* calculate cov(Rrs[ib],Rrs[iw])   */
                        float cov = 0;
                        F1[0] = 1 / t_sol[ib] / t_sen[ib] / radref[ib];
                        for (ix1 = 1; ix1 <= nbands_ac; ix1++)
                            F1[ix1] = derv_rhorc[ib][ix1 - 1];
                        F1[ix1] = derv_rh[ib];
                        F1[ix1 + 1] = derv_rhow_l[ib];
                        F1[ix1 + 2] = derv_taua_l[ib];

                        F2[0] = 1 / t_sol[iw] / t_sen[iw] / radref[iw];
                        for (ix1 = 1; ix1 <= nbands_ac; ix1++)
                            F2[ix1] = derv_rhorc[iw][ix1 - 1];
                        F2[ix1] = derv_rh[iw];
                        F2[ix1 + 1] = derv_rhow_l[iw];
                        F2[ix1 + 2] = derv_taua_l[iw];
                        cov += F1[0] * F2[0] *
                                       (corr_coef_rhot[ib * nwave + iw] * delta_Lt[ib] * radref[ib] *
                                        delta_Lt[iw] * radref[iw]);
                        // COV[0][0] = corr_coef_rhot[ib * nwave + iw] * delta_Lt[ib] * radref[ib] *
                        //             delta_Lt[iw] * radref[iw];
                        for (i = 0; i < nbands_ac; i++) {
                            ix1 = acbands_index[i];
                            cov += F1[0] * F2[i + 1] *  corr_nir[ib][i] * delta_Lt[ib] * radref[ib] * delta_Lt[ix1] * radref[ix1];
                            // COV[0][i + 1] =
                            //     corr_nir[ib][i] * delta_Lt[ib] * radref[ib] * delta_Lt[ix1] * radref[ix1];
                        }
                        // COV[0][nbands_ac + 1] = 0;
                        // COV[0][nbands_ac + 2] = 0;
                        // COV[0][nbands_ac + 3] = 0;

                        for (i = 0; i < nbands_ac; i++) {
                            ix2 = acbands_index[i];

                            cov+=F1[i + 1] * F2[0] * corr_nir[iw][i] * delta_Lt[iw] * radref[iw] * delta_Lt[ix2] * radref[ix2];
                            // COV[i + 1][0] =
                            //     corr_nir[iw][i] * delta_Lt[iw] * radref[iw] * delta_Lt[ix2] * radref[ix2];

                            for (j = 0; j < nbands_ac; j++) {
                                ix1 = acbands_index[j];
                                cov+=F1[i + 1] * F2[j + 1] * corr_nir[ix2][j] * delta_Lt[ix2] * radref[ix2] *
                                                    delta_Lt[ix1] * radref[ix1];
                                // COV[i + 1][j + 1] = corr_nir[ix2][j] * delta_Lt[ix2] * radref[ix2] *
                                //                     delta_Lt[ix1] * radref[ix1];
                            }
                            // COV[i + 1][nbands_ac + 1] = 0;
                            // COV[i + 1][nbands_ac + 2] = 0;
                            // COV[i + 1][nbands_ac + 3] = 0;
                        }

                        // for (ix2 = 1; ix2 < 4; ix2++)
                        //     for (ix1 = 0; ix1 < dim; ix1++)
                        //         COV[nbands_ac + ix2][ix1] = 0;
                        cov+=F1[nbands_ac + 1] * F2[nbands_ac + 1] *     uncertainty->drh[ip] * uncertainty->drh[ip];  
                        // COV[nbands_ac + 1][nbands_ac + 1] = pow(uncertainty->drh[ip], 2);
                        cov+=F1[nbands_ac + 2] * F2[nbands_ac + 2] *  dtLw_nir[aer_l] * radref[aer_l] *dtLw_nir[aer_l] * radref[aer_l];
                        // COV[nbands_ac + 2][nbands_ac + 2] = pow(dtLw_nir[aer_l] * radref[aer_l], 2);
                        cov+=F1[nbands_ac + 3] * F2[nbands_ac + 3] * last_dtaua_aer_l * last_dtaua_aer_l;
                        // COV[nbands_ac + 3][nbands_ac + 3] = pow(last_dtaua_aer_l, 2);

                        // for (ix1 = 0; ix1 < dim; ix1++) {
                        //     F1_temp[ix1] = 0.f;
                        // }

                        // for (ix2 = 0; ix2 < dim; ix2++)
                        //     for (ix1 = 0; ix1 < dim; ix1++) {
                        //         F1_temp[ix1] += F1[ix2] * COV[ix2][ix1];
                        //     }
                        covariance_matrix[ib * nwave + iw] = cov;
                        // for (ix1 = 0; ix1 < dim; ix1++)
                        //     covariance_matrix[ib * nwave + iw] += F1_temp[ix1] * F2[ix1];

                        covariance_matrix[ib * nwave + iw] *=
                            (brdf[ib] * brdf[iw] / (mu0 * mu0) / (fsol * fsol));
                    }
                }
            }

            /* Compute new estimated chlorophyll */
            if (want_nirLw) {
                refl_nir = Rrs[red];
                for (ib = aer_s; ib <= aer_l; ib += daer) {
                    last_tLw_nir[ib] = tLw_nir[ib];
                    if(uncertainty)
                        dlast_tLw_nir[ib] = dtLw_nir[ib];
                }
                for (ib = 0; ib < nwvis; ib++) {
                    Rrs[ib] = nLw[ib] / Fobar[ib];
                    if(uncertainty) {
                        for (iw = ib; iw < nwvis; iw++)
                            covariance_matrix[ib * nwave + iw] /= (Fobar[ib] * Fobar[iw]);
                        Rrs_unc[ib] = sqrt(covariance_matrix[ib * nwave + ib]);
                    }
                }
                chl = get_default_chl(l2rec, Rrs);

                // if we passed atmospheric correction but the spectral distribution of
                // Rrs is bogus (chl failed), assume this is a turbid-water case and
                // reseed iteration as if all 670 reflectance is from water.

                if (chl == badchl && iter_reset == 0 && iter < iter_max) {
                    chl = 10.0;
                    Rrs[red] = 1.0 * (Ltemp[red] - TLg[red]) / t_sol[red] / tg_sol[red] / mu0 / fsol / Fobar[red];
                    iter_reset = 1;
                    if(uncertainty)
                        *dchl=0.;
                }

                // if we already tried a reset, and still no convergence, force one last
                // pass with an assumption that all red radiance is water component, and
                // force iteration to end.  this will be flagged as atmospheric correction
                // failure, but a qualitatively useful retrieval may still result.

                if (chl == badchl && iter_reset == 1 && iter < iter_max) {
                    chl = 10.0;
                    Rrs[red] = 1.0 * (Ltemp[red] - TLg[red]) / t_sol[red] / tg_sol[red] / mu0 / fsol / Fobar[red];
                    iter = iter_max; // so iter will trigger maxiter flag and ATMFAIL
                    iter_reset = 2;
                    if(uncertainty)
                        *dchl=0.;
                }
            }

        } else {

            /* Aerosol determination failure */
            for (ib = 0; ib < nwave; ib++) {
                La [ib] = badval;
                tLw[ib] = badval;
                Lw[ib] = badval;
                nLw[ib] = badval;
                Rrs[ib] = badval;
            }
        }

        /* If NIR/SWIR switching, do secondary test for turbidity and reset if needed */
        if (iter == 1 && (aer_opt == AERRHSWIR)) {
            if (tindx >= 1.05 && status == 0) {
                for (ib = 0; ib < nwvis; ib++) {
                    Rrs[ib] = nLw[ib] / Fobar[ib];
                }
                chl = get_default_chl(l2rec, Rrs);
                //printf("Checking turbidity %d %f %f %f\n",ip,tindx,chl,nLw[nir_l]);
                if (chl > 0 && (chl <= 1.0 || nLw[nir_l] < 0.08)) {
                    iter_max = aer_iter_max;
                    aer_s = nir_s;
                    aer_l = nir_l;
                    daer = MAX(aer_l - aer_s, 1);
                    want_nirLw = 1;
                    //printf("Reverting to NIR %d %f %f %f\n",ip,tindx,chl,nLw[nir_l]);
                    goto NIRSWIR;
                } else
                    l1rec->flags[ip] |= TURBIDW;
            }
        }


        /* Shall we continue iterating */
        if (status != 0) {
            last_iter = 1;
        } else if (iter < iter_min) {
            last_iter = 0;
        } else if (want_nirLw && (fabs(refl_nir - last_refl_nir) < fabs(nir_chg * refl_nir) || refl_nir < 0.0)) {
            last_iter = 1;
        } else if (want_mumm || want_nirRrs) {
            last_iter = 1;
        } else if (iter > iter_max) {
            last_iter = 1;
        }

        last_refl_nir = refl_nir;
        if (aer_opt==AERRHSM) {//&& tLw_nir[aer_s]>0.
            for (ib = aer_s; ib <=  aer_l; ib ++) {
                if(sensorID==MODISA && ib<=aer_base && iter>2 && mbac_w[ib]!=0.){
                    mbac_w[ib] =  (iter*1.0/iter_max);
                    mbac_w[ib] = exp(-7*mbac_w[ib]);
                }
            }
        }

    } /* end of iteration loop */

    l2rec->num_iter[ip] = iter;


    /* If the atmospheric correction failed, we don't need to do more. */
    if (status != 0) {
        free(taur);
        free(tLw);
        free(rhown_nir);
        free(tLw_nir);
        free(last_tLw_nir);
        free(Ltemp);
        free(mbac_w);
        free(radref);

        if(uncertainty){
            free(dtLw);
            free(dtLw_nir);
            free(dlast_tLw_nir);
            free(dLtemp);
            free(last_dtaua);
            free(dtaua);
            free(derv_taua_l);
            free(derv_rhow_l);
            free(derv_rh);
            for (ib = 0; ib < nwave; ib++)
                free(derv_rhorc[ib]);
            free(derv_rhorc);

            free(F1);
            free(F2);
            free(F1_temp);
            for (ib = 0; ib < dim; ib++)
                free(COV[ib]);
            free(COV);

            for (ib = 0; ib < nwave; ib++)
                free(corr_nir[ib]);
            free(corr_nir);
        }
        return (status);

    }

    /* If we used a NIR Lw correction, we record the tLw as it was predicted. */
    if (want_nirLw || want_nirRrs || want_mumm) {
        for (ib = aer_s; ib <= aer_l; ib += daer) {
            tLw[ib] = tLw_nir[ib];
            Lw [ib] = tLw[ib] / t_sen[ib] * tg_sol[ib];
            nLw[ib] = Lw [ib] / t_sol[ib] / tg_sol[ib] / mu0 / fsol * brdf[ib];
            Rrs[ib] = nLw[ib] / Fobar[ib];
        }
    }

    /* Convert water-leaving radiances from band-averaged to band-centered.  */
    /* Switch mean solar irradiance from band-averaged to band centered also.*/
    if (l1_input->outband_opt >= 2) {
        nlw_outband(l1_input->evalmask, sensorID, wave, nwave, Lw, nLw, &l2rec->outband_correction[ip*nwave]);
        Fobar = l1file->Fonom;
    }

    /* Compute f/Q correction and apply to nLw */
    if (brdf_opt >= FOQMOREL) {
        ocbrdf(l2rec, ip, brdf_opt, wave, nwvis, solz, senz, delphi, ws, -1.0, nLw, Fobar, brdf);
        for (ib = 0; ib < nwvis; ib++) {
            nLw[ib] = nLw[ib] * brdf[ib];

            if(uncertainty) {
                Rrs_unc[ib] *= brdf[ib];
                for (iw = ib; iw < nwave; iw++)
                    covariance_matrix[ib * nwave + iw] *= (brdf[ib] * brdf[iw]);
            }
        }
    }

    /* Compute final Rrs */
    for (ib = 0; ib < nwave; ib++) {
        if (ib != aer_s && ib != aer_l) {
            Rrs[ib] = nLw[ib] / Fobar[ib];
            l2rec->Rrs[ip * nwave + ib] = Rrs[ib];
        }
    }

    if (uncertainty) {
        for (ib = nwvis; ib < nwave; ib++)
            for (iw = ib; iw < nwave; iw++)
                covariance_matrix[ib * nwave + iw] /= (Fobar[ib] * Fobar[iw]);

        for (ib = aer_s; ib < nwave; ib += daer) {
            if (Rrs_unc[ib] < 0)
                Rrs_unc[ib] = BAD_FLT;
            else {
                for (iw = ib; iw < nwave; iw++)
                    covariance_matrix[ib * nwave + iw] *= (brdf[ib] * brdf[iw]);
            }
        }
    }

    /* Compute final chl from final nLw (needed for flagging) */
    l2rec->chl[ip] = get_default_chl(l2rec, Rrs);

    /* if there is no lower bounding for aerosol selection, the error can't be quantified */
    if (uncertainty) {
        if (*aermodmin == *aermodmax || (*aermodrat2!=BAD_FLT && *aermodmin2 == *aermodmax2)) {
            for (ib = 0; ib < nwave; ib++) {
                Rrs_unc[ib] = BAD_FLT;
                for (iw = 0; iw < nwave; iw++)
                    covariance_matrix[ib * nwave + iw] = BAD_FLT;
            }
        } else {
            l2rec->chl_unc[ip] = *dchl;
        }
    }

    /*Determine Raman scattering contribution to Rrs*/
    run_raman_cor(l2rec, ip);

    free(taur);
    free(tLw);
    free(rhown_nir);
    free(tLw_nir);
    free(last_tLw_nir);
    free(Ltemp);
    free(mbac_w);
    free(radref);

    if(uncertainty){
        free(dtLw);
        free(dtLw_nir);
        free(dlast_tLw_nir);
        free(dLtemp);
        free(last_dtaua);
        free(dtaua);
        free(derv_taua_l);
        free(derv_rhow_l);
        free(derv_rh);
        for (ib = 0; ib < nwave; ib++)
            free(derv_rhorc[ib]);
        free(derv_rhorc);

        free(F1);
        free(F2);
        free(F1_temp);
        for (ib = 0; ib < dim; ib++)
            free(COV[ib]);
        free(COV);

        for (ib = 0; ib < nwave; ib++)
            free(corr_nir[ib]);
        free(corr_nir);
    }

    return (status);
}



