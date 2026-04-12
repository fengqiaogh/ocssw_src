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
    float derv_gc_ws; 
    float glint_coef_q,glint_coef_u;
    double nw;

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;
    uncertainty_t *uncertainty=l1rec->uncertainty;

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

    static float *Rrs_unc=NULL;
    float cov;

    double *tg_sol = &l1rec->tg_sol[ip * l1file->nbands];
    double *tg= &l1rec->tg[ip * l1file->nbands]; // double way absorption
    double *tg_sen = &l1rec->tg_sen[ip * l1file->nbands];

    float *dsensor=NULL, *dtg_sol=NULL;
    float *dtg_sen=NULL, *dvc=NULL,*last_derv_taua_rhorc=NULL;
    float dLt;

    float **derv_taua_rhorc; // derivative taua [nwave] to rhorc[nwave]
    float *dchl;
    static float *covariance_matrix=NULL;
    float F1[2], F2[2];
    static float *delta_Lt = NULL;

    float *taur,*radref;
    float *tLw;
    float *rhown_nir;
    float *tLw_nir, *last_tLw_nir;
    int32_t ib, ipb,iw, i, j, ix1, ix2;
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
    static int nbands_ac;
    static int *acbands_index=NULL;
    static float *derv_brdf_chl=NULL;
    float tmp,tmp_derv1,tmp_derv2=0.;
    float *last_derv_rhownir_rrs; //derivative of last_rhownir[nbands_ac] to rrs[nwave]
    float *last_derv_rhownir_chl; //derivative of last_rhownir[nbands_ac] to chl
    float *derv_rhownir_chl;
    float *derv_rhownir_rrs; // derivative of rhownir[nbands_ac] to rrs[nwave] from last iteration
    static int nbands_rhownir;  // No. of bands used in correction of non zero rhownir
    static int *bindex_rhownir; // index of the bands used in correction of non zero rhownir
    int rhownir_corr=0;            // 0: rhownir correction is not applied, 1: correction is applied
    static float *derv_outband_nlw=NULL;  // derivative of out-of-band correction to the nlw 

    if (firstCall == 1) {
        firstCall = 0;
        nwave = l1file->nbands;
        nbands_ac=input->nbands_ac;
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
                acbands_index[0]=windex(input->aer_wave_short,wave,nwave);
                acbands_index[1]=windex(input->aer_wave_long,wave,nwave);
                printf("NIR correction enabled --> for multi-scattering epsilon.\n");
                break;
            case AERRHSM:
            want_nirLw = 1; //This needs to be turned on, but a new SWIR water model is needed for it to work
                aer_iter_min = 1;
                aer_iter_max = input->aer_iter_max;
            daer=1;

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
            if ((delta_Lt = (float*)calloc(nwave, sizeof(float))) == NULL) {
                printf("-E- : Error allocating memory to delta_Lt\n");
                exit(FATAL_ERROR);
            }
            if ((derv_brdf_chl = (float*)calloc(nwave, sizeof(float))) == NULL) {
                printf("-E- : Error allocating memory to derv_brdf_chl\n");
                exit(FATAL_ERROR);
            }
            if ((covariance_matrix = (float*)calloc(nwave * nwave, sizeof(float))) == NULL) {
                printf("-E- : Error allocating memory to covariance_matrix\n");
                exit(FATAL_ERROR);
            }
            if ((Rrs_unc = (float*)calloc(nwave, sizeof(float))) == NULL) {
                printf("-E- : Error allocating memory to covariance_matrix\n");
                exit(FATAL_ERROR);
            }
            if ((derv_outband_nlw = (float*)calloc(nwave * nwave, sizeof(float))) == NULL) {
                printf("-E- : Error allocating memory to derv_outband_nlw\n");
                exit(FATAL_ERROR);
            }

            nbands_rhownir = uncertainty->nbands_rhownir;
            bindex_rhownir = uncertainty->bindx_rhownir;
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
        dchl=&uncertainty->dchl;
       // Rrs_unc = &l2rec->Rrs_unc[ip * l1file->nbands];
         uncertainty->dRrs=Rrs_unc;
       // uncertainty->covariance_matrix = covariance_matrix;

        ipb=ip*nwave;
        dsensor = &uncertainty->dsensor[ipb];
        dtg_sol = &uncertainty->dtg_sol[ipb];
        dtg_sen = &uncertainty->dtg_sen[ipb];
        dvc    = &uncertainty->dvc[ipb];
      
        if ((last_derv_taua_rhorc = (float *) calloc(nwave*nbands_ac, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to last_derv_taua_rhorc\n");
            exit(FATAL_ERROR);
        }

        if ((dLtemp = (float *) calloc(nwave, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to dLtemp\n");
            exit(FATAL_ERROR);
        }
        if ((derv_taua_rhorc = (float **)calloc(nwave, sizeof(float *))) == NULL) {
            printf("-E- : Error allocating memory to derv_taua_rhorc\n");
            exit(FATAL_ERROR);
        }
        if ((last_derv_rhownir_rrs = (float *) calloc(nwave*nbands_ac, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to last_derv_rhownir_rrs\n");
            exit(FATAL_ERROR);
        }
        if ((last_derv_rhownir_chl = (float *) calloc(nbands_ac, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to last_derv_rhownir_chl\n");
            exit(FATAL_ERROR);
        }
        if ((derv_rhownir_chl = (float *) calloc(nbands_ac, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to derv_rhownir_chl\n");
            exit(FATAL_ERROR);
        }
        if ((derv_rhownir_rrs = (float *) calloc(nbands_ac*nwave, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to derv_rhownir_rhorc\n");
            exit(FATAL_ERROR);
        }
        for (ib = 0; ib < nwave; ib++) {
            if ((derv_taua_rhorc[ib] = (float *)calloc(nwave, sizeof(float))) == NULL) {
                printf("-E- : Error allocating memory to derv_taua_rhorc[%d]\n", ib);
                exit(FATAL_ERROR);
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
            Rrs_unc[ib]=0.;
            for (i = 0; i < nwave; i++){
                covariance_matrix[ib * nwave + i] = 0.;
                if (i < nbands_ac)
                    last_derv_taua_rhorc[ib * nbands_ac + i] = 0.;
            }
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
            free(dLtemp);
            free(last_derv_taua_rhorc);

            for (ib = 0; ib < nwave; ib++)
                free(derv_taua_rhorc[ib]);
            free(derv_taua_rhorc);

            free(last_derv_rhownir_rrs);
            free(last_derv_rhownir_chl);
            free(derv_rhownir_chl);
            free(derv_rhownir_rrs);
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


        /* calculate uncertainty in Ltemp */
        if (uncertainty) {
            dLt = (1 /tg[ib]/ l1rec->polcor[ipb]);
            dLt = (dLt* dLt) * ( dsensor[ib]*dsensor[ib] );//+ dvc[ib]*dvc[ib]

            tmp=-l2rec->l1rec->Lt[ipb] /tg[ib]/l1rec->polcor[ipb]*(1/tg_sen[ipb]+uncertainty->derv_polcor_tgsen[ipb]/l1rec->polcor[ipb]);
            dLt += pow( tmp * dtg_sen[ib], 2);

            tmp=-l2rec->l1rec->Lt[ipb] /tg[ib]/l1rec->polcor[ipb]*(1/tg_sol[ipb]+uncertainty->derv_polcor_tgsol[ipb]/l1rec->polcor[ipb]);
            dLt += pow( tmp * dtg_sol[ib], 2);

            dLtemp[ib] = dLt;
        }
        radref[ib]=pi/(Fo[ib] * mu0);
    }


    /* Compute BRDF correction, if not dependent on radiance */
    if (brdf_opt < FOQMOREL && brdf_opt > NOBRDF)
        ocbrdf(l2rec, ip, brdf_opt, wave, nwvis, solz, senz, delphi, ws, -1.0, NULL, NULL, brdf,derv_brdf_chl);

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
        if (glint_opt == 3) {
            nw = get_nw(wave[ib]);
            getglint_iqu(senz, solz, delphi, ws, 0., &gc[ib], &glint_coef_q, &glint_coef_u, nw,&derv_gc_ws);
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
        }
        Rrs[green] = seed_green;
        Rrs[red ] = seed_red;
    }
    

    /* -------------------------------------------------------------------- */
    /* Begin iterations for aerosol with corrections for non-zero nLw(NIR) */
    /* -------------------------------------------------------------------- */

    if (aer_opt==AERRHSM ) {
        for(ib=0;ib<nbands_ac;ib++)
            mbac_w[acbands_index[ib]]=1.;
    }

    /*compute the uncertainty in TOA radiance */
    if (uncertainty) {
        for (ib = 0; ib < nwave; ib++){
            delta_Lt[ib] = sqrt(dLtemp[ib] + dvc[ib] * dvc[ib]);
            dvc[ib]=delta_Lt[ib];
        }
         for(ib=0;ib<nbands_ac;ib++){
            derv_rhownir_chl[ib]=0.;
            last_derv_rhownir_chl[ib]=0.;
         }           
    }

    while (!last_iter) {
        iter++;
        status = 0;

        /* Initialize tLw as surface + aerosol radiance */
        for (ib = 0; ib < nwave; ib++) {
            tLw[ib] = Ltemp[ib];
            
            if (uncertainty) {
               // uncertainty->derv_chl_rrs[ib]=0.;
                for (j = 0; j < nbands_ac; j++) {
                    ipb = ib * nbands_ac + j;

                    /* the max iteration is 2 for glint correction*/
                    if (iter <= 2)
                        last_derv_taua_rhorc[ipb] = derv_taua_rhorc[ib][acbands_index[j]];
                }
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

            if (iter>=2 && uncertainty) {
                for (ib = 0; ib < nbands_ac; ib++){
                    for (iw = 0; iw < nbands_rhownir; iw++) {
                        ipb=ib*nwave+bindex_rhownir[iw];
                        last_derv_rhownir_rrs[ipb] = uncertainty->derv_rhownir_rrs[ipb];
                    }
                    last_derv_rhownir_chl[ib] = uncertainty->derv_rhownir_chl[ib];
                }
            }
            rhownir_corr=0;

            ipb = ip*nwave;
            get_rhown_eval(input->fqfile, Rrs, wave, aer_s, aer_l, nwave, &l1rec->sw_a_avg[ipb], &l1rec->sw_bb_avg[ipb], chl, solz, senz, delphi, rhown_nir, l2rec,ip);

            for (ib = aer_s; ib <= aer_l; ib += daer) {

                /* Convert NIR reflectance to TOA W-L radiance */
                tLw_nir[ib] = rhown_nir[ib] / pi * Fo[ib] * mu0 * t_sol[ib] * t_sen[ib] / brdf[ib];

                /* Iteration damping */
                tLw_nir[ib] = (1.0 - df) * tLw_nir[ib] + df * last_tLw_nir[ib];

                /* Ramp-up ?*/
                if (want_ramp) {
                    if (chl > 0.0 && chl <= cbot) {
                        tLw_nir[ib] = 0.0;
                    }
                    else if ((chl > cbot) && (chl < ctop)) {
                        tLw_nir[ib] *= (cslp * chl + cint);
                    }
                }

                /* Remove estimated NIR water-leaving radiance */
                tLw[ib] -= tLw_nir[ib];
                if (tLw_nir[ib] != 0.0)
                    rhownir_corr = 1;
            }
            if (rhownir_corr && uncertainty) {
                for (ib = 0; ib < nbands_ac; ib++) {
                    i = acbands_index[ib];
                    tmp_derv1 = 1 / pi * Fo[i] * mu0 * t_sol[i] * t_sen[i] / brdf[i];

                    derv_rhownir_chl[ib]=(1.0-df)*uncertainty->derv_rhownir_chl[ib]+df*last_derv_rhownir_chl[ib];
                    derv_rhownir_chl[ib]*=tmp_derv1;
                    if ((chl > cbot) && (chl < ctop)) {
                        tmp = tLw_nir[i] / (cslp * chl + cint) * cslp;
                        //derv_rhownir_chl[ib]+=tmp;
                    }
                    derv_rhownir_chl[ib]*=radref[i];

                    for (iw = 0; iw < nbands_rhownir; iw++) {
                        ix2 = bindex_rhownir[iw];
                        ipb = ib * nwave + ix2;

                        derv_rhownir_rrs[ipb]= uncertainty->derv_rhownir_rrs[ipb]*(1.0 - df)+last_derv_rhownir_rrs[ipb] * df;
                        derv_rhownir_rrs[ipb] *= tmp_derv1;

                        /* contribution of rhorc through t_sol and t_sen*/
                        tmp=rhown_nir[i] / pi * Fo[i] * mu0 / brdf[i];
                       // derv_rhownir_rhorc[ipb]+= (tmp * t_sol[i]*derv_tsen_rhorc[i][ix2] );
                       // derv_rhownir_rhorc[ipb]+= (tmp * t_sen[i]*derv_tsol_rhorc[i][ix2] );  
                        
                        if ((chl > cbot) && (chl < ctop)){
                            derv_rhownir_rrs[ipb] *= (cslp * chl + cint);

                            /*contribution of rhorc through chl */
                            /* the hard coded ctop value cause artifacts in Rrs uncertainty*/
                            /*  remove this contribution temporarily.  CAUTION:!!!!!!*/
                            /*tmp=tLw_nir[i]/(cslp * chl + cint)*cslp;
                            tmp_derv2= derv_rhownir_rhorc[ipb];
                            for (j = 1; j < nbands_uncertainty; j++) {
                                ix1 = bindex_uncertainty[j];
                                if(uncertainty->derv_chl_rrs[ix1]!=0.)
                                    tmp_derv2+=tmp*uncertainty->derv_chl_rrs[ix1]/Fobar[ix1]*derv_nlw_rhorc[ix1][ix2];
                            }*/
                            //derv_rhownir_rhorc[ipb] = tmp_derv2;
                        }
                            
                        derv_rhownir_rrs[ipb] *= radref[i];
                    }
                }
            }
        }

        if (want_glintcorr && uncertainty) {
            for (ib = 0; ib < nwave; ib++){
                if (iter == 2)
                    uncertainty->derv_rhog_taua[ib] *= radref[ib];
                uncertainty->derv_rhog_rhorc_l[ib]*=(radref[ib]/radref[aer_l]);
            }
                
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
            }
            if (uncertainty){
                compute_uncertainty(l2rec, ip, aer_l, rhownir_corr, want_glintcorr, 1, radref, delta_Lt,
                                    derv_rhownir_chl, derv_rhownir_rrs, tLw, last_derv_taua_rhorc, derv_taua_rhorc, brdf,
                                    covariance_matrix);
                memcpy(uncertainty->covariance_matrix,covariance_matrix,nwave*nwave*sizeof(float));
                }

            /* Compute new estimated chlorophyll */
            if (want_nirLw) {
                refl_nir = Rrs[red];
                for (ib = aer_s; ib <= aer_l; ib += daer) {
                    last_tLw_nir[ib] = tLw_nir[ib];
                }
                for (ib = 0; ib < nwvis; ib++) {
                    Rrs[ib] = nLw[ib] / Fobar[ib];
                    if(uncertainty) {
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
                taua[ib]= badval;
            }
        }

        /* If NIR/SWIR switching, do secondary test for turbidity and reset if needed */
        if (iter == 1 && (aer_opt == AERRHSWIR)) {
            if (tindx >= 1.05 && status == 0) {
                for (ib = 0; ib < nwvis; ib++) {
                    Rrs[ib] = nLw[ib] / Fobar[ib];
                }
                chl = get_default_chl(l2rec, Rrs);
                if (chl > 0 && (chl <= 1.0 || nLw[nir_l] < 0.08)) {
                    iter_max = aer_iter_max;
                    aer_s = nir_s;
                    aer_l = nir_l;
                    daer = MAX(aer_l - aer_s, 1);
                    want_nirLw = 1;
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
            free(dLtemp);
            free(last_derv_taua_rhorc);

            for (ib = 0; ib < nwave; ib++) 
                free(derv_taua_rhorc[ib]);
            free(derv_taua_rhorc);

            free(last_derv_rhownir_rrs);
            free(derv_rhownir_rrs);
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

    if (uncertainty) {
        compute_uncertainty(l2rec, ip, aer_l, rhownir_corr, want_glintcorr, 2, radref, delta_Lt,
                            derv_rhownir_chl, derv_rhownir_rrs, tLw, last_derv_taua_rhorc, derv_taua_rhorc,
                            brdf, covariance_matrix);
        memcpy(uncertainty->covariance_matrix, covariance_matrix, nwave * nwave * sizeof(float));
    }

    /* Convert water-leaving radiances from band-averaged to band-centered.  */
    /* Switch mean solar irradiance from band-averaged to band centered also.*/ 
    float *outband_corr=&l2rec->outband_correction[ip*nwave];
    if (l1_input->outband_opt >= 2) { 

        Fobar = l1file->Fonom;
        if (uncertainty) { 
            for (ix1 = 0; ix1 < nwave; ix1++)
                for (ix2 = 0; ix2 < nwave; ix2++)
                    derv_outband_nlw[ix1 * nwave + ix2] =
                        covariance_matrix[ix1 * nwave + ix2] * Fobar[ix1] * Fobar[ix2];
        }

        nlw_outband(l1_input->evalmask, sensorID, wave, nwave, Lw, nLw, outband_corr,derv_outband_nlw,&tmp_derv2);

        if (uncertainty && tmp_derv2!=0.) {
             /*calculate covariance matrix accounting for the correction     */
            /* cov(Rrs1, Rrs2)=F1.Cov(X1,X2).F2                             */
            /*  F1: matrix of partial derivative of Rrs1 to X1(matrix)      */
            /*  F2: matrix of partial derivative of Rrs2 to X2(matrix)     */
            /*  COV(X1,X2):variance-covariance matrix of (X1,X2)           */
            /*  X1(nLw_unc[ib], ratio)  ratio is used in calculating outband_corr */
            /*  X2(nLw_unc[iw], ratio)                                     */
            for (ib = 0; ib < nwave; ib++) {
                F1[0]=outband_corr[ib];
                F1[1]=derv_outband_nlw[ib];
                for (iw = ib; iw < nwave; iw++) {
                    ipb=ib*nwave+iw;
                    tmp_derv1=Fobar[ib] * Fobar[iw];
                    F2[0] = outband_corr[iw];
                    F2[1] = derv_outband_nlw[iw];
                    
                    cov=F1[0]*F2[0]*covariance_matrix[ipb]*tmp_derv1;
                    cov+=F1[1]*F2[1]*tmp_derv2;
                    covariance_matrix[ipb]=cov/tmp_derv1;

                    if (iw == ib)
                        Rrs_unc[ib] = sqrt(covariance_matrix[ipb]);
                }
            }
        }

    }

    /* Compute f/Q correction and apply to nLw */
    if (brdf_opt >= FOQMOREL) {
        if (uncertainty) {
            tmp_derv1=*dchl;
            for (ib = 0; ib < nwave; ib++)
                for (iw = ib; iw < nwave; iw++) {
                    ipb = ib * nwave + iw;
                    derv_outband_nlw[ipb]=covariance_matrix[ipb];
                }
        }
        ocbrdf(l2rec, ip, brdf_opt, wave, nwvis, solz, senz, delphi, ws, -1.0, nLw, Fobar, brdf,
               derv_brdf_chl);
        for (ix1 = 0; ix1 < nwvis; ix1++) {
            nLw[ix1] = nLw[ix1] * brdf[ix1];
        }

        if (uncertainty && derv_brdf_chl[0] != 0) {
            /*calculate covariance matrix accounting for the correction     */
            /* cov(Rrs1, Rrs2)=F1.Cov(X1,X2).F2                             */
            /*  F1: matrix of partial derivative of Rrs1 to X1(matrix)      */
            /*  F2: matrix of partial derivative of Rrs2 to X2(matrix)     */
            /*  COV(X1,X2):variance-covariance matrix of (X1,X2)           */
            /*  X1(nLw_unc[ib], chl)  chl is used in calculating brdf     */
            /*  X2(nLw_unc[iw], chl)                                     */
            *dchl = tmp_derv1;
            for (ib = 0; ib < nwave; ib++) {
                F1[0] = brdf[ib];
                F1[1] = derv_brdf_chl[ib];
                for (iw = ib; iw < nwave; iw++) {
                    ipb = ib * nwave + iw;
                    tmp_derv1 = Fobar[ib] * Fobar[iw];
                    F2[0] = brdf[iw];
                    F2[1] = derv_brdf_chl[iw];

                    cov = F1[0] * F2[0] * derv_outband_nlw[ipb] * tmp_derv1;
                    cov += F1[1] * F2[1] * (*dchl * *dchl);
                    covariance_matrix[ipb] = cov / tmp_derv1;
                }
            }
        }
    }

    /* Compute final Rrs */
    for (ib = 0; ib < nwave; ib++) {
        if (ib != aer_s && ib != aer_l) {
            Rrs[ib] = nLw[ib] / Fobar[ib];
            l2rec->Rrs[ip * nwave + ib] = Rrs[ib];
            if (uncertainty)
                l2rec->Rrs_unc[ip*nwave+ib] = sqrt(covariance_matrix[ib * nwave + ib]);
        }
    }

    /* Compute final chl from final nLw (needed for flagging) */
    l2rec->chl[ip] = get_default_chl(l2rec, Rrs);

    /* if there is no lower or upper bounding aerosol, the uncertainty can't be quantified */
    if (uncertainty) {
        if (*aermodmin == *aermodmax || (*aermodrat2!=BAD_FLT && *aermodmin2 == *aermodmax2)) {
            for (ib = 0; ib < nwave; ib++) {
                l2rec->Rrs_unc[ip*nwave+ib] = BAD_FLT;
                for (iw = 0; iw < nwave; iw++)
                    covariance_matrix[ib * nwave + iw] = BAD_FLT;
            }
        } else {
            if (*dchl > 0 && l2rec->chl[ip]!=BAD_FLT)
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
        free(dLtemp);
        free(last_derv_taua_rhorc);
    
        for (ib = 0; ib < nwave; ib++)  
            free(derv_taua_rhorc[ib]);
        free(derv_taua_rhorc);

        free(last_derv_rhownir_rrs);
        free(last_derv_rhownir_chl);
        free(derv_rhownir_chl);
        free(derv_rhownir_rrs);
    }

    return (status);
}
