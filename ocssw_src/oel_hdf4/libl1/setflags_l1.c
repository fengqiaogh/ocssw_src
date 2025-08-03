#include "l1.h"

void l1_mask_set(l1str *l1rec, int32_t ip) {
    l1rec->mask[ip] = (l1_input->landmask && l1rec->land [ip]) ||
            (l1_input->bathmask && l1rec->swater [ip]) ||
            (l1_input->cloudmask && l1rec->cloud [ip]) ||
            (l1_input->glintmask && l1rec->glint [ip]) ||
            (l1_input->stlightmask && l1rec->stlight[ip]) ||
            (l1_input->satzenmask && l1rec->senzmax[ip]) ||
            (l1_input->sunzenmask && l1rec->solzmax[ip]) ||
            (l1_input->hiltmask && l1rec->hilt [ip]) ||
            (l1rec->filter[ip]) ||
            (l1rec->navfail[ip]);
}

int setflags(l1str *l1rec) {
    static int firstCall = 1;

    static int ib412;
    static int ib555;
    static int ib670;
    static int ib865;
    static int ibcloud;

    int32_t evalmask = l1_input->evalmask;
    float mu0;
    int32_t ip;
    float albedo;
    filehandle *l1file = l1rec->l1file;
    int32_t nwave = l1file->nbands;

    if (firstCall) {
        firstCall = 0;
        ib412 = windex(412., l1file->fwave, nwave);
        ib555 = windex(550., l1file->fwave, nwave);
        ib670 = windex(670., l1file->fwave, nwave);
        ib865 = windex(865., l1file->fwave, nwave);
        ibcloud = windex(l1_input->cloud_wave, l1file->fwave, nwave);

        printf("\nUsing %6.1f nm channel for cloud flagging over water.\n", l1file->fwave[ibcloud]);
        printf("Using %6.1f nm channel for cloud flagging over land.\n\n", l1file->fwave[ib412]);

    }

    /*                                                     */
    /* Now set flags for each pixel                        */
    /*                                                     */
    for (ip = 0; ip < l1rec->npix; ip++) {

        mu0 = cos(l1rec->solz[ip] / OEL_RADEG);

        /* Check view angle limits */
        if (l1rec->senz[ip] > l1_input->satzen) {
            l1rec->senzmax[ip] = ON;
        }
        if (l1rec->solz[ip] > l1_input->sunzen) {
            l1rec->solzmax[ip] = ON;
        }


        /* Check for glint */
        if (!l1rec->land[ip] && l1rec->solz[ip] < SOLZNIGHT) {
            if (l1rec->glint_coef[ip] > l1_input->glint) {
                l1rec->glint[ip] = ON;
                if(l1rec->glint_coef[ip] > l1_input->extreme_glint )
                    l1rec->hilt[ip]=ON;
            }
        }


        /* Check for clouds (daytime only) */
        if (l1rec->solz[ip] < SOLZNIGHT) {

            if (!l1rec->land[ip]) {
                l1rec->cloud_albedo[ip] = l1rec->rhos[ip * nwave + ibcloud]
                        - l1rec->TLg[ip * nwave + ibcloud] * OEL_PI / mu0 / l1rec->Fo[ibcloud];
                albedo = l1_input->albedo;
            } else {
                l1rec->cloud_albedo[ip] = l1rec->rhos[ip * nwave + ib412];
                albedo = l1_input->albedo * 4.0;
            }
            // If cloud mask file is there, use that 
            if ( l1_input->cld_msk_file[0]) {
                char cld_category;
                l1rec->cloud[ip] = get_sdps_cld_mask( l1rec, ip, &cld_category );
            }  else {
                //  use all other standard cloud data
                if ((evalmask & MODCLOUD) != 0) {
                    l1rec->cloud[ip] = get_cldmask(l1rec, ip);
                } else {
                    if (l1rec->cloud_albedo[ip] > albedo)
                        l1rec->cloud[ip] = ON;
                }
                /* Eval switch to enable spectral cloud test, recover bright water */
                if ((l1_input->cloud_eps > 0.0) != 0) {
                    float min_cld_rhos = 1.0;
                    float max_cld_rhos = 1.0;
                    float cld_rhos_ratio = 0.0;
                    if (l1rec->cloud[ip] == ON && !l1rec->land[ip] && l1rec->cloud_albedo[ip] < albedo * 10) {
                        min_cld_rhos = MIN(MIN(MIN(
                            l1rec->rhos[ip * nwave + ib865],
                            l1rec->rhos[ip * nwave + ib670]),
                            l1rec->rhos[ip * nwave + ib555]),
                            l1rec->rhos[ip * nwave + ib412]);
                        max_cld_rhos = MAX(MAX(MAX(
                            l1rec->rhos[ip * nwave + ib865],
                            l1rec->rhos[ip * nwave + ib670]),
                            l1rec->rhos[ip * nwave + ib555]),
                            l1rec->rhos[ip * nwave + ib412]);
                        if (min_cld_rhos > 0.0) {
                            cld_rhos_ratio = max_cld_rhos / min_cld_rhos;
                            if (cld_rhos_ratio > l1_input->cloud_eps && (evalmask & MODCLOUD) == 0)
                                l1rec->cloud[ip] = OFF;
                        }
                    }
                }
            }

        }

        /* Set masking */
        l1_mask_set(l1rec, ip);

    }

    return (1);
}

void setflagbits_l1(int level, l1str *l1rec, int32_t ipix) {
    int32_t npix, spix, epix;
    int32_t ip, ib, iw;
    int32_t nwave;

    if (l1rec == NULL) {
        printf("-E- %s line %d: attempt to set flags from NULL L1 record.\n",
                __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }

    npix = l1rec->npix;
    nwave = l1rec->l1file->nbands;

    if (ipix < 0) {
        spix = 0;
        epix = npix - 1;
    } else {
        spix = ipix;
        epix = ipix;
    }

    /* Level-0, Level-1A flag bits (require readl1) */
    if (level == 0) {
        for (ip = spix; ip <= epix; ip++) {
            if (l1rec->hilt [ip]) l1rec->flags[ip] |= HILT;
            if (l1rec->stlight[ip]) l1rec->flags[ip] |= STRAYLIGHT;
            if (l1rec->navfail[ip]) l1rec->flags[ip] |= NAVFAIL;
            if (l1rec->navwarn[ip]) l1rec->flags[ip] |= NAVWARN;
        }
    }
        /* Level-1B flag bits (require loadl1) */
    else if (level == 1) {
        for (ip = spix; ip <= epix; ip++) {
            if (l1rec->land [ip]) l1rec->flags[ip] |= LAND;
            if (l1rec->swater [ip]) l1rec->flags[ip] |= COASTZ;
            if (l1rec->cloud [ip]) l1rec->flags[ip] |= CLOUD;
            if (l1rec->ice [ip]) l1rec->flags[ip] |= SEAICE;
            if (l1rec->glint [ip]) l1rec->flags[ip] |= HIGLINT;
            if (l1rec->solzmax[ip]) l1rec->flags[ip] |= HISOLZEN;
            if (l1rec->senzmax[ip]) l1rec->flags[ip] |= HISATZEN;
            if (l1rec->filter [ip]) l1rec->flags[ip] |= FILTER;
            if (l1rec->glint_coef[ip] > GLINT_MIN)
                l1rec->flags[ip] |= MODGLINT;
            for (iw = 0; iw < nwave; iw++) {
                ib = iw;
                if (l1rec->dpol[ip * nwave + ib] > l1_input->hipol) {
                    l1rec->flags[ip] |= HIPOL;
                    break;
                }
            }
        }

    } else {
        printf("-E- %s line %d: attempt to set l1 flags at bogus level %d.\n",
                __FILE__, __LINE__, level);
        exit(FATAL_ERROR);
    }
}

