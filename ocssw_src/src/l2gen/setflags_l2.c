#include "l12_proto.h"

char isCoccolith(l2str *l2rec, int32_t ip);
char isTurbid(l2str *l2rec, int32_t ip);
float aerindex(l2str *l2rec, int32_t ip);
char isOptShallow(l2str *l2rec, int32_t ip);

enum abs_aer_option {
    RHOWINDEX = 1,
    TOTEXINCTION = 2
};


void setflagbits_l2(l2str *l2rec, int32_t ipix) {
    static int firstCall = 1;
    static int ib490;
    static int ib510;
    static int ib555;

    int32_t npix, spix, epix;
    int32_t ip;
    int32_t nwave; 

    l1str *l1rec = l2rec->l1rec;

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

    if (l2rec == NULL) {
        printf("-E- %s line %d: attempt to set flags from NULL L2 record.\n",
                __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }

    if (firstCall) {
        firstCall = 0;
        float* wave = l2rec->l1rec->l1file->fwave;
        ib490 = windex(490., wave, nwave);
        ib510 = windex(510., wave, nwave);
        ib555 = windex(550., wave, nwave);
    }

    for (ip = spix; ip <= epix; ip++) {
        if (l2rec->eps[ip] < input->epsmin ||
                l2rec->eps[ip] > input->epsmax || l2rec->chi2[ip]>5.) {
            l1rec->flags[ip] |= ATMWARN;
        }

        if (!l1rec->land[ip] && !l1rec->cloud[ip] && !l1rec->ice[ip]){
            if (isOptShallow(l2rec,ip))
                l1rec->flags[ip] |= OPSHAL;
        }

        if ((l2rec->Lw[ip * nwave + ib490] < 0.0) ||
                (l2rec->Lw[ip * nwave + ib510] < 0.0) ||
                (l2rec->Lw[ip * nwave + ib555] < 0.0)) {
            l1rec->flags[ip] |= ATMWARN;
        }

        if (l2rec->nLw[ip * nwave + ib555] < input->nlwmin)
            l1rec->flags[ip] |= LOWLW;

        if (isCoccolith(l2rec, ip))
            l1rec->flags[ip] |= COCCOLITH;

        if (input->aer_opt != AERWANGSWIR && input->aer_opt != AERRHSWIR) {
            if (isTurbid(l2rec, ip))
                l1rec->flags[ip] |= TURBIDW;
        }

        if (l2rec->chl[ip] == BAD_FLT)
            l1rec->flags[ip] |= CHLFAIL;

        else if (l2rec->chl[ip] > CHL_MAX || l2rec->chl[ip] < CHL_MIN)
            l1rec->flags[ip] |= CHLWARN;

        if (l2rec->num_iter[ip] > input->aer_iter_max) {
            l1rec->flags[ip] |= MAXAERITER;
            l1rec->flags[ip] |= ATMWARN;
        }

        if (input->absaer_opt > 0) {
            switch (input->absaer_opt) {
                case RHOWINDEX:
                    l2rec->aerindex[ip] = aerindex(l2rec, ip);
                    if (l2rec->aerindex[ip] < input->absaer)
                        l1rec->flags[ip] |= ABSAER;
                    break;
                case TOTEXINCTION:
                    if (l1rec->anc_aerosol) {
                        float aaod = l1rec->anc_aerosol->total_aerosol_ext[ip]
                                - l1rec->anc_aerosol->total_aerosol_scat[ip];
                        if (aaod > 0.01)
                            l1rec->flags[ip] |= ABSAER;
                    } else {
                        printf("Warning: Cannot apply GMAO total extinction absorbing aerosol test without\n");
                        printf("         anc_aerosol inputs defined\n");
                        printf("         Setting absaer_opt=0\n\n");
                        input->absaer_opt = 0;
                    }
                    break;
            }
        }
    }
}


#define Between(a,x,b)(x >= a && x <= b)

char isCoccolith(l2str *l2rec, int32_t ip) {
    static float firstCall = 1;
    static float c1, c2, c3, c4, c5, c6, c7, c8;
    static int ib443;
    static int ib510;
    static int ib555;
    static int ib670;

    int32_t nbands = l2rec->l1rec->l1file->nbands;

    if (firstCall) {
        firstCall = 0;
        float* wave = l2rec->l1rec->l1file->fwave;
        ib443 = windex(443., wave, nbands);
        ib510 = windex(510., wave, nbands);
        ib555 = windex(550., wave, nbands);
        ib670 = windex(670., wave, nbands);
        c1 = input->coccolith[0];
        c2 = input->coccolith[1];
        c3 = input->coccolith[2];
        c4 = input->coccolith[3];
        c5 = input->coccolith[4];
        c6 = input->coccolith[5];
        c7 = input->coccolith[6];
        c8 = input->coccolith[7];
    }

    float nLw443 = l2rec->nLw[ip * nbands + ib443];
    float nLw510 = l2rec->nLw[ip * nbands + ib510];
    float nLw555 = l2rec->nLw[ip * nbands + ib555];
    float La670 = l2rec->La [ip * nbands + ib670];

    if (nLw443 >= c1 &&
            nLw510 >= 0.0 &&
            nLw555 >= c2 &&
            La670 <= 1.1 &&
            Between(c3, nLw443 / nLw555, c4) &&
            Between(c5, nLw510 / nLw555, c6) &&
            Between(c7, nLw443 / nLw510, c8))

        return (1);
    else
        return (0);
}

char isTurbid(l2str *l2rec, int32_t ip) {

    static float Fo;
    static int firstCall = 1;
    static int ib670;

    if (firstCall) {
        firstCall = 0;
        ib670 = windex(670., l2rec->l1rec->l1file->fwave, l2rec->l1rec->l1file->nbands);
        if (l1_input->outband_opt >= 2)
            Fo = l2rec->l1rec->l1file->Fonom[ib670];
        else
            Fo = l2rec->l1rec->l1file->Fobar[ib670];
    }

    if (l2rec->nLw[ip * l2rec->l1rec->l1file->nbands + ib670] / Fo > 0.0012)
        return (1);
    else
        return (0);
}

char isOptShallow(l2str *l2rec, int32_t ip) {
    static float firstCall = 1;
    static int ib443;
    static int ib555;
    static int ib670;

    float* wave = l2rec->l1rec->l1file->fwave;

    if (firstCall) {
        firstCall = 0;
        ib443 = windex(443., wave, l2rec->l1rec->l1file->nbands);
        ib555 = windex(550., wave, l2rec->l1rec->l1file->nbands);
        ib670 = windex(670., wave, l2rec->l1rec->l1file->nbands);
    }

    //McKinna, L. & P.J. Werdell (2018). Approach for identifying optically 
    //shallow pixels when processing ocean-color imagery, Opt. Express, 26(22),
    //A915-A928, doi: 10.1364/OE.26.00A915
    
    float a670, a555, u670, u555, bbp670, bbp555, quasiC555, od555;
    
    float depth = 0. - l2rec->l1rec->dem[ip];
    
    float Rrs443 = l2rec->Rrs[ip * l2rec->l1rec->l1file->nbands + ib443];
    float Rrs555 = l2rec->Rrs[ip * l2rec->l1rec->l1file->nbands + ib555];
    float Rrs670 = l2rec->Rrs[ip * l2rec->l1rec->l1file->nbands + ib670];
    
    float rrsSub443 = Rrs443 / (0.52 + 1.7 * Rrs443);
    float rrsSub555 = Rrs555 / (0.52 + 1.7 * Rrs555);
    float rrsSub670 = Rrs670 / (0.52 + 1.7 * Rrs670);
    
    float aw670 = l2rec->l1rec->sw_a[ib670];
    float bbw670 = l2rec->l1rec->sw_bb[ib670];
    float bbw555 = l2rec->l1rec->sw_bb[ib555];

    //Check valid Rrs values
    if (Rrs443 < 0.0 || Rrs555 < 0.0 || Rrs670 < 0.0 ) {
        return 0;
    }
    
    //1. if water column greater than 40m, not flagged.
    if (depth > 40.0) {
        return 0;
    }
    
    //2. if water column shallower than 5m, flagged.
    if ( depth <= 5.0) {
        return 1;
    }     

    /*Use QAA methodology adopted by Barnes et al (2013) with reference*/ 
    /* wavelength set to 667nm*/ 
    
    a670 = aw670 + 0.07 * pow((rrsSub670 / rrsSub443), 1.1);
    u670 = (-0.089 + pow(0.089*0.089 + 4.0 * 0.125 * rrsSub670, 0.5)) / (2.0 * 0.125);
    u555 = (-0.089 + pow(0.089*0.089 + 4.0 * 0.125 * rrsSub555, 0.5)) / (2.0 * 0.125);

    bbp670 = ((u670 * a670) / (1.0 - u670)) - bbw670;

    bbp555 = bbp670 * pow((wave[ib670] / wave[ib555]), 1.03373);
    a555 = ((1.0 - u555)* (bbp555 + bbw555)) / u555;

    /*Estimate c555, the beam attenuation coefficient*/
    /*Assume the mean particulate backscatter ratio is 0.015*/
    quasiC555 = a555 + (bbp555/0.015) + (bbw555/0.5);

    /*Estimate the water column's optical depth at 547 nm*/
    od555 = quasiC555 * depth;
    
    /*Flag as optically shallow pixel if od(547) less than or equal to 20 */
    /*At this threshold, we assume the seafloor has a non-negligible effect*/ 
    /* on the the short wave spectral water-leaving reflectance signal*/
 
    if (od555 <= 20.0) {
        return 1;
    } else { 
        
        return 0;
    }

}

/* ---------------------------------------------------------------------------------------- */
/* aerindx() - compute aerosol index for absorbing aerosol test                             */

/* ---------------------------------------------------------------------------------------- */
float aerindex(l2str *l2rec, int32_t ip) {
    static int firstCall = 1;

    static int ib412;
    static int ib510;

    static int32_t mask = ATMFAIL | LAND | CHLFAIL;
    static int i510;
    static double poly_coeff[7];
    static double poly_coeff_412[7] = {0.73172938, -0.54565918, 0.20022312, -0.12036241, 0.11687968, -0.018732825, -0.0095574674};

    double poly_coeff_510[7] = {0.58543905, -0.013125745, -0.059568208, -0.016141849, 0.0035106655, 0.0012957265, 1.4235445e-05};
    double poly_coeff_531[7] = {0.51483158, 0.15893415, -0.051696975, -0.055474007, -0.0029635773, 0.0053882411, 0.0010106600};
    double poly_coeff_551[7] = {0.47507366, 0.25216739, -0.0096908094, -0.070882408, -0.012501495, 0.0061436085, 0.0015798359};
    double poly_coeff_555[7] = {0.43681192, 0.26663018, 0.016592559, -0.068132662, -0.015470602, 0.0051694309, 0.0015132129};
    //double poly_coeff_412_S[7] = {0.72884183, -0.54380414, 0.20225533, -0.12180914, 0.11508257, -0.017784535, -0.0095387145};
    //double poly_coeff_412_M[7] = {0.73461692, -0.54751422, 0.19819091, -0.11891567, 0.11867679, -0.019681115, -0.0095762203};
    float tLw_pred;
    float Lt_pred;
    float Lt_meas;
    float mu0;
    float index = 100.0, index510, index412;
    int32_t ipb, i;
    double logchl, lchl, nLw_510, nLw_412;


    /* load parameters */
    if (firstCall) {

        /* determine index of the bands for the absorbing aerosol analysis */
        ib412 = windex(412., l2rec->l1rec->l1file->fwave, l2rec->l1rec->l1file->nbands);
        ib510 = windex(510., l2rec->l1rec->l1file->fwave, l2rec->l1rec->l1file->nbands);
        i510 = l2rec->l1rec->l1file->iwave[ib510];

        switch (i510) {
        case 510:
            for (i = 0; i < 7; i++) poly_coeff[i] = poly_coeff_510[i];
            break;
        case 531:
            for (i = 0; i < 7; i++) poly_coeff[i] = poly_coeff_531[i];
            break;
        case 551:
            for (i = 0; i < 7; i++) poly_coeff[i] = poly_coeff_551[i];
            break;
        case 555:
            for (i = 0; i < 7; i++) poly_coeff[i] = poly_coeff_555[i];
            break;
        default:
            for (i = 0; i < 7; i++) poly_coeff[i] = poly_coeff_510[i];
            break;
        }
        /*
                switch (l2rec->sensorID) {
                    case SEAWIFS: 
                        for (i=0; i<7; i++) poly_coeff_412[i] = poly_coeff_412_S[i];
                        break;
                    case MODISA: 
                    case HMODISA: 
                        for (i=0; i<7; i++) poly_coeff_412[i] = poly_coeff_412_M[i];
                        break;
                    default:
                        break;
                }
         */
        firstCall = 0;

    }


    if ((l2rec->l1rec->flags[ip] & mask) != 0)
        return (index);

    if (l2rec->chl[ip] < 0.1)
        return (index);


    logchl = lchl = log10(l2rec->chl[ip]);
    nLw_510 = poly_coeff[0] + logchl * poly_coeff[1];
    nLw_412 = poly_coeff_412[0] + logchl * poly_coeff_412[1];
    for (i = 0; i < 5; i++) {
        logchl *= lchl;
        nLw_510 += poly_coeff[i + 2] * logchl;
        nLw_412 += poly_coeff_412[i + 2] * logchl;
    }

    /* Convert water-leaving radiances from band-centered to band-averaged. */
    /*
    if (input->outband_opt >= 2) {
          nLw_510 /= l2rec->f_outband[ip*l2rec->nbands + ib510];
          nLw_412 /= l2rec->f_outband[ip*l2rec->nbands + ib412];
    }
     */


    /* cos of solz */
    mu0 = cos(l2rec->l1rec->solz[ip] / RADEG);


    /* bring water-leaving radiance at 510nm to the TOA */
    ipb = ip * l2rec->l1rec->l1file->nbands + ib510;
    tLw_pred = (float) nLw_510 / l2rec->brdf[ipb]
            * l2rec->l1rec->t_sol[ipb]
            * l2rec->l1rec->t_sen[ipb]
            * l2rec->l1rec->tg_sol[ipb]
            * l2rec->l1rec->tg_sen[ipb]
            * l2rec->l1rec->polcor[ipb]
            * l2rec->l1rec->t_o2[ipb]
            * mu0
            * l2rec->l1rec->fsol;

    /* calculate TOA predicted radiance  */
    Lt_pred = tLw_pred
            + ((l2rec->l1rec->TLg[ipb]
            + l2rec->La[ipb])
            * l2rec->l1rec->t_o2[ipb]
            + l2rec->l1rec->Lr[ipb]
            + l2rec->l1rec->tLf[ipb])
            * l2rec->l1rec->polcor[ipb]
            * l2rec->l1rec->tg_sol[ipb]
            * l2rec->l1rec->tg_sen[ipb];

    /* get measured TOA radiance */
    Lt_meas = l2rec->l1rec->Lt[ipb];

    /* obtain percent difference between measured and predicted TOA radiance */
    index510 = 100.0 * (Lt_meas - Lt_pred) / Lt_meas;

    if (index510 >= 0.0)
        return (index);


    /* bring water-leaving radiance to the TOA */
    ipb = ip * l2rec->l1rec->l1file->nbands + ib412;
    tLw_pred = (float) nLw_412 / l2rec->brdf[ipb]
            * l2rec->l1rec->t_sol[ipb]
            * l2rec->l1rec->t_sen[ipb]
            * l2rec->l1rec->tg_sol[ipb]
            * l2rec->l1rec->tg_sen[ipb]
            * l2rec->l1rec->polcor[ipb]
            * l2rec->l1rec->t_o2[ipb]
            * mu0
            * l2rec->l1rec->fsol;

    /* calculate predicted TOA radiance  */
    Lt_pred = tLw_pred
            + ((l2rec->l1rec->TLg[ipb]
            + l2rec->La[ipb])
            * l2rec->l1rec->t_o2[ipb]
            + l2rec->l1rec->Lr[ipb]
            + l2rec->l1rec->tLf[ipb])
            * l2rec->l1rec->polcor[ipb]
            * l2rec->l1rec->tg_sol[ipb]
            * l2rec->l1rec->tg_sen[ipb];

    /* get measured TOA radiance */
    Lt_meas = l2rec->l1rec->Lt[ipb];

    /* obtain percent difference between measured and predicted TOA radiance */
    index412 = 100.0 * (Lt_meas - Lt_pred) / Lt_meas + 2.0;

    return (index412);

}


