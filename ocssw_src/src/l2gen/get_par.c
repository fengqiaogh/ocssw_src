/* --------------------------------------------------------------- */
/* get_par.c - computes photosynthetically active radiation        */
/*                                                                 */
/* Inputs:                                                         */
/*     l2rec - level-2 structure containing one complete scan      */
/*             after atmospheric correction.                       */
/* Outputs:                                                        */
/*     par   - Photosynthetically Active Radiation just above the  */
/*             surface from SeaWiFS level 1b instantaneous         */
/*             radiances at 412, 443, 490, 510, 555, and 670 nm.   */
/*                                                                 */
/* Algorithm Provided By: Robert Frouin,                           */
/*                        Scripps Institution of Oceanography      */
/*                                                                 */
/* Written By: B. A. Franz, SAIC GSC, SIMBIOS, 30 July 1999        */
/* Modified:                                                       */
/*  Fall 2013, Robert Lossing, SAIC - Adapted underlying code to C */
/* --------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>

#include "l12_proto.h"
#include "par_utils.h"
#include "smi_climatology.h"

static int *alloc_bindx(int nbands, int *bindx) {
    if ((bindx = (int *) calloc(nbands, sizeof(int))) == NULL) {
        printf("-E- : Error allocating memory to bindx in get_par\n");
        exit(FATAL_ERROR);
    }
    return bindx;
}

static int scan_processed = -1;
static float *parb, *parc, *par0, *para, *mu_est;

float *par_planar_a_inst, *par_planar_b_inst;

void get_par_scalar(l2str *l2rec, float par[]) {
    if (scan_processed != l2rec->l1rec->iscan || scan_processed == -1) {
        get_par2(l2rec, par);
    }

    for (size_t ip = 0; ip < l2rec->l1rec->npix; ip++) {
        calc_scalar_par_mean_cosine(l2rec, ip, para[ip], parc[ip], &par0[ip], &mu_est[ip]);
        par[ip] = par0[ip];
    }
};

void get_par_below_surface(l2str *l2rec, float par[]) {
    if (scan_processed != l2rec->l1rec->iscan || scan_processed == -1) {
        get_par2(l2rec, para);
    }
    for (size_t ip = 0; ip < l2rec->l1rec->npix; ip++)
        par[ip] = parb[ip];
};

void get_mu_cosine(l2str *l2rec, float mu[]) {
    get_par_scalar(l2rec, par0);
    for (size_t ip = 0; ip < l2rec->l1rec->npix; ip++) {
        mu[ip] = mu_est[ip];
        if (mu_est[ip] == BAD_FLT) {
            mu[ip] = BAD_FLT;
        }
    }
}

void get_par2(l2str *l2rec, float par[]) {
    if (scan_processed == l2rec->l1rec->iscan && scan_processed != -1) {
        for (size_t ip = 0; ip < l2rec->l1rec->npix; ip++)
            par[ip] = para[ip];
    }

    int32_t ip, ib, ipb, iw;
    float *Lt;
    float angst;
    float taua;

    static int32_t mask = SEAICE | LAND | HIGLINT | NAVFAIL;

    static float *lambda;
    static float *Fobar;
    static float *Taur;
    static float *kO3;
    static int *bindx;
    static int nwave;
    static int firstCall = TRUE;

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;

    if (firstCall) {
        firstCall = FALSE;
        parb = (float *) calloc(l1rec->npix, sizeof(float));
        parc = (float *) calloc(l1rec->npix, sizeof(float));
        par0 = (float *) calloc(l1rec->npix, sizeof(float));
        para = (float *) calloc(l1rec->npix, sizeof(float));
        mu_est = (float *) calloc(l1rec->npix, sizeof(float));
        par_planar_a_inst = (float *) calloc(l1rec->npix, sizeof(float));
        par_planar_b_inst = (float *) calloc(l1rec->npix, sizeof(float));
        int16_t year, day;
        double sec;
        unix2yds(l2rec->l1rec->scantime, &year, &day, &sec);
        if (strlen(input->cld_rad1) == 0) {
            printf("-Error-: rad1 file is not provided. Exiting ... \n ");
            exit(EXIT_FAILURE);
        }
        if (strlen(input->cld_rad2) == 0) {
            printf("-Error-: rad2 file is not provided. Exiting ... \n ");
            exit(EXIT_FAILURE);
        }
        if (strlen(input->cld_rad3) == 0) {
            printf("-Error-: rad3 file is not provided. Exiting ... \n ");
            exit(EXIT_FAILURE);
        }
        if (strlen(input->anc_aerosol1) == 0) {
            printf("-Error-: anc_aerosol1 file is not provided. Exiting ... \n ");
            exit(EXIT_FAILURE);
        }
        if (strlen(input->anc_aerosol2) == 0) {
            printf("-Error-: anc_aerosol2 file is not provided. Exiting ... \n ");
            exit(EXIT_FAILURE);
        }
        if (strlen(input->anc_aerosol3) == 0) {
            printf("-Error-: anc_aerosol3 file is not provided. Exiting ... \n ");
            exit(EXIT_FAILURE);
        }
        if (l2rec->l1rec->cld_rad == NULL) {
            printf("-Error-: cloud rad structure has not been allocated. Exiting ... \n ");
            exit(EXIT_FAILURE);
        }
        /* Initialize climatologies */
        smi_climatology_init(input->alphafile, day, ALPHA510);
        smi_climatology_init(input->tauafile, day, TAUA865);
        size_t total_waves = l1file->nbands;
        // size_t count = 0;
        size_t start_ib = total_waves;
        size_t end_ib = 0;
        switch (sensorId2InstrumentId(l1file->sensorID)) {
            case INSTRUMENT_MODIS:
                nwave = 3;
                bindx = alloc_bindx(nwave, bindx);
                bindx[0] = 2;
                bindx[1] = 6;
                bindx[2] = 7;
                break;
            case INSTRUMENT_SEAWIFS:
                nwave = 6;
                bindx = alloc_bindx(nwave, bindx);
                for (ib = 0; ib < nwave; ib++)
                    bindx[ib] = ib;
                break;
            case INSTRUMENT_VIIRS:
                nwave = 5;
                bindx = alloc_bindx(nwave, bindx);
                for (ib = 0; ib < nwave; ib++)
                    bindx[ib] = ib;
                break;
            case INSTRUMENT_MERIS:
                nwave = 7;
                bindx = alloc_bindx(nwave, bindx);
                for (ib = 0; ib < nwave; ib++)
                    bindx[ib] = ib;
                break;
            case INSTRUMENT_OCTS:
                nwave = 6;
                bindx = alloc_bindx(nwave, bindx);
                for (ib = 0; ib < nwave; ib++)
                    bindx[ib] = ib;
                break;
            case INSTRUMENT_OCI:

                for (ib = 0; ib < total_waves; ib++) {
                    if (l1file->fwave[ib] >= 400) {
                        start_ib = MIN(ib, start_ib);
                    }
                    if (l1file->fwave[ib] < 700) {
                        end_ib = MAX(ib, end_ib);
                    }
                }
                nwave = 20; // end_ib - start_ib + 1;  // to be changed
                bindx = alloc_bindx(nwave, bindx);
                size_t band_step = (end_ib - start_ib) / (nwave - 1);
                //printf("Wavelenght band is %d ", (int)(end_ib - start_ib));
                for (size_t ib = 0; ib < nwave; ib++) {
                    bindx[ib] = ib * band_step + start_ib;
                }

                break;
            default:
                printf("PAR not supported for this sensor (%d).\n", l1file->sensorID);
                exit(1);
                break;
        }
        if ((lambda = (float *) calloc(nwave, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to lambda in get_par\n");
            exit(FATAL_ERROR);
        }
        if ((Fobar = (float *) calloc(nwave, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to Fobar in get_par\n");
            exit(FATAL_ERROR);
        }
        if ((Taur = (float *) calloc(nwave, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to Taur in get_par\n");
            exit(FATAL_ERROR);
        }
        if ((kO3 = (float *) calloc(nwave, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to kO3 in get_par\n");
            exit(FATAL_ERROR);
        }

        /* Get band-pass dependent quantities */
        // consider passing in an array of band indexes to index
        // l2rec->(lambda,k03,Fobar,Taur) in calc_par
        for (iw = 0; iw < nwave; iw++) {
            ib = bindx[iw];
            lambda[iw] = l1file->fwave[ib];
            kO3[iw] = l1file->k_oz[ib];
            Fobar[iw] = l1file->Fobar[ib];
            Taur[iw] = l1file->Tau_r[ib];
        }
    }
    if ((Lt = (float *) calloc(nwave, sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to Lt in get_par\n");
        exit(FATAL_ERROR);
    }
    for (ip = 0; ip < l1rec->npix; ip++) {
        /* Grab radiances for this pixel */
        for (ib = 0; ib < nwave; ib++) {
            ipb = ip * l1file->nbands + bindx[ib];
            Lt[ib] = l1rec->Lt[ipb];
        }

        /* Skip pixel if masked */
        if (Lt[0] <= 0.0 || (l1rec->flags[ip] & mask) != 0 || l1rec->solz[ip] > 90.0) {
            par[ip] = BAD_FLT;
            para[ip] = par[ip];
            parb[ip] = BAD_FLT;
            mu_est[ip] = BAD_FLT;
            parc[ip] = BAD_FLT;
            par_planar_a_inst[ip] = BAD_FLT;
            par_planar_b_inst[ip] = BAD_FLT;
            l1rec->flags[ip] |= PRODFAIL;
            continue;
        }

        /* Get angstrom and AOT from climatology */
        angst = smi_climatology(l1rec->lon[ip], l1rec->lat[ip],
                                ALPHA510);  // ALPHA550
        taua = smi_climatology(l1rec->lon[ip], l1rec->lat[ip],
                               TAUA865);  // TAUA550
        switch (sensorId2InstrumentId(l1file->sensorID)) {
            case INSTRUMENT_SEAWIFS:
            case INSTRUMENT_MERIS:
            case INSTRUMENT_VIIRS:
            case INSTRUMENT_OCTS:
            case INSTRUMENT_OCI:
            case INSTRUMENT_MODIS:
                par[ip] = calc_par_impl_of_2023(l2rec, ip, nwave, Lt, taua, angst, lambda, Fobar, kO3, Taur,
                                                parb, parc);
                break;
            default:
                printf("PAR not supported for this sensor (%d).\n", l1file->sensorID);
                exit(1);
                break;
        }
        /* Convert to E/D/m^2 */
        if (par[ip] != BAD_FLT) {
            par[ip] *= EINSTEIN;
            parb[ip] *= EINSTEIN;
            parc[ip] *= EINSTEIN;
            par_planar_a_inst[ip] *= EINSTEIN;
            par_planar_b_inst[ip] *= EINSTEIN;
        } else {
            par[ip] = BAD_FLT;
            parb[ip] = BAD_FLT;
        }
        para[ip] = par[ip];
    }
    free(Lt);
    // process the scan
    scan_processed = l2rec->l1rec->iscan;
}

static int32_t ini_obs = 0;
static float observed_time;
static int32_t index_obs;

void get_taucld(l2str *l2rec, float taucld[]) {
    if (ini_obs == 0) {
        ini_obs = 1;
        int16_t year, month, mday;
        double sec;
        unix2ymds(l2rec->l1rec->scantime, &year, &month, &mday, &sec);
        observed_time = sec / 3600;
        size_t ntimes = l2rec->l1rec->cld_rad->ntimes;
        float min_diff = 73;
        for (size_t it = 0; it < ntimes; it++) {
            float diff = fabs(l2rec->l1rec->cld_rad->timecldrange[it] - observed_time);
            if (min_diff > diff) {
                min_diff = diff;
                index_obs = it;
            };
        }
    }
    size_t npix = l2rec->l1rec->npix;
    for (size_t ip = 0; ip < npix; ip++) {
        taucld[ip] = l2rec->l1rec->cld_rad->taucld[ip][index_obs];
    }
};

void get_clfr(l2str *l2rec, float clfr[]) {
    if (ini_obs == 0) {
        ini_obs = 1;
        int16_t year, month, mday;
        double sec;
        unix2ymds(l2rec->l1rec->scantime, &year, &month, &mday, &sec);
        observed_time = sec / 3600;
        size_t ntimes = l2rec->l1rec->cld_rad->ntimes;
        float min_diff = 73;
        for (size_t it = 0; it < ntimes; it++) {
            float diff = fabs(l2rec->l1rec->cld_rad->timecldrange[it] - observed_time);
            if (min_diff > diff) {
                min_diff = diff;
                index_obs = it;
            };
        }
    }
    size_t npix = l2rec->l1rec->npix;
    for (size_t ip = 0; ip < npix; ip++) {
        clfr[ip] = l2rec->l1rec->cld_rad->cfcld[ip][index_obs];
    }
};

void get_par(l2str *l2rec, float par[]) {
    int32_t ip, ib, ipb, iw;

    float *Lt;

    float angst;
    float taua;

    static int32_t mask = SEAICE | LAND | HIGLINT | NAVFAIL;

    static float *lambda;
    static float *Fobar;
    static float *Taur;
    static float *kO3;
    static int *bindx;
    static int nwave;
    static int firstCall = TRUE;

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;

    if (firstCall) {
        firstCall = FALSE;

        int16_t year, day;
        double sec;
        unix2yds(l2rec->l1rec->scantime, &year, &day, &sec);

        /* Initialize climatologies */
        smi_climatology_init(input->alphafile, day, ALPHA510);
        smi_climatology_init(input->tauafile, day, TAUA865);

        switch (sensorId2InstrumentId(l1file->sensorID)) {
            case INSTRUMENT_MODIS:
                nwave = 3;
                bindx = alloc_bindx(nwave, bindx);
                bindx[0] = 2;
                bindx[1] = 6;
                bindx[2] = 7;
                break;
            case INSTRUMENT_SEAWIFS:
                nwave = 6;
                bindx = alloc_bindx(nwave, bindx);
                for (ib = 0; ib < nwave; ib++)
                    bindx[ib] = ib;
                break;
            case INSTRUMENT_VIIRS:
                nwave = 5;
                bindx = alloc_bindx(nwave, bindx);
                for (ib = 0; ib < nwave; ib++)
                    bindx[ib] = ib;
                break;
            case INSTRUMENT_MERIS:
                nwave = 7;
                bindx = alloc_bindx(nwave, bindx);
                for (ib = 0; ib < nwave; ib++)
                    bindx[ib] = ib;
                break;
            case INSTRUMENT_OCTS:
                nwave = 6;
                bindx = alloc_bindx(nwave, bindx);
                for (ib = 0; ib < nwave; ib++)
                    bindx[ib] = ib;
                break;
            default:
                printf("PAR not supported for this sensor (%d).\n", l1file->sensorID);
                exit(1);
                break;
        }
        if ((lambda = (float *) calloc(nwave, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to lambda in get_par\n");
            exit(FATAL_ERROR);
        }
        if ((Fobar = (float *) calloc(nwave, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to Fobar in get_par\n");
            exit(FATAL_ERROR);
        }
        if ((Taur = (float *) calloc(nwave, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to Taur in get_par\n");
            exit(FATAL_ERROR);
        }
        if ((kO3 = (float *) calloc(nwave, sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to kO3 in get_par\n");
            exit(FATAL_ERROR);
        }

        /* Get band-pass dependent quantities */
        // consider passing in an array of band indexes to index l2rec->(lambda,k03,Fobar,Taur) in calc_par
        for (iw = 0; iw < nwave; iw++) {
            ib = bindx[iw];
            lambda[iw] = l1file->fwave[ib];
            kO3[iw] = l1file->k_oz[ib];
            Fobar[iw] = l1file->Fobar[ib];
            Taur[iw] = l1file->Tau_r[ib];
        }
    }
    if ((Lt = (float *) calloc(nwave, sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to Lt in get_par\n");
        exit(FATAL_ERROR);
    }

    for (ip = 0; ip < l1rec->npix; ip++) {
        /* Grab radiances for this pixel */
        for (ib = 0; ib < nwave; ib++) {
            ipb = ip * l1file->nbands + bindx[ib];
            Lt[ib] = l1rec->Lt[ipb];
        }

        /* Skip pixel if masked */
        if (Lt[0] <= 0.0 || (l1rec->flags[ip] & mask) != 0 || l1rec->solz[ip] > 90.0) {
            par[ip] = BAD_FLT;
            l1rec->flags[ip] |= PRODFAIL;
            continue;
        }

        /* Get angstrom and AOT from climatology */
        angst = smi_climatology(l1rec->lon[ip], l1rec->lat[ip], ALPHA510);
        taua = smi_climatology(l1rec->lon[ip], l1rec->lat[ip], TAUA865);

        switch (sensorId2InstrumentId(l1file->sensorID)) {
            case INSTRUMENT_SEAWIFS:
            case INSTRUMENT_MODIS:
            case INSTRUMENT_MERIS:
            case INSTRUMENT_VIIRS:
            case INSTRUMENT_OCTS:

                par[ip] = calc_par(l2rec, ip, nwave, Lt, taua, angst, lambda, Fobar, kO3, Taur);
                break;

            default:
                printf("PAR not supported for this sensor (%d).\n", l1file->sensorID);
                exit(1);
                break;
        }
        /* Convert to E/D/m^2 */
        if (par[ip] != BAD_FLT)
            par[ip] *= 1.193;
    }
    free(Lt);
}


void get_ipar2(l2str *l2rec, float ipar[]) {
    if (scan_processed != l2rec->l1rec->iscan || scan_processed == -1) {
        get_par2(l2rec, para);
    }
    for (size_t ip = 0; ip < l2rec->l1rec->npix; ip++) {
        if(par_planar_a_inst[ip]!=BAD_FLT)
            ipar[ip] = par_planar_a_inst[ip] * 1e6;
    }
}

void get_ipar_below_surface(l2str *l2rec, float ipar[]) {
    if (scan_processed != l2rec->l1rec->iscan || scan_processed == -1) {
        get_par2(l2rec, para);
    }
    for (size_t ip = 0; ip < l2rec->l1rec->npix; ip++) {
        if(par_planar_b_inst[ip]!=BAD_FLT)
            ipar[ip] = par_planar_b_inst[ip] * 1e6;
    }
}


void get_ipar_scalar(l2str *l2rec, float ipar[]) {
    if (scan_processed != l2rec->l1rec->iscan || scan_processed == -1) {
        get_par2(l2rec, para);
    }
    for (size_t ip = 0; ip < l2rec->l1rec->npix; ip++) {
        if (para[ip] == BAD_FLT) {
            ipar[ip] = BAD_FLT;
            continue;
        }
        calc_scalar_inst_par(l2rec, ip, par_planar_a_inst[ip], ipar + ip);
        ipar[ip]*=1e6;
    }
}
/*
 Subject:
 PAR routine
 Date:
 Fri, 23 Jul 1999 08:55:29 -0700
 From:
 Robert Frouin <rfrouin@ucsd.edu>
 To:
 chuck@seawifs, gfargion@simbios, gene@seawifs, wang@simbios, franz@seawifs
 CC:
 jmcpherson@ucsd.edu




 Greetings:

 A routine to compute daily PAR from SeaWiFS level 1b radiances is available
 at the following address: http://genius.ucsd.edu/~john/SeaWiFS_dir/ under
 the rubrique "PAR subroutine and test program".

 The routine requires as input year, month, day, time, latitude, longitude,
 SeaWiFS radiances in the first 6 spectral bands, solar zenith angle,
 viewing zenith angle, relative azimuth angle, aerosol optical thickness at
 865 nm, Angstrom coefficient, ozone amount, and surface pressure. Routine
 output is daily PAR.

 Thus a daily PAR value is computed for each instantaneous SeaWiFS
 observation, clear or cloudy. Diurnal variations are taken into account
 statistically. The algorithm is described succintly in the routine.

 During our discussion at GSFC in June, a first routine was supposed to be
 developed to provide a normalized cloud/surface albedo and then a second
 routine to compute daily PAR from the normalized albedo, and the second
 routine was to be applied when binning to the 9 km resolution. Now daily
 PAR is obtained using a single routine, which is more convenient.

 The binning to the 9 km resolution should be done as follows. First,
 weight-average the daily PAR estimates obtained from all SeaWiFS
 observations during the same day at each location (there might be several
 SeaWiFS observations of a surface target during the same day). The weight
 is the cosine of the sun zenith angle for the SeaWiFS observation. That is:

 PAR_avg = sum{cos[tetas(i)]*PAR(i)}/sum{cos[tetas(i)]}

 Second, simply average the values at all the locations within the 9 km bins.

 The routine requires aerosol data, ozone amount, surface pressure. If these
 parameters are missing (-999 or less), default values are used. If Eric
 Vermote's aerosol climatology (or Menghua's) is not available yet, please
 use default values for tests.

 At this time, the statistical diurnal function does not depend on latitude,
 longitude, and date, but will depend on these parameters in the second
 version of the code. Creating a date and location dependent function
 requires analysing several years of ISCCP data. We have the data, but a
 couple of weeks is needed to accomplish the task.

 I will present the algorithm at the SeaWiFS atmospheric correction meeting
 next week, and prepare a detailed documentation.

 Best, Robert.


 Robert Frouin
 Scripps Institution of Oceanography
 University of California San Diego
 9500 Gilman Drive
 La Jolla, CA 92093-0221
 Voice Tel.: 619/534-6243
 Fax Tel.: 619/534-7452


 */
