/* =========================================================== */
/* Module loadl1.c                                             */
/*                                                             */
/* Functions to fill a level-1b file with precomputed          */
/* atmospheric and masking data.                               */
/*                                                             */
/* Note: due to filtering capabilities of MSl12, it can not be */
/* assumed that the same l1rec pointer will be passed in later */
/* calls to this function, so all fields must be reloaded. In  */
/* addition, it is possible for MSl12 to process a L1B file    */
/* containing data from different time periods, so earth-sun   */
/* distance can not be assumed to be constant.                 */
/*                                                             */
/* Written By:                                                 */
/*                                                             */
/*     B. A. Franz                                             */
/*     SAIC General Sciences Corp.                             */
/*     NASA/SIMBIOS Project                                    */
/*     February 1999                                           */
/*                                                             */
/* Perturbed By:                                               */
/* E.Karakoylu (SAIC)                                          */
/* Summer 2015                                                 */
/* =========================================================== */

#include <stdio.h>
#include "l12_proto.h"
#include "smi_climatology.h"
#include "read_pixel_anc_file.h"
#include "anc_acq.h"
#include <libnav.h>
#include <xcal.h>
#include <smile.h>


int loadl1(filehandle *l1file, l1str *l1rec) {
    static double radeg = RADEG;
    static int32_t sensorID = -999;
    static float *aw;
    static float *bbw;
    int navfail_cnt = 0;

    int32_t ip, ipb, ib, iw, ix;
    double esdist;
    int32_t nbands = l1rec->l1file->nbands;

    double *rvs;
    double temp;
    int16_t year, day;
    double sec;
    unix2yds(l1rec->scantime, &year, &day, &sec);

    if (sensorID != l1file->sensorID) {

        sensorID = l1file->sensorID;
        aw = l1file->aw;
        bbw = l1file->bbw;


        printf("Loading land mask information from %s\n", input->land);
        if (land_mask_init() != 0) {
            printf("-E- %s : Unable to initialize land mask\n", __FILE__);
            exit(1);
        }

        printf("Loading DEM information from %s\n", input->demfile);
        if (input->dem_auxfile[0])
            printf("Loading auxiliary elevation file from %s\n", input->dem_auxfile);
        if (dem_init() != 0) {
            printf("-E- %s : Unable to initialize DEM \n", __FILE__);
            exit(1);
        }

        printf("Loading ice mask file from %s\n", input->icefile);
        if (ice_mask_init(input->icefile, (int) year,
                (int) day, input->ice_threshold) != 0) {
            printf("-E- %s : Unable to initialize ice mask\n", __FILE__);
            exit(1);
        }

    }

    /* Get correction for Earth-Sun distance and apply to Fo  */
    int32_t yr = (int32_t) year;
    int32_t dy = (int32_t) day;
    int32_t ms = (int32_t) (sec * 1.e3);
    esdist = esdist_(&yr, &dy, &ms);
    l1rec->fsol = pow(1.0 / esdist, 2);

    for (iw = 0; iw < nbands; iw++) {
        l1rec->Fo[iw] = l1file->Fobar[iw] * l1rec->fsol;
    }

    /* Apply vicarious cross-calibration gains */

    for (ix = 0; ix < l1_input->xcal_nwave; ix++) {
        if ((l1_input->xcal_opt[ix] & XCALRVS) != 0) {
            if ((ib = bindex_get(l1_input->xcal_wave[ix])) < 0) {
                printf("-E- %s line %d: xcal wavelength %f does not match sensor\n",
                        __FILE__, __LINE__, l1_input->xcal_wave[ix]);
                exit(1);
            };
            rvs = get_xcal(l1rec, XRVS, l1file->iwave[ib]);

            for (ip = 0; ip < l1rec->npix; ip++) {
                ipb = ip * nbands + ib;
                if (rvs == 0x0) {
                    l1rec->Lt[ipb] = -999.0;
                    continue;
                }
                if (l1rec->Lt[ipb] > 0.0 && l1rec->Lt[ipb] < 1000.0)
                    l1rec->Lt[ipb] /= rvs[ip];
            }
        }
    }

    for (ip = 0; ip < l1rec->npix; ip++) {
        /* Apply vicarious calibration */

        for (iw = 0; iw < nbands; iw++) {
            ipb = ip * nbands + iw;

            if (l1rec->Lt[ipb] > 0.0 && l1rec->Lt[ipb] < 1000.0) {
                l1rec->Lt[ipb] *= l1_input->gain[iw];
                l1rec->Lt[ipb] += l1_input->offset [iw];
            }
        }

        /*** Geolocation-based lookups ***/
        if (!l1rec->navfail[ip]) {

            /* Enforce longitude convention */
            if (l1rec->lon[ip] < -180.)
                l1rec->lon[ip] += 360.0;

            else if (l1rec->lon[ip] > 180.0)
                l1rec->lon[ip] -= 360.0;

           l1rec->dem[ip] = get_dem(l1rec->lat[ip], l1rec->lon[ip]);
            /* Get terrain height */
            if (input->proc_land) {
                if (get_height(l1rec, ip,
                        l1file->terrain_corrected) != 0) {
                    printf("-E- %s line %d: Error getting terrain height.\n",
                            __FILE__, __LINE__);
                    return (1);
                }
            } else
                l1rec->height[ip] = 0.0;

            /* Set land, bathymetry and ice flags */
            if (input->proc_ocean != 2 &&
                    l1file->format != FT_L3BIN ) {
                if( land_bath_mask(l1rec,ip) )
                    printf("-E- %s line %d: Error setting land and bathymetry masks.\n",
                            __FILE__, __LINE__);
            }

            if (!l1rec->land[ip] &&
                    ice_mask(l1rec->lon[ip], l1rec->lat[ip]) != 0) {
                l1rec->ice[ip] = ON;
            }
        } else {
            navfail_cnt++;
        }
        /*** end Geolocation-based lookups ***/

        /* Set sea surface temperature and salinity, and seawater optical properties */
        for (iw = 0; iw < nbands; iw++) {
            ipb = ip * nbands + iw;
            l1rec->sw_n [ipb] = 1.334;
            // center band
            l1rec->sw_a [ipb] = aw_spectra(l1file->fwave[iw], BANDW);
            l1rec->sw_bb[ipb] = bbw_spectra(l1file->fwave[iw], BANDW);
            // band-averaged
            l1rec->sw_a_avg [ipb] = aw [iw];
            l1rec->sw_bb_avg[ipb] = bbw[iw];
        }

        l1rec->sstref[ip] = BAD_FLT;
        l1rec->sssref[ip] = BAD_FLT;

        if (!l1rec->land[ip]) {
            float bbw_fac;
            l1rec->sstref[ip] = get_sstref(input->sstreftype, input->sstfile, l1rec, ip);
            l1rec->sssref[ip] = get_sssref(input->sssfile, l1rec->lon[ip], l1rec->lat[ip], (int) day);
            if (l1rec->sstref[ip] > BAD_FLT && l1rec->sssref[ip] > BAD_FLT && input->seawater_opt > 0) {
                for (iw = 0; iw < nbands; iw++) {
                    ipb = ip * nbands + iw;
                    l1rec->sw_n [ipb] = seawater_nsw(l1file->fwave[iw], l1rec->sstref[ip], l1rec->sssref[ip], NULL);
                    // scale bbw based on ratio of center-band model results for actual sea state versus
                    // conditions of Morel measurements used to derive center and band-averaged bbw
                    bbw_fac = seawater_bb(l1file->fwave[iw], l1rec->sstref[ip], l1rec->sssref[ip], 0.039) / seawater_bb(l1file->fwave[iw], 20.0, 38.4, 0.039);
                    l1rec->sw_bb[ipb] *= bbw_fac;
                    l1rec->sw_bb_avg[ipb] *= bbw_fac;
                }
            }
        }
        seawater_set(l1rec);

        /* Compute relative azimuth */
        /* CLASS AVHRR files contain relative azimuth so don't overwrite it */
        if (sensorID != AVHRR) {
            l1rec->delphi[ip] = l1rec->sena[ip] - 180.0 - l1rec->sola[ip];
        }
        if (l1rec->delphi[ip] < -180.)
            l1rec->delphi[ip] += 360.0;
        else if (l1rec->delphi[ip] > 180.0)
            l1rec->delphi[ip] -= 360.0;

        /* Precompute frequently used trig relations */
        l1rec->csolz[ip] = cos(l1rec->solz[ip] / radeg);
        l1rec->csenz[ip] = cos(l1rec->senz[ip] / radeg);

        /* Scattering angle */
        temp = sqrt((1.0 - l1rec->csenz[ip] * l1rec->csenz[ip])*(1.0 - l1rec->csolz[ip] * l1rec->csolz[ip]))
                * cos(l1rec->delphi[ip] / radeg);
        l1rec->scattang[ip] = acos(MAX(-l1rec->csenz[ip] * l1rec->csolz[ip] + temp, -1.0)) * radeg;

    }
    /* get derived quantities for band-dependent view angles */
    if (l1rec->geom_per_band != NULL)
        geom_per_band_deriv(l1rec);

    /* add ancillary data */
    // the navfail_cnt test is a kludge to prevent processing failures for scans that are entirely invalid.
    if (navfail_cnt != l1rec->npix) {
        if (setanc(l1rec) != 0)
            return (1);
    }  else {
        //  allocate the profile and cloud  storage if files are in use
        if( input->anc_profile1[0] != 0 )
            if (init_anc_add(l1rec) != 0)
                return 1;
        if( strlen(input->sfc_albedo ) )
            if (init_cld_dat(l1rec) != 0)
                return 1;
    }
    
    if (input->pixel_anc_file[0])
        read_pixel_anc_file(input->pixel_anc_file, l1rec);

    if (input->windspeed > -999)
        for (ip = 0; ip < l1rec->npix; ip++)
            l1rec->ws[ip] = input->windspeed;
    if (input->windangle > -999)
        for (ip = 0; ip < l1rec->npix; ip++)
            l1rec->wd[ip] = input->windangle;
    if (input->pressure > -999)
        for (ip = 0; ip < l1rec->npix; ip++)
            l1rec->pr[ip] = input->pressure;
    if (input->ozone > -999)
        for (ip = 0; ip < l1rec->npix; ip++)
            l1rec->oz[ip] = input->ozone;
    if (input->watervapor > -999)
        for (ip = 0; ip < l1rec->npix; ip++)
            l1rec->wv[ip] = input->watervapor;
    if (input->relhumid > -999)
        for (ip = 0; ip < l1rec->npix; ip++)
            l1rec->rh[ip] = input->relhumid;

    /* SWIM bathymetry */
    // for (ip = 0; ip < l1rec->npix; ip++)
    //     if (!l1rec->navfail[ip])
    //         l1rec->dem[ip] = get_elev(l1rec->lat[ip], l1rec->lon[ip]);

    /* add atmospheric cnomponents that do not depend on Lt */
    for (ip = 0; ip < l1rec->npix; ip++) {

        if (l1rec->is_l2){
            // l1rec->Lt is actually Rrs, so skip atmo, etc
        } else {
        /* ------------------------------------------------ */
        /* Ocean processing                                 */
        /* ------------------------------------------------ */
            if ((input->proc_ocean != 0) && !l1rec->land[ip] && !l1rec->navfail[ip]) {

                atmocor1(l1rec, ip);

                /* set polarization correction */
                polcor(l1rec, ip);

                /* add surface reflectance */
                get_rhos(l1rec, ip);
            }
                /* ------------------------------------------------ */
                /* Land Processing                                  */
                /* ------------------------------------------------ */
            else if (input->proc_land && l1rec->land[ip] && !l1rec->navfail[ip]) {
                atmocor1_land(l1rec, ip);
                get_rhos(l1rec, ip);
            }
                /* ------------------------------------------------ */
                /* General Processing                               */
                /* ------------------------------------------------ */
            // else {
            //     for (ib = 0; ib < nbands; ib++) {
            //         ipb = ip * nbands + ib;
            //         l1rec->Lr [ipb] = 0.0;
            //         l1rec->tg [ipb] = 1.0;
            //         l1rec->t_sol [ipb] = 1.0;
            //         l1rec->t_sen [ipb] = 1.0;
            //         l1rec->tg_sol[ipb] = 1.0;
            //         l1rec->tg_sen[ipb] = 1.0;
            //         l1rec->t_o2 [ipb] = 1.0;
            //         l1rec->t_h2o [ipb] = 1.0;
            //         l1rec->polcor [ipb] = 1.0;
            //     }
            // }
        }

    }
    /* set masks and flags */
    if (setflags(l1rec) == 0)
        return (1);

    setflagbits_l1(1, l1rec, -1);

    return (0);
}
