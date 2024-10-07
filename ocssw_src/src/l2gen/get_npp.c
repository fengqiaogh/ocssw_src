#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"
#include "get_npp.h"


double betasw_ZHH2009(double, double, double);
float opp_cafe( float PAR,
          float chl,
          float mld,
          float lat,
          int    yd,
          float aph_443,
          float adg_443,
          float bbp_443,
          float bbp_s,
          float sst );

/**
 * Calculate the Primary Productivity using 1 of 3 algorithms
 *  1) Behrenfeld-Falkowski: (BeFa)
 *  2) Behrenfeld-Falkowski  algorithm, but
             modifies the pb_opt function after Eppley (as
             implemented by Antoine and Morel): (Eppley)
    3) Primary Productivity using a chl:Carbon ratio.
                     This is a spectrally resolved version of the cbpm, using nine separate
                     wavelengths: (Updated CbPM)

      Adapted from Oregon State U.
      http://www.science.oregonstate.edu/ocean.productivity/

      by R. Healy at NASA
      January 2015
 */
void get_npp(l2str *l2rec, int prodnum, float prod[]) {
    static int32_t ib440;
    int32_t ip, ipb;
    float *kd, *par;
    static float32 *parin, *lat, *lon;
    static int nlat, nlon;
    float sst, chl, trise, tset, mld, bbp, zno3, irr;
    static int firstCall = 1, havefile;
    char *parfile = input->parfile;

    l1str *l1rec = l2rec->l1rec;

    int16_t year, day;
    double sec;
    unix2yds(l2rec->l1rec->scantime, &year, &day, &sec);

    char sdsname[H4_MAX_NC_NAME];
    int ncid, ndims;
    int32 sds_id;
    int status;
    nc_type rh_type; /* variable type */
    int dimids[H4_MAX_VAR_DIMS]; /* dimension IDs */
    int natts; /* number of attributes */
    size_t length;

    int nbands=l1rec->l1file->nbands;
    float *aph,*adg,*bbp_temp,*bbp_s;
    static int32_t ib443;

    if (firstCall && strcmp(parfile, "") != 0) {
        /* try netCDF first */
        if (nc_open(parfile, NC_NOWRITE, &ncid) == NC_NOERR) {


            strcpy(sdsname, "lat");

            status = nc_inq_varid(ncid, sdsname, &sds_id);
            if (status != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, sdsname, parfile);
                exit(1);
            }

            status = nc_inq_var(ncid, sds_id, 0, &rh_type, &ndims, dimids,
                    &natts);

            if (nc_inq_dimlen(ncid, dimids[0], &length) != NC_NOERR) {
                char name[H4_MAX_NC_NAME];
                nc_inq_dim(ncid, dimids[0], name, &length);
                fprintf(stderr,
                        "-E- %s line %d: could not get size of dimension \"%s\" in netCDF File.\n",
                        __FILE__, __LINE__, name);
                exit(1);
            }

            nlat = length;

            if ((lat = (float *) calloc(nlat, sizeof (float))) == NULL) {
                printf("-E- : Error allocating memory to tindx\n");
                exit(FATAL_ERROR);
            }

            if (nc_get_var(ncid, sds_id, lat) != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, sdsname, parfile);
                exit(1);
            }

            strcpy(sdsname, "lon");

            status = nc_inq_varid(ncid, sdsname, &sds_id);
            if (status != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, sdsname, parfile);
                exit(1);
            }

            status = nc_inq_var(ncid, sds_id, 0, &rh_type, &ndims, dimids,
                    &natts);
            if (status != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, sdsname, parfile);
                exit(1);
            }

            if (nc_inq_dimlen(ncid, dimids[0], &length) != NC_NOERR) {
                char name[H4_MAX_NC_NAME];
                nc_inq_dim(ncid, dimids[0], name, &length);
                fprintf(stderr,
                        "-E- %s line %d: could not get size of dimension \"%s\" in netCDF File.\n",
                        __FILE__, __LINE__, name);
                exit(1);
            }

            nlon = length;

            if ((lon = (float *) calloc(nlon, sizeof (float))) == NULL) {
                printf("-E- : Error allocating memory to tindx\n");
                exit(FATAL_ERROR);
            }

            if (nc_get_var(ncid, sds_id, lon) != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, sdsname, parfile);
                exit(1);
            }

            strcpy(sdsname, "par");

            status = nc_inq_varid(ncid, sdsname, &sds_id);
            if (status != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, sdsname, parfile);
                exit(1);
            }

            status = nc_inq_var(ncid, sds_id, 0, &rh_type, &ndims, dimids,
                    &natts);
            if (status != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, sdsname, parfile);
                exit(1);
            }
            if (ndims != 2) {
                fprintf(stderr, "-E- %s line %d:  Wrong number of dimensions for %s.  Need 2 got %d.\n",
                        __FILE__, __LINE__, sdsname, ndims);
                exit(1);

            }
            printf("PARFILE: %s nlat=%d nlon=%d\n", parfile, nlat, nlon);

            if ((parin = (float *) calloc(nlat * nlon, sizeof (float))) == NULL) {
                printf("-E- : Error allocating memory to tindx\n");
                exit(FATAL_ERROR);
            }

            if (nc_get_var_float(ncid, sds_id, parin) != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, sdsname, parfile);
                exit(1);
            }
            // apply  scale_factor and add_offset
            double add_offset = 0.0, scale_factor = 1.0;
            nc_get_att_double(ncid, sds_id, "add_offset", &add_offset);
            nc_get_att_double(ncid, sds_id, "scale_factor", &scale_factor);
            int i;
            for (i=0; i<nlat*nlon; i++){
                parin[i] *= scale_factor;
                parin[i] += add_offset;
            }

            havefile = 1;

        } else {
            fprintf(stderr, "-E- %s line %d:  Error opening parfile = %s.\n",
                    __FILE__, __LINE__, parfile);
            exit(1);
        }
    }

    if ((kd = (float *) calloc(l1rec->npix, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to par in get_opp\n");
        exit(FATAL_ERROR);
    }
    if ((par = (float *) calloc(l1rec->npix, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to par in get_opp\n");
        exit(FATAL_ERROR);
    }

    if (prodnum == CAT_npp_cbpm2) {
        if (input->iop_opt == IOPNONE) {
            printf("IOP-based Kd_lee product requires iop model selection (iop_opt).  ");
            printf("Using default model.\n");
            input->iop_opt = IOPDEFAULT;
            get_iops(l2rec, input->iop_opt);
        }
        Kd490_obpg(l2rec, kd);
    }

    /* get irradiance E/D/m^2*/

    if (havefile)
        get_par_clim(parin, lat, lon, nlat, nlon, l1rec->lat, l1rec->lon, l1rec->npix, par);
    else
        get_par(l2rec, par);

    for (ip = 0; ip < l1rec->npix; ip++) {
        sst = l1rec->sstref[ip];
        chl = l2rec->chl[ip];
        /*  Get the rise and set time for this day to get the number of hours of daylight:*/
        triseset(day, l1rec->lon[ip], l1rec->lat[ip], &trise, &tset);

        switch (prodnum) {
        case CAT_npp_mld:
            prod[ip] = get_mld(input->mldfile, l1rec->lon[ip], l1rec->lat[ip], day);
            break;
        case CAT_npp_zno3:
            prod[ip] = get_zno3(l1rec->lon[ip], l1rec->lat[ip], day);
            break;
        case CAT_npp_bbp:
            prod[ip] = l2rec->bb[l1rec->l1file->nbands * ip + ib440];
            break;
        case CAT_npp_par:
            prod[ip] = par[ip];
            break;
        default:

            if (par[ip] != BAD_FLT && chl > 0) {

                switch (prodnum) {
                case CAT_npp_vgpm:
                    if (sst > -2)
                        prod[ip] = npp_vgpm(chl, par[ip], sst, tset - trise);
                    else
                        l1rec->flags[ip] |= PRODFAIL;

                    break;
                case CAT_npp_eppley:
                    if (sst > -2)
                        prod[ip] = npp_eppley(chl, par[ip], sst, tset - trise);
                    else
                        l1rec->flags[ip] |= PRODFAIL;
                    //if (ip % 100 == 0)
                    //    printf("RJH: %f %f %f %f %f\n",prod[ip],chl, par[ip], sst, tset-trise);
                    break;
                case CAT_npp_cbpm2:
                    if (firstCall) {
                        ib440 = bindex_get(440);
                        if (ib440 < 0) {
                            printf("opp_cbpm2: incompatible sensor wavelengths (no 440 for backscatter coefficient).\n");
                            exit(1);
                        }
                    }
                    mld = get_mld(input->mldfile, l1rec->lon[ip], l1rec->lat[ip], day); //193
                    zno3 = get_zno3(l1rec->lon[ip], l1rec->lat[ip], day); //125
                    bbp = l2rec->bb[l1rec->l1file->nbands * ip + ib440];
                    irr = par[ip]; //53
                    //                    bbp = 0.005;
                    //                    mld = 50;
                    //                    zno3=125;
                    //                    irr=53;
                    //                    tset=19;
                    //                    trise=5;
                    //                    kdtmp=0.03; //kd[ip]
                    //                    chl = 0.1;
                    //                    if (mld == BAD_FLT) mld = 50;
                    //                    if (zno3 == BAD_FLT) zno3 = 50;
                    if (bbp > 0 && mld != BAD_FLT && zno3 != BAD_FLT) {
                        prod[ip] = npp_cbpm2(chl, bbp, irr, kd[ip], mld, zno3, tset - trise);
                        //if (ip % 1 == 0)
                        // printf("RJH: %d %d %f %f %f %f %f %f %f %f %d ",ip,l2rec->iscan,chl, irr, sst, bbp, kd[ip], mld, zno3, tset-trise, *l2rec->day);
                        //if (prod[ip] == 0) prod[ip] = BAD_FLT;
                        // if (ip % 1 == 0)
                        //   printf("%f\n",prod[ip]);
                    } else {
                        prod[ip] = BAD_FLT;
                        //l2rec->flags[ip] |= PRODFAIL;
                    }
                    break;
                case CAT_npp_cafe:
                    if(firstCall){
                        ib443=bindex_get(443);
                        if (ib443 < 0) {
                            printf("npp_cafe: incompatible sensor wavelengths (no 443 for aph, adg, and bbp).\n");
                            exit(1);
                        }
                    }
                    mld = get_mld(input->mldfile, l1rec->lon[ip], l1rec->lat[ip], day);

                    if (!giop_ran(l1rec->iscan)){
                          run_giop(l2rec);
                    }
                    aph=giop_get_aph_pointer();
                    adg=giop_get_adg_pointer();
                    bbp_temp=giop_get_bbp_pointer();
                    bbp_s=giop_get_bbp_s_pointer();
                    ipb=ip*nbands+ib443;

                    prod[ip]=opp_cafe(par[ip],chl,mld,l1rec->lat[ip],day,aph[ipb],adg[ipb],bbp_temp[ipb],bbp_s[ip],sst);
                    break;
                default:
                    printf("Error: %s : Unknown product specifier: %d\n", __FILE__, prodnum);
                    exit(1);
                    break;
                }
            } else {
                l1rec->flags[ip] |= PRODFAIL;
                prod[ip] = BAD_FLT;
            }
            break;
        }


    }

    free(kd);
    free(par);

    firstCall = 0;
}

void get_par_clim(float *parin, float *lat, float *lon, int nlat, int nlon, float *latp, float *lonp, int32_t npix, float *par) {

    int32_t i, j, n, incx, incy, startx, starty;
    float dx, dy;

    dy = *(lat + 1) - *(lat);
    dx = *(lon + 1) - *(lon);

    if (dx > 0) {
        incx = -1;
        startx = nlon - 1;
    } else {
        incx = 1;
        startx = 0;
    }
    if (dy > 0) {
        incy = -1;
        starty = nlat - 1;
    } else {
        incy = 1;
        starty = 0;
    }

    for (n = 0; n < npix; n++) {
        for (i = startx; i >= 0 && i < nlon && (fabs(lonp[n] - lon[i]) > fabs(dx)); i += incx);
        for (j = starty; j >= 0 && j < nlat && (fabs(latp[n] - lat[j]) > fabs(dy)); j += incy);
        //            printf("RJH: lat[%d]=%f latp[%d]=%f diff=%f dy=%f\n",j,lat[j],n,latp[n],fabs(latp[n] - lat[j]),dy);
        if (j < nlat && i < nlon && i > 0 && j > 0)
            par[n] = parin[j * nlon + i];
        else
            par[n] = BAD_FLT;

        //        printf("RJH: n=%d par=%f latp=%f lonp=%f lat[%d]=%f lon[%d]=%f\n",n,par[n], latp[n],lonp[n],j,lat[j],i,lon[i]);
    }
}

/**
!C--------------------------------------------------------------------------*\

   !Description:     opp_befa - computes daily primary productivity using
                     the Behrenfeld-Falkowski (BeFa) algorithm.  The BeFa
                     algorithm estimates productivity using surface chl
                     (mg m-3), surface irradiance (Einsteins m-2 d-1),
                     sea surface temperature (C).
             Pb_opt is modelled as a polynomial function of SST.

   !Input Parameters:
      @param[in] chl            Chlorophyll_a surface concentration in milligrams
                                chlorophyl per cubic meter
      @param[in] irr            Photosynthetically available radiation in Einsteins per
                                day per square meter
      @param[in] sst            Sea surface temperature in degrees Centigrade
      @param[in] dayL           Length day in decimal hours.

   !Output Parameters:
      @param[out]       Primary productivity in milligrams Carbon per square meter
                        per hour

   !Revision History:

      First programmed up by Monica Chen at Rutgers
      (1996)

      Revised by K. Turpie at NASA
      (August 1997)

      Maintained by Don Shea at NASA

      Now maintained by Robert O'Malley at Oregon State University
      (April, 2005 - present)

      Modified for inclusion in l2gen for OBPG by R. Healy at NASA
      (January 2015)

   !References and Credits

      Behrenfeld,M.J; Falkowski,P.G.; 1997. Photosynthetic Rates Derived
      from Satellite-Based Chlorophyll Concentration.  Limnology and
      Oceanography, Volume 42, Number 1

!END------------------------------------------------------------------------*\
 **/


float npp_vgpm(float chl, float irr, float sst, float dayL) {

    float chl_tot,
            z_eu,
            pb_opt,
            irrFunc,
            npp;


    /* Calculate euphotic depth (z_eu) with Morel's Case I model.            */
    /* Calculate chl_tot from Satellite Surface Chlorophyll Data.            */

    if (chl < 1.0)
        chl_tot = 38.0 * pow(chl, 0.425);
    else
        chl_tot = 40.2 * pow(chl, 0.507);


    z_eu = 200.0 * pow(chl_tot, (-.293));

    if (z_eu <= 102.0)
        z_eu = 568.2 * pow(chl_tot, (-.746));


    /* Calculate the Pb_opt from satellite sea surface temperature (sst).    */

    if (sst < -10.0)
        pb_opt = 0.00;
    else if (sst < -1.0)
        pb_opt = 1.13;
    else if (sst > 28.5)
        pb_opt = 4.00;
    else {
        pb_opt = 1.2956 + 2.749e-1 * sst + 6.17e-2 * pow(sst, 2) - 2.05e-2 * pow(sst, 3)
                + 2.462e-3 * pow(sst, 4) - 1.348e-4 * pow(sst, 5) + 3.4132e-6 * pow(sst, 6)
                - 3.27e-8 * pow(sst, 7);
    }


    /* calculate the irradiance function */

    irrFunc = 0.66125 * irr / (irr + 4.1);


    /* Return the primary production calculation.                            */

    npp = pb_opt * chl * dayL * irrFunc * z_eu;

    if (npp < 0)
	npp = BAD_FLT;

    return npp;
}

/*
!C--------------------------------------------------------------------------*\

   !Description:     opp_eppley - computes daily primary productivity using
                     the Behrenfeld-Falkowski (BeFa) algorithm, but
             modifies the pb_opt function after Eppley (as
             implemented by Antoine and Morel).  The BeFa
                     algorithm estimates productivity using surface chl
                     (mg m-3), surface irradiance (Einsteins m-2 d-1),
                     sea surface temperature (C), and day length (hours).
             Pb_opt is modelled as an exponential function of SST.

   !Input Parameters:
      @param[in] chl            Chlorophyll_a surface concentration in milligrams
                                chlorophyl per cubic meter
      @param[in] irr            Photosynthetically available radiation in Einsteins per
                                day per square meter
      @param[in] sst            Sea surface temperature in degrees Centigrade
      @param[in] dayL           Length day in decimal hours.

   !Output Parameters:
      @param[out]       Primary productivity in milligrams Carbon per square meter
                        per hour

   !Revision History:

      First programmed up by Monica Chen at Rutgers
      (1996)

      Revised by K. Turpie at NASA
      (August 1997)

      Maintained by Don Shea at NASA

      Now maintained by Robert O'Malley at Oregon State University
      (April, 2005 - present)

      Modified for inclusion in l2gen for OBPG by R. Healy at NASA
      (January 2015)

   !References and Credits

      Behrenfeld,M.J; Falkowski,P.G.; 1997. Photosynthetic Rates Derived
      from Satellite-Based Chlorophyll Concentration.  Limnology and
      Oceanography, Volume 42, Number 1

      Eppley, R.W.; 1972.  Temperature and Phytoplankton Growth in the Sea.
      Fishery Bulletin, Volume 79, Number 4

      Antoine, D.; Morel, A.; 1996.  Oceanic Primary Production
      1.  Adatptation of a Spectral Light-Photosynthesis Model
      in view of Application to Satellite Chlorophyll Observations

!END------------------------------------------------------------------------*\
 */

float npp_eppley(float chl,
        float irr,
        float sst,
        float dayL) {

    float chl_tot,
            z_eu,
            pb_opt,
            irrFunc,
            npp;


    /* Calculate euphotic depth (z_eu) with Morel's Case I model.            */
    /* Calculate chl_tot from Satellite Surface Chlorophyll Data.            */

    if (chl < 1.0)
        chl_tot = 38.0 * pow(chl, 0.425);
    else
        chl_tot = 40.2 * pow(chl, 0.507);


    z_eu = 200.0 * pow(chl_tot, (-.293));

    if (z_eu <= 102.0)
        z_eu = 568.2 * pow(chl_tot, (-.746));


    /* Calculate the Pb_opt from satellite sea surface temperature (sst).    */

    pb_opt = 1.54 * pow(10, 0.0275 * sst - 0.07);


    /* calculate the irradiance function */

    irrFunc = 0.66125 * irr / (irr + 4.1);


    /* Return the primary production calculation.                            */

    npp = pb_opt * chl * dayL * irrFunc * z_eu;

    if (npp < 0)
	npp = BAD_FLT;

    return npp;
}

/**

   !Description:     opp_cbpm2 - computes daily primary productivity using a chl:Carbon ratio.
                     This is a spectrally resolved version of the cbpm, using nine separate
                     wavelengths.  It is also depth resolved, integrating the effects from
                     the surface down to a fixed depth of 200 m.

                     The cbpm2 algorithm estimates productivity using chl (m-1), bbp (m-1),
             surface irradiance (Einsteins m-2 d-1), k490 (m-1), mld (m), zno3 (m)
                     and day length (hours).

Net primary productivity is carbon * growth rate, where carbon is proportional to particulate
backscatter

    carbon = 13000 * (bbp - 0.00035)

and growth rate is a function of nutrient and temperature stress (f(nut,T) and photoacclimation
(f(Ig))

    growth rate (u) = umax * f(nut,T) * f(Ig)

where:

    umax = 2

    f(nut,T) = ((Chl/C)sat - y0) / ((Chl/C)max - y0)

    f(Ig) = 1 - exp (-5 * Ig)


and:

    (Chl/C)sat = ratio of satellite observed chl and carbon (carbon from bbp)

    (Chl/C)max = 0.022 + (0.045-0.022) * exp (-3 * Ig)

    Ig = median mixed layer light level
       = surface irradiance * exp (-k(lambda) * MLD/2)

The above items are analyzed for nine separate wavelengths, and is vertically resolved to a depth
of 200 m.

For more details, please see the paper by Westberry, et al (2008)


   !Input Parameters:
      @param[in] chl            chlorophyll concentration
      @param[in] bbp            backscatter
      @param[in] irr            Photosynthetically available radiation in Einsteins per
                     day per square meter
      @param[in] k490           absorbence at 490nm
      @param[in] mld            mixing layer depth in meters
      @param[in] zno3           depth of the nitrocline
      @param[in] daylength      length of the day in decimal hours.

   !Output Parameters:
      @param[out]        Primary productivity in milligrams Carbon per square meter
                     per day

   !Dependencies:
      function austinPetzold_1986 ( float lambda, float K490 )

         given a reference k490 vlaue, determine k(lambda) for a specified lambda

         ref:
            Austin, R. W., and T. J. Petzold (1986), Spectral dependence of the diffuse
            attenuation coefficient of light in ocean waters, Opt. Eng., 25, 473 – 479

   !Revision History:

   08-16-2010 first release version (Robert O'Malley)
      [original code written in matlab by T. Westberry]

   01-05-2011   O'Malley
      add uMax trap on mu[m]
      correct z_eu determination

   !References and Credits

      Westberry, T. Behrenfeld, M.J., Siegel, D.A., and Boss, E.; 2008.  Carbon-based
      primary productivity modeling with vertically resolved photoacclimation.  Global
      Biogeochemical Cycles, Vol. 22, GB2024, doi:10.1029/2007GB003078

 */

double npp_cbpm2(double chl,
        double bbp,
        double irr,
        double k490,
        double mld,
        double zno3,
        double daylength) {

    double austinPetzold_1986(double, double);

    double uMax; /* max growth rate */
    double chlCarbonMax; /* max chl:carbon ration */
    double nutTempFunc; /* f(nut,T) */
    double chlCarbonSat; /* satalite chl:carbon ratio */
    double carbon; /* bbp converted to carbon */
    double IgFunc; /* f(Ig) */
    double IgFuncz; /* f(Ig) below the mixed layer depth */
    double z_eu; /* euphotic depth at 1% light level */
    double npp; /* net primary production */

    /* --------------------- */
    /*   spectral variables  */
    /* --------------------- */

    double lambda[] = {400, 412, 443, 490, 510, 555, 625, 670, 700};
    double parFraction[] = {0.0029, 0.0032, 0.0035, 0.0037, 0.0037, 0.0036, 0.0032, 0.0030, 0.0024};
    double X[] = {.11748, .122858, .107212, .07242, .05943, .03996, .04000, .05150, .03000};
    double e[] = {.64358, .653270, .673358, .68955, .68567, .64204, .64700, .69500, .60000};
    double Kw[] = {.01042, .007932, .009480, .01660, .03385, .06053, .28400, .43946, .62438};
    double Kd[9];
    double Kbio;
    double Kdif[9];

    double Klambda[9];
    double Eo[9];
    double Ez_mld[9];
    double par_mld;
    double delChlC;

    double y0;

    /* --------------------------- */
    /*   depth resolved variables  */
    /* --------------------------- */

    double z[200]; /* depths */
    double chl_C[200]; /* chl:c ratio */
    double chlz[200]; /* chl */
    double mu[200]; /* growth */
    double Ezlambda[9][200]; /* fraction of light at nine wavelengths */
    double parz[200]; /* total light */
    double prcnt[200]; /* percent light */
    double Cz[200]; /* carbon */
    double ppz[200]; /* npp */

    int i;
    int m;
    int mzeu;
    double r;
    double prcnt0;
    double prcnt1;
    double z0;
    double z1;
    double numerator;
    double denominator;
    double fraction;
    double deltaZ;

    if (irr <= 0.0) {
        return 0.0;
    }

    /* --------------------- */
    /*   initialize values   */
    /* --------------------- */

    z_eu = -9999; //  1.05.2011
    y0 = 0.0003; /* min  chl:c  when  mu = 0 */
    for (i = 0; i < 200; i++) {
        z[i] = (float) (i + 1);
    }
    r = 0.1;

    uMax = 2.0; /* after Banse (1991) */
    npp = 0.0;
    mzeu = 0.0;

    for (i = 0; i < 9; i++) {
        Klambda[i] = austinPetzold_1986(lambda[i], k490);
        Eo[i] = irr * parFraction[i];
        Ez_mld[i] = Eo[i] * 0.975 * exp(-Klambda[i] * mld / 2.0);
    }

    /* ----------------------------- */
    /*   reintegrate to get par at   */
    /*   depth ...                   */
    /*   do trapezoidal integration  */
    /* ----------------------------- */

    par_mld = 0.0;
    for (i = 0; i < 8; i++) {
        par_mld += (lambda[i + 1] - lambda[i])*(Ez_mld[i + 1] + Ez_mld[i]) / 2;
    }

    par_mld /= daylength;

    IgFunc = 1 - exp(-5.0 * par_mld);

    if (bbp < 0.00035)
        bbp = 0.00036;
    carbon = 13000.0 * (bbp - 0.00035);

    chlCarbonSat = chl / carbon;

    if (chlCarbonSat < y0) {
        chlCarbonSat = y0;
    }

    chlCarbonMax = 0.022 + (0.045 - 0.022) * exp(-3.0 * par_mld);
    delChlC = chlCarbonMax - chlCarbonSat;

    nutTempFunc = (chlCarbonSat - y0) / (chlCarbonMax - y0);

    /* ''''''''''''''''''''''''' */
    /*   calculate Kd offset     */
    /*   carry through to depth  */
    /*   non-chl attenuation     */
    /* ------------------------- */

    for (i = 0; i < 9; i++) {
        Kbio = X[i] * pow(chl, e[i]);
        Kd[i] = Kw[i] + Kbio;
        Kdif[i] = Klambda[i] - Kd[i];
    }

    /* ''''''''''''''''''''''''''''''''''' */
    /*   integrate down the water column   */
    /*   in one-meter steps                */
    /* ----------------------------------- */

    for (m = 0; m < 200; m++) {

        /* ---------------------------------------------- */
        /*   if you are in the mixed layer, do this way   */
        /* ---------------------------------------------- */

        if (z[m] < mld) {
            chl_C[m] = chlCarbonSat;
            chlz[m] = chl_C[m] * carbon;
            mu[m] = uMax * nutTempFunc * IgFunc;

            if (mu[m] > uMax) { //  1.05.2011
                mu[m] = uMax; //  1.05.2011
            } //  1.05.2011

            for (i = 0; i < 9; i++) {
                Ezlambda[i][m] = Eo[i]*0.975 * exp(-Klambda[i] * z[m]);
            }

            parz[m] = 0.0;
            for (i = 0; i < 8; i++) {
                parz[m] += (lambda[i + 1] - lambda[i])*(Ezlambda[i + 1][m] + Ezlambda[i][m]) / 2;
            }

            Cz[m] = carbon;

        } else {

            /* '''''''''''''''''''''''''''''''''''''''''''''''''''''''''' */
            /*   if below mixed layer must treat properties differently   */
            /* ---------------------------------------------------------- */

            for (i = 0; i < 9; i++) {
                Kbio = X[i] * pow(chlz[m - 1], e[i]); /*  after Morel & Maritorena (2001)  */
                Kd[i] = Kw[i] + Kbio;
                Kd[i] += Kdif[i];
                Ezlambda[i][m] = Ezlambda[i][m - 1] * exp(-Kd[i]*1.0);
            }

            parz[m] = 0.0;
            for (i = 0; i < 8; i++) {
                parz[m] += (lambda[i + 1] - lambda[i])*(Ezlambda[i + 1][m] + Ezlambda[i][m]) / 2;
            }

            deltaZ = zno3 - z[m];
            if (deltaZ < 0) {
                deltaZ = 0;
            }

            chl_C[m] = (0.022 + (0.045 - 0.022) * exp(-3.0 * parz[m] / daylength));
            chl_C[m] -= delChlC * (1 - exp(-0.075 * deltaZ));

            IgFuncz = 1 - exp(-5.0 * parz[m] / daylength);
            mu[m] = uMax * nutTempFunc * IgFuncz;

            if (mu[m] > uMax) { //  1.05.2011
                mu[m] = uMax; //  1.05.2011
            } //  1.05.2011

            if (mu[m - 1] >= r) {
                Cz[m] = carbon;
            } else {
                Cz[m] = carbon * mu[m - 1] / r;
            }

            chlz[m] = chl_C[m] * Cz[m];

        }

        prcnt[m] = parz[m] / (irr * 0.975);

        /*  track this to get to the euphotic depth  */

        if (prcnt[m] >= 0.01) {
            mzeu = m;
        } else {

            /* ''''''''''''''''''''''''''' */
            /*   now find 1% light depth   */
            /*   in case the user wants    */
            /*   to use this information   */
            /* --------------------------- */

            if (z_eu == -9999) { // 01.05.11
                prcnt0 = prcnt[mzeu];
                prcnt1 = prcnt[mzeu + 1];
                z0 = z[mzeu];
                z1 = z[mzeu + 1];
                numerator = prcnt0 - 0.01;
                denominator = prcnt0 - prcnt1;
                fraction = numerator / denominator;
                z_eu = z0 + (z1 - z0) * fraction;
            }
        }

        ppz[m] = mu[m] * Cz[m];

    }

    /* ------------------------------- */
    /*   do trapezoidal integration    */
    /*   from m = 0 to m = 200         */
    /* ------------------------------- */

    //  note:  186 m is the euphotic depth for pure water

    if (mzeu < 186) {
        npp = 0;
        for (i = 0; i < 199; i++) {
            npp += (z[i + 1] - z[i])*(ppz[i + 1] + ppz[i]) / 2;
        }
    } else {
        npp = BAD_FLT;
    }

    if (npp < 0)
	npp = BAD_FLT;

    return npp;
}

/* =================================================================  */

double austinPetzold_1986(double lambda,
        double K490) {

    double wave[] = {350, 360, 370, 380, 390, 400,
        410, 420, 430, 440, 450, 460, 470, 480, 490, 500,
        510, 520, 530, 540, 550, 560, 570, 580, 590, 600,
        610, 620, 630, 640, 650, 660, 670, 680, 690, 700};

    double M[] = {2.1442, 2.0504, 1.9610, 1.8772, 1.8009, 1.7383,
        1.7591, 1.6974, 1.6108, 1.5169, 1.4158, 1.3077, 1.1982, 1.0955, 1.0000, 0.9118,
        0.8310, 0.7578, 0.6924, 0.6350, 0.5860, 0.5457, 0.5146, 0.4935, 0.4840, 0.4903,
        0.5090, 0.5380, 0.6231, 0.7001, 0.7300, 0.7301, 0.7008, 0.6245, 0.4901, 0.2891};

    double Kdw[] = {0.0510, 0.0405, 0.0331, 0.0278, 0.0242, 0.0217,
        0.0200, 0.0189, 0.0182, 0.0178, 0.0176, 0.0176, 0.0179, 0.0193, 0.0224, 0.0280,
        0.0369, 0.0498, 0.0526, 0.0577, 0.0640, 0.0723, 0.0842, 0.1065, 0.1578, 0.2409,
        0.2892, 0.3124, 0.3296, 0.3290, 0.3559, 0.4105, 0.4278, 0.4521, 0.5116, 0.6514};

    double l0;
    double l1;
    double k0;
    double k1;
    double m0;
    double m1;
    double kdiff;
    double mdiff;
    double num;
    double den;
    double frac;
    double Kdw_l;
    double M_l;
    double Kd;

    int ref;
    int i;

    // -- INTERPOLATE TO WAVELENGTH OF INTEREST --  //

    for (i = 1; i < 36; i++) {
        if (wave[i] >= lambda) {
            l1 = wave[i];
            k1 = Kdw[i];
            m1 = M[i];
            l0 = wave[i - 1];
            k0 = Kdw[i - 1];
            m0 = M[i - 1];
            break;
        }
    }

    if(i==36) {
        printf("-E- %s:%d - lambda = %f, out of range\n", __FILE__, __LINE__, lambda);
        exit(EXIT_FAILURE);
    }

    num = lambda - l0;
    den = l1 - l0;
    frac = num / den;

    kdiff = k1 - k0;
    Kdw_l = k0 + frac*kdiff;

    mdiff = m1 - m0;
    M_l = m0 + frac*mdiff;


    // -- GET REFERENCE WAVELENGTH (=490 FOR NOW) AND APPLY MODEL -- //

    ref = 14;

    Kd = (M_l / M[ref]) * (K490 - Kdw[ref]) + Kdw_l;

    return Kd;

}


float opp_cafe( float PAR,
         float chl,
         float mld,
         float lat,
         int    yd,
         float aph_443,
         float adg_443,
         float bbp_443,
         float bbp_s,
         float sst ) {

/* ----------------------------------------------------------------------------------------------*/

/* Description: Computes daily net primary production following the Carbon, Absorption, Fluorescence
                Euphotic-resolving (CAFE) model. For a full description of the model, please see:

                Silsbe, G.M., M.J. Behrenfeld, K.H. Halsey, A.J. Milligan, and T.K. Westberry. 2016.
                The CAFE model. A net production model for global ocean phytoplankton. Global
                Biogeochemical Cycles. doi: 10.1002/2016GB005521.

   Model inputs:
    PAR     - Daily photosynthetic active radiation [mol photons m^-2 day^-1]
    chl     - Chlorophyll concentration, NASA OCI algorithm [mg m^-3]
    mld     - Mixed layer depth [m]
    lat     - Latitude [Degrees, north positive]
    yd      - Day of year
    aph_443 - Absorption due to phytoplankton at 443 nm, NASA GIOP algorithm [m-1]
    adg_443 - Absorption due to gelbstoff and detratial material at 443 nm, NASA GIOP model [m-1]
    bbp_443 - Particulate backscatter at 443 nm, NASA GIOP model [m-1]
    bbp_s   - Backscattering spectral parameter for GIOP model [dimensionless]
    sst     - Sea surface temperature [Degrees Celcuis]

  Modeling code below is divided into 10 steps:

    1) Declare all local variables
    2) Derive inherent optical properties (IOPs) at 10 nm increments from 400 to 700 nm
    3) Calculate the amount of energy in the eupotic zone that is absorbed by phytoplankton
    4) Derive the spectral attenunation coefficient and the attenuation coefficient of PAR
    5) Derive conversion factor (Eu) that converts downwelling planar irradiance to scalar
       irradiance
    6) Calculate the photoacclimation parameter (Ek) through depth
    7) Scale Ek to a spectrally explicit parameter (KPUR)
    8) Model the maximum quantum efficiency of net carbon fixation (phi_max)
    9) If the mixed layer depth (MLD) is shallower than the euphotic depth, then apply scaling
       factor to phytoplankton absorption beneath the MLD and adjust the energy absorbed by
       phytoplankton.
    10) Derive net phytoplankton production (NPP).
    11) Dependent function that calculates the wavelength, salinity and temperature dependence
        of the backscattering of pure water.
*/


/* ----------------------------------------------------------------------------------------------*/
/* Step 1: Declare all local variables */
/* ----------------------------------------------------------------------------------------------*/

  int w; /* Wavelength step 10 nm between 400 and 700 nm */
  int t; /* Time step, day is divided into 101 evenly spaced increments */
  int tsym; /* morning index for equivalent afternoon time */
  int z; /* Depth step, euphotic zone is divided into 101 evenly spaced increments */

  /* Inherent Optical Properties*/
  double aphi[31];    /* Absorption[wavelength] due to phytoplankton */
  double a[31];       /* Total absorption[wavelength] */
  double bb[31];      /* Total backscattering[wavelength] */
  double bbw[31];     /* Backscattering of pure water[wavelength] */

  /* Wavelength */
  double wv[] = { 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560,
          570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700 };

  /* Absorption due to pure water[wavelength] */
  double aw[] = { 0.00663, 0.00473, 0.00454, 0.00495, 0.00635, 0.00922, 0.00979, 0.0106, 0.0127,
          0.015, 0.0204, 0.0325, 0.0409, 0.0434, 0.0474, 0.0565, 0.0619, 0.0695,
          0.0896, 0.1351, 0.2224, 0.2644, 0.2755, 0.2916, 0.3108, 0.34, 0.41, 0.439, 0.465,
          0.516, 0.624 };

  /*Spectral shape of of aphi */
  double A_Bricaud[] = { 0.0241, 0.0287, 0.0328, 0.0359, 0.0378, 0.0350, 0.0328, 0.0309, 0.0281,
             0.0254, 0.0210, 0.0162, 0.0126, 0.0103, 0.0085, 0.0070, 0.0057, 0.0050,
             0.0051, 0.0054, 0.0052, 0.0055, 0.0061, 0.0066, 0.0071, 0.0078, 0.0108,
             0.0174, 0.0161, 0.0069, 0.0025 };

  double E_Bricaud[] = { 0.6877, 0.6834, 0.6664, 0.6478, 0.6266, 0.5993, 0.5961, 0.5970, 0.5890,
             0.6074, 0.6529, 0.7212, 0.7939, 0.8500, 0.9036, 0.9312, 0.9345, 0.9298,
             0.8933, 0.8589, 0.8410, 0.8548, 0.8704, 0.8638, 0.8524, 0.8155, 0.8233,
             0.8138, 0.8284, 0.9255, 1.0286 };

  /* Spectral shape of PAR */
  double PAR_spectrum[] = { 0.00227, 0.00218, 0.00239, 0.00189, 0.00297, 0.00348, 0.00345,
                0.00344, 0.00373, 0.00377, 0.00362, 0.00364, 0.00360, 0.00367,
                0.00354, 0.00368, 0.00354, 0.00357, 0.00363, 0.00332, 0.00358,
                0.00357, 0.00359, 0.00340, 0.00350, 0.00332, 0.00342, 0.00347,
                0.00342, 0.00290, 0.00314 };

  /* Parameters related to euphotic zone*/
  double decl;   /* Solar declination [rads] */
  double m0;     /* Coefficient to calculate kd (Lee et al. 2005) [m-1] */
  double DL;     /* Daylength [day] */
  double solzen; /* Solar zenith angle [rads] */
  double kd[31]; /* Downwelling attenuation coefficient of irradiance [m-1] */
  double zeu;    /* Euphotic depth [m] */
  double kdpar;  /* Downwelling attenuation coefficient of PAR [m-1] */

  /* Parameters related to light and absorbed photons through depth, time, and wavelength */
  double tseq[101];               /* Fractional time of day */
  double zseq[101];               /* Euphotic depth divided into 101 increments */
  double delz;                    /* Depth of zseq (i.e. zeu/101) */
  double absorbed_photons = 0.0;  /* Absorbed photons calculated in Step 3 */
  double PAR_noon[31];            /* PAR at solar noon[wavelength] [mol photons m-2 day-1 wv-] */
  double AP_tzw[101][101][31];    /* Absorbed photons[depth, time, wavelength] */
  double AP_tz[101][101];         /* Absorbed photons[depth, time] */
  double AP_tz2[101][101];         /* Absorbed photons[depth, time] */
  double AP_z[101];               /* Absorbed photons[depth] */
  double AP_z2[101];               /* Absorbed photons[depth] */
  double AP = {0.0};              /* Absorbed photons [mol m-2 day-1] */
  double AP2 = {0.0};              /* Absorbed photons [mol m-2 day-1] */
  double E_tzw[101][101][31];     /* Irradiance [depth, time, wavelength] */
  double E_tz[101][101];          /* Irradiance [depth, time] */
  double Eu;                      /* conversion factor (Eu) that converts downwelling planar
                                   irradiance to scalar irradiance */

  /* Parameters related to derivation of Ek*/
  double Ek[101];  /* Photoacclimation parameter[depth] (Behrenfeld et al. 2016) [mol m-2 day-1] */
  double IML;      /* Median irradiance in the mixed layer [umol m-2 day-1] */
  double Eg;       /* Growth irradiance at a specific depth */
  double Eg_mld;   /* Growth irradiance at a specific depth */

  /* Parameters used to convert Ek into the spectrally explicit equivalent KPUR */
  double KPUR[101];        /* Spectrally explicit Ek */
  double numerator;        /* Numerator to calculate spectral correction factor */
  double denominator;      /* Denominator to calculate spectral correction factor */
  double mean_aphi = 0.0;  /* Spectrally averaged absorption due to phytoplankton */

  /* Parameters used to calculate the maximum quantum yield of net carbon fixation */
  double phimax[101];  /* maximum quantum yield of net carbon fixation [mol C (mol photons)-1] */
  double phirange[2] = {0.018, 0.030};            /* range of phimax values */
  double Ekrange[2] = {150*0.086400, 10*0.0864};  /* relates phimax to Ek */
  double slope;                                   /* relates phimax to Ek */

  double aphi_fact[101]; /* Parameter to increase aphi beneath the MLD */

  /* Net primary production [mg C m-2 day-1] */
  double NPP_tz[101][101] = {{0.0}};
  double NPP_z[101] = {0.0};
  double NPP = {0.0};


  //FILE *fp;
  //  fp = fopen("test.txt", "a");

  /* ----------------------------------------------------------------------------------------------*/
  /* Step 2: Derive IOPs at 10 nm increments from 400 to 700 nm                                    */
  /*   Comment: Currently assume salinity is 32.5 ppmil                                            */
  /* ----------------------------------------------------------------------------------------------*/

  for ( w=0; w<31; w++ ){
    bbw[w]  = betasw_ZHH2009(wv[w], 32.5, sst) / 2;
    aphi[w] = aph_443 * A_Bricaud[w] * pow(chl, E_Bricaud[w]) / (0.03711 * pow(chl, 0.61479));
    a[w]    = aw[w] + aphi[w] + adg_443 * exp(-0.018 * (wv[w] - 443));
    bb[w]   = bbw[w] + bbp_443 * pow(443 / wv[w], bbp_s);
  }

  /* ----------------------------------------------------------------------------------------------*/
  /* Step 3: Calculate Absorbed Energy */
  /* ----------------------------------------------------------------------------------------------*/

  for (w=0; w<30; w++){
    absorbed_photons += 0.5 * 10 * (PAR_spectrum[w+1] * aphi[w+1]/a[w+1] +
                    PAR_spectrum[w] * aphi[w]/a[w]);
  }

  absorbed_photons *= PAR * 0.95;

  /* ----------------------------------------------------------------------------------------------*/
  /* Step 4: Derive Kd following Lee et al 2005, Eq. 11 and kdpar */
  /*         from Morel et al 2007, Eq. 9'                        */
  /*---------- ------------------------------------------------------------------------------------*/

  decl = 23.5 * cos (2 * M_PI * (yd - 172) / 365) * M_PI / 180;
  DL   = -1 * tan(lat*M_PI/180) * tan(decl);

  if (DL > 1)  {DL = 1;}   /* Check for daylength less than 0 hours */
  if (DL < -1) {DL = -1;}  /* Check for daylength greater than 24 hours */

  DL = acos(DL) / M_PI; /* Daylength in days */

  solzen = 90 - asin (sin(lat*M_PI/180) * sin(decl) - cos(lat*M_PI/180) * cos(decl) *
              cos(M_PI)) * 180 / M_PI;
  m0     = sqrt((1 + 0.005 * solzen) * (1 + 0.005 * solzen)); /* absolute value */

  for (w=0; w<31; w++){
    kd[w] =  m0 * a[w] + 4.18 * (1 - 0.52 * exp(-10.8 * a[w])) * bb[w];
  }

  kdpar = 0.0665 + (0.874 * kd[9]) - (0.00121 / kd[9]);
  zeu   = -1 * log (0.1 /(PAR * 0.95)) / kdpar;

  /* debug variable */
  //   oppDebugDepthEu = z_eu;


  /* ----------------------------------------------------------------------------------------------*/
  /* Step 5: Derive conversion factor (Eu) that converts downwelling planar irradiance to scalar   */
  /*   irradiance */
  /* ----------------------------------------------------------------------------------------------*/

  for (t=0; t<101; t++){
    ////  tseq[t] = (double)t / 100;
    zseq[t] = (double)t / 100 * ceil(zeu);
  }
  for (t=0; t<51; t++){
    tseq[t] = (double)t / 50;
  }
  delz = zseq[1] - zseq[0];

  for (w=0; w<31; w++){
    PAR_noon[w] = M_PI / 2 * PAR * 0.95  * PAR_spectrum[w];
  }

  ////for (t=0; t<101; t++){

  /* ==================================== */
  /*   light is symmetrical about noon    */
  /*   take advantage of that             */
  /*   to get rid of 1/2 the calculations */
  /* ==================================== */

  //  for (t=0; t<51; t++){
  for (t=0; t<26; t++){
    for (z=0; z<101; z++){
      for (w=0; w<31; w++){
    E_tzw[t][z][w]  = PAR_noon[w] * sin(M_PI*tseq[t]) * exp(-1*kd[w]*zseq[z]);
    AP_tzw[t][z][w] = E_tzw[t][z][w] * aphi[w];
      }
    }
  }
  for (t=26; t<51; t++){
    tsym = 25 - (t - 25);
    for (z=0; z<101; z++){
      for (w=0; w<31; w++){
    E_tzw[t][z][w]  = E_tzw[tsym][z][w];
    AP_tzw[t][z][w] = AP_tzw[tsym][z][w];
      }
    }
  }

  /* Integrate E and AP through wavelength*/
  for (z=0; z<101; z++){
    ////  for (t=0; t<101; t++){

  /* ==================================== */
  /*   light is symmetrical about noon    */
  /*   take advantage of that             */
  /*   to get rid of 1/2 the calculations */
  /* ==================================== */

    //    for (t=0; t<51; t++){
    for (t=0; t<26; t++){
      E_tz[t][z] = 0.0;
      AP_tz[t][z] = 0.0;
      for (w=0; w<30; w++){
    E_tz[t][z]  += 10 * (E_tzw[t][z][w+1] + E_tzw[t][z][w])/2;
    AP_tz[t][z] += 10 * (AP_tzw[t][z][w+1] + AP_tzw[t][z][w])/2;
      }
    }
    for (t=26; t<51; t++){
      tsym = 25 - (t - 25);
      for (w=0; w<30; w++){
    E_tz[t][z]  = E_tz[tsym][z];
    AP_tz[t][z] = AP_tz[tsym][z];
      }
    }
  }

  /* Integrate AP through time*/
  for (z=0; z<101; z++){
    AP_z[z] = 0.0;
    ////  for (t=0; t<100; t++){

  /* ==================================== */
  /*   light is symmetrical about noon    */
  /*   take advantage of that             */
  /*   to get rid of 1/2 the calculations */
  /* ==================================== */

    //    for (t=0; t<50; t++){
    for (t=0; t<25; t++){
      ////    AP_z[z] += 0.01 * (AP_tz[t+1][z] + AP_tz[t][z])/2;
      AP_z[z] += 0.02 * (AP_tz[t+1][z] + AP_tz[t][z])/2;
    }
    AP_z[z] = 2*AP_z[z];
 }

  /* Integrate AP through depth*/
  for (z=0; z<100; z++){
    AP += delz * (AP_z[z] + AP_z[z+1])/2;
  }

  /* Derive Upwelling Irradiance, absorbed energy is from Section 3 */
  Eu = absorbed_photons / AP;

  /* Multiply E_tz by Eu */

  ////for (t=0; t<101; t++){

  /* ==================================== */
  /*   light is symmetrical about noon    */
  /*   take advantage of that             */
  /*   to get rid of 1/2 the calculations */
  /* ==================================== */

  //  for (t=0; t<51; t++){
  for (t=0; t<26; t++){
    for (z=0; z<101; z++){
      E_tz[t][z] *= Eu;
      AP_tz[t][z] *= Eu;
    }
  }
  for (t=26; t<51; t++){
    tsym = 25 - (t - 25);
    for (z=0; z<101; z++){
      E_tz[t][z] =  E_tz[tsym][z];
      AP_tz[t][z] =  AP_tz[tsym][z];
    }
  }

  /* ----------------------------------------------------------------------------------------------*/
  /* Step 6: CALCULATE EK through depth (modified from Behrenfeld et al. 2016) */
  /* ----------------------------------------------------------------------------------------------*/

  IML   = (PAR * 0.95 / (DL* 24)) * exp(-0.5 * kdpar * mld);

  for (z=0; z<101; z++){
    Ek[z] = 19 * exp(0.038 * pow(PAR * 0.95 / (DL * 24), 0.45) / kdpar);
  }

  /* Modify Ek beneath the MLD */

  if (mld < zeu){
    Eg_mld = (PAR / DL)  * exp(-1*kdpar*mld);
    for (z=0; z<101; z++){
      Ek[z] *= (1 + exp(-0.15 * (0.95  * PAR / (DL * 24)))) / (1 + exp(-3 * IML));
      if (zseq[z] > mld){
    Eg = (PAR / DL)  * exp(-1*kdpar*zseq[z]);
    Ek[z] = 10 + (Ek[z] - 10) / (Eg_mld - 0.1) * (Eg -0.1);
      }
      if (Ek[z] < 10) {Ek[z] = 10;}
    }
  }

  for (z=0; z<101; z++){
    Ek[z] *= 0.0864;    /* 0.0864 #Convert to mol photons/m2/day */
  }

  /* debug variable */
  // oppDebugEk0 = Ek[0];


  /* ----------------------------------------------------------------------------------------------*/
  /* Step 7: Make Ek spectrally explicit (KPUR) using a spectral correction factor */
  /* ----------------------------------------------------------------------------------------------*/

  for (w=0; w<31; w++){
    mean_aphi += aphi[w];
  }
  mean_aphi /= 31;

  for (z=0; z<101; z++){
    numerator   = 0.0;
    denominator = 0.0;
    for (w=0; w<30; w++){
      ////    numerator +=   10 * 0.5 * (AP_tzw[49][z][w] + AP_tzw[49][z][w+1]);
      ////    denominator += 10 * 0.5 * (E_tzw[49][z][w] + E_tzw[49][z][w+1]);

      numerator +=   10 * 0.5 * (AP_tzw[24][z][w] + AP_tzw[24][z][w+1]);
      denominator += 10 * 0.5 * (E_tzw[24][z][w] + E_tzw[24][z][w+1]);
    }
    KPUR[z] = Ek[z] / (numerator / (denominator * mean_aphi) / 1.3);
  }

  /* ----------------------------------------------------------------------------------------------*/
  /* Step 8: Calculate phimax as a function of Ek */
  /* ----------------------------------------------------------------------------------------------*/

  slope = (phirange[1] - phirange[0]) / (Ekrange[1] - Ekrange[0]);

  for (z=0; z<101; z++){
    phimax[z] = phirange[1] + (Ek[z] - Ekrange[1]) * slope;
    if (phimax[z] > phirange[1]) {phimax[z] = phirange[1];}
    if (phimax[z] < phirange[0]) {phimax[z] = phirange[0];}
  }

  /* ----------------------------------------------------------------------------------------------*/
  /* Step 9: Scale aphi beneath the MLD to Ek */
  /* ----------------------------------------------------------------------------------------------*/

  if ( mld < zeu ) {
    for (z=0; z<101; z++){
      if (zseq[z] > mld){
    aphi_fact[z] = 1 + Ek[0]/Ek[z] * 0.15;
      } else {
    aphi_fact[z] = 1;
      }
    }

    //  for (z=21; z<23; z++){
    //    fprintf(fp, "z =  %d\n",z);
    //    fprintf(fp, " %lf",aphi_fact[z]);
    //    fprintf(fp,"\n");
    //  }

    /* Modify Irradiance and absorbed energy to account for any change to aphi */
    ////for (t=0; t<101; t++){

  /* ==================================== */
  /*   light is symmetrical about noon    */
  /*   take advantage of that             */
  /*   to get rid of 1/2 the calculations */
  /* ==================================== */

    //    for (t=0; t<51; t++){
    for (t=0; t<26; t++){
      for (w=0; w<31; w++){
    AP_tzw[t][0][w] = E_tzw[t][0][w] * aphi[w];
      }
      for (z=1; z<101; z++){
    for (w=0; w<31; w++){
      a[w]            = aw[w] + aphi[w] * aphi_fact[z] + adg_443 * exp(-0.018 * (wv[w] - 443));
      kd[w]           = m0 * a[w] + 4.18 * (1 - 0.52 * exp(-10.8 * a[w])) * bb[w];
      E_tzw[t][z][w]  = E_tzw[t][z-1][w] * exp(-1*kd[w]*delz);
      AP_tzw[t][z][w] = E_tzw[t][z][w] * aphi[w] * aphi_fact[z];
    }
      }
    }
    for (t=26; t<51; t++){
      tsym = 25 - (t - 25);
      for (w=0; w<31; w++){
    AP_tzw[t][0][w] = AP_tzw[tsym][0][w];
      }
      for (z=1; z<101; z++){
    for (w=0; w<31; w++){
      E_tzw[t][z][w]  = E_tzw[tsym][z][w];
      AP_tzw[t][z][w] = AP_tzw[tsym][z][w];
    }
      }
    }

    /* Integrate AP through wavelength and multiply by Eu*/
    for (z=0; z<101; z++){
      ////  for (t=0; t<101; t++){

  /* ==================================== */
  /*   light is symmetrical about noon    */
  /*   take advantage of that             */
  /*   to get rid of 1/2 the calculations */
  /* ==================================== */

      //      for (t=0; t<51; t++){
      for (t=0; t<26; t++){
    AP_tz2[t][z] = 0.0;
    for (w=0; w<30; w++){
      AP_tz2[t][z] += 10 * (AP_tzw[t][z][w+1] + AP_tzw[t][z][w])/2;
    }
    AP_tz2[t][z] *= Eu;
      }
      for (t=26; t<51; t++){
    tsym = 25 - (t - 25);
    AP_tz2[t][z] = AP_tz2[tsym][z];
      }
    }
    //    fprintf(fp, "\nEu =  %lf\n",Eu);
  } else {

    for (z=0; z<101; z++){
      ////  for (t=0; t<101; t++){
      for (t=0; t<51; t++){
    AP_tz2[t][z] = AP_tz[t][z];
      }
    }

  }


  /* ------------------------------------------------------------------------------------*/
  /* Step 9: CALCULATE NPP */
  /* ------------------------------------------------------------------------------------*/

  ////for (t=0; t<101; t++){

  /* ==================================== */
  /*   light is symmetrical about noon    */
  /*   take advantage of that             */
  /*   to get rid of 1/2 the calculations */
  /* ==================================== */

  //  for (t=0; t<51; t++){
  for (t=0; t<26; t++){
    for (z=0; z<101; z++){
      NPP_tz[t][z] = phimax[z] * AP_tz2[t][z] * tanh(KPUR[z] / E_tz[t][z]) * 12000;
    }
  }
  for (t=26; t<51; t++){
    tsym = 25 - (t - 25);
    for (z=0; z<101; z++){
      NPP_tz[t][z] = NPP_tz[tsym][z];
    }
  }


  /* Integrate NPP through time */
  for (z=0; z<101; z++){
    NPP_z[z] = 0.0;
    AP_z2[z] = 0.0;
    ////  for (t=0; t<100; t++){

  /* ==================================== */
  /*   light is symmetrical about noon    */
  /*   take advantage of that             */
  /*   to get rid of 1/2 the calculations */
  /* ==================================== */

    //    for (t=0; t<50; t++){
    for (t=0; t<25; t++){
      ////    NPP_z[z] += 0.01 * (NPP_tz[t+1][z] + NPP_tz[t][z])/2;
      NPP_z[z] += 0.02 * (NPP_tz[t+1][z] + NPP_tz[t][z])/2;
      AP_z2[z] += 0.02 * (AP_tz2[t+1][z] + AP_tz2[t][z])/2;
    }
    NPP_z[z] = 2*NPP_z[z];
    AP_z2[z] = 2*AP_z2[z];
  }

  /* Integrate NPP through depth */
  NPP = 0.0;
  AP2 = 0.0;
  for (z=0; z<100; z++){
    NPP += delz * (NPP_z[z] + NPP_z[z+1])/2;
    AP2 += delz * (AP_z2[z] + AP_z2[z+1])/2;
  }

  //  FILE *fp;
  //  fp = fopen("test.txt", "a");

  //  fprintf(fp, "\n %.14lf %.12lf %.13lf %.14lf\n", AP2, NPP, zeu, Ek[0]);
  //  fprintf(fp, " %lf %lf %lf\n",delz,Eu,absorbed_photons);

  //  for (z=0; z<100; z++){
  //    fprintf(fp, " %lf",AP_z2[z]);
  //  }

  //  fprintf(fp,"\n");


  //  fprintf(fp, " %.11lf %.11lf %.11lf\n",kd[0], kd[9], kdpar);
  //  fprintf(fp, " %.14lf %.14lf %.14lf\n",m0, a[9], bb[9]);
  //  fprintf(fp, " %.14lf %lf\n",bbw[9], wv[9]);
  //  fprintf(fp, " %.14lf %lf\n",aw[9], aphi[9]);

  //  for (w=0; w<30; w++){
  //    fprintf(fp, " %lf",aw[w] + aphi[w] + adg_443*exp(-0.018*(wv[w]-443)));
  //  }
  //  fprintf(fp,"\n");

  //  fprintf(fp, "adg_443:  %lf\n",adg_443);


  //  fprintf(fp,"      NPP     Ek     phimax    aphi_fact  KPUR    Ap_tz  zseq    index\n");
  //  for (z=0; z<101; z++){
  //    fprintf(fp, " %lf %lf %lf %lf %lf %lf %lf %lf %d\n",NPP_z[z],Ek[z],phimax[z],aphi_fact[z],KPUR[z],AP_tz[26][z],AP_tz2[26][z],zseq[z],z+1);
  //  }
  //  fprintf(fp,"\n");

  //  fclose(fp);

  return NPP;

}

double betasw_ZHH2009 (double lambda, double S, double Tc){

  /* function [theta,betasw,bsw,beta90sw]= betasw_ZHH2009(lambda,S,Tc,delta)
  Scattering by pure seawater: Effect of salinity Xiaodong Zhang, Lianbo Hu,
  and Ming-Xia He, Optics Express, 2009
  lambda (nm): wavelength
  Tc: temperature in degree Celsius, must be a scalar
  S: salinity, must be scalar
  delta: depolarization ratio, if not provided, default = 0.039 will be used.
  betasw: volume scattering at angles defined by theta. Its size is [x y],
  where x is the number of angles (x = length(theta)) and y is the number
  of wavelengths in lambda (y = length(lambda))
  beta90sw: volume scattering at 90 degree. Its size is [1 y]
  bw: total scattering coefficient. Its size is [1 y]
  for backscattering coefficients, divide total scattering by 2
  */

  /* values of the constants */
  double Na  = 6.0221417930e23 ;   /*  Avogadro's constant */
  double Kbz = 1.3806503e-23 ;     /*  Boltzmann constant */
  double Tk  = Tc + 273.15 ;       /*  Absolute tempearture */
  double M0  = 18e-3;              /*  Molecular weigth of water in kg/mol */
  double delta= 0.039;
  //double rad[181];
  double nsw;
  //double dnds;
  double IsoComp;
  double density;
  double dlnawds;
  double DFRI;
  double beta_df;
  double flu_con;
  double beta_cf;
  double beta90sw;
  double bsw;
  //double betasw;
  double n_air;
  double dnswds;

  //for (i=0; i<181; i++){
  //  rad[i] = i * M_PI/180;  /* angle in radian as a colum variable */
  //}

  /* nsw: absolute refractive index of seawater */
  /* dnds: partial derivative of seawater refractive index w.r.t. salinity */


  /* refractive index of air is from Ciddor (1996,Applied Optics) */
  n_air = 1.0 + (5792105.0 / (238.0185 - 1 / pow(lambda / 1e3, 2)) + 167917.0 /
                (57.362 - 1 / pow(lambda/1e3,2))) /1e8;

  /* refractive index of seawater is from Quan and Fry (1994, Applied Optics) */
  double n0 = 1.31405; double n1 = 1.779e-4 ; double n2 = -1.05e-6 ;
  double n3 = 1.6e-8 ; double n4 = -2.02e-6 ; double n5 = 15.868;
  double n6 = 0.01155; double n7 = -0.00423;  double n8 = -4382 ;
  double n9 = 1.1455e6;

  /* pure seawater */
  nsw    = n0 + (n1 + n2 * Tc + n3 * pow(Tc, 2)) * S + n4 * pow(Tc, 2) +
  (n5+n6*S+n7*Tc)/ lambda + n8 / pow(lambda, 2) + n9 / pow(lambda, 3);
  nsw   *= n_air;
  dnswds = (n1 + n2 * Tc + n3 * pow(Tc, 2) + n6 / lambda) * n_air;

  /* isothermal compressibility is from Lepple & Millero (1971,Deep Sea-Research), pages 10-11
  The error ~ ?0.004e-6 bar-1  */
  double kw;
  //double Btw_cal;
  //double Btw;
  double g0;
  double g1;
  double Ks;

  /* pure water secant bulk Millero (1980, Deep-sea Research) */
  kw = 19652.21 + 148.4206 * Tc - 2.327105 * pow(Tc, 2) + 1.360477e-2 * pow(Tc, 3) -
       5.155288e-5 * pow(Tc, 4);
  //Btw_cal = 1/kw;

  /* isothermal compressibility from Kell sound measurement in pure water */
  //Btw = (50.88630 + 0.717582 * Tc + 0.7819867e-3 * pow(Tc, 2) + 31.62214e-6 * pow(Tc, 3) -
  //       0.1323594e-6 * pow(Tc, 4) +
  //0.634575e-9 * pow(Tc, 5)) / (1 + 21.65928e-3 * Tc)*1e-6;

  /* seawater secant bulk */
  g0 = 54.6746 - 0.603459 * Tc + 1.09987e-2 * pow(Tc, 2) - 6.167e-5 * pow(Tc, 3);
  g1 = 7.944e-2 + 1.6483e-2 * Tc - 5.3009e-4 * pow(Tc, 2);

  Ks = kw + g0 * S + g1 * pow(S, 1.5);

  /* calculate seawater isothermal compressibility from the secant bulk */
  IsoComp = ( 1 / Ks * 1e-5); /* unit is pa */


  /* density of water and seawater,unit is Kg/m3, from UNESCO,38,1981

  density = density_sw(Tc, S);*/
  double a0 = 8.24493e-1; double a1 = -4.0899e-3; double a2 = 7.6438e-5;
  double a3 = -8.2467e-7; double a4 = 5.3875e-9; double a5 = -5.72466e-3;
  double a6 = 1.02270e-4; double a7 = -1.6546e-6; double a8 = 4.8314e-4;
  double b0 = 999.842594; double b1 = 6.793952e-2; double b2 = -9.09529e-3;
  double b3 = 1.001685e-4; double b4 = -1.120083e-6; double b5 = 6.536332e-9;

  /* density for pure water */
  density = b0 + b1 * Tc + b2 * pow(Tc, 2) + b3 * pow(Tc, 3) + b4 * pow(Tc, 4) + b5 * pow(Tc, 5);
  /* density for pure seawater */
  density += (( a0 + a1 * Tc + a2 * pow(Tc, 2) + a3 * pow(Tc, 3) + a4 * pow(Tc, 4)) * S +
  ( a5 + a6 * Tc + a7 * pow(Tc, 2)) * pow(S, 1.5) + a8 * pow(S, 2));

  /* water activity data of seawater is from Millero and Leung (1976,American
  Journal of Science,276,1035-1077). Table 19 was reproduced using
  Eq.(14,22,23,88,107) then were fitted to polynominal equation.
  dlnawds is partial derivative of natural logarithm of water activity
  w.r.t.salinity */

//  dlnawds = (-5.58651e-4 + 2.40452e-7 * Tc -3.12165e-9 * pow(Tc, 2) + 2.40808e-11 * pow(Tc, 3)) +
//            1.5*(1.79613e-5 -9.9422e-8 * Tc +2.08919e-9 * pow(Tc, 2) - 1.39872e-11 * pow(Tc, 3)) *
//            pow(S, 1.5) + 2 * (-2.31065e-6 - 1.37674e-9 * Tc -1.93316e-11 * pow(Tc, 2)) * S;

  dlnawds = (-5.58651e-4 + 2.40452e-7 * Tc -3.12165e-9 * pow(Tc, 2) + 2.40808e-11 * pow(Tc, 3)) +
            1.5*(1.79613e-5 -9.9422e-8 * Tc +2.08919e-9 * pow(Tc, 2) - 1.39872e-11 * pow(Tc, 3)) *
            pow(S, 0.5) + 2 * (-2.31065e-6 - 1.37674e-9 * Tc -1.93316e-11 * pow(Tc, 2)) * S;

  /* density derivative of refractive index from PMH model */
  double n_wat2;
  double base;
  //  double val;

  n_wat2 = pow(nsw, 2);
//  base   = (nsw/3 - 0.333333333333 / nsw);
  base   = (nsw/3 - (1.0/3.0) / nsw);
//  DFRI   = (n_wat2 - 1) * (1 + 0.666666666667 * (n_wat2 + 2) * pow(base, 2)) ;
  DFRI   = (n_wat2 - 1) * (1 + (2.0/3.0) * (n_wat2 + 2) * pow(base, 2)) ;
//  val = 1/3;
//  base   = (nsw/3) - (val/nsw);
//  DFRI   = (n_wat2 - 1) * (1 + (2/3) * (n_wat2 + 2) * pow(base, 2)) ;

  /* volume scattering at 90 degree due to the density fluctuation */
  beta_df = M_PI * M_PI / 2 * pow(lambda*1e-9, -4) * Kbz * Tk * IsoComp * pow(DFRI,2) *
            (6 + 6*delta) / (6-7*delta);

  /* volume scattering at 90 degree due to the concentration fluctuation */
  flu_con = S * M0 * pow(dnswds,2) / density / (-1 * dlnawds) / Na;
  beta_cf = 2 * M_PI * M_PI * pow(lambda*1e-9, -4) * pow(nsw, 2) * flu_con * (6 + 6 * delta) /
            (6 - 7 * delta);

  /* total volume scattering at 90 degree */
  beta90sw = beta_df + beta_cf;
  bsw      = 8 * M_PI / 3 * beta90sw * (2+delta) / (1+delta);


  //if ( lambda == 490 ) {
  //  FILE *fp2;
  //  fp2 = fopen("test2.txt", "a");
  //  fprintf(fp2, " %.9lf %.9lf %.9lf %.9lf\n",bsw/2, beta90sw, beta_df, beta_cf);
  //  fprintf(fp2, " %e %.9lf %e %.9lf\n",Kbz, Tk, IsoComp, DFRI);
  //  fprintf(fp2, " %.14lf %.14lf %.14lf \n",nsw,n_wat2,base);

//  fprintf(fp2, " %.9lf %e\n",nsw, flu_con);
//  fprintf(fp2, " %lf %lf %e %e %e %e\n",S, M0, dnswds, density, dlnawds, Na);
//  fprintf(fp2, " %lf\n",DFRI);
//}


  return(bsw);

}
