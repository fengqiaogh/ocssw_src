/*
 * get_atrem_corl1v3.c
 *
 *  ATmospheric REMoval (ATREM)
 *
 *  Created on: Feb 19, 2015
 *      Author: Rick Healy (richard.healy@nasa.gov)
 *      from ATREM Fortran Code by Bo Cai Gao
 *      Comments from Fortran Code unless noted otherwise
 *          Variable names may differ slightly from comments
 *          (mostly changed to lower case)
 *      nbands = NOBS in Fortran code
 *
 * 05/05/2016 - C stand alone (Version 3) completed for the most part - clean up continues...
 * 03/01/2016 - Began converting rest of code to C
 *              Netcdf files created
 * 01/15/2016 - Split the solar and sensor paths directly inside atrem.
 *               specify with atrem_splitpaths.  Only works with atrem_full on.
 *
 * 10/20/2015 - r.healy - added ATREM options to command line
 *                      atrem_opt - select gases for transmittance calculation
 *                      atrem_geom - turn on/off geometry recalculation
 *                                   0 = geometry calculated on error angle limit
 *                                   1 = recalculate every pixel
 *                      atrem_full - turn on/off full calculation (k-dist otherwise)
 *                      atrem_model - select atmospheric model (1-6)
 *                                  - 0= determine from latitude and day of year
 *
 *      *  Notes about water vapor VMRS and related quantities:                *
 *                                                                              *
 *     VAPVRT(60)   - a table containing 60 column vapor values (in unit of cm) *
 *                                                                              *
 *     VAP_SLANT(I) = VAPVRT(I) * 2.0, VAP_SLANT is a new table for containing  *
 *                    two-way total vapor amounts. Here the number "2" can be   *
 *                    changed to other numbers, e.g., 2.5, without major        *
 *                    effects on retrieved water vapor values.                  *
 *                                                                              *
 *     G_VAP(I = 1,..., NL) = true vapor geometric factor for each layer in     *
 *                    the model atmosphere (after adjusting for the elevated    *
 *                    surface.                                                  *
 *                                                                              *
 *     VMRM(I) = VMRM(I)*G_VAP(I). The VMRS are multiplied by the geometrical   *
 *                    factor. We can calculate the vapor transmittance on the   *
 *                    Sun-surface-sensor path by assuming a vertical path in    *
 *                    the model atmosphere with geometric-factor-adjusted VMRS. *
 *                                                                              *
 *     CLMVAP  = vertical column amount from ground to space in model atmosphere*
 *     CLMVAPP = vertical column amount from ground to aircraft or satellite    *
 *                    sensor in model atmosphere                                *
 *     Q       = 2.152E25 = # of molecules above the surface at one atmosphere  *
 *                    (in unit of  molecules/cm**2)                             *
 *                                                                              *
 *     VAP_SLANT_MDL= CLMVAP/COS(SOLZNI) + CLMVAPP/COS(OBSZNI) = total amount   *
 *                    of water vapor in the model atmosphere in the L-shaped    *
 *                    Sun-surface-plane ray path.                               *
 *                                                                              *
 *     G_VAP_EQUIV  = VAP_SLANT_MDL / CLMVAP = the "equivalent" geometrical     *
 *                    factor corresponding to the total slant vapor amount      *
 *                    VAP_SLANT_MDL and the column vapor amount CLMVAP.         *
 *                                                                              *
 *     SSH2O(I) (I = 1, ..., 60) - a pure scaling factor relative to the total  *
 *                    slant vapor amount of VAP_SLANT_MDL, and                  *
 *            SSH2O(I) = VAP_SLANT(I) / VAP_SLANT_MDL                           *
 *                                                                              *
 *     SH2O = one value of SSH2O(I). SH2O is used during generation of the      *
 *            look-up table.                                                    *
 *                                                                              *
 *     VAPTT  = VAP_SLANT_MDL*SH2O, is the absolute total vapor amount on the   *
 *                    L-shaped path corresponding to a spectrum stored in the   *
 *                    look-up table.                                            *
 *                                                                              *
 *     CLMWVP = 0.5*(VAPTTA+VAPTTB)/G_VAP_EQUIV, is the retrieved column water  *
 *                    vapor amount from imaging spectrometer data.              *
 ********************************************************************************
 ********************************************************************************
 ********************************************************************************!*
References: ********************************************************************!*                                         *
!*        Gao, B.-C., K. Heidebrecht, and A. F. H. Goetz, Derivation of scaled  *
!*           surface reflectances from AVIRIS data, Remote Sens. Env., 44,      *
!*           165-178, 1993.                            *
!*        Gao, B.-C., and C. O. Davis, Development of an operational algorithm  *
!*           for removing atmospheric effects from HYDICE and HSI data,         *
!*           in SPIE'96 Conference Proceedings, Vol. 2819, 45-55, 1996.         *
!*        Gao, B.-C., and A. F. H. Goetz, Column atmospheric water vapor and    *
!*           vegetation liquid water retrievals from airborne imaging          *
!*           spectrometer data, J. Geophys. Res., 95, 3549-3564, 1990.         *
!*        Goetz, A. F. H., and M. Herring, The high resolution imaging         *
!*           spectrometer (HIRIS) for Eos, IEEE Trans. Geosci. Remote Sens., 27,*
!*           136-144, 1989.                            *
!*        Goetz, A. F. H., G. Vane, J. Solomon, and B. N. Rock, Imaging        *
!*           spectrometry for Earth remote sensing, Science, 228, 1147-1153,1985*
!*        Green, R. O., and B.-C. Gao, A Proposed Update to the Solar Irradiance*
!*           Spectrum used in LOWTRAN and MODTRAN, in Summaries of the Fourth   *
!*           Annual JPL Airborne Geoscience Workshop, October 25-29, (Editor,   *
!*           R. O. Green), JPL Publ. 93-26, Vol. 1, pp. 81-84, Jet Propul. Lab, *
!*           Pasadena, Calif., 1993.                           *
!*        Kneizys, F. X., E. P. Shettle, L. W. Abreu, J. H. Chetwynd, G. P.     *
!*           Anderson, W. O. Gallery, J. E. A. Selby, and S. A. Clough, Users   *
!*           guide to LOWTRAN7, AFGL-TR-8-0177, Air Force Geophys. Lab.,        *
!*           Bedford, Mass., 1988.                         *
!*        Iqbal, M., An Introduction To Solar Radiation, Academic, San Diego,   *
!*           Calif., 1983.                             *
!*        Malkmus, W., Random Lorentz band model with exponential-tailed S line *
!*           intensity distribution function, J. Opt. Soc. Am., 57, 323-329,1967*
!*        Press, W. H., B. P. Flannery, S. A. Teukolsky, and W. T.  Vetterling, *
!*           Numerical Recipes-The ART of Scientific Computing, Cambridge       *
!*           University Press, 1986.                           *
!*        Rothman, L. S., et al., The HITRAN 2008 molecular spectroscopic       *
!*           database, JQSRT, 110, 533-572, 2009.                               *
!*        Solomon, S., R. W. Portmann, R. W. Sanders, J. S. Daniel, W. Madsen,  *
!*           B. Bartram, and E. G. Dutton, On the role of nitrogen dioxide in   *
!*           the absorption of solar radiation, J. Geophys. Res., 104,          *
!*           12047-12058, 1999.                                                 *
!*        Tanre, D., C. Deroo, P. Duhaut, M. Herman, J. J. Morcrette, J. Perbos,*
!*           and P. Y. Deschamps, Description of a computer code to simulate    *
!*           the satellite signal in the solar spectrum: the 5S code, Int.      *
!*           J. Remote Sens., 11, 659-668, 1990.                       *
!*        Tanre, D., C. Deroo, P. Duhaut, M. Herman, J. J. Morcrette, J. Perbos,*
!*           and P. Y. Deschamps, Simulation of the satellite signal in the     *
!*           solar spectrum (5S), Users' Guide (U. S. T. De Lille, 59655        *
!*           Villeneu d'Ascq, France: Laboratoire d'Optique Atmospherique),     *
!*      1986.                                  *
!*        Thuillier, G., et al., Solar irradiance reference spectra for two     *
!*           solar active levels, Adv. Space Res., 34, 256-261, 2004.           *
!*        Vane, G., R. O. Green, T. G. Chrien, H. T. Enmark, E. G. Hansen, and  *
!*           W. M. Porter, The Airborne Visible/Infrared Imaging Spectrometer,  *
!*           Remote Sens. Env., 44, 127-143, 1993.                              *
!*        Vane, G. (Ed), Airborne visible/infrared imaging spectrometer        *
!*      (AVIRIS), JPL Publ. 87-38, Jet Propul. Lab, Pasadena, Calif., 1987.*
!*        Vermote, E., D. Tanre, J. L. Deuze, M. Herman, and J. J. Morcrette,   *
!*           Second simulation of the satellite signal in the solar spectrum    *
!*           (6S), 6S User's Guide Version 1, NASA-GSFC, Greenbelt, Maryland,   *
!*           134 pages, 1994.                                                   *
!*                                         *
!********************************************************************************
 ********************************************************************************
 */

#include "atrem_corl1v3.h"
#include <sensorInfo.h>
#include <timeutils.h>
#include <math.h>
#define BUFSZ 100
#define MINDEGCHANGE 55 // Minimum change (degrees squared) in sum of zenith and azimuth angle squared
// for both solar and sensor angles before recalculating tran_table (transmittance table)

int get_atrem_cor(l1str *l1rec, int32_t ip, float *rhot, float *tg_tot, double *tg_sol, double *tg_sen) {
    int32_t sensorID = l1rec->l1file->sensorID;

    static paramstr P;
    int i, j, nb, modnum = input->atrem_model;
    static int firstCall = 1;
    int32_t nbands = l1rec->l1file->nbands;
    static double prev_ddeg = -999, prev_max_senz, prev_min_senz;
    static float **angle_limit, *ang_senz, *ang_solz;
    static int n_senz, n_solz, so_ang, se_ang;
    static int prev_modnum;
    float limitang, *anglelimit;

    //Initialize input paramters


    if (firstCall == 1) { //|| prev_dist > MINDEGCHANGE) {
        //INIT
        init_atrem(sensorID, &P, l1rec, nbands);
        prev_modnum = P.model;

        if (P.dogeom == 0) {
            printf("ATREM: Reading geometry limiting angles\n");
            if (get_angle_limits(&anglelimit, &ang_senz, &ang_solz, &n_senz, &n_solz)) {
                printf("-E- %s line %d : Error reading angle_limit file.\n",
                        __FILE__, __LINE__);
                exit(FATAL_ERROR);
            }

            if ((angle_limit = (float **) calloc(n_senz, sizeof (float *))) == NULL) {
                printf("-E- : Error allocating memory to tindx\n");
                exit(FATAL_ERROR);
            }

            for (j = 0; j < n_senz; j++) {
                if ((angle_limit[j] = (float *) calloc(n_solz, sizeof (float))) == NULL) {
                    printf("-E- : Error allocating memory to tindx\n");
                    exit(FATAL_ERROR);
                }
                for (i = 0; i < n_solz; i++) {
                    angle_limit[j][i] = *(anglelimit + j * (n_solz) + i);
                    //printf("RJH:2: angle_limit: %f %f %f\n",ang_solz[i],ang_senz[j],angle_limit[j][i]);

                }
            }

            prev_max_senz = l1rec->senz[ip];
            prev_min_senz = l1rec->senz[ip];
        }
    }

    if (P.dogeom == 0) {
        if (l1rec->senz[ip] > prev_max_senz) prev_max_senz = l1rec->senz[ip];
        if (l1rec->senz[ip] < prev_min_senz) prev_min_senz = l1rec->senz[ip];

        prev_ddeg = fabs(prev_max_senz - prev_min_senz);

        limitang = get_current_angle_limit(l1rec->senz[ip], l1rec->solz[ip], &se_ang, &so_ang, angle_limit, ang_senz, ang_solz, n_senz, n_solz);
        //printf("RJH: prev_max_senz=%f prev_min_senz=%f limitang[%d]=%f prev_ddeg=%f senz=%f solz=%f \n",prev_max_senz,prev_min_senz,ip,limitang,prev_ddeg,l1rec->senz[ip],l1rec->solz[ip]);
    }
    if (firstCall == 1 || prev_ddeg >= limitang || P.dogeom != 0) {

        //printf("Calculating Transmittance table for Atrem correction, ip=%d\n",ip);
        //start_time = now();
        //printf("\nBegin get_atrem_cor processing at %s\n\n", ydhmsf(start_time,'L'));


        geometry_l2gen_.splitpaths = 0;
        geometry_l2gen_.senzn_l2 = l1rec->senz[ip];
        geometry_l2gen_.senaz_l2 = l1rec->sena[ip];
        geometry_l2gen_.solzn_l2 = l1rec->solz[ip];
        getinput3_.vrto3 = l1rec->oz[ip];

        //geometry_(); //FORTRAN
        geometry();
        //printf("Processed geometry after %6.0f seconds\n",now()-start_time);

        //init_speccal_(); //FORTRAN
        init_spectral_calculations();
        //printf("Processed init_speccal after %6.0f seconds\n",now()-start_time);

        //tran_table_(); //FORTRAN
        tran_table();
        //printf("Processed tran_table after %6.0f seconds\n",now()-start_time);

        P.nh2o = init_speccal3_.nh2o;
        P.start_ndx[0] = init_speccal6_.ist1 - 1; // Fortran array index starts at 1, C at 0
        P.start_ndx[1] = init_speccal6_.ist2 - 1;
        P.end_ndx[0] = init_speccal6_.ied1 - 1;
        P.end_ndx[1] = init_speccal6_.ied2 - 1;
        P.start_p94 = init_speccal6_.istp94 - 1;
        P.end_p94 = init_speccal6_.iedp94 - 1;
        P.start_ndx[2] = init_speccal7_.ist3 - 1;
        P.start_ndx[3] = init_speccal7_.ist4 - 1;
        P.end_ndx[2] = init_speccal7_.ied3 - 1;
        P.end_ndx[3] = init_speccal7_.ied4 - 1;
        P.start_1p14 = init_speccal7_.ist1p14 - 1;
        P.end_1p14 = init_speccal7_.ied1p14 - 1;
        P.wt1 = init_speccal8_.wt1;
        P.wt2 = init_speccal8_.wt2;
        P.wt3 = init_speccal8_.wt3;
        P.wt4 = init_speccal8_.wt4;
        P.start2 = init_speccal10_.istrt2 - 1;
        P.end2 = init_speccal10_.iend2 - 1;
        P.ncv2 = init_speccal10_.ncv2;
        P.finst2 = init_speccal10_.finst2;
        P.natot = init_speccal11_.natot;
        P.nbtot = init_speccal11_.nbtot;
        P.nctot = init_speccal11_.nctot;
        P.ndtot = init_speccal11_.ndtot;
        P.g_vap_equiv = geometry3_.g_vap_equiv;
        P.r0p94 = tran_table1_.r0p94;
        P.r1p14 = tran_table1_.r1p14;
        P.vaptot = tran_table1_.vaptot;
        P.trntbl = tran_table1_.trntbl;


        //            printf("PARAM: nb: %d %d %d %d\n",P.nb1,P.nb2,P.nb3,P.nb4);
        //            printf("PARAM:nb: %d %d\n",P.nbp94,P.nb1p14);
        //            printf("PARAM:nh2o: %d\n",P.nh2o);
        //            printf("PARAM:nbands: %d\n",P.nbands);
        //            printf("PARAM:startndx: %d %d %d %d\n",P.start_ndx[0],P.start_ndx[1],P.start_ndx[2],P.start_ndx[3]);
        //            printf("PARAM:endndx: %d %d %d %d\n",P.end_ndx[0],P.end_ndx[1],P.end_ndx[2],P.end_ndx[3]);
        //            printf("PARAM:start_p94: %d\n",P.start_p94);
        //            printf("PARAM:end_p94: %d\n",P.end_p94);
        //            printf("PARAM:start_1p14: %d\n",P.start_1p14);
        //            printf("PARAM:end_1p14: %d\n",P.end_1p14);
        //            printf("PARAM:%d\n",P.start2);
        //            printf("PARAM:%d\n",P.end2);
        //
        //            printf("PARAM:ncv2: %d\n",P.ncv2);
        //            printf("PARAM:ntot: %d %d %d %d\n",P.natot,P.nbtot,P.nctot,P.ndtot);
        //            printf("PARAM:wt: %f %f %f %f\n",P.wt1,P.wt2,P.wt3,P.wt4);
        //            printf("PARAM:delta: %f %f\n",P.delta,P.delta2);
        //              printf("PARAM:g_vap_equiv: %f \n",P.g_vap_equiv);
        //            printf("PARAM:ip=%d\n",ip);
        //            printf("PARAM:prev_dist=%f\n",prev_ddeg);
        firstCall = 0;

        prev_max_senz = l1rec->senz[ip];
        prev_min_senz = l1rec->senz[ip];
    }

    //    if (l1rec->iscan != prevscan) {
    //        prevscan = l1rec->iscan;
    //    }
    if (modnum == 0) {
        int16_t year, day;
        double secs;
        unix2yds(l1rec->scantime, &year, &day, &secs);
        P.model = getModelNum(*l1rec->lat, day);
        if (P.model != prev_modnum) {
            nb = init_tpvmr_nc(P.model);
            if (nb <= 0) {
                printf("-E- %s line %d : Atmospheric data could not be initialized.\n",
                        __FILE__, __LINE__);
                exit(FATAL_ERROR);
            }
            //model_adj_(); //FORTRAN
            model_adjust();
            prev_modnum = P.model;
        }
    }
    l1rec->wv[ip] = get_atrem(tg_tot, rhot, &P); //returns water vapor(cm)

    // Calculate the separate path transmittances (solar and sensor) if selected
    if (input->atrem_splitpaths) {
        geometry_l2gen_.splitpaths = 1;
        geometry_l2gen_.water_vapor = l1rec->wv[ip]; //Misunderstanding.  This isn't used.
        geometry_l2gen_.ja = P.ja;
        geometry_l2gen_.jb = P.jb;
        geometry_l2gen_.f1a = P.f1a;
        geometry_l2gen_.f1b = P.f1b;
        geometry_l2gen_.f2a = P.f2a;
        geometry_l2gen_.f2b = P.f2b;

        // Call geometry with the indices from the tran_table
        //geometry_(); //FORTRAN
        geometry();

        // Don't need to call init_speccal since other trace gas transmittances don't depend on water vapor

        // recalculate transmittances based on separate solar/sensor paths using the tran_table indices from the first call
        //tran_table_(); //FORTRAN
        tran_table();

        for (i = 0; i < P.nbands; i++) {
            tg_sol[i] = tran_table_l2gen_.tg_sol[i];
            tg_sen[i] = tran_table_l2gen_.tg_sen[i];
        }

    }
    //printf("Water vapor[%d]=%f\n",ip,l1rec->wv[ip]);
    //printf("Processed get_atrem after %6.0f seconds\n",now()-start_time);


    return (0);
}

float get_atrem(float *tg_tot, float *rhot, paramstr *P) {

    double const1 = 0, const2 = 0, const3 = 0, const4 = 0, const5 = 0, const6 = 0;
    double ratio_94c, ratio_94co, ratio_114c, ratio_114co;
    int32_t i, ja, jb;
    float clmwvp;
    /*
     Arrays related to look-up table:
         VAPTOT: TOTAL SUN-SURFACE-SENSOR PATH WATER VAPOR IN UNITS OF CM
         R0P94 : Tc(0.94 um)/(WT1*Tc(0.86)+WT2*Tc(1.02))
         R1P14 : Tc(1.14 um)/(WT3*Tc(1.05)+WT4*Tc(1.23))

     Calculating 3-channel ratios from an observed spectrum, using a
     look-up table procedure to derive the column amount of water vapor
     and to find the associated simulated transmittance spectrum.
     */
    for (i = P->start_ndx[0]; i <= P->end_ndx[0]; i++) {
        const1 += rhot[i];
    }
    const1 /= P->nb1;

    for (i = P->start_ndx[1]; i <= P->end_ndx[1]; i++) {
        const2 += rhot[i];
    }
    const2 /= P->nb2;

    //      printf("const2=%f nb2=%d istr2=%d ind2=%d\n",const2,P->nb2,P->start_ndx[1],P->end_ndx[1]);

    for (i = P->start_p94; i <= P->end_p94; i++) {
        const3 += rhot[i];
    }
    const3 /= P->nbp94;

    ratio_94co = const3 / ((P->wt1 * const1) + (P->wt2 * const2));
    ratio_94c = ratio_94co;

    if (ratio_94co > 1.0) {
        const1 = 0.0;

        for (i = P->start_ndx[0]; i <= P->end_ndx[0]; i++) {
            const1 += (1.0 / rhot[i]);
        }
        const1 /= P->nb1;

        const2 = 0.0;
        for (i = P->start_ndx[1]; i <= P->end_ndx[1]; i++) {
            const2 += (1.0 / rhot[i]);
        }
        const2 /= P->nb2;
        const3 = 0.0;
        for (i = P->start_p94; i <= P->end_p94; i++) {
            const3 += (1.0 / rhot[i]);
        }
        const3 /= P->nbp94;

        ratio_94c = const3 / ((P->wt1 * const1) + (P->wt2 * const2));
    }

    debug_atrem.rp94 = ratio_94c;

    const4 = 0.0;
    for (i = P->start_ndx[2]; i <= P->end_ndx[2]; i++) {
        const4 += rhot[i];
    }
    const4 /= P->nb3;

    const5 = 0.0;
    for (i = P->start_ndx[3]; i <= P->end_ndx[3]; i++) {
        const5 += rhot[i];
    }
    const5 /= P->nb4;

    const6 = 0.0;
    for (i = P->start_1p14; i <= P->end_1p14; i++) {
        const6 += rhot[i];
    }
    const6 /= P->nb1p14;

    /* DEBUG
     *
     */
    debug_atrem.cst1 = const1;
    debug_atrem.cst2 = const2;
    debug_atrem.cst3 = const3;
    debug_atrem.cst4 = const4;
    debug_atrem.cst5 = const5;
    debug_atrem.cst6 = const6;

    ratio_114co = const6 / ((P->wt3 * const4) + (P->wt4 * const5));
    ratio_114c = ratio_114co;

    if (ratio_114co > 1.0) {

        const4 = 0.0;
        for (i = P->start_ndx[2]; i <= P->end_ndx[2]; i++) {
            const4 += (1.0 / rhot[i]);
        }
        const4 /= P->nb3;
        for (i = P->start_ndx[3]; i <= P->end_ndx[3]; i++) {
            const5 += (1.0 / rhot[i]);
        }
        const5 /= P->nb4;
        const6 = 0.0;
        for (i = P->start_1p14; i <= P->end_1p14; i++) {
            const6 += (1.0 / rhot[i]);
        }
        const6 /= P->nb1p14;
        ratio_114c = const6 / ((P->wt3 * const4) + (P->wt4 * const5));
    }

    debug_atrem.r1p14 = ratio_114c;

    double delta, deltab, fja, fjap1, fjb, fjbp1, vaptta, vapttb, specav, spec450;
    static double *speca, *specb;

    if (!speca) speca = (double*) allocateMemory(P->nbands * sizeof (double), "speca");
    if (!specb) specb = (double*) allocateMemory(P->nbands * sizeof (double), "specb");

    ja = P->nh2o / 2;
    ja = hunt(P->r0p94, P->nh2o, ratio_94c, ja);
    //      printf("xx[%d]=%f xx[%d]=%f JA=%d\n",ja-1,P->r0p94[ja-1],ja,P->r0p94[ja],ja);
    if (ja >= 0 && ja < P->nh2o) {
        delta = P->r0p94[ja + 1] - P->r0p94[ja];
        fja = (P->r0p94[ja + 1] - ratio_94c) / delta;
        fjap1 = (ratio_94c - P->r0p94[ja]) / delta;
        vaptta = fja * P->vaptot[ja] + fjap1 * P->vaptot[ja + 1];
        if (ratio_94co > 1.0) vaptta = -vaptta;
    } else {
        if (ja < 0) vaptta = P->vaptot[ja + 1];
        if (ja > P->nh2o) vaptta = P->vaptot[ja];
    }

    if (ratio_94co <= 1.0) {
        for (i = 0; i < P->nbands; i++) {
            if (ja >= 0 && ja < P->nh2o - 1) {
                speca[i] = fja * P->trntbl[ja][i] + fjap1 * P->trntbl[ja + 1][i];
            } else {
                if (ja < 0) speca[i] = P->trntbl[ja + 1][i];
                if (ja >= P->nh2o - 1) speca[i] = P->trntbl[ja][i];
            }
        }
    }

    if (ratio_94co > 1.0) {
        for (i = 0; i < P->nbands; i++) {
            if (ja >= 0 && ja < P->nh2o - 1) {
                speca[i] = 1.0 / (fja * P->trntbl[ja][i] + fjap1 * P->trntbl[ja + 1][i]);
            } else {
                if (ja < 0) speca[i] = 1.0 / P->trntbl[ja + 1][i];
                if (ja >= NH2OMAXM1) speca[i] = 1.0 / P->trntbl[ja][i];
            }
        }
    }

    jb = ja;

    jb = hunt(&P->r1p14[0], P->nh2o, ratio_114c, jb);
    //printf("RJH: 1p14 JB=%d\n",jb);

    debug_atrem.jac = ja;
    debug_atrem.jbc = jb;

    if (jb >= 0 && jb < P->nh2o - 1) {
        deltab = P->r1p14[jb + 1] - P->r1p14[jb];
        fjb = (P->r1p14[jb + 1] - ratio_114c) / deltab;
        fjbp1 = (ratio_114c - P->r1p14[jb]) / deltab;
        vapttb = fjb * P->vaptot[jb] + fjbp1 * P->vaptot[jb + 1];
        if (ratio_114co > 1.0) vapttb = -vapttb;
    } else {
        if (jb < 0) vapttb = P->vaptot[jb + 1];
        if (jb <= P->nh2o - 1) vapttb = P->vaptot[jb];
    }

    if (ratio_114co <= 1.0) {
        for (i = 0; i < P->nbands; i++) {
            if (jb >= 0 && jb < P->nh2o - 1) {
                specb[i] = fjb * P->trntbl[jb][i] + fjbp1 * P->trntbl[jb + 1][i];
            } else {
                if (jb < 0) specb[i] = P->trntbl[jb + 1][i];
                if (jb >= P->nh2o - 1) specb[i] = P->trntbl[jb][i];
            }
        }
    }
    if (ratio_114co > 1.0) {
        for (i = 0; i < P->nbands; i++) {
            if (jb >= 0 && jb < P->nh2o - 1) {
                specb[i] = 1.0 / (fjb * P->trntbl[jb][i] + fjbp1 * P->trntbl[jb + 1][i]);
            } else {
                if (jb < 0) specb[i] = 1.0 / P->trntbl[jb + 1][i];
                if (jb >= P->nh2o - 1) specb[i] = 1.0 / P->trntbl[jb][i];
            }
        }
    }

    clmwvp = 0.5 * (vaptta + vapttb) / P->g_vap_equiv;
    spec450 = 1.0 - 0.5 * (speca[P->idx450] + specb[P->idx450]); // should be 0.0

    //  Derivation of surface reflectances

    for (i = 0; i < P->nbands; i++) {
        specav = 0.5 * (speca[i] + specb[i]);
        tg_tot[i] = specav + spec450;
        //printf("RJH: Atrem: %d %f YY=%f %f %f\n",i,tg_tot[i],rhot[i],rhot[i]/specav,rhot[i]-rhot[i]/specav);
    }

    P->ja = ja;
    P->jb = jb;
    P->f1a = fja;
    P->f2a = fjap1;
    P->f1b = fjb;
    P->f2b = fjbp1;


    /*
    // Smooth the derived surface reflectance spectra
    // rhot = yy from atrem fortran code
    // tg_tot = vc from atrem fortran code
          int ii,j;
          double truncv;

         if (P->delta2 > P->delta) {

    //           * First, replace radiances <= 0 near the spectral overlapping parts of the
    //           * four AVIRIS spectrometers by radiances in the nearby AVIRIS' channels.

             for (i=P->natot-3; i< P->natot+2; i++) {
                 if (rhot[i] <= 0.0) rhot[i] = rhot[i-1];
             }
             for (i=P->nbtot-3; i< P->nbtot+2; i++) {
                 if (rhot[i] <= 0.0) rhot[i] = rhot[i-1];
             }
             for (i=P->nctot-3; i< P->nctot+2; i++) {
                 if (rhot[i] <= 0.0) rhot[i] = rhot[i-1];
             }
             for (i=P->ndtot-3; i< P->ndtot+2; i++) {
                 if (rhot[i] <= 0.0) rhot[i] = rhot[i-1];
             }
             for (i=P->start2-1; i<P->end2; i++){
                 truncv = 0.0;
                 ii = i - P->ncv2 - 1;
                 for (j=ii; j<i+P->ncv2; j++) {
                     truncv += rhot[j]*P->finst2[j-ii+1];
                 }
                 rhot[i] =  truncv;
             }
          }
     */

    return (clmwvp);

}

/*********************************************************************************
 *                                                                                *
 *  Name: HUNT                                                                    *
 *  Purpose: finds the element in array XX that is closest to value X.  Array AA  *
 *       must be monotonic, either increasing or decreasing.                      *
 *  Parameters: XX  - array to search                                             *
 *              N - number of elements in the array                               *
 *              X - element to search for closest match                           *
 *          JLO - index of the closest matching element                           *
 *  Algorithm: this subroutine was copied from Numerical Recipes                  *
 *  Modified for C and cleaned up by                                              *
 *    Richard Healy SAIC/NASA-GSFC 2/19/2015                                      *
 *    (Who didn't realize it was a Numerical Recipe's function until finished)    *
 *  Globals used: none                                                            *
 *  Global output: none                                                           *
 *  Return Codes: none                                                            *
 *  Special Considerations: none                                                  *
 *                                                                                *
 **********************************************************************************
 *
 */
int32_t hunt(float *xx, int32_t n, double x, int32_t jlo) {
    int32_t inc, jhi, jm;
    int ascnd;

    ascnd = (xx[n - 1] >= xx[0]);

    inc = 1;
    if (jlo >= 0 && jlo < n) {

        if ((x >= xx[jlo]) == ascnd) {
            jhi = jlo + inc;
            while (jhi < n && ((x >= xx[jhi]) == ascnd)) {
                jlo = jhi;
                inc += inc;
                jhi = jlo + inc;
            }
            if (jhi >= n)
                jhi = n;

        } else {
            jhi = jlo;
            jlo = jhi - inc;
            while (jlo >= 0 && ((x < xx[jlo]) == ascnd)) {
                jhi = jlo;
                inc += inc;
                jlo = jhi - inc;
            }
            if (jlo < 0)
                jlo = -1;

        }

    } else {
        jlo = 0;
        jhi = n + 1;
    }
    while (jhi - jlo != 1) {

        jm = (jhi + jlo) / 2;

        if ((x >= xx[jm]) == ascnd) {
            jlo = jm;
        } else {
            jhi = jm;
        }
    }

    if (x == xx[n - 1]) jlo = n - 2;
    if (x == xx[0]) jlo = 0;

    return (jlo);
}

int init_tpvmr(int model) {

    int i, nb;
    printf("ATREM: Initializing ATM for model number = %d\n", model);
    model--;
    if (model < 0 || model > MODELMAX) {
        printf("-E- %sline %d: Invalid ATM Model number\n Value must be between 1 and 7\n: get_atrem_cor3\n",
                __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }
    nb = tpvmr_init1_.tpvmr[0][model];

    for (i = 0; i < nb; i++) {
        getinput3_.h[i] = tpvmr_init1_.tpvmr[1 + (4 * i)][model];
        //Convert the atmospheric pressure from "mb" to "atm.":
        getinput3_.p[i] = tpvmr_init1_.tpvmr[2 + (4 * i)][model] / 1013.;
        getinput3_.t[i] = tpvmr_init1_.tpvmr[3 + (4 * i)][model];
        //Convert the VMR from the ppm unit in the model to absolute unit
        getinput3_.vmr[i] = tpvmr_init1_.tpvmr[4 + (4 * i)][model]*1.0E-06;
    }

    for (i = nb; i < MODELMAX; i++) {
        getinput3_.h[i] = 1000.;
        getinput3_.p[i] = 0.0;
        getinput3_.t[i] = 300;
        getinput3_.vmr[i] = 0.0;
    }
    getinput3_.nb = nb;
    getinput3_.nl = nb - 1;



    return nb;
}

void get_tpvmr(size_t layers, size_t models, int sds_id,
        char filename[FILENAME_MAX], char *varname, float* var_a) {
    size_t start[3], count[3];
    int retval;
    int xid;
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    count[0] = models;
    count[1] = layers;
    count[2] = 0;
    retval = nc_inq_varid(sds_id, varname, &xid);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, filename, varname);
        exit(FATAL_ERROR);
    }
    retval = nc_get_vara_float(sds_id, xid, start, count, var_a);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, filename, varname);
        exit(FATAL_ERROR);
    }
}

int init_tpvmr_nc(int model) {

    int i;
    char *filedir;
    char filename[FILENAME_MAX];
    char *tpvmr_file = "atrem_tpvmr.nc";
    int xid, yid, retval, sds_id;
    static size_t models, layers;
    static float *h, *p, *t, *vmr;

    if (!h) {

        if ((filedir = getenv("OCDATAROOT")) == NULL) {
            printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
            return (-1);
        }

        strcpy(filename, filedir);
        strcat(filename, "/common/");
        strcat(filename, tpvmr_file);
        strcat(filename, "\0");
        printf("ATREM: Initializing arrays for Atrem atmospheric models\n");

        retval = nc_open(filename, NC_NOWRITE, &sds_id);
        if (retval != NC_NOERR) {
            fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
                    __FILE__, __LINE__, filename);
            exit(FATAL_ERROR);
        }
        // Now read in the absorption coefficients
        retval = nc_inq_dimid(sds_id, "n_layers", &xid) \
               + nc_inq_dimlen(sds_id, xid, &layers) \
               + nc_inq_dimid(sds_id, "n_models", &yid)  \
               + nc_inq_dimlen(sds_id, yid, &models);

        if (retval) {
            fprintf(stderr, "-E- %s line %d: nc_inq_dimid(%s) failed.%d %d %d \n",
                    __FILE__, __LINE__, filename, (int) layers, (int) models, retval);
            exit(FATAL_ERROR);
        }

        h = (float *) calloc(layers*models, sizeof (float));
        if (!(h)) {
            fprintf(stderr, "-E- %s line %d: Failed to allocate memory for h.\n",
                    __FILE__, __LINE__);
            exit(FATAL_ERROR);
        }
        p = (float *) calloc(layers*models, sizeof (float));
        if (!(p)) {
            fprintf(stderr, "-E- %s line %d: Failed to allocate memory for p.\n",
                    __FILE__, __LINE__);
            exit(FATAL_ERROR);
        }
        t = (float *) calloc(layers*models, sizeof (float));
        if (!(t)) {
            fprintf(stderr, "-E- %s line %d: Failed to allocate memory for t.\n",
                    __FILE__, __LINE__);
            exit(FATAL_ERROR);
        }
        vmr = (float *) calloc(layers*models, sizeof (float));
        if (!(vmr)) {
            fprintf(stderr, "-E- %s line %d: Failed to allocate memory for vmr.\n",
                    __FILE__, __LINE__);
            exit(FATAL_ERROR);
        }
        get_tpvmr(layers, models, sds_id, filename, "h", h);
        get_tpvmr(layers, models, sds_id, filename, "p", p);
        get_tpvmr(layers, models, sds_id, filename, "t", t);
        get_tpvmr(layers, models, sds_id, filename, "vmr", vmr);
    }

    printf("ATREM: Initializing ATM for model number = %d\n", model);
    model--;

    if (model < 0 || model > MODELMAX) {
        printf("-E- %sline %d: Invalid ATM Model number\n Value must be between 1 and 7\n: get_atrem_cor3\n",
                __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }

    for (i = 0; i < layers; i++) {
        getinput3_.h[i] = h[model * layers + i];
        //Convert the atmospheric pressure from "mb" to "atm.":
        getinput3_.p[i] = p[model * layers + i] / 1013.;
        getinput3_.t[i] = t[model * layers + i] / 1013.;
        //Convert the VMR from the ppm unit in the model to absolute unit
        getinput3_.vmr[i] = vmr[model * layers + i]*1.0e-06;
    }

    for (i = layers; i < MODELMAX; i++) {
        getinput3_.h[i] = 1000.;
        getinput3_.p[i] = 0.0;
        getinput3_.t[i] = 300;
        getinput3_.vmr[i] = 0.0;
    }
    getinput3_.nb = layers;
    getinput3_.nl = layers - 1;



    return layers;
}

int getModelNum(float lat, int32_t day) {
    //  Determine which atmospheric model to use
    //     1 = tropical
    //     2 = mid latitude summer
    //     3 = mid latitude winter
    //     4 = subarctic summer
    //     5 = subarctic winter
    //     6 = US standard 1962
    int16 mon = (int) day / 31 + 1; // month of year (no need for perfection..at least according to the sea surface salinity reference algorithm)

    if (fabs(lat) < 30) return (1);
    else {
        switch (mon) {
        case 12:
        case 1:
        case 2:
            if (lat < 60 && lat > 30) return (3);
            if (lat < -30 && lat >-60) return (2);
            if (lat > 60) return (5);
            return (4);
            break;
        case 6:
        case 7:
        case 8:
            if (lat < 60 && lat > 30) return (2);
            if (lat < -30 && lat >-60) return (3);
            if (lat > 60) return (4);
            return (5);
            break;
        default:
            return (6);

        }


    }


}

int init_atrem(int32_t sensorID, paramstr *P, l1str *l1rec, int32_t nbands) {

    int32_t nb, atrem_opt = input->atrem_opt, atrem_full = input->atrem_full, atrem_geom = input->atrem_geom;
    int32_t atrem_splitpaths = input->atrem_splitpaths, atrem_model = input->atrem_model, gas_opt = input->gas_opt;
    float *fwhm, *xppp;
    float nwave;
    int i;
    char *filedir;
    char filename[FILENAME_MAX];
    int *model, flag_gas;

    model = (int *) calloc(1, sizeof (int));
    fwhm = (float *) calloc(1, sizeof (float));
    xppp = (float *) calloc(1, sizeof (float));

    if ((filedir = getenv("OCDATAROOT")) == NULL) {
        printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
        return (-1);
    }

    strcpy(filename, filedir);
    strcat(filename, "/");
    strcat(filename, sensorId2SensorDir(sensorID));
    strcat(filename, "/\0");

    //INIT - get_input_ must be called before any common block structures are filled in.
    //       It will also fill some of the parameters with default values
    //get_input_();

    strcpy(input_l2gen_.filename, filename);
    input_l2gen_.dln = strlen(filename);
    //printf("filedir=%s Filename=%s len=%d\n",filedir, filename,input_l2gen_.dln);

    if (!(fwhm = l1rec->l1file->fwhm)) {
        nwave = rdatreminfo(sensorID, 0, "fwhm", (void **) &fwhm);
    }

    getinput4_.fwhm = (float*) allocateMemory(nbands * sizeof (float), "fwhm");
    getinput4_.wavobs = (float*) allocateMemory(nbands * sizeof (float), "wavobs");

    for (i = 0; i < nbands; i++) {
        getinput4_.wavobs[i] = l1rec->l1file->fwave[i] / 1000.; // Convert nm to microns
        getinput4_.fwhm[i] = fwhm[i]; //fwhm should already be in microns
        //printf("->RJH:fwave(%d)=%f fwhm(%d) = %f \n",i+1,getinput4_.wavobs[i],i,fwhm[i]);
    }

    P->idx450 = bindex_get(450); // get the index nearest the 450nm wavelength
    if (P->idx450 < 0) {
        printf("-E- %s line %d : Unable to determine 450 nm index from spectrum.\n",
                __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }
    //    for (i=0;i<nbands;i++) {
    //        getinput4_.fwhm[i] = fwhm[i];
    // //       printf("->RJH:fwave(%d)=%f fwhm(%d) = %f\n",i+1,getinput4_.wavobs[i],i,fwhm[i]);
    //    }
    /*
    //  Determine which atmospheric model to use
    //  call to get_input sets h,p,t, and vmr arrays to default model 1
    //     1 = tropical
    //     2 = mid latitude summer
    //     3 = mid latitude winter
    //     4 = subarctic summer
    //     5 = subarctic winter
    //     6 = US standard 1962
    //     7 = user defined model

    //INIT
    //Get Atmospheric model to use
     */
    *model = 6;

    if (atrem_model == 0) {
        int16_t year, day;
        double secs;
        unix2yds(l1rec->scantime, &year, &day, &secs);
        *model = getModelNum(*l1rec->lat, day);
    } else {
        if (atrem_model <= 6) {
            *model = atrem_model;
            //   nwave = rdatreminfo(sensorID, 0, "model", (void **)&model);
        } else {
            printf("-E- %s line %d : Invalid atmospheric model, atrem_model = %d. Valid range is 0-6.\n",
                    __FILE__, __LINE__, atrem_model);
            exit(FATAL_ERROR);
        }
    }
    P->model = *model;
    nb = init_tpvmr_nc(P->model);

    if (nb <= 0) {
        printf("-E- %s line %d : Atmospheric data could not be initialized.\n",
                __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }

    getinput5_.nbands = nbands; //l1rec->nbands;
    getinput5_.dlt = 0; //used internally to the fortran routines
    getinput5_.dlt2 = 0; //to determine what to use for the width of the bands
    getinput5_.hsurf = 0; //mean surface elevation
    getinput14_.xppp = 700; //Plane/Satellite altitude (km)

    if (getinput5_.hsurf < 0 || getinput5_.hsurf > getinput3_.h[nb - 2]) {
        printf("-E- %s line %d : Elevation=%f must be less than the maximum elevation in the atmospheric model.\n",
                __FILE__, __LINE__, getinput5_.hsurf);
        exit(FATAL_ERROR);
    }

    nwave = rdatreminfo(sensorID, 0, "xppp", (void **) &xppp);

    if (getinput14_.xppp < getinput5_.hsurf) {
        printf("-E- %s line %d : Sensor altitude=%f must be greater than the bottom surface elevation.\n",
                __FILE__, __LINE__, getinput14_.xppp);
        exit(FATAL_ERROR);

    }
    /*
    //  Get ranges on curve for the two atmospheric windows surrounding the 1.14-um
    //    water vapor absorption feature and the center point of the 1.14-um water
    //    vapor absorption feature.  Enter:
    //         1. the midpoint of third window (0.6-2.5)
    //         2. number of points to average for third window (1-10)
    //         3. the midpoint of fourth window (0.6-2.5)
    //         4. number of points to average for fourth window (1-10)
    //         5. the midpoint of 1.14-um absorption feature (0.6-2.5)
    //         6. the number of points to average for the absorption feature (1-30)
     */
    // This is the default for HICO HS
    float *wndow1, *wndow2, *wp94c, *wndow3, *wndow4, *w1p14c;

    wndow1 = (float *) calloc(1, sizeof (float));
    wndow2 = (float *) calloc(1, sizeof (float));
    wndow3 = (float *) calloc(1, sizeof (float));
    wndow4 = (float *) calloc(1, sizeof (float));
    wp94c = (float *) calloc(1, sizeof (float));
    w1p14c = (float *) calloc(1, sizeof (float));

    getinput6_.wndow1 = 0.705;
    getinput6_.wndow2 = 0.745;
    getinput6_.wp94c = 0.725;
    getinput6_.wndow3 = 0.805;
    getinput6_.wndow4 = 0.845;
    getinput6_.w1p14c = 0.825;

    nwave = rdatreminfo(sensorID, 0, "window1", (void **) &wndow1);
    nwave = rdatreminfo(sensorID, 0, "window2", (void **) &wndow2);
    nwave = rdatreminfo(sensorID, 0, "window3", (void **) &wndow3);
    nwave = rdatreminfo(sensorID, 0, "window4", (void **) &wndow4);
    nwave = rdatreminfo(sensorID, 0, "wp94c", (void **) &wp94c);
    nwave = rdatreminfo(sensorID, 0, "w1p14c", (void **) &w1p14c);

    getinput6_.wndow1 = *wndow1;
    getinput6_.wndow2 = *wndow2;
    getinput6_.wp94c = *wp94c;
    getinput6_.wndow3 = *wndow3;
    getinput6_.wndow4 = *wndow4;
    getinput6_.w1p14c = *w1p14c;

    if (getinput6_.wndow1 < 0.6 || getinput6_.wndow1 > 2.5) {
        fprintf(stderr, "Invalid wavelength position for first atmospheric window " \
                        "in the .94-um region1:%f\nValid values are 0.6-2.5\n", getinput6_.wndow1);
        printf("-E- %s line %d : See stderr output for details.\n",
                __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }
    if (getinput6_.wndow2 < 0.6 || getinput6_.wndow2 > 2.5) {
        fprintf(stderr, "Invalid wavelength position for first atmospheric window " \
                        "in the .94-um region2:%f\nValid values are 0.6-2.5\n", getinput6_.wndow2);
        printf("-E- %s line %d : See stderr output for details.\n",
                __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }
    if (getinput6_.wp94c <= getinput6_.wndow1 || getinput6_.wp94c >= getinput6_.wndow2) {
        fprintf(stderr, "Invalid wavelength position for first atmospheric window " \
                        "in the .94-um region:%f\nValid range is: %f < value < %f \n", getinput6_.wp94c, getinput6_.wndow1, getinput6_.wndow2);
        printf("-E- %s line %d : See stderr output for details.\n",
                __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }
    if (getinput6_.wndow3 < 0.6 || getinput6_.wndow3 > 2.5) {
        fprintf(stderr, "Invalid wavelength position for first atmospheric window " \
                        "in the .94-um region3:%f\nValid values are 0.6-2.5\n", getinput6_.wndow3);
        printf("-E- %s line %d : See stderr output for details.\n",
                __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }
    if (getinput6_.wndow4 < 0.6 || getinput6_.wndow4 > 2.5) {
        fprintf(stderr, "Invalid wavelength position for first atmospheric window " \
                        "in the .94-um region4:%f\nValid values are 0.6-2.5\n", getinput6_.wndow4);
        printf("-E- %s line %d : See stderr output for details.\n",
                __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }
    if (getinput6_.w1p14c <= getinput6_.wndow3 || getinput6_.w1p14c >= getinput6_.wndow4) {
        fprintf(stderr, "Invalid wavelength position for first atmospheric window " \
                        "in the .94-um region:%f\nValid range is: %f < value < %f \n", getinput6_.w1p14c, getinput6_.wndow3, getinput6_.wndow4);
        printf("-E- %s line %d : See stderr output for details.\n",
                __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }

    // Set full_calc to 1 to turn on explicit, full calculation of
    // water vapor correction.  Setting value to 0 turns on speedy
    // h2o absorption calculation and sets all other gas options to 0

    //    int32_t *full_calc;
    //    full_calc = (int32_t *) calloc(1,sizeof(int32_t));
    //    if (atrem_full != 0) {
    //        *full_calc = (atrem_full < 0? 0:1);
    //    } else {
    //        nwave = rdatreminfo(sensorID, 0, "full_calc",    (void **) &full_calc);
    //    }

    /*
     * Allocate memory to global variables
     */

    init_speccal1_.tran_hi_others = (float*) allocateMemory(NP_HI * sizeof (float), "tran_hi_others");
    init_speccal12_.wavln_med = (float*) allocateMemory(NP_MED * sizeof (float), "wavln_med");
    init_speccal12_.wavln_std = (float*) allocateMemory(NP_STD * sizeof (float), "wavln_med");
    init_speccal13_.index_med = (int32_t*) allocateMemory(NP_MED * sizeof (int32_t), "index_med");
    init_speccal13_.wavln_med_index = (float*) allocateMemory(NP_MED * sizeof (float), "wavln_med_index");
    init_speccal15_.ncvhf = (int32_t*) allocateMemory(nbands * sizeof (int32_t), "ncvhf");

    init_speccal13_.tran_med_index = (float **) malloc(NH2OMAX * sizeof (float *));
    for (i = 0; i < NH2OMAX; i++)
        init_speccal13_.tran_med_index[i] = (float*) allocateMemory(NP_MED * sizeof (float), "tran_med_index");

    init_speccal15_.finstr = (float **) malloc(nbands * sizeof (float *));
    for (i = 0; i < nbands; i++)
        init_speccal15_.finstr[i] = (float*) allocateMemory(NINSTR_MAX * sizeof (float), "finstr");

    init_speccal16_.tran_o3_std = (float*) allocateMemory(NO3PT * sizeof (float), "tran_o3_std");
    init_speccal17_.tran_no2_std = (float*) allocateMemory(NO3PT * sizeof (float), "tran_no2_std");

    tran_table1_.vaptot = (float*) allocateMemory(NH2OMAX * sizeof (float), "vaptot");
    tran_table1_.r0p94 = (float*) allocateMemory(NH2OMAX * sizeof (float), "r0p94");
    tran_table1_.r1p14 = (float*) allocateMemory(NH2OMAX * sizeof (float), "r1p14");
    tran_table1_.trntblo = (float*) allocateMemory(nbands * sizeof (float), "trntblo");

    tran_table1_.trntbl = (float **) malloc(NH2OMAX * sizeof (float *));
    tran_table1_.tran_kd = (float **) malloc(NH2OMAX * sizeof (float *));
    tran_table1_.diff_tran = (float **) malloc(NH2OMAX * sizeof (float *));

    for (i = 0; i < NH2OMAX; i++) {
        tran_table1_.trntbl[i] = (float*) allocateMemory(nbands * sizeof (float), "trntbl");
        tran_table1_.tran_kd[i] = (float*) allocateMemory(nbands * sizeof (float), "tran_kd");
        tran_table1_.diff_tran[i] = (float*) allocateMemory(nbands * sizeof (float), "diff_tran");
    }

    tran_table_l2gen_.tg_sol = (float*) allocateMemory(nbands * sizeof (float), "tg_sol");
    tran_table_l2gen_.tg_sen = (float*) allocateMemory(nbands * sizeof (float), "tg_sen");
    tran_table_l2gen_.tg_solo = (float*) allocateMemory(nbands * sizeof (float), "tg_solo");
    tran_table_l2gen_.tg_seno = (float*) allocateMemory(nbands * sizeof (float), "tg_seno");


    if (atrem_splitpaths) {

        //allocate memory

        tran_tables_.tran_hi_sa [0] = (float*) allocateMemory(NP_HI * sizeof (float), "tran_hi_sa");
        tran_tables_.tran_hi_sa [1] = (float*) allocateMemory(NP_HI * sizeof (float), "tran_hi_sa");
        tran_tables_.tran_hi_sap1[0] = (float*) allocateMemory(NP_HI * sizeof (float), "tran_hi_sap1");
        tran_tables_.tran_hi_sap1[1] = (float*) allocateMemory(NP_HI * sizeof (float), "tran_hi_sap1");
        tran_tables_.tran_hi_sb [0] = (float*) allocateMemory(NP_HI * sizeof (float), "tran_hi_sb");
        tran_tables_.tran_hi_sb [1] = (float*) allocateMemory(NP_HI * sizeof (float), "tran_hi_sb");
        tran_tables_.tran_hi_sbp1[0] = (float*) allocateMemory(NP_HI * sizeof (float), "tran_hi_sbp1");
        tran_tables_.tran_hi_sbp1[1] = (float*) allocateMemory(NP_HI * sizeof (float), "tran_hi_sbp1");

        tran_tables1_.tran_med_index_sa_sol = (float*) allocateMemory(NP_MED * sizeof (float), "tran_med_index_sa_sol");
        tran_tables1_.tran_med_index_sa_sen = (float*) allocateMemory(NP_MED * sizeof (float), "tran_med_index_sa_sen");
        tran_tables1_.tran_med_index_sap1_sol = (float*) allocateMemory(NP_MED * sizeof (float), "tran_med_index_sap1_sol");
        tran_tables1_.tran_med_index_sap1_sen = (float*) allocateMemory(NP_MED * sizeof (float), "tran_med_index_sap1_sen");
        tran_tables1_.tran_med_index_sb_sol = (float*) allocateMemory(NP_MED * sizeof (float), "tran_med_index_sb_sol");
        tran_tables1_.tran_med_index_sb_sen = (float*) allocateMemory(NP_MED * sizeof (float), "tran_med_index_sb_sen");
        tran_tables1_.tran_med_index_sbp1_sol = (float*) allocateMemory(NP_MED * sizeof (float), "tran_med_index_sbp1_sol");
        tran_tables1_.tran_med_index_sbp1_sen = (float*) allocateMemory(NP_MED * sizeof (float), "tran_med_index_sbp1_sen");

        tran_tables1_.tran_med_sa_sol = (float*) allocateMemory(NP_MED * sizeof (float), "tran_med_sa_sol");
        tran_tables1_.tran_med_sa_sen = (float*) allocateMemory(NP_MED * sizeof (float), "tran_med_sa_sen");
        tran_tables1_.tran_med_sap1_sol = (float*) allocateMemory(NP_MED * sizeof (float), "tran_med_sap1_sol");
        tran_tables1_.tran_med_sap1_sen = (float*) allocateMemory(NP_MED * sizeof (float), "tran_med_sap1_sen");
        tran_tables1_.tran_med_sb_sol = (float*) allocateMemory(NP_MED * sizeof (float), "tran_med_sb_sol");
        tran_tables1_.tran_med_sb_sen = (float*) allocateMemory(NP_MED * sizeof (float), "tran_med_sb_sen");
        tran_tables1_.tran_med_sbp1_sol = (float*) allocateMemory(NP_MED * sizeof (float), "tran_med_sbp1_sol");
        tran_tables1_.tran_med_sbp1_sen = (float*) allocateMemory(NP_MED * sizeof (float), "tran_med_sbp1_sen");

        tran_tables1_.tran_std_sa_sol = (float*) allocateMemory(NP_STD * sizeof (float), "tran_std_sa_sol");
        tran_tables1_.tran_std_sa_sen = (float*) allocateMemory(NP_STD * sizeof (float), "tran_std_sa_sen");
        tran_tables1_.tran_std_sap1_sol = (float*) allocateMemory(NP_STD * sizeof (float), "tran_std_sap1_sol");
        tran_tables1_.tran_std_sap1_sen = (float*) allocateMemory(NP_STD * sizeof (float), "tran_std_sap1_sen");
        tran_tables1_.tran_std_sb_sol = (float*) allocateMemory(NP_STD * sizeof (float), "tran_std_sb_sol");
        tran_tables1_.tran_std_sb_sen = (float*) allocateMemory(NP_STD * sizeof (float), "tran_std_sb_sen");
        tran_tables1_.tran_std_sbp1_sol = (float*) allocateMemory(NP_STD * sizeof (float), "tran_std_sbp1_sol");
        tran_tables1_.tran_std_sbp1_sen = (float*) allocateMemory(NP_STD * sizeof (float), "tran_std_sbp1_sen");

        if (atrem_full == 0) {
            printf("ATREM: Turning on atrem_full because you selected atrem_splitpaths\n");
            atrem_full = 1;
        }
    }
    if (atrem_full > 0)
        printf("ATREM: Warning : full_calc !=0. Atrem will calculate transmittance table for every pixel\n");

    //    getinput5_.full_calc = *full_calc;

    getinput5_.full_calc = atrem_full;

    //    int32_t *dogeom;
    //    dogeom = (int32_t *) calloc(1,sizeof(int32_t));
    //    if (atrem_geom != 0) {
    //        *dogeom = (atrem_geom < 0? 0:1);
    //    } else {
    //        nwave = rdatreminfo(sensorID, 0, "dogeom",    (void **) &dogeom);
    //    }
    //
    //    P->dogeom        = *dogeom;

    if (atrem_geom > 0)
        printf("ATREM: Warning : dogeom !=0. Geometry will be calculated every pixel\n");

    P->dogeom = atrem_geom;
    //Default for HICO HS
    int32_t *nb1, *nb2, *nb3, *nb4, *nbp94, *nb1p14;
    nb1 = (int32_t *) calloc(1, sizeof (int32_t));
    nb2 = (int32_t *) calloc(1, sizeof (int32_t));
    nb3 = (int32_t *) calloc(1, sizeof (int32_t));
    nb4 = (int32_t *) calloc(1, sizeof (int32_t));
    nbp94 = (int32_t *) calloc(1, sizeof (int32_t));
    nb1p14 = (int32_t *) calloc(1, sizeof (int32_t));

    getinput7_.nb1 = 3;
    getinput7_.nb2 = 3;
    getinput7_.nb3 = 3;
    getinput7_.nb4 = 3;
    getinput7_.nbp94 = 5;
    getinput7_.nb1p14 = 5;

    nwave = rdatreminfo(sensorID, 0, "nb1", (void **) &nb1);
    nwave = rdatreminfo(sensorID, 0, "nb2", (void **) &nb2);
    nwave = rdatreminfo(sensorID, 0, "nb3", (void **) &nb3);
    nwave = rdatreminfo(sensorID, 0, "nb4", (void **) &nb4);
    nwave = rdatreminfo(sensorID, 0, "nbp94", (void **) &nbp94);
    nwave = rdatreminfo(sensorID, 0, "nb1p14", (void **) &nb1p14);

    getinput7_.nb1 = *nb1;
    getinput7_.nb2 = *nb2;
    getinput7_.nb3 = *nb3;
    getinput7_.nb4 = *nb4;
    getinput7_.nbp94 = *nbp94;
    getinput7_.nb1p14 = *nb1p14;

    if (getinput7_.nb1 < 1 || getinput7_.nb1 > 50) {
        fprintf(stderr, "Invalid number of channels for first wavelength position in " \
                "the .94-um region:%d\nValid values are 1-50\n", getinput7_.nb1);
        printf("-E- %s line %d : See stderr output for details.\n",
                __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }

    if (getinput7_.nb2 < 1 || getinput7_.nb2 > 50) {
        fprintf(stderr, "Invalid number of channels for first wavelength position in " \
                "the .94-um region:%d\nValid values are 1-50\n", getinput7_.nb2);
        printf("-E- %s line %d : See stderr output for details.\n",
                __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }
    if (getinput7_.nbp94 < 1 || getinput7_.nbp94 > 90) {
        fprintf(stderr, "Invalid number of channels for first wavelength position in " \
                "the .94-um region:%d\nValid values are 1-90\n", getinput7_.nbp94);
        printf("-E- %s line %d : See stderr output for details.\n",
                __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }
    if (getinput7_.nb3 < 1 || getinput7_.nb3 > 50) {
        fprintf(stderr, "Invalid number of channels for first wavelength position in " \
                "the .94-um region:%d\nValid values are 1-50\n", getinput7_.nb3);
        printf("-E- %s line %d : See stderr output for details.\n",
                __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }
    if (getinput7_.nb4 < 1 || getinput7_.nb4 > 50) {
        fprintf(stderr, "Invalid number of channels for first wavelength position in " \
                "the .94-um region:%d\nValid values are 1-50\n", getinput7_.nb4);
        printf("-E- %s line %d : See stderr output for details.\n",
                __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }
    if (getinput7_.nb1p14 < 1 || getinput7_.nb1p14 > 110) {
        fprintf(stderr, "Invalid number of channels for first wavelength position in " \
                "the .94-um region:%d\nValid values are 1-90\n", getinput7_.nbp94);
        printf("-E- %s line %d : See stderr output for details.\n",
                __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }

    //INIT
    int32_t *h2o, *co2, *o3, *n2o, *co, *ch4, *o2, *no2;
    h2o = (int32_t *) calloc(1, sizeof (int32_t));
    co2 = (int32_t *) calloc(1, sizeof (int32_t));
    o3 = (int32_t *) calloc(1, sizeof (int32_t));
    n2o = (int32_t *) calloc(1, sizeof (int32_t));
    co = (int32_t *) calloc(1, sizeof (int32_t));
    ch4 = (int32_t *) calloc(1, sizeof (int32_t));
    o2 = (int32_t *) calloc(1, sizeof (int32_t));
    no2 = (int32_t *) calloc(1, sizeof (int32_t));

    *h2o = 1;

    getinput1_.h2o = 1;
    getinput1_.co2 = 0;
    getinput1_.o3 = 0;
    getinput1_.n2o = 0;
    getinput1_.co = 0;
    getinput1_.ch4 = 0;
    getinput1_.o2 = 0;
    getinput1_.no2 = 0;

    //    nwave = rdatreminfo(sensorID, 0, "h2o", (void **) &h2o);
    nwave = rdatreminfo(sensorID, 0, "co2", (void **) &co2);
    nwave = rdatreminfo(sensorID, 0, "o3", (void **) &o3);
    nwave = rdatreminfo(sensorID, 0, "n2o", (void **) &n2o);
    nwave = rdatreminfo(sensorID, 0, "co", (void **) &co);
    nwave = rdatreminfo(sensorID, 0, "ch4", (void **) &ch4);
    nwave = rdatreminfo(sensorID, 0, "o2", (void **) &o2);
    nwave = rdatreminfo(sensorID, 0, "no2", (void **) &no2);

    getinput1_.co2 = *co2;
    getinput1_.o3 = *o3;
    getinput1_.n2o = *n2o;
    getinput1_.co = *co;
    getinput1_.ch4 = *ch4;
    getinput1_.o2 = *o2;
    getinput1_.no2 = *no2;

    //    The typical ozone amount is: 0.28-0.55 (atm-cm). The built-in NO2
    //    column amount is 5.0E+15 molecules/cm^2.
    if (atrem_opt > 0) {
        getinput1_.co2 = (atrem_opt & ATREM_CO2) > 0;
        getinput1_.o3 = (atrem_opt & ATREM_O3) > 0;
        getinput1_.n2o = (atrem_opt & ATREM_N2O) > 0;
        getinput1_.co = (atrem_opt & ATREM_CO) > 0;
        getinput1_.ch4 = (atrem_opt & ATREM_CH4) > 0;
        getinput1_.o2 = (atrem_opt & ATREM_O2) > 0;
        getinput1_.no2 = (atrem_opt & ATREM_NO2) > 0;

    }

    printf("ATREM: Gas Options:H2O:1 CO2:%d O3:%d N2O:%d CO:%d CH4:%d O2:%d NO2:%d\n",
            getinput1_.co2, getinput1_.o3, getinput1_.n2o, getinput1_.co, getinput1_.ch4, getinput1_.o2, getinput1_.no2);

    // Now check to make sure the same gas option isn't selected from the l2gen gas option.

    flag_gas = 0;
    if (gas_opt & H2O_BIT) {
        printf("ATREM: cannot be used with gas_opt bit mask=%d (H2O)\n", H2O_BIT);
        flag_gas = 1;
    }
    if (getinput1_.no2 == 1 && (gas_opt & NO2_BIT)) {
        printf("ATREM: cannot be used with gas_opt bit mask=%d (NO2)\n", NO2_BIT);
        flag_gas = 1;

    }
    if (getinput1_.co2 == 1 && (gas_opt & CO2_BIT)) {
        printf("ATREM: cannot be used with gas_opt bit mask=%d (CO2)\n", CO2_BIT);
        flag_gas = 1;
    }
    if (getinput1_.o3 == 1 && (gas_opt & O3_BIT)) {
        printf("ATREM: cannot be used with gas_opt bit mask=%d (O3)\n", O3_BIT);
        flag_gas = 1;
    }

    if (flag_gas) {
        printf("Error: Conflict using ATREM (gas_opt=16 bitmask) with atrem_opt=%d and gas_opt=%d.  " \
                "\nPlease resolve. Refer to command line options for atrem_opt and gas_opt.\n", atrem_opt, gas_opt);
        exit(1);
    }

    float *vrto3, *sno2;
    vrto3 = (float *) calloc(1, sizeof (float));
    sno2 = (float *) calloc(1, sizeof (float));

    getinput3_.vrto3 = 0.34; //total column ozone amount (atm-cm)
    getinput3_.sno2 = 1.0; //NO2 scaling factor (to 5.E+15 molecules/cm^2)

    nwave = rdatreminfo(sensorID, 0, "vrto3", (void **) &vrto3);
    nwave = rdatreminfo(sensorID, 0, "sno2", (void **) &sno2);

    getinput3_.vrto3 = *vrto3; //total column ozone amount (atm-cm)
    getinput3_.sno2 = *sno2; //NO2 scaling factor (to 5.E+15 molecules/cm^2)

    //model_adj_();   //FORTRAN
    model_adjust();

    P->nb1 = getinput7_.nb1;
    P->nb2 = getinput7_.nb2;
    P->nb3 = getinput7_.nb3;
    P->nb4 = getinput7_.nb4;
    P->nbp94 = getinput7_.nbp94;
    P->nb1p14 = getinput7_.nb1p14;
    P->nbands = getinput5_.nbands;
    P->delta = getinput5_.dlt;
    P->delta2 = getinput5_.dlt2;

    return nwave;
}

int get_angle_limits(float **anglelimit, float **insenz, float **insolz, int *n_senz, int *n_solz) {

    char filename[FILENAME_MAX];
    char *infile = "atrem_angle_limit.nc";
    char *filedir;

    char name[H4_MAX_NC_NAME];
    char sdsname[H4_MAX_NC_NAME];
    int ncid;
    int32 sds_id;
    int status;
    nc_type rh_type; /* variable type */
    int dimids[H4_MAX_VAR_DIMS]; /* dimension IDs */
    int ndims;
    int natts; /* number of attributes */
    size_t length;

    static float *senz, *solz;
    static float *angle_limit;

    if ((filedir = getenv("OCDATAROOT")) == NULL) {
        printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
        return (-1);
    }
    strcpy(filename, filedir);
    strcat(filename, "/common/");
    strcat(filename, infile);
    strcat(filename, "\0");

    printf("ATREM: Reading Angle_Limit FILE=%s\n", filename);

    if (nc_open(filename, NC_NOWRITE, &ncid) == NC_NOERR) {

        strcpy(sdsname, "senz");

        status = nc_inq_varid(ncid, sdsname, &sds_id);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                    __FILE__, __LINE__, sdsname, infile);
            exit(1);
        }

        status = nc_inq_var(ncid, sds_id, 0, &rh_type, &ndims, dimids,
                &natts);

        if (nc_inq_dimlen(ncid, dimids[0], &length) != NC_NOERR) {
            nc_inq_dim(ncid, dimids[0], name, &length);
            fprintf(stderr,
                    "-E- %s line %d: could not get size of dimension \"%s\" in netCDF File.\n",
                    __FILE__, __LINE__, name);
            exit(1);
        }

        *n_senz = length;

        if ((senz = (float *) calloc(*n_senz, sizeof (float))) == NULL) {
            printf("-E- : Error allocating memory to senz\n");
            exit(FATAL_ERROR);
        }

        if (nc_get_var(ncid, sds_id, senz) != NC_NOERR) {
            fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                    __FILE__, __LINE__, sdsname, infile);
            exit(1);
        }

        strcpy(sdsname, "solz");

        status = nc_inq_varid(ncid, sdsname, &sds_id);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                    __FILE__, __LINE__, sdsname, infile);
            exit(1);
        }

        status = nc_inq_var(ncid, sds_id, 0, &rh_type, &ndims, dimids,
                &natts);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                    __FILE__, __LINE__, sdsname, infile);
            exit(1);
        }

        if (nc_inq_dimlen(ncid, dimids[0], &length) != NC_NOERR) {
            nc_inq_dim(ncid, dimids[0], name, &length);
            fprintf(stderr,
                    "-E- %s line %d: could not get size of dimension \"%s\" in netCDF File.\n",
                    __FILE__, __LINE__, name);
            exit(1);
        }

        *n_solz = length;

        if ((solz = (float *) calloc(*n_solz, sizeof (float))) == NULL) {
            printf("-E- : Error allocating memory to solz\n");
            exit(FATAL_ERROR);
        }

        if (nc_get_var(ncid, sds_id, solz) != NC_NOERR) {
            fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                    __FILE__, __LINE__, sdsname, infile);
            exit(1);
        }

        strcpy(sdsname, "angle_limit");

        status = nc_inq_varid(ncid, sdsname, &sds_id);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                    __FILE__, __LINE__, sdsname, infile);
            exit(1);
        }

        status = nc_inq_var(ncid, sds_id, 0, &rh_type, &ndims, dimids,
                &natts);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                    __FILE__, __LINE__, sdsname, infile);
            exit(1);
        }
        if (ndims != 2) {
            fprintf(stderr, "-E- %s line %d:  Wrong number of dimensions for %s.  Need 2 got %d.\n",
                    __FILE__, __LINE__, sdsname, ndims);
            exit(1);

        }
        //            printf("ANGLE_LIMIT: %s n_senz=%d n_solz=%d\n",infile,*n_senz,*n_solz);

        if ((angle_limit = (float *) calloc((*n_senz) * (*n_solz), sizeof (float))) == NULL) {
            printf("-E- : Error allocating memory to angle_limit\n");
            exit(FATAL_ERROR);
        }

        //            for (i=0; i< *n_senz;i++) {
        //                if ( (anglelimit[i] = (float *)calloc((*n_solz) ,sizeof(float))) == NULL) {
        //                    printf("-E- : Error allocating memory to tindx\n");
        //                    exit(FATAL_ERROR);
        //                }
        //            }

        if (nc_get_var(ncid, sds_id, angle_limit) != NC_NOERR) {
            fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                    __FILE__, __LINE__, sdsname, infile);
            exit(1);
        }

        //            for (j=0;j<*n_senz;j++)
        //                for (i=0;i<*n_solz;i++) {
        //                   //anglelimit[j][i] = *(angle_limit+j*(*n_solz)+i);
        //                    printf("RJH: angle_limt: %f %f %f\n",solz[i],senz[j],*(angle_limit+j*(*n_solz)+i));
        //
        //                }

        *insenz = (float *) senz;
        *insolz = (float *) solz;
        *insolz = (float *) solz;
        *anglelimit = (float *) angle_limit;

    } else {
        fprintf(stderr, "-E- %s line %d:  Error opening infile = %s.\n",
                __FILE__, __LINE__, infile);
        return (1);
    }


    return (0);

}

float get_current_angle_limit(float insenz, float insolz, int *ii, int *jj, float **anglelimit, float *senz, float *solz, int n_senz, int n_solz) {
    int i, j;

    i = *ii;
    j = *jj;

    if (i < 0 || i >= n_senz) i = 0;
    if (j < 0 || j >= n_solz) j = 0;

    if (insenz > senz[i]) {
        while (i < n_senz && insenz > senz[i])
            i++;
        if (i >= n_senz) i = n_senz - 1;
    } else {
        while (i >= 0 && insenz < senz[i])
            i--;
        if (i < 0) i = 0;
    }

    if (insolz > solz[j]) {
        while (j < n_solz && insolz > solz[j])
            j++;
        if (j >= n_solz) j = n_solz - 1;
    } else {
        while (j >= 0 && insolz < solz[j])
            j--;
        if (j < 0) j = 0;
    }

    *ii = i;
    *jj = j;

    return (anglelimit[i][j]);
}

/********************************************************************************
!*                                                 *
!*  Name: FINDMATCH                                *
!*  Purpose: finds the closest match for ELEM in LIST                  *
!*  Parameters:  LIST - array of values to match ELEM to.  Elements is array    *
!*            should increase in value.                    *
!*  Algorithm: linearly compare ELEM to each element.  When ELEM is smaller     *
!*             than the LIST(I), then it is assumed to be closest to LIST(I-1)  *
!*             or LIST(I).  The one that has the smallest absolute difference   *
!*             to ELEM is returned as the closest match.                   *
!*  Globals used: none                                 *
!*  Global output: none                                *
!*  Return Codes: the closest matching element index                   *
!*  Special Considerations: none                               *
!*                                         *
!********************************************************************************/

int32_t findMatch(float *list, int32_t nbands, float elem) {
    int i;
    float diff1, diff2;

    for (i = 0; i < nbands; i++)
        if (list[i] > elem) break;

    diff1 = fabs(list[i - 1] - elem);
    diff2 = fabs(list[i] - elem);

    if (diff1 < diff2)
        return (i - 1);
    else
        return (i);

}

/*
!********************************************************************************
!*                                                 *
!*  Name: channelRatio                                *
!*  Purpose: Calculate 3-channel ratios.                           *
!*  Parameters: none                                   *
!*  Algorithm: The 0.94-um water vapor absorption channel is ratioed against    *
!*             the linear combination of two window channels near 0.86 and      *
!*             1.03 um to obtain one channel ratio for the 0.94-um band.           *
!*             Similar calculation is done for the 1.14-um water vapor band.    *
!*  Globals used: NB1,NB2,NBP94,NB3,NB4,NB1P14 - number of points used in       *
!*                          channel ratios for both the .94- and 1.14-um regions*
!*                IST1,IED1,IST2,IED2,ISTP94,IEDP94 - 3-channel ratioing        *
!*                          parameters for the 0.94-um water vapor band        *
!*                IST3,IED3,IST4,IED4,IST1P14,IED1P14 - 3-channel ratioing      *
!*                          parameters for the 1.14-um water vapor band.           *
!*       WT1,WT2,WT3,WT4,JA - Relative weights for the four window     *
!*                          channels used in channel-ratioing calculations. JA  *
!*             is an output parameter from a table searching       *
!*             routine.                        *
!*       TRNCAL -  Atmospheric transmittance spectra.              *
!*  Global output:R094,R114 - 3-channel ratio values for the 0.94- and 1.14-um  *
!*                          water vapor bands.                     *
!*  Return Codes: none                                 *
!*  Special Considerations: none                               *
!*                                         *
!********************************************************************************
 */
void channelRatio() {

    float const1[NH2OMAX], const2[NH2OMAX], const3[NH2OMAX], const4[NH2OMAX], const5[NH2OMAX], const6[NH2OMAX];
    int32_t nb1 = getinput7_.nb1, nb2 = getinput7_.nb2, nbp94 = getinput7_.nbp94, \
              nb3 = getinput7_.nb3, nb4 = getinput7_.nb4, nb1p14 = getinput7_.nb1p14;
    int32_t ist1 = init_speccal6_.ist1, ied1 = init_speccal6_.ied1, ist2 = init_speccal6_.ist2, \
              ied2 = init_speccal6_.ied2, istp94 = init_speccal6_.istp94, iedp94 = init_speccal6_.iedp94;
    int32_t ist3 = init_speccal7_.ist3, ied3 = init_speccal7_.ied3, ist4 = init_speccal7_.ist4, \
              ied4 = init_speccal7_.ied4, ist1p14 = init_speccal7_.ist1p14, ied1p14 = init_speccal7_.ied1p14;
    float wt1 = init_speccal8_.wt1, wt2 = init_speccal8_.wt2, wt3 = init_speccal8_.wt3, wt4 = init_speccal8_.wt4;

    int i, j;

    // Calculate average of spectra over window and water vapor absorption regions.

    for (j = 0; j < NH2OMAX; j++) {
        const1[j] = 0.0;
        const2[j] = 0.0;
        const3[j] = 0.0;
        const4[j] = 0.0;
        const5[j] = 0.0;
        const6[j] = 0.0;

        for (i = ist1 - 1; i < ied1; i++)
            const1[j] = const1[j] + tran_table1_.trntbl[j][i];

        const1[j] /= nb1;

        for (i = ist2 - 1; i < ied2; i++)
            const2[j] = const2[j] + tran_table1_.trntbl[j][i];
        const2[j] /= nb2;

        for (i = istp94 - 1; i < iedp94; i++)
            const3[j] = const3[j] + tran_table1_.trntbl[j][i];
        const3[j] /= nbp94;

        tran_table1_.r0p94[j] = const3[j] / (wt1 * const1[j] + wt2 * const2[j]);

        for (i = ist3 - 1; i < ied3; i++)
            const4[j] = const4[j] + tran_table1_.trntbl[j][i];
        const4[j] /= nb3;

        for (i = ist4 - 1; i < ied4; i++)
            const5[j] = const5[j] + tran_table1_.trntbl[j][i];
        const5[j] /= nb4;

        for (i = ist1p14 - 1; i < ied1p14; i++)
            const6[j] = const6[j] + tran_table1_.trntbl[j][i];
        const6[j] /= nb1p14;

        tran_table1_.r1p14[j] = const6[j] / (wt3 * const4[j] + wt4 * const5[j]);
    }


}

void kdistgasabs(float *kcdf, float *abscf, float*waveno, float *wavobs, int32_t nphi, int32_t n_layers, int32_t nbands) {
    //
    // R. Healy 3/2016 - converted from fortran kdist_gas_obs in atrem code which is based on Amir Ibrahim's
    //                   Matlab code
    // Construct k coefficients from gas transmission coefficients (from Amir Ibrahim's Matlab Code)
    //
    //     abscf[nlayers][np_hi]    (IN) :: the absorption coefficients
    //     waveno[np_hi]            (IN) :: wave number of the spectrum
    //     wavobs[nwave]            (IN) :: wavelengths of the instrument detector
    //     np_hi                    (IN) :: Number of high resolution spectral points,
    //                                      covering 3000 cm-1 to 18,000 cm-1 at
    //                                      point spacing of 0.05 cm-1
    //     nlayers                  (IN) :: number of atmospheric layers
    //     nwave                    (IN) :: number of wavelengths of the instrument
    //     kcdf[7][nwave][nlayers] (OUT):: K-Coefficients

    float **alayers, g7[7] = {0, 0.379, 0.6, 0.81, 0.9548, 0.9933, 1};
    float *diflam, *UV_lam, *IR_lam, *UV_dlam, *lam, *wn; //,*ir_dlam
    float *dwave, *dwn, *g_i, *k_i, *logk, *kint[7], **k7[7], *kk;
    int32_t **wavel_window, *binnum;
    float dlam, wmin, wmax, swav;
    float Q = {2.15199993E+25}; // Q=# of molecules above the surface at one atmosphere (molecules/cm**2)
    float AvN = {6.0225E+23};
    int kwavdn, kwavup, nbins;
    int32_t nsamp;
    int k, i, j, n, ndxwuv, ndxwir, ndxtot, iw;

    alayers = (float **) malloc(n_layers * sizeof (float *));
    for (i = 0; i < n_layers; i++) {
        if ((alayers[i] = (float *) malloc(nphi * sizeof (float))) == NULL) {
            printf("%s, %d - E - unable to allocate memory in kdistgasabs\n",
                    __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        for (n = 0; n < nphi; n++)
            alayers[i][n] = *(abscf + i * n_layers + n) * Q * 28.966 / AvN / 1.0e-6;
        //alayers[i][n] = abscf[i][n]*Q*28.966/ AvN / 1.0e-6;
    }

    diflam = (float *) malloc(nbands * sizeof (float));

    dlam = 9999.;
    wmin = 9999;
    wmax = -9999;
    for (i = 0; i < nbands - 1; i++) {
        diflam[i] = wavobs[i + 1] - wavobs[i];
        if (dlam > diflam[i]) dlam = diflam[i];
        if (wmin > wavobs[i]) wmin = wavobs[i];
        if (wmax < wavobs[i]) wmax = wavobs[i];
    }

    diflam[nbands - 1] = diflam[nbands - 2];

    if (wmin > wavobs[nbands - 1]) wmin = wavobs[nbands - 1];
    if (wmax < wavobs[nbands - 1]) wmax = wavobs[nbands - 1];

    ndxwuv = (wmin - 0.0001 - 0.3) / dlam + 1;
    ndxwir = (3.1 - (wmax + 0.0001)) / diflam[nbands - 1] + 1;
    ndxtot = nbands + ndxwuv + ndxwir;

    UV_lam = (float *) malloc(ndxwuv * sizeof (float));
    IR_lam = (float *) malloc(ndxwir * sizeof (float));
    UV_dlam = (float *) malloc(ndxwuv * sizeof (float));
    lam = (float *) malloc(ndxtot * sizeof (float));
    dwave = (float *) malloc(ndxtot * sizeof (float));
    dwn = (float *) malloc(ndxtot * sizeof (float));
    wn = (float *) malloc(ndxtot * sizeof (float));
    wavel_window = (int **) malloc(2 * sizeof (int *));
    binnum = (int32_t *) malloc(n_layers * sizeof (int32_t));

    for (i = 0; i < 2; i++)
        wavel_window[i] = (int *) malloc(ndxtot * sizeof (int));

    k = ndxtot - 1;
    swav = 0.3;

    for (i = 0; i < ndxwuv; i++) {
        UV_lam[i] = swav;
        UV_dlam[i] = diflam[0];
        swav = swav + dlam;
        lam[k] = UV_lam[i];
        dwn[k] = 10000. * dlam / pow(lam[k], 2);
        wn[k] = 10000. / lam[k];
        k--;
    }
    swav = wmin - 1;
    for (i = 0; i < nbands; i++) {
        lam[k] = swav;
        dwn[k] = 10000. * diflam[i] / pow(lam[k], 2);
        wn[k] = 10000. / lam[k];
        swav = swav + diflam[i];
        k--;
    }

    swav = wmax + 0.0001;

    for (i = 0; i < ndxwir; i++) {
        IR_lam[i] = swav;
        //ir_dlam[i] = diflam[nbands];
        swav = swav + diflam[nbands - 1];
        lam[k] = IR_lam[i];
        dwn[k] = 10000. * diflam[nbands - 1] / pow(lam[k], 2);
        wn[k] = 10000. / lam[k];
        k--;
    }

    for (i = 0; i < 7; i++) {
        kint[i] = (float *) malloc(ndxtot * sizeof (float));
        k7[i] = (float **) malloc(ndxtot * sizeof (float *));
        for (k = 0; k < ndxtot; k++)
            k7[i][k] = (float *) malloc(n_layers * sizeof (float));
    }

    iw = nbands - 1; // index for wavelngths for output array kcdf

    for (i = ndxwir; i < ndxtot - ndxwuv; i++) {
        kwavdn = -999;
        kwavup = -999;

        for (k = 0; k < nphi; k++) {
            if (kwavdn < 0 && waveno[k]> (wn[i] - dwn[i] / 2.0)) kwavdn = k;
            if (kwavup < 0 && waveno[k]> (wn[i] + dwn[i] / 2.0)) kwavup = k;
        }

        wavel_window[0][i] = kwavdn;
        wavel_window[1][i] = kwavup;

        if (kwavup < 0 || kwavdn < 0)
            nsamp = 1;
        else
            nsamp = kwavup - kwavdn + 1;

        nbins = nsamp;
        g_i = (float *) malloc((nbins + 1) * sizeof (float));
        k_i = (float *) malloc((nbins + 1) * sizeof (float));
        kk = (float *) malloc((nsamp + 1) * sizeof (float));
        logk = (float *) malloc((nbins + 1) * sizeof (float));

        for (k = 0; k < n_layers; k++) {
            if (kwavdn < 0 || kwavup < 0)
                for (j = 0; j < nsamp + 1; j++)
                    kk[j] = 0;
            else
                n = kwavdn - 1;
            for (j = 0; j < nsamp && n < nphi; j++) {
                kk[j] = alayers[k][n];
                n++;
            }
            binnum[k] = nbins;
            ecdf_(k_i, g_i, &binnum[k], kk, &nsamp);

            for (j = 0; j < nbins + 1; j++)
                logk[j] = log10(k_i[j]);

            for (j = 0; j < 7; j++) {

                kint[j][i] = linterp(g_i, logk, binnum[k], g7[j]);
                k7[j][i][k] = pow(10, kint[j][i]);

                if (i >= ndxwir && i < ndxtot - ndxwuv) {
                    if (isnan(k7[j][i][k])) k7[j][i][k] = 0;

                    //kcdf[j][iw][k] = k7[j][i][k];
                    //                    *(kcdf+j*nbands*7+iw*n_layers+k) = k7[j][i][k];
                    *(kcdf + j * nbands * n_layers + iw * n_layers + k) = k7[j][i][k];
                }
            }

        }

        free(logk);
        free(g_i);
        free(k_i);
        free(kk);

        if (i >= ndxwir && i < ndxtot - ndxwuv) iw--;

    }

    free(UV_lam);
    free(IR_lam);
    free(UV_dlam);
    //free ( ir_dlam );
    free(lam);
    free(dwave);
    free(dwn);
    free(wn);
    free(wavel_window[0]);
    free(wavel_window[1]);
    for (i = 0; i < 7; i++) {
        for (k = 0; k < ndxtot; k++)
            free(k7[i][k]);
        free(k7[i]);
        free(kint[i]);
    }
    for (i = 0; i < n_layers; i++)
        free(alayers[i]);
    free(alayers);
    free(binnum);
}

/*
 * Converted from Fortran by R. Healy (SAIC) 3/2016
!********************************************************************************
!*                                                                              *
!*  Name: model_adjust                                                          *
!*  Purpose: resets the bottom boundary of the input model if the surface       *
!*           elevation is greater than 0, and calculate the column water vapor  *
!*           amount in the selected model.                                      *
!*  Parameters: none                                                            *
!*  Algorithm: If the surface elevation > 0, the bottom layer temperature and   *
!*             water vapor volume mixing ratio are obtained through linear      *
!*             interpolation, while the bottom layer pressure is obtained       *
!*             through exponential interpolation.                               *
!*  Globals used:  H, T, P, VMR, NB, NL, HSURF - these values are adjusted if   *
!*                      HSURF > 0                                               *
!*  Global output:  CLMVAP - Column water vapor amount in unit of cm.           *
!*                       Q - Number of molecules above the surface at one       *
!*                           atmosphere in units of molecules/cm**2             *
!*  Return Codes: none                                                          *
!*  Special Considerations: none                                                *
!*                                                                              *
!********************************************************************************
 */

void model_adjust() {

    float *h = getinput3_.h, *t = getinput3_.t, *p = getinput3_.p, *vmr = getinput3_.vmr;
    int32_t nl = getinput3_.nl, nb = getinput3_.nb;
    float hsurf = getinput5_.hsurf;
    float xppp = getinput14_.xppp;

    static float hp[MODELMAX], tp[MODELMAX], pp[MODELMAX], vmrp[MODELMAX];

    float dhk, dhs, tsurf, vmrs, psurf, amtvrt, damtvrt, hplane;
    float q, tplane, pplane, vmrsp, dhss, dhkk, amtvrtp, damtvrtp;
    int32_t kk=-1, i, k_surf;

    q = 2.152e25;
    model_adj1_.q = q;

    //Determine index of H() such that H(index) < HSURF < H(index+1)

    for (i = 0; i < nb; i++)
        if (hsurf == h[i]) hsurf = h[i] + 0.0001;

    //       Determine index of H() such that H(index) < HSURF < H(index+1)
    locate_pos_(h, &nb, &hsurf, &k_surf);

    //    K_SURF is an index relative to the original model atmosphere (0 - 100 km)
    if (k_surf < 0) {
        fprintf(stderr, "***WARNING: Surface elevation smaller then lowest boundary of the model atmosphere.\n");
        k_surf = 0;
    } else {
        dhk = h[k_surf + 1] - h[k_surf];
        dhs = hsurf - h[k_surf];
        //linear interpolation for surface temperature (TSURF) and VMR (  VMRS)
        tsurf = t[k_surf]+(dhs / dhk)*(t[k_surf + 1] - t[k_surf]);
        vmrs = vmr[k_surf]+(dhs / dhk)*(vmr[k_surf + 1]);
        //exponential interpolation for surface pressure (PSURF)
        psurf = p[k_surf] * exp(-log(p[k_surf]));
        h[0] = hsurf;
        p[0] = psurf;
        t[0] = tsurf;
        vmr[0] = vmrs;

        nb -= k_surf;
        nl = nb - 1;

        for (i = 1; i < nb; i++) {
            h[i] = h[i + k_surf];
            p[i] = p[i + k_surf];
            t[i] = t[i + k_surf];
            vmr[i] = vmr[i + k_surf];
        }
        //Zero out pressures and VMRS of top atmospheric layers.
        for (i = nb; i < MODELMAX; i++) {
            h[i] = 1000;
            p[i] = 0.0;
            t[i] = 300.;
            vmr[i] = 0.0;
        }
    }
    amtvrt = 0.0;

    for (i = 0; i < nl; i++) {
        damtvrt = (q) *(p[i] - p[i + 1])*(vmr[i] + vmr[i + 1]) / 2.0;
        amtvrt += damtvrt;
    }

    model_adj1_.clmvap = amtvrt / 3.34e22;

    printf("ATREM: Model_adjust: Column vapor amount in model atmosphere from ground to space = %f cm\n", model_adj1_.clmvap);
    /*
         Setting the upward atmospheric path's T, P, and VMR profiles:

          1st duplicate the entire atmospheric profiles from the downward path
              to the upward path

     */

    for (i = 0; i < MODELMAX; i++) {
        hp[i] = h[i];
        tp[i] = t[i];
        pp[i] = p[i];
        vmrp[i] = vmr[i];
    }
    //   Set the highest plane altitude to the upper bound of model atmosphere
    hplane = xppp;

    if (hplane >= 100.0) hplane = 100.0 - 0.0001;
    //Do special processing if the plane height (HPLANE) is greater than HP(1)
    if (hplane > hp[0]) {
        // Reset Plane altitude HPLANE (= XPPP) to a larger value if
        // HPLANE.EQ.HP(I) to avoid possible problems in table
        // searching using LOCATE
        for (i = 0; i < MODELMAX; i++)
            if (hplane == hp[i]) hplane = hp[i] - 0.0001;
        //Determine index of HP() such that HP(index) < HPLANE < H(index+1)
        locate_pos_(hp, &nb, &hplane, &kk);

        if (kk < 0)
            printf("WARNING: Plane altitude less then lowest boundary of the model atmosphere.\n");
        else {
            dhkk = hp[kk + 1] - hp[kk];
            dhss = hplane - hp[kk];
            //linear interpolation for plane temperature (TPLANE) and VMR (  VMRSP)
            tplane = tp[kk] + (dhss / dhkk)*(tp[kk + 1] - tp[kk]);
            vmrsp = vmrp[kk] + (dhss / dhkk)*(vmrp[kk + 1] - vmrp[kk]);
            //exponential interpolation for plane pressure (PPLANE)
            pplane = pp[kk] * exp(-log(pp[kk] / pp[kk + 1]) * dhss / dhkk);
            hp[kk + 1] = hplane;
            pp[kk + 1] = pplane;
            tp[kk + 1] = tplane;
            vmrp[kk + 1] = vmrsp;
            //Zero out pressures and VMRP of top atmospheric layers
            if (kk < MODELMAX - 2) {
                for (i = kk + 2; i < MODELMAX; i++) {
                    hp[i] = 1000;
                    pp[i] = 0.0;
                    tp[i] = 300.;
                    vmrp[i] = 0.0;
                }
            }
        }
    }
    amtvrtp = 0.0;
    for (i = 0; i < kk; i++) {
        damtvrtp = (q) *(pp[i] - pp[i + 1])*(vmrp[i] + vmrp[i + 1]) / 2.0;
        amtvrtp += damtvrtp;
    }
    model_adj3_.clmvapp = amtvrtp / 3.34e22;
    printf("ATREM: Model_adjust: Column vapor amount below plane = %f cm\n", model_adj3_.clmvapp);

    //Indices and parameters for the plane layer
    //model_adj3_.k_plane = kk;model_adj4_.k_surf = k_surf; // For final C version
    //printf("ATREM: MODEL_ADJUST - be sure to change next line to model_adj3_.k_plane = kk;model_adj4_.k_surf = k_surf once final conversion to C is complete\n ");
    model_adj3_.k_plane = kk + 1;
    model_adj4_.k_surf = k_surf + 1;

    model_adj3_.dvap_plane = q * (pp[kk] - pp[kk + 1]) * (vmrp[kk] + vmrp[kk + 1]) / 2.0 / 3.34e22;
    model_adj3_.dvap_layer = q * (p[kk] - p[kk + 1]) * (vmr[kk] + vmr[kk + 1]) / 2.0 / 3.34e22;
    model_adj3_.dp_plane = pp[kk] - pp[kk + 1];
    model_adj3_.dp_layer = p[kk] - p[kk + 1];

}

/*
 * ********************************************************************************
!*                                                 *
!*  Name: LOCATE                                       *
!*  Purpose: given an array XX of length N, and given a value X, returns a value*
!*           J such that X is between XX(J) and XX(J+1).  XX must be monotonic, *
!*  Parameters: XX - monotonic array of values                     *
!*              N  - number of elements in XX                      *
!*              X  - value that will be matched to the XX array                *
!*              J  - index into the XX array where XX(J) <= X <= XX(J+1)           *
!*  Algorithm:  bisectional table searching, copied from Numerical Recipes.     *
!*  Globals used: none                                 *
!*  Global output: none                                *
!*  Return Codes: J=0 or J=N is returned to indicate that X is out of range     *
!*  Special Considerations: none                               *
!*                                         *
!********************************************************************************
 */
void locate_pos_(float *xx, int32_t *n1, float *x1, int32_t *jj) {
    int32_t j, ju, jl, jm, n = *n1;
    float x = *x1;

    jl = -1;
    ju = n;
    while ((ju - jl) > 1) {
        jm = (ju + jl) / 2;
        if ((xx[n - 1] > xx[0]) && (x > xx[jm]))
            jl = jm;
        else
            ju = jm;
    }
    if (x == xx[0])
        j = -1;
    else if (x == xx[n - 1])
        j = n - 2;
    else
        j = jl;

    *jj = j;
}

/* Converted from Fortran by R. Healy 3/2016
 *
!********************************************************************************
!*                                                                              *
!*  Name: GEOMETRY                                                              *
!*  Purpose: Calculates the solar and the observational geometric factors.      *
!*  Parameters: none                                                            *
!*  Algorithm: The solar geometry was obtained based on the latitude, longitude,*
!*             GMT time using programs written by W. Mankin at National Center  *
!*             for Atmospheric Research in Boulder, Colorado. The               *
!*             geometric factors for CO2, O3, N2O, CO, CH4, and O2 are based    *
!*             only on the solar and observational angles. Sixty artificial     *
!*             geometric factors for H2O are set up to produce a transmittance  *
!*             table for different atmospheric water vapor amounts.             *
!*  Globals used:  VRTO3    - Column O3 amount in units of atm-cm               *
!*      IMN,IDY,IYR,IH,IM,IS - time and date of data measurements               *
!*      XLATD,XLATM,XLATS,LATHEM    - Latitude of measured area                 *
!*      XLONGD,XLONGM,XLONGS,LNGHEM - Longitude of measured area                *
!*      CLMVAP - Column water vapor in unit of cm in the model atmosphere       *
!*  Global output:                                                              *
!*      SOLZNI,SOLAZ,OBSZNI,OBSPHI,IDAY - Solar zenith angle, solar azimuth     *
!*            angle, observational zenith angle, observational azimuth angle,   *
!*            and the day number in the year                                    *
!*      GCO2,GO3,GN2O,GCO,GCH4,GO2,SSH2O,GGEOM,TOTLO3 - Geometric factors for   *
!*            different gases, and total O3 amount in the sun-surface ray path. *
!*            The geometric factor is defined as: if the vertical column amount *
!*            of the gas is equal 1, then GGAS is equal to the total amount of  *
!*            the gas in the combined Sun-surface-sensor ray path.              *
!*  Return Codes: none                                                          *
!*  Special Considerations: none                                                *
!*                                                                              *
!********************************************************************************
 */
void geometry() {
    float vap_slant[NH2OMAX];
    int md[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334}, leap_year;
    float ggeom;
    /*
        The VAPVRT array contains 60 column water vapor values used for generating
         a table of transmittance spectra with different amount of water vapor.
         The values in VAPVRT is designed so that there is approximately 2% change
         in the .94-um H2O channel ratio for column water vapor in the range .4-6 cm.
     */
    float vapvrt[NH2OMAX] = {.00, .02, .06, .11, .16, .21, .26, .31, .36, .40,
        .43, .46, .50, .54, .58, .62, .66, .70, .75, .80,
        .86, .92, .98, 1.06, 1.14, 1.22, 1.3, 1.4, 1.5, 1.6,
        1.7, 1.8, 1.9, 2.05, 2.2, 2.35, 2.55, 2.75, 2.95, 3.2,
        3.5, 3.8, 4.1, 4.4, 4.7, 5.0, 5.3, 5.6, 6.0, 6.4,
        7.0, 7.7, 8.5, 9.4, 10.4, 11.6, 13.0, 15.0, 25.0, 50.};
    float hplane = getinput14_.xppp;
    float obszni = geometry_l2gen_.senzn_l2, obsphi = geometry_l2gen_.senaz_l2, solzni = geometry_l2gen_.solzn_l2;
    float clmvap = model_adj1_.clmvap, clmvapp = model_adj3_.clmvapp;
    float *ssh2o = geometry2_.ssh2o;
    float dvap_plane = model_adj3_.dvap_plane, dp_plane = model_adj3_.dp_plane, dvap_layer = model_adj3_.dvap_layer, dp_layer = model_adj3_.dp_layer;
    float mu, mu0;
    int i, k_plane = model_adj3_.k_plane, splitpaths = geometry_l2gen_.splitpaths;
    float vap_sol, vap_sen;
    static int firstCall = 1;

    /*

     VAP_SLANT is a new array for containing two-way total vapor amounts
         Here the number "2" can be changed to other numbers, e.g., 2.5,
         without any major effects in derived water vapor values from
         imaging spectrometer data.

     */

    for (i = 0; i < NH2OMAX; i++)
        vap_slant[i] = vapvrt[i] * (float) 2.0;

    obszni /= RAD_DEG;
    obsphi /= RAD_DEG;
    solzni /= RAD_DEG;

    mu0 = (float) 1. / cos(solzni);
    mu = (float) 1. / cos(obszni);
    ggeom = mu0 + mu;

    geometry5_.mu0 = mu0;
    geometry5_.mu = mu;
    geometry2_.ggeom = ggeom;
    geometry1_.obszni = obszni;
    geometry1_.obsphi = obsphi;
    geometry1_.solzni = solzni;

    if (want_verbose || firstCall) printf("ATREM: GGEOM =%f OBSZNI = %f  OBSPHI = %f solzni=%f degrees :: MU0=%f, MU = %f\n",
            ggeom, obszni, obsphi, solzni, mu0, mu);
    geometry2_.gco2 = ggeom;

    geometry2_.go3 = ggeom;
    if (hplane < 27.) geometry2_.go3 = ggeom - (float) 1. / cos(obszni);

    geometry2_.gn2o = ggeom;
    geometry2_.gco = ggeom;
    geometry2_.gch4 = ggeom;
    geometry2_.go2 = ggeom;

    geometry2_.totlo3 = geometry2_.go3 * getinput3_.vrto3;

    if (want_verbose || firstCall) printf("ATREM: TOTLO3 = %f %f\n", geometry2_.totlo3, getinput3_.vrto3);

    /*
     *  Initialize newly created geometrical factors for each atmospheric
       layers (here G_VAP and G_OTHER are true geometrical factors that
       account for the actual Sun-surface-plane path lengths in the
       model atmosphere):
     *
     */
    //printf("ATREM: GEOMETRY - be sure to change line to k_plane NOT k_plane-1 once final conversion to C is complete\n ");
    //---For layers below the plane layer---
    //for (i=0;i<k_plane;i++) {
    for (i = 0; i < k_plane - 1; i++) {
        geometry3_.g_vap[i] = ggeom;
        geometry3_.g_other[i] = ggeom;
    }

    //printf("ATREM: GEOMETRY - be sure to change line to k_plane+1 NOT k_plane once final conversion to C is complete\n ");
    /*---For layers above the plane layer */
    //     for (i=k_plane+1;i<MODELMAX;i++) {
    for (i = k_plane; i < MODELMAX; i++) {
        geometry3_.g_vap[i] = ggeom - (float) 1. / cos(obszni);
        geometry3_.g_other[i] = ggeom - (float) 1. / cos(obszni);
    }
    /*      ---Special treatment for the plane layer to take account the
               "shorter" upward path length
     */
    //printf("ATREM: GEOMETRY - be sure to change line to k_plane NOT k_plane-1 once final conversion to C is complete\n ");
    /*
        geometry3_.g_vap[k_plane]   = ggeom - 1./cos(obszni)  \
                                  + dvap_plane/dvap_layer/cos(obszni);
        geometry3_.g_other[k_plane] = ggeom - 1./cos(obszni)  \
                                  + dp_plane/dp_layer/cos(obszni);
     */
    geometry3_.g_vap[k_plane - 1] = ggeom - 1. / cos(obszni)
            + dvap_plane / dvap_layer / cos(obszni);
    geometry3_.g_other[k_plane - 1] = ggeom - 1. / cos(obszni)
            + dp_plane / dp_layer / cos(obszni);
    /*
         Calculate the water vapor SCALING factor relative to the total amount
            of water vapor in the model atmosphere in the L-shaped
            Sun-surface-plane ray path.
     */
    geometry4_.vap_slant_mdl = clmvap / cos(solzni) + clmvapp / cos(obszni);
    vap_sol = clmvapp*mu0;
    vap_sen = clmvapp*mu;
    if (want_verbose || firstCall) printf("ATREM: VAP_SLANT_MDL = %f cm\n", geometry4_.vap_slant_mdl);

    /*
         The "equivalent" geometrical factor corresponding to the total
            slant vapor amount of VAP_SLANT_MDL':
     */

    geometry3_.g_vap_equiv = geometry4_.vap_slant_mdl / clmvap;

    if (want_verbose || firstCall)
        printf("ATREM: G_VAP_EQUIV = %f clmvap=%f VAP_SOL=%f,VAP_SEN=%f \n", geometry3_.g_vap_equiv, clmvap,
            vap_sol, vap_sen);

    for (i = 0; i < NH2OMAX; i++) {
        ssh2o[i] = vap_slant[i] / geometry4_.vap_slant_mdl;
        //         if (geometry_l2gen_.water_vapor > 0) {
        if (splitpaths) {
            geometry5_.ssh2o_s[0][i] = vapvrt[i] / vap_sol;
            geometry5_.ssh2o_s[1][i] = vapvrt[i] / vap_sen;
            //             printf("SSH2O_S(1): %d %g %g | %g | %g \n",i,geometry5_.ssh2o_s[0][i],vapvrt[i],vap_sol,solzni*RAD_DEG);
            //             printf("SSH2O_S(2): %d %g %g | %g | %g \n",i,geometry5_.ssh2o_s[1][i],vapvrt[i],vap_sen,obszni*RAD_DEG);
        }
        //         } else if (splitpaths != 0) {
        //             printf("ATREM: Split paths is not working because WaterVapor is 0\n");
        //         }
    }
    /*
         Calculate the number of days that have passed in this year.  Take leap year
         into account.
     */
    geometry1_.day = md[getinput8_.imn] + getinput8_.idy;
    leap_year = getinput8_.iyr - (4 * (getinput8_.iyr / 4));
    if ((!leap_year) && (geometry1_.day > 59) && (getinput8_.imn != 2)) geometry1_.day++;
    firstCall = 0;
}

void get_abscf_data(int levels, int bands, int sds_id, char filename[FILENAME_MAX], float* abscf, char *varname) {
    size_t start[3], count[3];
    int retval;
    int xid;

    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    count[0] = levels;
    count[1] = bands;
    count[2] = 0;
    retval = nc_inq_varid(sds_id, varname, &xid);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, filename, varname);
        exit(FATAL_ERROR);
    }
    retval = nc_get_vara_float(sds_id, xid, start, count, abscf);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, filename, varname);
        exit(FATAL_ERROR);
    }

    return;
}

void apply_gas_trans(int32_t k_surf, int levels, int bands, float* abscf, float* dp, float* vmrm, float *tran_hi) {
    static float *sumcf;
    int32_t i, j;
    float Q = 2.152e25;

    if ((sumcf = (float *) calloc(bands, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to sumcf\n");
        exit(FATAL_ERROR);
    }

    //printf("ATREM: INIT_SPECCAL - be sure to change below line to i=k_surf and i - k_surf (NOT i - k_surf + 1) once final conversion to C is complete\n ");

    for (i = 0; i < bands; i++)
        sumcf[i] = 0.;

    for (i = k_surf - 1; i < levels; i++)
        for (j = 0; j < bands; j++)
            sumcf[j] -= abscf[bands * i + j] * dp[i - k_surf + 1] * vmrm[i - k_surf + 1];

    for (i = 0; i < bands; i++)
        tran_hi[i] *= exp(sumcf[i] * Q * (float) 28.966 / (float) 6.0225E+23 / (float) 1.0E-06);

    free(sumcf);
}

/* Ported from Fortran code by R. Healy 3/17/2016
 * Comments from Fortran Code - unless noted otherwise
!********************************************************************************
!*                                         *
!*  Name: INIT_SPECCAL                                 *
!*  Purpose: initialize global data for spectrum calculations.             *
!*  Parameters: none                                   *
!*  Algorithm: initialize data.                            *
!*  Globals used: AH2O, APH2O, BH2O, BPH2O, SODLT, SOGAM, O3CF - Band model     *
!*                             parameters for spectral calculations.            *
!*                WNDOW1, WNDOW2, WP94C, WNDOW3, WNDOW4, W1P14C - Center        *
!*                             positions of window and water vapor absorption   *
!*                             channels used in 3-channel ratio calculations.   *
!*                NB1, NB2, NBP94, NB3, NB4, NB1P14 - Number of narrow channels *
!*                             to form broader window and absorption channels.  *
!*  Global output:                                 *
!*     IH2OLQ,RLQAMT,NGASTT,NH2O,VSTART,VEND - Flag for including liquid           *
!*                     water, liquid water amount (cm), total number of gases   *
!*                     (typically 8), number of water vapor values, starting    *
!*                     and ending wavelengths in internal calculations.         *
!*     NO3PT,NCV,NCVHAF,NCVTOT,VMIN,ISTART,IEND  - Number of O3 abs. coef.     *
!*                     points, parameters for gaussian function and spectral    *
!*                     calculations.                           *
!*     ISTCAL,IEDCAL,DP,PM,TM,VMRM - Parameters for spectral calculations       *
!*     IST1,IED1,IST2,IED2,ISTP94,IEDP94 - 3-channel ratioing parameters    for    *
!*                     the 0.94-um water vapor band.                   *
!*     IST3,IED3,IST4,IED4,IST1P14,IED1P14 - 3-channel ratioing parameters for  *
!*                     the 1.14-um water vapor band.                   *
!*     WT1,WT2,WT3,WT4,JA - Relative weights for the four window channels       *
!*                     used in channel-ratioing calculations. JA is a          *
!*                     output parameter from a table searching routine.        *
!*     NCV2,NCVHF2,NCVTT2,ISTRT2,IEND2,FINST2 - Parameters for smoothing           *
!*                     output reflectance spectra.                         *
!*     NATOT,NBTOT,NCTOT,NDTOT - Number of channels for the four AVIRIS'           *
!*                     grating spectrometers (A, B, C, and D).             *
!*  Return Codes: None.                                *
!*  Special Considerations:  Some parameters may need to be fine-tuned.        *
!*                                         *
!********************************************************************************
!*  Notes about water vapor VMRS and related quantities:                   *
!*                                         *
!*     VAPVRT(60)   - a table containing 60 column vapor values (in unit of cm) *
!*                                         *
!*     VAP_SLANT(I) = VAPVRT(I) * 2.0, VAP_SLANT is a new table for containing  *
!*                    two-way total vapor amounts. Here the number "2" can be   *
!*                    changed to other numbers, e.g., 2.5, without major        *
!*                    effects on retrieved water vapor values.                  *
!*                                         *
!*     G_VAP(I = 1,..., NL) = true vapor geometric factor for each layer in     *
!*                    the model atmosphere (after adjusting for the elevated    *
!*                    surface.                                                  *
!*                                         *
!*     VMRM(I) = VMRM(I)*G_VAP(I). The VMRS are multiplied by the geometrical   *
!*                    factor. We can calculate the vapor transmittance on the   *
!*                    Sun-surface-sensor path by assuming a vertical path in    *
!*                    the model atmosphere with geometric-factor-adjusted VMRS. *
!*                                         *
!*     CLMVAP  = vertical column amount from ground to space in model atmosphere*
!*     CLMVAPP = vertical column amount from ground to aircraft or satellite    *
!*                    sensor in model atmosphere                                *
!*     Q       = 2.152E25 = # of molecules above the surface at one atmosphere  *
!*                    (in unit of  molecules/cm**2)                *
!*                                         *
!*     VAP_SLANT_MDL= CLMVAP/COS(SOLZNI) + CLMVAPP/COS(OBSZNI) = total amount   *
!*                    of water vapor in the model atmosphere in the L-shaped    *
!*                    Sun-surface-plane ray path.                               *
!*                                         *
!*     G_VAP_EQUIV  = VAP_SLANT_MDL / CLMVAP = the "equivalent" geometrical     *
!*                    factor corresponding to the total slant vapor amount      *
!*                    VAP_SLANT_MDL and the column vapor amount CLMVAP.         *
!*                                         *
!*     SSH2O(I) (I = 1, ..., 60) - a pure scaling factor relative to the total  *
!*                    slant vapor amount of VAP_SLANT_MDL, and                  *
!*           SSH2O(I) = VAP_SLANT(I) / VAP_SLANT_MDL               *
!*                                         *
!*     SH2O = one value of SSH2O(I). SH2O is used during generation of the      *
!*           look-up table.                            *
!*                                         *
!*     VAPTT  = VAP_SLANT_MDL*SH2O, is the absolute total vapor amount on the   *
!*                    L-shaped path corresponding to a spectrum stored in the   *
!*                    look-up table.                                       *
!*                                         *
!*     CLMWVP = 0.5*(VAPTTA+VAPTTB)/G_VAP_EQUIV, is the retrieved column water  *
!*                    vapor amount from imaging spectrometer data.         *
!********************************************************************************
 */
void init_spectral_calculations() {
    int32_t co2 = getinput1_.co2, o2 = getinput1_.o2, n2o = getinput1_.n2o, co = getinput1_.co, ch4 = getinput1_.ch4, no2 = getinput1_.no2, o3 = getinput1_.o3;
    float *fwhm = getinput4_.fwhm, *wavobs = getinput4_.wavobs, *dp = init_speccal5_.dp, *tm = init_speccal5_.tm, *pm = init_speccal5_.pm, *vmrm = init_speccal5_.vmrm;
    float *finst2 = init_speccal10_.finst2;
    int32_t nbands = getinput5_.nbands;
    ;
    float totlo3 = geometry2_.totlo3; //=no2cf_init1_.rno2cf;
    float dlt2 = getinput5_.dlt2, wavcv2, dwvavr, sumins;
    float sno2 = getinput3_.sno2, totno2;
    float go3 = geometry2_.go3, Q = model_adj1_.q;
    static float *abscf_co2, *abscf_o2, *abscf_n2o, *abscf_ch4, *abscf_co, const1, *o3cf, *rno2cf; //=o3cf_init1_.o3cf;
    static int firstCall = 1;
    int32_t i, j, k_surf = model_adj4_.k_surf, nl = getinput3_.nl;
    int32_t iwndw1, iwndw2, iwndw4, iwndw5, iwp94c, iw1p14c, nchnla, nchnlb, nchnlc, nchnld, nb1haf, nb2haf, nb3haf, nb4haf;
    char *abscf_file = "abscf_gas.nc";
    char *abscf_no2o3file = "abscf_no2o3.nc";
    static size_t ncf;

    char *filedir;
    char filename[FILENAME_MAX];
    size_t start[3], count[3];
    int xid, yid;
    int retval, sds_id;
    static size_t bands, levels;

    float vrtno2 = 5.0E+15, gno2, sclgas;
    float *sumcf;

    if (firstCall) {
        const1 = 4.0 * log(2); // gaussian curve area constant
        init_speccal3_.nh2o = NH2OMAX; // number of water vapor values in table
        //Wavelength of medium resolution spectrum FWHM=.2 nm, .56-3.1 um

        for (i = 0; i < NP_MED; i++) {
            init_speccal12_.wavln_med[i] = (float) VSTART + i * (float) DWAVLN;
            /*
             Note: The grids of WAVNO_HI do not match the grids of 10000./WAVLN_MED.
                   INDEX_MED is a specially designed index for finding close matches
                   between the two kinds of grids.
             */
            init_speccal13_.index_med[i] = ((float) 10000. / init_speccal12_.wavln_med[i] - (float) 3000.) / (float) DWAVNO + 1;
            /*
             Note:     WAVLN_MED_INDEX(I) is very close to WAVLN_MED(I),
                   and WAVLN_MED_INDEX(I) >= WAVLN_MED(I)


             */
            init_speccal13_.wavln_med_index[i] = (float) 10000. / ((init_speccal13_.index_med[i] - (float) 1.)*(float) DWAVNO + (float) 3000.);
        }
        //Wavelength of medium resolution spectrum FWHM=.2 nm, .3-3.1 um
        for (i = 0; i < NP_STD; i++)
            init_speccal12_.wavln_std[i] = (float) 0.3 + i * (float) DWAVLN;
        /*
         Initialize arrays for smoothing medium resolution spectrum (DLT_MED = 0.2 nm,
                 and point spacing DWAVLN = 0.0001 micron) to coarser spectral
                 resolution data from imaging spectrometers.

                          NCVHF(I) = ( FACDLT * FWHM(I) / DWAVLN + 1.)

         */
        for (i = 0; i < nbands; i++)
            init_speccal15_.ncvhf[i] = (float) FACDLT * fwhm[i] / (float) DWAVLN + (float) 1.0;

        // parameters and arrays to smooth output surface reflectance spectrum
        wavcv2 = FACDLT*dlt2;
        /*
         * Find the largest value in the FWHM array, and use this value in calculation
           of indices for smoothing output reflectance spectra. This smoothing
           algorithm should work well with grating spectrometers having nearly
           constant spectral resolutions, but not so well for prism spectrometers
           having variable spectral resolution.
         */
        dwvavr = fwhm[0];

        for (i = 1; i < nbands; i++)
            if (dwvavr < fwhm[i]) dwvavr = fwhm[i];

        init_speccal10_.ncv2 = wavcv2 / dwvavr;
        init_speccal10_.ncvhf2 = init_speccal10_.ncv2 + 1;
        init_speccal10_.ncvtt2 = init_speccal10_.ncv2 + 1;

        if (dlt2 != 0) {
            sumins = 0;
            for (i = init_speccal10_.ncvhf2 - 1; i < init_speccal10_.ncvtt2; i++) {
                finst2[i] = exp(-const1 * pow((i - init_speccal10_.ncvhf2 + 1) * dwvavr / dlt2, 2));
                sumins += finst2[i];
            }
            for (i = 0; i < init_speccal10_.ncvhf2 - 1; i++) {
                finst2[i] = finst2[init_speccal10_.ncvtt2 - i - 1];
                sumins += finst2[i];
            }
            sumins *= dwvavr;

            for (i = 0; i < init_speccal10_.ncvtt2; i++) {
                finst2[i] *= (dwvavr / sumins);
            }
        }

        init_speccal10_.istrt2 = init_speccal10_.ncvhf2;
        init_speccal10_.iend2 = nbands - init_speccal10_.ncvhf2;

        /*  number of channels of the four AVIRIS spectrometers.  These are used
            in removing null AVIRIS radiance values in the overlap portions of two
            adjacent spectrometers.

            Note by R.Healy - The Aviris correction isn't being used in l2gen yet. (3/2016)
            There are no sensor specific modifications being applied.  Atrem has only been
            tested on HICO data.
         */
        nchnla = 32;
        nchnlb = 64;
        nchnlc = 64;
        nchnld = 64;

        init_speccal11_.natot = nchnla;
        init_speccal11_.nbtot = nchnla + nchnlb;
        init_speccal11_.nctot = nchnla + nchnlb + nchnlc;
        init_speccal11_.ndtot = nchnla + nchnlb + nchnlc + nchnld;

        /* Resetting window wavelength positions and calculating weights for
            window and absorption channels used in 3-channel ratioing.
            Note that the C version of this function returns 0 based indices
         */
        iwndw1 = findMatch(wavobs, nbands, getinput6_.wndow1);
        iwndw2 = findMatch(wavobs, nbands, getinput6_.wndow2);

        getinput6_.wndow1 = wavobs[iwndw1];
        getinput6_.wndow2 = wavobs[iwndw2];

        if ((getinput7_.nb1 % 2) == 0) getinput7_.nb1 += 1;
        if ((getinput7_.nb2 % 2) == 0) getinput7_.nb2 += 1;

        nb1haf = (getinput7_.nb1 - 1) / 2;
        nb2haf = (getinput7_.nb2 - 1) / 2;

        // Adding one to all the indexes because of legacy fortran - subtracted in get_atrem_cor
        init_speccal6_.ist1 = iwndw1 - nb1haf + 1;
        init_speccal6_.ied1 = iwndw1 + nb1haf + 1;
        init_speccal6_.ist2 = iwndw2 - nb2haf + 1;
        init_speccal6_.ied2 = iwndw2 + nb2haf + 1;

        iwp94c = findMatch(wavobs, nbands, getinput6_.wp94c);
        getinput6_.wp94c = wavobs[iwp94c];

        if ((getinput7_.nbp94 % 2) == 0) getinput7_.nbp94 += 1;

        nb3haf = (getinput7_.nbp94 - 1) / 2;
        init_speccal6_.istp94 = iwp94c - nb3haf + 1;
        init_speccal6_.iedp94 = iwp94c + nb3haf + 1;

        init_speccal8_.wt1 = (getinput6_.wndow2 - getinput6_.wp94c) / (getinput6_.wndow2 - getinput6_.wndow1);
        init_speccal8_.wt2 = (getinput6_.wp94c - getinput6_.wndow1) / (getinput6_.wndow2 - getinput6_.wndow1);

        iwndw4 = findMatch(wavobs, nbands, getinput6_.wndow3);
        iwndw5 = findMatch(wavobs, nbands, getinput6_.wndow4);

        getinput6_.wndow3 = wavobs[iwndw4];
        getinput6_.wndow4 = wavobs[iwndw5];

        if ((getinput7_.nb3 % 2) == 0) getinput7_.nb3 += 1;
        if ((getinput7_.nb4 % 2) == 0) getinput7_.nb4 += 1;

        nb3haf = (getinput7_.nb3 - 1) / 2;
        nb4haf = (getinput7_.nb4 - 1) / 2;

        init_speccal7_.ist3 = iwndw4 - nb3haf + 1;
        init_speccal7_.ied3 = iwndw4 + nb3haf + 1;
        init_speccal7_.ist4 = iwndw5 - nb4haf + 1;
        init_speccal7_.ied4 = iwndw5 + nb4haf + 1;

        iw1p14c = findMatch(wavobs, nbands, getinput6_.w1p14c);
        getinput6_.w1p14c = wavobs[iw1p14c];

        if ((getinput7_.nb1p14 % 2) == 0) getinput7_.nb1p14 += 1;

        nb3haf = (getinput7_.nb1p14 - 1) / 2;
        init_speccal7_.ist1p14 = iw1p14c - nb3haf + 1;
        init_speccal7_.ied1p14 = iw1p14c + nb3haf + 1;

        init_speccal8_.wt3 = (getinput6_.wndow4 - getinput6_.w1p14c) / (getinput6_.wndow4 - getinput6_.wndow3);
        init_speccal8_.wt4 = (getinput6_.w1p14c - getinput6_.wndow3) / (getinput6_.wndow4 - getinput6_.wndow3);

        // Open the netcdf4 input file
        if ((filedir = getenv("OCDATAROOT")) == NULL) {
            printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
            exit(FATAL_ERROR);
        }

        strcpy(filename, filedir);
        strcat(filename, "/common/");
        strcat(filename, abscf_file);
        retval = nc_open(filename, NC_NOWRITE, &sds_id);
        if (retval != NC_NOERR) {
            fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
                    __FILE__, __LINE__, filename);
            exit(FATAL_ERROR);
        }
        // Now read in the absorption coefficients
        retval = nc_inq_dimid(sds_id, "level", &xid)
                || nc_inq_dimlen(sds_id, xid, &levels)
                || nc_inq_dimid(sds_id, "band", &yid)
                || nc_inq_dimlen(sds_id, yid, &bands);

        if (retval) {
            fprintf(stderr, "-E- %s line %d: nc_inq_dimid(%s) failed.%d %d %d \n",
                    __FILE__, __LINE__, filename, (int) levels, (int) bands, retval);
            exit(FATAL_ERROR);
        }

        abscf_o2 = (float *) calloc(levels*bands, sizeof (float));
        if (!(abscf_o2)) {
            fprintf(stderr, "-E- %s line %d: Failed to allocate memory for abscf_n2o.\n",
                    __FILE__, __LINE__);
            exit(FATAL_ERROR);
        }

        abscf_co2 = (float *) calloc(levels*bands, sizeof (float));
        if (!(abscf_co2)) {
            fprintf(stderr, "-E- %s line %d: Failed to allocate memory for abscf_co2.\n",
                    __FILE__, __LINE__);
            exit(FATAL_ERROR);
        }
        abscf_co = (float *) calloc(levels*bands, sizeof (float));
        if (!(abscf_co)) {
            fprintf(stderr, "-E- %s line %d: Failed to allocate memory for abscf_co.\n",
                    __FILE__, __LINE__);
            exit(FATAL_ERROR);
        }
        abscf_ch4 = (float *) calloc(levels*bands, sizeof (float));
        if (!(abscf_ch4)) {
            fprintf(stderr, "-E- %s line %d: Failed to allocate memory for abscf_ch4.\n",
                    __FILE__, __LINE__);
            exit(FATAL_ERROR);
        }
        abscf_n2o = (float *) calloc(levels*bands, sizeof (float));
        if (!(abscf_n2o)) {
            fprintf(stderr, "-E- %s line %d: Failed to allocate memory for abscf_n2o.\n",
                    __FILE__, __LINE__);
            exit(FATAL_ERROR);
        }

        get_abscf_data(levels, bands, sds_id, filename, abscf_o2, "abscf_o2");
        get_abscf_data(levels, bands, sds_id, filename, abscf_co2, "abscf_co2");
        get_abscf_data(levels, bands, sds_id, filename, abscf_co, "abscf_co");
        get_abscf_data(levels, bands, sds_id, filename, abscf_ch4, "abscf_ch4");
        get_abscf_data(levels, bands, sds_id, filename, abscf_n2o, "abscf_n2o");

        nc_close(sds_id);

        strcpy(filename, filedir);
        strcat(filename, "/common/");
        strcat(filename, abscf_no2o3file);
        strcat(filename, "\0");
        retval = nc_open(filename, NC_NOWRITE, &sds_id);
        if (retval != NC_NOERR) {
            fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
                    __FILE__, __LINE__, filename);
            exit(FATAL_ERROR);
        }
        // Now read in the absorption coefficients
        retval = nc_inq_dimid(sds_id, "ncf", &xid)
                || nc_inq_dimlen(sds_id, xid, &ncf);

        if (retval) {
            fprintf(stderr, "-E- %s line %d: nc_inq_dimid(%s) failed. %d %d \n",
                    __FILE__, __LINE__, filename, (int) ncf, retval);
            exit(FATAL_ERROR);
        }

        o3cf = (float *) calloc(ncf, sizeof (float));
        if (!(o3cf)) {
            fprintf(stderr, "-E- %s line %d: Failed to allocate memory for abscf_o3.\n",
                    __FILE__, __LINE__);
            exit(FATAL_ERROR);
        }
        rno2cf = (float *) calloc(ncf, sizeof (float));
        if (!(o3cf)) {
            fprintf(stderr, "-E- %s line %d: Failed to allocate memory for abscf_no2.\n",
                    __FILE__, __LINE__);
            exit(FATAL_ERROR);
        }

        start[0] = 0;
        start[1] = 0;
        start[2] = 0;
        count[0] = ncf;
        count[1] = 0;
        count[2] = 0;

        retval = nc_inq_varid(sds_id, "o3cf", &xid);

        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, filename, "o3cf");
            exit(FATAL_ERROR);
        }
        retval = nc_get_vara_float(sds_id, xid, start, count, o3cf);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, filename, "o3cf");
            exit(FATAL_ERROR);
        }

        retval = nc_inq_varid(sds_id, "no2cf", &xid);

        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, filename, "no2cf");
            exit(FATAL_ERROR);
        }
        retval = nc_get_vara_float(sds_id, xid, start, count, rno2cf);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, filename, "no2cf");
            exit(FATAL_ERROR);
        }


        nc_close(sds_id);

        firstCall = 0;
    }

    // Initialize the TRAN_HI array for high resolution spectrum:
    for (i = 0; i < bands; i++)
        init_speccal1_.tran_hi_others[i] = 1.0;

    /*
         Calculate medium resolution O3 transmittances (0.2 nm, 2-way) with
         a point spacing of 0.1 nm between 0.3 and 0.8 micron.
     */
    if (o3) {
        for (i = 0; i < NO3PT; i++) {
            init_speccal16_.tran_o3_std[i] = exp(-totlo3 * o3cf[i]);
            //printf("TRAN_O3_STD %d %g %f %f\n",i+1,init_speccal16_.tran_o3_std[i],totlo3,o3cf[i]);
        }
    } else {
        /* If NO2 is not intended to be included in total atmospheric gaseous
           transmittance calculations, assigning TRAN_NO2_STD = 1.0:
         */
        for (i = 0; i < NO3PT; i++)
            init_speccal16_.tran_o3_std[i] = 1.0;

    }
    /*
         Calculate medium resolution NO2 transmittances (0.2 nm, 2-way) with
         a point spacing of 0.1 nm between 0.3 and 0.8 micron.
     */

    vrtno2 = 5.0E+15;
    vrtno2 *= sno2;

    gno2 = go3;
    totno2 = gno2 * vrtno2;

    if (no2) {
        for (i = 0; i < NO2PT; i++) {
            init_speccal17_.tran_no2_std[i] = exp(-totno2 * rno2cf[i]);
            //printf("TRAN_NO2_STD %d %g %g %g\n",i+1,init_speccal17_.tran_no2_std[i],totno2,rno2cf[i]);
        }
    } else {
        /* If NO2 is not intended to be included in total atmospheric gaseous
           transmittance calculations, assigning TRAN_NO2_STD = 1.0:
         */
        for (i = 0; i < NO2PT; i++)
            init_speccal17_.tran_no2_std[i] = 1.0;
    }

    for (i = 0; i < nl; i++) {
        dp[i] = getinput3_.p[i] - getinput3_.p[i + 1];
        pm[i] = (getinput3_.p[i] + getinput3_.p[i + 1]) / (float) 2.0;
        tm[i] = (getinput3_.t[i] + getinput3_.t[i + 1]) / (float) 2.0;
    }

    /*

!C
!C Calculate high resolution transmittances (0.05 cm-1) of CO2, N2O, CO,
!C     CH4, and O2 in the 0.56 - 3.1 micron range, and save values for
!C     calculating total atmospheric transmittances later.
!C     Because water vapor amounts are allowed to vary,
!C     the high resolution water vapor transmittances are calculated
!C     in subroutines TRAN_TABLE and TRANCAL. TRAN_TABLE provides variable
!C     water vapor amounts, and calls TRANCAL for the calculation of
!C     corresponding vapor transmittance spectrum.
!C

     */

    // For CO2 transmittance calculation -
    if (co2) {
        /*
    --     On 2/7/2013 B.-C. Gao made the modification - Increased SCLCO2 i
    --     from 1.0 to 1.1 to reflect the fact that the CO2 VMR reached the
    --     2012 level of 390 ppmv.
    --
         */
        sclgas = 1.1;
        /*
                 Scale the VMRM by the two-way path geometrical factors. The geometric
                           factors, G_OTHER, varies with atmospheric layer number for
                           aircraft observational geometries.
         */
        for (i = 0; i < nl; i++)
            vmrm[i] = sclgas * (float) 355.0 * (float) 1.e-6 * geometry3_.g_other[i];

        apply_gas_trans(k_surf, levels, bands, abscf_co2, dp, vmrm, init_speccal1_.tran_hi_others);
        //        for (i = 0; i < bands; i++)
        //            if (init_speccal1_.tran_hi_others[i]<1.0) printf("CO2: %d tran_hi_others=%f\n",i, init_speccal1_.tran_hi_others[i]);

    }


    /*
        --------------------------------------------
               For N2O transmittance calculation.
     */
    if (n2o) {
        sclgas = 0.3;
        for (i = 0; i < nl; i++)
            vmrm[i] = sclgas * (float) 1.e-6 * geometry3_.g_other[i];

        apply_gas_trans(k_surf, levels, bands, abscf_n2o, dp, vmrm, init_speccal1_.tran_hi_others);
        //        for (i = 0; i < bands; i++)
        //            if (init_speccal1_.tran_hi_others[i]<1.0)printf("N2O: %d tran_hi_others=%f\n",i, init_speccal1_.tran_hi_others[i]);
    }
    /*
       --------------------------------------------
              For CO transmittance calculation.
     */
    if (co) {
        sclgas = 0.1;
        for (i = 0; i < nl; i++)
            vmrm[i] = sclgas * (float) 1.e-6 * geometry3_.g_other[i];

        apply_gas_trans(k_surf, levels, bands, abscf_co, dp, vmrm, init_speccal1_.tran_hi_others);
        //        for (i = 0; i < bands; i++)
        //            if (init_speccal1_.tran_hi_others[i]<1.0)printf("CO: %d tran_hi_others=%f\n",i, init_speccal1_.tran_hi_others[i]);
    }
    /*
--------------------------------------------
       For CH4 transmittance calculation.
       For assigning CH4 VMRM
 NOTE: The scaling factor of 0.8 for the CH4 VMRS was obtained by comparing
       transmittance spectra calculated using our program, which is based on
       the Malkmus narrow band spectral model, with a ratioed spectrum
       provided by G. C. Toon at Jet Propulsion Laboratory (JPL). The JPL
       ratio spectrum was obtained by ratioing a high resolution (0.005 cm-1)
       solar spectrum measured at ground level against a solar spectrum
       measured at an altitude of approximately 35 km with the same Fourier
       Transform spectrometer. The high resolution ratio spectrum was
       degraded to a resolution of 10 nm during our derivation of the
       scaling factor for the CH4 VMRS.

        sclch4=0.8;
     */
    if (ch4) {
        sclgas = 1.0;
        for (i = 0; i < nl; i++)
            vmrm[i] = sclgas * (float) 1.6 * (float) 1.0e-6 * geometry3_.g_other[i]; //I don't know why sclgas=1 and then it's multiplied by 1.6 - rjh

        apply_gas_trans(k_surf, levels, bands, abscf_ch4, dp, vmrm, init_speccal1_.tran_hi_others);
        //        for (i = 0; i < bands; i++)
        //            if (init_speccal1_.tran_hi_others[i]<1.0)printf("CH4: %d tran_hi_others=%f\n",i, init_speccal1_.tran_hi_others[i]);

    }
    if (o2) {
        /*
         * ***Modified by Bo-Cai Gao on 2/7/2013 to increase O2 absorption
          ---  coefficients by the factor SCL_O2 for wavelengths > 1.2 micron
          ---  in order to model properly the atmospheric O2 band centered
          ---  near 1.265 micron.
         *
         */
        sclgas = 0.21;

        if ((sumcf = (float *) calloc(bands, sizeof (float))) == NULL) {
            printf("-E- : Error allocating memory to sumcf\n");
            exit(FATAL_ERROR);
        }

        for (i = 0; i < nl; i++)
            vmrm[i] = sclgas * geometry3_.g_other[i];

        for (i = 0; i < bands; i++)
            sumcf[i] = 0.;
        /*
         * ***Modified by Bo-Cai Gao on 2/7/2013 to increase O2 absorption
            ---  coefficients by the factor SCL_O2 for wavelengths > 1.2 micron
            ---  i& < 1.3333 micron in order to model properly the atmospheric
            ---  O2 band centered near 1.265 micron.
         *
         */
        sclgas = 2.60;
        //printf("ATREM: INIT_SPECCAL - be sure to change below line to i=k_surf and i - k_surf (NOT i - k_surf + 1) once final conversion to C is complete\n ");

        for (i = k_surf - 1; i < levels; i++)
            for (j = 0; j < bands; j++)
                sumcf[j] -= abscf_o2[bands * i + j] * dp[i - k_surf + 1] * vmrm[i - k_surf + 1];

        for (i = 0; i < 9000; i++)
            init_speccal1_.tran_hi_others[i] *= exp(sumcf[i] * Q * (float) 28.966 / (float) 6.0225E+23 / (float) 1.0E-06);
        for (i = 9000; i < 106600; i++)
            init_speccal1_.tran_hi_others[i] *= exp(sumcf[i] * Q * sclgas * (float) 28.966 / (float) 6.0225E+23 / (float) 1.0E-06);
        for (i = 106600; i < NP_HI; i++)
            init_speccal1_.tran_hi_others[i] *= exp(sumcf[i] * Q * (float) 28.966 / (float) 6.0225E+23 / (float) 1.0E-06);
        //        for (i = 0; i < bands; i++) {
        //            printf("OX: %d tran_hi_others %f\n",i, init_speccal1_.tran_hi_others[i]);
        //            if (i==94923) {
        //                printf("OX: %d %f %f\n",i,exp(sumcf[i]*Q*sclgas* (float)28.966 / (float)6.0225E+23 /(float) 1.0E-06),sumcf[i]);
        //                for (k=k_surf-1; k<levels;k++)
        //                    printf("OX: ABSCF: %d %f %f %f\n",k,abscf_o2[bands * k + i],dp[k-k_surf+1],vmrm[k-k_surf+1]);            }
        //        }
        free(sumcf);
    }
}

/*
 * Converted to C 4/2016 - R. Healy
 * Note this is the merger of two subroutines (TRANCAL with TRAN_TABLE) in the original Atrem Fortran code.
 * The original code looped through the 60 Water Vapor values in TRAN_TABLE and called TRANCAL.
!********************************************************************************
!*                                             *
!*  Name: TRAN_TABLE                                   *
!*  Purpose: This subroutine generates a table consisting of 60 atmospheric     *
!*           transmittance spectra at the solar and observational               *
!*           geometry and with 60 column water vapor values. The table also     *
!*           includes the total amounts of column water vapor used in the       *
!*           calculations, and the 3-channel ratios calculated from the window  *
!*           and absorption channels in and around the 0.94- and 1.14-um water  *
!*           vapor bands.                              *
!*  Parameters: none                                   *
!*  Algorithm: For each of the 60 water vapor amounts, calculate the           *
!*             atmospheric transmittance, and save                 *
!*  Globals Used: NH2O   - number of column water vapor values             *
!*       VAPTT  - geometrically adjusted water vapor total         *
!*       R094   - channel ratio for .94 um region              *
!*       R114   - channel ratio for 1.14 um region             *
!*       TRNCAL - atmospheric transmittance spectra            *
!*  Global Output: VAPTOT() - array containing geometrically adjusted water     *
!*               vapor values                      *
!*        ROP94()  - array containing channel ratios for .94 um region *
!*        R1P14()  - array containing channel ratios for 1.14 um region*
!*        TRNTBL() - 2 dimensional array containing one transmittance  *
!*               spectrum for each column water vapor amount       *
!*  Return Codes: none                                                          *
!*  Special Considerations: none                                                *
!*                                         *
!********************************************************************************
!********************************************************************************
!*  TRANCAL combined with TRAN_TABLE by R. Healy 4/28/2015                      *
!*  Name: TRANCAL                                  *
!*  Purpose: This program calculates combined transmittances of H2O, CO2, O3,   *
!*      N2O, CO, CH4, and O2.                                                   *
!*  Parameters: none.                                  *
!*  Algorithm: The calculations were based on the line-by-line absorption       *
!*      parameters supplied by William R. Ridgway of NASA/GSFC.            *
!*  Global output:VAPTT  - geometrically adjusted water vapor total.           *
!*       R094   - channel ratio for 0.94 um region.            *
!*       R114   - channel ratio for 1.14 um region.            *
!*       TRNCAL - total transmittances of all gases that match the     *
!*                         resolutions of imaging spectrometers.               *
!*  Return Codes: none.                                *
!*  Special Considerations: The high resolution (0.05 cm-1) line-by-line        *
!*      absorption parameters cover the 0.555 - 3.33 micron spectral range      *
!*      (3000 - 18000 cm-1). The medium resolution ozone absorption             *
!*      coefficients covers the 0.3-0.8 um spectral range. The line-by-line     *
!*      high resolution spectra were first smoothed to medium resolution        *
!*      spectra (resolution = 0.2 nm, wavelength spacing = 0.1 nm) covering     *
!*      the 0.56 - 3.1 micron spectral region. The medium resolution spectra    *
!*      of O3 and other gases are combined (in SUBROUTINE TRAN_SMOOTH) to form  *
!*      a single medium resolution spectrum from 0.3 to 3.1 micron. This        *
!*      combined spectrum (medium resolution) is then smoothed to lower         *
!*      resolutions to match the resolutions of imaging spectrometers. The      *
!*      smoothing is also done in SUBROUTINE TRAN_SMOOTH.                       *
!*                                         *
!********************************************************************************
 */

void tran_table() {

    int32_t co2 = getinput1_.co2, o2 = getinput1_.o2, n2o = getinput1_.n2o, co = getinput1_.co, ch4 = getinput1_.ch4, no2 = getinput1_.no2, o3 = getinput1_.o3;
    float *wavobs = getinput4_.wavobs, *dp = init_speccal5_.dp, *vmrm = init_speccal5_.vmrm;
    float *vmr = getinput3_.vmr;
    int32_t nbands = getinput5_.nbands, nh2o = init_speccal3_.nh2o, fullcalc = getinput5_.full_calc;
    int32_t nl = getinput3_.nl, splitpaths = geometry_l2gen_.splitpaths, ja = geometry_l2gen_.ja, jb = geometry_l2gen_.jb;
    float Q = model_adj1_.q;
    float mu = geometry5_.mu, mu0 = geometry5_.mu0, *ssh2o = geometry2_.ssh2o, *g_vap = geometry3_.g_vap;
    float dvap_plane = model_adj3_.dvap_plane, dvap_layer = model_adj3_.dvap_layer;
    int32_t k_plane = model_adj3_.k_plane, k_surf = model_adj4_.k_surf;
    float *vaptot = tran_table1_.vaptot;
    int i, j, k, n;
    static int firstCall = 1;

    static float *abscf_h2o;
    static float *tran_hi, *wavno_hi;

    char *filedir;
    char *abscf_file = "abscf_gas.nc";
    char filename[FILENAME_MAX];
    size_t start[3], count[3];
    int xid, yid;
    int retval, sds_id;
    static size_t bands, levels;

    static int32_t n_g = 7, doSmooth;
    float tg[7] = {0.0, 0.379, 0.5106, 0.81, 0.9548, 0.9933, 1.0}, f1, f2;

    static float *tkcdf;
    static float *sum_kd, *sumcf, **sumcf_s;
    float vmrm_s[2][MODELMAX];

    if (firstCall) {

        if (fullcalc)
            printf("ATREM: **WARNING** Full_calc is on.  Processing will be extremely slow...\n");
        else
            printf("ATREM: Full_calc is off.  Processing using k-distribution method...\n");

        doSmooth = co2 || n2o || co || ch4 || o2 || o3 || no2;

        // Open the netcdf4 input file
        if ((filedir = getenv("OCDATAROOT")) == NULL) {
            printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
            exit(FATAL_ERROR);
        }

        strcpy(filename, filedir);
        strcat(filename, "/common/");
        //        strcpy(filename,input_l2gen_.filename);
        strcat(filename, abscf_file);
        strcat(filename, "\0");
        retval = nc_open(filename, NC_NOWRITE, &sds_id);
        if (retval != NC_NOERR) {
            fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
                    __FILE__, __LINE__, filename);
            exit(FATAL_ERROR);
        }
        // Now read in the absorption coefficients
        retval = nc_inq_dimid(sds_id, "level", &xid)
                || nc_inq_dimlen(sds_id, xid, &levels)
                || nc_inq_dimid(sds_id, "band", &yid)
                || nc_inq_dimlen(sds_id, yid, &bands);

        if (retval) {
            fprintf(stderr, "-E- %s line %d: nc_inq_dimid(%s) failed. %d %d %d \n",
                    __FILE__, __LINE__, filename, (int) levels, (int) bands, retval);
            exit(FATAL_ERROR);
        }

        if (!(abscf_h2o = (float *) malloc(levels * bands * sizeof (float)))) {
            printf("-E- : Error allocating memory to abscf_h2o\n");
            exit(FATAL_ERROR);
        }


        start[0] = 0;
        start[1] = 0;
        start[2] = 0;
        count[0] = levels;
        count[1] = bands;
        count[2] = 0;

        retval = nc_inq_varid(sds_id, "abscf_h2o", &xid);

        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, filename, "abscf_h2o");
            exit(FATAL_ERROR);
        }
        retval = nc_get_vara_float(sds_id, xid, start, count, abscf_h2o);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, filename, "abscf_h2o");
            exit(FATAL_ERROR);
        }

        if (!wavno_hi && !(wavno_hi = (float *) malloc(NP_HI * sizeof (float *)))) {
            printf("-E- : Error allocating memory to wavno_hi\n");
            exit(FATAL_ERROR);
        }

        retval = nc_inq_varid(sds_id, "waveno", &xid);

        count[0] = bands;
        count[1] = 0;

        retval = nc_get_vara_float(sds_id, xid, start, count, wavno_hi);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, filename, "wavno_hi");
            exit(FATAL_ERROR);
        }

        nc_close(sds_id);

        if (!(tkcdf = (float *) calloc(n_g * nbands * levels, sizeof (float)))) {
            printf("-E- : Error allocating memory to tkcdf\n");
            exit(FATAL_ERROR);
        }

        if (!(sum_kd = (float *) calloc(nbands, sizeof (float)))) {
            printf("-E- : Error allocating memory to sum_kd\n");
            exit(FATAL_ERROR);
        }

    }

    /*
        ! Moved from INIT_SPECAL because it seems more appropriate here. - rjh 01/04/2016
        !C Initial water vapor VMRs for repeated use in other subroutines and
        !C    adjust layered water vapor VMRM with geometrical factors.
     */
    //if (firstCall) printf("ATREM: TRAN_CAL - be sure to remove below 2 lines (k_surf--;k_plane--) once final conversion to C is complete\n ");
    k_plane--; // RJH - this is a relic of the fortran code - may want change this throughout the code later
    k_surf--; // Ditto

    for (i = 0; i < nl; i++) {
        vmrm[i] = (vmr[i] + vmr[i + 1]) / (float) 2.0;
        /*
        Scale the VMRM by the two-way path geometrical factors. The geometric
        factors, G_VAP, varies with atmospheric layer number for
        aircraft observational geometries.
         */
        if (splitpaths) { // perform the split path calculations
            vmrm_s[0][i] = vmrm[i] * mu0;
            if (i < k_plane)
                vmrm_s[1][i] = vmrm[i] * mu;
            else {
                if (i == k_plane)
                    vmrm_s[1][i] = vmrm[i] * mu * dvap_plane / dvap_layer;
                else
                    vmrm_s[1][i] = 0.;

            }
        }
        vmrm[i] *= g_vap[i];
    }

    if (firstCall || fullcalc) {
        // For each water vapor amount, calculate the geometrically adjusted water
        //     vapor amount, the channel ratios, and the transmittance spectrum.
        if (!sumcf && !(sumcf = (float *) calloc(NP_HI, sizeof (float)))) { //FREE
            printf("-E- : Error allocating memory to sumcf\n");
            exit(FATAL_ERROR);
        }
        if (splitpaths) {
            if (!sumcf_s) {
                if (!(sumcf_s = (float **) malloc(2 * sizeof (float*)))) { //FREE
                    printf("-E- : Error allocating memory to sumcf_s\n");
                    exit(FATAL_ERROR);
                }
                for (i = 0; i < 2; i++) {
                    if (!(sumcf_s[i] = (float *) calloc(NP_HI, sizeof (float)))) { //FREE
                        printf("-E- : Error allocating memory to sumcf_s\n");
                        exit(FATAL_ERROR);
                    }
                }
            }
        }
        for (i = 0; i < bands; i++) {
            sumcf[i] = 0;
            if (splitpaths) {
                sumcf_s[0][i] = 0;
                sumcf_s[1][i] = 0;
            }
            for (j = k_surf; j < levels; j++) {
                sumcf[i] -= abscf_h2o[bands * j + i] * dp[j - k_surf] * vmrm[j - k_surf];
                // printf("sumcf[%d]=%g abscf_h2o[%d][%d]=%g %f %f \n ",i,sumcf[i],j,i,abscf_h2o[bands*j+i],dp[j-k_surf],vmrm[j-k_surf]);
                if (splitpaths) {
                    sumcf_s[0][i] -= abscf_h2o[bands * j + i] * dp[j - k_surf] * vmrm_s[0][j - k_surf];
                    sumcf_s[1][i] -= abscf_h2o[bands * j + i] * dp[j - k_surf] * vmrm_s[1][j - k_surf];
                }
            }
        }

        if (!tran_hi && !(tran_hi = (float *) malloc(nh2o * NP_HI * sizeof (float *)))) {
            printf("-E- : Error allocating memory to tran_hi\n");
            exit(FATAL_ERROR);
        }

        for (i = 0; i < nh2o; i++) {
            for (j = 0; j < bands; j++) {
                tran_hi[i * NP_HI + j] = exp(sumcf[j] * ssh2o[i] * Q * (float) 28.966 / (float) 6.0225e+23 / (float) 1.0e-6);
            }

            // Total amount of water vapor (in unit of cm) corresponding to the spectrum.

            vaptot[i] = geometry4_.vap_slant_mdl * ssh2o[i];
        }

        if (splitpaths) {
            for (j = 0; j < 2; j++) {
                for (i = 0; i < bands; i++) {
                    //                    printf("TRAN_TABLE: SSH2O_S JA: %d %g  \n",  ja,  geometry5_.ssh2o_s[j][ja]);
                    //                    printf("TRAN_TABLE: SSH2O_S JA+1: %d %g  \n",ja+1,geometry5_.ssh2o_s[j][ja+1]);
                    //                    printf("TRAN_TABLE: SSH2O_S JB: %d %g  \n",  jb,  geometry5_.ssh2o_s[j][jb]);
                    //                    printf("TRAN_TABLE: SSH2O_S JB+1: %d %g  \n",jb+1,geometry5_.ssh2o_s[j][jb+1]);
                    tran_tables_.tran_hi_sa [j][i] = exp(sumcf_s[j][i] * geometry5_.ssh2o_s[j][ja - 1] * Q * (float) 28.966 / (float) 6.0225e+23 / (float) 1.0e-6);
                    tran_tables_.tran_hi_sap1[j][i] = exp(sumcf_s[j][i] * geometry5_.ssh2o_s[j][ja ] * Q * (float) 28.966 / (float) 6.0225e+23 / (float) 1.0e-6);
                    tran_tables_.tran_hi_sb [j][i] = exp(sumcf_s[j][i] * geometry5_.ssh2o_s[j][jb - 1] * Q * (float) 28.966 / (float) 6.0225e+23 / (float) 1.0e-6);
                    tran_tables_.tran_hi_sbp1[j][i] = exp(sumcf_s[j][i] * geometry5_.ssh2o_s[j][jb ] * Q * (float) 28.966 / (float) 6.0225e+23 / (float) 1.0e-6);
                }

            }
        }
    }

    if (!fullcalc) {
        // Calculate the k-distrubition coefficients

        if (firstCall) {
            //kdist_gas_abs_(tkcdf,abscf_h2o,&nphi,wavno_hi,wavobs,&nbands) ; //FORTRAN
            kdistgasabs(tkcdf, abscf_h2o, wavno_hi, wavobs, bands, levels, nbands); //C
        }
        // Perform the k-distribution calculation Trapezoidal integral for transmittance over the NBANDS bands

        for (i = 0; i < nh2o; i++) {
            for (j = 0; j < nbands; j++) tran_table1_.tran_kd[i][j] = 1.0;

            for (j = k_surf; j < levels; j++) {

                //tkcdf[n_g][nbands][levels]
                for (n = 0; n < nbands; n++) {
                    sum_kd[n] = 0.0;
                    f1 = exp(-tkcdf[0 * nbands * levels + levels * n + j] * dp[j - k_surf] * vmrm[j - k_surf] * ssh2o[i]);

                    for (k = 0; k < n_g - 1; k++) {
                        f2 = exp(-tkcdf[(k + 1) * nbands * levels + levels * n + j] * dp[j - k_surf] * vmrm[j - k_surf] * ssh2o[i]);
                        sum_kd[n] += (f1 + f2)*(tg[k + 1] - tg[k]) / (float) 2.0;
                        f1 = f2;
                    }
                    tran_table1_.tran_kd[i][n] *= sum_kd[n];
                }
            }
            // Total amount of water vapor (in unit of cm) corresponding to the spectrum.

            vaptot[i] = geometry4_.vap_slant_mdl * ssh2o[i];
        }
    }

    if (fullcalc || firstCall) {

        // Smooth the high resolution spectra to resolutions of measured spectrum
        //tran_smooth_(tran_hi); //FORTRAN
        tran_smooth(tran_hi);
        if (!fullcalc) {
            free(wavno_hi);
            free(tran_hi);
        }

    }

    for (i = 0; i < nbands; i++) {
        tran_table1_.trntblo[i] = 1.0;
        tran_table_l2gen_.tg_solo[i] = 1.0;
        tran_table_l2gen_.tg_seno[i] = 1.0;
    }

    //    if (doSmooth) tran_smooth_others_(); //FORTRAN
    if (doSmooth) tran_smooth_others();

    if (!fullcalc) {

        if (firstCall)
            for (k = 0; k < nbands; k++)
                printf("%d ATREM: FAST)TRNTBL=%f %f\n", k, tran_table1_.tran_kd[59][k] * tran_table1_.diff_tran[0][k], tran_table1_.trntblo[k]);

        for (i = 0; i < nh2o; i++) {
            for (k = 0; k < nbands; k++) {
                /*
                     Multiplying the sensor resolution water vapor transmittance with combined
                     transmittances of CO2, N2O, CO, CH4, and O2:
                 */
                tran_table1_.trntbl[i][k] = tran_table1_.tran_kd[i][k] * tran_table1_.diff_tran[i][k] * tran_table1_.trntblo[k];
                if (firstCall)
                    printf("ATREM: DIFF_TRAN: %d %d %f %f %f %f\n", i, k,
                        tran_table1_.trntbl[i][k], tran_table1_.trntblo[k], tran_table1_.tran_kd[i][k], tran_table1_.diff_tran[i][k]);
            }
        }

    } else {

        /*
             Multiplying the sensor resolution water vapor transmittance with combined
             transmittances of CO2, N2O, CO, CH4, and O2:
         */
        for (k = 0; k < nbands; k++) {
            for (i = 0; i < nh2o; i++) {
                tran_table1_.trntbl[i][k] *= tran_table1_.trntblo[k];
            }

            if (splitpaths) {
                tran_table_l2gen_.tg_sol[k] *= tran_table_l2gen_.tg_solo[k];
                tran_table_l2gen_.tg_sen[k] *= tran_table_l2gen_.tg_seno[k];
            }

        }
    }

    //chnlratio_(); // FORTRAN
    channelRatio(); // C version

    firstCall = 0;



}

/* Converted to C from Fortran by R. Healy (richard.healy@nasa.gov) 4/2016
 *
!********************************************************************************
!*                                                                              *
!*  Name: TRAN_SMOOTH                                                           *
!*  Purpose: This program is to smooth the line-by-line high resolution         *
!*           spectrum to lower resolution spectrum that matches the resolutions *
!*           of imaging spectrometer data.                                      *
!*  Parameters: none.                                                           *
!*  Algorithm: The smoothing is done in two stages. The 1st stage is to smooth  *
!*             the high resolution spectrum to medium resolution spectrum at a  *
!*             constant FWHM (0.2 nm) and a constant wavelength interval        *
!*             (0.1 nm). The 2nd stage smoothing is to smooth the medium        *
!*             resolution spectrum to resolutions of input imaging spectrometer *
!*             data.                                                            *
!*  Globals used:  The global variables used are contained in the file          *
!*                         "COMMONS_INC"                                        *
!*  Global output: TRNCAL - total transmittances of all gases that match the    *
!*                          resolutions of imaging spectrometers.               *
!*  Return Codes:  none.                                                        *
!*                                                                              *
!********************************************************************************
 */
void tran_smooth(float *tran_hi) {

    static int firstCall = 1;

    /*
     * Converted to C from Fortran by R. Healy (richard.healy@nasa.gov) 4/2016
     *
     First stage of smoothing - smooth line-by-line high resolution spectrum with
     over 300,000 points (point spacing of 0.05 cm-1) to medium resolution
     spectrum (resolution of 0.2 nm and point spacing of 0.1 nm) with about
     25,000 points.

     The smoothing of line-by-line spectrum is done in wavenumber domain. For
     a spectrum with a constant 0.2 nm resolution in wavelength domain, it has
     variable resolution in wavenumber domain. This effect is properly taken
     care of in the design of smoothing functions.

     Because the high resolution spectrum is in wavenumber units (cm-1), while
     the medium resolution spectrum is in wavelength units, the two kinds of
     grids do not automatically match. In order to match the grids, arrays
     of INDEX_MED and TRAN_MED_INDEX are specially designed. The desired
     medium resolution spectrum, TRAN_MED, at constant 0.1 nm point spacing
     and 0.2 nm resolution is obtained through linear interpolation of
     TRAN_MED_INDEX array.

     */
    int32_t i, j, k;
    float *wavln_med_index = init_speccal13_.wavln_med_index, *fwhm = getinput4_.fwhm, *wavln_std = init_speccal12_.wavln_std;
    float *wavobs = getinput4_.wavobs;
    static float **finstr_wavno, *fwhm_wavno;
    static int32_t ncvhf_wavno[NP_MED], *ia, ia_p1;
    int32_t ncvtot_wavno, ncvtot, nbands = getinput5_.nbands, *ncvhf = init_speccal15_.ncvhf, *index_med = init_speccal13_.index_med;

    float *tran_med_sa_sol = tran_tables1_.tran_med_sa_sol, *tran_med_sap1_sol = tran_tables1_.tran_med_sap1_sol;
    float *tran_med_sa_sen = tran_tables1_.tran_med_sa_sen, *tran_med_sap1_sen = tran_tables1_.tran_med_sap1_sen;
    float *tran_med_sb_sol = tran_tables1_.tran_med_sb_sol, *tran_med_sbp1_sol = tran_tables1_.tran_med_sbp1_sol;
    float *tran_med_sb_sen = tran_tables1_.tran_med_sb_sen, *tran_med_sbp1_sen = tran_tables1_.tran_med_sbp1_sen;
    float *tran_med_index_sa_sol = tran_tables1_.tran_med_index_sa_sol, *tran_med_index_sap1_sol = tran_tables1_.tran_med_index_sap1_sol;
    float *tran_med_index_sa_sen = tran_tables1_.tran_med_index_sa_sen, *tran_med_index_sap1_sen = tran_tables1_.tran_med_index_sap1_sen;
    float *tran_med_index_sb_sol = tran_tables1_.tran_med_index_sb_sol, *tran_med_index_sbp1_sol = tran_tables1_.tran_med_index_sbp1_sol;
    float *tran_med_index_sb_sen = tran_tables1_.tran_med_index_sb_sen, *tran_med_index_sbp1_sen = tran_tables1_.tran_med_index_sbp1_sen;
    float *tran_std_sa_sol = tran_tables1_.tran_std_sa_sol, *tran_std_sap1_sol = tran_tables1_.tran_std_sap1_sol;
    float *tran_std_sa_sen = tran_tables1_.tran_std_sa_sen, *tran_std_sap1_sen = tran_tables1_.tran_std_sap1_sen;
    float *tran_std_sb_sol = tran_tables1_.tran_std_sb_sol, *tran_std_sbp1_sol = tran_tables1_.tran_std_sbp1_sol;
    float *tran_std_sb_sen = tran_tables1_.tran_std_sb_sen, *tran_std_sbp1_sen = tran_tables1_.tran_std_sbp1_sen;
    float *tg_sol = tran_table_l2gen_.tg_sol, *tg_sen = tran_table_l2gen_.tg_sen;

    float tran_ia[NH2OMAX], tran_iap1[NH2OMAX];
    float tran_ia_sa_sol, tran_ia_sa_sen, tran_ia_sap1_sol, tran_ia_sap1_sen;
    float tran_iap1_sa_sol, tran_iap1_sa_sen, tran_iap1_sap1_sol, tran_iap1_sap1_sen;
    float tran_ia_sb_sol, tran_ia_sb_sen, tran_ia_sbp1_sol, tran_ia_sbp1_sen;
    float tran_iap1_sb_sol, tran_iap1_sb_sen, tran_iap1_sbp1_sol, tran_iap1_sbp1_sen;
    float sumins, dlt, fj, fjm1, dlt_ia, fia, fia_p1;
    float **tran_std, **tran_med;
    float tg_sol_a, tg_sen_a, tg_sol_ap1, tg_sen_ap1;
    float tg_sol_b, tg_sen_b, tg_sol_bp1, tg_sen_bp1;
    float f1a = geometry_l2gen_.f1a, f2a = geometry_l2gen_.f2a;
    float f1b = geometry_l2gen_.f1b, f2b = geometry_l2gen_.f2b;
    int splitpaths = geometry_l2gen_.splitpaths;
    int32_t np_std = NP_STD, ndx1, ndx2;

    if (firstCall) {
        if (!(fwhm_wavno = (float *) malloc(NP_MED * sizeof (float)))) {
            printf("-E- : Error allocating memory to fwhm_wavno\n");
            exit(FATAL_ERROR);
        }
        if (!(finstr_wavno = (float **) malloc(NP_MED * sizeof (float*)))) {
            printf("-E- : Error allocating memory to finstr_wavno\n");
            exit(FATAL_ERROR);
        }
        for (j = 0; j < NP_MED; j++) {
            if (!(finstr_wavno[j] = (float *) malloc(NINSTRF * sizeof (float)))) {
                printf("-E- : Error allocating memory to finstr_wavno\n");
                exit(FATAL_ERROR);
            }
            fwhm_wavno[j] = (float) 10000. * DLT_MED / (wavln_med_index[j] * wavln_med_index[j]);
            ncvhf_wavno[j] = (FACDLT * fwhm_wavno[j] / DWAVNO + 1.);
            ncvtot_wavno = 2 * ncvhf_wavno[j] - 1;

            sumins = 0;
            for (i = ncvhf_wavno[j] - 1; i < ncvtot_wavno; i++) {
                finstr_wavno[j][i] = exp(-CONST1 * pow((float) (i - ncvhf_wavno[j] + 1) * DWAVNO / fwhm_wavno[j], 2));
                sumins += finstr_wavno[j][i];
                //                if (j==0) printf("A:finstr_wavno sumins: %d %g %g %d\n",i,sumins,finstr_wavno[j][i],i-ncvhf_wavno[j]+1);
            }

            for (i = 0; i < ncvhf_wavno[j] - 1; i++) {
                finstr_wavno[j][i] = finstr_wavno[j][ncvtot_wavno - i - 1];
                sumins += finstr_wavno[j][i];
                //                if (j==0) printf("B:finstr_wavno sumins: %d %g %g %d\n",i,sumins,finstr_wavno[j][i],ncvtot_wavno-i-1);
            }

            sumins *= DWAVNO;

            for (i = 0; i < ncvtot_wavno; i++) {
                finstr_wavno[j][i] *= (DWAVNO / sumins);
                //                printf("RJH: finstr_wavno %d %d %g %g \n",i,j,finstr_wavno[j][i],sumins);
            }

        }

        ia = (int32_t*) allocateMemory(nbands * sizeof (int32_t), "ia");

        //  Index searching...
        //
        //Calculate instrumental response functions...


        for (j = 0; j < nbands; j++) {
            sumins = 0.0;
            ncvtot = 2 * ncvhf[j] - 1;
            for (i = ncvhf[j] - 1; i < ncvtot; i++) {
                init_speccal15_.finstr[j][i] = exp(-CONST1 * pow((float) (i - ncvhf[j] + 1) * DWAVLN / fwhm[j], 2));
                sumins += init_speccal15_.finstr[j][i];
                //                if (j==0) printf("c:A:finstr sumins: %d %g %g %d fwhm=%g %g %g %g %g %g\n",i,sumins,init_speccal15_.finstr[j][i],i-ncvhf[j]+1,fwhm[j],DWAVLN,
                //                        (float) (i-ncvhf[j]+1)*DWAVLN/fwhm[j],pow( (float) (i-ncvhf[j]+1)*DWAVLN/fwhm[j],2),-CONST1*pow( (float) (i-ncvhf[j]+1)*DWAVLN/fwhm[j],2),
                //                        exp(-CONST1*pow( (float) (i-ncvhf[j]+1)*DWAVLN/fwhm[j],2)));
            }

            for (i = 0; i < ncvhf[j] - 1; i++) {
                init_speccal15_.finstr[j][i] = init_speccal15_.finstr[j][ncvtot - i - 1];
                sumins += init_speccal15_.finstr[j][i];
                //                if (j==0) printf("c:B:finstr sumins: %d %g %g %d\n",i,sumins,init_speccal15_.finstr[j][i],ncvtot-i-1);
            }

            sumins *= DWAVLN;

            for (i = 0; i < ncvtot; i++)
                init_speccal15_.finstr[j][i] *= (DWAVLN / sumins);

            ia[j] = hunt(wavln_std, np_std, (double) wavobs[j], ia[j]);
        }

        firstCall = 0;
    }

    if (!(tran_std = (float **) malloc(NH2OMAX * sizeof (float *)))) {
        printf("-E- : Error allocating memory to tran_std\n");
        exit(FATAL_ERROR);
    }
    if (!(tran_med = (float **) malloc(NH2OMAX * sizeof (float *)))) {
        printf("-E- : Error allocating memory to tran_med\n");
        exit(FATAL_ERROR);
    }

    for (i = 0; i < NH2OMAX; i++) {
        if (!(tran_std[i] = (float *) calloc(NP_STD, sizeof (float *)))) {
            printf("-E- : Error allocating memory to tran_std\n");
            exit(FATAL_ERROR);
        }
        if (!(tran_med[i] = (float *) calloc(NP_MED, sizeof (float *)))) {
            printf("-E- : Error allocating memory to tran_med\n");
            exit(FATAL_ERROR);
        }
    }
    for (j = 0; j < NP_MED; j++) {
        for (i = 0; i < NH2OMAX; i++)
            init_speccal13_.tran_med_index[i][j] = 0.;
        if (splitpaths) {
            tran_med_index_sa_sen [j] = 0;
            tran_med_index_sap1_sen[j] = 0;
            tran_med_index_sb_sen [j] = 0;
            tran_med_index_sbp1_sen[j] = 0;
            tran_med_index_sa_sol [j] = 0;
            tran_med_index_sap1_sol[j] = 0;
            tran_med_index_sb_sol [j] = 0;
            tran_med_index_sbp1_sol[j] = 0;
        }
        ndx1 = index_med[j]-(ncvhf_wavno[j] - 1);
        ndx2 = index_med[j] + ncvhf_wavno[j] - 1;

        for (k = ndx1 - 1; k < ndx2; k++) {
            for (i = 0; i < NH2OMAX; i++)
                init_speccal13_.tran_med_index[i][j] += tran_hi[i * NP_HI + k] * finstr_wavno[j][k - index_med[j] + ncvhf_wavno[j]];


            //                if (i==29)
            //                    printf(" TRAN_MED_INDEX: %d %d <%g> %g %g %d %d %d\n",j,k,init_speccal13_.tran_med_index[i][j],tran_hi[i*NP_HI+k],
            //                                  finstr_wavno[j][k-index_med[j]+ncvhf_wavno[j]],
            //                            index_med[j],ncvhf_wavno[j],k-index_med[j]+ncvhf_wavno[j]);
            if (splitpaths) {
                //Note that the tran_hi_* arrays are 0 for solar path and 1 for sensor
                tran_med_index_sa_sol [j] += tran_tables_.tran_hi_sa [0][k] * finstr_wavno[j][k - index_med[j] + ncvhf_wavno[j]];
                tran_med_index_sa_sen [j] += tran_tables_.tran_hi_sa [1][k] * finstr_wavno[j][k - index_med[j] + ncvhf_wavno[j]];
                tran_med_index_sap1_sol[j] += tran_tables_.tran_hi_sap1[0][k] * finstr_wavno[j][k - index_med[j] + ncvhf_wavno[j]];
                tran_med_index_sap1_sen[j] += tran_tables_.tran_hi_sap1[1][k] * finstr_wavno[j][k - index_med[j] + ncvhf_wavno[j]];
                tran_med_index_sb_sol [j] += tran_tables_.tran_hi_sb [0][k] * finstr_wavno[j][k - index_med[j] + ncvhf_wavno[j]];
                tran_med_index_sb_sen [j] += tran_tables_.tran_hi_sb [1][k] * finstr_wavno[j][k - index_med[j] + ncvhf_wavno[j]];
                tran_med_index_sbp1_sol[j] += tran_tables_.tran_hi_sbp1[0][k] * finstr_wavno[j][k - index_med[j] + ncvhf_wavno[j]];
                tran_med_index_sbp1_sen[j] += tran_tables_.tran_hi_sbp1[1][k] * finstr_wavno[j][k - index_med[j] + ncvhf_wavno[j]];
                //                if (j==8204) printf("SPLIT: TRAN_MED_INDEX: %d %d <%g> %g %g %d %d %d %g\n",j,k,tran_med_index_sa_sol[j],tran_tables_.tran_hi_sa[0][k],
                //                        finstr_wavno[j][k-index_med[j]+ncvhf_wavno[j]],
                //                        index_med[j],ncvhf_wavno[j],k-index_med[j]+ncvhf_wavno[j],
                //                        tran_tables_.tran_hi_sa  [0][k]*finstr_wavno[j][k-index_med[j]+ncvhf_wavno[j]]);
            }
        }
    }
    /*
     *
     Linear interpolation to get TRAN_MED from TRAN_MED_INDEX:
     (Note that WAVLN_MED_INDEX(J) >= WAVLN_MED(J)    )
     *
     */
    for (i = 0; i < NH2OMAX; i++) {
        tran_med[i][0] = init_speccal13_.tran_med_index[i][0];
        tran_med[i][NP_MED - 1] = init_speccal13_.tran_med_index[i][NP_MED - 1];
    }
    if (splitpaths) {
        tran_med_sa_sol [0] = tran_med_index_sa_sol [0];
        tran_med_sa_sol [NP_MED - 1] = tran_med_index_sa_sol [NP_MED - 1];
        tran_med_sa_sen [0] = tran_med_index_sa_sen [0];
        tran_med_sa_sen [NP_MED - 1] = tran_med_index_sa_sen [NP_MED - 1];

        tran_med_sap1_sol[0] = tran_med_index_sap1_sol[0];
        tran_med_sap1_sol[NP_MED - 1] = tran_med_index_sap1_sol[NP_MED - 1];
        tran_med_sap1_sen[0] = tran_med_index_sap1_sen[0];
        tran_med_sap1_sen[NP_MED - 1] = tran_med_index_sap1_sen[NP_MED - 1];

        tran_med_sb_sol [0] = tran_med_index_sb_sol [0];
        tran_med_sb_sol [NP_MED - 1] = tran_med_index_sb_sol [NP_MED - 1];
        tran_med_sb_sen [0] = tran_med_index_sb_sen [0];
        tran_med_sb_sen [NP_MED - 1] = tran_med_index_sb_sen [NP_MED - 1];

        tran_med_sbp1_sol[0] = tran_med_index_sbp1_sol[0];
        tran_med_sbp1_sol[NP_MED - 1] = tran_med_index_sbp1_sol[NP_MED - 1];
        tran_med_sbp1_sen[0] = tran_med_index_sbp1_sen[0];
        tran_med_sbp1_sen[NP_MED - 1] = tran_med_index_sbp1_sen[NP_MED - 1];
    }

    for (j = 1; j < NP_MED - 1; j++) {
        if (wavln_med_index[j] <= init_speccal12_.wavln_med[j]) {
            for (i = 0; i < NH2OMAX; i++)
                tran_med[i][j] = init_speccal13_.tran_med_index[i][j];
            //                printf(" TRAN_MED_WAVLN: %d %d <%g> <%g> %g %g\n",j,29,tran_med[29][j],init_speccal13_.tran_med_index[29][j],
            //                        wavln_med_index[j],init_speccal12_.wavln_med[j]);
            if (splitpaths) {
                tran_med_sa_sol [j] = tran_med_index_sa_sol [j];
                tran_med_sa_sen [j] = tran_med_index_sa_sen [j];
                tran_med_sap1_sol [j] = tran_med_index_sap1_sol [j];
                tran_med_sap1_sen [j] = tran_med_index_sap1_sen [j];
                tran_med_sb_sol [j] = tran_med_index_sb_sol [j];
                tran_med_sb_sen [j] = tran_med_index_sb_sen [j];
                tran_med_sbp1_sol [j] = tran_med_index_sbp1_sol [j];
                tran_med_sbp1_sen [j] = tran_med_index_sbp1_sen [j];
            }
        } else {
            dlt = wavln_med_index[j] - wavln_med_index[j - 1];
            fjm1 = (wavln_med_index[j] - init_speccal12_.wavln_med[j]) / dlt;
            fj = (init_speccal12_.wavln_med [j] - wavln_med_index[j - 1]) / dlt;
            for (i = 0; i < NH2OMAX; i++) {
                tran_med[i][j] = fjm1 * init_speccal13_.tran_med_index[i][j - 1] + fj * init_speccal13_.tran_med_index[i][j];
                //                if (i==29)
                //                    printf(" TRAN_MED_INTERP: %d %d <%g> <%g> %g %g %g\n",j,i,init_speccal13_.tran_med_index[i][j-1],init_speccal13_.tran_med_index[i][j],
                //                            fjm1,fj,dlt);
            }
            if (splitpaths) {
                tran_med_sa_sol [j] = fjm1 * tran_med_index_sa_sol [j - 1] + fj * tran_med_index_sa_sol [j];
                tran_med_sa_sen [j] = fjm1 * tran_med_index_sa_sen [j - 1] + fj * tran_med_index_sa_sen [j];
                tran_med_sap1_sol [j] = fjm1 * tran_med_index_sap1_sol[j - 1] + fj * tran_med_index_sap1_sol[j];
                tran_med_sap1_sen [j] = fjm1 * tran_med_index_sap1_sen[j - 1] + fj * tran_med_index_sap1_sen[j];
                tran_med_sb_sol [j] = fjm1 * tran_med_index_sb_sol [j - 1] + fj * tran_med_index_sb_sol [j];
                tran_med_sb_sen [j] = fjm1 * tran_med_index_sb_sen [j - 1] + fj * tran_med_index_sb_sen [j];
                tran_med_sbp1_sol [j] = fjm1 * tran_med_index_sbp1_sol[j - 1] + fj * tran_med_index_sbp1_sol[j];
                tran_med_sbp1_sen [j] = fjm1 * tran_med_index_sbp1_sen[j - 1] + fj * tran_med_index_sbp1_sen[j];
                //                printf("SPLIT: TRAN_MED_SA_SOL= %d <%g> <%g %g> %g %g\n",j,tran_med_sa_sol[j],
                //                        fjm1,fj,tran_med_index_sa_sol[j-1],tran_med_index_sa_sol[j]);
            }
        }
    }

    /*
     * --- Here multiplying O3 and NO2 spectra and other spectrum at medium resolution:
     *
     */
    for (j = 0; j < NH2OMAX; j++) {
        for (i = 0; i < NPSHIF; i++)
            tran_std[j][i] = 1.0;
        for (i = NPSHIF; i < NP_STD; i++)
            tran_std[j][i] = tran_med[j][i - NPSHIF];

    }
    if (splitpaths) {
        for (i = 0; i < NP_STD; i++) {
            tran_std_sa_sol [i] = 1.0;
            tran_std_sa_sen [i] = 1.0;
            tran_std_sap1_sol[i] = 1.0;
            tran_std_sap1_sen[i] = 1.0;
            tran_std_sb_sol [i] = 1.0;
            tran_std_sb_sen [i] = 1.0;
            tran_std_sbp1_sol[i] = 1.0;
            tran_std_sbp1_sen[i] = 1.0;
        }
        for (i = NPSHIF; i < NP_STD; i++) {
            tran_std_sa_sol [i] *= tran_med_sa_sol [i - NPSHIF];
            tran_std_sa_sen [i] *= tran_med_sa_sen [i - NPSHIF];
            tran_std_sap1_sol[i] *= tran_med_sap1_sol[i - NPSHIF];
            tran_std_sap1_sen[i] *= tran_med_sap1_sen[i - NPSHIF];
            tran_std_sb_sol [i] *= tran_med_sb_sol [i - NPSHIF];
            tran_std_sb_sen [i] *= tran_med_sb_sen [i - NPSHIF];
            tran_std_sbp1_sol[i] *= tran_med_sbp1_sol[i - NPSHIF];
            tran_std_sbp1_sen[i] *= tran_med_sbp1_sen[i - NPSHIF];
            //                printf("%d TRAN_STD_SA_SOL= %g %g\n",i,tran_std_sa_sol[i],tran_med_sa_sol[i-NPSHIF]);

        }
    }
    /*
     * The 2nd stage of smoothing - smooth the medium resolution spectrum (resolution
       of 0.2 nm and point spacing of 0.1 nm) with about 25,000 points to match
       the coarser and variable resolution spectrum from imaging spectrometers.

       Initialize some index parameters:

     */

    for (j = 0; j < nbands; j++) {

        for (i = 0; i < NH2OMAX; i++) {
            tran_table1_.trntbl[i][j] = 0.0;
            tran_ia[i] = 0.0;
            tran_iap1[i] = 0.0;
        }

        if (splitpaths) {
            tg_sol[j] = 0.0;
            tg_sen[j] = 0.0;
            tran_ia_sa_sol = 0.0;
            tran_ia_sa_sen = 0.0;
            tran_ia_sap1_sol = 0.0;
            tran_ia_sap1_sen = 0.0;
            tran_iap1_sa_sol = 0.0;
            tran_iap1_sa_sen = 0.0;
            tran_iap1_sap1_sol = 0.0;
            tran_iap1_sap1_sen = 0.0;
            tran_ia_sb_sol = 0.0;
            tran_ia_sb_sen = 0.0;
            tran_ia_sbp1_sol = 0.0;
            tran_ia_sbp1_sen = 0.0;
            tran_iap1_sb_sol = 0.0;
            tran_iap1_sb_sen = 0.0;
            tran_iap1_sbp1_sol = 0.0;
            tran_iap1_sbp1_sen = 0.0;

        }

        // Smoothing ...
        for (k = ia[j]-(ncvhf[j] - 1); k < ia[j] + ncvhf[j]; k++) {
            for (i = 0; i < NH2OMAX; i++) {
                tran_ia[i] += tran_std[i][k] * init_speccal15_.finstr[j][k - ia[j] + ncvhf[j] - 1];
                //                if(i==29) printf("TRAN_IA: %d %d %d %g %g %d %g %d\n",j,k,ia[j],tran_ia[i],tran_std[i][k],k-ia[j]+ncvhf[j]-1,init_speccal15_.finstr[j][k-ia[j]+ncvhf[j]-1],ncvhf[j]);
            }
            if (splitpaths) {
                tran_ia_sa_sol += tran_std_sa_sol [k] * init_speccal15_.finstr[j][k - ia[j] + ncvhf[j] - 1];
                tran_ia_sa_sen += tran_std_sa_sen [k] * init_speccal15_.finstr[j][k - ia[j] + ncvhf[j] - 1];
                tran_ia_sap1_sol += tran_std_sap1_sol[k] * init_speccal15_.finstr[j][k - ia[j] + ncvhf[j] - 1];
                tran_ia_sap1_sen += tran_std_sap1_sen[k] * init_speccal15_.finstr[j][k - ia[j] + ncvhf[j] - 1];
                tran_ia_sb_sol += tran_std_sb_sol [k] * init_speccal15_.finstr[j][k - ia[j] + ncvhf[j] - 1];
                tran_ia_sb_sen += tran_std_sb_sen [k] * init_speccal15_.finstr[j][k - ia[j] + ncvhf[j] - 1];
                tran_ia_sbp1_sol += tran_std_sbp1_sol[k] * init_speccal15_.finstr[j][k - ia[j] + ncvhf[j] - 1];
                tran_ia_sbp1_sen += tran_std_sbp1_sen[k] * init_speccal15_.finstr[j][k - ia[j] + ncvhf[j] - 1];
                //                printf("%d %d TRAN_IA_SA_SOL= %g %g %g %d\n",j,k,tran_ia_sa_sol, tran_std_sa_sol  [k],init_speccal15_.finstr[j][k-ia[j]+ncvhf[j]-1],ia[j]);
            }
        }

        ia_p1 = ia[j] + 1;
        for (k = ia_p1 - (ncvhf[j] - 1); k < ia_p1 + ncvhf[j]; k++) {
            for (i = 0; i < NH2OMAX; i++) {
                tran_iap1[i] += tran_std[i][k] * init_speccal15_.finstr[j][k - ia_p1 + ncvhf[j] - 1];
            }
            if (splitpaths) {
                tran_iap1_sa_sol += tran_std_sa_sol [k] * init_speccal15_.finstr[j][k - ia_p1 + ncvhf[j] - 1];
                tran_iap1_sa_sen += tran_std_sa_sen [k] * init_speccal15_.finstr[j][k - ia_p1 + ncvhf[j] - 1];
                tran_iap1_sap1_sol += tran_std_sap1_sol[k] * init_speccal15_.finstr[j][k - ia_p1 + ncvhf[j] - 1];
                tran_iap1_sap1_sen += tran_std_sap1_sen[k] * init_speccal15_.finstr[j][k - ia_p1 + ncvhf[j] - 1];
                tran_iap1_sb_sol += tran_std_sb_sol [k] * init_speccal15_.finstr[j][k - ia_p1 + ncvhf[j] - 1];
                tran_iap1_sb_sen += tran_std_sb_sen [k] * init_speccal15_.finstr[j][k - ia_p1 + ncvhf[j] - 1];
                tran_iap1_sbp1_sol += tran_std_sbp1_sol[k] * init_speccal15_.finstr[j][k - ia_p1 + ncvhf[j] - 1];
                tran_iap1_sbp1_sen += tran_std_sbp1_sen[k] * init_speccal15_.finstr[j][k - ia_p1 + ncvhf[j] - 1];
                //                printf("%d %d TRAN_IAP1_SA_SOL= %g %g %g\n",j,k,tran_iap1_sa_sol, tran_std_sa_sol  [k],init_speccal15_.finstr[j][k-ia[j]+ncvhf[j]-1]);
            }
        }

        // Linear interpolation to get TRNCAL from TRAN_IA and TRAN_IAP1:
        dlt_ia = wavln_std[ia_p1] - wavln_std[ia[j]];
        fia = (wavln_std[ia_p1] - wavobs[j]) / dlt_ia;
        fia_p1 = (float) 1.0 - fia;

        for (i = 0; i < NH2OMAX; i++)
            tran_table1_.trntbl[i][j] = fia * tran_ia[i] + fia_p1 * tran_iap1[i];

        if (splitpaths) {
            tg_sol_a = fia * tran_ia_sa_sol + fia_p1*tran_iap1_sa_sol;
            tg_sol_ap1 = fia * tran_ia_sap1_sol + fia_p1*tran_iap1_sap1_sol;
            tg_sen_a = fia * tran_ia_sa_sen + fia_p1*tran_iap1_sa_sen;
            tg_sen_ap1 = fia * tran_ia_sap1_sen + fia_p1*tran_iap1_sap1_sen;

            tg_sol_b = fia * tran_ia_sb_sol + fia_p1*tran_iap1_sb_sol;
            tg_sol_bp1 = fia * tran_ia_sbp1_sol + fia_p1*tran_iap1_sbp1_sol;
            tg_sen_b = fia * tran_ia_sb_sen + fia_p1*tran_iap1_sb_sen;
            tg_sen_bp1 = fia * tran_ia_sbp1_sen + fia_p1*tran_iap1_sbp1_sen;

            tg_sol[j] = (float) 0.5 * (f1a * tg_sol_a + f2a * tg_sol_ap1 +
                    f1b * tg_sol_b + f2b * tg_sol_bp1);
            tg_sen[j] = (float) 0.5 * (f1a * tg_sen_a + f2a * tg_sen_ap1 +
                    f1b * tg_sen_b + f2b * tg_sen_bp1);
            //            printf("%d %d TG_SOL_SA= %g %g %g %g %g %g %g <%g %g %g %g>\n",j,ia[j],tg_sol[j],tg_sol_a,tg_sol_ap1,
            //                    tran_ia_sa_sol,tran_iap1_sa_sol,
            //                    tran_ia_sap1_sol,tran_iap1_sap1_sol,
            //                    f1a,f2a,f1b,f2b);
            //            printf("%d %d TG_SOL_SB= %g %g %g %g %g %g %g <%g %g %g %g>\n",j,ia[j],tg_sol[j],tg_sol_b,tg_sol_bp1,
            //                    tran_ia_sb_sol,tran_iap1_sb_sol,
            //                    tran_ia_sbp1_sol,tran_iap1_sbp1_sol,
            //                    f1a,f2a,f1b,f2b);
        }

        if (!getinput5_.full_calc) {
            for (k = 0; k < NH2OMAX; k++) {
                tran_table1_.diff_tran[k][j] = (tran_table1_.trntbl[k][j] - tran_table1_.tran_kd[k][j]) /
                        tran_table1_.tran_kd[k][j]+(float) 1.0;
                //                      printf("RJH: TRAN2: %d %d %f %f %f\n",j,k,tran_table1_.tran_kd[k][j],tran_table1_.trntbl[k][j],tran_table1_.diff_tran[k][j]);
            }
        }
    } //nbands

    for (i = 0; i < NH2OMAX; i++) {
        free(tran_std[i]);
        free(tran_med[i]);
    }
    free(tran_std);
    free(tran_med);
}

/* Converted to C from Fortran by R. Healy (richard.healy@nasa.gov) 4/2016
 *
!********************************************************************************
!*                                                                              *
!*  Name: TRAN_SMOOTH_OTHERS                                                    *
!*  Purpose: This program is to smooth the line-by-line high resolution         *
!*           spectrum to lower resolution spectrum that matches the resolutions *
!*           of imaging spectrometer data - for other absorbing gases           *
!*  Parameters: none.                                                           *
!*  Algorithm: The smoothing is done in two stages. The 1st stage is to smooth  *
!*             the high resolution spectrum to medium resolution spectrum at a  *
!*             constant FWHM (0.2 nm) and a constant wavelength interval        *
!*             (0.1 nm). The 2nd stage smoothing is to smooth the medium        *
!*             resolution spectrum to resolutions of input imaging spectrometer *
!*             data.                                                            *
!*  Globals used:  The global variables used are contained in the file          *
!*                         "COMMONS_INC"                                        *
!*  Global output: TRNCAL - total transmittances of all gases that match the    *
!*                          resolutions of imaging spectrometers.               *
!*  Return Codes:  none.                                                        *
!*                                                                              *
!********************************************************************************
 */
void tran_smooth_others() {

    static int firstCall = 1;

    /*
     * Converted to C from Fortran by R. Healy (richard.healy@nasa.gov) 4/2016
     *
     First stage of smoothing - smooth line-by-line high resolution spectrum with
     over 300,000 points (point spacing of 0.05 cm-1) to medium resolution
     spectrum (resolution of 0.2 nm and point spacing of 0.1 nm) with about
     25,000 points.

     The smoothing of line-by-line spectrum is done in wavenumber domain. For
     a spectrum with a constant 0.2 nm resolution in wavelength domain, it has
     variable resolution in wavenumber domain. This effect is properly taken
     care of in the design of smoothing functions.

     Because the high resolution spectrum is in wavenumber units (cm-1), while
     the medium resolution spectrum is in wavelength units, the two kinds of
     grids do not automatically match. In order to match the grids, arrays
     of INDEX_MED and TRAN_MED_INDEX are specially designed. The desired
     medium resolution spectrum, TRAN_MED, at constant 0.1 nm point spacing
     and 0.2 nm resolution is obtained through linear interpolation of
     TRAN_MED_INDEX array.

     */
    int32_t i, j, k;
    float *wavln_med_index = init_speccal13_.wavln_med_index, *fwhm = getinput4_.fwhm, *wavln_std = init_speccal12_.wavln_std;
    float *wavobs = getinput4_.wavobs;
    static float **finstr_wavno, *fwhm_wavno;
    static int32_t *ncvhf_wavno, *ia, ia_p1;
    int32_t ncvtot_wavno, ncvtot, nbands = getinput5_.nbands, *ncvhf = init_speccal15_.ncvhf, *index_med = init_speccal13_.index_med;

    float *trntblo = tran_table1_.trntblo;
    float *tg_solo = tran_table_l2gen_.tg_solo, *tg_seno = tran_table_l2gen_.tg_seno;

    float *tran_hi_others = init_speccal1_.tran_hi_others;
    float tran_ia, tran_iap1;
    float tran_ias_sol, tran_ias_sen, tran_iap1s_sol, tran_iap1s_sen;
    float sumins, dlt, fj, fjm1, dlt_ia, fia, fia_p1;
    //   Transmittance of medium resolution data
    float *tran_std_o, *tran_med_o, *tran_med_index_o;
    float *tran_std_os_sol, *tran_std_os_sen;
    float *tran_med_os_sol, *tran_med_os_sen, *tran_med_index_os_sol, *tran_med_index_os_sen;
    int splitpaths = geometry_l2gen_.splitpaths;
    float airmass, airm, mu = geometry5_.mu, mu0 = geometry5_.mu0;
    float tran_hi_others_sol, tran_hi_others_sen;


    int32_t np_std = NP_STD, ndx1, ndx2;

    if (firstCall) {
        if (!(fwhm_wavno = (float *) malloc(NP_MED * sizeof (float)))) {
            printf("-E- : Error allocating memory to fwhm_wavno\n");
            exit(FATAL_ERROR);
        }
        if (!(finstr_wavno = (float **) malloc(NP_MED * sizeof (float*)))) {
            printf("-E- : Error allocating memory to finstr_wavno\n");
            exit(FATAL_ERROR);
        }
        ncvhf_wavno = (int32_t*) allocateMemory(NP_MED * sizeof (int32_t), "ncvhf_wavno");

        for (j = 0; j < NP_MED; j++) {
            if (!(finstr_wavno[j] = (float *) malloc(NINSTRF * sizeof (float)))) {
                printf("-E- : Error allocating memory to finstr_wavno\n");
                exit(FATAL_ERROR);
            }

            fwhm_wavno[j] = (float) 10000. * DLT_MED / (wavln_med_index[j] * wavln_med_index[j]);
            ncvhf_wavno[j] = (FACDLT * fwhm_wavno[j] / DWAVNO + 1.);
            ncvtot_wavno = 2 * ncvhf_wavno[j] - 1;

            sumins = 0;
            for (i = ncvhf_wavno[j] - 1; i < ncvtot_wavno; i++) {
                finstr_wavno[j][i] = exp(-CONST1 * pow((float) (i - ncvhf_wavno[j] + 1) * DWAVNO / fwhm_wavno[j], 2));
                sumins += finstr_wavno[j][i];
                //                if (j==0) printf("A:finstr_wavno sumins: %d %g %g %d\n",i,sumins,finstr_wavno[j][i],i-ncvhf_wavno[j]+1);
            }

            for (i = 0; i < ncvhf_wavno[j] - 1; i++) {
                finstr_wavno[j][i] = finstr_wavno[j][ncvtot_wavno - i - 1];
                sumins += finstr_wavno[j][i];
                //                if (j==0) printf("B:finstr_wavno sumins: %d %g %g %d\n",i,sumins,finstr_wavno[j][i],ncvtot_wavno-i-1);
            }

            sumins *= DWAVNO;

            for (i = 0; i < ncvtot_wavno; i++) {
                finstr_wavno[j][i] *= (DWAVNO / sumins);
                //                printf("RJH: finstr_wavno %d %d %g %g \n",i,j,finstr_wavno[j][i],sumins);
            }

        }

        ia = (int32_t*) allocateMemory(nbands * sizeof (int32_t), "ia");
        //Calculate instrumental response functions...

        for (j = 0; j < nbands; j++) {
            sumins = 0.0;
            ncvtot = 2 * ncvhf[j] - 1;
            //            if (j==0) printf("NCVTOT: %d %d\n",ncvtot,ncvhf[j]);
            for (i = ncvhf[j] - 1; i < ncvtot; i++) {
                init_speccal15_.finstr[j][i] = exp(-CONST1 * pow((float) (i - ncvhf[j] + 1) * DWAVLN / fwhm[j], 2));
                sumins += init_speccal15_.finstr[j][i];
                //                if (j==0) printf("c:A:finstr sumins: %d %g %g %d fwhm=%g %g %g %g %g %g\n",i,sumins,init_speccal15_.finstr[j][i],i-ncvhf[j]+1,fwhm[j],DWAVLN,
                //                        (float) (i-ncvhf[j]+1)*DWAVLN/fwhm[j],pow( (float) (i-ncvhf[j]+1)*DWAVLN/fwhm[j],2),-CONST1*pow( (float) (i-ncvhf[j]+1)*DWAVLN/fwhm[j],2),
                //                        exp(-CONST1*pow( (float) (i-ncvhf[j]+1)*DWAVLN/fwhm[j],2)));
            }

            for (i = 0; i < ncvhf[j] - 1; i++) {
                init_speccal15_.finstr[j][i] = init_speccal15_.finstr[j][ncvtot - i - 1];
                sumins += init_speccal15_.finstr[j][i];
                //                if (j==0) printf("c:B:finstr sumins: %d %g %g %d\n",i,sumins,init_speccal15_.finstr[j][i],ncvtot-i-1);
            }

            sumins *= DWAVLN;

            for (i = 0; i < ncvtot; i++)
                init_speccal15_.finstr[j][i] *= (DWAVLN / sumins);

            ia[j] = hunt(wavln_std, np_std, (double) wavobs[j], ia[j]);
        }

        firstCall = 0;
    }

    airmass = 1.0 / mu0 + 1.0 / mu;

    /*
     *  !!!High resolution transmittances for CO2, N2O, CO, CH4, and O2 should
         also be calculated for wavelength start = 0.56 micron.
     *
     */

    if (!(tran_std_o = (float *) calloc(NP_STD, sizeof (float)))) {
        printf("-E- : Error allocating memory to tran_std_o\n");
        exit(FATAL_ERROR);
    }
    if (!(tran_std_os_sol = (float *) calloc(NP_STD, sizeof (float)))) {
        printf("-E- : Error allocating memory to tran_std_os_sol\n");
        exit(FATAL_ERROR);
    }
    if (!(tran_std_os_sen = (float *) calloc(NP_STD, sizeof (float)))) {
        printf("-E- : Error allocating memory to tran_std_os_sen\n");
        exit(FATAL_ERROR);
    }
    if (!(tran_med_o = (float *) calloc(NP_MED, sizeof (float)))) {
        printf("-E- : Error allocating memory to tran_med_o\n");
        exit(FATAL_ERROR);
    }
    if (!(tran_med_index_o = (float *) calloc(NP_MED, sizeof (float)))) {
        printf("-E- : Error allocating memory to tran_med_index_o\n");
        exit(FATAL_ERROR);
    }
    if (!(tran_med_os_sol = (float *) calloc(NP_MED, sizeof (float)))) {
        printf("-E- : Error allocating memory to tran_med_os_sol\n");
        exit(FATAL_ERROR);
    }
    if (!(tran_med_os_sen = (float *) calloc(NP_MED, sizeof (float)))) {
        printf("-E- : Error allocating memory to tran_med_os_sen\n");
        exit(FATAL_ERROR);
    }
    if (!(tran_med_index_os_sol = (float *) calloc(NP_MED, sizeof (float)))) {
        printf("-E- : Error allocating memory to tran_med_index_os_sol\n");
        exit(FATAL_ERROR);
    }
    if (!(tran_med_index_os_sen = (float *) calloc(NP_MED, sizeof (float)))) {
        printf("-E- : Error allocating memory to tran_med_index_os_sen\n");
        exit(FATAL_ERROR);
    }

    for (j = 0; j < NP_MED; j++) {
        ndx1 = index_med[j]-(ncvhf_wavno[j] - 1);
        ndx2 = index_med[j] + ncvhf_wavno[j] - 1;

        //printf(" NDX1,NDX2= %d %d %d %d\n",ndx1,ndx2,index_med[j],ncvhf_wavno[j]);
        for (k = ndx1 - 1; k < ndx2; k++) {
            tran_med_index_o[j] += tran_hi_others[k] * finstr_wavno[j][k - index_med[j] + ncvhf_wavno[j]];

            //                    printf(" TRAN_MED_INDEX_O: %d %d <%g> %g %g %d %d %d %g\n",j,k,tran_med_index_o[j],tran_hi_others[k],
            //                            finstr_wavno[j][k-index_med[j]+ncvhf_wavno[j]],
            //                            index_med[j],ncvhf_wavno[j],k-index_med[j]+ncvhf_wavno[j],tran_hi_others[k]*finstr_wavno[j][k-index_med[j]+ncvhf_wavno[j]]);
            if (splitpaths) {
                //Note that the tran_hi_* arrays are 0 for solar path and 1 for sensor
                airm = -(float) log(tran_hi_others[k]) / airmass;
                tran_hi_others_sol = (float) exp(-airm / mu0);
                tran_hi_others_sen = (float) exp(-airm / mu);
                tran_med_index_os_sol[j] += tran_hi_others_sol * finstr_wavno[j][k - index_med[j] + ncvhf_wavno[j]];
                tran_med_index_os_sen[j] += tran_hi_others_sen * finstr_wavno[j][k - index_med[j] + ncvhf_wavno[j]];

            }
        }
    }
    /*
     *
     Linear interpolation to get TRAN_MED from TRAN_MED_INDEX:
     (Note that WAVLN_MED_INDEX(J) >= WAVLN_MED(J)    )
     *
     */
    tran_med_o[0] = tran_med_index_o[0];
    tran_med_o[NP_MED - 1] = tran_med_index_o[NP_MED - 1];

    if (splitpaths) {
        tran_med_os_sol[0] = tran_med_index_os_sol[0];
        tran_med_os_sol[NP_MED - 1] = tran_med_index_os_sol[NP_MED - 1];
        tran_med_os_sen[0] = tran_med_index_os_sen[0];
        tran_med_os_sen[NP_MED - 1] = tran_med_index_os_sen[NP_MED - 1];
    }

    for (j = 1; j < NP_MED - 1; j++) {
        if (wavln_med_index[j] <= init_speccal12_.wavln_med[j]) {
            tran_med_o[j] = tran_med_index_o[j];
            //                printf(" TRAN_MED_WAVLN: %d <%g> <%g> %g %g\n",j,tran_med_o[j],tran_med_index_o][j],
            //                        wavln_med_index[j],init_speccal12_.wavln_med[j]);
            if (splitpaths) {
                tran_med_os_sol[j] = tran_med_index_os_sol[j];
                tran_med_os_sen[j] = tran_med_index_os_sen[j];
            }
        } else {
            dlt = wavln_med_index[j] - wavln_med_index[j - 1];
            fjm1 = (wavln_med_index[j] - init_speccal12_.wavln_med[j]) / dlt;
            fj = (init_speccal12_.wavln_med [j] - wavln_med_index[j - 1]) / dlt;
            tran_med_o[j] = fjm1 * tran_med_index_o[j - 1] + fj * tran_med_index_o[j];
            //                    printf(" TRAN_MED_INTERP: %d <%g> <%g> %g %g %g\n",j,tran_med_index_o[j-1],tran_med_index_o[j],
            //                            fjm1,fj,dlt);
            if (splitpaths) {
                tran_med_os_sol[j] = fjm1 * tran_med_index_os_sol[j - 1] + fj * tran_med_index_os_sol[j];
                tran_med_os_sen[j] = fjm1 * tran_med_index_os_sen[j - 1] + fj * tran_med_index_os_sen[j];
                //                printf(" TRAN_MED_INTERP_OS: %d <%g> <%g> %g %g %g\n",j,tran_med_index_os_sol[j-1],tran_med_index_os_sol[j],
                //                          fjm1,fj,dlt);
            }
        }
    }

    /*
     * --- Here multiplying O3 and NO2 spectra and other spectrum at medium resolution:
     *
     */

    for (i = 0; i < NP_STD; i++) tran_std_o[i] = 1.0;

    if (splitpaths) {
        for (i = 0; i < NP_STD; i++) {
            tran_std_os_sol[i] = 1.0;
            tran_std_os_sen[i] = 1.0;
        }
    }

    for (i = 0; i < NO3PT; i++) tran_std_o[i] = init_speccal16_.tran_o3_std[i] * init_speccal17_.tran_no2_std[i];

    if (splitpaths) {
        for (i = 0; i < NO3PT; i++) {
            airm = -(float) log(tran_std_o[i]) / airmass;
            tran_std_os_sol[i] = (float) exp(-airm / mu0);
            tran_std_os_sen[i] = (float) exp(-airm / mu);
        }
    }

    for (i = NPSHIF; i < NP_STD; i++) tran_std_o[i] *= tran_med_o[i - NPSHIF];

    if (splitpaths) {
        for (i = NPSHIF; i < NP_STD; i++) {
            tran_std_os_sol[i] *= tran_med_os_sol[i - NPSHIF];
            tran_std_os_sen[i] *= tran_med_os_sen[i - NPSHIF];
            //                printf("%d TRAN_STD_OS_SOL= %g %g\n",i,tran_std_os_sol[i],tran_med_os_sol[i-NPSHIF]);
            //                printf("%d TRAN_STD_OS_SEN= %g %g\n",i,tran_std_os_sen[i],tran_med_os_sen[i-NPSHIF]);
        }
    }

    /*
     * The 2nd stage of smoothing - smooth the medium resolution spectrum (resolution
       of 0.2 nm and point spacing of 0.1 nm) with about 25,000 points to match
       the coarser and variable resolution spectrum from imaging spectrometers.

       Initialize some index parameters:

     */

    for (j = 0; j < nbands; j++) {

        tran_ia = 0.0;
        tran_iap1 = 0.0;

        if (splitpaths) {
            tran_ias_sol = 0.0;
            tran_iap1s_sol = 0.0;
            tran_ias_sen = 0.0;
            tran_iap1s_sen = 0.0;
        }

        // Smoothing ...
        for (k = ia[j]-(ncvhf[j] - 1); k < ia[j] + ncvhf[j]; k++) {
            tran_ia += tran_std_o[k] * init_speccal15_.finstr[j][k - ia[j] + ncvhf[j] - 1];
            //                printf("TRAN_IA: %d %d %d %g %g %d %g %d\n",j,k,ia[j],tran_ia,tran_std_o[k],k-ia[j]+ncvhf[j]-1,init_speccal15_.finstr[j][k-ia[j]+ncvhf[j]-1],ncvhf[j]);

            if (splitpaths) {
                tran_ias_sol += tran_std_os_sol [k] * init_speccal15_.finstr[j][k - ia[j] + ncvhf[j] - 1];
                tran_ias_sen += tran_std_os_sen [k] * init_speccal15_.finstr[j][k - ia[j] + ncvhf[j] - 1];
                //                                printf("TRAN_IAS_SOL: %d %d %g %g %d %g %d %d\n",j,k,tran_ias_sol,tran_std_os_sol[k],k-ia[j]+ncvhf[j]-1,mu0,ia[j],ncvhf[j]);
                //                                printf("TRAN_IAS_SEN: %d %d %g %g %d %g %d %d\n",j,k,tran_ias_sen,tran_std_os_sen[k],k-ia[j]+ncvhf[j]-1,mu,ia[j],ncvhf[j]);
            }
        }

        ia_p1 = ia[j] + 1;
        for (k = ia_p1 - (ncvhf[j] - 1); k < ia_p1 + ncvhf[j]; k++) {
            tran_iap1 += tran_std_o[k] * init_speccal15_.finstr[j][k - ia_p1 + ncvhf[j] - 1];
            if (splitpaths) {
                tran_iap1s_sol += tran_std_os_sol [k] * init_speccal15_.finstr[j][k - ia_p1 + ncvhf[j] - 1];
                tran_iap1s_sen += tran_std_os_sen [k] * init_speccal15_.finstr[j][k - ia_p1 + ncvhf[j] - 1];
            }
        }

        // Linear interpolation to get TRNCAL from TRAN_IA and TRAN_IAP1:
        dlt_ia = wavln_std[ia_p1] - wavln_std[ia[j]];
        fia = (wavln_std[ia_p1] - wavobs[j]) / dlt_ia;
        fia_p1 = (float) 1.0 - fia;

        trntblo[j] = fia * tran_ia + fia_p1*tran_iap1;

        if (splitpaths) {

            tg_solo[j] = fia * tran_ias_sol + fia_p1*tran_iap1s_sol;
            tg_seno[j] = fia * tran_ias_sen + fia_p1*tran_iap1s_sen;
            //            printf("TG_SOLO: %d %g %g %g %g %g\n",j,tg_solo[j],tran_ias_sol,tran_iap1s_sol,fia,fia_p1);
            //            printf("TG_SENO: %d %g %g %g %g %g\n",j,tg_seno[j],tran_ias_sen,tran_iap1s_sen,fia,fia_p1);
        }
    }

    free(tran_std_o);
    free(tran_med_o);
    free(tran_std_os_sol);
    free(tran_std_os_sen);
    free(tran_med_index_o);
    free(tran_med_os_sol);
    free(tran_med_os_sen);
    free(tran_med_index_os_sol);
    free(tran_med_index_os_sen);
}

//THE END!!!
