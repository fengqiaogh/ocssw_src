/*
 * get_atrem_cor2.c
 *
 *  Created on: Feb 19, 2015
 *      Author: Rick Healy (richard.healy@nasa.gov)
 *
 *      *  Notes about water vapor VMRS and related quantities:                *
 *                                          *
 *     VAPVRT(60)   - a table containing 60 column vapor values (in unit of cm) *
 *                                          *
 *     VAP_SLANT(I) = VAPVRT(I) * 2.0, VAP_SLANT is a new table for containing  *
 *                    two-way total vapor amounts. Here the number "2" can be   *
 *                    changed to other numbers, e.g., 2.5, without major        *
 *                    effects on retrieved water vapor values.                  *
 *                                          *
 *     G_VAP(I = 1,..., NL) = true vapor geometric factor for each layer in     *
 *                    the model atmosphere (after adjusting for the elevated    *
 *                    surface.                                                  *
 *                                          *
 *     VMRM(I) = VMRM(I)*G_VAP(I). The VMRS are multiplied by the geometrical   *
 *                    factor. We can calculate the vapor transmittance on the   *
 *                    Sun-surface-sensor path by assuming a vertical path in    *
 *                    the model atmosphere with geometric-factor-adjusted VMRS. *
 *                                          *
 *     CLMVAP  = vertical column amount from ground to space in model atmosphere*
 *     CLMVAPP = vertical column amount from ground to aircraft or satellite    *
 *                    sensor in model atmosphere                                *
 *     Q       = 2.152E25 = # of molecules above the surface at one atmosphere  *
 *                    (in unit of  molecules/cm**2)                 *
 *                                          *
 *     VAP_SLANT_MDL= CLMVAP/COS(SOLZNI) + CLMVAPP/COS(OBSZNI) = total amount   *
 *                    of water vapor in the model atmosphere in the L-shaped    *
 *                    Sun-surface-plane ray path.                               *
 *                                          *
 *     G_VAP_EQUIV  = VAP_SLANT_MDL / CLMVAP = the "equivalent" geometrical     *
 *                    factor corresponding to the total slant vapor amount      *
 *                    VAP_SLANT_MDL and the column vapor amount CLMVAP.         *
 *                                          *
 *     SSH2O(I) (I = 1, ..., 60) - a pure scaling factor relative to the total  *
 *                    slant vapor amount of VAP_SLANT_MDL, and                  *
 *            SSH2O(I) = VAP_SLANT(I) / VAP_SLANT_MDL               *
 *                                          *
 *     SH2O = one value of SSH2O(I). SH2O is used during generation of the      *
 *            look-up table.                            *
 *                                          *
 *     VAPTT  = VAP_SLANT_MDL*SH2O, is the absolute total vapor amount on the   *
 *                    L-shaped path corresponding to a spectrum stored in the   *
 *                    look-up table.                                        *
 *                                          *
 *     CLMWVP = 0.5*(VAPTTA+VAPTTB)/G_VAP_EQUIV, is the retrieved column water  *
 *                    vapor amount from imaging spectrometer data.          *
 ********************************************************************************
 *
 *
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
 *                                  !*  References:                                        *
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
 *
 */

/* compile this code with:
 gcc -O1 -Wpadded -Wpacked -malign-double -mpreferred-stack-boundary=8 -o get_atrem_cor3   \
 get_atrem_cor3.c rdatreminfo.c atrem_app_refl_plus_gas_removal_for_l2gen3.o cubeio.o     \
 tpvmr_init.o solar_irr_PC.o bndprms.o -lgfortran \
 -L/disk01/home/rhealy/ocssw/build/cbuild/src/libgenutils -lgenutils \
 -L/disk01/home/rhealy/ocssw/build/lib3/src/netcdf/netcdf-fortran-4.2/fortran/.libs \
 -lnetcdff -L/disk01/home/rhealy/ocssw/build/lib3/src/netcdf/netcdf-4.3.1.1/liblib/.libs -lnetcdf \
 -L/disk01/home/rhealy/ocssw/build/lib3/src/hdf5/hdf5-1.8.10-patch1/src/.libs -lhdf5 \
 -L/disk01/home/rhealy/ocssw/build/lib3/src/hdf5/hdf5-1.8.10-patch1/hl/src/.libs -lhdf5_hl \
 -L/disk01/home/rhealy/ocssw/build/lib3/src/netcdf/netcdf-4.3.1.1/libsrc4/.libs -lnetcdf4  -lz \
 -L/disk01/home/rhealy/ocssw/build/lib3/src/hdf5/hdf5-1.8.10-patch1/src/.libs -lhdf5

 *   compile fortran routines:
    gfortran -malign-double -c atrem_app_refl_plus_gas_removal_for_l2gen2.f90 cubeio.f90 \
              bndprms.f solar_irr_PC.f tpvmr_init.f

---- Previous version without netcdf -----

  gcc -Wpadded -Wpacked -malign-double -mpreferred-stack-boundary=8 -o get_atrem_cor \
    get_atrem_cor.c atrem_app_refl_plus_gas_removal_for_l2gen.o cubeio.o \
    tpvmr_init.o solar_irr_PC.o bndprms.o -lgfortran -L/disk01/home/rhealy/ocssw/build/cbuild/src/libgenutils -lgenutils
 *
 *   compile fortran routines:
    gfortran -malign-double  -I/disk01/home/rhealy/ocssw/build/lib3/src/netcdf/netcdf-fortran-4.2/fortran -c atrem_app_refl_plus_gas_removal_for_l2gen3.f90 cubeio.f90 \
              bndprms.f solar_irr_PC.f tpvmr_init.f

              Example:
              get_atrem_cor2 < input/test_input.txt  > cor2_test90b.out
 */
#include "atrem_corl1.h"
#include <sensorDefs.h>
#include <timeutils.h>
#include <math.h>
#define BUFSZ 100
#define MINDEGCHANGE 55 // Minimum change (degrees squared) in sum of zenith and azimuth angle squared
// for both solar and sensor angles before recalculating tran_table (transmittance table)

int get_atrem_cor(int32_t sensorID, l1str *l1rec, int32_t ip, float *rhot, float *tg_tot, float *tg_sol, float *tg_sen) {

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
            printf("Reading geometry limiting angles for Atrem\n");
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

        geometry_();
        //printf("Processed geometry after %6.0f seconds\n",now()-start_time);

        init_speccal_();
        //printf("Processed init_speccal after %6.0f seconds\n",now()-start_time);

        tran_table_();
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


        //            printf("nb: %d %d %d %d\n",P.nb1,P.nb2,P.nb3,P.nb4);
        //            printf("nb: %d %d\n",P.nbp94,P.nb1p14);
        //            printf("nh2o: %d\n",P.nh2o);
        //            printf("nobs: %d\n",P.nobs);
        //            printf("startndx: %d %d %d %d\n",P.start_ndx[0],P.start_ndx[1],P.start_ndx[2],P.start_ndx[3]);
        //            printf("endndx: %d %d %d %d\n",P.end_ndx[0],P.end_ndx[1],P.end_ndx[2],P.end_ndx[3]);
        //            printf("start_p94: %d\n",P.start_p94);
        //            printf("end_p94: %d\n",P.end_p94);
        //            printf("start_1p14: %d\n",P.start_1p14);
        //            printf("end_1p14: %d\n",P.end_1p14);
        //            printf("%d\n",P.start2);
        //            printf("%d\n",P.end2);
        //
        //            printf("ncv2: %d\n",P.ncv2);
        //            printf("ntot: %d %d %d %d\n",P.natot,P.nbtot,P.nctot,P.ndtot);
        //            printf("wt: %f %f %f %f\n",P.wt1,P.wt2,P.wt3,P.wt4);
        //            printf("delta: %f %f\n",P.delta,P.delta2);
        //              printf("g_vap_equiv: %f \n",P.g_vap_equiv);
        //            printf("ip=%d\n",ip);
        //            printf("prev_dist=%f\n",prev_ddeg);
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
            nb = init_tpvmr(P.model);
            if (nb <= 0) {
                printf("-E- %s line %d : Atmospheric data could not be initialized.\n",
                        __FILE__, __LINE__);
                exit(FATAL_ERROR);
            }
            model_adj_();
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
        geometry_();

        // Don't need to call init_speccal since other trace gas transmittances don't depend on water vapor

        // recalculate transmittances based on separate solar/sensor paths using the tran_table indices from the first call
        tran_table_();

        for (i = 0; i < P.nobs; i++) {
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

    double delta, deltab, fja, fjap1, fjb, fjbp1, vaptta, vapttb;
    double speca[NBANDS], specb[NBANDS], specav, spec450;
    ja = P->nh2o / 2;
    ja = hunt(P->r0p94, P->nh2o, ratio_94c, ja);
    //      printf("xx[%d]=%f xx[%d]=%f JA=%d\n",ja-1,P->r0p94[ja-1],ja,P->r0p94[ja],ja);
    if (ja >= 0 && ja < NH2OMAX) {
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
        for (i = 0; i < P->nobs; i++) {
            if (ja >= 0 && ja < NH2OMAXM1) {
                speca[i] = fja * P->trntbl[ja][i] + fjap1 * P->trntbl[ja + 1][i];
            } else {
                if (ja < 0) speca[i] = P->trntbl[ja + 1][i];
                if (ja >= NH2OMAXM1) speca[i] = P->trntbl[ja][i];
            }
        }
    }

    if (ratio_94co > 1.0) {
        for (i = 0; i < P->nobs; i++) {
            if (ja >= 0 && ja < NH2OMAXM1) {
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

    if (jb >= 0 && jb < NH2OMAXM1) {
        deltab = P->r1p14[jb + 1] - P->r1p14[jb];
        fjb = (P->r1p14[jb + 1] - ratio_114c) / deltab;
        fjbp1 = (ratio_114c - P->r1p14[jb]) / deltab;
        vapttb = fjb * P->vaptot[jb] + fjbp1 * P->vaptot[jb + 1];
        if (ratio_114co > 1.0) vapttb = -vapttb;
    } else {
        if (jb < 0) vapttb = P->vaptot[jb + 1];
        if (jb <= NH2OMAXM1) vapttb = P->vaptot[jb];
    }

    if (ratio_114co <= 1.0) {
        for (i = 0; i < P->nobs; i++) {
            if (jb >= 0 && jb < NH2OMAXM1) {
                specb[i] = fjb * P->trntbl[jb][i] + fjbp1 * P->trntbl[jb + 1][i];
            } else {
                if (jb < 0) specb[i] = P->trntbl[jb + 1][i];
                if (jb >= NH2OMAXM1) specb[i] = P->trntbl[jb][i];
            }
        }
    }
    if (ratio_114co > 1.0) {
        for (i = 0; i < P->nobs; i++) {
            if (jb >= 0 && jb < NH2OMAXM1) {
                specb[i] = 1.0 / (fjb * P->trntbl[jb][i] + fjbp1 * P->trntbl[jb + 1][i]);
            } else {
                if (jb < 0) specb[i] = 1.0 / P->trntbl[jb + 1][i];
                if (jb >= NH2OMAXM1) specb[i] = 1.0 / P->trntbl[jb][i];
            }
        }
    }

    clmwvp = 0.5 * (vaptta + vapttb) / P->g_vap_equiv;
    spec450 = 1.0 - 0.5 * (speca[P->idx450] + specb[P->idx450]); // should be 0.0

    //  Derivation of surface reflectances

    for (i = 0; i < P->nobs; i++) {
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
    printf("RJH: Initializing ATM for model number = %d\n", model);
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

    model = (int *) calloc(1, sizeof (model));
    fwhm = (float *) calloc(1, sizeof (fwhm));
    xppp = (float *) calloc(1, sizeof (xppp));

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
    get_input_();

    strcpy(input_l2gen_.filename, filename);
    input_l2gen_.dln = strlen(filename);
    //printf("filedir=%s Filename=%s len=%d\n",filedir, filename,input_l2gen_.dln);

    if (!(fwhm = l1rec->l1file->fwhm)) {
        nwave = rdatreminfo(sensorID, 0, "fwhm", (void **) &fwhm);
    }
    for (i = 0; i < nbands; i++) {
        getinput4_.wavobs[i] = l1rec->l1file->fwave[i] / 1000.; // Convert nm to microns
        getinput4_.fwhm[i] = fwhm[i];
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
    nb = init_tpvmr(P->model);

    if (nb <= 0) {
        printf("-E- %s line %d : Atmospheric data could not be initialized.\n",
                __FILE__, __LINE__);
        exit(FATAL_ERROR);
    }

    getinput5_.nobs = nbands; //l1rec->nbands;
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
    float *wndow1 = 0, *wndow2 = 0, *wp94c = 0, *wndow3 = 0, *wndow4 = 0, *w1p14c = 0;

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

    if (atrem_splitpaths > 0 && atrem_full == 0) {
        printf("ATREM: Turning on atrem_full because you selected atrem_splitpaths\n");
        atrem_full = 1;
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

    printf("ATREM Gas Options:H2O:1 CO2:%d O3:%d N2O:%d CO:%d CH4:%d O2:%d NO2:%d\n",
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

    model_adj_();

    P->nb1 = getinput7_.nb1;
    P->nb2 = getinput7_.nb2;
    P->nb3 = getinput7_.nb3;
    P->nb4 = getinput7_.nb4;
    P->nbp94 = getinput7_.nbp94;
    P->nb1p14 = getinput7_.nb1p14;
    P->nobs = getinput5_.nobs;
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
    int ncid, ndims;
    int32 sds_id;
    int status;
    nc_type rh_type; /* variable type */
    int dimids[H4_MAX_VAR_DIMS]; /* dimension IDs */
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

    printf("ATREM Angle_Limit FILE=%s\n", filename);

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

int findMatch_(float *list, int *nobs, float *elem) {
    int i;
    float diff1, diff2;

    for (i = 0; i<*nobs && list[i] <= *elem; i++)

        diff1 = fabs(list[i - 1] - *elem);
    diff2 = fabs(list[i] - *elem);

    if (diff1 < diff2)
        return (i - 1);
    else
        return (i);

}
