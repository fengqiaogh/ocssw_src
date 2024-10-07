/*
 * get_atrem_cor2.c
 *
 *  Created on: Feb 19, 2015
 *      Author: rhealy
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
 */

/* compile this code with:
 gcc -O1 -Wpadded -Wpacked -malign-double -mpreferred-stack-boundary=8 -o get_atrem_cor3   \
 get_atrem_cor3.c atrem_app_refl_plus_gas_removal_for_l2gen3.o cubeio.o     \
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
#include "atrem_cor.h"
#define BUFSZ 100

int main(int argc, char *argv[]) {

    FILE *fp;
    char *fname;
    char line[BUFSZ];
    paramstr P;
    double *YY;
    int cnt = 0, i, j, k, l;
    double start_time, tot_time = 0;

    //   char *f_rad = "atrem_yy_pix2line3_before.txt";
    char *f_rad = "atrem_yy_pix512line2000_before.txt";
    if (argc == 2) {
        f_rad = argv[1];
        //    }else{
        //        fprintf(stderr,"%s: Please supply a file name\n",argv[0]);
        //        exit(1);
    }
    //    P.r0p94 = (float *) calloc(TBLMAX,sizeof(float));
    //    P.r1p14 = (float *) calloc(TBLMAX,sizeof(float));
    //    P.finst2 = (float *) calloc(FINSTMAX,sizeof(float));
    //    P.vaptot = (float *) calloc(TBLMAX,sizeof(float));
    //    P.trntbl = (float *) calloc(TBLMAX,sizeof(t_array));

    /* call fortran routines and transfer parameters over to C struct */
    start_time = now();
    printf("\nBegin get_atrem_cor processing at %s\n\n", ydhmsf(start_time, 'L'));

    get_input_();
    printf("Processed get_input after %6.0f seconds\n", now() - start_time);
    model_adj_();
    printf("Processed model_adj after %6.0f seconds\n", now() - start_time);
    geometry_();
    printf("Processed geometry after %6.0f seconds\n", now() - start_time);
    init_speccal_();
    printf("Processed init_speccal after %6.0f seconds\n", now() - start_time);
    solar_irr_pc_();
    printf("Processed solar_irr_pc after %6.0f seconds\n", now() - start_time);
    tran_table_();
    printf("Processed tran_table after %6.0f seconds\n", now() - start_time);

    P.nb1 = getinput7_.nb1;
    P.nb2 = getinput7_.nb2;
    P.nb3 = getinput7_.nb3;
    P.nb4 = getinput7_.nb4;
    P.nbp94 = getinput7_.nbp94;
    P.nb1p14 = getinput7_.nb1p14;
    P.nh2o = init_speccal3_.nh2o;
    P.nobs = getinput5_.nobs;
    P.delta = getinput5_.dlt;
    P.delta2 = getinput5_.dlt2;
    P.start_ndx[0] = init_speccal6_.ist1 - 1;
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
    P.start2 = init_speccal10_.istrt2 - 1;
    P.end2 = init_speccal10_.iend2 - 1;
    P.ncv2 = init_speccal10_.ncv2;
    P.natot = init_speccal11_.natot;
    P.nbtot = init_speccal11_.nbtot;
    P.nctot = init_speccal11_.nctot;
    P.ndtot = init_speccal11_.ndtot;
    P.wt1 = init_speccal8_.wt1;
    P.wt2 = init_speccal8_.wt2;
    P.wt3 = init_speccal8_.wt3;
    P.wt4 = init_speccal8_.wt4;
    P.g_vap_equiv = geometry3_.g_vap;
    P.r0p94 = tran_table1_.r0p94;
    P.r1p14 = tran_table1_.r1p14;
    P.vaptot = tran_table1_.vaptot;
    P.trntbl = tran_table1_.trntbl;
    P.finst2 = init_speccal10_.finst2;

    printf("nb: %d %d %d %d\n", P.nb1, P.nb2, P.nb3, P.nb4);
    printf("nb: %d %d\n", P.nbp94, P.nb1p14);
    printf("nh2o: %d\n", P.nh2o);
    printf("nobs: %d\n", P.nobs);
    printf("startndx: %d %d %d %d\n", P.start_ndx[0], P.start_ndx[1], P.start_ndx[2], P.start_ndx[3]);
    printf("endndx: %d %d %d %d\n", P.end_ndx[0], P.end_ndx[1], P.end_ndx[2], P.end_ndx[3]);
    printf("start_p94: %d\n", P.start_p94);
    printf("end_p94: %d\n", P.end_p94);
    printf("start_1p14: %d\n", P.start_1p14);
    printf("end_1p14: %d\n", P.end_1p14);
    printf("%d\n", P.start2);
    printf("%d\n", P.end2);

    printf("ncv2: %d\n", P.ncv2);
    printf("ntot: %d %d %d %d\n", P.natot, P.nbtot, P.nctot, P.ndtot);
    printf("wt: %f %f %f %f\n", P.wt1, P.wt2, P.wt3, P.wt4);
    printf("delta: %f %f\n", P.delta, P.delta2);
    printf("g_vap: %f \n", P.g_vap_equiv);



    double *YYp, cst1, cst2, cst3, cst4, cst5, cst6, rp94, r1p14;
    int ii, ll, kk;

    YY = (double *) calloc(P.nobs, sizeof (double));
    YYp = (double *) calloc(P.nobs, sizeof (double));
    if ((fp = fopen(f_rad, "r")) == NULL) {
        fprintf(stderr, "Unable to open %s\n", f_rad);
        exit(1);
    }
    int32_t ja, jb;

    for (l = 1999; l < 2000; l++) {
        for (k = 511; k < 512; k++) {

            for (i = 0; i < P.nobs && fgets(line, BUFSZ, fp); i++) {
                sscanf(line, "%d %d %d %lf", &ii, &kk, &ll, &YY[i]);
                YYp[i] = YY[i];
            }
            start_time = now();
            cnt = get_atrem(YY, P, &rp94, &r1p14, &cst1, &cst2, &cst3, &cst4, &cst5, &cst6, &ja, &jb);
            tot_time += now() - start_time;
            for (i = 0; i < P.nobs; i++) {
                printf("%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d \n", i + 1, kk, ll, YYp[i], YY[i], rp94, r1p14, cst1, cst2, cst3, cst4, cst5, cst6, ja, jb);
            }
        }

    }
    printf("Processed %d get_atrem calls after %6.0f seconds\n", 2000 * 512, tot_time);
    return (0);
}

int get_atrem(double *yy, paramstr P, double *rp94, double *r1p14, double *cst1, double *cst2, double *cst3, double *cst4, double *cst5, double *cst6, int32_t *jac, int32_t *jbc) {

    double const1 = 0, const2 = 0, const3 = 0, const4 = 0, const5 = 0, const6 = 0;
    double ratio_94c, ratio_94co, ratio_114c, ratio_114co;
    int32_t i, ja, jb;
    double clmwvp;
    /*
     Arrays related to look-up table:
         VAPTOT: TOTAL SUN-SURFACE-SENSOR PATH WATER VAPOR IN UNITS OF CM
         R0P94 : Tc(0.94 um)/(WT1*Tc(0.86)+WT2*Tc(1.02))
         R1P14 : Tc(1.14 um)/(WT3*Tc(1.05)+WT4*Tc(1.23))

     Calculating 3-channel ratios from an observed spectrum, using a
     look-up table procedure to derive the column amount of water vapor
     and to find the associated simulated transmittance spectrum.
     */
    for (i = P.start_ndx[0]; i <= P.end_ndx[0]; i++) {
        const1 += yy[i];
    }
    const1 /= P.nb1;

    for (i = P.start_ndx[1]; i <= P.end_ndx[1]; i++) {
        const2 += yy[i];
    }
    const2 /= P.nb2;

    //      printf("const2=%f nb2=%d istr2=%d ind2=%d\n",const2,P.nb2,P.start_ndx[1],P.end_ndx[1]);

    for (i = P.start_p94; i <= P.end_p94; i++) {
        const3 += yy[i];
    }
    const3 /= P.nbp94;

    ratio_94co = const3 / ((P.wt1 * const1) + (P.wt2 * const2));
    ratio_94c = ratio_94co;

    if (ratio_94co > 1.0) {
        const1 = 0.0;

        for (i = P.start_ndx[0]; i <= P.end_ndx[0]; i++) {
            const1 += (1.0 / yy[i]);
        }
        const1 /= P.nb1;

        const2 = 0.0;
        for (i = P.start_ndx[1]; i <= P.end_ndx[1]; i++) {
            const2 += (1.0 / yy[i]);
        }
        const2 /= P.nb2;
        const3 = 0.0;
        for (i = P.start_p94; i <= P.end_p94; i++) {
            const3 += (1.0 / yy[i]);
        }
        const3 /= P.nbp94;

        ratio_94c = const3 / ((P.wt1 * const1) + (P.wt2 * const2));
    }

    *rp94 = ratio_94c;

    const4 = 0.0;
    for (i = P.start_ndx[2]; i <= P.end_ndx[2]; i++) {
        const4 += yy[i];
    }
    const4 /= P.nb3;

    const5 = 0.0;
    for (i = P.start_ndx[3]; i <= P.end_ndx[3]; i++) {
        const5 += yy[i];
    }
    const5 /= P.nb4;

    const6 = 0.0;
    for (i = P.start_1p14; i <= P.end_1p14; i++) {
        const6 += yy[i];
    }
    const6 /= P.nb1p14;

    /* DEBUG
     *
     */
    *cst1 = const1;
    *cst2 = const2;
    *cst3 = const3;
    *cst4 = const4;
    *cst5 = const5;
    *cst6 = const6;

    ratio_114co = const6 / ((P.wt3 * const4) + (P.wt4 * const5));
    ratio_114c = ratio_114co;

    if (ratio_114co > 1.0) {

        const4 = 0.0;
        for (i = P.start_ndx[2]; i <= P.end_ndx[2]; i++) {
            const4 += (1.0 / yy[i]);
        }
        const4 /= P.nb3;
        for (i = P.start_ndx[3]; i <= P.end_ndx[3]; i++) {
            const5 += (1.0 / yy[i]);
        }
        const5 /= P.nb4;
        const6 = 0.0;
        for (i = P.start_1p14; i <= P.end_1p14; i++) {
            const6 += (1.0 / yy[i]);
        }
        const6 /= P.nb1p14;
        ratio_114c = const6 / ((P.wt3 * const4) + (P.wt4 * const5));
    }

    *r1p14 = ratio_114c;

    double delta, deltab, fja, fjap1, fjb, fjbp1, vaptta, vapttb;
    double speca[NBANDS], specb[NBANDS], specav;
    ja = P.nh2o / 2;
    ja = hunt(P.r0p94, P.nh2o, ratio_94c, ja);
    if (ja >= 0 && ja < TBLMAX) {
        delta = P.r0p94[ja + 1] - P.r0p94[ja];
        fja = (P.r0p94[ja + 1] - ratio_94c) / delta;
        fjap1 = (ratio_94c - P.r0p94[ja]) / delta;
        vaptta = fja * P.vaptot[ja] + fjap1 * P.vaptot[ja + 1];
        if (ratio_94co > 1.0) vaptta = -vaptta;
    } else {
        if (ja < 0) vaptta = P.vaptot[ja + 1];
        if (ja > P.nh2o) vaptta = P.vaptot[ja];
    }

    if (ratio_94co <= 1.0) {
        for (i = 0; i < P.nobs; i++) {
            if (ja >= 0 && ja < TBLMAXM1) {
                speca[i] = fja * P.trntbl[ja][i] + fjap1 * P.trntbl[ja + 1][i];
            } else {
                if (ja < 0) speca[i] = P.trntbl[ja + 1][i];
                if (ja >= TBLMAXM1) speca[i] = P.trntbl[ja][i];
            }
        }
    }

    if (ratio_94co > 1.0) {
        for (i = 0; i < P.nobs; i++) {
            if (ja >= 0 && ja < TBLMAXM1) {
                speca[i] = 1.0 / (fja * P.trntbl[ja][i] + fjap1 * P.trntbl[ja + 1][i]);
            } else {
                if (ja < 0) speca[i] = 1.0 / P.trntbl[ja + 1][i];
                if (ja >= TBLMAXM1) speca[i] = 1.0 / P.trntbl[ja][i];
            }
        }
    }

    jb = ja;

    jb = hunt(&P.r1p14[0], P.nh2o, ratio_114c, jb);
    *jac = ja;
    *jbc = jb;

    if (jb >= 0 && jb < TBLMAXM1) {
        deltab = P.r1p14[jb + 1] - P.r1p14[jb];
        fjb = (P.r1p14[jb + 1] - ratio_114c) / deltab;
        fjbp1 = (ratio_114c - P.r1p14[jb]) / deltab;
        vapttb = fjb * P.vaptot[jb] + fjbp1 * P.vaptot[jb + 1];
        if (ratio_114co > 1.0) vapttb = -vapttb;
    } else {
        if (jb < 0) vapttb = P.vaptot[jb + 1];
        if (jb <= TBLMAXM1) vapttb = P.vaptot[jb];
    }

    if (ratio_114co <= 1.0) {
        for (i = 0; i < P.nobs; i++) {
            if (jb >= 0 && jb < TBLMAXM1) {
                specb[i] = fjb * P.trntbl[jb][i] + fjbp1 * P.trntbl[jb + 1][i];
            } else {
                if (jb < 0) specb[i] = P.trntbl[jb + 1][i];
                if (jb >= TBLMAXM1) specb[i] = P.trntbl[jb][i];
            }
        }
    }
    if (ratio_114co > 1.0) {
        for (i = 0; i < P.nobs; i++) {
            if (jb >= 0 && jb < TBLMAXM1) {
                specb[i] = 1.0 / (fjb * P.trntbl[jb][i] + fjbp1 * P.trntbl[jb + 1][i]);
            } else {
                if (jb < 0) specb[i] = 1.0 / P.trntbl[jb + 1][i];
                if (jb >= TBLMAXM1) specb[i] = 1.0 / P.trntbl[jb][i];
            }
        }
    }

    clmwvp = 0.5 * (vaptta + vapttb) / P.g_vap_equiv;


    //  Derivation of surface reflectances

    for (i = 0; i < P.nobs; i++) {
        specav = 0.5 * (speca[i] + specb[i]);
        yy[i] /= specav;
    }




    // Smooth the derived surface reflectance spectra

    int ii, j;
    double truncv;

    if (P.delta2 > P.delta) {
        /*
         * First, replace radiances <= 0 near the spectral overlapping parts of the
         * four AVIRIS spectrometers by radiances in the nearby AVIRIS' channels.
         */
        for (i = P.natot - 3; i < P.natot + 2; i++) {
            if (yy[i] <= 0.0) yy[i] = yy[i - 1];
        }
        for (i = P.nbtot - 3; i < P.nbtot + 2; i++) {
            if (yy[i] <= 0.0) yy[i] = yy[i - 1];
        }
        for (i = P.nctot - 3; i < P.nctot + 2; i++) {
            if (yy[i] <= 0.0) yy[i] = yy[i - 1];
        }
        for (i = P.ndtot - 3; i < P.ndtot + 2; i++) {
            if (yy[i] <= 0.0) yy[i] = yy[i - 1];
        }
        for (i = P.start2 - 1; i < P.end2; i++) {
            truncv = 0.0;
            ii = i - P.ncv2 - 1;
            for (j = ii; j < i + P.ncv2; j++) {
                truncv += yy[j] * P.finst2[j - ii + 1];
            }
            yy[i] = truncv;
        }
    }


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

    if (jlo >= 0 && jlo < n) {

        inc = 1;
        if ((x >= xx[jlo]) == ascnd) {
            jhi = jlo + inc;
            while (jhi < n && (x >= xx[jhi]) == ascnd) {
                if (jhi >= n) {
                    jhi = n;
                    break;
                } else if ((x >= xx[jhi]) == ascnd) {
                    jlo = jhi;
                    inc += inc;
                    jhi = jlo + inc;
                }
            }

        } else {
            jhi = jlo;
            jlo = jhi - inc;
            while (jlo >= 0 && (x < xx[jlo]) == ascnd) {
                if (jlo < 0) {
                    jlo = -1;
                    break;
                } else if ((x < xx[jlo]) == ascnd) {
                    jhi = jlo;
                    inc += inc;
                    jlo = jhi - inc;
                }
            }

        }

    } else {
        jlo = -1;
        jhi = n;
    }

    while (jhi - jlo != 1) {

        jm = (jhi + jlo) / 2;

        if ((x == xx[jm]) == ascnd) {
            jlo = jm;
        } else {
            jhi = jm;
        }
    }

    if (x == xx[n - 1]) jlo = n - 2;
    if (x == xx[0]) jlo = 0;

    return (jlo);
}

