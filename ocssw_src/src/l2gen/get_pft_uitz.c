/*Uitz algorithm for determining the distribution of phytoplankton communities based upon surface chlorophyll
 *
 * This algorithm utilizes the near-surface chlorophyll a concentration as derived from satellite ocean color observation to infer the the phytoplankton functional type composition.
 *
 * Citation:
 * Uitz,J.,H. Claustre, A. Morel, and S.B. Hooker (2006), Vertical distributin of phytoplankton communities in open ocean: An assessment based on surface chlorophyll, J. Geophys. Res., 111, C08005, doi:10.1029/2005JC003207.
 *
 * This algorithm was developed by Robert Lossing  in October 2013 */

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "l12_proto.h"

void calc_pft_uitz(float mld, float lat, float chl, float *fm, float *fn, float *fp) {

    /*Check Chlorophyll and mixed layer depth (mld) inputs and pixel fill values with ERROR CODES (2 if chl<0 and 5 if mld <=0) */

    if (chl < 0) {
        *fm = BAD_FLT;
        *fn = BAD_FLT;
        *fp = BAD_FLT;
        return;
    }
    if (mld <= 0) {
        *fm = BAD_FLT;
        *fn = BAD_FLT;
        *fp = BAD_FLT;
        return;
    }
    float f;
    /*Estimate the euphotic depth (Zeu) using chl Eq. 10 of Morel et al., RSE 111, 69-88, 2007 */

    float x, a, y, zeu;
    x = log(chl);
    a = 0.0186 * pow(x, 3);
    y = 1.524 - 0.436 * x - 0.0145 * pow(x, 2) + a;
    zeu = pow(10, y);

    /*Zm equals the mixed layer depth (mld) */
    float zm;
    zm = mld;

    /*Determine if waters are mixed or stratified. */
    if (zeu / zm > 1) {
        /* Determination of coeffiecients and variables for equation 7.
         *
         * Explanation of variable names:
         * trophcat = Trophic categories defined with respect to the chlorophyll a concentration within the surface layer, [Chla]surf
         *  chpcb,ps,pcmax,pzeta,pzetamax,pdeltamax = parameter placeholder value for equation 7 */

        int trophcathi, trophcatlo;

        /* Values of the five parameters (Table 5) to be used in equation (7) bbtained for the average dimensionless vertical profiles of
         * Chla,Micro-Chla,Nano-Chla, and Pico-Chla, for each trophic class of stratified water (S1 to S9) (Table 5)*/

        // If the variable has a number at the end of it,it is related to a PFT size class.1=chla, 2= micro, 3= nano, 4 = pico. //

        //Table 5 Parameters for Chla
        float cb1 [9] = {0.471, 0.533, 0.428, 0.570, 0.611, 0.390, 0.569, 0.835, 0.188};
        float s1 [9] = {0.135, 0.172, 0.138, 0.173, 0.214, 0.109, 0.183, 0.298, 0};
        float cmax1 [9] = {1.572, 1.194, 1.015, 0.766, 0.676, 0.788, 0.608, 0.382, 0.885};
        float zetamax1 [9] = {0.969, 0.921, 0.905, 0.814, 0.663, 0.521, 0.452, 0.512, 0.378};
        float deltazeta1 [9] = {0.393, 0.435, 0.630, 0.586, 0.539, 0.681, 0.744, 0.625, 1.081};
        //Table 5 Parameters for Micro-Chla
        float cb2[9] = {0.036, 0.071, 0.076, 0.071, 0.145, 0.173, 0.237, 0.331, 0.891};
        float s2 [9] = {0.020, 0.020, 0.021, 0.021, 0.050, 0.044, 0.077, 0.105, 0.302};
        float cmax2 [9] = {0.122, 0.173, 0.126, 0.160, 0.163, 0.161, 0.158, 0.278, 0};
        float zetamax2 [9] = {1.012, 0.885, 0.835, 0.776, 0.700, 0.600, 0.521, 0.451, 0.277};
        float deltazeta2 [9] = {0.532, 0.406, 0.424, 0.546, 0.479, 0.508, 0.543, 0.746, 1.014};
        //Table 5 Parameters for Nano-Chla
        float cb3 [9] = {0.138, 0.129, 0.142, 0.192, 0.188, 0.331, 0.201, 0.227, 0.171};
        float s3 [9] = {0.033, 0.014, 0, 0.037, 0.055, 0.132, 0.084, 0.081, 0};
        float cmax3 [9] = {0.764, 0.589, 0.463, 0.400, 0.418, 0.294, 0.350, 0.198, 0.088};
        float zetamax3 [9] = {0.980, 0.899, 0.872, 0.782, 0.650, 0.501, 0.402, 0.181, 0.375};
        float deltazeta3 [9] = {0.451, 0.454, 0.526, 0.535, 0.640, 0.516, 0724, 0.690, 0.352};
        //Table 5 Parameters for Pico-Chla
        float cb4 [9] = {0.222, 0.242, 0.254, 0.271, 0.159, 0.176, 0.009, 0.094, 0.051};
        float s4 [9] = {0.114, 0.109, 0.099, 0.100, 0.052, 0.071, 0, 0.040, 0.023};
        float cmax4 [9] = {0.906, 0.627, 0.437, 0.255, 0.176, 0.129, 0.251, 0.109, 0};
        float zetamax4 [9] = {0.970, 0.977, 0.969, 0.858, 0.574, 0.458, 0.239, 0.187, 0.052};
        float deltazeta4 [9] = {0.352, 0.427, 0.634, 0.637, 0.650, 0.626, 0.943, 0.618, 0.417};

        //interpolation
        int i;
        float f;

        //for (i = 0; i<9 ;i++)

        /*Trophic Categories for Stratified Waters were Defined With Respect to the Chlorophyll a Concentration Within the Surface Layer,
         * [Chla]surf, and the Associated Parameters. Stratified Waters were separated into 9 trophic categories"trophcat". */

        float avgchlasurf_stratified [9] = {.032, 0.062, 0.098, 0.158, 0.244, 0.347, 0.540, 1.235, 2.953};


        // Average Chla surface mg m-3 values from Table 3 for stratified waters.
        if (chl <= avgchlasurf_stratified[0]) {
            trophcathi = 0;
            trophcatlo = 0;
        } else if (chl > avgchlasurf_stratified[0] && chl <= avgchlasurf_stratified[1]) {
            trophcatlo = 0;
            trophcathi = 1;
        } else if (chl >= avgchlasurf_stratified[1] && chl <= avgchlasurf_stratified[2]) {
            trophcatlo = 1;
            trophcathi = 2;
        } else if (chl >= avgchlasurf_stratified[2] && chl <= avgchlasurf_stratified[3]) {
            trophcatlo = 2;
            trophcathi = 3;
        } else if (chl >= avgchlasurf_stratified[3] && chl <= avgchlasurf_stratified[4]) {
            trophcatlo = 3;
            trophcathi = 4;
        } else if (chl >= avgchlasurf_stratified[4] && chl <= avgchlasurf_stratified[5]) {
            trophcatlo = 4;
            trophcathi = 5;
        } else if (chl >= avgchlasurf_stratified[5] && chl <= avgchlasurf_stratified[6]) {
            trophcatlo = 5;
            trophcathi = 6;
        } else if (chl >= avgchlasurf_stratified[6] && chl <= avgchlasurf_stratified[7]) {
            trophcatlo = 6;
            trophcathi = 7;
        } else if (chl >= avgchlasurf_stratified[7] && chl < avgchlasurf_stratified[8]) {
            trophcatlo = 7;
            trophcathi = 8;
        } else {
            trophcatlo = 8;
            trophcathi = 8;
        }
        // Returns the proper lower and upper Table 5 coeffecients necessary for Equation 7 using the calculated Trophic Category / Stratified Class

        /* weighting scheme uses Kd490 as input into Gordon and Clark 1980. Estimate Kd490 using Chl Eq. 8 of Morel et al., RSE 111, 69-88, 2007. */

        float kd;
        kd = 0.0166 + 0.0773 * pow(chl, 0.6715);


        /* weighting scheme */
        //float od [11] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
        float od; //Optical Depth
        float wt; //Weighting Scheme
        float sumczetachlahi, sumczetamicrohi, sumczetananohi, sumczetapicohi, sumwt;
        float sumczetachlalo, sumczetamicrolo, sumczetananolo, sumczetapicolo;
        float increment;


        sumczetachlahi = 0;
        sumczetamicrohi = 0;
        sumczetananohi = 0;
        sumczetapicohi = 0;
        sumczetachlalo = 0;
        sumczetamicrolo = 0;
        sumczetananolo = 0;
        sumczetapicolo = 0;
        increment = .1;
        sumwt = 0;
        od = 0;

        for (i = 0; i < 11; i++) {

            wt = pow(exp(-kd * zeu * od), 2);
            sumwt += wt;

            /* Equation 7: The upper and lower PFTs fractional values are calulated to provide the upper/lower interpolation boundaries */
            //czeta1[i] = cb1[trophcatlo]-s1[trophcatlo]*od + cmax1[trophcatlo] * exp(-(pow(od-zetamax1[trophcatlo])/deltazeta1[trophcatlo]),2) *wt;
            sumczetachlalo += (cb1[trophcatlo] - s1[trophcatlo] * od + cmax1[trophcatlo] * exp(-pow((od - zetamax1[trophcatlo]) / deltazeta1[trophcatlo], 2))) * wt;
            //Micro -Chla
            sumczetamicrolo += (cb2[trophcatlo] - s2[trophcatlo] * od + cmax2[trophcatlo] * exp(-pow((od - zetamax2[trophcatlo]) / deltazeta2[trophcatlo], 2))) * wt;
            // Nano-Chla
            sumczetananolo += (cb3[trophcatlo] - s3[trophcatlo] * od + cmax3[trophcatlo] * exp(-pow((od - zetamax3[trophcatlo]) / deltazeta3[trophcatlo], 2))) * wt;
            // Pico-Chla
            sumczetapicolo += (cb4[trophcatlo] - s4[trophcatlo] * od + cmax4[trophcatlo] * exp(-pow((od - zetamax4[trophcatlo]) / deltazeta4[trophcatlo], 2))) * wt;

            //czeta1[i] = cb1[trophcathi]-s1[trophcathi]*od + cmax1[trophcathi] * exp(-(pow(od-zetamax1[trophcathi])/deltazeta1[trophcathi]),2) *wt;
            sumczetachlahi += (cb1[trophcathi] - s1[trophcathi] * od + cmax1[trophcathi] * exp(-pow((od - zetamax1[trophcathi]) / deltazeta1[trophcathi], 2))) * wt;
            //Micro -Chla
            sumczetamicrohi += (cb2[trophcathi] - s2[trophcathi] * od + cmax2[trophcathi] * exp(-pow((od - zetamax2[trophcathi]) / deltazeta2[trophcathi], 2))) * wt;
            // Nano-Chla
            sumczetananohi += (cb3[trophcathi] - s3[trophcathi] * od + cmax3[trophcathi] * exp(-pow((od - zetamax3[trophcathi]) / deltazeta3[trophcathi], 2))) * wt;
            // Pico-Chla
            sumczetapicohi += (cb4[trophcathi] - s4[trophcathi] * od + cmax4[trophcathi] * exp(-pow((od - zetamax4[trophcathi]) / deltazeta4[trophcathi], 2))) * wt;

            od += increment;
        }
        //float chla;
        float chlmhi, chlnhi, chlphi;
        float chlmlo, chlnlo, chlplo;
        float chlm, chln, chlp;

        //returned PFT chlorophyll concentration values mg/m-3
        //chla = sumczetachla/sumwt;
        chlmlo = sumczetamicrolo / sumwt;
        chlnlo = sumczetananolo / sumwt;
        chlplo = sumczetapicolo / sumwt;

        chlmhi = sumczetamicrohi / sumwt;
        chlnhi = sumczetananohi / sumwt;
        chlphi = sumczetapicohi / sumwt;

        //Final Interpolation of Calculated Low and High PFT chlorophyll concentration (mg m-3) values

        if (trophcathi == 0 && trophcatlo == 0) {
            f = 1;
        }
        else if (trophcathi == 8 && trophcatlo == 8) {
            f = 1;
        }
        else {
            f = (chl - avgchlasurf_stratified[trophcatlo]) / (avgchlasurf_stratified[trophcathi] - avgchlasurf_stratified[trophcatlo]);
        }
        chlm = ((1 - f) * chlmlo) + (f * chlmhi);
        chln = ((1 - f) * chlnlo) + (f * chlnhi);
        chlp = ((1 - f) * chlplo) + (f * chlphi);

        if (chlm > 1)
            chlm = 1;

        if (chln > 1)
            chln = 1;

        if (chlp > 1)
            chlp = 1;

        if (chlm < 0)
            chlm = 0;

        if (chln < 0)
            chln = 0;

        if (chlp < 0)
            chlp = 0;

        /* fractional contribution of micro, nano, and picoplankton */

        *fm = chlm / (chlm + chln + chlp);
        *fn = chln / (chlm + chln + chlp);
        *fp = chlp / (chlm + chln + chlp);

    }
        /* If Zeu/Zm is NOT greater than 1 then the waters are mixed */
    else {

        if (lat >= -60) {
            /* The following statements determine the corresponding PFTs fractional proportions in GLOBALLY MIXED WATERS (southern mixed waters excluded). */
            int trophcathi, trophcatlo;

            float global_mixed_waters_fmicro [5] = {.180, .241, .281, .522, .909};
            float global_mixed_waters_fnano [5] = {.507, .498, .572, .381, .053};
            float global_mixed_waters_fpico [5] = {.312, .262, .147, .097, .037};

            /* The following procedure determines the high and low fractional values needed to interpolate the average fractional proportion of each PFTs class.
             * The mixed waters into five trophic classes based upon average surface chlorophyll a (mg m-3) */

            float avgchlasurf_global_mixed [5] = {0.234, 0.593, 0.891, 1.540, 7.964};

            if (chl <= avgchlasurf_global_mixed[0]) {
                trophcathi = 0;
                trophcatlo = 0;
            } else if (chl > avgchlasurf_global_mixed[0] && chl <= avgchlasurf_global_mixed[1]) {
                trophcatlo = 0;
                trophcathi = 1;
            } else if (chl > avgchlasurf_global_mixed[1] && chl <= avgchlasurf_global_mixed[2]) {
                trophcatlo = 1;
                trophcathi = 2;
            } else if (chl > avgchlasurf_global_mixed[2] && chl <= avgchlasurf_global_mixed[3]) {
                trophcatlo = 2;
                trophcathi = 3;
            } else if (chl > avgchlasurf_global_mixed[3] && chl <= avgchlasurf_global_mixed[4]) {
                trophcatlo = 3;
                trophcathi = 4;
            } else {
                trophcatlo = 4;
                trophcathi = 4;
            }

            //Final Interpolation of Calculated Low and High PFT chlorophyll concentration (mg m-3) values
            if (trophcathi == 0 && trophcatlo == 0) {
                f = 1;
            }
            else if (trophcathi == 4 && trophcatlo == 4) {
                f = 1;
            }
            else {
                f = (chl - avgchlasurf_global_mixed[trophcatlo]) / (avgchlasurf_global_mixed[trophcathi] - avgchlasurf_global_mixed[trophcatlo]);
            }

            float fmicro, fnano, fpico;

            fmicro = ((1 - f) * global_mixed_waters_fmicro[trophcatlo]) + (f * global_mixed_waters_fmicro[trophcathi]);
            fnano = ((1 - f) * global_mixed_waters_fnano[trophcatlo]) + (f * global_mixed_waters_fnano[trophcathi]);
            fpico = ((1 - f) * global_mixed_waters_fpico[trophcatlo]) + (f * global_mixed_waters_fpico[trophcathi]);

            if (fmicro > 1)
                fmicro = 1;
            if (fnano > 1)
                fnano = 1;
            if (fpico > 1)
                fpico = 1;

            if (fmicro < 0)
                fmicro = 0;
            if (fnano < 0)
                fnano = 0;
            if (fpico < 0)
                fpico = 0;

            /* fractional contribution of micro, nano, and picoplankton */

            *fm = fmicro;
            *fn = fnano;
            *fp = fpico;


        } else {
            /*For southern mixed waters */
            /* The following procedure determines the high and low fractional values needed to interpolate the average fractional proportion of each PFTs class.
             * The mixed waters are separated into five trophic classes based upon average surface chlorophyll a (mg m-3) */
            int trophcathi, trophcatlo;
            float avgchlasurf_southern_mixed [5] = {0.345, 0.605, 0.889, 1.956, 5.755};

            if (chl <= avgchlasurf_southern_mixed[0]) {
                trophcathi = 0;
                trophcatlo = 0;
            } else if (chl > avgchlasurf_southern_mixed[0] && chl <= avgchlasurf_southern_mixed[1]) {
                trophcatlo = 0;
                trophcathi = 1;
            } else if (chl > avgchlasurf_southern_mixed[1] && chl <= avgchlasurf_southern_mixed[2]) {
                trophcatlo = 1;
                trophcathi = 2;
            } else if (chl > avgchlasurf_southern_mixed[2] && chl <= avgchlasurf_southern_mixed[3]) {
                trophcatlo = 2;
                trophcathi = 3;
            } else if (chl > avgchlasurf_southern_mixed[3] && chl <= avgchlasurf_southern_mixed[4]) {
                trophcatlo = 3;
                trophcathi = 4;
            } else {
                trophcatlo = 4;
                trophcathi = 4;
            }

            //Final Interpolation of Calculated Low and High PFT chlorophyll concentration (mg m-3) values


            float southern_mixed_waters_fmicro [5] = {.534, .507, .446, .409, .137};
            float southern_mixed_waters_fnano [5] = {.360, .441, .507, .566, .853};
            float southern_mixed_waters_fpico [5] = {.106, .052, .047, .025, .010};

            float f;

            if (trophcathi == 0 && trophcatlo == 0) {
                f = 1;
            }
            else if (trophcathi == 4 && trophcatlo == 4) {
                f = 1;
            }
            else {
                f = (chl - avgchlasurf_southern_mixed[trophcatlo]) / (avgchlasurf_southern_mixed[trophcathi] - avgchlasurf_southern_mixed[trophcatlo]);
            }

            float fmicro, fnano, fpico;

            fmicro = ((1 - f) * southern_mixed_waters_fmicro[trophcatlo]) + (f * southern_mixed_waters_fmicro[trophcathi]);
            fnano = ((1 - f) * southern_mixed_waters_fnano[trophcatlo]) + (f * southern_mixed_waters_fnano[trophcathi]);
            fpico = ((1 - f) * southern_mixed_waters_fpico[trophcatlo]) + (f * southern_mixed_waters_fpico[trophcathi]);

            if (fmicro > 1)
                fmicro = 1;
            if (fnano > 1)
                fnano = 1;
            if (fpico > 1)
                fpico = 1;

            if (fmicro < 0)
                fmicro = 0;
            if (fnano < 0)
                fnano = 0;
            if (fpico < 0)
                fpico = 0;

            /* fractional contribution of micro, nano, and picoplankton */

            *fm = fmicro;
            *fn = fnano;
            *fp = fpico;
        }
    }
}

void get_pft_uitz(l2str *l2rec, l2prodstr *p, float prod[]) {
    int i;
    float mld, fm, fn, fp;
    int16_t year, day;
    double sec;
    unix2yds(l2rec->l1rec->scantime, &year, &day, &sec);

    l1str *l1rec = l2rec->l1rec;

    for (i = 0; i < l1rec->npix; i++) {
        mld = get_mld(input->mldfile, l1rec->lon[i], l1rec->lat[i], day);
        if (l2rec->chl[i] == BAD_FLT) {
            prod[i] = BAD_FLT;
            continue;
        }

        switch (p->cat_ix) {
        case CAT_microplankton_uitz:
            calc_pft_uitz(mld, l1rec->lat[i], l2rec->chl[i], &prod[i], &fn, &fp);
            break;
        case CAT_nanoplankton_uitz:
            calc_pft_uitz(mld, l1rec->lat[i], l2rec->chl[i], &fm, &prod[i], &fp);
            break;
        case CAT_picoplankton_uitz:
            calc_pft_uitz(mld, l1rec->lat[i], l2rec->chl[i], &fm, &fn, &prod[i]);
            break;

        default:
            printf("calc_chla_uitz can not produce product %s\n", p->algorithm_id);
            exit(1);
        }
        if (isnan(prod[i])) {
            prod[i] = BAD_FLT;
            l1rec->flags[i] |= PRODFAIL;
        }
    }
}

