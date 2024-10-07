#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"

#define MINWVDIFF 10
#define REFWVL 443

typedef struct _context {
    float *wvl_i;
    float *tar_rrs;
    float *aw_i, *bbw_i, *aBric_i, *bBric_i;
    float *aw_o, *bbw_o, *aBric_o, *bBric_o;
    float *aphRef, *adgRef, *bbpRef;
    int *sh_strt_idx, *num_wvl_i;
} ccstr;

float conv_rrs_to_555(float Rrs, float wave, float uRrs_in, float *uRrs_out ) {
    float sw, a1, b1, a2, b2;
    float temp_Rrs;

    if (fabs(wave - 555) > 2) {
        if (fabs(wave - 547) <= 2) {
            sw = 0.001723;
            a1 = 0.986;
            b1 = 0.081495;
            a2 = 1.031;
            b2 = 0.000216;
        } else if (fabs(wave - 550) <= 2) {
            sw = 0.001597;
            a1 = 0.988;
            b1 = 0.062195;
            a2 = 1.014;
            b2 = 0.000128;
        } else if (fabs(wave - 560) <= 2) {
            sw = 0.001148;
            a1 = 1.023;
            b1 = -0.103624;
            a2 = 0.979;
            b2 = -0.000121;
        } else if (fabs(wave - 565) <= 2) {
            sw = 0.000891;
            a1 = 1.039;
            b1 = -0.183044;
            a2 = 0.971;
            b2 = -0.000170;
        } else {
            printf("-E- %s line %d: Unable to convert Rrs at %f to 555nm.\n", __FILE__, __LINE__, wave);
            exit(1);
        }

        if (Rrs < sw){
        	temp_Rrs=Rrs;
        	Rrs = pow(10.0, a1 * log10(Rrs) - b1);
        	if(uRrs_out) {
        	    *uRrs_out=Rrs*a1/temp_Rrs*uRrs_in;
            }
        } else{
        	Rrs = a2 * Rrs - b2;
        	if(uRrs_out) {
        	    *uRrs_out=a2*uRrs_in;
            }
        }
    }

    return (Rrs);
}

/*
int windex(float wave, float twave[], int ntwave)
{
    int   iw, index;
    float wdiff;
    float wdiffmin = 99999.;

    for (iw=0; iw<ntwave; iw++) {

        // break on exact match 
        if (twave[iw] == wave) {
            index = iw;
            break;
        }

        // look for closest 
        wdiff = fabs(twave[iw]-wave);
        if (wdiff < wdiffmin) {
            wdiffmin = wdiff;
            index = iw;
        }
    }

    return(index);
}
 */
void grabBricaud(float wvl[], int nwvl, float* a_bricaud, float* b_bricaud) {

    float bricaud_L_A_B[][17] ={
        {410, 412, 413, 440, 443, 488, 490, 510, 530, 531, 547, 550, 555, 560, 665, 667, 670},
        {0.0313, 0.0323, 0.032775, 0.0403, 0.0394, 0.0279, 0.0274, 0.018, 0.0117, 0.0115, 0.00845,
            0.008, 0.007, 0.0062, 0.0152, 0.01685, 0.0189},
        {0.283, 0.286, 0.28775, 0.332, 0.3435, 0.369, 0.361, 0.260, 0.139, 0.134, 0.0625,
            0.052, 0.0315, 0.016, 0.134, 0.140, 0.149}
    };

    int wvlToSrch = 17;
    int kb = 0, ko = 0;
    int srchCounter = 0;

    /* if (wvl[nwvl-1] > bricaud_L_A_B[0][wvlToSrch - 1]){
             printf ("-E- %s line %d: Bricaud coefficient not available for band %f :(\n",
                     __FILE__,__LINE__,wvl[ko]);
             exit(1);
         }*/
    while (ko < nwvl && srchCounter < wvlToSrch) {

        if (wvl[ko] == bricaud_L_A_B[0][kb]) {
            a_bricaud[ko] = bricaud_L_A_B[1][kb];
            b_bricaud[ko] = bricaud_L_A_B[2][kb];
            ko++;
            srchCounter = 0;
        } else if (wvl[ko] < bricaud_L_A_B[0][kb]) {
            kb--;
            srchCounter += 1;
        } else {
            kb++;
            srchCounter += 1;
        }
    }
}

void grabAwBw(float wvl[], int nwvl, float aw[], float bw[]) {
    int i;
    for (i = 0; i < nwvl; i++) {
        aw[i] = aw_spectra(wvl[i], 11);
        bw[i] = bbw_spectra(wvl[i], 11);
    }
    //    float pwLAwBw[][17] = 
    //            {
    //                {410,412,413,440,443,488,490,510,530,531,547,550,555,560,665,667,670},
    //                {0.00473,0.00455056,0.00449607,0.00635,0.00706914,0.0145167,0.015,0.0325,0.0434,
    //                    0.0439153,0.0531686,0.0565,0.0596,0.0619,0.429,0.434888,0.439},
    //                {0.0059145,0.00579201,0.00573196,0.00437075,0.00424592,0.00281659,0.00276835,0.00233870,
    //                    0.00198953,0.00197385,0.00174280,0.00170333,0.00163999,0.00157958,0.000772104,
    //                    0.000762543,0.000748479}
    //            };
    //
    //    int wvlToSrch = 17;
    //    int kb=0,ko=0;
    //    int srchCounter = 0;
    //    
    //    /*if (wvl[nwvl-1] > pwLAwBw[0][wvlToSrch-1]){
    //        printf ("-E-: purewater coefficient not found for band %f >:(\n",wvl[ko]);
    //           exit(1);
    //        }
    //    */
    //    while (ko<nwvl && srchCounter<wvlToSrch){
    //        if (wvl[ko] == pwLAwBw[0][kb]){
    //            aw[ko] = pwLAwBw[1][kb];
    //            bw[ko] = pwLAwBw[2][kb];
    //            ko++;
    //            srchCounter = 0;
    //        }
    //        
    //        else if (wvl[ko] < pwLAwBw[0][kb]){
    //            kb--;
    //            srchCounter+=1;
    //        }
    //        else{ 
    //            kb++;
    //            srchCounter+=1;
    //        }
    //    }
}

void calcQAA443(float ab_surf_rrs[], float wvl[], int nwvl, ccstr *ctxt) {
    int ib;
    float *blw_surf_rrs, *u_ptr;
    float c1 = 0.52, c2 = 1.7;
    float g0 = 0.089, g1 = 0.125;
    float a0, bbp0, bbp0_exp, a0_exp;
    float ratio_aph, slope_adg, ratio_adg;
    static int firstCall = 1;
    static int rIdx = -1, gIdx = -1, bgIdx = -1, refIdx = -1, purIdx = -1;

    blw_surf_rrs = (float*) malloc(nwvl * sizeof (float));
    u_ptr = (float*) malloc(nwvl * sizeof (float));

    // get some indices:
    if (firstCall) {
        rIdx = windex(670, wvl, nwvl); // red (~670) band,  idx=5 in Melin's code
        gIdx = windex(550, wvl, nwvl); // green (~555)   ,  idx=4 ---------------
        bgIdx = windex(490, wvl, nwvl); // blue-green (~490),idx=2 ---------------
        refIdx = windex(REFWVL, wvl, nwvl); //ref (~443),        idx=1 ---------------
        purIdx = windex(412, wvl, nwvl); // purple (~412),    idx=0 ---------------
        firstCall = 0;
    }
    float bbp[refIdx + 1], a[refIdx + 1];
    // BEG ----- red band correction
    float up_667 = 20 * pow(ab_surf_rrs[gIdx], 1.5);
    float lw_667 = 0.9 * pow(ab_surf_rrs[gIdx], 1.7);
    if (ab_surf_rrs[rIdx] > up_667 || ab_surf_rrs[rIdx] < lw_667) {
        ab_surf_rrs[rIdx] = 1.27 * pow(ab_surf_rrs[gIdx], 1.47) +
                0.00018 * pow(ab_surf_rrs[gIdx] / (ab_surf_rrs[bgIdx]), -3.19);

    }
    // END ----- red band correciton

    for (ib = 0; ib < nwvl; ib++) {

        *(blw_surf_rrs + ib) = ab_surf_rrs[ib] / (c1 + c2 * ab_surf_rrs[ib]);
        *(u_ptr + ib) = (-g0 + pow(pow(g0, 2) + 4 * g1 *
                *(blw_surf_rrs + ib), 0.5)) / (2 * g1);
    }
    a0_exp = log10(
            (*(blw_surf_rrs + refIdx) + *(blw_surf_rrs + bgIdx)) /
            (*(blw_surf_rrs + gIdx) + 5 * (*(blw_surf_rrs + rIdx) /
            *(blw_surf_rrs + bgIdx)) *
            *(blw_surf_rrs + rIdx))
            );
    a0 = *(ctxt->aw_i + gIdx) + pow(10.0, (-1.146 - (1.366 * a0_exp) -(0.469 * pow(a0_exp, 2))));
    bbp0 = *(u_ptr + gIdx) * a0 / (1 - *(u_ptr + gIdx)) - *(ctxt->bbw_i + gIdx);
    bbp0_exp = 2 * (1 - 1.2 * exp(-0.9 * *(blw_surf_rrs + refIdx) / *(blw_surf_rrs + gIdx)));
    ratio_aph = 0.74 + (0.2 / (0.8 + *(blw_surf_rrs + refIdx) / *(blw_surf_rrs + gIdx)));
    slope_adg = 0.015 + 0.002 / (0.6 + *(blw_surf_rrs + refIdx) / *(blw_surf_rrs + gIdx));
    ratio_adg = exp(slope_adg * (wvl[refIdx] - wvl[purIdx]));
    for (ib = 0; ib <= refIdx; ib++) {

        bbp[ib] = bbp0 * pow((wvl[gIdx] / wvl[ib]), bbp0_exp);
        a[ib] = (1 - *(u_ptr + ib)) * (bbp[ib] + *(ctxt->bbw_i + ib)) / *(u_ptr + ib);
    }
    *(ctxt->adgRef) = ((a[purIdx] - ratio_aph * a[refIdx]) - (*(ctxt->aw_i + purIdx) - ratio_aph * *(ctxt->aw_i + refIdx))) /
            (ratio_adg - ratio_aph);
    *(ctxt->aphRef) = a[refIdx] - *(ctxt->aw_i + refIdx) - *(ctxt->adgRef);
    *(ctxt->bbpRef) = bbp[refIdx];
}

/* 
 *Interpolation for cases with more than 2 input bands used 
 */

float invDistInterp(ccstr *ctxt, float xout) {
    /*
    float invDistInterp(float xin[],float yin[],int ninout,float xout){*/
    float wt, interpOut;
    float wtSum = 0.0;
    float prodSum = 0.0;
    int in;
    for (in = 0; in<*ctxt->num_wvl_i; in++) {
        wt = 1.0 / fabs(ctxt->wvl_i[in] - xout);
        prodSum += wt * ctxt->tar_rrs[in];
        wtSum += wt;

    }
    interpOut = prodSum / wtSum;
    return interpOut;
}

/* 
 * Function to identify band or bands (at most 2) to be used 
 */

void idInputBand(float inputWvl[], float targetWvl, ccstr *ctxt, int nin) {
    /* Function to match the target Wvl with one or two input Wvl. 
     * Two wavelengths are deemed necessary if the target wvl and the closest input wvl
     * are more than 10nm apart. If this is the case and the input wvl is smaller than the 
     * target wvl, then the next higher wvl is also recruited. Similarly, if the input wvl is greater 
     * than the target wvl, the next smaller wvl is also recruited for shifting.
     */
    int closeIdx, kb;
    int nout = 2;
    int idx[2];
    closeIdx = windex(targetWvl, inputWvl, nin);
    if (((inputWvl[closeIdx] - targetWvl) > MINWVDIFF) && (closeIdx > 0)) {
        idx[0] = closeIdx - 1;
        idx[1] = closeIdx;
    } else if (((inputWvl[closeIdx] - targetWvl) < -MINWVDIFF) && (closeIdx < (nin - 1))) {
        idx[0] = closeIdx;
        idx[1] = closeIdx + 1;
    } else {
        nout = 1;
        idx[0] = closeIdx;
    }
    ctxt->sh_strt_idx = (int*) malloc(nout * sizeof (int));
    ctxt->wvl_i = (float*) malloc(nout * sizeof (float));
    ctxt->num_wvl_i = (int*) malloc(sizeof (int));
    *(ctxt->num_wvl_i) = nout;
    //*startWvlIdx = (int*)malloc(*nout*sizeof(int));
    for (kb = 0; kb < nout; kb++) {
        *(ctxt->sh_strt_idx + kb) = idx[kb];
        *(ctxt->wvl_i + kb) = inputWvl[idx[kb]];
    }
    //*(*(startWvlIdx)+kb) = idx[kb];

}

/* 
 * Actual band shifting function
 */

void shiftBand(float inputWvl[], float inputRrs[], float chla, int numBands, float tarWvl, ccstr* ctxt)
 {
    //float* wvlIn,*rrsIn;
    float wvlIn, rrsIn, refwvl;
    float s = 0.015, g0 = 0.089, g1 = 0.125, c1 = -0.52, c2 = 1.7;
    static int gIdx = -1, refIdx = -1; //,firstCall=1;
    float b2g, sdg, yy; //,chla;
    float ll_i, aph_i, adg_i, bbp_i;
    float a_tot_i, bb_tot_i, qaa_fwd_i, qaa_rrs_bbw_i, qaa_rrs_aw_i;
    float ll_o, aph_o, adg_o, bbp_o;
    float a_tot_o, bb_tot_o, qaa_fwd_o, qaa_rrs_bbw_o, qaa_rrs_aw_o, correc_factor;
    int kb, nin, idx;


    //if(firstCall){
    gIdx = windex(550, inputWvl, numBands);
    refIdx = windex(REFWVL, inputWvl, numBands);
    //firstCall=0;
    refwvl = inputWvl[refIdx];
    //}
    nin = *(ctxt->num_wvl_i);

    b2g = inputRrs[refIdx] / inputRrs[gIdx];
    sdg = s + 0.002 / (0.6 + b2g);
    yy = 2 * (1 - 1.2 * exp(-0.9 * b2g));

    //    chla = pow((*(ctxt->aphRef) / *(ctxt->aBric_i+refIdx)),1/(1-*(ctxt->bBric_i+refIdx)));
    for (kb = 0; kb < nin; kb++) {
        idx = *(ctxt->sh_strt_idx + kb);
        wvlIn = inputWvl[idx];
        rrsIn = inputRrs[idx];

        ll_i = wvlIn - refwvl;
        aph_i = aph_bricaud(wvlIn, chla); //pow( (*(ctxt->aBric_i+idx) * chla),*(ctxt->bBric_i+idx));
        adg_i = *(ctxt->adgRef) * exp(-sdg * ll_i);
        bbp_i = *(ctxt->bbpRef) * pow((443 / wvlIn), yy);
        a_tot_i = aph_i + adg_i + *(ctxt->aw_i + idx);
        bb_tot_i = bbp_i + *(ctxt->bbw_i + idx);
        qaa_fwd_i = bb_tot_i / (a_tot_i + bb_tot_i);
        qaa_rrs_bbw_i = (g0 + g1 * qaa_fwd_i) * qaa_fwd_i;
        qaa_rrs_aw_i = (c1 * qaa_rrs_bbw_i) / ((c2 * qaa_rrs_bbw_i) - 1);

        ll_o = tarWvl - REFWVL;
        aph_o = aph_bricaud(tarWvl, chla); //aph_o = pow( ( *(ctxt->aBric_o) * chla),*(ctxt->bBric_o));
        adg_o = *(ctxt->adgRef) * exp(-sdg * ll_o);
        bbp_o = *(ctxt->bbpRef) * pow((443 / tarWvl), yy);
        a_tot_o = aph_o + adg_o + *(ctxt->aw_o);
        bb_tot_o = bbp_o + *(ctxt->bbw_o);
        qaa_fwd_o = bb_tot_o / (a_tot_o + bb_tot_o);
        qaa_rrs_bbw_o = (g0 + g1 * qaa_fwd_o) * qaa_fwd_o;
        qaa_rrs_aw_o = (c1 * qaa_rrs_bbw_o) / ((c2 * qaa_rrs_bbw_o) - 1);
        correc_factor = qaa_rrs_aw_o / qaa_rrs_aw_i;
        *(ctxt->tar_rrs + kb) = rrsIn * correc_factor;
    }

}

void deMalloc(ccstr* ctxt) {
    free(ctxt->aw_i);
    free(ctxt->bbw_i);
    free(ctxt->aBric_i);
    free(ctxt->bBric_i);
    free(ctxt->aw_o);
    free(ctxt->bbw_o);
    free(ctxt->aBric_o);
    free(ctxt->bBric_o);
    free(ctxt->aphRef);
    free(ctxt->adgRef);
    free(ctxt->bbpRef);
    free(ctxt->tar_rrs);
    free(ctxt->num_wvl_i);
    free(ctxt->wvl_i);
    free(ctxt->sh_strt_idx);

}

float bioBandShift(float wvl[], float rrs[], float chla, int nw, float tarWvl) {
    ccstr sh_ctxt; //band-shifting context structure
    float result;
    idInputBand(wvl, tarWvl, &sh_ctxt, nw);
    sh_ctxt.aw_i = (float*) malloc(nw * sizeof (float));
    sh_ctxt.bbw_i = (float*) malloc(nw * sizeof (float));
    sh_ctxt.aBric_i = (float*) malloc(nw * sizeof (float));
    sh_ctxt.bBric_i = (float*) malloc(nw * sizeof (float));
    sh_ctxt.aw_o = (float*) malloc(sizeof (float));
    sh_ctxt.bbw_o = (float*) malloc(sizeof (float));
    sh_ctxt.aBric_o = (float*) malloc(sizeof (float));
    sh_ctxt.bBric_o = (float*) malloc(sizeof (float));
    sh_ctxt.aphRef = (float*) malloc(sizeof (float));
    sh_ctxt.adgRef = (float*) malloc(sizeof (float));
    sh_ctxt.bbpRef = (float*) malloc(sizeof (float));
    sh_ctxt.tar_rrs = (float*) malloc(sizeof (float));

    grabAwBw(wvl, nw, sh_ctxt.aw_i, sh_ctxt.bbw_i);
    grabAwBw(&tarWvl, 1, sh_ctxt.aw_o, sh_ctxt.bbw_o);

    calcQAA443(rrs, wvl, nw, &sh_ctxt);
    shiftBand(wvl, rrs, chla, nw, tarWvl, &sh_ctxt);
    if (*sh_ctxt.num_wvl_i > 1)
        result = invDistInterp(&sh_ctxt, tarWvl);

    else
        result = *(sh_ctxt.tar_rrs);
    deMalloc(&sh_ctxt);
    return result;
}
