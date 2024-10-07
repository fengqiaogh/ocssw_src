/*
 * @file bioOptBandShift.c
 * @Brief Bio-optical band shift based on Melin & Sclep (2015).
 * @author Erdem Karakoylu
 * @author Ocean Biology Processing Group, Nasa Goddard Space Flight Center
 * 
 * Version History:
 * V. 0.1:  - hardwired iop algorithm and LUT to get aw/bbw, adg(443),aph(443),bbp(443)
 *          - aph(443) computed w/ Bricaud 95 Coeffs. 
 * V. 0.2:  - Removed  hardwired qaa & and aw, bbw code
 *          - Added interface (qaa_iops_4_bshift) to get_qaa.c; returns adg(443),aph(443),bbp(443). 
 *          - Added access to aph_bricaud_1995 in aph.c and get_aw_bbw  
 * V. 0.3:  - Replaced aph_bricaud_1995 by aph_bricaud_1998 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "l12_proto.h"
#define REFWVL 443.0

typedef struct _context {
    int *ibvc_bin;
    float greenband;
    float *wvl_i, *wvl_o;
    float *aw_i, *bbw_i;
    float *aw_o, *bbw_o;
    float aw_ref, bbw_ref, wvl_ref;
    int nBandPairs;
} ccstr;

float BandShift(ccstr*, l2str*, int32_t, float);
void CheckNULL(void*, char*, int, char*);
float ToAboveWater(float);
float ToBelowWater(float);
float Get_bb(l2str*, ccstr*, int32_t, int, char);
float Get_a(l2str*, ccstr*, int32_t, int, char);
float BandInterp(ccstr*, float*, int, int);
void Setup(l2str*, ccstr*);

void qaa_iops_4_bshift(l2str*, float*, float*);
//void pml_iops_4_bshift(l2str*,float*,float*);
// void gsm_iops_4_bshift(l2str*,float*,float*);
float aph_bricaud_1998(float, float);

static int nband_vc = 7;
static float vc_bands[7] = {412., 443., 490., 510., 531., 555., 670.};
static int current_iscan = -1;
static int current_npix = -1;
static float *sh_prod_arr, *adg_ref, *bbp_ref;

void bioOptBandShift(l2str *l2rec, l2prodstr *p, float prod[]) {

    int32_t ibvc = p->prod_ix;
    int32_t ib, ip, ipb;
    static ccstr vcc;
    static int firstRun = 1;
    if (firstRun) {
        Setup(l2rec, &vcc);
        firstRun = 0;
    }
    /* sh_prod_arr is a 2D array accessed as a simple pointer of dims npix rows, and nband_vc cols.
       If band shift not necessary, corresponding cols filled w/ Rrs already in l2rec. 
       Occasional segment below -- happens only when switching to new row 
       avoids recalculation for multiple products otherwise*/
    if (l2rec->l1rec->iscan != current_iscan) {
        current_iscan = l2rec->l1rec->iscan;
        current_npix = l2rec->l1rec->npix;
        if (adg_ref != NULL)
            free(adg_ref);
        CheckNULL(adg_ref = (float*) malloc(current_npix * sizeof (float)),
                "adg_ref", __LINE__, __FILE__);
        if (bbp_ref != NULL)
            free(bbp_ref);
        CheckNULL(bbp_ref = (float*) malloc(current_npix * sizeof (float)),
                "bbp_ref", __LINE__, __FILE__);
        if (sh_prod_arr != NULL)
            free(sh_prod_arr);
        CheckNULL(sh_prod_arr = (float*) malloc(current_npix * nband_vc * sizeof (float))
                , "sh_prod_arr", __LINE__, __FILE__);
        qaa_iops_4_bshift(l2rec, adg_ref, bbp_ref);
        //pml_iops_4_bshift(l2rec,adg_ref,bbp_ref);
        //gsm_iops_4_bshift(l2rec,adg_ref,bbp_ref);
        for (ip = 0; ip < current_npix; ip++) {
            for (ib = 0; ib < nband_vc; ib++) {
                // fill sh_prod_arr, by bandshifting if necessary, by copying l2rec rrs otherwise. 
                if (vcc.ibvc_bin[ib]) {
                    *(sh_prod_arr + ip * nband_vc + ib) = BandShift(&vcc, l2rec, ip, vc_bands[ib]);
                } else {

                    ipb = ip * l2rec->l1rec->l1file->nbands + bindex_get((int32_t) vc_bands[ib]);
                    *(sh_prod_arr + ip * nband_vc + ib) = l2rec->Rrs[ipb];
                }
            }
        }
    }

    for (ip = 0; ip < l2rec->l1rec->npix; ip++)
        prod[ip] = *(sh_prod_arr + ip * nband_vc + ibvc);

}

void Setup(l2str *l2rec, ccstr *vccp) {
    int kb, iflag;

    float lambda_i[4][5] = {
        {510., 555., -1., -1., -1.}, //SWF
        {488., 488., 531., 547., 667.}, //MODA
        {510., 560., 560., 665., -1.}, //MER
        {488., 488., 555., 555., 670.} //VIIRS
    };

    float lambda_o[4][5] = {
        {531., 531., -1., -1., -1.}, //SWF
        {490., 510., 510., 555., 670.}, //MODA
        {531., 531., 555., 670., -1.}, //MER
        {490., 510., 510., 531., 670.} //VIIRS
    };

    int ibvc_bin[4][7] = {//index of band needing bandshift for each sensor
        {0, 0, 0, 0, 1, 0, 0}, //SWF
        {0, 0, 1, 1, 0, 1, 1}, //MODA
        {0, 0, 0, 0, 1, 1, 1}, //MER
        {0, 0, 1, 1, 1, 0, 1} //VIIRS
    };
    vccp->wvl_ref = REFWVL;
    switch (l2rec->l1rec->l1file->sensorID) {
    case SEAWIFS:
    {
        vccp->nBandPairs = 2; // 2 input bands will be used to compute 531.
        vccp->greenband = 555;
        iflag = 0;
        break;
    }

    case MODISA:
    {
        vccp->nBandPairs = 5; //inp bands: 488,488,531,547,667, resp
        vccp->greenband = 547;
        iflag = 1;
        break;
    }

    case MERIS:
    {
        vccp->nBandPairs = 5;
        vccp->greenband = 560;
        iflag = 2;
        break;
    }
    case VIIRSN:
    {
        vccp->nBandPairs = 5;
        vccp->greenband = 555;
        iflag = 3;
        break;
    }
    }

    CheckNULL(vccp->ibvc_bin = (int*) malloc(nband_vc * sizeof (int)), "vccp->ibvc_bin", __LINE__, __FILE__);
    CheckNULL(vccp->wvl_i = (float*) malloc(vccp->nBandPairs * sizeof (float)), "vccp->wvl_i", __LINE__, __FILE__);
    CheckNULL(vccp->wvl_o = (float*) malloc(vccp->nBandPairs * sizeof (float)), "vccp->wvl_o", __LINE__, __FILE__);
    CheckNULL(vccp->aw_i = (float*) malloc(vccp->nBandPairs * sizeof (float)), "vccp->aw_i", __LINE__, __FILE__);
    CheckNULL(vccp->aw_o = (float*) malloc(vccp->nBandPairs * sizeof (float)), "vccp->aw_o", __LINE__, __FILE__);
    CheckNULL(vccp->bbw_i = (float*) malloc(vccp->nBandPairs * sizeof (float)), "vccp->aw_all", __LINE__, __FILE__);
    CheckNULL(vccp->bbw_o = (float*) malloc(vccp->nBandPairs * sizeof (float)), "vccp->bbw_o", __LINE__, __FILE__);

    for (kb = 0; kb < vccp->nBandPairs; kb++) {
        vccp->wvl_i[kb] = lambda_i[iflag][kb];
        vccp->wvl_o[kb] = lambda_o[iflag][kb];

    }
    for (kb = 0; kb < nband_vc; kb++)
        vccp->ibvc_bin[kb] = ibvc_bin[iflag][kb];

    get_aw_bbw(l2rec, vccp->wvl_i, vccp->nBandPairs, vccp->aw_i, vccp->bbw_i);
    get_aw_bbw(l2rec, vccp->wvl_o, vccp->nBandPairs, vccp->aw_o, vccp->bbw_o);
    get_aw_bbw(l2rec, &vccp->wvl_ref, 1, &vccp->aw_ref, &vccp->bbw_ref);
}

float BandInterp(ccstr* vcc, float* rrs_sh, int idx_0, int num_in) {
    /* Inverse distance interpolation
     * E.g. for bandshift from SeaWiFS, rrs(510) and rrs(555) are used to calculate
     * two separate rrs(531). These are used here to calculate a final rrs(531).
     * Max of 2 input bands is assumed.
     */

    float *wts, wts_sum = 0;
    float rrs_sh_interp = 0;
    int idx, ib;
    CheckNULL(wts = (float*) malloc(num_in * sizeof (float)), "wts", __LINE__, __FILE__);

    for (idx = 0; idx < num_in; idx++) {
        ib = idx_0 + idx;
        wts[idx] = 1.0 / fabs(vcc->wvl_o[ib] - vcc->wvl_i[ib]);
        wts_sum += wts[idx];
    }

    for (idx = 0; idx < num_in; idx++) {
        rrs_sh_interp += (wts[idx] * rrs_sh[idx]) / (wts_sum);

    }
    free(wts);
    return rrs_sh_interp;
}

float BandShift(ccstr* vccp, l2str *l2rec, int32_t ip, float target_band) {
    float g0 = 0.089, g1 = 0.1245;
    float rrs_in, *rrs_e_f, rrs_f_i, rrs_f_f, rrs_e_temp, Rrs_e_f;
    int ipb, ib, idx, t_num = 0, t_start;
    float bb_i, bb_f;
    float a_i, a_f, u_i, u_f;

    // locate target_band
    for (ib = 0; ib < vccp->nBandPairs; ib++) {
        if (vccp->wvl_o[ib] == target_band) {
            if (t_num == 0)
                t_start = ib;
            t_num++;
        }
    }

    if (t_num == 0) {
        printf("\n -E- line %d Unable to find target band,%f => t_num = 0 :( \n", __LINE__, target_band);
        exit(FATAL_ERROR);
    }

    CheckNULL(rrs_e_f = (float*) malloc(t_num * sizeof (float)), "rrs_e_f", __LINE__, __FILE__);

    for (idx = 0; idx < t_num; idx++) {
        ib = t_start + idx;
        ipb = ip * l2rec->l1rec->l1file->nbands + bindex_get((int32_t) vccp->wvl_i[ib]);
        rrs_in = ToBelowWater(l2rec->Rrs[ipb]);
        // get IOPs
        bb_i = Get_bb(l2rec, vccp, ip, ib, 'i');
        bb_f = Get_bb(l2rec, vccp, ip, ib, 'f');
        a_i = Get_a(l2rec, vccp, ip, ib, 'i');
        a_f = Get_a(l2rec, vccp, ip, ib, 'f');
        // get u
        u_i = bb_i / (a_i + bb_i);
        u_f = bb_f / (a_f + bb_f);
        // rrs (input & output)  from forward model
        rrs_f_i = g0 * u_i + g1 * pow(u_i, 2);
        rrs_f_f = g0 * u_f + g1 * pow(u_f, 2);
        // estimate shifted rrs
        rrs_e_f[idx] = rrs_in * rrs_f_f / rrs_f_i;

    }
    if (t_num > 1) {//run inverse distance weighted interpolation if > 1 input band
        rrs_e_temp = BandInterp(vccp, rrs_e_f, t_start, t_num);
        Rrs_e_f = ToAboveWater(rrs_e_temp);
    } else {
        Rrs_e_f = ToAboveWater(rrs_e_f[0]);
    }
    free(rrs_e_f);
    return Rrs_e_f;
}

void CheckNULL(void *ptr, char* ptrName, int lineNum, char* filename) {
    if (ptr == NULL) {
        printf("\n --E-- pointer %s at line %d is %s not allocated :(\n", ptrName, lineNum, filename);
        exit(FATAL_ERROR);
    }
}

float Get_bb(l2str *l2rec, ccstr *vccp, int32_t ip, int vccidx, char key) {
    int refIdx, greenIdx, ipb;
    float rrsRef, rrsGreen, eta, wavl;
    float bbw, bbp, bb;

    refIdx = bindex_get((int32_t) REFWVL);
    ipb = ip * l2rec->l1rec->l1file->nbands + refIdx;
    rrsRef = ToBelowWater(l2rec->Rrs[ipb]);

    greenIdx = bindex_get(vccp->greenband);
    ipb = ip * l2rec->l1rec->l1file->nbands + greenIdx;
    rrsGreen = ToBelowWater(l2rec->Rrs[ipb]);

    if (key == 'i') {
        wavl = vccp->wvl_i[vccidx];
        bbw = vccp->bbw_i[vccidx];
    } else {
        wavl = vccp->wvl_o[vccidx];
        bbw = vccp->bbw_o[vccidx];
    }
    eta = 2. * (1 - 1.2 * exp(-0.9 * rrsRef / rrsGreen));
    bbp = bbp_ref[ip] * pow((REFWVL / wavl), eta);
    bb = bbp + bbw;
    return (bb);
}

float Get_a(l2str *l2rec, ccstr *vccp, int32_t ip, int vccidx, char key) {

    int refIdx, greenIdx, ipb;
    float rrsRef, rrsGreen, wavl, s;
    float aw, aph, adg, a;

    refIdx = bindex_get((int32_t) REFWVL);
    ipb = ip * l2rec->l1rec->l1file->nbands + refIdx;
    rrsRef = ToBelowWater(l2rec->Rrs[ipb]);

    greenIdx = bindex_get((int32_t) vccp->greenband);
    ipb = ip * l2rec->l1rec->l1file->nbands + greenIdx;
    rrsGreen = ToBelowWater(l2rec->Rrs[ipb]);

    if (key == 'i') {
        wavl = vccp->wvl_i[vccidx];
        aw = vccp->aw_i[vccidx];
    } else {
        wavl = vccp->wvl_o[vccidx];
        aw = vccp->aw_o[vccidx];
    }
    s = 0.015 + 0.002 / (0.6 + (rrsRef / rrsGreen));
    adg = adg_ref[ip] * exp(-s * (wavl - REFWVL));
    aph = aph_bricaud_1998(wavl, l2rec->chl[ip]);
    a = aph + adg + aw;
    return (a);
}

float ToBelowWater(float rrsAbove) {
    float rrsBelow;
    rrsBelow = rrsAbove / (0.52 + 1.7 * rrsAbove);
    return rrsBelow;
}

float ToAboveWater(float rrsBelow) {
    float rrsAbove;
    rrsAbove = 0.52 * rrsBelow / (1 - 1.7 * rrsBelow);
    return rrsAbove;
}
