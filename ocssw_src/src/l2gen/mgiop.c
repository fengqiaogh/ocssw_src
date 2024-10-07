#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"
#include "chl.h"
#include <stdio.h>
#include "giop.h"

static float badval = BAD_FLT;

void get_mgiop(l2str *l2rec, l2prodstr *p, float prod[]) {
    static int firstCall = 1;
    int prodID = p->cat_ix;
    int ib = p->prod_ix;

    l1str *l1rec = l2rec->l1rec;
    int32 nwave = l1rec->l1file->nbands;
    int32 npix = l1rec->npix;

    int32_t ip, ipb, ibb;
    int isdg, ieta, ipar;
    int npar = 3;

    float def_sdg[] = {0.01, 0.015, 0.02};
    float def_eta[] = {0.0, 0.33, 0.67, 1.0};

    static int32_t *nval;
    static float *schl;
    static float *sbbp;
    static float *sadg;
    static float *saph;
    static float *snoc;
    static float *sdia;
    static float **spar;

    // pointers from giop.c
    static float *chl;
    static float *adg;
    static float *bbp;
    static float *aph;
    static float **fit_par;

    if (firstCall) {

        firstCall = 0;

        // get number of eigenvalues (= # aph + 1 eta + 1 Sdg)
        //npar = table_column_count(input->giop_aph_file)+1;  // memory leak here?
        if (strcmp(input->giop_aph_file, "aph_nas"))
            npar = 4;

        if ((aph = calloc(npix * nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for NAS.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((adg = calloc(npix * nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for NAS.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((bbp = calloc(npix * nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for NAS.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((chl = calloc(npix, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for NAS.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((fit_par = allocate2d_float(npix, npar)) == NULL) {
            printf("-E- %s line %d : error allocating memory for NAS.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        if ((schl = calloc(npix, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for NAS.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((snoc = calloc(npix, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for NAS.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((sdia = calloc(npix, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for NAS.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((nval = calloc(npix, sizeof (double))) == NULL) {
            printf("-E- %s line %d : error allocating memory for NAS.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((sbbp = calloc(npix * nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for NAS.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((sadg = calloc(npix * nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for NAS.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((saph = calloc(npix * nwave, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for NAS.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((spar = allocate2d_float(npix, npar)) == NULL) {
            printf("-E- %s line %d : error allocating memory for NAS.\n",
                    __FILE__, __LINE__);
            exit(1);
        }


        // set bbp spectral shape to user defined power-law
        // set aph spectral shape to tabulated
        //input->giop_bbp_opt=1;
        //input->giop_aph_opt=0;

        // use NAS-specific tabulated aph
        //parse_file_name("$OCDATAROOT/common/aph_nas_early.txt", tmp_file);
        //strcpy(input->giop_aph_file, tmp_file);

        // further restrict rrsdiff_max and enable iteration
        //input->giop_rrs_diff=0.1;
        input->giop_iterate = 1;
    }

    // initialize arrays and counters

    for (ip = 0; ip < npix; ip++) {

        nval[ip] = 0;
        schl[ip] = 0.0;
        sdia[ip] = 0.0;
        snoc[ip] = 0.0;
        chl[ip] = badval;

        ipb = ip*nwave;
        for (ibb = 0; ibb < nwave; ibb++) {
            sbbp[ipb + ibb] = 0.0;
            saph[ipb + ibb] = 0.0;
            sadg[ipb + ibb] = 0.0;
            bbp[ipb + ibb] = badval;
            aph[ipb + ibb] = badval;
            adg[ipb + ibb] = badval;
        }

        for (ipar = 0; ipar < npar; ipar++) {
            fit_par[ip][ipar] = badval;
            spar[ip][ipar] = 0;
        }
    }


    // loop through three Sdg and four eta

    for (isdg = 0; isdg < 3; isdg++) for (ieta = 0; ieta < 4; ieta++) {

            input->giop_bbp_s = def_eta[ieta];
            input->giop_adg_s = def_sdg[isdg];

            run_giop(l2rec);

            chl = giop_get_chl_pointer();
            adg = giop_get_adg_pointer();
            aph = giop_get_aph_pointer();
            bbp = giop_get_bbp_pointer();
            fit_par = giop_get_fitpar_pointer();

            for (ip = 0; ip < l1rec->npix; ip++) {

                ipb = ip * nwave + ib;

                if (chl[ip] > 0.005 && chl[ip] < 200.0) {

                    nval[ip]++;

                    schl[ip] += chl[ip];
                    sdia[ip] += fit_par[ip][0];
                    snoc[ip] += fit_par[ip][1];

                    sadg[ipb] += adg[ipb];
                    saph[ipb] += aph[ipb];
                    sbbp[ipb] += bbp[ipb];

                    for (ipar = 0; ipar < npar; ipar++) {
                        spar[ip][ipar] += fit_par[ip][ipar];
                    }

                }
            }
        }


    for (ip = 0; ip < npix; ip++) {

        // flag and skip if pixel already masked

        if (l1rec->mask[ip]) {
            l1rec->flags[ip] |= PRODFAIL;
            continue;
        }

        ipb = ip * nwave + ib;

        switch (prodID) {

        case CAT_aph_mgiop:
            prod[ip] = (float) saph[ipb] / nval[ip];
            break;

        case CAT_adg_mgiop:
            prod[ip] = (float) sadg[ipb] / nval[ip];
            break;

        case CAT_bbp_mgiop:
            prod[ip] = (float) sbbp[ipb] / nval[ip];
            break;

        case CAT_chl_mgiop:
            prod[ip] = (float) schl[ip] / nval[ip];
            break;

        case CAT_npix_mgiop:
            prod[ip] = (int) nval[ip];
            break;

        case CAT_crat_mgiop:
            prod[ip] = (float) fabs(snoc[ip] / nval[ip]) / fabs(sdia[ip] / nval[ip]);
            break;

        case CAT_fitpar_mgiop:
            prod[ip] = (float) spar[ip][ib] / nval[ip];
            break;

        default:
            printf("-E- %s line %d : erroneous product ID %d passed to GIOP.\n",
                    __FILE__, __LINE__, prodID);
            exit(1);
        }
    }

    return;

}

