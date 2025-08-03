/*---------------------------------------------------------------------*/
/* calcite.c -  get calcium carbonate concentration.                   */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/* Outputs:                                                            */
/*     caco3 - calcium carbonate concentration, per pixel        .     */
/*                                                                     */
/* Written by: W. Robinson, GSC, 7 Jun 2000.                           */
/*             S. Bailey, OCDPG, July 2004, conversion to C.           */
/*             B. Franz, OCDPG, Sep 2004, sensor generalization and    */
/*                 implementation of 2-Band algorithm.                 */
/*                                                                     */
/*             2014: Standardized to use common table for 2-band alg   */
/*             and adjust green nLw as needed for sensor.              */
/*                                                                     */
/*             2014: Changed bbstar from 4 to 1.628.                   */
/*                                                                     */
/*---------------------------------------------------------------------*/

#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"
#include "l2_flags.h"

#define BAD_CACO3 BAD_FLT

static int32_t  caco3_msk = LAND | HIGLINT | CLOUD | HILT;
static float bbstr = 1.628;
static float caco3min = 1e-5;     // BCB - was 1.18e-5
static float caco3hi   = 0.0005;  // BCB - was 0.003
static float fixedbbstar = 1.28;

/* --------------------------------------------------------------------- */
/* calcite_3b() - calcium carbonate concentration from 3-Band algorith.. */
/*                                                                       */
/* Gordon, H.R. Boynton, G.C., Balch, W.M., Groom, S.B., Harbour, D.S.,  */
/* Smyth, T.J., Retrieval of Coccolithophore Calcite Concentration from  */
/* SeaWiFS Imagery, GRL, 28, 8, 1587-1590.                               */
/*                                                                       */
/* --------------------------------------------------------------------- */

float calcite_3b(l2str *l2rec, int32_t ip) {
    static int firstCall = 1;
    static int maxiter = 10;
    static float ftrans = 6.179; /* (1/.298)*(1/.543)   */

    static float wave[3] = {670., 760., 870.}; /* approx. wavelengths */
    static int bx [3];
    static float aw [3];
    static float bbw [3];
    static float bbc [3];
    static float t [3];
    static float b68diff;
    static float b78diff;
    static float fw1, fw2;

    static float oobswf[3][8] = {
        {0.000313529, 0.000770558, 0.00152194, 0.000155573,
            0.00116455, 0.0, 0.000445433, 0.000124172},
        {0.000201709, 6.96143e-05, 7.00147e-06, 2.28957e-07,
            4.17788e-05, 0.00159814, 0.0, 0.00536827},
        {0.000463807, 8.54003e-05, 2.47401e-05, 0.000755890,
            0.00587073, 0.00021686, 0.0111331, 0.0}
    };

    int32_t ipb, ib, i;
    float *rhoaw;
    float rho[3];
    float bbc_cclth, r8_cclth, aeps_cclth;
    float bbcinit, bbctol;
    int numiter;
    int32_t nwave, status = 0;
    float *awptr, *bbwptr;
    float caco3;
    float bbw546;
    float newbbstar; // BCB - new, dynamic bb*

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;
    nwave = l1file->nbands;

    if (firstCall) {

        /* save coeffs for the three bands, resolve actual sensor wave */
        for (i = 0; i < 3; i++) {
            bx [i] = windex(wave[i], l1file->fwave, nwave);
            wave[i] = l1file->fwave[bx[i]];

            bbc [i] = 0.0;
        }

        b68diff = wave[2] - wave[0];
        b78diff = wave[1] - wave[0];

        /* spectral dependence of bbc */
        fw1 = pow(wave[0] / wave[1], 1.35);
        fw2 = pow(wave[0] / wave[2], 1.35);

        firstCall = 0;
    }

    /* set aw & bbw */

    ipb = ip*nwave;
    awptr = &l1rec->sw_a_avg[ipb];
    bbwptr = &l1rec->sw_bb_avg[ipb];

    for (i = 0; i < 3; i++) {
        aw [i] = awptr [bx[i]];
        bbw [i] = bbwptr[bx[i]];
    }
    bbw546 = seawater_bb(546.0, l1rec->sstref[ip], l1rec->sssref[ip], 0.039);

    status = 0;
    numiter = 0;
    bbctol = 100.;
    bbcinit = 0.00085;
    caco3 = BAD_CACO3;
    bbc[0] = 0.000;

    /* skip pixel if already masked (this should not include ATMFAIL) */
    if ((l1rec->flags[ip] & caco3_msk) != 0) {
        return (caco3);
    }

    if ((rhoaw = (float *) calloc(nwave, sizeof (float))) == NULL) {
        printf("-E- : Error allocating memory to rhoaw\n");
        exit(FATAL_ERROR);
    }

    /* compute the aerosol/water reflectance (include out-of-band correction for SeaWiFS) */
    if (l1file->sensorID == SEAWIFS) {
        for (ib = 0; ib < nwave; ib++) {
            ipb = nwave * ip + ib;
            rhoaw[ib] = ((l2rec->l1rec->Lt[ipb] / l1rec->tg_sol[ipb] / l1rec->tg_sen[ipb] - l1rec->tLf[ipb]
                    - l1rec->Lr[ipb]) / l1rec->t_o2[ipb] - l1rec->TLg[ipb]) * OEL_PI / l1rec->Fo[ib] / l1rec->csolz[ip];
        }
        for (i = 0; i < 3; i++) {
            rho[i] = rhoaw[bx[i]];
            for (ib = 0; ib < nwave; ib++) {
                rho[i] -= rhoaw[ib] * oobswf[i][ib];
                if (rho[i] <= 0.0) status = 1;
            }
        }
    } else {
        for (i = 0; i < 3; i++) {
            ib = bx[i];
            ipb = nwave * ip + ib;
            rho[i] = ((l2rec->l1rec->Lt[ipb] / l1rec->tg_sol[ipb] / l1rec->tg_sen[ipb] - l1rec->tLf[ipb]
                    - l1rec->Lr[ipb]) / l1rec->t_o2[ipb] - l1rec->TLg[ipb]) * OEL_PI / l1rec->Fo[ib] / l1rec->csolz[ip];
            if (rho[i] <= 0.0) status = 1;
        }
    }

    /* skip pixel on negative surface reflectance */
    if (status != 0) {
        free(rhoaw);
        return (caco3);
    }

    /* compute total transmittance */
    for (i = 0; i < 3; i++) {
        ipb = nwave * ip + bx[i];
        t[i] = l1rec->tg_sol[ipb] * l1rec->tg_sen[ipb] * l1rec->t_sol[ipb] * l1rec->t_sen[ipb];
    }

    /* compute backscatter at 546 nm */
    while (bbctol > 5. && numiter < maxiter) {

        numiter++;

        bbc[1] = bbc[0] * fw1;
        bbc[2] = bbc[0] * fw2;

        /* reflectance at longest wavelength */
        r8_cclth = rho[2] - (bbw[2] + bbc[2]) / (aw[2] + bbw[2] + bbc[2]) / ftrans * t[2];

        if ((r8_cclth > 0.09) || (r8_cclth < 0)) {
            status = 1;
            bbc[0] = 0;
            break;
        }

        /* atmospheric epsilon at two longest wavelengths */
        aeps_cclth = log((rho[1] - (bbw[1] + bbc[1]) / (aw[1] + bbw[1] + bbc[1]) / ftrans * t[1]) / r8_cclth) / b78diff;

        if (aeps_cclth > 0.4) {
            status = 1;
            bbc[0] = 0;
            break;
        }

        /* --------------- */
        bbc[0] = (rho[0] - r8_cclth * exp(aeps_cclth * b68diff)) / t[0] * (aw[0] + bbw[0] + bbc[0]) * ftrans - bbw[0];

        if ((bbc[0] <= 0) || isnan(bbc[0])) {
            status = 1;
            bbc[0] = 0;
            break;
        }

        bbctol = fabs((bbcinit - bbc[0]) / bbcinit)*100.;
        bbcinit = bbc[0];
    }


    if (status == 0) {

        bbc_cclth = bbc[0] / pow((546. / wave[0]), 1.35) - bbw546;
        if (bbc_cclth > 0) {
            // convert to calcite in moles/m^3 (Balch 2005)
            caco3 = bbc_cclth / bbstr;
            newbbstar = pow(10,(0.2007476*log10(caco3)*log10(caco3) + 1.033187*log10(caco3) + 1.069821));

            if (newbbstar < fixedbbstar) {
                newbbstar = fixedbbstar;
            }

            caco3 = caco3*bbstr/newbbstar;
        }
    }

    free(rhoaw);

    return (caco3);
}


/* --------------------------------------------------------------------- */
/* calcite_2b() - calcium carbonate concentration from 2-Band algorith.. */
/*                                                                       */
/* Gordon, H.R. and Balch, W.M., MODIS Detached Coccolith Concentration  */
/* Algorithm Theoretical Basis Document,  April 30, 1999                 */
/*                                                                       */
/* --------------------------------------------------------------------- */

#define N443 490
#define N550 456

float calcite_2b(l2str *l2rec, int32_t ip) {
    static int firstCall = 1;
    static int bandShift = 0;

    static float* t443;
    static float* t550;
    typedef float tbb_t[N550];
    static tbb_t* tbb;
    static int32_t ib443;
    static int32_t ib550;

    float caco3;
    int i443 = 0, i550 = 0;
    float x, y, a, b;
    int32_t nwave, i, nc;
    float x443, x550;
    float bb1, bb2, bb;
    float newbbstar; // BCB - new, dynamic bb*

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;
    nwave = l1file->nbands;

    if (firstCall) {

        t443 = (float*) allocateMemory(N443 * sizeof (float), "t443");
        t550 = (float*) allocateMemory(N550 * sizeof (float), "t550");
        tbb = (tbb_t*) allocateMemory(N443 * sizeof (tbb_t), "tbb");

        FILE *fp;
        char filename[FILENAME_MAX];
        char line [80];

        strcpy(filename, input->picfile);
        if (strlen(filename) == 0) {
            printf("-E- %s line %d: No picfile specified.\n", __FILE__, __LINE__);
            exit(1);
        }
        printf("Loading PIC 2-band algorithm table %s\n", filename);

        ib443 = windex(443., l1file->fwave, nwave);
        ib550 = windex(550., l1file->fwave, nwave);

        if (strstr(filename, "common") != NULL) {
            printf("Assuming PIC table is for 443nm and 555nm.\n");
            bandShift = 1;
            ib443 = bindex_get(443);
            ib550 = bindex_get_555(l1file->sensorID);
            if (ib443 < 0 || ib550 < 0) {
                printf("-E- %s line %d: required bands not available PIC\n",
                        __FILE__, __LINE__);
                exit(1);
            }
        }

        if ((fp = fopen(filename, "r")) == NULL) {
            fprintf(stderr, "-E- %s line %d: unable to open %s for reading\n",
                    __FILE__, __LINE__, filename);
            exit(1);
        }

        // Skip comment lines
        nc = 0;
        while (fgets(line, 80, fp)) {
            if (line[0] != '#' && line[0] != '\n') {
                break;
            }
            nc++;
        }
        rewind(fp);
        for (i = 0; i < nc; i++)
            fgets(line, 80, fp);

        // Load table

        for (i443 = 0; i443 < N443; i443++)
            for (i550 = 0; i550 < N550; i550++) {
                fscanf(fp, "%f %f %f %f\n", &x, &y, &a, &b);
                t443[i443] = x;
                t550[i550] = y;
                tbb [i443][i550] = a;
            }

        firstCall = 0;
    }

    caco3 = BAD_CACO3;

    // skip pixel if already masked (this includes ATMFAIL)
    if (l1rec->mask[ip]) {
        return (caco3);
    }

    x443 = l2rec->nLw[ip * nwave + ib443];
    x550 = l2rec->nLw[ip * nwave + ib550];

    // if required radiances are negative, fail; BCB - correct comment
    if (x550 <= 0.0) {
        return (caco3);
    }
    if (x443 <= 0.0) {
        return (caco3);
    }

    // adjust nLw to 555 based on Rrs555/Rrs5xx ratio
    if (bandShift) {
        float Rrs555 = conv_rrs_to_555(l2rec->Rrs[ip * nwave + ib550], l1file->fwave[ib550],-99, NULL);
        x550 *= Rrs555 / l2rec->Rrs[ip * nwave + ib550];
    }

    // locate bounding table indices
    for (i = 0; i < N443; i++) {
        if (x443 < t443[i]) {
            i443 = i;
            break;
        }
    }
    if (x443 >= t443[N443-1]) // BCB - Explicit check for high side of the table
      i443 = N443;
    
    for (i=0; i<N550; i++) {
        if (x550 < t550[i]) {
            i550 = i;
            break;
        }
    }
    if (x550 >= t550[N550-1]) // BCB - Explicit check for high side of the table
      i550 = N550;

    // radiances less than table entries, fail and don't call 3band; BCB - change failure
    if (i443 <=0) {  // BCB - is nLw443 is lower than any entry in table?
      return(-1.0);  // BCB - yes, set failure mode
    }
    if (i550 <=0) {  // BCB - is nLw550 is lower than any entry in table?
      return(-1.0);  // BCB - yes, set failure mode
    }

    // radiances greater than table entries, fail
    if (i443 >= N443 || i550 >= N550) {
        return (caco3);
    }

    // radiances associated with missing table entries, fail
    if (tbb[i443 - 1][i550 - 1] > 998.9 || tbb[i443 ][i550 - 1] > 998.9 ||
            tbb[i443 - 1][i550 ] > 998.9 || tbb[i443 ][i550 ] > 998.9) {
        return (caco3);
    }

    // interpolate to get bb(546) 
    bb1 = tbb[i443 - 1][i550 - 1] + (x443 - t443[i443 - 1])*
            (tbb[i443][i550 - 1] - tbb[i443 - 1][i550 - 1]) / (t443[i443] - t443[i443 - 1]);

    bb2 = tbb[i443 - 1][i550 ] + (x443 - t443[i443 - 1])*
            (tbb[i443][i550 ] - tbb[i443 - 1][i550 ]) / (t443[i443] - t443[i443 - 1]);

    bb = bb1 + (x550 - t550[i550 - 1]) * (bb2 - bb1) / (t550[i550] - t550[i550 - 1]);

    // convert to calcite in moles/m^3 (Balch 2005)
    caco3 = bb / bbstr;

    newbbstar = pow(10,(0.2007476*log10(caco3)*log10(caco3) + 1.033187*log10(caco3) + 1.069821));

    if (newbbstar < fixedbbstar) {
      newbbstar = fixedbbstar;
    }

    caco3 = caco3*bbstr/newbbstar;

    //    if (caco3 < caco3min ) {
    //       caco3 = caco3min;
    //}

    return(caco3);
}



/* --------------------------------------------------------------------- */
/* calcite_c() - calcium carbonate concentration (combined algorithm)    */
/* --------------------------------------------------------------------- */
float calcite_c(l2str *l2rec, int32_t ip) {
    float caco3 = BAD_CACO3;
    int32_t shallow;
    int32_t turbid;
    int32_t shallowDepth = 30;          // BCB - how shallow is too shallow?

    turbid  = ((l2rec->l1rec->flags[ip] & TURBIDW) != 0);    // BCB - is the TURBID flag set?
    shallow = (abs(l2rec->l1rec->dem[ip]) < shallowDepth);  // BCB - use depth, rather than SHALLOW flag


    if (turbid & shallow) {
      caco3 = BAD_CACO3;                // BCB - if it's shallow and turbid, we don't believe it.
    } else {
        caco3 = calcite_2b(l2rec,ip);     // BCB - calculate 2 band value
        if (caco3 < 0.0) {                // BCB - did it fail?
            if (caco3 < -2.0) {             // BCB - yes, how did it fail?
                caco3 = calcite_3b(l2rec,ip); // BCB - normal 2 band failure, call the 3 band
                if (caco3 < caco3hi) {        // BCB - is 3 band value believable?
                    caco3 = BAD_CACO3;          // BCB - nope, fail.
                }
            } else {
                caco3 = BAD_CACO3;            // BCB - 2 band failed out of the low side of the table, fail
            }
        }
    }

    // if valid value
    if (caco3 > BAD_CACO3) {
        caco3 = MAX(caco3, caco3min);
    }

    return (caco3);
}

/* --------------------------------------------------------------------- */
/* calcite_ci2() - calcium carbonate concentration - CI2    algorithm    */
/*                reference: doi:10.1002/2017JC013146                    */

/* --------------------------------------------------------------------- */
float calcite_ci2(l2str *l2rec, int32_t ip) {
    float caco3 = BAD_CACO3;
    int nwave = l2rec->l1rec->l1file->nbands;
    static int32_t ibRed;
    static int32_t ibGreen;
    static int firstCall = 1;

    if (firstCall) {
        ibRed = windex(667., l2rec->l1rec->l1file->fwave, nwave);
        ibGreen = windex(550., l2rec->l1rec->l1file->fwave, nwave);
        firstCall = 0;
    }

    float RrsRed = l2rec->Rrs[ip * nwave + ibRed];
    float RrsGreen = l2rec->Rrs[ip * nwave + ibGreen];
    if (RrsRed >= 0.0 && RrsGreen >= 0.0) {
        caco3 = 1.3055 * (RrsGreen - RrsRed) - 0.00188;    // SRP 09/21 changed from previous caco3 = 0.4579 * (RrsGreen - RrsRed) - 0.0006;
    }
    if (caco3<0) {
      caco3=BAD_CACO3;
    }
    return (caco3);
}

/* --------------------------------------------------------------------- */
/* calcite_ciNIR() - calcium carbonate concentration -                   */
/*                    CI748 or CI869 algorithm                           */
/*                reference: doi:10.1002/2017JC013146                    */

/* --------------------------------------------------------------------- */
float calcite_ciNIR(l2str *l2rec, int32_t ip, int32_t NIR) {
    float caco3 = BAD_CACO3;
    float CI = BAD_FLT;
    int nwave = l2rec->l1rec->l1file->nbands;
    static int32_t ibNIR;
    static int32_t ibRed;
    static int32_t ibGreen;
    static float wvlratio;
    static int firstCall = 1;

    if (firstCall) {
        if (NIR == 869) {
            ibNIR = windex(869., l2rec->l1rec->l1file->fwave, nwave);

        } else {
            ibNIR = windex(748., l2rec->l1rec->l1file->fwave, nwave);
        }
        ibRed = windex(667., l2rec->l1rec->l1file->fwave, nwave);
        ibGreen = windex(550., l2rec->l1rec->l1file->fwave, nwave);
        wvlratio = (l2rec->l1rec->l1file->fwave[ibRed] -
                l2rec->l1rec->l1file->fwave[ibGreen]) /
                (l2rec->l1rec->l1file->fwave[ibNIR] -
                l2rec->l1rec->l1file->fwave[ibGreen]);
        firstCall = 0;
    }

    float RrsNIR = l2rec->Rrs[ip * nwave + ibNIR];
    if (RrsNIR == BAD_FLT)
        RrsNIR = 0.0;
    float RrsRed = l2rec->Rrs[ip * nwave + ibRed];
    float RrsGreen = l2rec->Rrs[ip * nwave + ibGreen];
    if (RrsRed >= 0.0 && RrsGreen >= 0.0) {
        CI = RrsRed - (RrsGreen + (wvlratio * (RrsNIR - RrsGreen)));
        if (NIR == 869) {
            caco3 = -0.8013 * CI - 0.00076;
        } else {
            caco3 = -1.3764 * CI - 0.00071;
        }
    }
    if (caco3<0) {
      caco3=BAD_CACO3;
    }

    return (caco3);
}

/* ------------------------------------------------------------------- */
/* calcite() - l2_hdf_generic interface for calcite (pic)              */
/* ------------------------------------------------------------------- */
void calcite(l2str *l2rec, l2prodstr *p, float prod[]) {
    int32_t ip;

    for (ip = 0; ip < l2rec->l1rec->npix; ip++) {
        switch (p->cat_ix) {
        case CAT_calcite:
            prod[ip] = calcite_c(l2rec, ip);
            break;
        case CAT_calcite_2b:
            prod[ip] = calcite_2b(l2rec, ip);
            break;
        case CAT_calcite_3b:
            prod[ip] = calcite_3b(l2rec, ip);
            break;
        case CAT_calcite_ci2:
            prod[ip] = calcite_ci2(l2rec, ip);
            break;
        case CAT_calcite_ci748:
            prod[ip] = calcite_ciNIR(l2rec, ip,748);
            break;
        case CAT_calcite_ci869:
            prod[ip] = calcite_ciNIR(l2rec, ip,869);
            break;
        default:
            printf("Error: %s : Unknown product specifier: %d\n", __FILE__, p->cat_ix);
            exit(1);
            break;
        }

        if (prod[ip] == BAD_CACO3)
            l2rec->l1rec->flags[ip] |= PRODFAIL;
    }

    return;
}
